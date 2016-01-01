// small-radiosity by Magnus Burenius
#include <algorithm>
#include <functional>
#include <numeric>
#include <random>
#include <valarray>
#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;
// Fast settings:
const auto S            = 20;		    // Number of patches per rect is S*S.			// SCENE
const auto PHOTONS      = 10000;	    // 300000;//30000;	// Number of particles used for each color component	// COMPUTE RADIOSITY
const auto DRAW_SAMPLES = 20;		    // Number of gaussian samples used for drawing.	// RENDERING
const auto WIDTH        = 250.0;		// Width of image		// RENDERING
const auto HEIGHT       = 250.0;		// Height of image		// RENDERING
// Slow settings:
//const auto S            = 64;//30;    // Number of patches per rect is S*S.
//const auto PHOTONS      = 100000000;  // Number of particles used for each color component
//const auto DRAW_SAMPLES = 100;        // Number of gaussian samples used for drawing.
//const auto WIDTH        = 4096.0;     // Width of image
//const auto HEIGHT       = 4096.0;     // Height of image
auto g = bind(      normal_distribution<double>(), mt19937());
auto u = bind(uniform_real_distribution<double>(), mt19937());
// Definition of vector type:
using vec = valarray<double>;
double dot(vec a, vec b)
{
    return (a * b).sum();
}
vec normalized(vec v)
{
    return v / sqrt(dot(v, v));
}
vec crossProduct(vec x, vec y)
{
    return vec{x[1]*y[2]-x[2]*y[1], x[2]*y[0]-x[0]*y[2], x[0]*y[1]-x[1]*y[0]};
}
// Used to describe the geometrical properties of a rectangular surface:
struct Rect
{
    vec p;      // A corner in the rectangle.
    vec x;      // Vector going from p to the right nearby corner.
    vec y;      // Vector going from p to the top nearby corner.
    vec n;      // Normal of the plane. The cross product of x and y.
    vec xn;     // normalized(x).
    vec yn;     // normalized(y).
    vec nn;     // normalized(n).
    double a;   // Area of rectangle. Given by length of x cross y.
    Rect(vec P, vec X, vec Y) : p{P}, x{X}, y{Y}, n{crossProduct(x, y)},
        xn{normalized(x)}, yn{normalized(y)}, nn{normalized(n)}, a{sqrt(dot(n, n))} {}
};
const auto L = 555.0; // Length of Cornell Box side.
struct Scene // Define the Cornell box:
{
    static const auto NUM_RECTANGLES = 15;
    static const auto NUM_PATCHES = NUM_RECTANGLES * S * S;

    Rect rectangles[NUM_RECTANGLES] = {
        Rect{ vec{ 0, 0, 0 }, vec{  0, 0, L }, vec{ L, 0, 0 } },    // Floor
        Rect{ vec{ 0, L, 0 }, vec{  L, 0, 0 }, vec{ 0, 0, L } },    // Ceiling
        Rect{ vec{ L, 0, L }, vec{ -L, 0, 0 }, vec{ 0, L, 0 } },	    // Back wall
        Rect{ vec{ 0, 0, 0 }, vec{  0, L, 0 }, vec{ 0, 0, L } },    // Right wall
        Rect{ vec{ L, 0, 0 }, vec{  0, 0, L }, vec{ 0, L, 0 } },    // Left wall
        // Tall block:
        Rect{ vec{ 314, 330, 454 }, vec{  158, 0,  -49 }, vec{ -49,   0, -158 } },
        Rect{ vec{ 314,   0, 454 }, vec{  158, 0,  -49 }, vec{   0, 330,    0 } },
        Rect{ vec{ 423,   0, 247 }, vec{ -158, 0,   49 }, vec{   0, 330,    0 } },
        Rect{ vec{ 265,   0, 296 }, vec{   49, 0,  158 }, vec{   0, 330,    0 } },
        Rect{ vec{ 472,   0, 405 }, vec{  -49, 0, -158 }, vec{   0, 330,    0 } },
        // Short block:
        Rect{ vec{  81, 165, 223 }, vec{  158, 0,   49 }, vec{  49,   0, -158 } },
        Rect{ vec{  81,   0, 223 }, vec{  158, 0,   49 }, vec{   0, 165,    0 } },
        Rect{ vec{ 288,   0, 114 }, vec{ -158, 0,  -49 }, vec{   0, 165,    0 } },
        Rect{ vec{ 130,   0,  65 }, vec{  -49, 0,  158 }, vec{   0, 165,    0 } },
        Rect{ vec{ 239,   0, 272 }, vec{   49, 0, -158 }, vec{   0, 165,    0 } } };

    // The color quantities for each patch:
    // TODO: rename?
    valarray<vec> R{ .75 * vec{ 1, 1, 1 }, NUM_PATCHES }; // Reflectance for each patch
    valarray<vec> B{       vec{ 0, 0, 0 }, NUM_PATCHES }; // Radiosity for each patch

    vec lightPos   =    L * vec{ .5, .8, .5 };
    vec lightPower = 1e11 * vec{  1,  1,  1 };

    Scene()
    {
        fill(&R[3 * S * S], &R[4 * S * S], vec{ .25, .75, .25 }); // Color right wall
        fill(&R[4 * S * S], &R[5 * S * S], vec{ .75, .25, .25 }); // Color left wall
    }
};
// Information about an intersection between a ray and a Rect:
struct Intersection
{
    double distance;
    double u; // 0-1 coordinate of the intersection within the rectangle.
    double v; // 0-1 coordinate of the intersection within the rectangle.
    vec position;
    int rectangleIndex;
    Intersection() : distance{numeric_limits<double>::max()} {}
    operator bool() const { return distance < numeric_limits<double>::max(); }
};
// Compute the first intersection along a ray:
Intersection ComputeIntersection(const Scene& scene, vec start, vec dir)
{
    const auto e = 0.0001;
    start += e * dir;
    dir = normalized(dir);
    auto firstIntersection = Intersection();
    for (auto r = 0; r < Scene::NUM_RECTANGLES; ++r)
    {
        const auto& rectangle = scene.rectangles[r];
        auto i = Intersection();
        i.distance = dot(rectangle.p - start, rectangle.nn) / dot(dir, rectangle.nn);
        if (i.distance < 0 || firstIntersection.distance < i.distance)
            continue;
        i.rectangleIndex = r;
        i.position = start + i.distance * dir;
        const auto p = i.position - rectangle.p;
        i.u = dot(p, rectangle.x) / dot(rectangle.x, rectangle.x);
        i.v = dot(p, rectangle.y) / dot(rectangle.y, rectangle.y);
        if (0 - e < i.u && i.u < 1 + e && 0 - e < i.v && i.v < 1 + e)
            firstIntersection = i;
    }
    return firstIntersection;
}
int RectToPatch(int rect_index, int rect_patch_u, int rect_path_v)
{
    return rect_patch_u + rect_path_v * S + rect_index * S * S;
}
// Get patch index of an intersection, with an optional u and v offset:
int Patch(const Intersection& i, double du = 0, double dv = 0)
{
    const auto rect_path_u = max(min(int((i.u + du) * S), S - 1), 0);
    const auto rect_path_v = max(min(int((i.v + dv) * S), S - 1), 0);
    return RectToPatch(i.rectangleIndex, rect_path_u, rect_path_v);
}
// Get an "interpolated" value for the radiosity using gaussian samples:
vec Radiosity(const Scene& scene, const Intersection& i)
{
    auto r = vec{0, 0, 0};
    if (!i)
        return r;
    //return B[Patch(i)]; // No interpolation
    // Gaussian filter:
    const auto DRAW_STD = 0.5 / S; // STD of gaussian filter.
    for (auto s = 0; s < DRAW_SAMPLES; ++s)
        r += scene.B[Patch(i, DRAW_STD * g(), DRAW_STD * g())];
    return r / double(DRAW_SAMPLES);
}
vec RandomDiffuseReflectionDirection(vec tangent1, vec tangent2, vec normal)
{
    const auto a = u();
    const auto dir = sqrt(a) * normalized(vec{g(), g(), 0}) + vec{0, 0, sqrt(1 - a)};
    return tangent1 * dir[0] + tangent2 * dir[1] + normal * dir[2];
}
vec RandomDirection()
{
    return normalized(vec{ g(), g(), g() });
}
bool PhotonIsAbsorbed(double reflectance)
{
    return u() > reflectance;
}
const auto NUM_COLOR_CHANNELS = 3;
// Compute the radiosity of the scene by bouncing around photons.
void ComputeRadiosity(Scene& scene)
{
    for (auto c = 0; c < NUM_COLOR_CHANNELS; ++c)
    {
        for (auto p = 0; p < PHOTONS; ++p)
        {
            auto ray_dir = RandomDirection();
            auto ray_start = scene.lightPos;
            auto i = ComputeIntersection(scene, ray_start, ray_dir);
            while (i && !PhotonIsAbsorbed(scene.R[Patch(i)][c]))
            {
                const auto& r = scene.rectangles[i.rectangleIndex];
                scene.B[Patch(i)][c] += scene.lightPower[c] / PHOTONS * S * S / r.a;
                ray_dir = RandomDiffuseReflectionDirection(r.xn, r.yn, r.nn);
                ray_start = i.position;
                i = ComputeIntersection(scene, ray_start, ray_dir);
            }
        }
    }
}
// Do gamma correction and truncation of color:
int ScreenColor(double c)
{
    return min(int(pow(c, 1 / 2.2)), 255);
}
// Render the image to a ppm-file:
void Render(const Scene& scene, vec cameraPos, double focalLength, const char* filename)
{
    ofstream file(filename);
    file << "P3\n" << int(WIDTH) << " " << int(HEIGHT) << "\n255\n";
    for (auto y = 0.0; y < HEIGHT; ++y)
    {
        for (auto x = 0.0; x < WIDTH; ++x)
        {
            const auto rayDir = vec{WIDTH / 2 - x, HEIGHT / 2 - y, focalLength};
            const auto color = Radiosity(scene, ComputeIntersection(scene, cameraPos, rayDir));
            for (auto c = 0; c < NUM_COLOR_CHANNELS; ++c)
                file << ScreenColor(color[c]) << " ";
        }
    }
}
void SaveLightmaps(const Scene& scene, const char* fileNameStart)
{
    for (auto r = 0; r < Scene::NUM_RECTANGLES; ++r)
    {
        auto ss = stringstream();
        ss << fileNameStart << r << ".ppm";
        auto file = ofstream(ss.str());
        file << "P3" << endl << S << " " << S << endl << 255 << endl;
        for (auto y = 0; y < S; ++y)
        {
            for (auto x = 0; x < S; ++x)
            {
                const auto color = scene.B[RectToPatch(r, x, y)];
                for (auto c = 0; c < NUM_COLOR_CHANNELS; ++c)
                    file << ScreenColor(color[c]) << " ";
            }
        }
    }
}
// Here we go:
int main()
{
    auto scene = Scene();
    cout << "Computing radiosity." << endl;
    ComputeRadiosity(scene);
    cout << "Saving lightmaps to files." << endl;
    SaveLightmaps(scene, "lightmap");
    cout << "Rendering image to file." << endl;
    const auto focalLength = 1.4 * WIDTH;
    const auto cameraPos = L * vec{.5, .5, -1.4};
    Render(scene, cameraPos, focalLength, "image.ppm");
}


//void LoadLightmaps()
//{
//	for(size_t q=0; q<Q;++q)
//	{
//		stringstream ss;
//		ss << "lightmaps30/lightmap" << q << ".ppm";
//		ifstream f( ss.str() );
//		string temp;
//		size_t w,h,color;
//		f >> temp >> w >> h >> color;
//		for( size_t y=0; y<S; ++y )
//		{
//			for( size_t x=0; x<S; ++x )
//			{
//				vec color = vec{0,0,0);
//				for( int c=0; c<3; ++c )
//					f >> color[c];
//				B[QuadToPatch(q,x,y)] = 700.*color; // N.B. truncated over 255
//			}
//		}
//	}
//}



//void LoadLightmaps()
//{
//	for(size_t q=0; q<Q;++q)
//	{
//		stringstream ss;
//		ss << "lightmap" << q << ".ppm";
//		ifstream f( ss.str() );
//		string temp;
//		size_t w,h,color;
//		f >> temp >> w >> h >> color;
//		for( size_t y=0; y<S; ++y )
//		{
//			for( size_t x=0; x<S; ++x )
//			{
//				vec color = vec{0,0,0);
//				for( int c=0; c<3; ++c )
//					f >> color[c];
//				radiosity[QuadToPatch(q,x,y)] = color; // N.B. truncated over 255
//			}
//		}
//	}
//}
