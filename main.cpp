// small-radiosity by Magnus Burenius
#include <algorithm>
#include <functional>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <valarray>
using namespace std;
// Fast settings:
const auto S            = 20;		    // Number of patches per rect is S*S.
const auto PHOTONS      = 10000;	    // Number of particles used for each color component.
const auto DRAW_SAMPLES = 20;		    // Number of gaussian samples used for drawing.
const auto WIDTH        = 250.0;		// Width of image.
const auto HEIGHT       = 250.0;		// Height of image.
// Slow settings:
//const auto S            = 64;         // Number of patches per rect is S*S.
//const auto PHOTONS      = 100000000;  // Number of particles used for each color component.
//const auto DRAW_SAMPLES = 100;        // Number of gaussian samples used for drawing.
//const auto WIDTH        = 4096.0;     // Width of image.
//const auto HEIGHT       = 4096.0;     // Height of image.
auto g = bind(      normal_distribution<double>(), mt19937());
auto u = bind(uniform_real_distribution<double>(), mt19937());
// Definition of vector type:
using vec = valarray<double>;
double dot(vec a, vec b)
{
    return (a * b).sum();
}
double norm(vec a)
{
    return sqrt(dot(a, a));
}
vec normalized(vec a)
{
    return a / norm(a);
}
vec crossProduct(vec a, vec b)
{
    return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
}
const auto NUM_COLOR_CHANNELS = 3;
const auto BLACK = vec{ 0.00, 0.00, 0.00 };
const auto WHITE = vec{ 0.75, 0.75, 0.75 };
const auto RED   = vec{ 0.75, 0.25, 0.25 };
const auto GREEN = vec{ 0.25, 0.75, 0.25 };
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
        xn{normalized(x)}, yn{normalized(y)}, nn{normalized(n)}, a{ norm(n)} {}
};
// Define the scene which is a cubic room with two blocks inside:
const auto L = 555.0; // Length of the side of the room.
struct Scene
{
    static const auto NUM_RECTANGLES = 15;
    static const auto NUM_PATCHES = NUM_RECTANGLES * S * S;

    static int rectToPatch(int rect_index, int rect_patch_u, int rect_path_v)
    {
        return rect_patch_u + rect_path_v * S + rect_index * S * S;
    }

    Rect rectangles[NUM_RECTANGLES] = {
        Rect{ { 0, 0, 0 }, {  0, 0, L }, { L, 0, 0 } }, // Floor
        Rect{ { 0, L, 0 }, {  L, 0, 0 }, { 0, 0, L } }, // Ceiling
        Rect{ { L, 0, L }, { -L, 0, 0 }, { 0, L, 0 } }, // Back wall
        Rect{ { 0, 0, 0 }, {  0, L, 0 }, { 0, 0, L } }, // Right wall
        Rect{ { L, 0, 0 }, {  0, 0, L }, { 0, L, 0 } }, // Left wall
        // Tall block:
        Rect{ { 314, 330, 454 }, {  158, 0,  -49 }, { -49,   0, -158 } },
        Rect{ { 314,   0, 454 }, {  158, 0,  -49 }, {   0, 330,    0 } },
        Rect{ { 423,   0, 247 }, { -158, 0,   49 }, {   0, 330,    0 } },
        Rect{ { 265,   0, 296 }, {   49, 0,  158 }, {   0, 330,    0 } },
        Rect{ { 472,   0, 405 }, {  -49, 0, -158 }, {   0, 330,    0 } },
        // Short block:
        Rect{ {  81, 165, 223 }, {  158, 0,   49 }, {  49,   0, -158 } },
        Rect{ {  81,   0, 223 }, {  158, 0,   49 }, {   0, 165,    0 } },
        Rect{ { 288,   0, 114 }, { -158, 0,  -49 }, {   0, 165,    0 } },
        Rect{ { 130,   0,  65 }, {  -49, 0,  158 }, {   0, 165,    0 } },
        Rect{ { 239,   0, 272 }, {   49, 0, -158 }, {   0, 165,    0 } } };

    // The color quantities for each patch:
    valarray<vec> reflectance{ WHITE, NUM_PATCHES }; // Reflectance for each patch
    valarray<vec>   radiosity{ BLACK, NUM_PATCHES }; // Radiosity for each patch

    Scene()
    {
        fill(&reflectance[rectToPatch(3, 0, 0)], &reflectance[rectToPatch(4, 0, 0)], GREEN); // Color right wall.
        fill(&reflectance[rectToPatch(4, 0, 0)], &reflectance[rectToPatch(5, 0, 0)], RED); // Color left wall.
    }
};
struct Light
{
    vec position = L * vec{ .5, .8, .5 };
    vec power = 1e11 * vec{ 1,  1,  1 };
};
struct Camera
{
    vec position = L * vec{ .5, .5, -1.4 };
    double focalLength = 1.4 * WIDTH;
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
// Find the first intersection along a ray:
Intersection findIntersection(const Scene& scene, vec start, vec dir)
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
        const vec p = i.position - rectangle.p;
        i.u = dot(p, rectangle.x) / dot(rectangle.x, rectangle.x);
        i.v = dot(p, rectangle.y) / dot(rectangle.y, rectangle.y);
        if (0 - e < i.u && i.u < 1 + e && 0 - e < i.v && i.v < 1 + e)
            firstIntersection = i;
    }
    return firstIntersection;
}
// Get patch index of an intersection, with an optional u and v offset:
int patchIndex(const Intersection& i, double du = 0, double dv = 0)
{
    const auto rect_path_u = max(min(int((i.u + du) * S), S - 1), 0);
    const auto rect_path_v = max(min(int((i.v + dv) * S), S - 1), 0);
    return Scene::rectToPatch(i.rectangleIndex, rect_path_u, rect_path_v);
}
// Get an "interpolated" value for the radiosity using gaussian samples:
vec sampleRadiosity(const Scene& scene, const Intersection& i)
{
    auto r = BLACK;
    if (!i)
        return r;
    //return B[patchIndex(i)]; // No interpolation.
    // Gaussian filter:
    const auto DRAW_STD = 0.5 / S; // STD of gaussian filter.
    for (auto s = 0; s < DRAW_SAMPLES; ++s)
        r += scene.radiosity[patchIndex(i, DRAW_STD * g(), DRAW_STD * g())];
    return r / double(DRAW_SAMPLES);
}
vec randomDirection()
{
    return normalized({ g(), g(), g() });
}
vec randomDiffuseReflectionDirection(vec tangent1, vec tangent2, vec normal)
{
    const auto a = u();
    const vec dir = sqrt(a) * normalized({g(), g(), 0}) + vec{0, 0, sqrt(1 - a)};
    return tangent1 * dir[0] + tangent2 * dir[1] + normal * dir[2];
}
bool isPhotonAbsorbed(double reflectance)
{
    return u() > reflectance;
}
// Compute the radiosity of the scene by bouncing around photons.
void illuminateScene(Scene& scene, const Light& light)
{
    for (auto c = 0; c < NUM_COLOR_CHANNELS; ++c)
    {
        for (auto p = 0; p < PHOTONS; ++p)
        {
            auto ray_dir = randomDirection();
            auto ray_start = light.position;
            auto i = findIntersection(scene, ray_start, ray_dir);
            while (i && !isPhotonAbsorbed(scene.reflectance[patchIndex(i)][c]))
            {
                const auto& r = scene.rectangles[i.rectangleIndex];
                scene.radiosity[patchIndex(i)][c] += light.power[c] / PHOTONS * S * S / r.a;
                ray_dir = randomDiffuseReflectionDirection(r.xn, r.yn, r.nn);
                ray_start = i.position;
                i = findIntersection(scene, ray_start, ray_dir);
            }
        }
    }
}
// Do gamma correction and truncation of color:
int screenColor(double c)
{
    return min(int(pow(c, 1 / 2.2)), 255);
}
// Render the image to a ppm-file:
void renderCameraImage(const Scene& scene, const Camera& camera, const char* filename)
{
    ofstream file(filename);
    file << "P3" << endl << int(WIDTH) << " " << int(HEIGHT) << endl << "255" << endl;
    for (auto y = 0.0; y < HEIGHT; ++y)
    {
        for (auto x = 0.0; x < WIDTH; ++x)
        {
            const auto rayDir = vec{WIDTH / 2 - x, HEIGHT / 2 - y, camera.focalLength};
            const auto intersection = findIntersection(scene, camera.position, rayDir);
            const auto radiosity = sampleRadiosity(scene, intersection);
            for (auto c = 0; c < NUM_COLOR_CHANNELS; ++c)
                file << screenColor(radiosity[c]) << " ";
        }
    }
}
void saveLightmaps(const Scene& scene)
{
    for (auto r = 0; r < Scene::NUM_RECTANGLES; ++r)
    {
        const auto filename = "lightmap" + std::to_string(r) +  ".ppm";
        ofstream file(filename);
        file << "P3" << endl << S << " " << S << endl << 255 << endl;
        for (auto y = 0; y < S; ++y)
        {
            for (auto x = 0; x < S; ++x)
            {
                const auto radiosity = scene.radiosity[scene.rectToPatch(r, x, y)];
                for (auto c = 0; c < NUM_COLOR_CHANNELS; ++c)
                    file << screenColor(radiosity[c]) << " ";
            }
        }
    }
}
int main()
{
    auto scene = Scene{};
    auto light = Light{};
    auto camera = Camera{};
    cout << "Illuminating the scene." << endl;
    illuminateScene(scene, light);
    cout << "Saving lightmaps to files." << endl;
    saveLightmaps(scene);
    cout << "Rendering camera image to file." << endl;
    renderCameraImage(scene, camera, "cameraImage.ppm");
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
