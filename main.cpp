// SmallRadiosity by Magnus Burenius
#include <iostream>	// Not needed for minimal
#include <fstream>
#include <math.h>
#include <valarray>
#include <algorithm>
#include <numeric>
#include <random>
#include <functional>
#include <sstream>
#include <vector>
using namespace std;
//// Fast settings:
const double WIDTH			= 250;		// Width of image		// RENDERING
const double HEIGHT			= 250;		// Height of image		// RENDERING
const size_t PHOTONS		= 10000;	// 300000;//30000;	// Number of particles used for each color component	// COMPUTE RADIOSITY
const size_t DRAW_SAMPLES	= 20;		// Number of gaussian samples used for drawing.	// RENDERING
const size_t RECTS			= 15;		// Number of rectangles.						// SCENE
const size_t S				= 20;		// Number of patches per rect is S*S.			// SCENE
const size_t PATCHES		= RECTS*S*S;// Number of patches in total.					// SCENE
const double DRAW_STD		= 0.5/S;	// STD of gaussian filter						// RENDERING (this could be a hidden parameter)
// Slow settings:
//const double WIDTH			= 4096;	// Width of image
//const double HEIGHT			= 4096;	// Height of image
//const size_t PHOTONS		= 100000000;// Number of particles used for each color component
//const size_t DRAW_SAMPLES	= 100;		// Number of gaussian samples used for drawing.
//const size_t RECTS			= 15;		// Number of rectangles.
//const size_t S				= 64;//30;		// Number of patches per rect is S*S.
//const size_t PATCHES		= RECTS*S*S;// Number of patches in total.
//const double DRAW_STD		= 0.5/S;	// STD of gaussian filter
auto g = bind(      normal_distribution<double>(), mt19937());
auto u = bind(uniform_real_distribution<double>(), mt19937());
// Definition of vector type:
using vec = valarray<double>;
double	dot( const vec& a, const vec& b )	{ return (a*b).sum();			}
vec		normalize( const vec& v )			{ return v / sqrt( dot(v,v) );	}
// Used to describe the geometrical properties of a rectangular surface:
struct Rect
{
	vec p;		// A corner in the rectangle.
	vec x;		// Vector going from p to the right nearby corner.
	vec y;		// Vector going from p to the top nearby corner.
	vec n;		// Normal of the plane. The cross product of x and y.
	vec xn;		// normalize(x)
	vec yn;		// normalize(y)
	vec nn;		// normalize(n)
	double a;	// Area of rectangle. Given by length of x cross y.
	Rect(vec P, vec X, vec Y) : p{ P }, x{ X }, y{ Y },
		n{ vec{ x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0] } },
		xn{ normalize(x) }, yn{ normalize(y) }, nn{ normalize(n) }, a{ sqrt(dot(n, n)) } {}
};

struct Scene
{
	vector<Rect>	rectangles;
	vector<vec>		reflectance;
	vector<vec>		radiosity;
	vec				lightPos;
	vec				lighPower;
};

// The color quantities for each patch:
valarray<vec> R(.75*vec{ 1, 1, 1 }, PATCHES);	// Reflectance for each patch
valarray<vec> B(	vec{ 0, 0, 0 }, PATCHES);	// Radiosity for each patch

// Define the Cornell box:
const double L = 555;		// Length of Cornell Box side.
Rect rects[RECTS] = {
	Rect{ vec{ 0, 0, 0 }, vec{  0, 0, L }, vec{ L, 0, 0 } },	// Floor
	Rect{ vec{ 0, L, 0 }, vec{  L, 0, 0 }, vec{ 0, 0, L } },	// Ceiling
	Rect{ vec{ L, 0, L }, vec{ -L, 0, 0 }, vec{ 0, L, 0 } },	// Back wall
	Rect{ vec{ 0, 0, 0 }, vec{  0, L, 0 }, vec{ 0, 0, L } },	// Right wall
	Rect{ vec{ L, 0, 0 }, vec{  0, 0, L }, vec{ 0, L, 0 } },	// Left wall
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
// Define light:
const vec lightPos  =    L * vec{ .5, .8,   .5 };
const vec lighPower = 1e11 * vec{  1,  1,    1 };
// Information about an intersection between a ray and a Rect:
struct Intersection
{
	double dist,u,v;	// Distance to intersection and its u,v coordinates.
	vec pos;			// Intersection point.
	size_t rectIndex;	// Index of intersecting Rect.
	Intersection() : dist{numeric_limits<double>::max()} {}
	operator bool() const { return dist < numeric_limits<double>::max(); }
};
// Compute the first intersections along a ray:
Intersection ComputeIntersection( vec start, vec dir )
{
	const double e = 0.0001;
	start += e*dir;
	dir = normalize( dir );
	Intersection firstIntersection;
	//firstIntersection.dist = numeric_limits<double>::max();
	for( size_t r=0; r<RECTS; ++r )
	{	// Check quad ray intersection:
		Intersection i;
		i.dist = dot(rects[r].p-start,rects[r].n) / dot( dir, rects[r].n );
		if( i.dist < 0 || firstIntersection.dist < i.dist )
			continue;
		i.rectIndex = r;
		i.pos = start + i.dist*dir;
		vec p = i.pos - rects[r].p;
		i.u = dot( p, rects[r].x ) / dot(rects[r].x, rects[r].x );
		i.v = dot( p, rects[r].y ) / dot(rects[r].y, rects[r].y );
		if(  0-e < i.u && i.u < 1+e && 0-e < i.v && i.v < 1+e )
			firstIntersection = i;
	}
	return firstIntersection;
}
// Get patch index of an intersection, with an optional u,v offset:
size_t Patch( const Intersection& i, double du=0, double dv=0 )
{
	size_t ui = min( size_t((i.u+du)*S), S-1 );
	size_t vi = min( size_t((i.v+dv)*S), S-1 );
	return ui + vi*S + i.rectIndex*S*S;
}
// Get an "interpolated" value for the radiosity using gaussian samples:
vec Radiosity( const Intersection& i )
{
	vec r = vec{ 0, 0, 0 };
	if( !i ) //.infinity() )
		return r;
	//return B[ Patch(i) ]; // No interpolation
	// Gaussian filter:
	for( size_t s=0; s<DRAW_SAMPLES; ++s )
		r += B[ Patch(i,DRAW_STD*g(),DRAW_STD*g()) ];
	return r / double(DRAW_SAMPLES);
}
// Third component is assumed to be in the direction of the normal:
vec RandomDiffuseReflectionDir()
{
	double a = u();
	return sqrt(a)*normalize(vec{g(),g(),0}) + sqrt(1-a)*vec{0,0,1};
}

void ComputeRadiosity()
{
	for( size_t c=0; c<3; ++c)
	{
		for( size_t p=0; p<PHOTONS; ++p )	// Generate lines/rays:
		{
			// Set initial photon position and direction:
			vec start = lightPos;
			vec dir = normalize(vec{ g(), g(), g() });
			Intersection i = ComputeIntersection( start, dir );
			while( i && u() < R[Patch(i)][c] )
			{
				const Rect& r = rects[ i.rectIndex ];
				B[ Patch(i) ][c] += lighPower[c]/PHOTONS*S*S/r.a;
				start = i.pos;
				vec d = RandomDiffuseReflectionDir();
				dir = r.xn*d[0] + r.yn*d[1] + r.nn*d[2];
				i = ComputeIntersection( start, dir );
			}
		}
	}
}
// Do gamma correction and truncation of color:
int ScreenColor( double c )
{
	return min( int(pow(c,1/2.2)), 255 ); 
}
// Render the image to a ppm-file:
void Render( vec cameraPos, double focalLength, const char* filename )
{
	ofstream file( filename );
	file << "P3\n" << int(WIDTH) << " " << int(HEIGHT) << "\n255\n";
	for( double y=0; y<HEIGHT; ++y )
	{
		for( double x=0; x<WIDTH; ++x )
		{
			vec pixelDir = vec{ WIDTH / 2 - x, HEIGHT / 2 - y, focalLength };
			vec color = Radiosity( ComputeIntersection( cameraPos, pixelDir ) );
			for( size_t c=0; c<3; ++c )
				file << ScreenColor( color[c] ) << " ";
		}
	}
}
size_t RectToPatch(size_t r, size_t i, size_t j)
{
	return i + j*S + r*S*S;
}
void SaveLightmaps( const char* fileNameStart )
{
	for (size_t r = 0; r<RECTS; ++r)
	{
		stringstream ss;
		ss << fileNameStart << r << ".ppm";
		ofstream file(ss.str());
		file << "P3\n" << S << " " << S << "\n255\n";
		for (size_t y = 0; y<S; ++y)
		{
			for (size_t x = 0; x<S; ++x)
			{
				vec color = B[RectToPatch(r, x, y)];
				for (int c = 0; c<3; ++c)
					file << ScreenColor(color[c]) << " ";
			}
		}
	}
}
// Here we go:
int main()
{
	fill( &R[3*S*S], &R[4*S*S], vec{.25,.75,.25} );	// Color right wall
	fill( &R[4*S*S], &R[5*S*S], vec{.75,.25,.25} );	// Color left wall
	cout << "Computing radiosity." << endl;
	ComputeRadiosity();
	cout << "Saving lightmaps to files." << endl;
	SaveLightmaps("lightmap");
	cout << "Rendering image to file." << endl; // Rendering image from a 
	double focalLength = 1.4*WIDTH;
	vec cameraPos = L*vec{.5,.5,-1.4};
	Render( cameraPos, focalLength, "image.ppm" );
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
