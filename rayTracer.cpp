/*
Author:  Anton Medvdev, amedvedev2013@my.fit.edu
Course:  CSE 5280, Spring 2017
Project: Proj 03, Ray Tracing
*/

// imports 
#include <fstream>
#include <cmath>

//Vector Struct
struct Vector {
	
  double x,y,z;
  
  Vector(double x, double y, double z) : x(x), y(y), z(z) {}
  
  //Vector Operations
  Vector operator + (const Vector& v) const { return Vector(x+v.x, y+v.y, z+v.z); }
  Vector operator - (const Vector& v) const { return Vector(x-v.x, y-v.y, z-v.z); }
  Vector operator * (double d) const { return Vector(x*d, y*d, z*d); }
  Vector operator / (double d) const { return Vector(x/d, y/d, z/d); }
  
  //norm
  Vector norm() const {
    double mg = sqrt(x*x + y*y + z*z);
    return Vector(x/mg,y/mg,z/mg);
  }
};

//Ray Struct
struct Ray {
  Vector o,d;
  Ray(const Vector& o, const Vector& d) : o(o), d(d) {}
};

//Vector dot product
inline double dot(const Vector& a, const Vector& b) {
  return (a.x*b.x + a.y*b.y + a.z*b.z);
}

//Sphere
struct Sphere {
  Vector c;
  double r;
  Sphere(const Vector& c, double r) : c(c), r(r) {}
  Vector getNormal(const Vector& pi) const { return (pi - c) / r; }
  bool intersect(const Ray& ray, double &t) const {
    const Vector o = ray.o;
    const Vector d = ray.d;
    const Vector oc = o - c;
    const double b = 2 * dot(oc, d);
    const double c = dot(oc, oc) - r*r;
    double disc = b*b - 4 * c;
    if (disc < 1e-4) return false;
    disc = sqrt(disc);
    const double t0 = -b - disc;
    const double t1 = -b + disc;
    t = (t0 < t1) ? t0 : t1;
    return true;
  }
};


void clamp255(Vector& col) {
  col.x = (col.x > 255) ? 255 : (col.x < 0) ? 0 : col.x;
  col.y = (col.y > 255) ? 255 : (col.y < 0) ? 0 : col.y;
  col.z = (col.z > 255) ? 255 : (col.z < 0) ? 0 : col.z;
}

int main() {

  //Height and Width
  const int H = 500;
  const int W = 500;

 
  const Vector white(255, 255, 255);
  const Vector greyShade(192, 192, 192);
  const Vector blueShade(65,36,250);

  const Sphere sphere(Vector(W*0.5, H*0.5, 50), 50);
 
  const Sphere light(Vector(0, 0, 50), 1);

  std::ofstream out("out.ppm");
  out << "P3\n" << W << ' ' << H << ' ' << "255\n";

  double t;
  
  Vector pix_col(greyShade);

  // For each pixel 
  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      
      const Ray ray(Vector(x,y,0),Vector(0,0,1));
      
      pix_col = greyShade;
           
      if ((sphere.intersect(ray, t))) {
      
        const Vector pi = ray.o + ray.d*t;
        const Vector L = light.c - pi;
        
        const Vector N = sphere.getNormal(pi);
        const double dt = dot(L.norm(), N.norm());

		//
        pix_col = (blueShade + white*dt) * 0.5;
        clamp255(pix_col);
        
      }
      out << (int)pix_col.x << ' '
          << (int)pix_col.y << ' '
          << (int)pix_col.z << '\n';
    }
  }
}