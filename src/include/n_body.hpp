#include <iostream>
#include <ctime>
#include "../include/math_3d.h"
#include "../include/distance.hpp"
#include "../include/pipeline.hpp"

#define TRY(command, str_error) try {                                   \
    if (command) {                                                      \
        throw std::runtime_error(str_error);                                 \
    }                                                                   \
} catch (const std::exception& e) {                                     \
    std::cerr << e.what() << std::endl;                                 \
    std::cerr << "Run with --help for more information." << std::endl;  \
    abort();                                                            \
}

#if defined(__GNUC__) || defined(__clang__)
#define likely(x)   __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)
#else
#define likely(x)   (x)
#define unlikely(x) (x)
#endif


#define NDIM 3

#define Max(a, b)  (a * (a >= b) + b * (a < b))

#define Sqr(x) ((x) * (x))
        
#define VProd(v) ((v).x * (v).y * (v).z)

#define DO_MOL(n) for (int i = 0; i < n; i++)

#define VSCopy(v2, s1, v1)      \
        (v2).x = (s1) * (v1).x, \
        (v2).y = (s1) * (v1).y, \
        (v2).z = (s1) * (v1).z 

#define PropEst(v) \
v.sum1, v.sum2

#define RandR() \
(static_cast<double>(rand()) / (double)RAND_MAX)

typedef struct {
    double val, sum1, sum2;
    void PropZero() {sum1 = 0.0f; sum2 = 0.0f;};
    void PropAccum() {sum1 += val; sum2 += Sqr(val);};
    void PropAvg(int num) {sum1 /= num; sum2 = sqrt(Max(sum2 / num - Sqr(sum1), 0.0f));};
} Prop;

typedef struct {
    Vector3f *p;
    Vector3f *f;
    Vector3f *v;
    double *m;
} DataMol;

// void init_particles ();
// void MoveBody();
void SingleStep ();
void init();

/*
#define VScale(v, s) \
(v).x *= s, \
(v).y *= s, \
(v).z *= s

#define VAdd(v1, v2, v3)  \
(v1).x = (v2).x + (v3).x, \
(v1).y = (v2).y + (v3).y, \
(v1).z = (v2).z + (v3).z

#define VSub(v1, v2, v3)  \
(v1).x = (v2).x - (v3).x, \
(v1).y = (v2).y - (v3).y, \
(v1).z = (v2).z - (v3).z

#define VMul(v1, v2, v3)  \
(v1).x = (v2).x * (v3).x, \
(v1).y = (v2).y * (v3).y, \
(v1).z = (v2).z * (v3).z

#define VDiv(v1, v2, v3)  \
(v1).x = (v2).x / (v3).x, \
(v1).y = (v2).y / (v3).y, \
(v1).z = (v2).z / (v3).z

#define VDot(v1, v2) \
((v1).x * (v2).x + (v1).y * (v2).y + (v1).z * (v2).z)

#define VSAdd(v1, v2, s3, v3) \
(v1).x = (v2).x + (s3) * (v3).x, \
(v1).y = (v2).y + (s3) * (v3).y, \
(v1).z = (v2).z + (s3) * (v3).z

#define VSet(v, sx, sy, sz) \
(v).x = sx, \
(v).y = sy, \
(v).z = sz

#define VSetAll(v, s) VSet (v, s, s, s)
#define VZero(v) VSetAll (v, 0) 
#define VVAdd(v1, v2) VAdd (v1, v1, v2)
#define VVSAdd(v1, s2, v2) VSAdd (v1, v1, s2, v2)
#define VLenSq(v) VDot (v, v)

#define VWrap(v, t) \
if (v.t >= 0.5 * region.t) v.t -= region.t; \
else if (v.t < -0.5 * region.t) v.t += region.t

#define VWrapAll(v) \
{                   \
VWrap (v, x);       \
VWrap (v, y);       \
VWrap (v, z);       \
}
*/

// void ApplyBoundaryCond ()
// {
//     DO_MOL(nMol) {
//         if (unlikely(Mol.p[i].x < -region.x)) Mol.p[i].x = fmod(region.x * 2 - Mol.p[i].x, 2 * region.x) * 0.5 + 0.5 * region.x;
//         if (unlikely(Mol.p[i].x > region.x)) Mol.p[i].x = fmod(-region.x * 2 + Mol.p[i].x, 2 * region.x) * 0.5 - 0.5 * region.x;

//         if (unlikely(Mol.p[i].y < -region.y)) Mol.p[i].y = fmod(region.y * 2 - Mol.p[i].y, 2 * region.y) * 0.5 + 0.5 * region.y;
//         if (unlikely(Mol.p[i].y > region.y)) Mol.p[i].y = fmod(-region.y * 2 + Mol.p[i].y, 2 * region.y) * 0.5 - 0.5 * region.y;

//         if (unlikely(Mol.p[i].z < -region.z)) Mol.p[i].z = fmod(region.z * 2 - Mol.p[i].z, 2 * region.z) * 0.5 + 0.5 * region.z;
//         if (unlikely(Mol.p[i].z > region.z)) Mol.p[i].z = fmod(-region.z * 2 + Mol.p[i].z, 2 * region.z) * 0.5 - 0.5 * region.z;
//     }
// }

// void MoveParticles(double dt)
// {
//     DO_MOL(nMol) {
//         Vector3f dv(
//             Mol.f[i].x / Mol.m[i] * dt,
//             Mol.f[i].y / Mol.m[i] * dt,
//             Mol.f[i].z / Mol.m[i] * dt
//         );
        
//         Vector3f dp(
//             (Mol.v[i].x + dv.x * 0.5) * dt,
//             (Mol.v[i].y + dv.y * 0.5) * dt,
//             (Mol.v[i].z + dv.z * 0.5) * dt
//         );

//         Mol.v[i] = Mol.v[i].VAdd(dv);
//         Mol.p[i] = Mol.p[i].VAdd(dp);

//         if (unlikely(Mol.p[i].x < -region.x)) Mol.p[i].x = fmod(region.x * 2 - Mol.p[i].x, 2 * region.x) * 0.5 + 0.5 * region.x;
//         if (unlikely(Mol.p[i].x > region.x)) Mol.p[i].x = fmod(-region.x * 2 + Mol.p[i].x, 2 * region.x) * 0.5 - 0.5 * region.x;

//         if (unlikely(Mol.p[i].y < -region.y)) Mol.p[i].y = fmod(region.y * 2 - Mol.p[i].y, 2 * region.y) * 0.5 + 0.5 * region.y;
//         if (unlikely(Mol.p[i].y > region.y)) Mol.p[i].y = fmod(-region.y * 2 + Mol.p[i].y, 2 * region.y) * 0.5 - 0.5 * region.y;

//         if (unlikely(Mol.p[i].z < -region.z)) Mol.p[i].z = fmod(region.z * 2 - Mol.p[i].z, 2 * region.z) * 0.5 + 0.5 * region.z;
//         if (unlikely(Mol.p[i].z > region.z)) Mol.p[i].z = fmod(-region.z * 2 + Mol.p[i].z, 2 * region.z) * 0.5 - 0.5 * region.z;


//         // Mol.p[i].x = (Mol.p[i].x >= -region.x) * Mol.p[i].x + (Mol.p[i].x < -region.x) * (2 * region.x - Mol.p[i].x);  
//         // Mol.p[i].x = (Mol.p[i].x <= region.x)  * Mol.p[i].x + (Mol.p[i].x > region.x) * (-2 * region.x + Mol.p[i].x); 
//         // Mol.p[i].y = (Mol.p[i].y >= -region.y) * Mol.p[i].y + (Mol.p[i].y < -region.y) * (2 * region.y - Mol.p[i].y); 
//         // Mol.p[i].y = (Mol.p[i].y <= region.y)  * Mol.p[i].y + (Mol.p[i].y > region.y) * (-2 * region.y + Mol.p[i].y); 
//         // Mol.p[i].z = (Mol.p[i].z >= -region.z) * Mol.p[i].z + (Mol.p[i].z < -region.z) * (2 * region.z - Mol.p[i].z); 
//         // Mol.p[i].z = (Mol.p[i].z <= region.z)  * Mol.p[i].z + (Mol.p[i].z > region.z) * (-2 * region.z + Mol.p[i].z); 

//         // Mol.f[i].x = Mol.f[i].y = Mol.f[i].z = 0;

//         distances[i].index = i;
//         distances[i].dist = Mol.p[i].Distance(pipeline.camera.Params.WorldPos);
//     }
// }