#pragma once

#include <iostream>
#include <ctime>
#include <string>
#include <vector>
#include "../include/math_3d.h"

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

enum SIDE {
    RIGHT,
    LEFT,
    TOP,
    BOTTOM,
    FRONT,
    BACK
};

typedef struct {
    Vector3f *p;
    Vector3f *f;
    Vector3f *v;
    double *m;
} DataMol;

typedef struct OneMolStruct {
    Vector3f p;
    Vector3f f;
    Vector3f v;
    double m;
    OneMolStruct (Vector3f _p, Vector3f _f, Vector3f _v, double _m)
    {
        p = _p;
        f = _f;
        v = _v;
        m = _m;
    };
} OneMol;

typedef struct {
    std::vector<OneMol> *esc[6];
    std::vector<int> *n[6];
} Escapees;

typedef struct {
    double val, sum1, sum2;
    void PropZero() {sum1 = 0.0f; sum2 = 0.0f;};
    void PropAccum() {sum1 += val; sum2 += Sqr(val);};
    void PropAvg(int num) {sum1 /= num; sum2 = sqrt(Max(sum2 / num - Sqr(sum1), 0.0f));};
} Prop;

void CloseFile(int rank);
void OpenFile(int rank);
void WriteParams(int commsize);

void WriteSystem ();
void CloseSystemFile();
void OpenSystemFile();
void AccumProps (int icode);
void EvalProps ();
void GetvSum (int n);

Escapees FindEscapees(int n, const Vector3i &crank, const Vector3i &dims, const Vector3f &center);
void SingleStep (int n);
Vector3f GetCenter(const Vector3i &crank, const Vector3i &dims);
bool AllocArrays (int commsize);
void SetupJob (int rank, int n, Vector3f center);
void SetParams (const Vector3i &dims);
