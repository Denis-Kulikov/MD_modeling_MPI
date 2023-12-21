#include <iostream>
#include <ctime>
#include "../include/math_3d.h"
#include "../include/distance.hpp"
#include "../include/pipeline.hpp"

#define TRY(command, str_error) try {                                   \
    if (command) {                                                      \
        throw std::runtime_error(str_error);                            \
    }                                                                   \
} catch (const std::exception& e) {                                     \
    std::cerr << e.what() << std::endl;                                 \
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

void SingleStep (int first, int last);
bool AllocArrays ();
void SetupJob ();
void SetParams ();
void init();
