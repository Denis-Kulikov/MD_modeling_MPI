#include <ctime>
#include "../include/math_3d.h"
#include "../include/distance.hpp"
#include "../include/pipeline.hpp"

// const double G = 6.67e-11;
const double G = 10;

extern int n;
extern Vector3f *p;
extern Vector3f *f;
extern Vector3f *v;
extern double *m;
extern struct distance_by_index *distances;
extern Pipeline pipeline;

extern double width_space;
extern double height_space;
extern double length_space;

void move_particles(double dt)
{
    double boundary_offset = 0.001;

    for (int i = 0; i < n; i++) {
        Vector3f dv(
            f[i].x / m[i] * dt,
            f[i].y / m[i] * dt,
            f[i].z / m[i] * dt
        );
        
        Vector3f dp(
            (v[i].x + dv.x * 0.5) * dt,
            (v[i].y + dv.y * 0.5) * dt,
            (v[i].z + dv.z * 0.5) * dt
        );

        if (p[i].x < -width_space) {
            p[i].x = -width_space + boundary_offset; 
            v[i].x = -v[i].x; 
        }
        if (p[i].x > width_space) {
            p[i].x = width_space - boundary_offset;
            v[i].x = -v[i].x;
        }

        if (p[i].y < -height_space) {
            p[i].y = -height_space + boundary_offset;
            v[i].y = -v[i].y;
        }
        if (p[i].y > height_space) {
            p[i].y = height_space - boundary_offset;
            v[i].y = -v[i].y;
        }

        if (p[i].z < -length_space) {
            p[i].z = -length_space + boundary_offset;
            v[i].z = -v[i].z;
        }
        if (p[i].z > length_space) {
            p[i].z = length_space - boundary_offset;
            v[i].z = -v[i].z;
        }


        v[i] = VectorAdd(v[i], dv);
        p[i] = VectorAdd(p[i], dp);

        f[i].x = f[i].y = f[i].z = 0;

        distances[i].index = i;
        distances[i].dist = p[i].distance(pipeline.camera.Params.WorldPos);
    }
}

void calculate_forces(double size)
{
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            Vector3f dir = VectorSubtract(p[j], p[i]);
            double dist = p[i].distance(p[j]);
            double mag = (G * m[i] * m[j]) / powf(dist, 2);

            if (dist < size * 2) {
                Vector3f normal = VectorNormalize(VectorSubtract(p[i], p[j]));

                double v1n = VectorDot(v[i], normal);
                double v2n = VectorDot(v[j], normal);

                v[i] = VectorSubtract(v[i], VectorScale(normal, 2.0f * v1n));
                v[j] = VectorSubtract(v[j], VectorScale(normal, 2.0f * v2n));
            }

            f[i] = VectorAdd(f[i], VectorScale(VectorNormalize(dir), mag));
            f[j] = VectorSubtract(f[j], VectorScale(VectorNormalize(dir), mag));
        }
    }
}

void init_partiecle ()
{
    srand(time(NULL));
    
    for (int i = 0; i < n; i++) {
        p[i].x = (i % 2 == 0 ? -1 : 1) * (i / 2 % 2 == 0 ? 1 : -1) * width_space / 2.0;
        p[i].y = (i / 4 % 2 == 0 ? 1 : -1) * height_space / 2.0;
        p[i].z = (i / 8 % 2 == 0 ? 1 : -1) * length_space / 2.0;

        p[i].x += (rand() / (double)RAND_MAX - 0.5) * 2 / 2 * width_space;
        p[i].y += (rand() / (double)RAND_MAX - 0.5) * 2 / 2 * height_space;
        p[i].z += (rand() / (double)RAND_MAX - 0.5) * 2 / 2 * length_space;

        v[i].x = rand() / (double)RAND_MAX - G;
        v[i].y = rand() / (double)RAND_MAX - G;
        v[i].z = rand() / (double)RAND_MAX - G;

        m[i] = rand() / (double)RAND_MAX * 0.8 + 0.2;
        
        f[i].x = f[i].y = f[i].z = 0;

        distances[i].index = i;
        distances[i].dist = p[i].distance(pipeline.camera.Params.WorldPos);
    }
}

void move_body(double size)
{
    double dt = 1e-5;
    calculate_forces(size); 
    move_particles(dt);
}
