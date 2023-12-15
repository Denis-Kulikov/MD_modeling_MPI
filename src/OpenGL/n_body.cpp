#include "../include/n_body.hpp"

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

const double epsilon = 5.0;

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

double lj_potential(double r, double sigma) {
    if (r == 0.0) {
        return 0.0;
    }
    double term1 = 4.0 * epsilon * (pow(sigma / r, 12));
    double term2 = 4.0 * epsilon * (pow(sigma / r, 6));
    return term1 - term2;
}

void calculate_forces(double radius)
{
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            Vector3f dir = VectorSubtract(p[j], p[i]);
            double dist = p[i].distance(p[j]);

            if (dist < radius * 2) {
                Vector3f normal = VectorNormalize(VectorSubtract(p[i], p[j]));

                double v1n = VectorDot(v[i], normal);
                double v2n = VectorDot(v[j], normal);

                v[i] = VectorSubtract(v[i], VectorScale(normal, 2.0f * v1n));
                v[j] = VectorSubtract(v[j], VectorScale(normal, 2.0f * v2n));
            }

            double lj_energy = lj_potential(dist, radius);

            Vector3f force = VectorScale(VectorNormalize(dir), lj_energy);

            f[i] = VectorAdd(f[i], force);
            f[j] = VectorSubtract(f[j], force);
        }
    }
}

void init_partiecle ()
{
    srand(time(NULL));
    
    for (int i = 0; i < n; i++) {
        p[i].x = (i % 2 == 0 ? -1 : 1) * (i / 2 % 2 == 0 ? 1 : -1) * width_space / 2.0 +
                 (rand() / (double)RAND_MAX - 0.5)  * width_space * 0.9999;
        p[i].y = (i / 4 % 2 == 0 ? 1 : -1) * height_space / 2.0 +
                 (rand() / (double)RAND_MAX - 0.5) * height_space * 0.9999;
        p[i].z = (i / 8 % 2 == 0 ? 1 : -1) * length_space / 2.0 + 
                 (rand() / (double)RAND_MAX - 0.5) * length_space * 0.9999;

        v[i].x = 0.05 * (i % 2 == 0 ? -1 : 1) * (i / 2 % 2 == 0 ? 1 : -1) * rand() / (double)RAND_MAX;
        v[i].y = 0.05 * (i / 4 % 2 == 0 ? 1 : -1) * rand() / (double)RAND_MAX;
        v[i].z = 0.05 * (i / 8 % 2 == 0 ? 1 : -1) * rand() / (double)RAND_MAX;

        m[i] = rand() / (double)RAND_MAX * 0.8 + 0.2;
        
        f[i].x = f[i].y = f[i].z = 0;

        distances[i].index = i;
        distances[i].dist = p[i].distance(pipeline.camera.Params.WorldPos);
    }
}

void move_body(double radius)
{
    double dt = 1e-2;
    calculate_forces(radius); 
    move_particles(dt);
}
