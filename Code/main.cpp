#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <algorithm>
using namespace std;

const int N = 4;
const int NS = 6;
const double EPS = 1e-12;

struct Vec3 {
    double x, y, z;

    Vec3 operator+(const Vec3& b) const { return {x+b.x, y+b.y, z+b.z}; }
    Vec3 operator-(const Vec3& b) const { return {x-b.x, y-b.y, z-b.z}; }
    Vec3 operator*(double s) const { return {x*s, y*s, z*s}; }
    Vec3 operator/(double s) const { return {x/s, y/s, z/s}; }
    Vec3& operator+=(const Vec3& b) { x+=b.x; y+=b.y; z+=b.z; return *this; }
    Vec3& operator-=(const Vec3& b) { x-=b.x; y-=b.y; z-=b.z; return *this; }
    Vec3& operator*=(double s) { x*=s; y*=s; z*=s; return *this; }
};

Vec3 operator*(double s, const Vec3& v) { return v*s; }
double dot(const Vec3& a, const Vec3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};
}
double norm2(const Vec3& v) { return dot(v,v); }
double norm(const Vec3& v) { return sqrt(norm2(v)); }
Vec3 normalized(const Vec3& v) {
    double n = norm(v);
    return (n > EPS) ? v/n : Vec3{0,0,0};
}

struct Particle {
    Vec3 r;      // position
    Vec3 v;      // velocity
    Vec3 f;      // total spring force
    double m;    // mass
};

struct Spring {
    int i, j;
};

struct Box {
    double minCoord;
    double maxCoord;
};

struct Sphere {
    bool enabled;
    Vec3 center;
    double radius;
};

struct CollisionStats {
    long long wallHits = 0;
    long long sphereHits = 0;
    double firstWallTime = -1.0;
    double firstSphereTime = -1.0;
};

struct EnergyData {
    double kinetic = 0.0;
    double spring = 0.0;
    double comKinetic = 0.0;
    double rotational = 0.0;
    double vibrationKinetic = 0.0;
    double vibrationTotal = 0.0;
    double total = 0.0;
};

Vec3 centerOfMass(const Particle p[]) {
    double M = 0.0;
    Vec3 R{0,0,0};
    for (int i=0;i<N;i++) {
        M += p[i].m;
        R += p[i].r * p[i].m;
    }
    return R / M;
}

Vec3 centerOfMassVelocity(const Particle p[]) {
    double M = 0.0;
    Vec3 V{0,0,0};
    for (int i=0;i<N;i++) {
        M += p[i].m;
        V += p[i].v * p[i].m;
    }
    return V / M;
}

void calculateForces(Particle p[], const Spring springs[], double k, double a) {
    for (int i=0;i<N;i++) p[i].f = {0,0,0};

    for (int s=0;s<NS;s++) {
        int i = springs[s].i;
        int j = springs[s].j;

        Vec3 d = p[j].r - p[i].r;
        double L = norm(d);
        if (L < EPS) continue;

        Vec3 force = d * (k * (L-a) / L);
        p[i].f += force;
        p[j].f -= force;
    }
}

double springEnergy(const Particle p[], const Spring springs[], double k, double a) {
    double E = 0.0;
    for (int s=0;s<NS;s++) {
        Vec3 d = p[springs[s].j].r - p[springs[s].i].r;
        double ext = norm(d) - a;
        E += 0.5 * k * ext * ext;
    }
    return E;
}

double kineticEnergy(const Particle p[]) {
    double E = 0.0;
    for (int i=0;i<N;i++) E += 0.5*p[i].m*norm2(p[i].v);
    return E;
}

bool solve3x3(double A[3][3], const Vec3& b, Vec3& x) {
    double M[3][4] = {
        {A[0][0], A[0][1], A[0][2], b.x},
        {A[1][0], A[1][1], A[1][2], b.y},
        {A[2][0], A[2][1], A[2][2], b.z}
    };

    for (int col=0; col<3; col++) {
        int pivot = col;
        for (int row=col+1; row<3; row++)
            if (fabs(M[row][col]) > fabs(M[pivot][col])) pivot = row;

        if (fabs(M[pivot][col]) < 1e-14) return false;
        if (pivot != col)
            for (int c=col;c<4;c++) swap(M[col][c], M[pivot][c]);

        double div = M[col][col];
        for (int c=col;c<4;c++) M[col][c] /= div;

        for (int row=0; row<3; row++) {
            if (row == col) continue;
            double factor = M[row][col];
            for (int c=col;c<4;c++) M[row][c] -= factor*M[col][c];
        }
    }

    x = {M[0][3], M[1][3], M[2][3]};
    return true;
}

EnergyData calculateEnergies(const Particle p[], const Spring springs[], double k, double a) {
    EnergyData E;
    E.kinetic = kineticEnergy(p);
    E.spring = springEnergy(p, springs, k, a);

    double M = 0.0;
    for (int i=0;i<N;i++) M += p[i].m;

    Vec3 R = centerOfMass(p);
    Vec3 V = centerOfMassVelocity(p);
    E.comKinetic = 0.5 * M * norm2(V);

    double I[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    Vec3 L{0,0,0};

    for (int i=0;i<N;i++) {
        Vec3 q = p[i].r - R;
        Vec3 u = p[i].v - V;
        double q2 = norm2(q);
        double m = p[i].m;

        I[0][0] += m*(q2 - q.x*q.x);
        I[1][1] += m*(q2 - q.y*q.y);
        I[2][2] += m*(q2 - q.z*q.z);
        I[0][1] -= m*q.x*q.y; I[1][0] = I[0][1];
        I[0][2] -= m*q.x*q.z; I[2][0] = I[0][2];
        I[1][2] -= m*q.y*q.z; I[2][1] = I[1][2];

        L += cross(q, u) * m;
    }

    Vec3 omega{0,0,0};
    if (solve3x3(I, L, omega)) {
        E.rotational = 0.5 * dot(omega, L);
    } else {
        E.rotational = 0.0;
        omega = {0,0,0};
    }

    E.vibrationKinetic = 0.0;
    for (int i=0;i<N;i++) {
        Vec3 q = p[i].r - R;
        Vec3 u = p[i].v - V;
        Vec3 uVib = u - cross(omega, q);
        E.vibrationKinetic += 0.5*p[i].m*norm2(uVib);
    }

    E.vibrationTotal = E.vibrationKinetic + E.spring;
    E.total = E.kinetic + E.spring;
    return E;
}

// Returns earliest sphere collision time in [0, maxTime].
// Particle is treated as a point and sphere as perfectly rigid/elastic.
bool sphereCollisionTime(const Particle& p, const Sphere& sphere,
                         double maxTime, double& tHit) {
    if (!sphere.enabled) return false;

    Vec3 r = p.r - sphere.center;
    double R = sphere.radius;
    double dist = norm(r);

    // If numerically on the surface and moving inward, collide immediately.
    if (fabs(dist - R) < 1e-10 && dot(r, p.v) < 0.0) {
        tHit = 0.0;
        return true;
    }

    double A = norm2(p.v);
    if (A < EPS) return false;

    double B = 2.0 * dot(r, p.v);
    double C = norm2(r) - R*R;
    double D = B*B - 4.0*A*C;
    if (D < 0.0) return false;

    double sqrtD = sqrt(max(0.0, D));
    double t1 = (-B - sqrtD)/(2.0*A);
    double t2 = (-B + sqrtD)/(2.0*A);

    double best = 1e100;
    if (t1 > 1e-12 && t1 <= maxTime + 1e-12) best = t1;
    if (t2 > 1e-12 && t2 <= maxTime + 1e-12) best = min(best, t2);

    if (best == 1e100) return false;
    tHit = best;
    return true;
}

void driftParticleWithCollisions(Particle& p, double dt,
                                 const Box& box, const Sphere& sphere,
                                 double stepStartTime, CollisionStats& stats) {
    double remaining = dt;
    double elapsed = 0.0;

    // Protect against any tiny numerical penetration of the sphere.
    if (sphere.enabled) {
        Vec3 r = p.r - sphere.center;
        double d = norm(r);
        if (d < sphere.radius) {
            Vec3 n = (d > EPS) ? r/d : Vec3{1,0,0};
            p.r = sphere.center + n*(sphere.radius + 1e-12);
            double vn = dot(p.v, n);
            if (vn < 0.0) p.v -= n*(2.0*vn);
        }
    }

    for (int bounce=0; bounce<12 && remaining > 1e-14; bounce++) {
        double tBest = remaining + 1.0;
        int hitType = 0;  // 0 none, 1 xwall, 2 ywall, 3 zwall, 4 sphere

        // X wall
        if (p.v.x > EPS) {
            double t = (box.maxCoord - p.r.x)/p.v.x;
            if (t >= -1e-14 && t <= remaining + 1e-14 && t < tBest) { tBest=max(0.0,t); hitType=1; }
        } else if (p.v.x < -EPS) {
            double t = (box.minCoord - p.r.x)/p.v.x;
            if (t >= -1e-14 && t <= remaining + 1e-14 && t < tBest) { tBest=max(0.0,t); hitType=1; }
        }

        // Y wall
        if (p.v.y > EPS) {
            double t = (box.maxCoord - p.r.y)/p.v.y;
            if (t >= -1e-14 && t <= remaining + 1e-14 && t < tBest) { tBest=max(0.0,t); hitType=2; }
        } else if (p.v.y < -EPS) {
            double t = (box.minCoord - p.r.y)/p.v.y;
            if (t >= -1e-14 && t <= remaining + 1e-14 && t < tBest) { tBest=max(0.0,t); hitType=2; }
        }

        // Z wall
        if (p.v.z > EPS) {
            double t = (box.maxCoord - p.r.z)/p.v.z;
            if (t >= -1e-14 && t <= remaining + 1e-14 && t < tBest) { tBest=max(0.0,t); hitType=3; }
        } else if (p.v.z < -EPS) {
            double t = (box.minCoord - p.r.z)/p.v.z;
            if (t >= -1e-14 && t <= remaining + 1e-14 && t < tBest) { tBest=max(0.0,t); hitType=3; }
        }

        // Sphere
        double tSphere;
        if (sphereCollisionTime(p, sphere, remaining, tSphere) && tSphere < tBest) {
            tBest = tSphere;
            hitType = 4;
        }

        if (hitType == 0 || tBest > remaining) {
            p.r += p.v * remaining;
            elapsed += remaining;
            remaining = 0.0;
            break;
        }

        // Move exactly to collision.
        p.r += p.v * tBest;
        elapsed += tBest;
        remaining -= tBest;

        double hitTime = stepStartTime + elapsed;

        if (hitType == 1) {
            p.v.x = -p.v.x;
            stats.wallHits++;
            if (stats.firstWallTime < 0.0) stats.firstWallTime = hitTime;
        } else if (hitType == 2) {
            p.v.y = -p.v.y;
            stats.wallHits++;
            if (stats.firstWallTime < 0.0) stats.firstWallTime = hitTime;
        } else if (hitType == 3) {
            p.v.z = -p.v.z;
            stats.wallHits++;
            if (stats.firstWallTime < 0.0) stats.firstWallTime = hitTime;
        } else if (hitType == 4) {
            Vec3 n = normalized(p.r - sphere.center);
            double vn = dot(p.v, n);
            p.v -= n * (2.0*vn);
            p.r += n * 1e-12; // tiny outward nudge
            stats.sphereHits++;
            if (stats.firstSphereTime < 0.0) stats.firstSphereTime = hitTime;
        }

        // Tiny forward move is not needed; the event logic and nudge avoid re-hits.
        if (remaining < 1e-14) break;
    }

    // Safety clamp against roundoff at box boundaries.
    p.r.x = min(box.maxCoord, max(box.minCoord, p.r.x));
    p.r.y = min(box.maxCoord, max(box.minCoord, p.r.y));
    p.r.z = min(box.maxCoord, max(box.minCoord, p.r.z));
}

void initializeTetrahedron(Particle p[], double a, double m,
                           const Vec3& comPosition, const Vec3& velocity) {
    double s = a/(2.0*sqrt(2.0));
    Vec3 q[N] = {
        { s,  s,  s},
        { s, -s, -s},
        {-s,  s, -s},
        {-s, -s,  s}
    };

    for (int i=0;i<N;i++) {
        p[i].r = comPosition + q[i];
        p[i].v = velocity;
        p[i].f = {0,0,0};
        p[i].m = m;
    }
}

void writeHeader(ofstream& out) {
    out << "time";
    for (int i=0;i<N;i++) out << ",p"<<i<<"x,p"<<i<<"y,p"<<i<<"z";
    out << ",com_x,com_y,com_z"
        << ",kinetic,spring,com_kinetic,rotational,vibration_kinetic,vibration_total,total,relative_energy_error\n";
}

void writeRow(ofstream& out, double t, const Particle p[],
              const EnergyData& E, double E0) {
    Vec3 R = centerOfMass(p);
    out << setprecision(15) << t;
    for (int i=0;i<N;i++) out << ',' << p[i].r.x << ',' << p[i].r.y << ',' << p[i].r.z;
    double relErr = (fabs(E0) > EPS) ? fabs(E.total-E0)/fabs(E0) : 0.0;
    out << ',' << R.x << ',' << R.y << ',' << R.z
        << ',' << E.kinetic
        << ',' << E.spring
        << ',' << E.comKinetic
        << ',' << E.rotational
        << ',' << E.vibrationKinetic
        << ',' << E.vibrationTotal
        << ',' << E.total
        << ',' << relErr << '\n';
}

struct RunSummary {
    double initialEnergy = 0.0;
    double finalEnergy = 0.0;
    double maxRelativeError = 0.0;
    double minComEnergy = 1e100;
    double maxVibrationEnergy = 0.0;
    double maxRotationalEnergy = 0.0;
    double timeOfMaxVibration = 0.0;
    CollisionStats collisions;
};

RunSummary runSimulation(const string& name,
                         const Vec3& initialCOM,
                         const Vec3& initialVelocity,
                         bool sphereEnabled,
                         double dt,
                         double totalTime,
                         int outputEvery) {
    const double a = 1.0;
    const double k = 10.0;
    const double m = 1.0;
    const Box box{-3.0, 3.0};
    const Sphere sphere{sphereEnabled, {0,0,0}, 0.75};
    const Spring springs[NS] = {
        {0,1},{0,2},{0,3},{1,2},{1,3},{2,3}
    };

    Particle p[N];
    initializeTetrahedron(p, a, m, initialCOM, initialVelocity);
    calculateForces(p, springs, k, a);

    EnergyData E = calculateEnergies(p, springs, k, a);
    double E0 = E.total;

    RunSummary summary;
    summary.initialEnergy = E0;
    summary.minComEnergy = E.comKinetic;
    summary.maxVibrationEnergy = E.vibrationTotal;
    summary.maxRotationalEnergy = E.rotational;

    ofstream out(name + ".csv");
    if (!out) {
        cerr << "Could not open output file: " << name << ".csv\n";
        return summary;
    }
    writeHeader(out);
    writeRow(out, 0.0, p, E, E0);

    int steps = static_cast<int>(llround(totalTime/dt));

    for (int step=0; step<steps; step++) {
        double t = step*dt;

        // Velocity Verlet: half kick
        for (int i=0;i<N;i++) p[i].v += p[i].f * (0.5*dt/p[i].m);

        // Drift for dt, resolving exact hard-wall / hard-sphere events
        for (int i=0;i<N;i++)
            driftParticleWithCollisions(p[i], dt, box, sphere, t, summary.collisions);

        // New spring forces
        calculateForces(p, springs, k, a);

        // Velocity Verlet: second half kick
        for (int i=0;i<N;i++) p[i].v += p[i].f * (0.5*dt/p[i].m);

        double timeNow = (step+1)*dt;
        E = calculateEnergies(p, springs, k, a);
        double relErr = fabs(E.total-E0)/fabs(E0);
        summary.maxRelativeError = max(summary.maxRelativeError, relErr);
        summary.minComEnergy = min(summary.minComEnergy, E.comKinetic);
        if (E.vibrationTotal > summary.maxVibrationEnergy) {
            summary.maxVibrationEnergy = E.vibrationTotal;
            summary.timeOfMaxVibration = timeNow;
        }
        summary.maxRotationalEnergy = max(summary.maxRotationalEnergy, E.rotational);

        if ((step+1)%outputEvery == 0 || step == steps-1)
            writeRow(out, timeNow, p, E, E0);
    }

    summary.finalEnergy = E.total;
    return summary;
}

void printSummary(const string& title, const RunSummary& s) {
    cout << "\n=== " << title << " ===\n";
    cout << setprecision(10);
    cout << "Initial total energy       : " << s.initialEnergy << '\n';
    cout << "Final total energy         : " << s.finalEnergy << '\n';
    cout << "Maximum relative error     : " << s.maxRelativeError << '\n';
    cout << "Minimum COM kinetic energy : " << s.minComEnergy << '\n';
    cout << "Maximum vibration energy   : " << s.maxVibrationEnergy << '\n';
    cout << "Time of max vibration      : " << s.timeOfMaxVibration << '\n';
    cout << "Maximum rotational energy  : " << s.maxRotationalEnergy << '\n';
    cout << "Wall hits                  : " << s.collisions.wallHits << '\n';
    cout << "Sphere hits                : " << s.collisions.sphereHits << '\n';
    if (s.collisions.firstWallTime >= 0.0)
        cout << "First wall collision time  : " << s.collisions.firstWallTime << '\n';
    if (s.collisions.firstSphereTime >= 0.0)
        cout << "First sphere collision time: " << s.collisions.firstSphereTime << '\n';
}

int main() {
    const double dt = 0.0001;
    const int outputEvery = 50;

    // Experiment 1: tetrahedron bouncing from box walls.
    RunSummary wall = runSimulation(
        "wall_experiment",
        {-1.5, 0.0, 0.0},
        { 1.2, 0.0, 0.0},
        false,
        dt,
        8.0,
        outputEvery
    );

    // Experiment 2: tetrahedron hits an off-centre rigid sphere at box centre.
    RunSummary sphere = runSimulation(
        "sphere_experiment",
        {-2.0, 0.45, 0.25},
        { 1.0, 0.0, 0.0},
        true,
        dt,
        4.0,
        outputEvery
    );

    printSummary("WALL EXPERIMENT", wall);
    printSummary("SPHERE EXPERIMENT", sphere);

    cout << "\nCSV files created:\n";
    cout << "  wall_experiment.csv\n";
    cout << "  sphere_experiment.csv\n";
    cout << "\nUpload these CSV files here and we can check the results and make the plots/report.\n";

    return 0;
}