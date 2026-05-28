#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "Vec3.hpp"
#include "MieTable.hpp"

#include <cmath>
#include <vector>
#include <iostream>
#include <cstdlib>

constexpr double PI = 3.14159265358979323846;
constexpr double EPS = 1e-6;

struct Ray {
    Vec3 origin;
    Vec3 dir;
};

struct Segment {
    Vec3 a;
    Vec3 b;
    float r;
    float g;
    float bcol;
};

std::vector<Segment> segments;
MieTable mieTable;

double sphereRadius = 1.0;
double dropletRadiusUm = 100.0;

int internalReflections = 1;

double cameraDistance = 6.0;
double cameraYaw = 25.0;
double cameraPitch = 20.0;

int lastMouseX = 0;
int lastMouseY = 0;
bool dragging = false;

Vec3 reflectVec(const Vec3& I, const Vec3& N) {
    return normalize(I - N * (2.0 * dot(I, N)));
}

bool refractVec(const Vec3& I, const Vec3& N, double eta, Vec3& T) {
    double cosi = -dot(N, I);
    double sint2 = eta * eta * (1.0 - cosi * cosi);

    if (sint2 > 1.0) {
        return false;
    }

    double cost = std::sqrt(1.0 - sint2);
    T = normalize(I * eta + N * (eta * cosi - cost));
    return true;
}

bool intersectSphere(const Ray& ray, double radius, double& tHit) {
    Vec3 oc = ray.origin;

    double A = dot(ray.dir, ray.dir);
    double B = 2.0 * dot(oc, ray.dir);
    double C = dot(oc, oc) - radius * radius;

    double D = B * B - 4.0 * A * C;

    if (D < 0.0) {
        return false;
    }

    double sqrtD = std::sqrt(D);

    double t1 = (-B - sqrtD) / (2.0 * A);
    double t2 = (-B + sqrtD) / (2.0 * A);

    if (t1 > EPS) {
        tHit = t1;
        return true;
    }

    if (t2 > EPS) {
        tHit = t2;
        return true;
    }

    return false;
}

double waterRefractiveIndex(double wavelength_nm) {
    double lambda_um = wavelength_nm / 1000.0;
    return 1.322 + 0.003 / (lambda_um * lambda_um);
}

void wavelengthToRGB(double wavelength_nm, float& r, float& g, float& b) {
    r = 0.0f;
    g = 0.0f;
    b = 0.0f;

    if (wavelength_nm >= 380.0 && wavelength_nm < 440.0) {
        r = static_cast<float>(-(wavelength_nm - 440.0) / (440.0 - 380.0));
        g = 0.0f;
        b = 1.0f;
    } else if (wavelength_nm < 490.0) {
        r = 0.0f;
        g = static_cast<float>((wavelength_nm - 440.0) / (490.0 - 440.0));
        b = 1.0f;
    } else if (wavelength_nm < 510.0) {
        r = 0.0f;
        g = 1.0f;
        b = static_cast<float>(-(wavelength_nm - 510.0) / (510.0 - 490.0));
    } else if (wavelength_nm < 580.0) {
        r = static_cast<float>((wavelength_nm - 510.0) / (580.0 - 510.0));
        g = 1.0f;
        b = 0.0f;
    } else if (wavelength_nm < 645.0) {
        r = 1.0f;
        g = static_cast<float>(-(wavelength_nm - 645.0) / (645.0 - 580.0));
        b = 0.0f;
    } else if (wavelength_nm <= 750.0) {
        r = 1.0f;
        g = 0.0f;
        b = 0.0f;
    }

    double factor = 1.0;

    if (wavelength_nm >= 380.0 && wavelength_nm < 420.0) {
        factor = 0.3 + 0.7 * (wavelength_nm - 380.0) / (420.0 - 380.0);
    } else if (wavelength_nm > 700.0 && wavelength_nm <= 750.0) {
        factor = 0.3 + 0.7 * (750.0 - wavelength_nm) / (750.0 - 700.0);
    }

    r = static_cast<float>(r * factor);
    g = static_cast<float>(g * factor);
    b = static_cast<float>(b * factor);
}

double scatteringAngleDeg(const Vec3& incoming, const Vec3& outgoing) {
    Vec3 a = normalize(incoming);
    Vec3 b = normalize(outgoing);

    double c = clampDouble(dot(a, b), -1.0, 1.0);
    return std::acos(c) * 180.0 / PI;
}

void addSegment(const Vec3& a, const Vec3& b, float r, float g, float colb) {
    segments.push_back({a, b, r, g, colb});
}

bool traceDropRay(
    double impactParameter,
    double wavelength_nm,
    int reflections,
    float cr,
    float cg,
    float cb
) {
    double nAir = 1.0;
    double nWater = waterRefractiveIndex(wavelength_nm);

    Vec3 incomingDir = {1.0, 0.0, 0.0};

    Ray ray;
    ray.origin = {-3.0, impactParameter, 0.0};
    ray.dir = incomingDir;

    double t = 0.0;

    if (!intersectSphere(ray, sphereRadius, t)) {
        return false;
    }

    Vec3 p1 = ray.origin + ray.dir * t;
    Vec3 n1 = normalize(p1);

    addSegment(ray.origin, p1, cr, cg, cb);

    Vec3 insideDir;

    if (!refractVec(ray.dir, n1, nAir / nWater, insideDir)) {
        return false;
    }

    Vec3 currentPoint = p1 + insideDir * EPS;
    Vec3 currentDir = insideDir;

    for (int k = 0; k < reflections; ++k) {
        Ray innerRay{currentPoint, currentDir};

        double tInner = 0.0;

        if (!intersectSphere(innerRay, sphereRadius, tInner)) {
            return false;
        }

        Vec3 p = innerRay.origin + innerRay.dir * tInner;
        Vec3 n = normalize(p);

        addSegment(currentPoint, p, cr, cg, cb);

        currentDir = reflectVec(currentDir, n);
        currentPoint = p + currentDir * EPS;
    }

    Ray exitRay{currentPoint, currentDir};

    double tExit = 0.0;

    if (!intersectSphere(exitRay, sphereRadius, tExit)) {
        return false;
    }

    Vec3 pExit = exitRay.origin + exitRay.dir * tExit;
    Vec3 nExit = normalize(pExit);

    addSegment(currentPoint, pExit, cr, cg, cb);

    Vec3 outsideDir;

    if (!refractVec(currentDir, nExit * (-1.0), nWater / nAir, outsideDir)) {
        outsideDir = reflectVec(currentDir, nExit);
    }

    double theta = scatteringAngleDeg(incomingDir, outsideDir);

    double mieWeight = mieTable.lookupNearest(
    dropletRadiusUm,
    wavelength_nm,
    theta
    );

    mieWeight =
    std::log(1.0 + 20.0*mieWeight) /
    std::log(21.0);


    float wr = static_cast<float>(cr * mieWeight);
    float wg = static_cast<float>(cg * mieWeight);
    float wb = static_cast<float>(cb * mieWeight);

    Vec3 pFar = pExit + outsideDir * 3.0;
    addSegment(pExit, pFar, wr, wg, wb);

    return true;
}

void rebuildRays() {
    segments.clear();

    std::vector<double> wavelengths = {
        430.0, 470.0, 510.0, 550.0, 590.0, 630.0, 680.0
    };

    for (double lambda : wavelengths) {
        float r = 0.0f;
        float g = 0.0f;
        float b = 0.0f;

        wavelengthToRGB(lambda, r, g, b);

        for (double p = -0.85; p <= 0.85; p += 0.17) {
            traceDropRay(p, lambda, internalReflections, r, g, b);
        }
    }
}

void drawSphere() {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glColor4f(0.55f, 0.80f, 1.0f, 0.22f);
    glutSolidSphere(sphereRadius, 64, 64);

    glDisable(GL_BLEND);

    glColor3f(0.75f, 0.95f, 1.0f);
    glutWireSphere(sphereRadius, 32, 32);
}

void drawAxes() {
    glLineWidth(1.5f);
    glBegin(GL_LINES);

    glColor3f(1.0f, 0.1f, 0.1f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(2.0f, 0.0f, 0.0f);

    glColor3f(0.1f, 1.0f, 0.1f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 2.0f, 0.0f);

    glColor3f(0.1f, 0.3f, 1.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 2.0f);

    glEnd();
}

void drawRays() {
    glLineWidth(2.0f);
    glBegin(GL_LINES);

    for (const auto& s : segments) {
        glColor3f(s.r, s.g, s.bcol);
        glVertex3d(s.a.x, s.a.y, s.a.z);
        glVertex3d(s.b.x, s.b.y, s.b.z);
    }

    glEnd();
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glTranslated(0.0, 0.0, -cameraDistance);
    glRotated(cameraPitch, 1.0, 0.0, 0.0);
    glRotated(cameraYaw, 0.0, 1.0, 0.0);

    drawAxes();
    drawSphere();
    drawRays();

    glutSwapBuffers();
}

void reshape(int w, int h) {
    if (h == 0) {
        h = 1;
    }

    glViewport(0, 0, w, h);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    double aspect = static_cast<double>(w) / static_cast<double>(h);
    gluPerspective(45.0, aspect, 0.1, 100.0);

    glMatrixMode(GL_MODELVIEW);
}

void keyboard(unsigned char key, int, int) {
    if (key == 27) {
        std::exit(0);
    }

    if (key == 'w' || key == 'W') {
        cameraDistance -= 0.3;
        if (cameraDistance < 2.0) {
            cameraDistance = 2.0;
        }
    }

    if (key == 's' || key == 'S') {
        cameraDistance += 0.3;
    }

    if (key == '1') {
        internalReflections = 1;
        std::cout << "Mode: primary rainbow, 1 internal reflection\n";
        rebuildRays();
    }

    if (key == '2') {
        internalReflections = 2;
        std::cout << "Mode: secondary rainbow, 2 internal reflections\n";
        rebuildRays();
    }

    if (key == 'r' || key == 'R') {
        cameraDistance = 6.0;
        cameraYaw = 25.0;
        cameraPitch = 20.0;
    }

    glutPostRedisplay();
}

void mouse(int button, int state, int x, int y) {
    if (button == GLUT_LEFT_BUTTON) {
        dragging = (state == GLUT_DOWN);
        lastMouseX = x;
        lastMouseY = y;
    }
}

void motion(int x, int y) {
    if (!dragging) {
        return;
    }

    int dx = x - lastMouseX;
    int dy = y - lastMouseY;

    cameraYaw += dx * 0.5;
    cameraPitch += dy * 0.5;

    lastMouseX = x;
    lastMouseY = y;

    glutPostRedisplay();
}

void initOpenGL() {
    glClearColor(0.02f, 0.02f, 0.035f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LINE_SMOOTH);
    glShadeModel(GL_SMOOTH);
}

int main(int argc, char** argv) {
    bool loaded = mieTable.loadCSV("data/mie_table.csv");

    if (!loaded) {
        std::cout << "Mie table not found. Running with fallback intensity = 1.\n";
        std::cout << "Run first: python tools/generate_mie_table.py\n";
    }

    rebuildRays();

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(1000, 750);
    glutCreateWindow("Single Drop Rainbow with Mie Lookup Table");

    initOpenGL();

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);

    std::cout << "Controls:\n";
    std::cout << "Mouse drag - rotate camera\n";
    std::cout << "W/S        - zoom\n";
    std::cout << "1          - primary rainbow\n";
    std::cout << "2          - secondary rainbow\n";
    std::cout << "R          - reset camera\n";
    std::cout << "Esc        - exit\n";

    glutMainLoop();

    return 0;
}



