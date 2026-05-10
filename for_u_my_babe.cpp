/*

everything for u my babe

установка пакетов Arch Linux: 
    sudo pacman -S --needed base-devel sdl2 mesa glu

компи:
    g++ rainbow_arch_wayland.cpp -std=c++17 $(sdl2-config --cflags --libs) -lGL -lGLU -O2 -o RainbowSimulation

запуск ч/з wayland
    SDL_VIDEODRIVER=wayland ./RainbowSimulation

При запуске в терминале должно быть:
    Current SDL video driver: wayland


Controls:
        W/S/A/D      move observer forward/back/left/right
        Q/E          move observer down/up
        Mouse drag   rotate camera
        Arrow keys   move light direction
        + / -        change droplet radius
        [ / ]        change wavelength in monochrome mode
        1            geometric optics mode
        2            wave approximation mode
        3            Mie-like approximation mode
        B            white light / single wavelength
        V            cycle visualization mode
        R            reset scene
        L            show/hide ray paths
        N            show/hide normals
        P            pause animation
        8            primary + secondary rainbow
        9            primary rainbow only
        0            no internal reflections
        Esc          exit

Notes:
    This version does not use GLUT. It uses SDL2 for native Wayland window/input.
    Rendering uses old OpenGL immediate mode for simplicity and portability in a student project.
    Exact Mie theory is replaced by a qualitative Mie-like approximation.


*/







#include <SDL2/SDL.h>
#include <GL/gl.h>
#include <GL/glu.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>

const double PI = 3.14159265358979323846;

static double clampd(double x, double a, double b) {
    return x < a ? a : (x > b ? b : x);
}

static double rad(double d) { return d * PI / 180.0; }

struct Vec3 {
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double X, double Y, double Z) : x(X), y(Y), z(Z) {}

    Vec3 operator+(const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3 operator-(const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    Vec3 operator*(double k) const { return Vec3(x * k, y * k, z * k); }
    Vec3 operator/(double k) const { return Vec3(x / k, y / k, z / k); }
    Vec3& operator+=(const Vec3& v) { x += v.x; y += v.y; z += v.z; return *this; }
    Vec3& operator-=(const Vec3& v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
};

static Vec3 operator*(double k, const Vec3& v) { return v * k; }

static double dot(const Vec3& a, const Vec3& b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

static Vec3 cross(const Vec3& a, const Vec3& b) {
    return Vec3(a.y*b.z - a.z*b.y,
                a.z*b.x - a.x*b.z,
                a.x*b.y - a.y*b.x);
}

static double length(const Vec3& v) {
    return std::sqrt(dot(v, v));
}

static Vec3 normalize(const Vec3& v) {
    double l = length(v);
    if (l < 1e-12) return Vec3(0, 0, 0);
    return v / l;
}

static Vec3 reflect(const Vec3& I, const Vec3& N) {
    return normalize(I - N * (2.0 * dot(I, N)));
}

static bool refract(const Vec3& I, const Vec3& N, double n1, double n2, Vec3& T) {
    Vec3 i = normalize(I);
    Vec3 n = normalize(N);
    double cosi = clampd(dot(i, n), -1.0, 1.0);
    double eta = n1 / n2;

    if (cosi > 0.0) {
        n = n * -1.0;
        cosi = -cosi;
    }

    double k = 1.0 - eta * eta * (1.0 - cosi * cosi);
    if (k < 0.0) return false;

    T = normalize(i * eta + n * (eta * (-cosi) - std::sqrt(k)));
    return true;
}

struct Color {
    double r, g, b;
    Color() : r(1), g(1), b(1) {}
    Color(double R, double G, double B) : r(R), g(G), b(B) {}
};

static void glColor(const Color& c, double alpha = 1.0) {
    glColor4d(clampd(c.r, 0, 1), clampd(c.g, 0, 1), clampd(c.b, 0, 1), clampd(alpha, 0, 1));
}

static Color wavelengthToRGB(double wavelength) {
    double R = 0.0, G = 0.0, B = 0.0;
    double wl = wavelength;

    if (wl >= 380 && wl < 440) { R = -(wl - 440) / (440 - 380); G = 0.0; B = 1.0; }
    else if (wl >= 440 && wl < 490) { R = 0.0; G = (wl - 440) / (490 - 440); B = 1.0; }
    else if (wl >= 490 && wl < 510) { R = 0.0; G = 1.0; B = -(wl - 510) / (510 - 490); }
    else if (wl >= 510 && wl < 580) { R = (wl - 510) / (580 - 510); G = 1.0; B = 0.0; }
    else if (wl >= 580 && wl < 645) { R = 1.0; G = -(wl - 645) / (645 - 580); B = 0.0; }
    else if (wl >= 645 && wl <= 780) { R = 1.0; G = 0.0; B = 0.0; }

    double factor;
    if (wl >= 380 && wl < 420) factor = 0.3 + 0.7 * (wl - 380) / (420 - 380);
    else if (wl >= 420 && wl <= 700) factor = 1.0;
    else if (wl > 700 && wl <= 780) factor = 0.3 + 0.7 * (780 - wl) / (780 - 700);
    else factor = 0.0;

    const double gamma = 0.8;
    return Color(std::pow(R * factor, gamma),
                 std::pow(G * factor, gamma),
                 std::pow(B * factor, gamma));
}

static double waterRefractiveIndex(double wavelengthNm) {
    double t = clampd((wavelengthNm - 380.0) / (780.0 - 380.0), 0.0, 1.0);
    return 1.3435 - 0.0125 * t;
}

struct RayPath {
    std::vector<Vec3> points;
    Color color;
    bool valid;
    RayPath() : valid(false) {}
};

struct Settings {
    int width = 1280;
    int height = 800;
    int model = 0;
    int viewMode = 0;
    bool whiteLight = true;
    bool showRays = true;
    bool showNormals = true;
    bool paused = false;
    double wavelength = 550.0;
    double dropletRadius = 1.0;
    double rainDensity = 1.0;
    double sizeSpread = 0.10;
    int internalReflections = 2;
    double fov = 60.0;
} S;

struct CameraState {
    Vec3 pos = Vec3(0.0, 1.2, 8.0);
    double yaw = 180.0;
    double pitch = -8.0;
} Cam;

Vec3 lightDir = normalize(Vec3(-0.9, -0.25, -0.35));
bool keyState[SDL_NUM_SCANCODES] = {false};
bool mouseDown = false;
int lastMouseX = 0;
int lastMouseY = 0;
double animTime = 0.0;

static Vec3 cameraForward() {
    double cy = std::cos(rad(Cam.yaw));
    double sy = std::sin(rad(Cam.yaw));
    double cp = std::cos(rad(Cam.pitch));
    double sp = std::sin(rad(Cam.pitch));
    return normalize(Vec3(sy * cp, sp, cy * cp));
}

static Vec3 cameraRight() {
    return normalize(cross(cameraForward(), Vec3(0, 1, 0)));
}

static Vec3 cameraUp() {
    return normalize(cross(cameraRight(), cameraForward()));
}

static bool raySphere(const Vec3& O, const Vec3& D, const Vec3& C, double R, double& t0, double& t1) {
    Vec3 L = O - C;
    double a = dot(D, D);
    double b = 2.0 * dot(D, L);
    double c = dot(L, L) - R*R;
    double disc = b*b - 4*a*c;
    if (disc < 0.0) return false;
    double s = std::sqrt(disc);
    t0 = (-b - s) / (2*a);
    t1 = (-b + s) / (2*a);
    return true;
}

static RayPath traceDropRay(double impact, double wavelength, int reflections) {
    RayPath path;
    path.color = wavelengthToRGB(wavelength);

    Vec3 C(0, 0, 0);
    double R = S.dropletRadius;
    Vec3 D = normalize(Vec3(0, 0, 1));
    Vec3 O(impact * R, 0, -4.0 * R);
    double t0, t1;
    if (!raySphere(O, D, C, R, t0, t1)) return path;

    Vec3 P = O + D * t0;
    Vec3 N = normalize(P - C);
    Vec3 inside;
    if (!refract(D, N, 1.0, waterRefractiveIndex(wavelength), inside)) return path;

    path.points.push_back(O);
    path.points.push_back(P);

    Vec3 currentP = P + inside * 1e-5;
    Vec3 currentD = inside;

    for (int k = 0; k < reflections; ++k) {
        double a0, a1;
        if (!raySphere(currentP, currentD, C, R, a0, a1)) return path;
        Vec3 hit = currentP + currentD * a1;
        path.points.push_back(hit);
        Vec3 n = normalize(hit - C);
        currentD = reflect(currentD, n);
        currentP = hit + currentD * 1e-5;
    }

    double b0, b1;
    if (!raySphere(currentP, currentD, C, R, b0, b1)) return path;
    Vec3 exitP = currentP + currentD * b1;
    path.points.push_back(exitP);

    Vec3 exitN = normalize(exitP - C);
    Vec3 out;
    if (!refract(currentD, exitN * -1.0, waterRefractiveIndex(wavelength), 1.0, out)) {
        out = reflect(currentD, exitN);
    }
    path.points.push_back(exitP + out * 3.0 * R);
    path.valid = true;
    return path;
}

static double rainbowAnglePrimary(double wavelength) {
    double t = clampd((wavelength - 380.0) / 400.0, 0.0, 1.0);
    return 40.2 + 2.2 * t;
}

static double rainbowAngleSecondary(double wavelength) {
    double t = clampd((wavelength - 380.0) / 400.0, 0.0, 1.0);
    return 54.5 - 3.6 * t;
}

static double gaussian(double x, double sigma) {
    return std::exp(-(x*x)/(2.0*sigma*sigma));
}

static double modelIntensity(double angleDeg, double centerDeg, double wavelength, int order) {
    double width = order == 1 ? 0.28 : 0.45;
    double base = gaussian(angleDeg - centerDeg, width);

    if (S.model == 0) return base;

    if (S.model == 1) {
        double radiusFactor = clampd(1.6 / S.dropletRadius, 0.2, 3.5);
        double phase = (angleDeg - centerDeg) * (10.0 + 4.0 * radiusFactor);
        double fringes = 0.62 + 0.38 * std::cos(phase);
        double supernumerary = gaussian(angleDeg - (centerDeg - 1.1), 2.2) * fringes;
        return 0.75 * base + 0.45 * supernumerary * (1.0 - clampd(S.sizeSpread * 3.0, 0.0, 0.9));
    }

    double broadening = clampd(0.5 / S.dropletRadius, 0.0, 1.3);
    double soft = gaussian(angleDeg - centerDeg, width + broadening);
    double ripple = 0.85 + 0.15 * std::cos((angleDeg - centerDeg) * 18.0 + wavelength * 0.03);
    return soft * ripple;
}

static void drawLine(const Vec3& a, const Vec3& b, const Color& c, double width = 2.0, double alpha = 1.0) {
    glLineWidth((GLfloat)width);
    glColor(c, alpha);
    glBegin(GL_LINES);
    glVertex3d(a.x, a.y, a.z);
    glVertex3d(b.x, b.y, b.z);
    glEnd();
    glLineWidth(1.0f);
}

static void drawGrid(double size, int lines) {
    glColor4d(0.25, 0.25, 0.28, 0.65);
    glBegin(GL_LINES);
    for (int i = -lines; i <= lines; ++i) {
        double p = size * i / lines;
        glVertex3d(p, 0, -size); glVertex3d(p, 0, size);
        glVertex3d(-size, 0, p); glVertex3d(size, 0, p);
    }
    glEnd();
}

static void drawWireCube(const Vec3& center, const Vec3& scale, const Color& c, double alpha) {
    glPushMatrix();
    glTranslated(center.x, center.y, center.z);
    glScaled(scale.x, scale.y, scale.z);
    glColor(c, alpha);
    glBegin(GL_LINES);
    double v[8][3] = {
        {-0.5,-0.5,-0.5},{0.5,-0.5,-0.5},{0.5,0.5,-0.5},{-0.5,0.5,-0.5},
        {-0.5,-0.5,0.5},{0.5,-0.5,0.5},{0.5,0.5,0.5},{-0.5,0.5,0.5}
    };
    int e[12][2] = {{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7}};
    for (int i = 0; i < 12; ++i) {
        glVertex3dv(v[e[i][0]]);
        glVertex3dv(v[e[i][1]]);
    }
    glEnd();
    glPopMatrix();
}

static void drawSphereMesh(const Vec3& center, double radius, const Color& c, double alpha, bool wire) {
    glPushMatrix();
    glTranslated(center.x, center.y, center.z);
    glColor(c, alpha);

    int slices = 36;
    int stacks = 18;
    if (wire) {
        glLineWidth(1.0f);
        for (int i = 1; i < stacks; ++i) {
            double th = -PI/2.0 + PI * i / stacks;
            double y = radius * std::sin(th);
            double rr = radius * std::cos(th);
            glBegin(GL_LINE_LOOP);
            for (int j = 0; j < slices; ++j) {
                double ph = 2.0 * PI * j / slices;
                glVertex3d(rr * std::cos(ph), y, rr * std::sin(ph));
            }
            glEnd();
        }
        for (int j = 0; j < slices; j += 3) {
            double ph = 2.0 * PI * j / slices;
            glBegin(GL_LINE_STRIP);
            for (int i = 0; i <= stacks; ++i) {
                double th = -PI/2.0 + PI * i / stacks;
                glVertex3d(radius * std::cos(th) * std::cos(ph),
                           radius * std::sin(th),
                           radius * std::cos(th) * std::sin(ph));
            }
            glEnd();
        }
    } else {
        for (int i = 0; i < stacks; ++i) {
            double th1 = -PI/2.0 + PI * i / stacks;
            double th2 = -PI/2.0 + PI * (i + 1) / stacks;
            glBegin(GL_QUAD_STRIP);
            for (int j = 0; j <= slices; ++j) {
                double ph = 2.0 * PI * j / slices;
                glVertex3d(radius * std::cos(th1) * std::cos(ph),
                           radius * std::sin(th1),
                           radius * std::cos(th1) * std::sin(ph));
                glVertex3d(radius * std::cos(th2) * std::cos(ph),
                           radius * std::sin(th2),
                           radius * std::cos(th2) * std::sin(ph));
            }
            glEnd();
        }
    }
    glPopMatrix();
}

static void drawRainVolume() {
    drawWireCube(Vec3(0, 2.4, -7.0), Vec3(8.0, 4.0, 2.5), Color(0.35, 0.55, 1.0), 0.20);

    glPointSize(3.0f);
    glBegin(GL_POINTS);
    for (int i = 0; i < 360; ++i) {
        double a = i * 12.9898;
        double b = i * 78.233;
        double x = std::fmod(std::sin(a) * 43758.5453, 1.0);
        double y = std::fmod(std::sin(b) * 24634.6345, 1.0);
        double z = std::fmod(std::sin(a+b) * 12345.1234, 1.0);
        if (x < 0) x += 1; if (y < 0) y += 1; if (z < 0) z += 1;
        glColor4d(0.55, 0.72, 1.0, 0.25);
        glVertex3d((x - 0.5) * 8.0, 0.5 + y * 4.0, -7.0 + (z - 0.5) * 2.5);
    }
    glEnd();
}

static void drawSingleDropMode() {
    Vec3 C(0, 1.3, -4.0);
    drawSphereMesh(C, S.dropletRadius, Color(0.45, 0.75, 1.0), 0.22, false);
    drawSphereMesh(C, S.dropletRadius, Color(0.8, 0.95, 1.0), 0.90, true);

    std::vector<double> wavelengths;
    if (S.whiteLight) {
        double vals[] = {650, 590, 530, 480, 420};
        wavelengths.assign(vals, vals + 5);
    } else {
        wavelengths.push_back(S.wavelength);
    }

    for (size_t wi = 0; wi < wavelengths.size(); ++wi) {
        double impact = -0.65 + wi * 0.065;
        RayPath rp = traceDropRay(impact, wavelengths[wi], S.internalReflections);
        if (!rp.valid) continue;
        for (size_t i = 1; i < rp.points.size(); ++i) {
            drawLine(rp.points[i-1] + C, rp.points[i] + C, rp.color, 3.0, 0.95);
        }
    }

    if (S.showNormals) {
        Vec3 p = C + Vec3(-0.65*S.dropletRadius, 0, -std::sqrt(1.0-0.65*0.65)*S.dropletRadius);
        Vec3 n = normalize(p - C);
        drawLine(p, p + n * 0.8, Color(1,1,1), 2.0, 0.9);
    }
}

static void basisAround(const Vec3& axis, Vec3& u, Vec3& v) {
    Vec3 a = std::fabs(axis.y) < 0.9 ? Vec3(0,1,0) : Vec3(1,0,0);
    u = normalize(cross(axis, a));
    v = normalize(cross(u, axis));
}

static void drawRainbowRing(const Vec3& origin, const Vec3& antiSolar, const Vec3& u, const Vec3& v,
                            double angleDeg, double dist, const Color& c, double alpha, double width) {
    double theta = rad(angleDeg);
    glLineWidth((GLfloat)width);
    glColor(c, alpha);
    glBegin(GL_LINE_LOOP);
    for (int i = 0; i < 360; ++i) {
        double phi = rad(i);
        Vec3 dir = normalize(antiSolar * std::cos(theta) +
                             (u * std::cos(phi) + v * std::sin(phi)) * std::sin(theta));
        Vec3 p = origin + dir * dist;
        glVertex3d(p.x, p.y, p.z);
    }
    glEnd();
    glLineWidth(1.0f);
}

static void drawRainbowCones() {
    Vec3 antiSolar = normalize(lightDir * -1.0);
    Vec3 origin = Cam.pos;
    Vec3 u, v;
    basisAround(antiSolar, u, v);
    double dist = 18.0;

    std::vector<double> wavelengths;
    if (S.whiteLight) {
        for (int i = 0; i < 32; ++i) wavelengths.push_back(380.0 + i * (400.0 / 31.0));
    } else {
        wavelengths.push_back(S.wavelength);
    }

    for (double wl : wavelengths) {
        Color c = wavelengthToRGB(wl);
        double primary = rainbowAnglePrimary(wl);
        double secondary = rainbowAngleSecondary(wl);

        if (S.model == 0) {
            drawRainbowRing(origin, antiSolar, u, v, primary, dist, c, 0.95, 3.0);
            if (S.internalReflections >= 2) drawRainbowRing(origin, antiSolar, u, v, secondary, dist, c, 0.35, 2.0);
        }
        else if (S.model == 1) {
            drawRainbowRing(origin, antiSolar, u, v, primary, dist, c, 0.90, 3.0);
            for (int fringe = 1; fringe <= 5; ++fringe) {
                drawRainbowRing(origin, antiSolar, u, v, primary - 0.55 * fringe, dist, c, 0.45 / fringe, 2.0);
            }
            if (S.internalReflections >= 2) drawRainbowRing(origin, antiSolar, u, v, secondary, dist, c, 0.25, 2.0);
        }
        else {
            for (int layer = -4; layer <= 4; ++layer) {
                double offset = layer * 0.22;
                double alpha = 0.38 * gaussian(offset, 0.55);
                double ripple = 0.75 + 0.25 * std::cos(wl * 0.04 + layer * 1.7);
                drawRainbowRing(origin, antiSolar, u, v, primary + offset, dist, c, alpha * ripple, 3.0);
            }
            if (S.internalReflections >= 2) {
                for (int layer = -3; layer <= 3; ++layer) {
                    double offset = layer * 0.28;
                    drawRainbowRing(origin, antiSolar, u, v, secondary + offset, dist, c, 0.10 * gaussian(offset, 0.7), 3.0);
                }
            }
        }
    }

    if (S.internalReflections >= 2) {
        glLineWidth(10.0f);
        for (double angle = 45.0; angle < 50.0; angle += 0.7) {
            drawRainbowRing(origin, antiSolar, u, v, angle, dist, Color(0,0,0), 0.22, 10.0);
        }
        glLineWidth(1.0f);
    }

    drawLine(origin, origin + antiSolar * 5.0, Color(1.0, 1.0, 1.0), 2.0, 0.8);
}

static void drawLight() {
    Vec3 sunPos = Cam.pos - lightDir * 10.0 + Vec3(0, 4, 0);
    drawSphereMesh(sunPos, 0.35, Color(1.0, 0.88, 0.35), 1.0, false);
    for (int i = -2; i <= 2; ++i) {
        Vec3 offset(i * 0.35, 0, 0);
        drawLine(sunPos + offset, sunPos + offset + lightDir * 4.0, Color(1.0, 0.9, 0.45), 1.5, 0.55);
    }
}

static void drawAngularGraph() {
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0, S.width, 0, S.height);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glDisable(GL_DEPTH_TEST);

    int x0 = 80, y0 = 90, w = S.width - 160, h = 250;
    glColor4d(1,1,1,0.8);
    glBegin(GL_LINE_LOOP);
    glVertex2i(x0, y0); glVertex2i(x0+w, y0); glVertex2i(x0+w, y0+h); glVertex2i(x0, y0+h);
    glEnd();

    int count = S.whiteLight ? 12 : 1;
    for (int wlIndex = 0; wlIndex < count; ++wlIndex) {
        double wl = S.whiteLight ? 390.0 + wlIndex * (380.0 / 11.0) : S.wavelength;
        Color c = wavelengthToRGB(wl);
        glColor(c, 0.95);
        glBegin(GL_LINE_STRIP);
        for (int px = 0; px <= w; ++px) {
            double angle = 35.0 + 25.0 * px / w;
            double I = modelIntensity(angle, rainbowAnglePrimary(wl), wl, 1);
            if (S.internalReflections >= 2) I += 0.45 * modelIntensity(angle, rainbowAngleSecondary(wl), wl, 2);
            I = clampd(I, 0.0, 1.4);
            double yy = y0 + h * clampd(I / 1.4, 0, 1);
            glVertex2d(x0 + px, yy);
        }
        glEnd();
    }

    glEnable(GL_DEPTH_TEST);
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
}

static std::string modelName() {
    if (S.model == 0) return "geometric optics";
    if (S.model == 1) return "wave approximation";
    return "Mie-like approximation";
}

static std::string viewName() {
    if (S.viewMode == 0) return "full 3D scene";
    if (S.viewMode == 1) return "single droplet ray tracing";
    if (S.viewMode == 2) return "observer rainbow image";
    return "angular intensity graph";
}

static void printStatus(bool force = false) {
    static Uint32 last = 0;
    Uint32 now = SDL_GetTicks();
    if (!force && now - last < 1200) return;
    last = now;
    std::printf("model=%s | view=%s | observer=(%.2f, %.2f, %.2f) | R=%.2f | wl=%.0f nm | n=%.4f | %s | reflections=%d\n",
                modelName().c_str(), viewName().c_str(),
                Cam.pos.x, Cam.pos.y, Cam.pos.z,
                S.dropletRadius, S.wavelength, waterRefractiveIndex(S.wavelength),
                S.whiteLight ? "white" : "mono", S.internalReflections);
    std::fflush(stdout);
}

static void setCameraProjection() {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(S.fov, (double)S.width / (double)S.height, 0.05, 200.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    Vec3 f = cameraForward();
    Vec3 u = cameraUp();
    Vec3 center = Cam.pos + f;
    gluLookAt(Cam.pos.x, Cam.pos.y, Cam.pos.z,
              center.x, center.y, center.z,
              u.x, u.y, u.z);
}

static void display(SDL_Window* window) {
    glViewport(0, 0, S.width, S.height);
    glClearColor(0.025f, 0.028f, 0.04f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    setCameraProjection();

    if (S.viewMode == 0) {
        drawGrid(12.0, 24);
        drawLight();
        drawRainVolume();
        drawRainbowCones();
        if (S.showRays) drawSingleDropMode();
    }
    else if (S.viewMode == 1) {
        drawGrid(8.0, 16);
        drawLight();
        drawSingleDropMode();
    }
    else if (S.viewMode == 2) {
        drawGrid(12.0, 24);
        drawRainbowCones();
    }
    else if (S.viewMode == 3) {
        drawAngularGraph();
    }

    SDL_GL_SwapWindow(window);
}

static void resetScene() {
    Cam.pos = Vec3(0.0, 1.2, 8.0);
    Cam.yaw = 180.0;
    Cam.pitch = -8.0;
    lightDir = normalize(Vec3(-0.9, -0.25, -0.35));
    S.dropletRadius = 1.0;
    S.wavelength = 550.0;
    S.model = 0;
    S.viewMode = 0;
    S.whiteLight = true;
    S.showRays = true;
    S.showNormals = true;
    S.internalReflections = 2;
}

static void updateMotion(double dt) {
    double speed = 4.0 * dt;
    Vec3 f = cameraForward();
    Vec3 r = cameraRight();
    if (keyState[SDL_SCANCODE_W]) Cam.pos += f * speed;
    if (keyState[SDL_SCANCODE_S]) Cam.pos -= f * speed;
    if (keyState[SDL_SCANCODE_A]) Cam.pos -= r * speed;
    if (keyState[SDL_SCANCODE_D]) Cam.pos += r * speed;
    if (keyState[SDL_SCANCODE_Q]) Cam.pos.y -= speed;
    if (keyState[SDL_SCANCODE_E]) Cam.pos.y += speed;
}

static void handleSpecialKey(SDL_Keycode key) {
    Vec3 horizontalRight = normalize(cross(lightDir, Vec3(0,1,0)));
    if (length(horizontalRight) < 1e-6) horizontalRight = Vec3(1,0,0);
    if (key == SDLK_LEFT)  lightDir = normalize(lightDir - horizontalRight * 0.05);
    if (key == SDLK_RIGHT) lightDir = normalize(lightDir + horizontalRight * 0.05);
    if (key == SDLK_UP)    lightDir = normalize(lightDir + Vec3(0, 0.04, 0));
    if (key == SDLK_DOWN)  lightDir = normalize(lightDir - Vec3(0, 0.04, 0));
}

static void handleKeyDown(SDL_Keycode key) {
    switch (key) {
        case SDLK_ESCAPE: std::exit(0); break;
        case SDLK_PLUS:
        case SDLK_EQUALS:
        case SDLK_KP_PLUS: S.dropletRadius = clampd(S.dropletRadius + 0.05, 0.15, 4.0); printStatus(true); break;
        case SDLK_MINUS:
        case SDLK_UNDERSCORE:
        case SDLK_KP_MINUS: S.dropletRadius = clampd(S.dropletRadius - 0.05, 0.15, 4.0); printStatus(true); break;
        case SDLK_LEFTBRACKET: S.wavelength = clampd(S.wavelength - 5.0, 380.0, 780.0); printStatus(true); break;
        case SDLK_RIGHTBRACKET: S.wavelength = clampd(S.wavelength + 5.0, 380.0, 780.0); printStatus(true); break;
        case SDLK_1: S.model = 0; printStatus(true); break;
        case SDLK_2: S.model = 1; printStatus(true); break;
        case SDLK_3: S.model = 2; printStatus(true); break;
        case SDLK_b: S.whiteLight = !S.whiteLight; printStatus(true); break;
        case SDLK_v: S.viewMode = (S.viewMode + 1) % 4; printStatus(true); break;
        case SDLK_l: S.showRays = !S.showRays; printStatus(true); break;
        case SDLK_n: S.showNormals = !S.showNormals; printStatus(true); break;
        case SDLK_p: S.paused = !S.paused; printStatus(true); break;
        case SDLK_r: resetScene(); printStatus(true); break;
        case SDLK_0: S.internalReflections = 0; printStatus(true); break;
        case SDLK_9: S.internalReflections = 1; printStatus(true); break;
        case SDLK_8: S.internalReflections = 2; printStatus(true); break;
        case SDLK_LEFT:
        case SDLK_RIGHT:
        case SDLK_UP:
        case SDLK_DOWN: handleSpecialKey(key); printStatus(true); break;
        default: break;
    }
}

int main(int argc, char** argv) {
    (void)argc;
    (void)argv;

    resetScene();

    if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER) != 0) {
        std::fprintf(stderr, "SDL_Init failed: %s\n", SDL_GetError());
        return 1;
    }

    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 2);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 1);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);

    SDL_Window* window = SDL_CreateWindow(
        "Rainbow Simulation — Arch Linux Wayland SDL2 OpenGL",
        SDL_WINDOWPOS_CENTERED,
        SDL_WINDOWPOS_CENTERED,
        S.width,
        S.height,
        SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE | SDL_WINDOW_SHOWN
    );

    if (!window) {
        std::fprintf(stderr, "SDL_CreateWindow failed: %s\n", SDL_GetError());
        SDL_Quit();
        return 1;
    }

    SDL_GLContext ctx = SDL_GL_CreateContext(window);
    if (!ctx) {
        std::fprintf(stderr, "SDL_GL_CreateContext failed: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    SDL_GL_SetSwapInterval(1);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    const char* driver = SDL_GetCurrentVideoDriver();
    std::printf("Current SDL video driver: %s\n", driver ? driver : "unknown");
    std::printf("OpenGL renderer: %s\n", glGetString(GL_RENDERER));
    std::printf("OpenGL version: %s\n", glGetString(GL_VERSION));
    std::printf("Controls: WASD/QE move, mouse drag rotate, arrows light, 1/2/3 model, V view, B white/mono, +/- radius, [] wavelength, Esc quit.\n");
    printStatus(true);

    bool running = true;
    Uint64 lastCounter = SDL_GetPerformanceCounter();

    while (running) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) running = false;
            else if (e.type == SDL_WINDOWEVENT && e.window.event == SDL_WINDOWEVENT_SIZE_CHANGED) {
                S.width = std::max(1, e.window.data1);
                S.height = std::max(1, e.window.data2);
                glViewport(0, 0, S.width, S.height);
            }
            else if (e.type == SDL_KEYDOWN && !e.key.repeat) {
                if (e.key.keysym.scancode < SDL_NUM_SCANCODES) keyState[e.key.keysym.scancode] = true;
                handleKeyDown(e.key.keysym.sym);
            }
            else if (e.type == SDL_KEYUP) {
                if (e.key.keysym.scancode < SDL_NUM_SCANCODES) keyState[e.key.keysym.scancode] = false;
            }
            else if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_LEFT) {
                mouseDown = true;
                lastMouseX = e.button.x;
                lastMouseY = e.button.y;
            }
            else if (e.type == SDL_MOUSEBUTTONUP && e.button.button == SDL_BUTTON_LEFT) {
                mouseDown = false;
            }
            else if (e.type == SDL_MOUSEMOTION && mouseDown) {
                int dx = e.motion.x - lastMouseX;
                int dy = e.motion.y - lastMouseY;
                lastMouseX = e.motion.x;
                lastMouseY = e.motion.y;
                Cam.yaw += dx * 0.25;
                Cam.pitch -= dy * 0.25;
                Cam.pitch = clampd(Cam.pitch, -89.0, 89.0);
            }
        }

        Uint64 nowCounter = SDL_GetPerformanceCounter();
        double dt = double(nowCounter - lastCounter) / double(SDL_GetPerformanceFrequency());
        lastCounter = nowCounter;
        dt = clampd(dt, 0.0, 0.05);

        updateMotion(dt);
        if (!S.paused) animTime += dt;

        display(window);
        printStatus(false);
    }

    SDL_GL_DeleteContext(ctx);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
