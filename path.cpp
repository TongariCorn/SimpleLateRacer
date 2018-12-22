#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"

#include <array>
#include <initializer_list>
#include <iostream>
#include <cmath>
#include <memory>
#include <vector>
#include <limits>
#include <optional>

constexpr double Pi = std::ldexp(std::acos(0.0), 1);

inline std::optional<std::pair<float, float>> quadratic(float _a, float _b, float _c) {
    double d = (double)_b * (double)_b - 4 * (double)_a * (double)_c;
    if (d < 0) return std::nullopt;
    double rootD = std::sqrt(d);
    std::pair<float, float> result = std::make_pair<float, float>((-_b-rootD)/(2*_a), (-_b+rootD)/(2*_a));
    if (result.second < result.first) {
        float temp = result.first;
        result.first = result.second;
        result.second = temp;
    }
    return result;
}

// Vector
template<class T, size_t N>
class Vector {
    std::array<T, N> c;

public:
    Vector(T _t = 0.0) { c.fill(_t); }
    Vector(const std::initializer_list<T> _list) {
        size_t n = 0;
        for (auto&& iter = begin(_list); iter != end(_list); iter++) {
            c[n] = *iter;
            if (++n == N) return;
        }
    }
    Vector(const Vector& _v) : c(_v.c) {}
    Vector(const Vector&& _v) : c(std::move(_v.c)) {}

    Vector operator+(const Vector& _v) const {
        Vector v;
        for (size_t i = 0; i < N; i++) v.c[i] = c.at(i) + _v.c.at(i);
        return v;
    }

    Vector operator-() const {
        Vector v;
        for (size_t i = 0; i < N; i++) v.c[i] = -c.at(i);
        return v;
    }

    Vector operator-(const Vector& _v) const {
        Vector v;
        for (size_t i = 0; i < N; i++) v.c[i] = c.at(i) - _v.c.at(i);
        return v;
    }

    Vector operator*(T _t) const {
        Vector v;
        for (size_t i = 0; i < N; i++) v.c[i] = c.at(i) * _t;
        return v;
    }

    Vector operator*(const Vector& _v) const {
        Vector v;
        for (size_t i = 0; i < N; i++) v.c[i] = c.at(i) * _v.c.at(i);
        return v;
    }

    Vector operator/(T _t) const {
        Vector v;
        for (size_t i = 0; i < N; i++) v.c[i] = c.at(i) / _t;
        return v;
    }

    Vector& operator=(const Vector& _v) {
        c = _v.c;
        return *this;
    }

    Vector& operator=(Vector&& _v) {
        c = std::move(_v.c);
        return *this;
    }

    bool operator==(const Vector& _v) const {
        for (size_t i = 0; i < N; i++) {
            if (c.at(i) != _v.c.at(i)) return false;
        }
        return true;
    }

    T& operator[](size_t i) { return c[i]; }
    T at(size_t i) const { return c.at(i); }

    T length() const {
        return std::sqrt(dot(*this));
    }

    Vector normalize() const {
        Vector v = *this;
        return v / v.length();
    }

    T dot(const Vector& _v) const {
        T sum = 0.0;
        for (size_t i = 0; i < N; i++) sum += c.at(i) * _v.c.at(i);
        return sum;
    }

    void print(std::ostream& _ostream) const {
        _ostream<<"Vector<"<<N<<"> { ";
        for (size_t i = 0; i < N; i++) {
            _ostream<<c.at(i)<<", ";
        }
        _ostream<<"}"<<std::endl;
    }
};

using Vector3f = Vector<float, 3>;

Vector3f cross(const Vector3f& _v1, const Vector3f& _v2) {
    return Vector3f({_v1.at(1) * _v2.at(2) - _v1.at(2) * _v2.at(1),
                     _v1.at(2) * _v2.at(0) - _v1.at(0) * _v2.at(2),
                     _v1.at(0) * _v2.at(1) - _v1.at(1) * _v2.at(0)});
}


//  Ray
class Ray {
    Vector3f d, o;

public:
    Ray() = default;
    Ray(const Vector3f& _o, const Vector3f& _d) : o(_o), d(_d) {}

    Vector3f operator()(float _t) const { return o + d * _t; }
    Vector3f getOrigin() const { return o; }
    Vector3f getDir() const { return d; }
    
    void print(std::ostream& _ostream) const {
        _ostream<<"Ray { origin: { "<<o.at(0)<<", "<<o.at(1)<<", "<<o.at(2)<<"}, dir: { "<<d.at(0)<<", "<<d.at(1)<<", "<<d.at(2)<<"} }"<<std::endl;
    }
};


// Matrix
template<class T, size_t M, size_t N>
class Matrix {
    Vector<Vector<T, N>, M> m;

public:
    Matrix(T _t = 0.0) : m(_t) {}
    Matrix(const Matrix& _m) : m(_m.m) {}
    Matrix(Matrix&& _m) : m(std::move(_m.m)) {}
    Matrix(const std::initializer_list<Vector<T, N>> _list) : m(_list) {}

    Matrix<T, N, M> transpose() const {
        Matrix<T, N, M> t_m;
        for (size_t i = 0; i < M; i++) {
            for (size_t j = 0; j < N; j++) {
                t_m[j][i] = m.at(i).at(j);
            }
        }
        return t_m;
    }

    template<size_t O>
    Matrix<T, M, O> operator*(const Matrix<T, N, O>& _m) const {
        Matrix<T, O, N> t_m = _m.transpose();
        Matrix<T, M, O> new_m;
        for (size_t i = 0; i < M; i++) {
            for (size_t j = 0; j < O; j++) {
                new_m[i][j] = m.at(i).dot(t_m.at(j));
            }
        }
        return new_m;
    }

    Vector<T, M> operator*(const Vector<T, N>& _v) const {
        Vector<T, M> v;
        for (size_t i = 0; i < M; i++) v[i] = m.at(i).dot(_v);
        return v;
    }

    Vector<T, N>& operator[](size_t _i) { return m[_i]; }
    Vector<T, N> at(size_t _i) const { return m.at(_i); }

    void print(std::ostream& _ostream) const {                                                                                                                                                                       
        _ostream<<"Matrix<"<<M<<", "<<N<<"> {"<<std::endl;
        for (int i = 0; i < M; i++) {
            _ostream<<"    ";
            for (int j = 0; j < N - 1; j++) {
                _ostream<<m.at(i).at(j)<<", ";
            }
            if (i != M - 1) _ostream<<m.at(i).at(N - 1)<<", "<<std::endl;
            else _ostream<<m.at(M - 1).at(N - 1)<<std::endl;
        }
        _ostream<<"}"<<std::endl;
    }

};

using Mat4x4f = Matrix<float, 4, 4>;

// Transform
class Transform {
    Mat4x4f m, inv_m;

public:
    Transform() = default;
    Transform(const Mat4x4f& _m, const Mat4x4f& _inv_m) : m(_m), inv_m(_inv_m) {}
    Transform(Mat4x4f&& _m, Mat4x4f&& _inv_m) : m(std::move(_m)), inv_m(std::move(_inv_m)) {}
    Transform operator*(const Transform& _t) const { return Transform(m * _t.m, _t.inv_m * inv_m); }

    Vector3f act(const Vector3f& _v, float w) const {
        Vector<float, 4> v{_v.at(0), _v.at(1), _v.at(2), w};
        v = m * v;
        return Vector3f{v.at(0), v.at(1), v.at(2)};
    }

    Vector3f inv(const Vector3f& _v, float w) const {
        Vector<float, 4> v{_v.at(0), _v.at(1), _v.at(2), w};
        v = inv_m * v;
        return Vector3f{v.at(0), v.at(1), v.at(2)};
    }

    Ray act(const Ray& _r) const {
        return Ray(act(_r.getOrigin(), 1.0), act(_r.getDir(), 0.0));
    }

    Ray inv(const Ray& _r) const {
        return Ray(inv(_r.getOrigin(), 1.0), inv(_r.getDir(), 0.0));
    }

    void print(std::ostream& _ostream) const {
        _ostream<<"Transform {"<<std::endl;
        _ostream<<"    m :"<<std::endl;
        m.print(_ostream);
        _ostream<<"    inv :"<<std::endl;
        inv_m.print(_ostream);
        _ostream<<"}"<<std::endl;
    }
};

inline Mat4x4f identity() {
    return Mat4x4f{ {1.0, 0.0, 0.0, 0.0},
                    {0.0, 1.0, 0.0, 0.0},
                    {0.0, 0.0, 1.0, 0.0},
                    {0.0, 0.0, 0.0, 1.0}};
}

inline Transform translate(const Vector3f& _delta) {
    Mat4x4f m(identity());
    m[0][3] = _delta.at(0);
    m[1][3] = _delta.at(1);
    m[2][3] = _delta.at(2);

    Mat4x4f inv_m(identity());
    inv_m[0][3] = -_delta.at(0);
    inv_m[1][3] = -_delta.at(1);
    inv_m[2][3] = -_delta.at(2);
    return Transform(m, inv_m);
}    

inline Transform scale(float _x, float _y, float _z) {
    Mat4x4f m;
    m[0][0] = _x;
    m[1][1] = _y;
    m[2][2] = _z;
    m[3][3] = 1.0f;
    
    Mat4x4f inv_m;
    inv_m[0][0] = 1.0f / _x;
    inv_m[1][1] = 1.0f / _y;
    inv_m[2][2] = 1.0f / _z;
    inv_m[3][3] = 1.0f;
    return Transform(m, inv_m);
}
        
inline Transform rotateX(float _theta) {
    Mat4x4f m;
    m[0][0] = 1.0f;
    m[1][1] = std::cos(_theta);
    m[1][2] = -std::sin(_theta);
    m[2][1] = std::sin(_theta);   
    m[2][2] = std::cos(_theta);
    m[3][3] = 1.0f;
        
    return Transform(m, m.transpose());
}

inline Transform rotateY(float _theta) {
    Mat4x4f m;
    m[0][0] = std::cos(_theta);
    m[0][2] = std::sin(_theta);
    m[1][1] = 1.0f;
    m[2][0] = -std::sin(_theta);
    m[2][2] = std::cos(_theta);
    m[3][3] = 1.0f;
    
    return Transform(m, m.transpose());
}
            
inline Transform rotateZ(float _theta) {
    Mat4x4f m;
    m[0][0] = std::cos(_theta);
    m[0][1] = -std::sin(_theta);
    m[1][0] = std::sin(_theta);
    m[1][1] = std::cos(_theta);
    m[2][2] = 1.0f;
    m[3][3] = 1.0f;
    
    return Transform(m, m.transpose());
}


class BRDF;

// SurfaceInteraction
struct SurfaceInteraction {
    std::unique_ptr<BRDF> brdf;
    Vector3f point;
    Vector3f normal;
    Vector3f wout;
    Vector3f tangent;

    void actedBy(Transform _t) {
        point = _t.act(point, 1.0);
        normal = _t.act(normal, 0.0).normalize();
        wout = _t.act(wout, 0.0).normalize();
        tangent = _t.act(tangent, 0.0).normalize();
    }

    Vector3f worldToLocal(const Vector3f& _v) const {
        return Vector3f{_v.dot(tangent), _v.dot(normal), _v.dot(cross(tangent, normal))};
    }
};


// Shape
class Shape {
public:
    Shape() = default;
    virtual ~Shape() = default;
    virtual bool intersect(const Ray&, float*, SurfaceInteraction*) const = 0;
};

class Plane : public Shape {
    float width, height;
public:
    Plane(float _w = 1.0, float _h = 1.0) : width(_w), height(_h) {}
    ~Plane() = default;
    bool intersect(const Ray& _ray, float* _t, SurfaceInteraction* _si) const {
        Vector3f o = _ray.getOrigin();
        Vector3f d = _ray.getDir();
        if (d.at(2) == 0.0) return false;

        float t = std::abs(o.at(2)) / d.at(2);
        if (t < 0.0000000001) return false;
        Vector3f point = _ray(t);
        if (point.at(0) < -(width/2.0) || (width/2.0) < point.at(0)
                || point.at(1) < -(height/2.0) || (height/2.0) < point.at(1)) return false;
        *_t = t;

        _si->point = point;
        _si->normal = Vector3f{0.0f, 0.0f, 1.0f} * (o.at(2) / std::abs(o.at(2)));
        _si->tangent = Vector3f{1.0f, 0.0f, 0.0f};
        _si->wout = -d.normalize();
        return true;
    }
};

// |O + tD| = r, (dx^2+dy^2+dz^2)t^2 + 2(dxox+dyoy+dzoz)t + ox^2+oy^2+oz^2-r^2=0
class Sphere : public Shape {
    float r;
public:
    Sphere(float _r = 1.0) : r(_r) {}
    ~Sphere() = default;
    bool intersect(const Ray& _ray, float* _t, SurfaceInteraction* _si) const {
        Vector3f o = _ray.getOrigin();
        Vector3f d = _ray.getDir();
        float a = d.at(0) * d.at(0) + d.at(1) * d.at(1) + d.at(2) * d.at(2);
        float b = 2.0 * (d.at(0) * o.at(0) + d.at(1) * o.at(1) + d.at(2) * o.at(2));
        float c = o.at(0) * o.at(0) + o.at(1) * o.at(1) + o.at(2) * o.at(2) - r * r;
        auto r_op = quadratic(a, b, c);
        if (!r_op) return false;

        auto ts = r_op.value();
        if (ts.second < 0.0000000001f) return false;
        *_t = ts.first;

        _si->point = _ray(ts.first);
        _si->normal = _si->point;
        if (_si->normal.at(0) == 0.0 && _si->normal.at(2) == 0.0) _si->tangent = Vector3f{1.0f, 0.0f, 0.0f};
        else _si->tangent = Vector3f{_si->normal.at(2), 0.0f, -_si->normal.at(0)}.normalize();
        _si->wout = -d.normalize();
        return true;
    }
};


using RGB = Vector3f;

// BRDF
class BRDF {
public:
    BRDF() = default;
    virtual ~BRDF() = default;
    virtual RGB f(const Vector3f& _wo, const Vector3f& _wi) const = 0;
};

class LambertianReflection : public BRDF {
    const RGB color;

public:
    LambertianReflection(const RGB& _color) : color(_color) {}
    ~LambertianReflection() = default;
    RGB f(const Vector3f& _wo, const Vector3f& _wi) const {
        return color / Pi;
    }
};


// Primitive
class Primitive {
    std::shared_ptr<Shape> shape;
    Transform worldToLocal;
    RGB color;

public:
    Primitive(const std::shared_ptr<Shape> _shape, const Transform& _t, const RGB& _c) : shape(_shape), worldToLocal(_t), color(_c) {}
    bool intersect(const Ray& _ray, float* _t, SurfaceInteraction* _si) const {
        Ray r = worldToLocal.inv(_ray);

        bool result = false;
        if (shape) result = shape->intersect(r, _t, _si);
        if (result) {
            _si->actedBy(worldToLocal);
        }
        _si->brdf.reset(new LambertianReflection(color));
        return result;
    }
};


// Light
class Light {
public:
    Light() = default;
    virtual ~Light() = default;
    virtual RGB sampleLi(const Vector3f& _p, Ray* _wi, float* _pdf) const = 0;
};

class PointLight : public Light {
    Vector3f pos;
    RGB color;

public:
    PointLight(const Vector3f& _p, const RGB& _c) : pos(_p), color(_c) {}
    ~PointLight() {}
    RGB sampleLi(const Vector3f& _p, Ray* _wi, float* _pdf) const {
        *_wi = Ray(pos, _p - pos);
        *_pdf = 1.0;
        float length = _wi->getDir().length();
        return color / (length * length);
    }
};


// Scene
class Scene {
    std::shared_ptr<Light> light;
    std::vector<std::shared_ptr<Primitive>> primitives;

public:
    Scene(std::shared_ptr<Light> _light) : light(_light) {}

    bool intersect(const Ray& _ray, SurfaceInteraction* _si, float* _t) const {
        float result = false;
        float min_t = std::numeric_limits<float>::max();
        for (auto&& primitive : primitives) {
            float t = std::numeric_limits<float>::max();
            SurfaceInteraction si;
            if (primitive->intersect(_ray, &t, &si)) {
                result = true;
                if (t < min_t) {
                    min_t = t;
                    *_si = std::move(si);
                }
            }
        }
        *_t = min_t;
        return result;
    }
    void addPrimitive(const std::shared_ptr<Primitive> _primitive) { primitives.push_back(_primitive); }
    const std::shared_ptr<Light> getLight() const { return light; }
};

// Integrator
class Integrator {
public:
    virtual ~Integrator() = default;
    void render(const Scene& _scene) const {
        const int width = 512;
        const int height = 512;
        char pixels[width * height * 3];
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                float x = ((float)j * 2.0 - (float)width) / (float)width * 2.0;
                float y = ((float)i * 2.0 - (float)height) / (float)height * 2.0;
                Ray ray(Vector3f{0.0f, 0.0f, -4.0f}, Vector3f{x, -y, 2.4f}.normalize());
                RGB rgb = li(ray, _scene);

                int pos = (i * width + j) * 3;
                pixels[pos]   = (unsigned char)(255.f * std::min(rgb.at(0), 1.0f));
                pixels[pos+1] = (unsigned char)(255.f * std::min(rgb.at(1), 1.0f));
                pixels[pos+2] = (unsigned char)(255.f * std::min(rgb.at(2), 1.0f));
            }
        }
        stbi_write_png("out.png", width, height, STBI_rgb, pixels, 0);
    }
    
    RGB estimateDirect(const SurfaceInteraction& _si, const std::shared_ptr<Light> _light, const Scene& _scene) const {
        RGB ld{0.0, 0.0, 0.0};

        Ray wi;
        float lightPdf;
        RGB l = _light->sampleLi(_si.point, &wi, &lightPdf);
        if (lightPdf > 0.0f) {
            float t;
            SurfaceInteraction si;
            if (!(_scene.intersect(wi, &si, &t) && t <= 0.999)) {
                auto l_wout = _si.worldToLocal(_si.wout);
                auto l_wi = _si.worldToLocal((-wi.getDir().normalize()));
                RGB f = _si.brdf->f(l_wout, l_wi) * std::abs(l_wi.dot(Vector3f{0.0, 1.0, 0.0}));
                ld = f * l / lightPdf;
            }
        }

        return ld;
    }

    virtual RGB li(const Ray&, const Scene&) const = 0;
};

class DirectIntegrator : public Integrator {
    public:
    RGB li(const Ray& _ray, const Scene& _scene) const {
        SurfaceInteraction si;
        float t;
        if (_scene.intersect(_ray, &si, &t)) {
            auto light = _scene.getLight();
            if (light) return estimateDirect(si, light, _scene);
        }
        return RGB{0.0, 0.0, 0.0};
    }
};

int main() {
    Transform plane1_wtl = translate({0.0f, 0.0f, 3.5f});
    std::shared_ptr<Shape> plane1(new Plane(8.0, 8.0));
    std::shared_ptr<Primitive> p_plane1(new Primitive(plane1, plane1_wtl, RGB{1.0, 1.0, 1.0}));
    
    Transform plane2_wtl = translate({-3.5f, 0.0f, 0.0f}) * rotateY(-Pi / 2.0f);
    std::shared_ptr<Shape> plane2(new Plane(8.0, 8.0));
    std::shared_ptr<Primitive> p_plane2(new Primitive(plane2, plane2_wtl, RGB{1.0, 0.0, 0.0}));

    Transform plane3_wtl = translate({3.5f, 0.0f, 0.0f}) * rotateY(Pi / 2.0f);
    std::shared_ptr<Shape> plane3(new Plane(8.0, 8.0));
    std::shared_ptr<Primitive> p_plane3(new Primitive(plane3, plane3_wtl, RGB{0.0, 1.0, 0.0}));

    Transform plane4_wtl = translate({0.0f, 3.5f, 0.0f}) * rotateX(-Pi / 2.0f);
    std::shared_ptr<Shape> plane4(new Plane(8.0, 8.0));
    std::shared_ptr<Primitive> p_plane4(new Primitive(plane4, plane4_wtl, RGB{1.0, 1.0, 1.0}));

    Transform plane5_wtl = translate({0.0f, -3.5f, 0.0f}) * rotateX(Pi / 2.0f);
    std::shared_ptr<Shape> plane5(new Plane(8.0, 8.0));
    std::shared_ptr<Primitive> p_plane5(new Primitive(plane5, plane5_wtl, RGB{1.0, 1.0, 1.0}));

    Transform sphere1_wtl = translate({-1.5f, -2.8f, 1.5f});
    std::shared_ptr<Shape> sphere1(new Sphere(0.7));
    std::shared_ptr<Primitive> p_sphere1(new Primitive(sphere1, sphere1_wtl, RGB{1.0, 1.0, 0.0}));

    Transform sphere2_wtl = translate({2.0f, -2.5f, 1.0f});
    std::shared_ptr<Shape> sphere2(new Sphere(1.0));
    std::shared_ptr<Primitive> p_sphere2(new Primitive(sphere2, sphere2_wtl, RGB{0.6, 0.7, 0.8}));

    std::shared_ptr<Light> light(new PointLight(Vector3f{0.0f, 3.4f, 0.0f}, RGB{40.0f, 35.0f, 30.0f}));
    Scene scene(light);
    scene.addPrimitive(p_plane1);
    scene.addPrimitive(p_plane2);
    scene.addPrimitive(p_plane3);
    scene.addPrimitive(p_plane4);
    scene.addPrimitive(p_plane5);
    scene.addPrimitive(p_sphere1);
    scene.addPrimitive(p_sphere2);

    DirectIntegrator integrator;
    integrator.render(scene);

    return 0;
}
