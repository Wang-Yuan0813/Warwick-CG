// CG_Day1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <string>
#define SQ(x) ((x)*(x))
template<typename T>
static T lerp(const T a, const T b, float t) {
    return a * (1.0f - t) + (b * t);
}

template <typename T>
T max(T a, T b) {
    return (a > b) ? a : b;
}
template <typename T>
T min(T a, T b) {
    return (a < b) ? a : b;
}

class Vec3 {
public:
    union {
        float v[3];
        struct { float x, y, z; };
    };
    Vec3(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}
    Vec3() : x(0), y(0), z(0) {}
    Vec3 operator+(const Vec3& pVec) const {
        return Vec3(v[0] + pVec.v[0], v[1] + pVec.v[1], v[2] + pVec.v[2]);
    }
    Vec3 operator+(const float val) const {
        return Vec3(v[0] + val, v[1] + val, v[2] + val);
    }
    Vec3 operator-(const Vec3& pVec) const {
        return Vec3(v[0] - pVec.v[0], v[1] - pVec.v[1], v[2] - pVec.v[2]);
    }
    Vec3 operator-(const float val) const {
        return Vec3(v[0] - val, v[1] - val, v[2] - val);
    }
    Vec3 operator-() const {
        return Vec3(-v[0], -v[1], -v[2]);
    }
    Vec3& operator+=(const Vec3& pVec) {
        v[0] += pVec.v[0];
        v[1] += pVec.v[1];
        v[2] += pVec.v[2];
        return *this;
    }
    Vec3& operator+=(const float val) {
        v[0] += val;
        v[1] += val;
        v[2] += val;
        return *this;
    }
    Vec3& operator-=(const Vec3& pVec) {
        v[0] -= pVec.v[0];
        v[1] -= pVec.v[1];
        v[2] -= pVec.v[2];
        return *this;
    }
    Vec3& operator-=(const float val) {
        v[0] -= val;
        v[1] -= val;
        v[2] -= val;
        return *this;
    }
    Vec3 operator*(const Vec3& pVec) const {
        return Vec3(v[0] * pVec.v[0], v[1] * pVec.v[1], v[2] * pVec.v[2]);
    }
    Vec3 operator*(const float val) const {
        return Vec3(v[0] * val, v[1] * val, v[2] * val);
    }
    Vec3& operator*=(const Vec3& pVec) {
        v[0] *= pVec.v[0];
        v[1] *= pVec.v[1];
        v[2] *= pVec.v[2];
        return *this;
    }
    Vec3& operator*=(const float val) {
        v[0] *= val;
        v[1] *= val;
        v[2] *= val;
        return *this;
    }
    Vec3 operator/(const Vec3& pVec) const {
        return Vec3(v[0] / pVec.v[0], v[1] / pVec.v[1], v[2] / pVec.v[2]);
    }
    Vec3 operator/(const float val) const {
        return Vec3(v[0] / val, v[1] / val, v[2] / val);
    }
    Vec3& operator/=(const Vec3& pVec) {
        v[0] /= pVec.v[0];
        v[1] /= pVec.v[1];
        v[2] /= pVec.v[2];
        return *this;
    }
    Vec3& operator/=(const float val) {
        v[0] /= val;
        v[1] /= val;
        v[2] /= val;
        return *this;
    }
    float& operator[](const int index) {
        return v[index];
    }
    float length() const {
        return sqrtf(SQ(v[0]) + SQ(v[1]) + SQ(v[2]));
    }
    float lenthSquare() const {
        return SQ(v[0]) + SQ(v[1]) + SQ(v[2]);
    }
    Vec3 normalize() {
        float len = 1.0f / length();
        return Vec3(x * len, y * len, z * len);
    }
    float normalize_GetLength() {
        float _length = length();
        float len = 1.0f / _length;
        v[0] *= len; v[1] *= len; v[2] *= len;
        return _length;
    }
    float dot(const Vec3& pVec) const {
        return v[0] * pVec.v[0] + v[1] * pVec.v[1] + v[2] * pVec.v[2];
    }
    Vec3 cross(const Vec3& v1) {
        return Vec3(v1.v[1] * v[2] - v1.v[2] * v[1],
            v1.v[2] * v[0] - v1.v[0] * v[2],
            v1.v[0] * v[1] - v1.v[1] * v[0]);
    }
    void print() const {
        std::cout << '(' << v[0] << ',' << v[1] << ',' << v[2] << ')' << std::endl;
    }
    std::string print2str() const {
        return std::string('(' + std::to_string(v[0]) + ',' + std::to_string(v[1]) + ',' + std::to_string(v[2]) + ')');
    }
    Vec3 Max(const Vec3& v1, const Vec3& v2)
    {
        return Vec3(max(v1.v[0], v2.v[0]),
            max(v1.v[1], v2.v[1]),
            max(v1.v[2], v2.v[2]));
    }
    Vec3 Min(const Vec3& v1, const Vec3& v2)
    {
        return Vec3(min(v1.v[0], v2.v[0]),
            min(v1.v[1], v2.v[1]),
            min(v1.v[2], v2.v[2]));
    }
    float Max() const
    {
        return max(x, max(y, z));
    }
    float Min() const
    {
        return min(x, min(y, z));
    }
};
//vector functions
float Dot(const Vec3& v1, const Vec3& v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}
class Vec4 {
public:
    union {
        float v[4];
        struct { float x, y, z, w; };
    };

    Vec4() : x(0), y(0), z(0), w(0) {}
    Vec4(float _x, float _y, float _z, float _w) : x(_x), y(_y), z(_z), w(_w) {}

    Vec4 operator+(const Vec4& pVec) const {
        return Vec4(v[0] + pVec.v[0], v[1] + pVec.v[1], v[2] + pVec.v[2], v[3] + pVec.v[3]);
    }
    Vec4 operator+(const float val) const {
        return Vec4(v[0] + val, v[1] + val, v[2] + val, v[3] + val);
    }
    Vec4 operator-(const Vec4& pVec) const {
        return Vec4(v[0] - pVec.v[0], v[1] - pVec.v[1], v[2] - pVec.v[2], v[3] - pVec.v[3]);
    }
    Vec4 operator-(const float val) const {
        return Vec4(v[0] - val, v[1] - val, v[2] - val, v[3] - val);
    }
    Vec4 operator-() const {
        return Vec4(-v[0], -v[1], -v[2], -v[3]);
    }
    Vec4& operator+=(const Vec4& pVec) {
        v[0] += pVec.v[0];
        v[1] += pVec.v[1];
        v[2] += pVec.v[2];
        v[3] += pVec.v[3];
        return *this;
    }
    Vec4& operator+=(const float val) {
        v[0] += val;
        v[1] += val;
        v[2] += val;
        v[3] += val;
        return *this;
    }
    Vec4& operator-=(const Vec4& pVec) {
        v[0] -= pVec.v[0];
        v[1] -= pVec.v[1];
        v[2] -= pVec.v[2];
        v[3] -= pVec.v[3];
        return *this;
    }
    Vec4& operator-=(const float val) {
        v[0] -= val;
        v[1] -= val;
        v[2] -= val;
        v[3] -= val;
        return *this;
    }
    Vec4 operator*(const Vec4& pVec) const {
        return Vec4(v[0] * pVec.v[0], v[1] * pVec.v[1], v[2] * pVec.v[2], v[3] * pVec.v[3]);
    }
    Vec4 operator*(const float val) const {
        return Vec4(v[0] * val, v[1] * val, v[2] * val, v[3] * val);
    }
    Vec4& operator*=(const Vec4& pVec) {
        v[0] *= pVec.v[0];
        v[1] *= pVec.v[1];
        v[2] *= pVec.v[2];
        v[3] *= pVec.v[3];
        return *this;
    }
    Vec4& operator*=(const float val) {
        v[0] *= val;
        v[1] *= val;
        v[2] *= val;
        v[3] *= val;
        return *this;
    }
    Vec4 operator/(const Vec4& pVec) const {
        return Vec4(v[0] / pVec.v[0], v[1] / pVec.v[1], v[2] / pVec.v[2], v[3] / pVec.v[3]);
    }
    Vec4 operator/(const float val) const {
        return Vec4(v[0] / val, v[1] / val, v[2] / val, v[3] / val);
    }
    Vec4& operator/=(const Vec4& pVec) {
        v[0] /= pVec.v[0];
        v[1] /= pVec.v[1];
        v[2] /= pVec.v[2];
        v[3] /= pVec.v[3];
        return *this;
    }
    Vec4& operator/=(const float val) {
        v[0] /= val;
        v[1] /= val;
        v[2] /= val;
        v[3] /= val;
        return *this;
    }
    float& operator[](const int index) {
        return v[index];
    }
    float length() const {
        return sqrtf(SQ(v[0]) + SQ(v[1]) + SQ(v[2]) + SQ(v[3]));
    }
    float lenthSquare() const {
        return SQ(v[0]) + SQ(v[1]) + SQ(v[2]) + SQ(v[3]);
    }
    Vec4 normalize() const {
        float L = length();
        if (L == 0.0f) return Vec4(0, 0, 0, 0);
        float len = 1.0f / L;
        return Vec4(x * len, y * len, z * len, w * len);
    }
    float normalizeGetLength() {
        float _length = length();
        if (_length == 0.0f) return 0.0f;
        float len = 1.0f / _length;
        v[0] *= len; v[1] *= len; v[2] *= len; v[3] *= len;
        return _length;
    }
    float dot(const Vec4& pVec) const {
        return v[0] * pVec.v[0] + v[1] * pVec.v[1] + v[2] * pVec.v[2] + v[3] * pVec.v[3];
    }
    // Cross uses xyz, w set to 0
    Vec4 cross(const Vec4& v1) const {
        return Vec4(v1.v[1] * v[2] - v1.v[2] * v[1],
            v1.v[2] * v[0] - v1.v[0] * v[2],
            v1.v[0] * v[1] - v1.v[1] * v[0],
            0.0f);
    }
    void print() const {
        std::cout << '(' << v[0] << ',' << v[1] << ',' << v[2] << ',' << v[3] << ')' << std::endl;
    }
    std::string print2str() const {
        return std::string('(' + std::to_string(v[0]) + ',' + std::to_string(v[1]) + ',' +
            std::to_string(v[2]) + ',' + std::to_string(v[3]) + ')');
    }
    Vec4 Max(const Vec4& v1, const Vec4& v2)
    {
        return Vec4(max(v1.v[0], v2.v[0]),
            max(v1.v[1], v2.v[1]),
            max(v1.v[2], v2.v[2]),
            max(v1.v[3], v2.v[3]));
    }
    Vec4 Min(const Vec4& v1, const Vec4& v2)
    {
        return Vec4(min(v1.v[0], v2.v[0]),
            min(v1.v[1], v2.v[1]),
            min(v1.v[2], v2.v[2]),
            min(v1.v[3], v2.v[3]));
    }
    float Max() const
    {
        return max(x, max(y, max(z, w)));
    }
    float Min() const
    {
        return min(x, min(y, min(z, w)));
    }
};

class Matrix {
public:
    union {
        float a[4][4];
        float m[16];
    };
    Matrix() {
        memset(m, 0, 16 * sizeof(float));
        a[0][0] = 1; a[1][1] = 1; a[2][2] = 1; a[3][3] = 1;
    }
    Matrix(float* _m) {
        m[0] = _m[0];   m[1] = _m[1];   m[2] = _m[2];   m[3] = _m[3];
        m[4] = _m[4];   m[5] = _m[5];   m[6] = _m[6];   m[7] = _m[7];
        m[8] = _m[8];   m[9] = _m[9];   m[10] = _m[10]; m[11] = _m[11];
        m[12] = _m[12]; m[13] = _m[13]; m[14] = _m[14]; m[15] = _m[15];
    }
    void print() const {
        for (int i = 0; i < 4; i++) {
            std::cout << '|' << '\t';
            for (int j = 0; j < 4; j++) {
                std::cout << a[i][j] << '\t';
            }
            std::cout << '|' << std::endl;
        }
    }
    Vec4 mul(const Vec4& v) {
        return Vec4(
            v.x * m[0] + v.y * m[1] + v.z * m[2] + v.w * m[3],
            v.x * m[4] + v.y * m[5] + v.z * m[6] + v.w * m[7],
            v.x * m[8] + v.y * m[9] + v.z * m[10] + v.w * m[11],
            v.x * m[12] + v.y * m[13] + v.z * m[14] + v.w * m[15]);
    }
    Vec3 mulPoint(const Vec3& v) {
        Vec3 v1 = Vec3(
            (v.x * m[0] + v.y * m[1] + v.z * m[2]) + m[3],
            (v.x * m[4] + v.y * m[5] + v.z * m[6]) + m[7],
            (v.x * m[8] + v.y * m[9] + v.z * m[10]) + m[11]);
    }
    Vec3 mulVec(const Vec3& v) {
        return Vec3(
            (v.x * m[0] + v.y * m[1] + v.z * m[2]),
            (v.x * m[4] + v.y * m[5] + v.z * m[6]),
            (v.x * m[8] + v.y * m[9] + v.z * m[10]));
    }
    Matrix mul(const Matrix& matrix) const {
        Matrix ret;
        ret.m[0] = m[0] * matrix.m[0] + m[1] * matrix.m[4] + m[2] * matrix.m[8] + m[3] * matrix.m[12];
        ret.m[1] = m[0] * matrix.m[1] + m[1] * matrix.m[5] + m[2] * matrix.m[9] + m[3] * matrix.m[13];
        ret.m[2] = m[0] * matrix.m[2] + m[1] * matrix.m[6] + m[2] * matrix.m[10] + m[3] * matrix.m[14];
        ret.m[3] = m[0] * matrix.m[3] + m[1] * matrix.m[7] + m[2] * matrix.m[11] + m[3] * matrix.m[15];
        ret.m[4] = m[4] * matrix.m[0] + m[5] * matrix.m[4] + m[6] * matrix.m[8] + m[7] * matrix.m[12];
        ret.m[5] = m[4] * matrix.m[1] + m[5] * matrix.m[5] + m[6] * matrix.m[9] + m[7] * matrix.m[13];
        ret.m[6] = m[4] * matrix.m[2] + m[5] * matrix.m[6] + m[6] * matrix.m[10] + m[7] * matrix.m[14];
        ret.m[7] = m[4] * matrix.m[3] + m[5] * matrix.m[7] + m[6] * matrix.m[11] + m[7] * matrix.m[15];
        ret.m[8] = m[8] * matrix.m[0] + m[9] * matrix.m[4] + m[10] * matrix.m[8] + m[11] * matrix.m[12];
        ret.m[9] = m[8] * matrix.m[1] + m[9] * matrix.m[5] + m[10] * matrix.m[9] + m[11] * matrix.m[13];
        ret.m[10] = m[8] * matrix.m[2] + m[9] * matrix.m[6] + m[10] * matrix.m[10] + m[11] * matrix.m[14];
        ret.m[11] = m[8] * matrix.m[3] + m[9] * matrix.m[7] + m[10] * matrix.m[11] + m[11] * matrix.m[15];
        ret.m[12] = m[12] * matrix.m[0] + m[13] * matrix.m[4] + m[14] * matrix.m[8] + m[15] * matrix.m[12];
        ret.m[13] = m[12] * matrix.m[1] + m[13] * matrix.m[5] + m[14] * matrix.m[9] + m[15] * matrix.m[13];
        ret.m[14] = m[12] * matrix.m[2] + m[13] * matrix.m[6] + m[14] * matrix.m[10] + m[15] * matrix.m[14];
        ret.m[15] = m[12] * matrix.m[3] + m[13] * matrix.m[7] + m[14] * matrix.m[11] + m[15] * matrix.m[15];
        return ret;
    }
    float determinant() {
        Matrix inv;
        inv[0] = m[5] * m[10] * m[15] - m[5] * m[11] * m[14] - m[9] * m[6] * m[15] + m[9] * m[7] * m[14] + m[13] * m[6] * m[11] - m[13] * m[7] * m[10];
        inv[4] = -m[4] * m[10] * m[15] + m[4] * m[11] * m[14] + m[8] * m[6] * m[15] - m[8] * m[7] * m[14] - m[12] * m[6] * m[11] + m[12] * m[7] * m[10];
        inv[8] = m[4] * m[9] * m[15] - m[4] * m[11] * m[13] - m[8] * m[5] * m[15] + m[8] * m[7] * m[13] + m[12] * m[5] * m[11] - m[12] * m[7] * m[9];
        inv[12] = -m[4] * m[9] * m[14] + m[4] * m[10] * m[13] + m[8] * m[5] * m[14] - m[8] * m[6] * m[13] - m[12] * m[5] * m[10] + m[12] * m[6] * m[9];
        inv[1] = -m[1] * m[10] * m[15] + m[1] * m[11] * m[14] + m[9] * m[2] * m[15] - m[9] * m[3] * m[14] - m[13] * m[2] * m[11] + m[13] * m[3] * m[10];
        inv[5] = m[0] * m[10] * m[15] - m[0] * m[11] * m[14] - m[8] * m[2] * m[15] + m[8] * m[3] * m[14] + m[12] * m[2] * m[11] - m[12] * m[3] * m[10];
        inv[9] = -m[0] * m[9] * m[15] + m[0] * m[11] * m[13] + m[8] * m[1] * m[15] - m[8] * m[3] * m[13] - m[12] * m[1] * m[11] + m[12] * m[3] * m[9];
        inv[13] = m[0] * m[9] * m[14] - m[0] * m[10] * m[13] - m[8] * m[1] * m[14] + m[8] * m[2] * m[13] + m[12] * m[1] * m[10] - m[12] * m[2] * m[9];
        inv[2] = m[1] * m[6] * m[15] - m[1] * m[7] * m[14] - m[5] * m[2] * m[15] + m[5] * m[3] * m[14] + m[13] * m[2] * m[7] - m[13] * m[3] * m[6];
        inv[6] = -m[0] * m[6] * m[15] + m[0] * m[7] * m[14] + m[4] * m[2] * m[15] - m[4] * m[3] * m[14] - m[12] * m[2] * m[7] + m[12] * m[3] * m[6];
        inv[10] = m[0] * m[5] * m[15] - m[0] * m[7] * m[13] - m[4] * m[1] * m[15] + m[4] * m[3] * m[13] + m[12] * m[1] * m[7] - m[12] * m[3] * m[5];
        inv[14] = -m[0] * m[5] * m[14] + m[0] * m[6] * m[13] + m[4] * m[1] * m[14] - m[4] * m[2] * m[13] - m[12] * m[1] * m[6] + m[12] * m[2] * m[5];
        inv[3] = -m[1] * m[6] * m[11] + m[1] * m[7] * m[10] + m[5] * m[2] * m[11] - m[5] * m[3] * m[10] - m[9] * m[2] * m[7] + m[9] * m[3] * m[6];
        inv[7] = m[0] * m[6] * m[11] - m[0] * m[7] * m[10] - m[4] * m[2] * m[11] + m[4] * m[3] * m[10] + m[8] * m[2] * m[7] - m[8] * m[3] * m[6];
        inv[11] = -m[0] * m[5] * m[11] + m[0] * m[7] * m[9] + m[4] * m[1] * m[11] - m[4] * m[3] * m[9] - m[8] * m[1] * m[7] + m[8] * m[3] * m[5];
        inv[15] = m[0] * m[5] * m[10] - m[0] * m[6] * m[9] - m[4] * m[1] * m[10] + m[4] * m[2] * m[9] + m[8] * m[1] * m[6] - m[8] * m[2] * m[5];
        float det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];
        return det;
    }
    Matrix invert() const {
        Matrix inv;
        inv[0] = m[5] * m[10] * m[15] - m[5] * m[11] * m[14] - m[9] * m[6] * m[15] + m[9] * m[7] * m[14] + m[13] * m[6] * m[11] - m[13] * m[7] * m[10];
        inv[4] = -m[4] * m[10] * m[15] + m[4] * m[11] * m[14] + m[8] * m[6] * m[15] - m[8] * m[7] * m[14] - m[12] * m[6] * m[11] + m[12] * m[7] * m[10];
        inv[8] = m[4] * m[9] * m[15] - m[4] * m[11] * m[13] - m[8] * m[5] * m[15] + m[8] * m[7] * m[13] + m[12] * m[5] * m[11] - m[12] * m[7] * m[9];
        inv[12] = -m[4] * m[9] * m[14] + m[4] * m[10] * m[13] + m[8] * m[5] * m[14] - m[8] * m[6] * m[13] - m[12] * m[5] * m[10] + m[12] * m[6] * m[9];
        inv[1] = -m[1] * m[10] * m[15] + m[1] * m[11] * m[14] + m[9] * m[2] * m[15] - m[9] * m[3] * m[14] - m[13] * m[2] * m[11] + m[13] * m[3] * m[10];
        inv[5] = m[0] * m[10] * m[15] - m[0] * m[11] * m[14] - m[8] * m[2] * m[15] + m[8] * m[3] * m[14] + m[12] * m[2] * m[11] - m[12] * m[3] * m[10];
        inv[9] = -m[0] * m[9] * m[15] + m[0] * m[11] * m[13] + m[8] * m[1] * m[15] - m[8] * m[3] * m[13] - m[12] * m[1] * m[11] + m[12] * m[3] * m[9];
        inv[13] = m[0] * m[9] * m[14] - m[0] * m[10] * m[13] - m[8] * m[1] * m[14] + m[8] * m[2] * m[13] + m[12] * m[1] * m[10] - m[12] * m[2] * m[9];
        inv[2] = m[1] * m[6] * m[15] - m[1] * m[7] * m[14] - m[5] * m[2] * m[15] + m[5] * m[3] * m[14] + m[13] * m[2] * m[7] - m[13] * m[3] * m[6];
        inv[6] = -m[0] * m[6] * m[15] + m[0] * m[7] * m[14] + m[4] * m[2] * m[15] - m[4] * m[3] * m[14] - m[12] * m[2] * m[7] + m[12] * m[3] * m[6];
        inv[10] = m[0] * m[5] * m[15] - m[0] * m[7] * m[13] - m[4] * m[1] * m[15] + m[4] * m[3] * m[13] + m[12] * m[1] * m[7] - m[12] * m[3] * m[5];
        inv[14] = -m[0] * m[5] * m[14] + m[0] * m[6] * m[13] + m[4] * m[1] * m[14] - m[4] * m[2] * m[13] - m[12] * m[1] * m[6] + m[12] * m[2] * m[5];
        inv[3] = -m[1] * m[6] * m[11] + m[1] * m[7] * m[10] + m[5] * m[2] * m[11] - m[5] * m[3] * m[10] - m[9] * m[2] * m[7] + m[9] * m[3] * m[6];
        inv[7] = m[0] * m[6] * m[11] - m[0] * m[7] * m[10] - m[4] * m[2] * m[11] + m[4] * m[3] * m[10] + m[8] * m[2] * m[7] - m[8] * m[3] * m[6];
        inv[11] = -m[0] * m[5] * m[11] + m[0] * m[7] * m[9] + m[4] * m[1] * m[11] - m[4] * m[3] * m[9] - m[8] * m[1] * m[7] + m[8] * m[3] * m[5];
        inv[15] = m[0] * m[5] * m[10] - m[0] * m[6] * m[9] - m[4] * m[1] * m[10] + m[4] * m[2] * m[9] + m[8] * m[1] * m[6] - m[8] * m[2] * m[5];
        float det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];
        if (det == 0) {
            // Handle this case 
            std::cout << "det = 0!!" << std::endl;
        }
        det = 1.0f / det;
        for (int i = 0; i < 16; i++) {
            inv[i] = inv[i] * det;
        }
        return inv;
    }
    Matrix transpose() const {
        Matrix inv;
        inv[0] = m[0];  inv[1] = m[4];  inv[2] = m[8];  inv[3] = m[12];
        inv[4] = m[1];  inv[5] = m[5];  inv[6] = m[9];  inv[7] = m[13];
        inv[8] = m[2];  inv[9] = m[6];  inv[10] = m[10]; inv[11] = m[14];
        inv[12] = m[3]; inv[13] = m[7]; inv[14] = m[11]; inv[15] = m[15];
        return inv;
    }
    float& operator[](const int index) {
        return m[index];
    }
    Matrix operator*(const Matrix& matrix) {
        return mul(matrix);
    }
};
//translation
Vec4 Translation(const Vec4& v, float x, float y, float z) {
    Matrix m;
    m.m[3] = x; m.m[7] = y; m.m[11] = z;
    return m.mul(v);
}
//rotation
Vec4 RotateX(const Vec4& v, float theta) {
    return Vec4(
        (v.x + v.y * 0 + v.z * 0),
        (v.y * cos(theta) - v.z * sin(theta)),
        (v.y * sin(theta) + v.z * cos(theta)),
        v.w);
}
Vec4 RotateY(const Vec4& v, float theta) {
    return Vec4(
        (v.x * cos(theta) + v.y * 0 + v.z * sin(theta)),
        (v.x * 0 + v.y * 1 + v.z * 0),
        (-v.x * sin(theta) + v.y * 0 + v.z * cos(theta)),
        v.w);
}
Vec4 RotateZ(const Vec4& v, float theta) {
    return Vec4(
        (v.x * cos(theta) - v.y * sin(theta) + v.z * 0),
        (v.x * sin(theta) + v.y * cos(theta) + v.z * 0),
        (v.x * 0 + v.y * 0 + v.z * 1),
        v.w);
}
//scaling
Vec4 Scaling(const Vec4& v, float x, float y, float z) {
    Matrix m;
    m.m[0] = x; m.m[5] = y; m.m[10] = z;
    return m.mul(v);
}
class SphericalCoordinates {
public:
    float theta;
    float phi;
    float r;
    SphericalCoordinates(Vec3 v, bool isZup = true) {//z up, for shading
        if (isZup) {
            theta = acosf(v.z);
            phi = atan2f(v.y, v.x);
            r = 1;
        }
        else {
            theta = acosf(v.y);
            phi = atan2f(v.z, v.x);
            r = 1;
        }
    }
    Vec3 convertZup() const {
        return Vec3(r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta));
    }
    Vec3 convertYup() const {
        return Vec3(r * sin(theta) * cos(phi), r * cos(theta), r * sin(theta) * sin(phi));
    }
};


class Quaternions {
public:
    union {//q = d + ai + bj + ck
        struct {
            float a;
            float b;
            float c;
            float d;
        };
        float q[4];
    };
    Quaternions() : a(0), b(0), c(0), d(0) {}
    Quaternions(float _a, float _b, float _c, float _d) : a(_a), b(_b), c(_c), d(_d) {}
    Matrix toMatrix() const {
        float aa = a * a, ab = a * b, ac = a * c;
        float bb = b * b, bc = b * c, cc = c * c;
        float da = d * a, db = d * b, dc = d * c;
        Matrix m;
        m[0] = 1 - 2 * (bb + cc); m[1] = 2 * (ab - dc); m[2] = 2 * (ac + db); m[3] = 0;
        m[4] = 2 * (ab + dc); m[5] = 1 - 2 * (aa + cc); m[6] = 2 * (bc - da); m[7] = 0;
        m[8] = 2 * (ac - db); m[9] = 2 * (bc + da); m[10] = 1 - 2 * (aa + bb); m[11] = 0;
        m[12] = m[13] = m[14] = 0; m[15] = 1;
        return m;
    }
    float magnitude() const {
        return sqrtf(SQ(a) + SQ(b) + SQ(c) + SQ(d));
    }
    float magnitudeSquare() const {
        return SQ(a) + SQ(b) + SQ(c) + SQ(d);
    }
    void normalize() {
        float m = magnitude();
        *this = *this / m;
    }
    Quaternions conjugate() const {
        return Quaternions(-a, -b, -c, d);
    }
    Quaternions inverse() const {
        float m2 = magnitudeSquare();
        return Quaternions(-a / m2, -b / m2, -c / m2, d / m2);
    }
    Quaternions slerp(const Quaternions& q, float t) const {
        Quaternions res;
        if (t < 0 || t > 1) {
            std::cout << "slerp error, invalid t!!" << std::endl;
            return res;
        }
        float theta = findAngle(q);
        std::cout << cosf(theta) << std::endl;
        float q1_co = sinf(theta * (1 - t)) / sinf(theta);
        float q2_co = sinf(theta * t) / sinf(theta);
        res.a = q1_co * a + q2_co * q.a;
        res.b = q1_co * b + q2_co * q.b;
        res.c = q1_co * c + q2_co * q.c;
        res.d = q1_co * d + q2_co * q.d;
        return res;
    }
    Vec3 pointRotate(const Vec3& v) {
        Quaternions p(v.v[0], v.v[1], v.v[2], 0);
        p = *this * p * conjugate();
        return Vec3(p.a, p.b, p.c);
    }
    float dot(const Quaternions& q) const {
        return a * q.a + b * q.b + c * q.c + d * q.d;
    }
    float findAngle(const Quaternions& q) const {
        return acosf(dot(q));
    }
    Quaternions operator+(const Quaternions& q) const {
        return Quaternions(a + q.a, b + q.b, c + q.c, d + q.d);
    }
    Quaternions operator+(const float val) const {
        return Quaternions(a + val, b + val, c + val, d + val);
    }
    Quaternions operator-(const Quaternions& q) const {
        return Quaternions(a - q.a, b - q.b, c - q.c, d - q.d);
    }
    Quaternions operator-(const float val) const {
        return Quaternions(a - val, b - val, c - val, d - val);
    }
    Quaternions operator-() const {
        return Quaternions(-a, -b, -c, -d);
    }
    Quaternions operator*(const Quaternions& q) const {
        return Quaternions(d * q.a + a * q.d + b * q.c - c * q.b,
            d * q.b - a * q.c + b * q.d + c * q.a,
            d * q.c + a * q.b - b * q.a + c * q.d,
            d * q.d - a * q.a - b * q.b - c * q.c);
    }
    Quaternions operator/(const float val) const {
        return Quaternions(a / val, b / val, c / val, d / val);
    }
    void print() const {
        std::cout << '(' << a << ',' << b << ',' << c << ',' << d << ')' << std::endl;
    }
};
class Colour {
public:
    float r, g, b, a;
    Colour(float _r, float _g, float _b, float _a) : r(_r), g(_g), b(_b), a(_a) {}
    Colour(unsigned char _r, unsigned char _g, unsigned char _b, unsigned char _a)
        : r(_r / 255.0f), g(_g / 255.0f), b(_b / 255.0f), a(_a / 255.0f) {
    }
    Colour operator+(const Colour& c) const {
        return Colour(c.r + r, c.g + g, c.b + b, c.a + a);
    }
    Colour operator*(const Colour& c) const {
        return Colour(c.r * r, c.g * g, c.b * b, c.a * a);
    }
    Colour operator*(const float _a) const {
        return Colour(r * _a, g * _a, b * _a, a * _a);
    }
    Colour operator/(const float _a) const {
        return Colour(r / _a, g / _a, b / _a, a / _a);
    }

};

class ShadingFrame {
public:
    union {
        Vec3 v[3];
        struct { Vec3 i, j, k; };
        float m[3][3];
    };
    ShadingFrame(Vec3 _v) {
        _v = _v.normalize();
        //create vectors linearly independent
        Vec3 v1, v2;
        do {
            v1 = randomVector();
            v2 = randomVector();
        } while (fabs(determinant3x3(_v, v1, v2)) < 1e-9);
        //Gram-Schmidt
        v[0] = _v;
        v[1] = v1 - (v[0] * (v1.dot(v[0]) / v[0].dot(v[0])));
        v[2] = v2 - (v[0] * (v2.dot(v[0]) / v[0].dot(v[0])));
        v[2] = v[2] - (v[1] * (v[2].dot(v[1]) / v[1].dot(v[1])));
        //normlize        
        v[1] = v[1].normalize();
        v[2] = v[2].normalize();
    }
    float determinant3x3(Vec3& v0, Vec3& v1, Vec3& v2) {
        return v0[0] * (v1[1] * v2[2] - v1[2] * v2[1])
            - v0[1] * (v1[0] * v2[2] - v1[2] * v2[0])
            + v0[2] * (v1[0] * v2[1] - v1[1] * v2[0]);
    }
    Vec3 randomVector() {
        return Vec3(rand() % 10 + 1, rand() % 10 + 1, rand() % 10 + 1);
    }
    void print() const {
        std::cout << "orthonormal basis:" << std::endl;
        v[0].print();
        v[1].print();
        v[2].print();
    }
};
//above undone
int main() {
    srand((unsigned)time(NULL));
    ShadingFrame s(Vec3(1, 2, 3));
    s.print();
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            std::cout << s.m[i][j] << '\t';
        }
        std::cout << std::endl;
    }
    /*float a[4][4];
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            a[i][j] = rand() % 10;
    Matrix m(a[0]);
    m.print();
    std::cout << "----------" << std::endl;
    m.transpose().print();
    float a[4][4];
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            a[i][j] = rand() % 10;
    Matrix m(a[0]);
    m.print();

    ShadingFrame s(Vec3(1,2,3));
    s.print();

    Quaternions a(1, 2, 3, 4);
    Quaternions b(2, 3, 4, 5);
    a.normalize();
    b.normalize();
    a.print();
    b.print();
    std::cout << "----------" << std::endl;
    Quaternions c = a.slerp(b, 0.5);
    c.print();
    std::cout << "----------" << std::endl;
    Vec3 v(4, 5, 6);
    a.pointRotate(v).print();*/

    return 0;
}