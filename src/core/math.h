/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Computes barycentric coordinates.
 */
template<class T>
inline T barycentric(const T& a, const T& b, const T& c, const float u, const float v) {
    return a * (1 - u - v) + b * u + c * v;
}

/**
 * Restricts a value to a given interval.
 */
template<class T>
inline T clamp(T v, T min, T max) {
    return std::min(std::max(v, min), max);
}

/**
 * Checks if vector is zero.
 */
inline bool isZero(const v3f v) {
    return glm::dot(v, v) < Epsilon;
}

/**
 * Generates coordinate system.
 */
inline void coordinateSystem(const v3f& a, v3f& b, v3f& c) {
    if (std::abs(a.x) > std::abs(a.y)) {
        float invLen = 1.f / std::sqrt(a.x * a.x + a.z * a.z);
        c = v3f(a.z * invLen, 0.f, -a.x * invLen);
    } else {
        float invLen = 1.f / std::sqrt(a.y * a.y + a.z * a.z);
        c = v3f(0.f, a.z * invLen, -a.y * invLen);
    }
    b = glm::cross(c, a);
}

/**
 * Converts RGB value to luminance.
 */
inline float getLuminance(const v3f& rgb) {
    return glm::dot(rgb, v3f(0.212671f, 0.715160f, 0.072169f));
}

/**
 * Pseudo-random sampler (Mersenne Twister 19937) structure.
 */
struct Sampler {
    std::mt19937 g;
    std::uniform_real_distribution<float> d;
    explicit Sampler(int seed) {
        g = std::mt19937(seed);
        d = std::uniform_real_distribution<float>(0.f, 1.f);
    }
    float next() { return d(g); }
    p2f next2D() { return {d(g), d(g)}; }
    void setSeed(int seed) {
        g.seed(seed);
        d.reset();
    }
};

/**
 * 1D discrete distribution.
 */
struct Distribution1D {
    std::vector<float> cdf{0};
    bool isNormalized = false;

    inline void add(float pdfVal) {
        cdf.push_back(cdf.back() + pdfVal);
    }

    size_t size() {
        return cdf.size() - 1;
    }

    float normalize() {
        float sum = cdf.back();
        for (float& v : cdf) {
            v /= sum;
        }
        isNormalized = true;
        return sum;
    }

    inline float pdf(size_t i) const {
        assert(isNormalized);
        return cdf[i + 1] - cdf[i];
    }

    int sample(float sample) const {
        assert(isNormalized);
        const auto it = std::upper_bound(cdf.begin(), cdf.end(), sample);
        return clamp(int(distance(cdf.begin(), it)) - 1, 0, int(cdf.size()) - 2);
    }
};


/**
 * Warping functions.
 */
namespace Warp {


inline v3f squareToUniformSphere(const p2f& sample) {
    float z = 1 - (2 * sample.x);
    float r = std::sqrt(std::max(0.f, 1.f - (z * z)));
    float phi = 2 * M_PI * sample.y;

    v3f v(r * cos(phi), r * sin(phi), z);

    return v;
}

inline float squareToUniformSpherePdf() {
    float pdf = INV_FOURPI;

    return pdf;
}

inline v3f squareToUniformHemisphere(const p2f& sample) {
	float z = sample.x;
	float r = std::sqrt(std::max(0.f, 1.f - (z * z)));
    float phi = 2 * M_PI * sample.y;

	v3f v(r * cos(phi), r * sin(phi), z);

    return v;
}

inline float squareToUniformHemispherePdf(const v3f& v) {
    float pdf = INV_TWOPI;

    return pdf;
}

inline v2f squareToUniformDiskConcentric(const p2f& sample) {
    float r = std::sqrt(sample.x);
    float theta = 2 * M_PI * sample.y;

    v2f v(r * cos(theta), r * sin(theta));
    return v;
}

inline v3f squareToCosineHemisphere(const p2f& sample) {
    p2f d = squareToUniformDiskConcentric(sample);
    float z = std::sqrt(std::max(0.f, 1.f - d.x * d.x - d.y * d.y));

    v3f v(d.x, d.y, z);
    return v;
}

inline float squareToCosineHemispherePdf(const v3f& v) {
    float cosTheta = v.z; //TODO: look into if this is correct or not
    float pdf = cosTheta * INV_PI;
    return pdf;
}

inline v3f squareToPhongLobe(const p2f& sample, float n) {

    //p2f d = squareToUniformDiskConcentric(sample);
    //float z = std::sqrt(std::max(0.f, 1.f - d.x * d.x - d.y * d.y));

    //v3f v(d.x, d.y, z);

    float alpha = std::acos(std::pow(sample.x, 1.f/(n+1.f)));
    float phi = 2.f*M_PI*sample.y;
    v3f v(std::sin(alpha)*std::cos(phi), std::sin(alpha)*std::sin(phi), std::cos(alpha));

    //p2f d = squareToUniformDiskConcentric(sample);
    //float sintheta = sqrt(-d.x * (d.x - 2));
    //float phi = 2 * M_PI * d.y;
    //float z = std::sqrt(std::max(0.f, 1.f - d.x * d.x - d.y * d.y)); //TODO: compare with others to make sure correct

    //v3f v(sintheta * cos(phi), sintheta * sin(phi), z);  //could be that z = cos phi instead
    return v;
}

inline float squareToPhongLobePdf(const v3f& v, float n) {
    //v3f z(0,0,1);
    //float alpha = glm::dot(v, z);
    //float cosAlpha = cos(alpha);
    float power = std::pow(v.z, n);
    float pdf = (n+2.f)*INV_TWOPI*power;

    //float pdf = acos(pow(sqrt(v.x), v.y + 1));

    //float pdf = INV_TWOPI;

    //float alpha = acos(1/(n+1));
    //float phi = 2*M_PI*sample.y;
    //v3f v(sin(alpha)*cos(phi), sin(alpha)*sin(phi), cos(alpha));

    return pdf;
}

}

TR_NAMESPACE_END