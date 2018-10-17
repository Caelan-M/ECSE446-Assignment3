/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once
#include <random>

TR_NAMESPACE_BEGIN

/**
 * Reflection occlusion integrator
 */
struct ROIntegrator : Integrator {

    float m_exponent;

    explicit ROIntegrator(const Scene& scene) : Integrator(scene) {
        m_exponent = scene.config.integratorSettings.ro.exponent;
    }

    inline v3f reflect(const v3f& d) const {
        return v3f(-d.x, -d.y, d.z);
    }

    v3f render(const Ray& ray, Sampler& sampler) const override {
        v3f Li(0.f);

        SurfaceInteraction surfInt;

        if(scene.bvh -> intersect(ray, surfInt)){

            v3f sampleDir = Warp::squareToPhongLobe(sampler.next2D(), m_exponent);
            float pdf = Warp::squareToPhongLobePdf(sampleDir, m_exponent);

            v3f reflectDir = reflect(surfInt.wo);

            surfInt.wi = glm::toMat4(glm::quat(v3f(0.0f, 0.0f, 1.0f), reflectDir)) * v4f(sampleDir, 1.f);
            surfInt.wi = glm::normalize(surfInt.wi);

            surfInt.wi = surfInt.frameNs.toWorld(surfInt.wi);
            float cosFact = dot(surfInt.wi, surfInt.frameNs.n);

            SurfaceInteraction i;
            Ray sampleRay(surfInt.p, surfInt.wi);
            if(!scene.bvh -> intersect(sampleRay, i))
                Li =  (m_exponent + 2) * INV_TWOPI * v3f(1.f) * std::max(0.f, cosFact) * std::pow(std::max(sampleDir.z, 0.f), m_exponent) / (pdf);
        }
        return Li;
    }
};

TR_NAMESPACE_END