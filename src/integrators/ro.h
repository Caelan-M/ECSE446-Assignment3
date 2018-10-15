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
	    // TODO: Implement this

        SurfaceInteraction surfInt;
//        SurfaceInteraction i;

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


//        if(scene.bvh->intersect(ray, surfInt)) {
//
//            glm::vec3 sampleDir = Warp::squareToPhongLobe(sampler.next2D(), m_exponent);
//
//            glm::vec3 reflectDir = reflect(surfInt.wo);
//            surfInt.wi = glm::toMat4(glm::quat(v3f(0.f,0.f,1.f), reflectDir)) * v4f(sampleDir, 1.f);
//            surfInt.wi = glm::normalize(surfInt.wi);
//
//            float cosFact4 = dot(surfInt.frameNs.n, surfInt.wi);
//            float cosFact1 = surfInt.wi.z;
//            surfInt.wi = surfInt.frameNs.toWorld(surfInt.wi);
//            float cosFact = surfInt.wi.z;
//            //glm::vec3 reflectDir = reflect(surfInt.wo);
//            //surfInt.wi = glm::normalize(glm::mat4(glm::quat(v3f(0.f,0.f,1.f), reflectDir)) * v4f(surfInt.wi, 1.f));
//
//
//
//            //surfInt.wi = surfInt.frameNs.toWorld(surfInt.wi);
//
//            //Ray sampleRay(surfInt.p, surfInt.wi, Epsilon, scene.aabb.getBSphere().radius * 0.5);
//            Ray sampleRay(surfInt.p, surfInt.wi);
//            if(!scene.bvh->intersect(sampleRay, i)) {
//
//                float cosFact5 = dot(surfInt.frameNs.n, surfInt.wi);
//                float pdf = Warp::squareToPhongLobePdf(sampleDir, m_exponent);
//                float power = std::pow(std::max(0.f, sampleDir.z), m_exponent);
//
//                Li = v3f(1.f) * (m_exponent + 2.f) * INV_TWOPI / pdf * power;// * std::max(0.f, cosFact5);
//
//            }
//        }

        return Li;
    }
};

TR_NAMESPACE_END