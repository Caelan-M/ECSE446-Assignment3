/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Ambient occlusion integrator
 */
struct AOIntegrator : Integrator {
    explicit AOIntegrator(const Scene& scene) : Integrator(scene) { }

    v3f render(const Ray& ray, Sampler& sampler) const override {
        v3f Li(0.f);
	    // TODO: Implement this


        SurfaceInteraction surfInt;
        SurfaceInteraction i;

        if(scene.bvh->intersect(ray, surfInt)) {
            //glm::vec3 sampleDir = glm::normalize(Warp::squareToCosineHemisphere(sampler.next2D()));
            //glm::vec3 sampleDir = glm::normalize(Warp::squareToUniformHemisphere(sampler.next2D()));
            glm::vec3 sampleDir = glm::normalize(Warp::squareToUniformSphere(sampler.next2D()));

            surfInt.wi = surfInt.frameNs.toWorld(sampleDir);

            Ray sampleRay(surfInt.p, surfInt.wi, Epsilon, scene.aabb.getBSphere().radius * 0.5);
            if(!scene.bvh->intersect(sampleRay, i)) {
                float cosFact = dot(surfInt.frameNs.n, surfInt.wi);

                Li = v3f(1.f) / Warp::squareToUniformSpherePdf() / M_PI * std::max(0.f, cosFact);
                //Li = v3f(1.f) / Warp::squareToUniformHemispherePdf(sampleDir) / M_PI * std::max(0.f, cosFact);
                //Li = v3f(1.f) / Warp::squareToCosineHemispherePdf(sampleDir) / M_PI * std::max(0.f, cosFact);

            }
        }
        return Li;
    }
};

TR_NAMESPACE_END