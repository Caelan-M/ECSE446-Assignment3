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
        v3f position = scene.getFirstLightPosition();
        v3f intensity = scene.getFirstLightIntensity();

        if(scene.bvh->intersect(ray, surfInt)) {
            //float cosFact = dot(surfInt.wi, surfInt.frameNs.n);
            glm::vec3 sampleDir = glm::normalize(Warp::squareToCosineHemisphere(sampler.next2D()));
            //glm::vec3 sampleDir = glm::normalize(Warp::squareToUniformHemisphere(sampler.next2D()));
            //glm::vec3 sampleDir = glm::normalize(Warp::squareToUniformSphere(sampler.next2D()));
            //glm::vec3 sampleDir = glm::normalize(Warp::squareToPhongLobe(sampler.next2D()));
            //Ray sampleRay(surfInt.p, sampleDir, Epsilon, scene.aabb.getBSphere().radius - Epsilon);

            //glm::vec3 sampleDir = Warp::squareToUniformSphere(sampler.next2D());
            //glm::vec3 worldDir = surfInt.frameNs.toWorld(sampleDir);
            surfInt.wi = surfInt.frameNs.toWorld(sampleDir);

            Ray sampleRay(surfInt.p, surfInt.wi, Epsilon, scene.aabb.getBSphere().radius * 0.5);
            if(!scene.bvh->intersect(sampleRay, i)) {
                //float dist2hit = glm::distance(surfInt.p, i.p);
                //if(dist2hit < scene.aabb.getBSphere().radius)
                //if(cosFact > 0)
                //Li = v3f(1.f);
                //float cosFact = dot(surfInt.frameNs.n, surfInt.wi);
                //surfInt.wi = glm::normalize(position - surfInt.p);
                //surfInt.wi = surfInt.frameNs.toLocal(surfInt.wi);

                float cosFact = surfInt.frameNs.cosTheta(surfInt.wi);
                float cosFact1 = dot(surfInt.frameNs.n, surfInt.wi);
                float cosFact2 = dot(surfInt.frameNs.toWorld(surfInt.frameNs.n), surfInt.wi);
                float cosFact3 = dot(surfInt.frameNs.n, sampleDir);
                float cosFact4 = dot(surfInt.frameNs.toWorld(surfInt.frameNs.n), surfInt.wi);
                float cosFact5 = dot(surfInt.frameNs.n, surfInt.wi);

                //float correctSide = dot(worldDir, surfInt.frameNs.n); //this check is for sphere to make sure it polls on the corrct side of the surface
                //if(cosFact4 > 0)
                //Li = v3f(1.f) / Warp::squareToUniformSpherePdf() / M_PI * std::max(0.f, cosFact5);//* abs(cosFact1);//*cosFact;
                //Li = v3f(1.f) / Warp::squareToUniformHemispherePdf(sampleDir) / M_PI * std::max(0.f, cosFact5);
                Li = v3f(1.f) / Warp::squareToCosineHemispherePdf(sampleDir) / M_PI * std::max(0.f, cosFact5);
                //Li = v3f(1.f) / Warp::squareToPhongLobePdf(sampleDir) / M_PI * std::max(0.f, cosFact5);

                //Li = v3f(1.f)/Warp::squareToUniformHemispherePdf();
                //Li = v3f(1.f)/Warp::squareToCosineHemispherePdf(surfInt.wi);
                //Li = v3f(1.f)/Warp::squareToPhongLobePdf(surfInt.wi);
            }
        }

//        Your next task is to implement three (3) Monte Carlo estimators for ambient occlusion according the following
//        importance sampling schemes: uniform spherical direction sampling, uniform hemispherical direction sampling,
//        and cosine-weighted hemispherical direction sampling. You can find the API you need to implement for the
//        AOIntegrator in src/integrators/ao.h, and the algorithm should roughly follow the structure below:
//
//        1. Intersect your eye rays with the scene geometry.

//        2. If there's an intersection i, solve for the appropriate Monte Carlo AO estimate at this shade point:
//        you can sample directions in your MC estimator using the sampling routines you developed earlier.

//        3. When evaluating the visibility in the AO integrand of your MC estimator, take care when specifying
//        relevant positions and directions; remember, all intersection routines expect world-space coordinates.
//        Here, you will need to compute the parameters of a shadow ray based on i.p and i.wi.

//        4. When computing the contribution of each MC sample to the final integral estimate, don't forget to
//        evaluate all the integrand terms and to divide by the appropriate PDF evaluation.
//
//        Note that you will have to limit the length of the shadow rays using a scene bounding sphere heuristic.
//        Think about what would happen if you rendered an indoor scene with AO: any shadow ray you trace would
//        hit the scene and the final image would be completely black. To avoid this problem, set the maximum
//        shadow ray length to half of the bounding sphere radius. Use the scene.aabb.getBSphere() to retrieve
//        this sphere and its corresponding radius.
//
//        If you did not implement multiple pixel samples in A1, you can control the MC sampling rate as a
//        parameter of the integrator, looping and averaging over your samples directly in the integrator.
//        In this case, you will shade rays passing through the pixel centers, but you will sample multiple
//        directions per pixel in AOIntegrator (and ROIntegrator).

        return Li;
    }
};

TR_NAMESPACE_END