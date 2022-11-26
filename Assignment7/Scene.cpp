//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"


void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum){
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}

// Implementation of Path Tracing
Vector3f Scene::castRayRecursive(const Ray &ray, int depth) const
{
    // TO DO Implement Path Tracing Algorithm here
    // pseudo-code of Path Tracing Algorithm
    // shade(p, wo)
    //      sampleLight(inter, pdf_light)
    //      Get x, ws, NN, emit from inter
    //      Shoot a ray from p to x
    //      If the ray is not blocked in the middle
    //          L_dir = emit * eval(wo, ws, N) * dot(ws, N) * dot(ws, NN) / |x-p|^2 / pdf_light
    //      L_indir = 0.0
    //      Test Russian Roulette with probability RussianRoulette
    //      wi = sample(wo, N)
    //      Trace a ray r(p, wi)
    //      If ray r hit a non-emitting object at q
    //          L_indir = shade(q, wi) * eval(wo, wi, N) * dot(wi, N) / pdf(wo, wi, N) / RussianRoulette
    //      Return L_dir + L_indir
    float epsilon = 0.0005;

    Intersection hitObject = this->intersect(ray);

    if (!hitObject.happened) return Vector3f(0.0, 0.0, 0.0);

    if (hitObject.m->hasEmission()) return hitObject.m->getEmission();

    // wo: direction from object point to eye position
    Vector3f wo = -ray.direction;
    Vector3f no = hitObject.normal;
    
    Intersection srcLight;
    float pdf_light;
    this->sampleLight(srcLight, pdf_light);

    // ws: direction from light source to object point
    Vector3f objectToLight = srcLight.coords - hitObject.coords;
    Vector3f ws = objectToLight.normalized();
    Vector3f ns = srcLight.normal;

    Intersection blockObject = this->intersect(Ray(hitObject.coords, ws, 0.0));

    Vector3f L_dir(0.0, 0.0, 0.0);
    if (blockObject.distance - objectToLight.norm() > -epsilon) {
        float distanceSquare = std::pow(objectToLight.norm(), 2);
        float objectCos = std::max(0.0f, dotProduct(no, ws));
        float LightCos = std::max(0.0f, dotProduct(ns, -ws));
        L_dir = srcLight.emit 
            * hitObject.m->eval(ws, wo, no) 
                * objectCos * LightCos / distanceSquare / pdf_light;
    }

    Vector3f L_indir(0.0, 0.0, 0.0);
    if (get_random_float() < RussianRoulette) {
        Vector3f wi = normalize(hitObject.m->sample(wo, no));
        
        float pdf_hemi = hitObject.m->pdf(wi, wo, no);
        if (pdf_hemi > epsilon) {
            Ray nextRay(hitObject.coords, wi);
            Intersection nextObject = this->intersect(nextRay);

            if (nextObject.happened && !nextObject.m->hasEmission()) {
                float nextCos = std::max(0.0f, dotProduct(wi, no));
                L_indir = 
                    this->castRay(nextRay, depth + 1)
                        * hitObject.m->eval(wi, wo, no) 
                            * nextCos / pdf_hemi / RussianRoulette;
            }
        }
    }
    return L_dir + L_indir;
}


// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    return castRayRecursive(ray, 0);
}