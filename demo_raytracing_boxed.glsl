#define PI 			3.141592653589793
#define DIFFUSE 1
#define METAL 2
#define GLASS 3
precision highp float;

struct Material {
    int materialType;
    vec3 albedo; // diffuse reflection (r,g,b)
    float fuzz;
    float refractionIndex;
};

struct HitRecord {
    vec3 p;
    float t;
    vec3 normal;
    Material material;
};

struct Sphere {
    vec3 center;
    float radius;
    Material material;
};


struct Ray {
    vec3 origin;
    vec3 direction;
};

struct PointLight {
    vec3 position;
};

vec2 randState;


float rand2D()
{
    randState.x = fract(sin(dot(randState.xy + iTime, vec2(12.9898, 78.233))) * 43758.5453);
    randState.y = fract(sin(dot(randState.xy + iTime, vec2(12.9898, 78.233))) * 43758.5453);;

        return randState.x;
}

/**
 * Approximating the contribution of the Fresnel factor in the specular reflection
 *
 * @see https://en.wikipedia.org/wiki/Schlick%27s_approximation
 *
 * @param cosine the cosine of the angle between the direction from which the 
 * incident light is coming and the normal of the interface.
 * @param refractionIndex refraction index of the interface
 */
float schlick(float cosine, float refractionIndex) {
    float r0 = (1.0 - refractionIndex) / (1.0 + refractionIndex);
    r0 = r0 * r0;
    return r0 + (1.0 - r0) * pow((1.0 - cosine), 5.0);
}


// random direction in unit sphere (for lambert brdf)
vec3 random_in_unit_sphere()
{
    float phi = 2.0 * PI * rand2D();
    float cosTheta = 2.0 * rand2D() - 1.0;
    float u = rand2D();

    float theta = acos(cosTheta);
    float r = pow(u, 1.0 / 3.0);

    float x = r * sin(theta) * cos(phi);
    float y = r * sin(theta) * sin(phi);
    float z = r * cos(theta);

    return vec3(x, y, z);
}


bool intersect(Ray ray, Sphere sphere, float t_min, float t_max, out HitRecord record) {

    // return false;

    vec3 origin = ray.origin;
    vec3 direction = ray.direction;
    vec3 center = sphere.center;
    vec3 oc = origin - center;
    float a = dot(direction, direction);
    float b = dot(direction, oc);
    float c = dot(oc, oc) - sphere.radius * sphere.radius;
    float discriminant = b * b - a * c;

    bool isHit = false;
    if (discriminant >= 0.0) {

        float temp = (-b - sqrt(discriminant)) / a;
        if (temp > t_min && temp < t_max) {
            isHit = true;
            record.t = temp;
            record.p = ray.origin + temp * ray.direction;
            record.normal = normalize(record.p - sphere.center);
            record.material = sphere.material;
            return isHit;
        }

        temp = (-b + sqrt(discriminant)) / a;
        if (temp > t_min && temp < t_max) {
            isHit = true;          
            record.t = temp;
            record.p = ray.origin + temp * ray.direction;
            record.normal = normalize(record.p - sphere.center);
            record.material = sphere.material;
            return isHit;
        }
    }      

    return isHit;
}



Sphere[] spheres = Sphere[10](
    Sphere(vec3(-2.0, 2.0, 1.0), 1.0, Material(GLASS, vec3(0.1, 0.7, 0.1), 0.0, 1.0)),
    Sphere(vec3(-4.0, 1.0, 0.8), 1.0, Material(DIFFUSE, vec3(0.4, 0.2, 0.1), 0.0, 1.0)),
    Sphere(vec3(1.0, 1.0, 0.0), 1.0, Material(METAL, vec3(0.5, 0.5, 0.5), 0.0, 1.0)),
    Sphere(vec3(2.0, 2.8, 0.0), 1.0, Material(METAL, vec3(0.5, 0.5, 0.5), 0.0, 1.0)),
    Sphere(vec3(3.0, 1.0, 0.0), 1.0, Material(METAL, vec3(0.5, 0.5, 0.5), 0.0, 1.0)),
    Sphere(vec3(0.0, -10000.0, 0.0), 10000.0, Material(DIFFUSE, vec3(0.6, 0.6, 0.6), 0.0, 1.0)), // BOTTOM
    Sphere(vec3(0.0, 10010.0, 0.0), 10000.0, Material(DIFFUSE, vec3(0.6, 0.6, 0.6), 0.0, 1.0)), // TOP
    Sphere(vec3(-10006.0, 0.0, 0.0), 10000.0, Material(DIFFUSE, vec3(0.6, 0.5, 1.0), 0.0, 1.0)), // LEFT
    Sphere(vec3(10006.0, 0.0, 0.0), 10000.0, Material(DIFFUSE, vec3(1.0, 0.5, 0.6), 0.0, 1.0)), // RIGHT
    Sphere(vec3(0.0, 0.0, -10005.0), 10000.0, Material(DIFFUSE, vec3(0.9, 0.9, 0.9), 0.0, 1.0)) // BACK
);  

PointLight[] pointLights = PointLight[1](
    PointLight(vec3(0.0, 3.0, 1.0))
);

/*
vec3 computeLight(HitRecord record) {
    
    vec3 I = vec3(0.0);

    for (int i = 0; i < pointLights.length(); i++) {
        vec3 L = pointsLights[i].position - record.p;
        float dist2Light = length(L);
        
        for (int i = 0; i < spheres.length(); i++) {
            
        }
    }
    return vec3(0.0);
}
*/

bool intersectWorld(Ray ray, float t_min, float t_max, out HitRecord record) {
    bool isHit = false;

    for (int i = 0; i < spheres.length(); i++) {
        bool isCurrHit = false;
        HitRecord currentRecord;
        isCurrHit = intersect(ray, spheres[i], t_min, t_max, currentRecord);
        if (isCurrHit) {
            isHit = true;
            t_max = currentRecord.t;
            record = currentRecord;
        }
    }

    return isHit;
}



bool scatter(inout Ray ray, HitRecord record) {
    int materialType = record.material.materialType;
    
    if (materialType == METAL) {
        vec3 reflected = reflect(ray.direction, record.normal);
        float fuzz = record.material.fuzz > 1.0 ? 1.0 : record.material.fuzz < 0.0 ? 0.0 : record.material.fuzz;
        if (fuzz > 0.0) {
            ray = Ray(record.p, normalize(reflected) + fuzz * random_in_unit_sphere());
        } else {
            ray = Ray(record.p, reflected);
        }
        return (dot(ray.direction, record.normal) > 0.0);
        
    } else if (materialType == GLASS) {

        vec3 outward_normal;
        vec3 reflected = reflect(ray.direction, record.normal);
        float ni_over_nt;
        float reflected_prob;
        float cosine;
        float refractionIndex = record.material.refractionIndex;

        if (dot(ray.direction, record.normal) > 0.0) {
            outward_normal = - record.normal;
            ni_over_nt = record.material.refractionIndex;
            cosine = refractionIndex * dot(normalize(ray.direction), record.normal);
        } else {
            outward_normal = record.normal;
            ni_over_nt = 1.0 / record.material.refractionIndex;
            cosine = -refractionIndex * dot(normalize(ray.direction), record.normal);
        }

        vec3 refracted = refract(normalize(ray.direction), outward_normal, ni_over_nt);

        if (refracted.x != 0.0 && refracted.y != 0.0 && refracted.z != 0.0) {
            reflected_prob = schlick(cosine, refractionIndex);
            // ray = Ray(record.p, refracted);
        } else {
            reflected_prob  = 1.0;
            // ray = Ray(record.p, reflected);
        }

        if (rand2D() < reflected_prob) {
            ray = Ray(record.p, reflected);          
        } else {
            ray = Ray(record.p, refracted);
        }

        return true;

    } else {
        ray = Ray(record.p, record.normal + random_in_unit_sphere());
        return true;
    }
}



#define MAX_DEPTH 4

vec3 shootRay(Ray ray, float t_min, float t_max) {
    HitRecord record;

    int hitCounts = 0;
    bool isHit = intersectWorld(ray, t_min, t_max, record);
    vec3 scale = vec3(1.0, 1.0, 1.0);

    while(hitCounts < MAX_DEPTH) {
        hitCounts++;
        bool isScatter = scatter(ray, record);

        if (!isScatter) return vec3(0.0, 0.0, 0.0);

        scale *= record.material.albedo;
        isHit = intersectWorld(ray, t_min, t_max, record);   
    }

    float t = (normalize(ray.direction).y + 1.0) * 0.5;
    vec3 color = (1.0 - t) * vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 1.0);

    return scale * color;
}


struct Camera {
    vec3 center;
    vec3 lower_left;
    vec3 hori;
    vec3 vert;
};

Camera getCamera(vec3 lookfrom, vec3 lookat, vec3 vup, float vfov, float aspect) {
    vec3 u, v, w;
    Camera camera;
    float theta = vfov * PI / 180.0;
    float half_height = tan(theta / 2.0);
    float half_width = aspect * half_height;
    camera.center = lookfrom;
    w = normalize(lookfrom - lookat);
    u = normalize(cross(vup, w));
    v = cross(w, u);
    camera.lower_left = camera.center - half_width * u - half_height * v - w;
    camera.hori = 2.0 * half_width * u;
    camera.vert = 2.0 * half_height * v;
    return camera;
}

Ray getRay(Camera camera, float u, float v) {
    return Ray(camera.center, camera.lower_left + u * camera.hori + v * camera.vert - camera.center);
}


void mainImage( out vec4 frag_color, in vec2 frag_coord )
{
    if (ivec2(frag_coord) == ivec2(0)) {
        frag_color = iResolution.xyxy;
        return;
    }
    float ratio = iResolution.x / iResolution.y;

    Camera closeUpReflection = getCamera(vec3(2.0, 2.0, 2.0), vec3(2.0, 1.5, 0.0), vec3(0.0, 1.0, 0.0), 90.0, ratio);
    Camera wholeScene = getCamera(vec3(3.0, 2.0, 2.0), vec3(2.0, 1.5, 0.0), vec3(0.0, 1.0, 0.0), 90.0, ratio);
    Camera wholeSceneFront = getCamera(vec3(0.0, 2.0, 5.0), vec3(0.0, 2.0, 0.0), vec3(0.0, 1.0, 0.0), 90.0, ratio);
    
    Camera camera = wholeSceneFront;

    vec3 color = vec3(0.0, 0.0, 0.0);

    randState = frag_coord.xy / iResolution.xy;
    
    float u = float(frag_coord.x + rand2D()) / float(iResolution.x);
    float v = float(frag_coord.y + rand2D()) / float(iResolution.y);

    Ray ray = getRay(camera, u, v);

    color = shootRay(ray, 0.01, 1000.0);

    if (texelFetch(iChannel0, ivec2(0),0).xy == iResolution.xy) {        
        frag_color = vec4(color,1) + texelFetch(iChannel0, ivec2(frag_coord), 0);
    } else {        
        frag_color = vec4(color,1);
    }

}