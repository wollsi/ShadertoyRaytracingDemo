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

bool intersect(Ray ray, Sphere s, float t_min, float t_max, out HitRecord rec) {

    vec3 origin = ray.origin;
    vec3 direction = ray.direction;
    vec3 center = s.center;
    vec3 oc = origin - center;
    float a = dot(direction, direction);
    float b = dot(direction, oc);
    float c = dot(oc, oc) - s.radius * s.radius;
    float discriminant = b * b - a * c;

    bool isHit = false;
    if (discriminant > 0.0) {

        float temp = (-b - sqrt(discriminant)) / a;
        if (temp > t_min && temp < t_max) {
            isHit = true;
            rec.t = temp;
            rec.p = ray.origin + temp * ray.direction;
            rec.normal = normalize(rec.p - s.center);
            rec.material = s.material;
            return isHit;
        }

        temp = (-b + sqrt(discriminant)) / a;
        if (temp > t_min && temp < t_max) {
            isHit = true;          
            rec.t = temp;
            rec.p = ray.origin + temp * ray.direction;
            rec.normal = normalize(rec.p - s.center);
            rec.material = s.material;
            return isHit;
        }
    }      

    return isHit;
}


Sphere[] spheres = Sphere[6](
    Sphere(vec3(-2.0, 2.0, 1.0), 1.0, Material(GLASS, vec3(0.1, 0.7, 0.1), 0.0, 1.0)),
    Sphere(vec3(-4.0, 1.0, 0.8), 1.0, Material(DIFFUSE, vec3(0.4, 0.2, 0.1), 0.0, 1.0)),
    Sphere(vec3(1.0, 1.0, 0.0), 1.0, Material(METAL, vec3(0.5, 0.5, 0.5), 0.0, 1.0)),
    Sphere(vec3(2.0, 2.8, 0.0), 1.0, Material(METAL, vec3(0.5, 0.5, 0.5), 0.0, 1.0)),
    Sphere(vec3(3.0, 1.0, 0.0), 1.0, Material(METAL, vec3(0.5, 0.5, 0.5), 0.0, 1.0)),
    Sphere(vec3(0.0, -1000.0, 0.0), 1000.0, Material(DIFFUSE, vec3(0.5, 0.5, 0.5), 0.0, 1.0))
);  

bool intersectWorld(Ray ray, float t_min, float t_max, out HitRecord rec) {
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


bool scatter(inout Ray ray, HitRecord rec) {
    int materialType = rec.material.materialType;
    
    if (materialType == METAL) {
        vec3 reflected = reflect(ray.direction, rec.normal);
        float fuzz = rec.material.fuzz > 1.0 ? 1.0 : rec.material.fuzz < 0.0 ? 0.0 : rec.material.fuzz;
        if (fuzz > 0.0) {
            ray = Ray(rec.p, normalize(reflected) + fuzz * random_in_unit_sphere());
        } else {
            ray = Ray(rec.p, reflected);
        }
        return (dot(ray.direction, rec.normal) > 0.0);
        
    } else if (materialType == GLASS) {

        vec3 outward_normal;
        vec3 reflected = reflect(ray.direction, rec.normal);
        float ni_over_nt;
        float reflected_prob;
        float cosine;
        float refractionIndex = rec.material.refractionIndex;

        if (dot(ray.direction, rec.normal) > 0.0) {
            outward_normal = - rec.normal;
            ni_over_nt = rec.material.refractionIndex;
            cosine = refractionIndex * dot(normalize(ray.direction), rec.normal);
        } else {
            outward_normal = rec.normal;
            ni_over_nt = 1.0 / rec.material.refractionIndex;
            cosine = -refractionIndex * dot(normalize(ray.direction), rec.normal);
        }

        vec3 refracted = refract(normalize(ray.direction), outward_normal, ni_over_nt);

        if (refracted.x != 0.0 && refracted.y != 0.0 && refracted.z != 0.0) {
            reflected_prob = schlick(cosine, refractionIndex);
            // ray = Ray(rec.p, refracted);
        } else {
            reflected_prob  = 1.0;
            // ray = Ray(rec.p, reflected);
        }

        if (rand2D() < reflected_prob) {
            ray = Ray(rec.p, reflected);          
        } else {
            ray = Ray(rec.p, refracted);
        }

        return true;

    } else {
        ray = Ray(rec.p, rec.normal + random_in_unit_sphere());
        return true;
    }
}

#define MAX_DEPTH 4

vec3 shootRay(Ray ray, float t_min, float t_max) {
    HitRecord rec;

    int hitCounts = 0;
    bool isHit = intersectWorld(ray, t_min, t_max, rec);
    vec3 scale = vec3(1.0, 1.0, 1.0);

    while(isHit && hitCounts < MAX_DEPTH) {
        hitCounts++;
        bool isScatter = scatter(ray, rec);

        if (!isScatter) return vec3(0.0, 0.0, 0.0);

        scale *= rec.material.albedo;
        isHit = intersectWorld(ray, t_min, t_max, rec);   
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
    
    Camera camera = wholeScene;

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