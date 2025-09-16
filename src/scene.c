// #include "softr.c"
// #include "softr-avx.c"
// #include "softr-tight.c"
#include "softr-tight-avx.c"

static Vtx cube_verts[] = {
    { -1, -1, -1 },
    { -1, -1,  1 },
    { -1,  1, -1 },
    
    { -1,  1,  1 },
    {  1, -1, -1 },
    {  1, -1,  1 },
    {  1,  1, -1 },
    {  1,  1,  1 },
};
static Vtx cube_colours[] = {
    { 0.2f, 0.2f, 0.2f },
    { 0, 0, 1 },
    { 0, 1, 0 },
    { 0, 1, 1 },
    { 1, 0, 0 },
    { 1, 0, 1 },
    { 1, 1, 0 },
    { 1, 1, 1 },
};
static Vtx cube_verts_tform[countof(cube_verts)];
static U16 cube_indices[] = {
    0, 2, 1,
    1, 2, 3,
    2, 7, 3,
    3, 7, 1,
    7, 5, 1, // TRI FOR BENCHMARKING
    1, 5, 0,
    5, 4, 0,
    0, 4, 2,
    4, 6, 2,
    2, 6, 7,
    6, 4, 7,
    7, 4, 5,
};

static struct timespec timer;
void timer_start(void) {
    clock_gettime(CLOCK_MONOTONIC, &timer);
}

F64 timer_elapsed_us(void) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    struct timespec diff = now;
    diff.tv_sec -= timer.tv_sec;
    diff.tv_nsec -= timer.tv_nsec;
    F64 time_us = (F64)diff.tv_sec * 1000000.0 + (F64)diff.tv_nsec / 1000.0;
    return time_us;
}
void timer_elapsed(const char *label) {
    printf("%s: %fus\n", label, timer_elapsed_us());
}

static Depth depth;

void RenderScene(Tex tex) {
    static F32 i = 1.f;
    i += 1;
    
    if (tex.width != depth.width || tex.height != depth.height) {
        depth.width = tex.width;
        depth.height = tex.height;
        USize new_size = tex.width * tex.height * sizeof(depth.pixels);
        if (depth.pixels)
            depth.pixels = realloc(depth.pixels, new_size);
        else
            depth.pixels = malloc(new_size);
    }
    memset(depth.pixels, 0, depth.width*depth.height*sizeof(depth.pixels));
    
    for (USize y = 0; y < tex.height; ++y) {
        for (USize x = 0; x < tex.width; ++x) {
            SetPixel(tex, x, y, (RGBA) {{ 0, 0, 0, 255 }});
        }
    }
    
    F32 t = i * 3.141592654f * 2.f / 1000.f;
    F32 cos_t = cosf(t);
    F32 sin_t = sinf(t);
    F32 scale = 0.4f; 
    
    static F64 sum = 0.0;
    static F64 count = 0.0;
    
    F32 x_t = 0.0;
    F32 y_t = 0.0;
    // for (F32 x_t = -1.0; x_t <= 1.0; x_t += 0.3f) {
    //     for (F32 y_t = -1.0; y_t <= 1.0; y_t += 0.3f) {
            for (USize v_i = 0; v_i < countof(cube_verts); ++v_i) {
                Vtx v = cube_verts[v_i];
                
                // scale
                v.x *= scale;
                v.y *= scale;
                v.z *= scale;
                
                // rotate around x axis
                {
                    F32 y = v.y;
                    F32 z = v.z;
                    v.y = y*cos_t - z*sin_t;
                    v.z = y*sin_t + z*cos_t;
                }
                
                // rotate around y axis
                {
                    F32 x = v.x;
                    F32 z = v.z;
                    v.x = x*cos_t - z*sin_t;
                    v.z = x*sin_t + z*cos_t;
                }
                
                // rotate around z axis
                {
                    F32 x = v.x;
                    F32 y = v.y;
                    v.x = x*cos_t - y*sin_t;
                    v.y = x*sin_t + y*cos_t;
                }
                
                // translate
                v.x += x_t;
                v.y += y_t;
                
                // scale z from -1..1 to 0..1
                v.z = v.z * 0.5f + 0.5f;
                
                F32 aspect = (F32)tex.width / (F32)tex.height;
                v.x /= aspect;
                
                cube_verts_tform[v_i] = v;
            }
            
            timer_start();
            for (USize idx_i = 0; idx_i < countof(cube_indices); idx_i += 3) {
                RenderTri(tex, depth, cube_verts_tform, cube_indices+idx_i, cube_colours);
            }
            sum += timer_elapsed_us();
    //     }
    // }
    count++;
    printf("RenderTris: %f\n", sum / count);
    
    // timer_elapsed("RenderTris");
}
