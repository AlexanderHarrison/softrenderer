

typedef struct Tex {
    U8 *pixels;
    USize width, height, pitch;
} Tex;

typedef union RGBA {
    struct {
        U8 r, g, b, a;
    };
    U32 u32;
} RGBA;

typedef struct Vtx {
    F32 x, y, z;
} Vtx;

typedef __m256 F32_8;
typedef __m256i U32_8;

typedef struct Vtx_8 {
    F32_8 x, y, z;
} Vtx_8;

static inline void SetPixel(Tex tex, USize x, USize y, RGBA colour) {
    *(U32*)(tex.pixels + y*tex.pitch + x*4) = bswap_32(colour.u32);
}

static inline F32_8 CalcArea_8(F32_8 x0, F32_8 y0, F32_8 x1, F32_8 y1, F32_8 x2, F32_8 y2) {
    F32_8 x01 = _mm256_sub_ps(x0, x1);
    F32_8 y20 = _mm256_sub_ps(y2, y0);
    F32_8 x02 = _mm256_sub_ps(x0, x2);
    F32_8 y10 = _mm256_sub_ps(y1, y0);
    
    F32_8 x01_y20 = _mm256_mul_ps(x01, y20);
    F32_8 x02_y10 = _mm256_mul_ps(x02, y10);
    
    return _mm256_sub_ps(x01_y20, x02_y10);
}

static inline Vtx_8 CalcBarycentric_8(F32_8 x, F32_8 y, Vtx *v1, Vtx *v2, Vtx *v3, F32 area_recip) {
    F32_8 v1_x = _mm256_set1_ps(v1->x);
    F32_8 v1_y = _mm256_set1_ps(v1->y);
    F32_8 v2_x = _mm256_set1_ps(v2->x);
    F32_8 v2_y = _mm256_set1_ps(v2->y);
    F32_8 v3_x = _mm256_set1_ps(v3->x);
    F32_8 v3_y = _mm256_set1_ps(v3->y);
    // return (Vtx_8) { v1_x, v1_y, v2_x };

    F32_8 a23 = CalcArea_8(x, y, v2_x, v2_y, v3_x, v3_y);
    F32_8 a31 = CalcArea_8(v1_x, v1_y, x, y, v3_x, v3_y);
    F32_8 a12 = CalcArea_8(v1_x, v1_y, v2_x, v2_y, x, y);
    
    F32_8 area_recip_8 = _mm256_set1_ps(area_recip);
    a23 = _mm256_mul_ps(a23, area_recip_8);
    a31 = _mm256_mul_ps(a31, area_recip_8);
    a12 = _mm256_mul_ps(a12, area_recip_8);
    return (Vtx_8) { a23, a31, a12 };
}

static inline F32_8 TriInterpolate_8(Vtx_8 bary, F32 f1, F32 f2, F32 f3) {
    F32_8 a1 = _mm256_mul_ps(bary.x, _mm256_set1_ps(f1));
    F32_8 a2 = _mm256_mul_ps(bary.y, _mm256_set1_ps(f2));
    F32_8 a3 = _mm256_mul_ps(bary.z, _mm256_set1_ps(f3));
    F32_8 sum = _mm256_add_ps(a1, a2);
    sum = _mm256_add_ps(sum, a3);
    return sum;
}

// returned U32 values will be 0 if inside, 1 if outside.
static inline U32_8 OutsideTri_8(Vtx_8 bary) {
    F32_8 o = _mm256_or_ps(bary.x, bary.y);
    o = _mm256_or_ps(o, bary.z);
    return (U32_8)o;
}

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

static F32 *depth;
USize depth_width = 0;
USize depth_height = 0;

__attribute__((noinline))
static void RenderTri_8(Tex tex, Vtx *vtx, U16 *idx, Vtx *vtx_colour) {
    F32 w = (F32)tex.width;
    F32 h = (F32)tex.height;
    F32 wr = 1.f / w;
    F32 hr = 1.f / h;
    
    USize v1_i = (USize)idx[0];
    USize v2_i = (USize)idx[1];
    USize v3_i = (USize)idx[2];
    
    Vtx *v1 = vtx + v1_i;
    Vtx *v2 = vtx + v2_i;
    Vtx *v3 = vtx + v3_i;
    
    F32 area = (v1->x - v2->x)*(v3->y - v1->y) - (v1->x - v3->x)*(v2->y - v1->y);
    if (-1e-5f < area && area < 1e-5f)
        return;
    
    if (area < 0.f)
        return;
    F32 area_recip = 1.f / area;
    
    F32 x_min = 1.f;
    F32 y_min = 1.f;
    F32 x_max = -1.f;
    F32 y_max = -1.f;
    
    if (v1->x < x_min) x_min = v1->x;
    if (v1->x > x_max) x_max = v1->x;
    if (v1->y < y_min) y_min = v1->y;
    if (v1->y > y_max) y_max = v1->y;
    
    if (v2->x < x_min) x_min = v2->x;
    if (v2->x > x_max) x_max = v2->x;
    if (v2->y < y_min) y_min = v2->y;
    if (v2->y > y_max) y_max = v2->y;
    
    if (v3->x < x_min) x_min = v3->x;
    if (v3->x > x_max) x_max = v3->x;
    if (v3->y < y_min) y_min = v3->y;
    if (v3->y > y_max) y_max = v3->y;
    
    F32 x_inc = 2.f * wr;
    F32 y_inc = 2.f * hr;
    
    ISize x_i_start = (ISize)((x_min + 1.f) * (w * 0.5f));
    ISize y_i_start = (ISize)((y_min + 1.f) * (h * 0.5f));
    ISize x_i_end = (ISize)((x_max + 1.f) * (w * 0.5f));
    ISize y_i_end = (ISize)((y_max + 1.f) * (h * 0.5f));
    
    if (x_i_start < 0) x_i_start = 0;
    if (y_i_start < 0) y_i_start = 0;
    if (x_i_end > (ISize)tex.width) x_i_end = (ISize)tex.width;
    if (y_i_end > (ISize)tex.height) y_i_end = (ISize)tex.height;
    
    F32 x_start_val = wr - 1.f;
    F32 y_start_val = hr - 1.f;
    
    F32 vc[12] = {
        vtx_colour[v1_i].x, vtx_colour[v2_i].x, vtx_colour[v3_i].x, // r
        vtx_colour[v1_i].y, vtx_colour[v2_i].y, vtx_colour[v3_i].y, // g
        vtx_colour[v1_i].z, vtx_colour[v2_i].z, vtx_colour[v3_i].z, // b
    };
    for (int i = 0; i < 9; ++i)
        vc[i] *= 255.f;
        
        
    // TODO: add viewport oob checks
    // TODO: reduce oob calculations somehow

    for (ISize y_i = y_i_start; y_i < y_i_end; y_i += 1) {
        // TODO: oob checks
        for (ISize x_i = x_i_start; x_i < x_i_end; x_i += 8) {
            F32_8 x = _mm256_set1_ps((F32)x_i);
            F32_8 y = _mm256_set1_ps((F32)y_i);
            F32_8 ints = _mm256_set_ps(7, 6, 5, 4, 3, 2, 1, 0);
            x = _mm256_add_ps(x, ints);
            
            F32_8 x_incs = _mm256_set1_ps(x_inc);
            F32_8 y_incs = _mm256_set1_ps(y_inc);
            x = _mm256_mul_ps(x, x_incs);
            y = _mm256_mul_ps(y, y_incs);
            
            F32_8 x_starts = _mm256_set1_ps(x_start_val);
            F32_8 y_starts = _mm256_set1_ps(y_start_val);
            x = _mm256_add_ps(x, x_starts);
            y = _mm256_add_ps(y, y_starts);
            
            Vtx_8 bary = CalcBarycentric_8(x, y, v1, v2, v3, area_recip);
            U32_8 tri_outside = OutsideTri_8(bary);
            
            int mask = _mm256_movemask_ps((F32_8)tri_outside);
            if (mask == 0xff)
                continue;
            
            F32 *depth_row = &depth[(USize)y_i*depth_width + (USize)x_i];
            F32_8 cur_depth = _mm256_loadu_ps(depth_row); // TODO: masked load based on oob vector
            
            F32_8 z = TriInterpolate_8(bary, v1->z, v2->z, v3->z);
            U32_8 z_good = (U32_8)_mm256_cmp_ps(z, cur_depth, _CMP_GE_OQ);
            U32_8 write_pixel = (U32_8)_mm256_andnot_ps((F32_8)tri_outside, (F32_8)z_good);
            
            mask = _mm256_movemask_ps((F32_8)write_pixel);
            if (mask == 0)
                continue;
                
            _mm256_maskstore_ps(depth_row, write_pixel, z);
            
            // TODO: oob and here
            
            F32_8 r = TriInterpolate_8(bary, vc[0], vc[1], vc[2]);
            F32_8 g = TriInterpolate_8(bary, vc[3], vc[4], vc[5]);
            F32_8 b = TriInterpolate_8(bary, vc[6], vc[7], vc[8]);
            
            r = (F32_8)_mm256_cvtps_epi32(r);
            g = (F32_8)_mm256_cvtps_epi32(g);
            b = (F32_8)_mm256_cvtps_epi32(b);
            F32_8 a = (F32_8)_mm256_set1_epi32(255);
            
            // rg     = [0r0r 0r0r 0g0g 0g0g 0r0r 0r0r 0g0g 0g0g]
            // ba     = [0b0b 0b0b 0a0a 0a0a 0b0b 0b0b 0a0a 0a0a]
            // rgba   = [rrrr gggg bbbb aaaa rrrr gggg bbbb aaaa]
            // colour = [rgba rgba rgba rgba rgba rgba rgba rgba]
            U32_8 rg = _mm256_packus_epi32((U32_8)r, (U32_8)g);
            U32_8 ba = _mm256_packus_epi32((U32_8)b, (U32_8)a);
            U32_8 rgba = _mm256_packus_epi16(rg, ba);
            
            #define R(A) A+0, A+4, A+8, A+12
            U32_8 shuffle = _mm256_set_epi8(R(19), R(18), R(17), R(16), R(3), R(2), R(1), R(0));
            F32_8 colour = (F32_8)_mm256_shuffle_epi8(rgba, shuffle);
            
            F32 *target = (F32*)(tex.pixels + (USize)y_i*tex.pitch + (USize)x_i*4);
            _mm256_maskstore_ps(target, write_pixel, colour);
        }
    }
}

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

void SoftRender(Tex tex) {
    static F32 i = 100;
    i += 1;
    
    if (tex.width != depth_width || tex.height != depth_height) {
        depth_width = tex.width;
        depth_height = tex.height;
        USize new_size = tex.width * tex.height * sizeof(*depth);
        if (depth)
            depth = realloc(depth, new_size);
        else
            depth = malloc(new_size);
    }
    memset(depth, 0, depth_width*depth_height*sizeof(*depth));
    
    for (USize y = 0; y < tex.height; ++y) {
        for (USize x = 0; x < tex.width; ++x) {
            SetPixel(tex, x, y, (RGBA) {{ 0, 0, 0, 255 }});
        }
    }
    
    F32 t = i * 3.141592654f * 2.f / 1000.f;
    F32 cos_t = cosf(t);
    F32 sin_t = sinf(t);
    
    static F64 sum = 0.0;
    static F64 count = 0.0;
    
    for (F32 x_t = -1.0; x_t <= 1.0; x_t += 0.3f) {
        for (F32 y_t = -1.0; y_t <= 1.0; y_t += 0.3f) {
            for (USize v_i = 0; v_i < countof(cube_verts); ++v_i) {
                Vtx v = cube_verts[v_i];
                
                // scale
                v.x *= 0.1f;
                v.y *= 0.1f;
                v.z *= 0.1f;
                
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
                RenderTri_8(tex, cube_verts_tform, cube_indices+idx_i, cube_colours);
            }
            sum += timer_elapsed_us();
        }
    }
    count++;
    printf("RenderTris: %f\n", sum / count);
    
    // timer_elapsed("RenderTris");
}
