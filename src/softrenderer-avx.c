

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

static inline U32 FloatBits(float f) {
    U32 i;
    memcpy(&i, &f, sizeof(U32));
    return i;
}

static inline float CalcArea(F32 x0, F32 y0, F32 x1, F32 y1, F32 x2, F32 y2) {
    return (x0 - x1)*(y2 - y0) - (x0 - x2)*(y1 - y0);
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

static inline Vtx CalcBarycentric(F32 x, F32 y, Vtx *v1, Vtx *v2, Vtx *v3, F32 area_recip) {
    F32 a23 = CalcArea(x, y, v2->x, v2->y, v3->x, v3->y);
    F32 a31 = CalcArea(v1->x, v1->y, x, y, v3->x, v3->y);
    F32 a12 = CalcArea(v1->x, v1->y, v2->x, v2->y, x, y);
    
    Vtx bary = {
        a23 * area_recip,
        a31 * area_recip,
        a12 * area_recip,
    };
    
    return bary;
}

static inline Vtx_8 CalcBarycentric_8(F32_8 x, F32_8 y, Vtx *v1, Vtx *v2, Vtx *v3, F32 area_recip) {
    F32_8 v1_x = _mm256_set1_ps(v1->x);
    F32_8 v1_y = _mm256_set1_ps(v1->y);
    F32_8 v2_x = _mm256_set1_ps(v2->x);
    F32_8 v2_y = _mm256_set1_ps(v2->y);
    F32_8 v3_x = _mm256_set1_ps(v3->x);
    F32_8 v3_y = _mm256_set1_ps(v3->y);

    F32_8 a23 = CalcArea_8(x, y, v2_x, v2_y, v3_x, v3_y);
    F32_8 a31 = CalcArea_8(v1_x, v1_y, x, y, v3_x, v3_y);
    F32_8 a12 = CalcArea_8(v1_x, v1_y, v2_x, v2_y, x, y);
    
    F32_8 area_recip_8 = _mm256_set1_ps(area_recip);
    a23 = _mm256_mul_ps(a23, area_recip_8);
    a31 = _mm256_mul_ps(a31, area_recip_8);
    a12 = _mm256_mul_ps(a12, area_recip_8);
    return (Vtx_8) { a23, a31, a12 };
}

static inline F32 TriInterpolate(Vtx *bary, F32 f1, F32 f2, F32 f3) {
    F32 z1 = f1*bary->x;
    F32 z2 = f2*bary->y;
    F32 z3 = f3*bary->z;
    F32 sum = z1 + z2 + z3;
    return sum;
}

static inline F32_8 TriInterpolate_8(Vtx_8 bary, F32 f1, F32 f2, F32 f3) {
    F32_8 a1 = _mm256_mul_ps(bary.x, _mm256_set1_ps(f1));
    F32_8 a2 = _mm256_mul_ps(bary.y, _mm256_set1_ps(f2));
    F32_8 a3 = _mm256_mul_ps(bary.z, _mm256_set1_ps(f3));
    F32_8 sum = _mm256_add_ps(a1, a2);
    sum = _mm256_add_ps(sum, a3);
    return sum;
}

static inline bool InsideTri(Vtx *bary) {
    U32 b1 = FloatBits(bary->x);
    U32 b2 = FloatBits(bary->y); 
    U32 b3 = FloatBits(bary->z);
    
    // U32 ands = b1 & b2 & b3;
    // return (ands >> 31) == 1;
    U32 ors = b1 | b2 | b3;
    return (ors >> 31) == 0;
    // U32 xors = ands ^ ors;
    // U32 xor_sign = xors >> 31;
    // return xor_sign == 0;
}

// returned U32 values will be 0 if inside, 1 if outside.
static inline U32_8 InsideTri_8(Vtx_8 bary) {
    U32_8 x = (U32_8)bary.x;
    U32_8 y = (U32_8)bary.y;
    U32_8 z = (U32_8)bary.z;
    // U32_8 x = _mm256_srli_epi32((U32_8)bary.x), 31);
    // U32_8 y = _mm256_srli_epi32((U32_8)bary.y, 31);
    // U32_8 z = _mm256_srli_epi32((U32_8)bary.z, 31);
    
    F32_8 o = _mm256_or_ps((F32_8)x, (F32_8)y);
    o = _mm256_or_ps(o, (F32_8)z);
    
    // U32_8 all_ones = _mm256_set1_epi32(-1);
    // o = _mm256_xor_ps(o, (F32_8)all_ones);
    
    return (U32_8)o;
}

static struct timespec timer;
void timer_start(void) {
    clock_gettime(CLOCK_MONOTONIC, &timer);
}
void timer_elapsed(const char *label) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    struct timespec diff = now;
    diff.tv_sec -= timer.tv_sec;
    diff.tv_nsec -= timer.tv_nsec;
    
    double time_us = (double)diff.tv_sec * 1000000.0 + (double)diff.tv_nsec / 1000.0;
    printf("%s: %fus\n", label, time_us);
}

static F32 *depth;
USize depth_width = 0;
USize depth_height = 0;

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
    if (-1e-7f < area && area < 1e-7f)
        return;
    
    if (area > 0.f)
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
    
    F32_8 ints = _mm256_set_ps(0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f); // TODO: check if load every iter is better?
    
    for (ISize y_i = y_i_start; y_i < y_i_end; y_i += 1) {
        // TODO: oob checks
        for (ISize x_i = x_i_start; x_i < x_i_end; x_i += 8) {
            F32_8 x = _mm256_set1_ps((F32)x_i);
            F32_8 y = _mm256_set1_ps((F32)y_i);
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
            U32_8 tri_outside = InsideTri_8(bary);
            // U32_8 tri_outside = _mm256_set1_epi32(0);
            F32_8 z = TriInterpolate_8(bary, v1->z, v2->z, v3->z);
            F32_8 cur_depth = _mm256_loadu_ps(&depth[(USize)y_i*depth_width + (USize)x_i]);
            
            U32_8 z_good = (U32_8)_mm256_cmp_ps(z, cur_depth, _CMP_GE_OQ);
            U32_8 write_pixel = (U32_8)_mm256_andnot_ps((F32_8)tri_outside, (F32_8)z_good);
            
            F32_8 r = TriInterpolate_8(bary, vtx_colour[v1_i].x, vtx_colour[v2_i].x, vtx_colour[v3_i].x);
            F32_8 g = TriInterpolate_8(bary, vtx_colour[v1_i].y, vtx_colour[v2_i].y, vtx_colour[v3_i].y);
            F32_8 b = TriInterpolate_8(bary, vtx_colour[v1_i].z, vtx_colour[v2_i].z, vtx_colour[v3_i].z);
            
            F32_8 f255 = _mm256_set1_ps(255.f);
            r = _mm256_mul_ps(r, f255);
            g = _mm256_mul_ps(g, f255);
            b = _mm256_mul_ps(b, f255);
            
            r = (F32_8)_mm256_cvtps_epi32(r);
            g = (F32_8)_mm256_cvtps_epi32(g);
            b = (F32_8)_mm256_cvtps_epi32(b);
            
            r = (F32_8)_mm256_slli_epi32((U32_8)r, 24);
            g = (F32_8)_mm256_slli_epi32((U32_8)g, 16);
            b = (F32_8)_mm256_slli_epi32((U32_8)b, 8);
            F32_8 a = (F32_8)_mm256_set1_epi32(255);
            
            F32_8 rg = _mm256_or_ps(r, g);
            F32_8 ba = _mm256_or_ps(b, a);
            F32_8 colour = _mm256_or_ps(rg, ba);
            
            F32 *target = (F32*)(tex.pixels + (USize)y_i*tex.pitch + (USize)x_i*4);
            _mm256_maskstore_ps(target, write_pixel, colour);
            
                // U32_8 x = _mm256_srli_epi32((U32_8)bary.x), 31);
            
                // F32 *pixel_depth = &depth[(USize)y_i*depth_width + (USize)x_i];
                // if (z >= *pixel_depth) {
                //     *pixel_depth = z;
                    
                //     F32 r = TriInterpolate(&bary, vtx_colour[v1_i].x, vtx_colour[v2_i].x, vtx_colour[v3_i].x);
                //     F32 g = TriInterpolate(&bary, vtx_colour[v1_i].y, vtx_colour[v2_i].y, vtx_colour[v3_i].y);
                //     F32 b = TriInterpolate(&bary, vtx_colour[v1_i].z, vtx_colour[v2_i].z, vtx_colour[v3_i].z);
                    
                //     F32 f = 255.f;
                //     RGBA colour = {{
                //         (U8)(r * f),
                //         (U8)(g * f),
                //         (U8)(b * f),
                //         255.f,
                //     }};
                    
                //     SetPixel(tex, (USize)x_i, (USize)y_i, colour);
                // }
            // }
        }
    }
}

static void RenderTri(Tex tex, Vtx *vtx, U16 *idx, Vtx *vtx_colour) {
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
    if (-1e-7f < area && area < 1e-7f)
        return;
    
    if (area > 0.f)
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
    
    F32 x_start = wr - 1.f;
    F32 y_start = hr - 1.f;
    
    for (ISize y_i = y_i_start; y_i < y_i_end; y_i += 1) {
        for (ISize x_i = x_i_start; x_i < x_i_end; x_i += 1) {
            F32 x = x_start + x_inc*(F32)x_i;
            F32 y = y_start + y_inc*(F32)y_i;
            
            Vtx bary = CalcBarycentric(x, y, v1, v2, v3, area_recip);
            
            if (InsideTri(&bary)) {
                F32 z = TriInterpolate(&bary, v1->z, v2->z, v3->z);
                
                F32 *pixel_depth = &depth[(USize)y_i*depth_width + (USize)x_i];
                if (z >= *pixel_depth) {
                    *pixel_depth = z;
                    
                    F32 r = TriInterpolate(&bary, vtx_colour[v1_i].x, vtx_colour[v2_i].x, vtx_colour[v3_i].x);
                    F32 g = TriInterpolate(&bary, vtx_colour[v1_i].y, vtx_colour[v2_i].y, vtx_colour[v3_i].y);
                    F32 b = TriInterpolate(&bary, vtx_colour[v1_i].z, vtx_colour[v2_i].z, vtx_colour[v3_i].z);
                    
                    F32 f = 255.f;
                    RGBA colour = {{
                        (U8)(r * f),
                        (U8)(g * f),
                        (U8)(b * f),
                        255.f,
                    }};
                    
                    SetPixel(tex, (USize)x_i, (USize)y_i, colour);
                }
            }
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
    0, 1, 2,
    1, 3, 2,
    2, 3, 7,
    3, 1, 7,
    7, 1, 5,
    1, 0, 5,
    5, 0, 4,
    0, 2, 4,
    4, 2, 6,
    2, 7, 6,
    6, 7, 4,
    7, 5, 4,
};

void SoftRender(Tex tex) {
    static F32 i = 0;
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
    
    timer_start();
    
    F32 t = i * 3.141592654f * 2.f / 1000.f;
    F32 cos_t = cosf(t);
    F32 sin_t = sinf(t);
    
    for (USize v_i = 0; v_i < countof(cube_verts); ++v_i) {
        Vtx v = cube_verts[v_i];
        
        // scale
        v.x *= 0.4f;
        v.y *= 0.4f;
        v.z *= 0.4f;
        
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
        
        // scale z from -1..1 to 0..1
        v.z = v.z * 0.5f + 0.5f;
        
        F32 aspect = (F32)tex.width / (F32)tex.height;
        v.x /= aspect;
        
        cube_verts_tform[v_i] = v;
    }
    
    for (USize idx_i = 0; idx_i < countof(cube_indices); idx_i += 3) {
        // RenderTri(tex, cube_verts_tform, cube_indices+idx_i, cube_colours);
        RenderTri_8(tex, cube_verts_tform, cube_indices+idx_i, cube_colours);
    }
    
    timer_elapsed("RenderTris");
}
