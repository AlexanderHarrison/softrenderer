

typedef struct Tex {
    U8 *pixels;
    USize width, height, pitch;
} Tex;

typedef struct Depth {
    F32 *pixels;
    USize width, height, pitch;
} Depth;

typedef union RGBA {
    struct {
        U8 r, g, b, a;
    };
    U32 u32;
} RGBA;

typedef struct Vtx {
    F32 x, y, z;
} Vtx;

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

static inline F32 TriInterpolate(Vtx *bary, F32 f1, F32 f2, F32 f3) {
    F32 z1 = f1*bary->x;
    F32 z2 = f2*bary->y;
    F32 z3 = f3*bary->z;
    F32 sum = z1 + z2 + z3;
    return sum;
}

static inline bool InsideTri(Vtx *bary) {
    U32 b1 = FloatBits(bary->x);
    U32 b2 = FloatBits(bary->y); 
    U32 b3 = FloatBits(bary->z);
    U32 ors = b1 | b2 | b3;
    return (ors >> 31) == 0;
}

__attribute__((noinline))
static void RenderTri(Tex tex, Depth depth, Vtx *vtx, U16 *idx, Vtx *vtx_colour) {
    F32 w = (F32)tex.width;
    F32 h = (F32)tex.height;
    F32 wr = 1.f / w;
    F32 hr = 1.f / h;
    
    F32 x_min = 1.f;
    F32 y_min = 1.f;
    F32 x_max = -1.f;
    F32 y_max = -1.f;
    
    USize v1_i = (USize)idx[0];
    USize v2_i = (USize)idx[1];
    USize v3_i = (USize)idx[2];
    
    Vtx *v1 = vtx + v1_i;
    Vtx *v2 = vtx + v2_i;
    Vtx *v3 = vtx + v3_i;
    
    F32 area = (v1->x - v2->x)*(v3->y - v1->y) - (v1->x - v3->x)*(v2->y - v1->y);
    if (-1e-7f < area && area < 1e-7f)
        return;
    
    if (area < 0.f)
        return;
    F32 area_recip = 1.f / area;

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
                
                F32 *pixel_depth = &depth.pixels[(USize)y_i*depth.width + (USize)x_i];
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
