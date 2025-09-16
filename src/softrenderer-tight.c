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

static inline void Swap(USize *a, USize *b) {
    USize tmp = *a;
    *a = *b;
    *b = tmp;
}

__attribute__((noinline))
static void RenderTri(Tex tex, Depth depth, Vtx *vtx, U16 *idx, Vtx *vtx_colour) {
    USize v1_i = (USize)idx[0];
    USize v2_i = (USize)idx[1];
    USize v3_i = (USize)idx[2];
    
    F32 area;
    {
        F32 v1x = vtx[v1_i].x; F32 v1y = vtx[v1_i].y;
        F32 v2x = vtx[v2_i].x; F32 v2y = vtx[v2_i].y;
        F32 v3x = vtx[v3_i].x; F32 v3y = vtx[v3_i].y;
        area = (v1x - v2x)*(v3y - v1y) - (v1x - v3x)*(v2y - v1y);
    }
    
    if (-1e-7f < area && area < 1e-7f)
        return;
    if (area < 0.f) // TODO: check
        return;

    F32 area_recip = 1.f / area;
    
    // order v1, v2, v3 by ascending y values
    if (vtx[v1_i].y > vtx[v2_i].y) { Swap(&v1_i, &v2_i); area_recip = -area_recip; }
    if (vtx[v2_i].y > vtx[v3_i].y) { Swap(&v2_i, &v3_i); area_recip = -area_recip; }
    if (vtx[v1_i].y > vtx[v2_i].y) { Swap(&v1_i, &v2_i); area_recip = -area_recip; }
    
    F32 w = (F32)tex.width;
    F32 h = (F32)tex.height;
    
    ISize y_i_a = (ISize)((vtx[v1_i].y + 1.f) * (h * 0.5f));
    ISize y_i_b = (ISize)((vtx[v2_i].y + 1.f) * (h * 0.5f));
    ISize y_i_c = (ISize)((vtx[v3_i].y + 1.f) * (h * 0.5f));
    
    if (y_i_a < 0) y_i_a = 0;
    if (y_i_b < 0) y_i_b = 0;
    if (y_i_c < 0) y_i_c = 0;
    if (y_i_a > (ISize)tex.height) y_i_a = (ISize)tex.height;
    if (y_i_b > (ISize)tex.height) y_i_b = (ISize)tex.height;
    if (y_i_c > (ISize)tex.height) y_i_c = (ISize)tex.height;
    
    F32 wr = 1.f / w;
    F32 hr = 1.f / h;
    F32 x_inc = 2.f * wr;
    F32 y_inc = 2.f * hr;
    F32 x_start = wr - 1.f;
    F32 y_start = hr - 1.f;
    
    F32 v1x = vtx[v1_i].x; F32 v1y = vtx[v1_i].y; F32 v1z = vtx[v1_i].z;
    F32 v2x = vtx[v2_i].x; F32 v2y = vtx[v2_i].y; F32 v2z = vtx[v2_i].z;
    F32 v3x = vtx[v3_i].x; F32 v3y = vtx[v3_i].y; F32 v3z = vtx[v3_i].z;
    
    Vtx *vc1 = vtx_colour + v1_i;
    Vtx *vc2 = vtx_colour + v2_i;
    Vtx *vc3 = vtx_colour + v3_i;
    
    F32 slope_a = (v2x - v1x)/(v2y - v1y);
    F32 slope_b = (v3x - v1x)/(v3y - v1y);
    
    for (ISize y_i = y_i_a; y_i < y_i_b; y_i += 1) {
        F32 y = y_start + y_inc*(F32)y_i;

        F32 tri_x_1 = v1x + slope_a*(y - v1y);
        F32 tri_x_2 = v1x + slope_b*(y - v1y);
        
        F32 tri_x_a = tri_x_1 < tri_x_2 ? tri_x_1 : tri_x_2;
        F32 tri_x_b = tri_x_1 > tri_x_2 ? tri_x_1 : tri_x_2;

        ISize x_i_a = (ISize)((tri_x_a + 1.f) * (w * 0.5f));
        ISize x_i_b = (ISize)((tri_x_b + 1.f) * (w * 0.5f)) + 1;
    
        if (x_i_a < 0) x_i_a = 0;
        if (x_i_b < 0) x_i_b = 0;
        if (x_i_a > (ISize)tex.width) x_i_a = (ISize)tex.width;
        if (x_i_b > (ISize)tex.width) x_i_b = (ISize)tex.width;
        
        for (ISize x_i = x_i_a; x_i < x_i_b; x_i += 1) {
            F32 x = x_start + x_inc*(F32)x_i;
            
            F32 a23 = CalcArea(x, y, v2x, v2y, v3x, v3y);
            F32 a31 = CalcArea(v1x, v1y, x, y, v3x, v3y);
            F32 a12 = CalcArea(v1x, v1y, v2x, v2y, x, y);
            Vtx bary = { a23 * area_recip, a31 * area_recip, a12 * area_recip };
            
            // prevent pathological triangles
            if (!InsideTri(&bary)) continue;
            
            F32 z = TriInterpolate(&bary, v1z, v2z, v3z);
            F32 *pixel_depth = &depth.pixels[(USize)y_i*depth.width + (USize)x_i];
            if (z >= *pixel_depth) {
                *pixel_depth = z;
                
                F32 r = TriInterpolate(&bary, vc1->x, vc2->x, vc3->x);
                F32 g = TriInterpolate(&bary, vc1->y, vc2->y, vc3->y);
                F32 b = TriInterpolate(&bary, vc1->z, vc2->z, vc3->z);
                
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
    
    slope_a = (v2x - v3x)/(v2y - v3y);
    slope_b = (v1x - v3x)/(v1y - v3y);
    
    for (ISize y_i = y_i_b; y_i < y_i_c; y_i += 1) {
        F32 y = y_start + y_inc*(F32)y_i;

        F32 tri_x_1 = v3x + slope_a*(y - v3y);
        F32 tri_x_2 = v3x + slope_b*(y - v3y);
        
        F32 tri_x_a = tri_x_1 < tri_x_2 ? tri_x_1 : tri_x_2;
        F32 tri_x_b = tri_x_1 > tri_x_2 ? tri_x_1 : tri_x_2;

        ISize x_i_a = (ISize)((tri_x_a + 1.f) * (w * 0.5f));
        ISize x_i_b = (ISize)((tri_x_b + 1.f) * (w * 0.5f)) + 1;
    
        if (x_i_a < 0) x_i_a = 0;
        if (x_i_b < 0) x_i_b = 0;
        if (x_i_a > (ISize)tex.width) x_i_a = (ISize)tex.width;
        if (x_i_b > (ISize)tex.width) x_i_b = (ISize)tex.width;
        
        for (ISize x_i = x_i_a; x_i < x_i_b; x_i += 1) {
            F32 x = x_start + x_inc*(F32)x_i;
            
            F32 a23 = CalcArea(x, y, v2x, v2y, v3x, v3y);
            F32 a31 = CalcArea(v1x, v1y, x, y, v3x, v3y);
            F32 a12 = CalcArea(v1x, v1y, v2x, v2y, x, y);
            Vtx bary = { a23 * area_recip, a31 * area_recip, a12 * area_recip };
            
            // prevent pathological triangles
            if (!InsideTri(&bary)) continue;
            
            F32 z = TriInterpolate(&bary, v1z, v2z, v3z);
            F32 *pixel_depth = &depth.pixels[(USize)y_i*depth.width + (USize)x_i];
            if (z >= *pixel_depth) {
                *pixel_depth = z;
                
                F32 r = TriInterpolate(&bary, vc1->x, vc2->x, vc3->x);
                F32 g = TriInterpolate(&bary, vc1->y, vc2->y, vc3->y);
                F32 b = TriInterpolate(&bary, vc1->z, vc2->z, vc3->z);
                
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
