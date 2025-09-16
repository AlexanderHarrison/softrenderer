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
    
    if (-1e-5f < area && area < 1e-5f)
        return;
    if (area < 0.f)
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
    
    F32 vc[12] = {
        vtx_colour[v1_i].x, vtx_colour[v2_i].x, vtx_colour[v3_i].x, // r
        vtx_colour[v1_i].y, vtx_colour[v2_i].y, vtx_colour[v3_i].y, // g
        vtx_colour[v1_i].z, vtx_colour[v2_i].z, vtx_colour[v3_i].z, // b
    };
    for (int i = 0; i < 9; ++i)
        vc[i] *= 255.f;

    F32 slope_a = (v2x - v1x)/(v2y - v1y);
    F32 slope_b = (v3x - v1x)/(v3y - v1y);
    
    if (fabs(slope_a) < 10000.f && fabs(slope_b) < 10000.f) {
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
            
            for (ISize x_i = x_i_a; x_i < x_i_b; x_i += 8) {
                F32_8 x = _mm256_set1_ps((F32)x_i);
                F32_8 ints = _mm256_set_ps(7, 6, 5, 4, 3, 2, 1, 0);
                x = _mm256_add_ps(x, ints);
    
                F32_8 x_i_end_ymm = _mm256_set1_ps((F32)x_i_b);
                U32_8 in_bounds = (U32_8)_mm256_cmp_ps(x, x_i_end_ymm, _CMP_LT_OQ);
    
                x = _mm256_mul_ps(x, _mm256_set1_ps(x_inc));
                x = _mm256_add_ps(x, _mm256_set1_ps(x_start));
                
                F32_8 vw1_x = _mm256_set1_ps(v1x);
                F32_8 vw1_y = _mm256_set1_ps(v1y);
                F32_8 vw2_x = _mm256_set1_ps(v2x);
                F32_8 vw2_y = _mm256_set1_ps(v2y);
                F32_8 vw3_x = _mm256_set1_ps(v3x);
                F32_8 vw3_y = _mm256_set1_ps(v3y);
                F32_8 yw = _mm256_set1_ps(y);
                F32_8 a23 = CalcArea_8(x, yw, vw2_x, vw2_y, vw3_x, vw3_y);
                F32_8 a31 = CalcArea_8(vw1_x, vw1_y, x, yw, vw3_x, vw3_y);
                F32_8 a12 = CalcArea_8(vw1_x, vw1_y, vw2_x, vw2_y, x, yw);
                
                F32_8 area_recip_w = _mm256_set1_ps(area_recip);
                Vtx_8 bary = {
                    _mm256_mul_ps(a23, area_recip_w),
                    _mm256_mul_ps(a31, area_recip_w),
                    _mm256_mul_ps(a12, area_recip_w),
                };

                U32_8 tri_outside = OutsideTri_8(bary);
                U32_8 write_pixel = (U32_8)_mm256_andnot_ps((F32_8)tri_outside, (F32_8)in_bounds);

                int mask = _mm256_movemask_ps((F32_8)write_pixel);
                if (mask == 0)
                    continue;

                F32 *depth_row = &depth.pixels[(USize)y_i*depth.width + (USize)x_i];
                F32_8 cur_depth = _mm256_maskload_ps(depth_row, in_bounds);

                F32_8 z = TriInterpolate_8(bary, v1z, v2z, v3z);
                U32_8 z_good = (U32_8)_mm256_cmp_ps(z, cur_depth, _CMP_GE_OQ);
                write_pixel = (U32_8)_mm256_and_ps((F32_8)write_pixel, (F32_8)z_good);
            
                mask = _mm256_movemask_ps((F32_8)write_pixel);
                if (mask == 0)
                    continue;

                _mm256_maskstore_ps(depth_row, write_pixel, z);
                
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
    
    slope_a = (v2x - v3x)/(v2y - v3y);
    slope_b = (v1x - v3x)/(v1y - v3y);
    
    if (fabs(slope_a) < 10000.f && fabs(slope_b) < 10000.f) {
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
            
            for (ISize x_i = x_i_a; x_i < x_i_b; x_i += 8) {
                F32_8 x = _mm256_set1_ps((F32)x_i);
                F32_8 ints = _mm256_set_ps(7, 6, 5, 4, 3, 2, 1, 0);
                x = _mm256_add_ps(x, ints);
    
                F32_8 x_i_end_ymm = _mm256_set1_ps((F32)x_i_b);
                U32_8 in_bounds = (U32_8)_mm256_cmp_ps(x, x_i_end_ymm, _CMP_LT_OQ);
    
                x = _mm256_mul_ps(x, _mm256_set1_ps(x_inc));
                x = _mm256_add_ps(x, _mm256_set1_ps(x_start));
                
                F32_8 vw1_x = _mm256_set1_ps(v1x);
                F32_8 vw1_y = _mm256_set1_ps(v1y);
                F32_8 vw2_x = _mm256_set1_ps(v2x);
                F32_8 vw2_y = _mm256_set1_ps(v2y);
                F32_8 vw3_x = _mm256_set1_ps(v3x);
                F32_8 vw3_y = _mm256_set1_ps(v3y);
                F32_8 yw = _mm256_set1_ps(y);
                F32_8 a23 = CalcArea_8(x, yw, vw2_x, vw2_y, vw3_x, vw3_y);
                F32_8 a31 = CalcArea_8(vw1_x, vw1_y, x, yw, vw3_x, vw3_y);
                F32_8 a12 = CalcArea_8(vw1_x, vw1_y, vw2_x, vw2_y, x, yw);
                
                F32_8 area_recip_w = _mm256_set1_ps(area_recip);
                Vtx_8 bary = {
                    _mm256_mul_ps(a23, area_recip_w),
                    _mm256_mul_ps(a31, area_recip_w),
                    _mm256_mul_ps(a12, area_recip_w),
                };
                
                U32_8 tri_outside = OutsideTri_8(bary);
                U32_8 write_pixel = (U32_8)_mm256_andnot_ps((F32_8)tri_outside, (F32_8)in_bounds);

                int mask = _mm256_movemask_ps((F32_8)write_pixel);
                if (mask == 0)
                    continue;

                F32 *depth_row = &depth.pixels[(USize)y_i*depth.width + (USize)x_i];
                F32_8 cur_depth = _mm256_maskload_ps(depth_row, in_bounds);

                F32_8 z = TriInterpolate_8(bary, v1z, v2z, v3z);
                U32_8 z_good = (U32_8)_mm256_cmp_ps(z, cur_depth, _CMP_GE_OQ);
                write_pixel = (U32_8)_mm256_and_ps((F32_8)write_pixel, (F32_8)z_good);
            
                mask = _mm256_movemask_ps((F32_8)write_pixel);
                if (mask == 0)
                    continue;

                _mm256_maskstore_ps(depth_row, write_pixel, z);
                
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
    
    /*slope_a = (v2x - v3x)/(v2y - v3y);
    slope_b = (v1x - v3x)/(v1y - v3y);
    
    if (fabs(slope_a) < 10000.f && fabs(slope_b) < 10000.f) {
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
            
            for (ISize x_i = x_i_a; x_i < x_i_b; x_i += 8) {
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
    }*/
}
