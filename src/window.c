#define SDL_MAIN_USE_CALLBACKS 1
#include "SDL3/SDL.h"
#include "SDL3/SDL_main.h"

#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <byteswap.h>
#include <math.h>
#include <immintrin.h>

#define countof(A) ((sizeof(A)/sizeof(*(A))))

typedef float F32;
typedef double F64;
typedef size_t USize;
typedef uint64_t U64;
typedef uint32_t U32;
typedef uint16_t U16;
typedef uint8_t U8;
typedef ssize_t ISize;
typedef int64_t I64;
typedef int32_t I32;
typedef int16_t I16;
typedef int8_t I8;

// #include "softrenderer.c"
#include "softrenderer-avx.c"

static SDL_Window *window;
static SDL_Renderer *renderer;
static SDL_Texture *rgba;

// static USize window_width = 640;
// static USize window_height = 480;
static USize window_width = 1920;
static USize window_height = 1080;

#define Expect_SDL(A) do {\
    if (!(A)) {\
        fprintf(stderr, "error in '" #A "' : %s\n", SDL_GetError());\
        return SDL_APP_FAILURE;\
    }\
} while (0)

SDL_AppResult SDL_AppInit(void **appstate, int argc, char *argv[])
{
    (void)appstate;
    (void)argc;
    (void)argv;

    SDL_SetAppMetadata("Software Renderer", "1.0", "com.example.software-renderer");

    Expect_SDL(
        SDL_Init(SDL_INIT_VIDEO)
    );
    Expect_SDL(
        SDL_CreateWindowAndRenderer("software renderer", (int)window_width, (int)window_height, 0, &window, &renderer)
    );
    Expect_SDL(
        SDL_SetRenderVSync(renderer, true)
    );
    
    rgba = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, (int)window_width, (int)window_height);
    Expect_SDL(rgba != NULL);
    
    return SDL_APP_CONTINUE;
}

SDL_AppResult SDL_AppEvent(void *appstate, SDL_Event *event)
{
    (void)appstate;
    
    if (event->type == SDL_EVENT_QUIT)
        return SDL_APP_SUCCESS;
    return SDL_APP_CONTINUE;
}

SDL_AppResult SDL_AppIterate(void *appstate)
{
    (void)appstate;
    
    U8 *pixels;
    int pitch;
    Expect_SDL(
        SDL_LockTexture(rgba, NULL, (void**)&pixels, &pitch)
    );
    
    SoftRender((Tex) { pixels, window_width, window_height, (USize)pitch });
    
    SDL_UnlockTexture(rgba);
    
    Expect_SDL(
        SDL_RenderTexture(renderer, rgba, NULL, NULL)
    );
    
    SDL_RenderPresent(renderer);
 
    return SDL_APP_CONTINUE;
}

void SDL_AppQuit(void *appstate, SDL_AppResult result)
{
    (void)appstate;
    (void)result;
}
