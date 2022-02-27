#include "scenes.hpp"

#include <raylib.h>
#include <raymath.h>
#include <rlgl.h>
#include <glm/mat4x4.hpp>
#include <glm/ext.hpp>

#if defined(PLATFORM_WEB)
    #include <emscripten/emscripten.h>
#endif

#include <iostream>

orni::GameState g_gameState;
std::function<void(orni::GameState&)> g_sceneFunc;

void load_default_scene()
{
     g_sceneFunc = orni::gen_test_scene_c();
}

void update_draw_frame()
{
    if (IsKeyDown(KEY_ONE))
    {
        g_sceneFunc = orni::gen_test_scene_a();
    }
    else if (IsKeyDown(KEY_TWO))
    {
        g_sceneFunc = orni::gen_test_scene_b();
    }
    else if (IsKeyDown(KEY_THREE))
    {
        g_sceneFunc = orni::gen_test_scene_b_b();
    }
    else if (IsKeyDown(KEY_FOUR))
    {
        g_sceneFunc = orni::gen_test_scene_c();
    }
    g_sceneFunc(g_gameState);
}

int main(int argc, char** argv)
{
    const int screenWidth = 1280;
    const int screenHeight = 720;

    SetAudioStreamBufferSizeDefault(14400 / 60 * 6);
    SetConfigFlags(FLAG_WINDOW_TRANSPARENT);
    InitWindow(screenWidth, screenHeight, "Nice");

    // default resources
    static Font defaultFont = LoadFontEx("Tuffy-Bold.ttf", 20, 0, 0);
    SetTextureFilter(defaultFont.texture, TEXTURE_FILTER_BILINEAR);
    g_gameState.m_pFont = &defaultFont;

    load_default_scene();

#if defined(PLATFORM_WEB)
    emscripten_set_main_loop(update_draw_frame, 0, 1);
#else

    SetTargetFPS(60);

    while (!WindowShouldClose())
    {
        update_draw_frame();
    }

#endif

    return 0;
}
