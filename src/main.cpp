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


void update_draw_frame()
{
    g_sceneFunc(g_gameState);
}

int main(int argc, char** argv)
{
    const int screenWidth = 800;
    const int screenHeight = 600;

    SetAudioStreamBufferSizeDefault(14400 / 60 * 6);
    SetConfigFlags(FLAG_WINDOW_TRANSPARENT);
    InitWindow(screenWidth, screenHeight, "Nice");

    g_sceneFunc = orni::gen_test_scene_a();

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
