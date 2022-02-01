

#include <raylib.h>

#if defined(PLATFORM_WEB)
    #include <emscripten/emscripten.h>
#endif

void UpdateDrawFrame()
{
    BeginDrawing();

    ClearBackground(CLITERAL(Color){ 0, 200, 0, 50  });

    EndDrawing();

}

int main(int argc, char** argv)
{
    const int screenWidth = 256;
    const int screenHeight = 256;
    const int screenScale = 2;

    SetAudioStreamBufferSizeDefault(14400 / 60 * 6);
    SetConfigFlags(FLAG_WINDOW_TRANSPARENT);
    InitWindow(screenWidth * screenScale, screenHeight * screenScale, "Nice");

    #if defined(PLATFORM_WEB)
    emscripten_set_main_loop(UpdateDrawFrame, 60, 1);
    #else
    
    SetTargetFPS(60);

    
    while (!WindowShouldClose())
    {
        UpdateDrawFrame();
    }
    #endif

    return 0;
}
