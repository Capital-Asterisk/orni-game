#include "scenes.hpp"
#include "salads.hpp"

#include <raylib.h>
#include <raymath.h>
#include <rlgl.h>

// todo
#if defined(PLATFORM_DESKTOP)
    #define GLSL_VERSION            330
#else   // PLATFORM_RPI, PLATFORM_ANDROID, PLATFORM_WEB
    #define GLSL_VERSION            100
#endif
#define RLIGHTS_IMPLEMENTATION
extern "C"
{
#include "rlights.h"
}

#include <glm/gtx/transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/matrix_interpolation.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtx/vector_angle.hpp>

#include <iostream>
#include <memory>
#include <unordered_map>
#include <sstream>
#include <string_view>


#if defined(PLATFORM_WEB)
    #include <emscripten/emscripten.h>
#endif

namespace orni
{



struct MainGameScene : ThreadedLoop
{
    ~MainGameScene()
    {
        m_running = false;
        m_syncFramesCv.notify_all();

        // LOL
        for (Image &rImage : m_gltfRayImage)
        {
            UnloadImage(rImage);
        }
        for (Texture &rTex : m_gltfRayTextures)
        {
            UnloadTexture(rTex);
        }
        for (Material &rMat : m_gltfRayMaterials)
        {
            RL_FREE(rMat.maps);
        }
        UnloadShader(m_shaderLit);
        UnloadMesh(m_cube);
        UnloadRenderTexture(m_ui);
        UnloadMaterial(m_mat);

        for (std::unique_ptr<SaladModel>& pSalad : m_salads)
        {
            if (!bool(pSalad))
            {
                continue;
            }
            rlUnloadVertexArray(pSalad->m_rayMesh.vaoId);
        }


        m_updater.join();
    }

    Salads_t                m_salads;
    FrogDyn                 m_frogs;

    Characters_t            m_characters;

    Inputs                  m_inputs;
    ToolGrab                m_toolGrab;
    ToolGrabRotater         m_toolGrabRotate;
    std::vector<GrabDisplay> m_toolGrabDisplay;
    bool                    m_mouseMoved;

    tinygltf::Model         m_gltf;
    std::vector<Image>      m_gltfRayImage;
    std::vector<Texture>    m_gltfRayTextures;
    std::vector<Material>   m_gltfRayMaterials;

    Camera3D                m_camera;
    glm::vec3               m_cameraOffset;

    Mesh                    m_cube;
    Material                m_mat;
    Shader                  m_shaderLit;
    std::array<Light, 4>    m_lights;

    RenderTexture2D         m_ui;

    int                     m_cstSteps = 17;
    float                   m_extPercent = 0.5f;

    float                   m_camDist{3.0f};
    float                   m_camYaw{0.78539f};
    float                   m_camPitch{0.0f};
    float                   m_time{0.0f};
    bool                    m_gravity{true};
};

inline bool g_printEmotions = false;

static void update_scene(MainGameScene &rScene, GameState &rGame)
{

    float const delta = 1.0f / 60.0f;
    CharB &rChar = rScene.m_characters.begin()->second;

    if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT))
    {
        Vector2 const mouseDelta = GetMouseDelta();
        float const pi = glm::pi<float>();
        float const sensitivity = 0.5f;

        rScene.m_camPitch = glm::clamp(rScene.m_camPitch - mouseDelta.y * delta * sensitivity, -pi/2.1f, pi/2.1f) ;
        rScene.m_camYaw = glm::mod(rScene.m_camYaw - mouseDelta.x * delta * sensitivity, pi * 2.0f) ;
    }

    {
        float const sensitivity = 0.2f;
        rScene.m_camDist = glm::clamp(rScene.m_camDist - rScene.m_camDist * GetMouseWheelMove() * sensitivity, 1.0f, 12.0f );
    }

    if (!IsKeyDown(KEY_LEFT_CONTROL))
    { }
    if (IsKeyPressed(KEY_E))
    {
        g_printEmotions = !g_printEmotions;
    }

    if (IsKeyPressed(KEY_ONE))
    {
        rScene.m_inputs.m_selected = 0;
    }
    else if (IsKeyPressed(KEY_TWO))
    {
        rScene.m_inputs.m_selected = 1;
    }

    glm::vec3 avgPos{0.0f};
    float totalMass = 0.0f;
    // draw frogs
    for (frog_id_t id = 0; id < rScene.m_frogs.m_ids.capacity(); id ++)
    {
        if (!rScene.m_frogs.m_ids.exists(id))
        {
            continue;
        }
        avgPos += glm::vec3(rScene.m_frogs.m_tf[id][3]) * rScene.m_frogs.m_mass[id];
        totalMass += rScene.m_frogs.m_mass[id];
    }

    avgPos /= totalMass;

    glm::vec3 &rCamPos = reinterpret_cast<glm::vec3&>(rScene.m_camera.position);
    glm::vec3 &rCamTgt = reinterpret_cast<glm::vec3&>(rScene.m_camera.target);
    glm::vec3 const diff = avgPos + rScene.m_cameraOffset - rCamTgt;
    // smoothly move camera to center of mass

    //rCamPos += diff * 0.5f;
    rCamTgt += diff * 0.01f;

    rCamPos = glm::quat(glm::vec3{rScene.m_camPitch, rScene.m_camYaw, 0.0f}) * glm::vec3{0.0f, 0.0f, rScene.m_camDist} + rCamTgt;

    update_inputs_rl(rScene.m_camera, rScene.m_inputs);

    rScene.m_inputs.m_lazor.m_salad = -1;
    if (update_tool_grab(rScene.m_salads, rChar.m_wetJoints, rScene.m_frogs, rScene.m_inputs, rScene.m_toolGrab))
    { }
    else if (update_tool_grab_rotate(rScene.m_toolGrabRotate, rScene.m_toolGrab, rScene.m_inputs, rScene.m_frogs, rScene.m_camera, delta))
    { }
    else if (offset_camera_lazor(rScene.m_salads, rScene.m_inputs, rScene.m_mouseMoved, rScene.m_cameraOffset, avgPos))
    { }
    update_tool_grab_pos(rScene.m_camera, rScene.m_inputs, rScene.m_frogs, rScene.m_toolGrab, delta);


    update_apples(rChar.m_apples, rChar.m_joints);

    rChar.m_eyeL.m_texturePos = calc_eye_pos(rChar.m_apples.m_dataOut[rChar.m_eyeL.m_apple], reinterpret_cast<glm::vec3&>(rScene.m_camera.position));
    rChar.m_eyeR.m_texturePos = calc_eye_pos(rChar.m_apples.m_dataOut[rChar.m_eyeR.m_apple], reinterpret_cast<glm::vec3&>(rScene.m_camera.position));

    float pushMag = (float(IsKeyDown(KEY_D)) - float(IsKeyDown(KEY_A))) * delta * 5.0f;

    glm::vec3 const push = glm::vec3{
                glm::sin(rScene.m_camYaw + 1.5708f),
                0.0f,
                glm::cos(rScene.m_camYaw + 1.5708f)} * pushMag;


    if (IsKeyPressed(KEY_G))
    {
        rScene.m_gravity = !rScene.m_gravity;
    }

    if (IsKeyDown(KEY_LEFT_ALT))
    {
        for (int i = 0; i < rScene.m_frogs.m_ids.size(); i ++)
        {
            rScene.m_frogs.m_vel[i].m_ang *= 0.50f;
        }
    }

    for (int i = 0; i < rScene.m_frogs.m_ids.capacity(); i ++)
    {
        // apply gravity
        rScene.m_frogs.m_extImp[i].m_lin.y -= 9.81f * delta * rScene.m_frogs.m_mass[i] * rScene.m_gravity;

        rScene.m_frogs.m_extImp[i].m_lin += push;

        auto p = glm::vec3(rScene.m_frogs.m_tf[i][3]);
    }

    float const cstPercent = 1.0f - rScene.m_extPercent;

    apply_ext_forces(rScene.m_frogs, delta * rScene.m_extPercent);

    // repeat constrain forces a few times
    for (int j = 0; j < rScene.m_cstSteps; j++)
    {
        float smldelta = delta * cstPercent / float(rScene.m_cstSteps);

        apply_baits(rScene.m_frogs, {64.0f, 0.25f, 16.0f, 0.5f, 8.0f, 0.2f, 8.0f, 0.25f}, smldelta);

        // beak-head spring workaround
        glm::vec3 const axis = rScene.m_frogs.m_tf[rChar.m_frogBeak][0];
        rScene.m_frogs.m_cstImp[rChar.m_frogBeak].m_ang -= axis * smldelta * 2.0f;
        rScene.m_frogs.m_cstImp[rChar.m_frogHead].m_ang += axis * smldelta * 2.0f;

        // stable head
        rScene.m_frogs.m_vel[rChar.m_frogHead].m_ang = frogdyn::oppose(rScene.m_frogs.m_vel[rChar.m_frogHead].m_ang, smldelta * 10.0f);

        apply_cst_forces(rScene.m_frogs, smldelta);

        calc_balls_pos(rScene.m_frogs);

        calc_frog_collisions(rScene.m_frogs);
    }

    // apply drag
    for (int i = 0; i < rScene.m_frogs.m_ids.size(); i ++)
    {
        rScene.m_frogs.m_vel[i].m_ang = frogdyn::oppose(rScene.m_frogs.m_vel[i].m_ang, delta * 0.2f);
        rScene.m_frogs.m_vel[i].m_lin = frogdyn::oppose(rScene.m_frogs.m_vel[i].m_lin, delta * 0.2f);
        //rScene.m_frogs.m_vel[i].m_lin *= 0.998f;
    }

    update_hoppers(rChar.m_wetJoints.m_hoppers, rScene.m_frogs, rChar.m_joints.m_nodeTf.data());

    meshdeform::calculate_joint_transforms(
            glm::mat4x4(1.0f),
            rChar.m_joints.m_pInverseBindIn,
            rChar.m_joints.m_nodeTf.data(),
            0,
            rChar.m_joints.m_jointTf.size(),
            rChar.m_joints.m_jointTf.data());

    for (int i = 0; i < 4; i ++)
    {

        meshdeform::apply_vertex_transform(
                rScene.m_characters[0].m_joints,
                rScene.m_salads[i]->m_spookM,
                rScene.m_salads[i]->m_tgt,
                rScene.m_salads[i]->m_pPosIn,
                rScene.m_salads[i]->m_pNrmIn,
                0,
                rScene.m_salads[i]->m_rayMesh.vertexCount,
                rScene.m_salads[i]->m_Pos.data(),
                rScene.m_salads[i]->m_Nrm.data());


    }

    rChar.m_soul.m_emPissed += (glm::length(rScene.m_frogs.m_vel[rChar.m_frogHead].m_ang) > 4.0f) * 0.2f * delta;
    rChar.m_soul.m_emWorry = glm::max(0.0f, rScene.m_toolGrab.m_grabs.size() - 7.0f) * 0.1f + rScene.m_toolGrab.m_grabs.size() * (!g_limits) * 0.25f;
    update_expressions(rChar.m_soul, delta);

    // breathing
    rScene.m_frogs.m_scale[rChar.m_frogBelly] = 1.0f + 0.1f * breath_cycle(rChar.m_soul.m_breathCycle);
}

static void update_scene_loop(MainGameScene &rScene, GameState &rGame)
{
    while (rScene.m_running.load())
    {
        {
            // Wait for render to finish
            {
                std::unique_lock<std::mutex> lock(rScene.m_syncFramesMtx);
                // wait until frame finishes rendering
                rScene.m_syncFramesCv.wait(lock, [&rScene]  { return rScene.m_framesRendered == rScene.m_framesUpdated || (!rScene.m_running.load()); }  );
            }

            //std::lock_guard<std::mutex>(rScene.m_updaterBusy);
            //std::cout << "doing update...\n";
            update_scene(rScene, rGame);

            //std::this_thread::sleep_for(std::chrono::milliseconds(1000));

            // Update finished
            {
                std::unique_lock<std::mutex> lock(rScene.m_syncFramesMtx);
                ++rScene.m_framesUpdated;
            }
            rScene.m_syncFramesCv.notify_all();

        }
    }
}

static void draw_scene(MainGameScene &rScene, GameState &rGame, float threading)
{
    // wait for new frame to be updated
    if (threading)
    {
        std::unique_lock<std::mutex> lock(rScene.m_syncFramesMtx);
        rScene.m_syncFramesCv.wait(lock, [&rScene]  { return rScene.m_framesRendered < rScene.m_framesUpdated; }  );
    }

    // acquire variables and upload stuff needed to render, this 'completes' rendering

    CharB const &rChar = rScene.m_characters.begin()->second;
    Soul const soul = rChar.m_soul;
    glm::vec2 const eyePosR = rChar.m_eyeR.m_texturePos;
    glm::vec2 const eyePosL = rChar.m_eyeL.m_texturePos;
    Inputs const inputs = rScene.m_inputs;
    Camera const cam = rScene.m_camera;


    SetShaderValue(rScene.m_shaderLit, rScene.m_shaderLit.locs[SHADER_LOC_VECTOR_VIEW], &rScene.m_camera.position.x, SHADER_UNIFORM_VEC3);
    UpdateLightValues(rScene.m_shaderLit, rScene.m_lights[0]);
    UpdateLightValues(rScene.m_shaderLit, rScene.m_lights[1]);
    UpdateLightValues(rScene.m_shaderLit, rScene.m_lights[2]);

    update_grab_displays(rScene.m_toolGrab, rScene.m_frogs, rScene.m_toolGrabDisplay);

    for (int i = 0; i < 4; i ++)
    {
        rlUpdateVertexBuffer(rScene.m_salads[i]->m_rayMesh.vboId[0], rScene.m_salads[i]->m_rayMesh.vertices, rScene.m_salads[i]->m_rayMesh.vertexCount*3*sizeof(float), 0);    // Update vertex position
        rlUpdateVertexBuffer(rScene.m_salads[i]->m_rayMesh.vboId[2], rScene.m_salads[i]->m_rayMesh.normals, rScene.m_salads[i]->m_rayMesh.vertexCount*3*sizeof(float), 0);     // Update vertex normals
    }

    if (threading)
    {
        std::unique_lock<std::mutex> lock(rScene.m_syncFramesMtx);
        ++rScene.m_framesRendered;
        lock.unlock();
        rScene.m_syncFramesCv.notify_all(); // update of next frame can resume
    }


    BeginDrawing();


        ClearBackground(Color{ 10, 10, 10, 100 });

        // draw eyes


        BeginTextureMode(rChar.m_eyeTexture);

            ClearBackground(Color{ 0, 0, 0, 0 });

            int const g = 128;
            int gx;
            int gy;

            // blink for 5 frames
            if (soul.m_blinkCdn < 5.0f/60.0f)
            {
                gx = 0; gy = 2;
            }
            else
            {
                if (soul.m_emWorry > 0.6f)
                {
                    gx = 4; gy = 1;
                }
                else if (soul.m_emWorry > 0.3f)
                {
                    gx = 3; gy = 1;
                }
                else if (soul.m_emPissed > 0.5f)
                {
                    gx = 1; gy = 1;
                }
                else
                {
                    gx = 0; gy = 1;
                }
            }

            float ox = gx * g;
            float oy = gy * g * 2;
            float wy = oy + g;
            Rectangle const rectR{0, g, g, g};
            Rectangle const rectL{g, g, g, g};


            // white stuff 1
            DrawTexturePro(rChar.m_eyeSheet, {ox, wy, -g, -g}, rectR, {}, 0.0f, WHITE);
            DrawTexturePro(rChar.m_eyeSheet, {ox, wy, g, -g}, rectL, {}, 0.0f, WHITE);

            // irises
            draw_iris(rChar.m_eyeSheet, 0, eyePosR);
            draw_iris(rChar.m_eyeSheet, 1, eyePosL);

            // erase overlap
            BeginBlendMode(BLEND_CUSTOM);
                rlSetBlendFactors(0, 770, 32774);
                DrawTexturePro(rChar.m_eyeSheet, {ox, wy, -g, -g}, rectR, {}, 0.0f, WHITE);
                DrawTexturePro(rChar.m_eyeSheet, {ox, wy, g, -g}, rectL, {}, 0.0f, WHITE);
            EndBlendMode();

            // outline
            DrawTexturePro(rChar.m_eyeSheet, {ox, oy, -g, -g}, rectR, {}, 0.0f, WHITE);
            DrawTexturePro(rChar.m_eyeSheet, {ox, oy, g, -g}, rectL, {}, 0.0f, WHITE);
        EndTextureMode();

        BeginMode3D(cam);
            BeginBlendMode(BLEND_ALPHA);
                DrawGrid(20, 1.0f);

                rlDisableBackfaceCulling();
                DrawModel(rScene.m_salads[0]->m_rayModel, Vector3{0.0f, 0.0f, 0.0f}, 1.0f, WHITE);
                DrawModel(rScene.m_salads[2]->m_rayModel, Vector3{0.0f, 0.0f, 0.0f}, 1.0f, WHITE);
                DrawModel(rScene.m_salads[3]->m_rayModel, Vector3{0.0f, 0.0f, 0.0f}, 1.0f, WHITE);
                DrawModel(rScene.m_salads[1]->m_rayModel, Vector3{0.0f, 0.0f, 0.0f}, 1.0f, WHITE);

                //glm::vec3 lol = rScene.m_inputs.m_mouseOrig + rScene.m_inputs.m_mouseDir * rScene.m_inputs.m_lazor.m_mcray.m_dist;
                //DrawSphereWires(reinterpret_cast<Vector3&>(lol), 0.01f, 4, 4, GREEN);

            EndBlendMode();
        EndMode3D();

        BeginTextureMode(rScene.m_ui);
            ClearBackground(Color{ 0, 0, 0, 0 });

            BeginBlendMode(BLEND_ALPHA);

            for (GrabDisplay const& grab : rScene.m_toolGrabDisplay)
            {
                if (grab.m_visible)
                {
                    unsigned char alpha = glm::clamp(350.0f - glm::distance(grab.m_screenPos, inputs.m_mousePos) * 0.5f, 100.0f, 220.0f);
                    Color color = grab.m_selected ? Color{100, 255, 100, alpha} : Color{255, 255, 255, alpha};
                    DrawCircle(grab.m_screenPos.x, grab.m_screenPos.y, ToolGrab::smc_radius, color);
                    DrawCircleLines(grab.m_screenPos.x, grab.m_screenPos.y, ToolGrab::smc_radius, color);
                }
                //DrawSphereWires(reinterpret_cast<Vector3 const&>(grab.m_pos), 0.1, 5, 6, Color{255, 255, 255, 255});
            }

            EndBlendMode();

            BeginMode3D(cam);



#if 0
                // wait first to avoid race condition when reading frogs
                if (threading)
                {
                    std::unique_lock<std::mutex> lock(rScene.m_syncFramesMtx);
                    rScene.m_syncFramesCv.wait(lock, [&rScene]  { return rScene.m_framesRendered < rScene.m_framesUpdated; }  );
                }

                // draw frogs
                for (frog_id_t id = 0; id < rScene.m_frogs.m_ids.capacity(); id ++)
                {
                    if (!rScene.m_frogs.m_ids.exists(id))
                    {
                        continue;
                    }

                    auto const &tf = rScene.m_frogs.m_tf[id];
                    auto transposed = glm::transpose(tf);

                    auto tip = glm::vec3(tf[3] - tf[1] * 0.5f);
                    auto tail = glm::vec3(tf[3] + tf[1] * 0.5f);

                    DrawMesh(rScene.m_cube, rScene.m_mat, reinterpret_cast<Matrix&>(transposed));

                    // draw the balls if the frog has balls
                    if (rScene.m_frogs.m_ballPos.contains(id))
                    {
                        // loop throught the balls
                        for (auto& ball : rScene.m_frogs.m_ballPos[id])
                        {
                            DrawSphereWires(reinterpret_cast<Vector3&>(ball.m_pos), ball.m_radius, 5, 6, Color{255, 255, 255, 255});
                        }
                    }
                }

                // draw baits
                for (int i = 0; i < rScene.m_frogs.m_baits.size(); i ++)
                {
                    auto const &rBait = rScene.m_frogs.m_baits[i];
                    glm::vec3 posA, velA, posB, velB;

                    if (rBait.m_a.m_id == -1)
                    {
                        // world-anchored
                        posA = rBait.m_a.m_offset;
                        velA = glm::vec3{0.0f};
                    }
                    else
                    {
                        glm::vec3 const offsetRotated(rScene.m_frogs.m_tf[rBait.m_a.m_id] * glm::vec4(rBait.m_a.m_offset, 0.0f));
                        posA = offsetRotated + glm::vec3(rScene.m_frogs.m_tf[rBait.m_a.m_id][3]);
                    }

                    glm::vec3 const offsetRotated(rScene.m_frogs.m_tf[rBait.m_b.m_id] * glm::vec4(rBait.m_b.m_offset, 0.0f));
                    posB = offsetRotated + glm::vec3(rScene.m_frogs.m_tf[rBait.m_b.m_id][3]);

                    DrawSphereWires(reinterpret_cast<Vector3&>(posA), 0.01f, 2, 4, Color{255, 0, 0, 255});
                    DrawSphere(reinterpret_cast<Vector3&>(posB), 0.005f, Color{0, 0, 255, 255});


                    if (rBait.m_a.m_id != -1)
                    {
                        Vector3 origA = reinterpret_cast<Vector3&>(rScene.m_frogs.m_tf[rBait.m_a.m_id][3]);
                        DrawLine3D(origA, reinterpret_cast<Vector3&>(posA), Color{255, 0, 0, 255});
                    }
                    Vector3 origB = reinterpret_cast<Vector3&>(rScene.m_frogs.m_tf[rBait.m_b.m_id][3]);
                    DrawLine3D(origB, reinterpret_cast<Vector3&>(posB), Color{0, 0, 255, 255});

                }

#endif
            EndMode3D();

            std::ostringstream readouts;
            readouts << "";
            if (g_printEmotions)
            {
                readouts << "Pissed:  " << int(soul.m_emPissed * 100) << "/100\n";
                readouts << "Worried: " << int(soul.m_emWorry * 100) << "/100\n";
                readouts << "Shy:     " << int(soul.m_emShy * 100) << "/100\n";
                readouts << "BadIdea: " << int(soul.m_emBadIdea * 100) << "/100\n";
            }

            DrawTextEx(*rGame.m_pFont, readouts.str().c_str(), Vector2{10.0, 100.0}, 20, 0, WHITE);

        EndTextureMode();

        auto &rTex = rScene.m_ui.texture;
        DrawTextureRec(rTex, Rectangle{0.0f, 0.0f, float(rTex.width), -float(rTex.height) * 1.0f}, Vector2{0.0f, 0.0f}, Color{255, 255, 255, 255});
        EndDrawing();
}

SceneFunc_t gen_the_game_scene(GameState &rGame)
{
    std::shared_ptr<MainGameScene> pScene = std::make_shared<MainGameScene>();
    MainGameScene &rScene = *pScene;

    std::string err;
    std::string warn;

    tinygltf::TinyGLTF loader;

    rScene.m_toolGrab.m_id = 0;
    rScene.m_toolGrabRotate.m_id = 1;
    rScene.m_inputs.m_selected = 0;

    rScene.m_gltfRayImage.resize(10);
    loader.SetImageLoader([] (
            tinygltf::Image *img, const int index,
            std::string *pErr, std::string *pWarn, int foo, int bar,
            const unsigned char *pData, int size, void *pUser) -> bool
    {
        auto *pImages = reinterpret_cast<std::vector<Image>*>(pUser);
        pImages->at(index) = LoadImageFromMemory(".png", pData, size);
        return true;
    }, &rScene.m_gltfRayImage);

    //image, image_idx, err, warn, 0, 0, &img.at(0),
    //static_cast<int>(img.size()), load_image_user_data);

    loader.LoadASCIIFromFile(&rScene.m_gltf, &err, &warn, "resources/gltf/salad0.gltf");

    // funny patch LOL!!!
    for (auto &rNode : rScene.m_gltf.nodes)
    {
        if (rNode.skin != -1)
        {
            rNode.skin = 0;
        }
    }

    // load textures
    rScene.m_gltfRayTextures.reserve(rScene.m_gltf.textures.size());
    for (auto const &tex : rScene.m_gltf.textures)
    {
        rScene.m_gltfRayTextures.push_back(LoadTextureFromImage(rScene.m_gltfRayImage[tex.source]));
    }

    // load materials
    rScene.m_gltfRayMaterials.reserve(rScene.m_gltf.materials.size());
    for (auto const &mat : rScene.m_gltf.materials)
    {
        Material &rMat = rScene.m_gltfRayMaterials.emplace_back(LoadMaterialDefault());
        if (auto it = mat.values.find("baseColorTexture"); it != mat.values.end())
        {
            int texId = it->second.json_double_value.at("index");
            SetMaterialTexture(&rMat, MATERIAL_MAP_DIFFUSE, rScene.m_gltfRayTextures[texId]);
        }
    }

    rScene.m_shaderLit = LoadShader(TextFormat("resources/shaders/glsl%i/base_lighting.vs", GLSL_VERSION),
                                   TextFormat("resources/shaders/glsl%i/lighting.fs", GLSL_VERSION));
    rScene.m_shaderLit.locs[SHADER_LOC_VECTOR_VIEW] = GetShaderLocation(rScene.m_shaderLit, "viewPos");
    int ambientLoc = GetShaderLocation(rScene.m_shaderLit, "ambient");
    glm::vec4 wo{ 7.5f, 7.0f, 9.0f, 10.0f };
    SetShaderValue(rScene.m_shaderLit, ambientLoc, &wo[0], SHADER_UNIFORM_VEC4);

    {
    Light& light = rScene.m_lights[0];
    light.type = LIGHT_DIRECTIONAL;
    light.color = {130, 110, 65, 255};
    light.position = {1.0f, 1.0f, 1.0f};
    light.target = {0.0f, 0.0f, 0.0f};
    light.enabled = true;
    light.enabledLoc    = GetShaderLocation(rScene.m_shaderLit, "lights[0].enabled");
    light.typeLoc       = GetShaderLocation(rScene.m_shaderLit, "lights[0].type");
    light.posLoc        = GetShaderLocation(rScene.m_shaderLit, "lights[0].position");
    light.targetLoc     = GetShaderLocation(rScene.m_shaderLit, "lights[0].target");
    light.colorLoc      = GetShaderLocation(rScene.m_shaderLit, "lights[0].color");
    }

    {
    Light& light = rScene.m_lights[1];
    light.type = LIGHT_DIRECTIONAL;
    light.color = {70, 50, 20, 255};
    light.position = {0.0f, 0.0f, -1.0f};
    light.target = {0.0f, 0.0f, 0.0f};
    light.enabled = true;
    light.enabledLoc    = GetShaderLocation(rScene.m_shaderLit, "lights[1].enabled");
    light.typeLoc       = GetShaderLocation(rScene.m_shaderLit, "lights[1].type");
    light.posLoc        = GetShaderLocation(rScene.m_shaderLit, "lights[1].position");
    light.targetLoc     = GetShaderLocation(rScene.m_shaderLit, "lights[1].target");
    light.colorLoc      = GetShaderLocation(rScene.m_shaderLit, "lights[1].color");
    }

    {
    Light& light = rScene.m_lights[2];
    light.type = LIGHT_DIRECTIONAL;
    light.color = {10, 16, 30, 255};
    light.position = {0.0f, -1.0f, 0.0f};
    light.target = {0.0f, 0.0f, 0.0f};
    light.enabled = true;
    light.enabledLoc    = GetShaderLocation(rScene.m_shaderLit, "lights[2].enabled");
    light.typeLoc       = GetShaderLocation(rScene.m_shaderLit, "lights[2].type");
    light.posLoc        = GetShaderLocation(rScene.m_shaderLit, "lights[2].position");
    light.targetLoc     = GetShaderLocation(rScene.m_shaderLit, "lights[2].target");
    light.colorLoc      = GetShaderLocation(rScene.m_shaderLit, "lights[2].color");
    }

    rScene.m_salads.reserve(10);
    metal_pipe(rScene.m_gltf, 0, rScene.m_characters, rScene.m_salads, rScene.m_gltfRayMaterials, rScene.m_frogs);


    for (Material &rMat : rScene.m_gltfRayMaterials)
    {
        rMat.shader = rScene.m_shaderLit;
    }

    rScene.m_cube = GenMeshCube(0.04545f, 0.1f, 0.01136f);
    rScene.m_mat = LoadMaterialDefault();


    rScene.m_ui = LoadRenderTexture(GetScreenWidth(), GetScreenHeight());


    rScene.m_camera.target = Vector3{ 0.0f, 6.9f, 0.0f };
    rScene.m_camera.position = Vector3{ 0.0f, 6.9f, 2.0f };
    rScene.m_camera.up = Vector3{ 0.0f, 1.0f, 0.0f };
    rScene.m_camera.fovy = 40.0f;
    rScene.m_camera.projection = CAMERA_PERSPECTIVE;

    // setup default pose

    for (frog_id_t id = 0; id < rScene.m_frogs.m_ids.capacity(); id ++)
    {
        if (!rScene.m_frogs.m_ids.exists(id))
        {
            continue;
        }
        rScene.m_frogs.m_tf[id][3] += glm::vec4(0.0f, 6.9f, 0.0f, 0.0f);
    }

    ToolGrab::Grab &rGrab = rScene.m_toolGrab.m_grabs.emplace_back();
    frog_id_t const frog = rScene.m_characters.begin()->second.m_frogTailTip;

    glm::vec3 const down{0.0f, -1.0f, 0.0f};
    glm::vec3 const fwd{0.0f, 0.0f, 1.0f};
    glm::vec3 const hit = rScene.m_frogs.m_tf[frog][3];
    rGrab.baitId = rScene.m_frogs.m_baitNextId++;

    FrogDyn::Bait &rBait = rScene.m_frogs.m_baits.emplace_back(
                FrogDyn::Bait{ {-1,  hit, down, fwd},
                               {frog, {0.0f, 0.0f, 0.0f}, down, fwd} });
    rBait.m_id = rGrab.baitId;
    if (g_limits)
    {
        rBait.m_forceLim = 15.2f * 10.0f * 4.0f;
    }

    // Start!

    float multithread = true;

    if (multithread)
    {
        rScene.m_running = true;
        rScene.m_updater = std::thread(&update_scene_loop, std::ref(rScene), std::ref(rGame));

        return [pScene = std::move(pScene)] (GameState &rGame) -> void
        {
            draw_scene(*pScene, rGame, true);
        };
    }
    else
    {
        rScene.m_running = true;

        return [pScene = std::move(pScene)] (GameState &rGame) -> void
        {
            update_scene(*pScene, rGame);

            draw_scene(*pScene, rGame, false);
        };
    }

}

} // namespace orni
