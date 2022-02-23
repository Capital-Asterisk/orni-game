#include "scenes.hpp"
#include "frogdyn.hpp"

#include <raylib.h>
#include <raymath.h>
#include <rlgl.h>

#include <glm/mat4x4.hpp>
#include <glm/gtx/string_cast.hpp>

#include <memory>
#include <iostream>


using namespace orni;

using namespace frogdyn;

struct TestSceneB
{
    Camera3D m_camera;

    Material m_mat;

    RenderTexture2D m_ui;

    Mesh m_cube;

    FrogDyn m_frogs;

    float m_time = 0.0f;

    float m_camDist = 16.0f;
    float m_camAngle = 0.0f;

    int m_cstSteps = 2;
    float m_extPercent = 0.5f;
};



static void draw_scene(TestSceneB &rScene)
{
    auto &t = rScene.m_time;

    float delta = 1.0f / 60.0f;
    t += GetFrameTime();

    rScene.m_frogs.m_baits[0].m_b.m_offset.y += (float(IsKeyDown(KEY_S)) - float(IsKeyDown(KEY_W))) * delta * 10.0f;
    rScene.m_camDist += (float(IsKeyDown(KEY_Z)) - float(IsKeyDown(KEY_X))) * delta * 10.0f;

    //std::cout << "a: " << rScene.m_frogs.m_baits[0].m_b.m_offset.y << "\n";

    rScene.m_camAngle += (float(IsKeyDown(KEY_RIGHT)) - float(IsKeyDown(KEY_LEFT))) * delta * 3.14159f;

    rScene.m_camera.position = Vector3{ std::sin(rScene.m_camAngle) * rScene.m_camDist, 3.0f, std::cos(rScene.m_camAngle) * rScene.m_camDist };

    BeginDrawing();

        ClearBackground(Color{ 10, 10, 10, 240 });

        BeginMode3D(rScene.m_camera);


            float pushMag = (float(IsKeyDown(KEY_DOWN)) - float(IsKeyDown(KEY_UP))) * delta * 20.0f;

            glm::vec3 const push = glm::vec3{
                        glm::sin(rScene.m_camAngle + 1.5708f),
                        0.0f,
                        glm::cos(rScene.m_camAngle + 1.5708f)} * pushMag;

            if (IsKeyDown(KEY_SPACE))
            {
                for (int i = 0; i < rScene.m_frogs.m_ids.size(); i ++)
                {
                    rScene.m_frogs.m_vel[i].m_lin *= 0.50f;
                }
            }


            if (IsKeyDown(KEY_LEFT_ALT))
            {
                for (int i = 0; i < rScene.m_frogs.m_ids.size(); i ++)
                {
                    rScene.m_frogs.m_vel[i].m_ang *= 0.50f;
                }
            }

            for (int i = 0; i < rScene.m_frogs.m_ids.size(); i ++)
            {
                // apply gravity
                rScene.m_frogs.m_extImp[i].m_lin.y -= 9.81f * delta * rScene.m_frogs.m_mass[i];

                rScene.m_frogs.m_extImp[i].m_lin += push;

                auto p = glm::vec3(rScene.m_frogs.m_tf[i][3]);
                DrawLine3D(Vector3{p.x, p.y, p.z}, Vector3{p.x + push.x * 5.0f, p.y + push.y * 5.0f, p.z + push.z * 5.0f}, Color{0, 255, 0, 255});
            }

            float const cstPercent = 1.0f - rScene.m_extPercent;

            apply_ext_forces(rScene.m_frogs, delta * rScene.m_extPercent);

            // repeat constrain forces a few times
            for (int j = 0; j < rScene.m_cstSteps; j++)
            {
                float smldelta = delta * cstPercent / float(rScene.m_cstSteps);

                apply_baits(rScene.m_frogs, {16.0f, 0.1f}, smldelta);

                apply_cst_forces(rScene.m_frogs, smldelta);

                calc_balls_pos(rScene.m_frogs);

                calc_frog_collisions(rScene.m_frogs);
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

                auto tf = rScene.m_frogs.m_tf[id];
                auto transposed = glm::transpose(tf);

                auto tip = glm::vec3(tf[3] - tf[1] * 0.5f);
                auto tail = glm::vec3(tf[3] + tf[1] * 0.5f);

//                DrawCylinderWiresEx(
//                        reinterpret_cast<Vector3&>(tip),
//                        reinterpret_cast<Vector3&>(tail),
//                        0.2f, 0.2f, 4, Color{255, 255, 255, 255});
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

            avgPos /= totalMass;

            // draw AABBs
            for (frog_id_t id : rScene.m_frogs.m_canCollide)
            {
                BoundingBox box{ reinterpret_cast<Vector3&>(rScene.m_frogs.m_aabb[id].m_min),
                                 reinterpret_cast<Vector3&>(rScene.m_frogs.m_aabb[id].m_max)};
                DrawBoundingBox(box, Color{0, 255, 0, 255});
            }

            // draw center of mass
            DrawSphereWires(reinterpret_cast<Vector3&>(avgPos), 0.3f, 2, 4, Color{255, 255, 0, 255});

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

                DrawSphereWires(reinterpret_cast<Vector3&>(posA), 0.1f, 2, 4, Color{255, 0, 0, 255});
                DrawSphere(reinterpret_cast<Vector3&>(posB), 0.05f, Color{0, 0, 255, 255});


                if (rBait.m_a.m_id != -1)
                {
                    Vector3 origA = reinterpret_cast<Vector3&>(rScene.m_frogs.m_tf[rBait.m_a.m_id][3]);
                    DrawLine3D(origA, reinterpret_cast<Vector3&>(posA), Color{255, 0, 0, 255});
                }
                Vector3 origB = reinterpret_cast<Vector3&>(rScene.m_frogs.m_tf[rBait.m_b.m_id][3]);
                DrawLine3D(origB, reinterpret_cast<Vector3&>(posB), Color{0, 0, 255, 255});

            }

            DrawGrid(10, 1.0f);

        EndMode3D();

        BeginTextureMode(rScene.m_ui);
            ClearBackground(Color{ 0, 0, 0, 0 });
            BeginMode3D(rScene.m_camera);

            EndMode3D();
        EndTextureMode();

        auto &rTex = rScene.m_ui.texture;
        DrawTextureRec(rTex, Rectangle{0.0f, 0.0f, float(rTex.width), -float(rTex.height) * 1.0f}, Vector2{0.0f, 0.0f}, Color{255, 255, 255, 255});
    EndDrawing();
}




SceneFunc_t orni::gen_test_scene_b()
{
    std::shared_ptr<TestSceneB> pScene = std::make_shared<TestSceneB>();
    TestSceneB &rScene = *pScene;

    rScene.m_extPercent = 0.7f;
    rScene.m_cstSteps = 8;

    rScene.m_mat = LoadMaterialDefault();
    rScene.m_cube = GenMeshCube(0.04545f * 10.0f, 0.1f * 10.0f, 0.01136f * 10.0f);
    rScene.m_ui = LoadRenderTexture(GetScreenWidth(), GetScreenHeight());

    add_frog(rScene.m_frogs, glm::translate(glm::mat4x4{1.0f}, {0.0f, 1.0f, 0.0f}), 0.1f );

    add_frog(rScene.m_frogs, glm::translate(glm::mat4x4{1.0f}, {0.0f, 1.0f, 0.0f}), 0.1f );

    add_frog(rScene.m_frogs, glm::translate(glm::mat4x4{1.0f}, {0.0f, 1.0f, 0.0f}), 0.1f );

    add_frog(rScene.m_frogs, glm::translate(glm::mat4x4{1.0f}, {0.0f, 1.0f, 0.0f}), 1.5f );

    glm::vec3 const down{0.0f, -1.0f, 0.0f};
    glm::vec3 const fwd{0.0f, 0.0f, 1.0f};

    // global constrained
    rScene.m_frogs.m_baits.emplace_back(
                FrogDyn::Bait{ {-1,  {0.0f, 7.0f, 0.0f}, down, fwd},
                               {0,   {0.0f, 1.0f, 0.0f}, down, fwd} });

    rScene.m_frogs.m_baits.emplace_back(
                FrogDyn::Bait{ {0,  {0.0f, 0.0f, 0.0f}, down, fwd},
                               {1,   {0.0f, 2.0f, 0.0f}, down, fwd} });

    rScene.m_frogs.m_baits.emplace_back(
                FrogDyn::Bait{ {1,   {0.0f, -1.0f, 0.0f}, down, fwd},
                               {2,   {0.0f,  1.0f, 0.0f}, down, fwd} });

    rScene.m_frogs.m_baits.emplace_back(
                FrogDyn::Bait{ {2,   {0.0f, -1.0f, 0.0f}, down, fwd},
                               {3,   {0.0f,  1.0f, 0.0f}, down, fwd} });

    rScene.m_camera.target = Vector3{ 0.0f, 2.0f, 0.0f };
    rScene.m_camera.up = Vector3{ 0.0f, 1.0f, 0.0f };
    rScene.m_camera.fovy = 50.0f;
    rScene.m_camera.projection = CAMERA_PERSPECTIVE;

    return [pScene = std::move(pScene)] (GameState &rGame) -> void
    {
        draw_scene(*pScene);
    };
}

SceneFunc_t orni::gen_test_scene_b_b()
{
    std::shared_ptr<TestSceneB> pScene = std::make_shared<TestSceneB>();
    TestSceneB &rScene = *pScene;

    rScene.m_cstSteps = 24;

    rScene.m_mat = LoadMaterialDefault();
    rScene.m_cube = GenMeshCube(0.04545f * 10.0f, 0.1f * 10.0f, 0.01136f * 10.0f);
    rScene.m_ui = LoadRenderTexture(GetScreenWidth(), GetScreenHeight());

    add_frog(rScene.m_frogs, glm::translate(glm::mat4x4{1.0f}, {1.0f, 1.0f, 0.0f}), 6.0f );

    add_frog(rScene.m_frogs, glm::translate(glm::mat4x4{1.0f}, {2.0f, 1.0f, 0.0f}), 0.2f );

    add_frog(rScene.m_frogs, glm::translate(glm::mat4x4{1.0f}, {3.0f, 1.0f, 0.0f}), 5.0f );

    add_frog(rScene.m_frogs, glm::translate(glm::mat4x4{1.0f}, {4.0f, 1.0f, 0.0f}), 4.0f );

    //rScene.m_frogs.m_vel[2].m_ang.y = 10.0f;

    // add balls
    rScene.m_frogs.m_balls.resize_ids(20);
    rScene.m_frogs.m_balls.resize_data(20);
    rScene.m_frogs.m_ballPos.resize_ids(20);
    rScene.m_frogs.m_ballPos.resize_data(20);

    rScene.m_frogs.m_canCollide = {0, 1, 2, 3};

    rScene.m_frogs.m_balls.emplace(0, { {{0.0f, 0.0f, 0.0f}, 0.4f, 1, 1}, {{0.0f, 0.5f, 0.0f}, 0.4f, 1, 1} } );
    rScene.m_frogs.m_ballPos.emplace(0, 2);

    rScene.m_frogs.m_balls.emplace(2, { {{0.0f, 0.0f, 0.0f}, 0.7f, 1, 1}, {{0.0f, -0.5f, 0.0f}, 0.7f, 1, 1} });
    rScene.m_frogs.m_ballPos.emplace(2, 2);

    rScene.m_frogs.m_balls.emplace(1, { {{0.0f, 0.0f, 0.0f}, 0.7f, 1, 1} });
    rScene.m_frogs.m_ballPos.emplace(1, 1);

    rScene.m_frogs.m_balls.emplace(3, { {{0.0f, 0.0f, 0.0f}, 0.7f, 1, 1} });
    rScene.m_frogs.m_ballPos.emplace(3, 1);

    glm::vec3 const down{0.0f, -1.0f, 0.0f};
    glm::vec3 const fwd{0.0f, 0.0f, 1.0f};

    rScene.m_frogs.m_baits.emplace_back(
                FrogDyn::Bait{ {-1,  {1.0f, 5.0f, 0.0f}, down, fwd},
                               {0,   {0.0f, 2.0f, 0.0f}, down, fwd} });

    rScene.m_frogs.m_baits.emplace_back(
                FrogDyn::Bait{ {0,  {0.0f, -2.0f, 0.0f}, down, fwd},
                               {1,   {0.0f, 2.0f, 0.0f}, down, fwd} });

    rScene.m_frogs.m_baits.emplace_back(
                FrogDyn::Bait{ {-1,  {-1.0f, 5.0f, 0.0f}, down, fwd},
                               {2,   {0.0f,  2.0f, 0.0f}, down, fwd} });

    rScene.m_frogs.m_baits.emplace_back(
                FrogDyn::Bait{ {0,   {0.0f, -2.0f, 0.0f}, down, fwd},
                               {3,   {0.0f,  2.0f, 0.0f}, down, fwd} });

    rScene.m_camera.target = Vector3{ 0.0f, 2.0f, 0.0f };
    rScene.m_camera.up = Vector3{ 0.0f, 1.0f, 0.0f };
    rScene.m_camera.fovy = 50.0f;
    rScene.m_camera.projection = CAMERA_PERSPECTIVE;

    return [pScene = std::move(pScene)] (GameState &rGame) -> void
    {
        draw_scene(*pScene);
    };
}

