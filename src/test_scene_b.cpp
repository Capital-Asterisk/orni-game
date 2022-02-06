#include "scenes.hpp"

#include <raylib.h>
#include <raymath.h>
#include <rlgl.h>

#include <glm/mat4x4.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtx/matrix_interpolation.hpp>

#include <memory>
#include <iostream>


using namespace orni;

struct TestSceneB
{
    Camera3D m_camera;

    Material m_mat;

    RenderTexture2D m_ui;

    Mesh m_cube;

    struct Insect
    {
        int m_id;
        glm::vec3 m_offset;
    };

    struct Bait
    {
        Insect m_a, m_b;
    };

    std::vector<glm::mat4x4> m_frog;
    std::vector<glm::vec3> m_extLinImp;
    std::vector<glm::vec3> m_extAngImp;
    std::vector<glm::vec3> m_cstLinImp;
    std::vector<glm::vec3> m_cstAngImp;
    std::vector<glm::vec3> m_linVel;
    std::vector<glm::vec3> m_angVel;
    std::vector<float> m_mass;

    std::vector<Bait> m_baits;


    float m_time = 0.0f;

    float m_camAngle = 0.0f;
};

void apply_impulse_at()
{

}

static void draw_scene(TestSceneB &rScene)
{
    auto &t = rScene.m_time;

    float delta = 1.0f / 60.0f;
    t += GetFrameTime();

    rScene.m_baits[0].m_b.m_offset.y += (float(IsKeyDown(KEY_S)) - float(IsKeyDown(KEY_W))) * delta * 5.0f;

    std::cout << "a: " << rScene.m_baits[0].m_b.m_offset.y << "\n";

    rScene.m_camAngle += (float(IsKeyDown(KEY_RIGHT)) - float(IsKeyDown(KEY_LEFT))) * delta * 3.14159f;

    rScene.m_camera.position = Vector3{ std::sin(rScene.m_camAngle) * 16.0f, 3.0f, std::cos(rScene.m_camAngle) * 16.0f };



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
                for (int i = 0; i < rScene.m_frog.size(); i ++)
                {
                    rScene.m_linVel[i] *= 0.50f;
                }
            }

            for (int i = 0; i < rScene.m_frog.size(); i ++)
            {
                // apply gravity
                rScene.m_extLinImp[i].y -= 9.81f * delta * rScene.m_mass[i];

                rScene.m_extLinImp[i] += push;

                auto p = glm::vec3(rScene.m_frog[i][3]);
                DrawLine3D(Vector3{p.x, p.y, p.z}, Vector3{p.x + push.x * 5.0f, p.y + push.y * 5.0f, p.z + push.z * 5.0f}, Color{0, 255, 0, 255});
            }

            float deltaA = delta / 2.0f;

            // apply external forces
            for (int i = 0; i < rScene.m_frog.size(); i ++)
            {
                rScene.m_linVel[i] += (rScene.m_extLinImp[i]) / rScene.m_mass[i];
                rScene.m_angVel[i] += (rScene.m_extAngImp[i]) / rScene.m_mass[i];

                rScene.m_frog[i][3] += glm::vec4(rScene.m_linVel[i], 0.0f) * deltaA;
                auto const& w = rScene.m_angVel[i];

                float const mag = glm::length(w);
                if (mag > 0.001f)
                {
                    auto const saved = rScene.m_frog[i][3];
                    rScene.m_frog[i][3] = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
                    rScene.m_frog[i] = glm::axisAngleMatrix(w / mag, mag * deltaA) * rScene.m_frog[i];
                    rScene.m_frog[i][3] = saved;
                }
            }

            float deltaB = delta / 8.0f;

            // repeat constrain forces a few times
            for (int j = 0; j < 8; j++)
            {

                // calculate baits
                for (int i = 0; i < rScene.m_baits.size(); i ++)
                {
                    auto const &rBait = rScene.m_baits[i];
                    glm::vec3 posA, velA, posB, velB;
                    float massA, massB;

                    if (rBait.m_a.m_id == -1)
                    {
                        // world-anchored
                        posA = rBait.m_a.m_offset;
                        velA = glm::vec3{0.0f};
                        massA = 9999999.0f;
                    }
                    else
                    {
                        glm::vec3 const offsetRotated(rScene.m_frog[rBait.m_a.m_id] * glm::vec4(rBait.m_a.m_offset, 0.0f));
                        posA = offsetRotated + glm::vec3(rScene.m_frog[rBait.m_a.m_id][3]);

                        glm::vec3 const velTan = glm::cross(rScene.m_angVel[rBait.m_a.m_id], offsetRotated);
                        velA = rScene.m_linVel[rBait.m_a.m_id] + velTan;
                        massA = rScene.m_mass[rBait.m_a.m_id];
                    }

                    glm::vec3 const offsetRotated(rScene.m_frog[rBait.m_b.m_id] * glm::vec4(rBait.m_b.m_offset, 0.0f));
                    posB = offsetRotated + glm::vec3(rScene.m_frog[rBait.m_b.m_id][3]);

                    glm::vec3 const velTan = glm::cross(rScene.m_angVel[rBait.m_b.m_id], offsetRotated);
                    velB = rScene.m_linVel[rBait.m_b.m_id] + velTan;

                    massB = rScene.m_mass[rBait.m_b.m_id];

                    float const minMass = glm::min(massA, massB);

                    glm::vec3 const     posRel  = posB - posA;
                    glm::vec3 const     velRel  = velB - velA;
                    float const         r       = glm::max(glm::length(posRel), 0.001f);
                    glm::vec3 const     dir     = posRel / r;

                    float const         velDot  = glm::dot(velRel, dir);

                    // PD control lol what
                    float const         p   = r * (16.0f + minMass * 2.0f);
                    float const         d   = velDot;
                    glm::vec3 const     pd  = -dir * p;

                    float const lenB = glm::length(rBait.m_b.m_offset);
                    glm::vec3 const     damp = -velRel * 0.25f;

                    glm::vec3 const totalLinImp = (pd + damp) * minMass;
                    glm::vec3 const totalAngImp = (glm::cross(offsetRotated, totalLinImp)) * minMass;

                    if (rBait.m_a.m_id != -1)
                    {
                        rScene.m_cstLinImp[rBait.m_a.m_id] -= totalLinImp;
                        rScene.m_cstAngImp[rBait.m_a.m_id] -= totalAngImp / (lenB+lenB*lenB);

                    }

                    rScene.m_cstLinImp[rBait.m_b.m_id] += totalLinImp;
                    rScene.m_cstAngImp[rBait.m_b.m_id] += totalAngImp / (lenB+lenB*lenB);

                }

                // apply constrain forces
                for (int i = 0; i < rScene.m_frog.size(); i ++)
                {

                    rScene.m_linVel[i] += (rScene.m_extLinImp[i] + rScene.m_cstLinImp[i]) / rScene.m_mass[i];
                    rScene.m_angVel[i] += (rScene.m_extAngImp[i] + rScene.m_cstAngImp[i]) / rScene.m_mass[i];

                    rScene.m_frog[i][3] += glm::vec4(rScene.m_linVel[i], 0.0f) * deltaB;
                    auto const& w = rScene.m_angVel[i];

                    float const mag = glm::length(w);
                    if (mag > 0.001f)
                    {
                        auto const saved = rScene.m_frog[i][3];
                        rScene.m_frog[i][3] = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
                        rScene.m_frog[i] = glm::axisAngleMatrix(w / mag, mag * deltaB) * rScene.m_frog[i];
                        rScene.m_frog[i][3] = saved;
                    }

                    rScene.m_cstLinImp[i] = glm::vec3(0.0f);
                    rScene.m_cstAngImp[i] = glm::vec3(0.0f);


                }
            }

            // clear external forces
            for (int i = 0; i < rScene.m_frog.size(); i ++)
            {

                rScene.m_extLinImp[i] = glm::vec3(0.0f);
                rScene.m_extAngImp[i] = glm::vec3(0.0f);
            }

            glm::vec3 avgPos{0.0f};
            float totalMass = 0.0f;

            // draw frogs
            for (int i = 0; i < rScene.m_frog.size(); i ++)
            {
                avgPos += glm::vec3(rScene.m_frog[i][3]) * rScene.m_mass[i];
                totalMass += rScene.m_mass[i];

                auto transposed = glm::transpose(rScene.m_frog[i]);

                DrawMesh(rScene.m_cube, rScene.m_mat, reinterpret_cast<Matrix&>(transposed));
            }

            avgPos /= totalMass;

            // draw center of mass
            DrawSphereWires(reinterpret_cast<Vector3&>(avgPos), 0.3f, 2, 4, Color{255, 255, 0, 255});

            // draw baits
            for (int i = 0; i < rScene.m_baits.size(); i ++)
            {
                auto const &rBait = rScene.m_baits[i];
                glm::vec3 posA, velA, posB, velB;

                if (rBait.m_a.m_id == -1)
                {
                    // world-anchored
                    posA = rBait.m_a.m_offset;
                    velA = glm::vec3{0.0f};
                }
                else
                {
                    glm::vec3 const offsetRotated(rScene.m_frog[rBait.m_a.m_id] * glm::vec4(rBait.m_a.m_offset, 0.0f));
                    posA = offsetRotated + glm::vec3(rScene.m_frog[rBait.m_a.m_id][3]);
                }

                glm::vec3 const offsetRotated(rScene.m_frog[rBait.m_b.m_id] * glm::vec4(rBait.m_b.m_offset, 0.0f));
                posB = offsetRotated + glm::vec3(rScene.m_frog[rBait.m_b.m_id][3]);

                DrawSphereWires(reinterpret_cast<Vector3&>(posA), 0.1f, 2, 4, Color{255, 0, 0, 255});
                DrawSphere(reinterpret_cast<Vector3&>(posB), 0.05f, Color{0, 0, 255, 255});


                if (rBait.m_a.m_id != -1)
                {
                    Vector3 origA = reinterpret_cast<Vector3&>(rScene.m_frog[rBait.m_a.m_id][3]);
                    DrawLine3D(origA, reinterpret_cast<Vector3&>(posA), Color{255, 0, 0, 255});
                }
                Vector3 origB = reinterpret_cast<Vector3&>(rScene.m_frog[rBait.m_b.m_id][3]);
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

    rScene.m_mat = LoadMaterialDefault();
    rScene.m_cube = GenMeshCube(0.5f, 0.5f, 0.5f);
    rScene.m_ui = LoadRenderTexture(GetScreenWidth(), GetScreenHeight());

    rScene.m_frog       .resize(4);
    rScene.m_extLinImp  .resize(4);
    rScene.m_extAngImp  .resize(4);
    rScene.m_cstLinImp  .resize(4);
    rScene.m_cstAngImp  .resize(4);
    rScene.m_linVel     .resize(4);
    rScene.m_angVel     .resize(4);
    rScene.m_mass       .resize(4);

    rScene.m_frog[0] = glm::mat4x4(1.0f);
    rScene.m_frog[0][3].y = 1.0f;
    rScene.m_mass[0] = 4.0f;

    rScene.m_frog[1] = glm::mat4x4(1.0f);
    rScene.m_frog[1][3].y = 1.0f;
    rScene.m_mass[1] = 5.0f;

    rScene.m_frog[2] = glm::mat4x4(1.0f);
    rScene.m_frog[2][3].y = 1.0f;
    rScene.m_mass[2] = 0.1f;

    rScene.m_frog[3] = glm::mat4x4(1.0f);
    rScene.m_frog[3][3].y = 1.0f;
    rScene.m_mass[3] = 0.1f;

    rScene.m_baits.resize(4);

    // global constrained
    rScene.m_baits[0].m_a.m_id = -1;
    rScene.m_baits[0].m_a.m_offset.y = 7;
    rScene.m_baits[0].m_b.m_id = 0;
    rScene.m_baits[0].m_b.m_offset.y = 3;

    // floating
//    rScene.m_baits[0].m_a.m_id = 3;
//    rScene.m_baits[0].m_a.m_offset.y = -3;
//    rScene.m_baits[0].m_b.m_id = 0;
//    rScene.m_baits[0].m_b.m_offset.y = 3;

    rScene.m_baits[1].m_a.m_id = 0;
    rScene.m_baits[1].m_a.m_offset.y = -1;
    rScene.m_baits[1].m_b.m_id = 1;
    rScene.m_baits[1].m_b.m_offset.y = 1;

    rScene.m_baits[2].m_a.m_id = 1;
    rScene.m_baits[2].m_a.m_offset.y = -0.1;
    rScene.m_baits[2].m_b.m_id = 2;
    rScene.m_baits[2].m_b.m_offset.y = 0.1;

    rScene.m_baits[3].m_a.m_id = 2;
    rScene.m_baits[3].m_a.m_offset.y = -0.1;
    rScene.m_baits[3].m_b.m_id = 3;
    rScene.m_baits[3].m_b.m_offset.y = 0.1;

    rScene.m_camera.target = Vector3{ 0.0f, 2.0f, 0.0f };
    rScene.m_camera.up = Vector3{ 0.0f, 1.0f, 0.0f };
    rScene.m_camera.fovy = 50.0f;
    rScene.m_camera.projection = CAMERA_PERSPECTIVE;

    return [pScene = std::move(pScene)] (GameState &rGame) -> void
    {
        draw_scene(*pScene);
    };
}
