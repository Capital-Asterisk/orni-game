#include "scenes.hpp"

#include <longeron/containers/intarray_multimap.hpp>
#include <longeron/id_management/registry.hpp>

#include <raylib.h>
#include <raymath.h>
#include <rlgl.h>

#include <glm/mat4x4.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtx/norm.hpp>
#include <glm/gtx/matrix_interpolation.hpp>

#include <memory>
#include <iostream>


using namespace orni;




using frog_id_t = int;

struct FrogDyn
{
    struct Ball
    {
        glm::vec3   m_pos;
        float       m_radius;
        int         m_bitType;
        int         m_bitCollide;
    };

    struct BallOut
    {
        glm::vec3   m_pos;
        float       m_radius;
    };

    struct Bait
    {
        struct Insect
        {
            frog_id_t   m_id;
            glm::vec3   m_offset;
        };

        Insect m_a, m_b;
    };

    struct AABB
    {
        glm::vec3 m_min, m_max;
    };

    struct CollisionCheck
    {
        frog_id_t m_a, m_b;
    };

    struct LinAng
    {
        glm::vec3 m_lin;
        glm::vec3 m_ang;
    };

    lgrn::IdRegistry<frog_id_t> m_ids;
    std::vector<glm::mat4x4>    m_tf;
    std::vector<LinAng>         m_extImp;
    std::vector<LinAng>         m_cstImp;
    std::vector<LinAng>         m_vel;
    std::vector<float>          m_mass;
    std::vector<float>          m_scale;

    std::vector<AABB> m_aabb;

    std::vector<frog_id_t> m_canCollide;
    lgrn::IntArrayMultiMap<frog_id_t, Ball>     m_balls;
    lgrn::IntArrayMultiMap<frog_id_t, BallOut>  m_ballPos;

    //std::vector<Ball> m_ballData;

    std::vector<Bait> m_baits;
};

void apply_baits(FrogDyn &rDyn)
{
    // calculate baits
    for (int i = 0; i < rDyn.m_baits.size(); i ++)
    {
        auto const &rBait = rDyn.m_baits[i];
        glm::vec3 posA, velA, posB, velB, offsetRotA, offsetRotB;
        float massA, massB;

        if (rBait.m_a.m_id == -1)
        {
            // world-anchored
            offsetRotA = rBait.m_a.m_offset;
            posA = rBait.m_a.m_offset;
            velA = glm::vec3{0.0f};
            massA = 9999999.0f;
        }
        else
        {
            offsetRotA = rDyn.m_tf[rBait.m_a.m_id] * glm::vec4(rBait.m_a.m_offset, 0.0f);
            posA = offsetRotA + glm::vec3(rDyn.m_tf[rBait.m_a.m_id][3]);

            glm::vec3 const velTan = glm::cross(rDyn.m_vel[rBait.m_a.m_id].m_ang, offsetRotA);
            velA = rDyn.m_vel[rBait.m_a.m_id].m_lin + velTan;
            massA = rDyn.m_mass[rBait.m_a.m_id];
        }

        offsetRotB = (rDyn.m_tf[rBait.m_b.m_id] * glm::vec4(rBait.m_b.m_offset, 0.0f));
        posB = offsetRotB + glm::vec3(rDyn.m_tf[rBait.m_b.m_id][3]);

        glm::vec3 const velTan = glm::cross(rDyn.m_vel[rBait.m_b.m_id].m_ang, offsetRotB);
        velB = rDyn.m_vel[rBait.m_b.m_id].m_lin + velTan;

        massB = rDyn.m_mass[rBait.m_b.m_id];



        glm::vec3 const     posRel  = posB - posA;
        glm::vec3 const     velRel  = velB - velA;

        // PD control lol what
        glm::vec3 const     p   = posRel * 32.0f;
        glm::vec3 const     d   = velRel * 0.5f;

        float const lenBSq = glm::length2(rBait.m_b.m_offset);
        float const lenASq = glm::length2(rBait.m_a.m_offset);

        float const minMass = glm::min(massA, massB);
        float const maxLen  = glm::max(lenASq, lenBSq);

        glm::vec3 const totalLinImp = (p + d) * minMass;
        glm::vec3 const totalAngImp =  (p + d) * minMass * 1.0f / (maxLen + 1.0f);

        if (rBait.m_a.m_id != -1)
        {
            rDyn.m_cstImp[rBait.m_a.m_id].m_lin += totalLinImp;
            rDyn.m_cstImp[rBait.m_a.m_id].m_ang += glm::cross(offsetRotA, totalAngImp);

        }

        rDyn.m_cstImp[rBait.m_b.m_id].m_lin -= totalLinImp;
        rDyn.m_cstImp[rBait.m_b.m_id].m_ang -= glm::cross(offsetRotB, totalAngImp);
    }
}

void apply_ext_forces(FrogDyn &rDyn, float delta)
{
    // apply external forces
    for (frog_id_t id = 0; id < rDyn.m_ids.capacity(); id ++)
    {
        if (!rDyn.m_ids.exists(id))
        {
            continue;
        }

        rDyn.m_vel[id].m_lin += (rDyn.m_extImp[id].m_lin) / rDyn.m_mass[id];
        rDyn.m_vel[id].m_ang += (rDyn.m_extImp[id].m_ang) / rDyn.m_mass[id];

        rDyn.m_tf[id][3] += glm::vec4(rDyn.m_vel[id].m_lin, 0.0f) * delta;
        auto const& w = rDyn.m_vel[id].m_ang;

        float const mag = glm::length(w);
        if (mag > 0.001f)
        {
            auto const saved = rDyn.m_tf[id][3];
            rDyn.m_tf[id][3] = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
            rDyn.m_tf[id] = glm::axisAngleMatrix(w / mag, mag * delta) * rDyn.m_tf[id];
            rDyn.m_tf[id][3] = saved;
        }

        rDyn.m_extImp[id].m_ang = glm::vec3(0.0f);
        rDyn.m_extImp[id].m_lin = glm::vec3(0.0f);
    }
}


void apply_cst_forces(FrogDyn &rDyn, float delta)
{
    // apply external forces
    for (frog_id_t id = 0; id < rDyn.m_ids.capacity(); id ++)
    {
        if (!rDyn.m_ids.exists(id))
        {
            continue;
        }

        rDyn.m_vel[id].m_lin += (rDyn.m_cstImp[id].m_lin) / rDyn.m_mass[id];
        rDyn.m_vel[id].m_ang += (rDyn.m_cstImp[id].m_ang) / rDyn.m_mass[id];

        rDyn.m_tf[id][3] += glm::vec4(rDyn.m_vel[id].m_lin, 0.0f) * delta;
        auto const& w = rDyn.m_vel[id].m_ang;

        float const mag = glm::length(w);
        if (mag > 0.001f)
        {
            auto const saved = rDyn.m_tf[id][3];
            rDyn.m_tf[id][3] = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
            rDyn.m_tf[id] = glm::axisAngleMatrix(w / mag, mag * delta) * rDyn.m_tf[id];
            rDyn.m_tf[id][3] = saved;
        }

        rDyn.m_cstImp[id].m_ang = glm::vec3(0.0f);
        rDyn.m_cstImp[id].m_lin = glm::vec3(0.0f);
    }
}

//void clear_ext_forces(FrogDyn &rDyn)
//{
//    for (frog_id_t id = 0; id < rDyn.m_ids.capacity(); id ++)
//    {
//        if (!rDyn.m_ids.exists(id))
//        {
//            continue;
//        }

//        rDyn.m_extImp[id].m_ang = glm::vec3(0.0f);
//        rDyn.m_extImp[id].m_lin = glm::vec3(0.0f);
//    }
//}

void calc_balls_pos(FrogDyn &rDyn)
{
    for (frog_id_t id : rDyn.m_canCollide)
    {
        auto const& tf      = rDyn.m_tf[id];
        auto &rAABB         = rDyn.m_aabb[id];
        auto const ballsIn  = rDyn.m_balls[id];
        auto posOut         = rDyn.m_ballPos[id];

        rAABB.m_min = glm::vec3{9999, 9999, 9999};
        rAABB.m_max = glm::vec3{-9999, -9999, -9999};

        for (unsigned int i = 0; i < ballsIn.size(); i ++)
        {
            auto const &ball = ballsIn[i];
            auto &rOut = posOut[i];
            rOut.m_pos = tf * glm::vec4(ball.m_pos, 1.0f);
            rOut.m_radius = ball.m_radius * rDyn.m_scale[id];

            rAABB.m_max = glm::max(rOut.m_pos + rOut.m_radius, rAABB.m_max);
            rAABB.m_min = glm::min(rOut.m_pos - rOut.m_radius, rAABB.m_min);
        }
    }
}

constexpr bool aabb_intersect(FrogDyn::AABB const& a, FrogDyn::AABB const& b) noexcept
{
    return     (a.m_min.x < b.m_max.x) && (a.m_max.x > b.m_min.x)
            && (a.m_min.y < b.m_max.y) && (a.m_max.y > b.m_min.y)
            && (a.m_min.z < b.m_max.z) && (a.m_max.z > b.m_min.z);
}

using BallCollisions_t = lgrn::IntArrayMultiMap<frog_id_t, FrogDyn::BallOut>::Span;

struct BallContact
{
    glm::vec3   m_pos{0.0f};
    glm::vec3   m_nrm{0.0f};
    int         m_count{0};
};

BallContact calc_ball_collisions(BallCollisions_t a, BallCollisions_t b)
{
    BallContact contact;
    // test each ball A against each ballB
    for (auto const& ballA : a)
    {
        for (auto const& ballB : b)
        {
            float const radC = ballA.m_radius + ballB.m_radius;
            float const distSq = glm::distance2(ballA.m_pos, ballB.m_pos);
            if (distSq < (radC * radC))
            {
                contact.m_count ++;
                contact.m_pos += (ballA.m_pos * ballB.m_radius + ballB.m_pos * ballA.m_radius) / radC;
                contact.m_nrm += (ballB.m_pos - ballA.m_pos) / glm::sqrt(distSq);
            }
        }
    }
    return contact;
}

void calc_frog_collisions(FrogDyn &rDyn)
{
    //   0  1  2  3  4  5
    // 0
    // 1 x
    // 2 x  x
    // 3 x  x  x
    // 4 x  x  x  x
    // 5 x  x  x  x  x

    for (int i = 1; i < rDyn.m_canCollide.size(); i ++)
    {
        for (int j = 0; j < i; j ++)
        {
            frog_id_t const a = rDyn.m_canCollide[i];
            frog_id_t const b = rDyn.m_canCollide[j];

            if (aabb_intersect(rDyn.m_aabb[a], rDyn.m_aabb[b]))
            {
                // AABBs intersect, now start checking balls
                BallContact contact = calc_ball_collisions(rDyn.m_ballPos[a], rDyn.m_ballPos[b]);
                if (contact.m_count != 0)
                {
                    glm::vec3 const pos = contact.m_pos / float(contact.m_count);
                    glm::vec3 const nrm = contact.m_nrm / float(contact.m_count);

                    glm::vec3 const posnrm = pos + nrm;

                    DrawSphereWires(reinterpret_cast<Vector3 const&>(pos), 0.1f, 5, 6, Color{255, 0, 255, 255});
                    DrawLine3D(reinterpret_cast<Vector3 const&>(pos), reinterpret_cast<Vector3 const&>(posnrm), Color{0, 255, 0, 255});
                }
            }
        }
    }
}

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

    rScene.m_frogs.m_baits[0].m_b.m_offset.y += (float(IsKeyDown(KEY_S)) - float(IsKeyDown(KEY_W))) * delta * 1.0f;
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

                apply_baits(rScene.m_frogs);

                apply_cst_forces(rScene.m_frogs, delta * cstPercent / float(rScene.m_cstSteps));

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

                auto tf =(rScene.m_frogs.m_tf[id]);

                auto tip = glm::vec3(tf[3] - tf[1] * 0.5f);
                auto tail = glm::vec3(tf[3] + tf[1] * 0.5f);

                DrawCylinderWiresEx(
                        reinterpret_cast<Vector3&>(tip),
                        reinterpret_cast<Vector3&>(tail),
                        0.2f, 0.2f, 4, Color{255, 255, 255, 255});
                //DrawMesh(rScene.m_cube, rScene.m_mat, reinterpret_cast<Matrix&>(transposed));

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



frog_id_t add_frog(FrogDyn &rDyn, glm::mat4x4 const& tf, float mass)
{
    frog_id_t const id = rDyn.m_ids.create();
    std::size_t const capacity = rDyn.m_ids.capacity();
    rDyn.m_cstImp   .resize(capacity);
    rDyn.m_extImp   .resize(capacity);
    rDyn.m_vel      .resize(capacity);
    rDyn.m_mass     .resize(capacity);
    rDyn.m_tf       .resize(capacity);
    rDyn.m_scale    .resize(capacity);
    rDyn.m_aabb     .resize(capacity);
    rDyn.m_balls    .resize_ids(capacity);

    rDyn.m_scale[id]    = 1.0f;
    rDyn.m_mass[id]     = mass;
    rDyn.m_tf[id]       = tf;

    return id;
}


SceneFunc_t orni::gen_test_scene_b()
{
    std::shared_ptr<TestSceneB> pScene = std::make_shared<TestSceneB>();
    TestSceneB &rScene = *pScene;

    rScene.m_extPercent = 0.7f;
    rScene.m_cstSteps = 8;

    rScene.m_mat = LoadMaterialDefault();
    rScene.m_cube = GenMeshCube(0.5f, 0.5f, 0.5f);
    rScene.m_ui = LoadRenderTexture(GetScreenWidth(), GetScreenHeight());

    add_frog(rScene.m_frogs, glm::translate(glm::mat4x4{1.0f}, {0.0f, 1.0f, 0.0f}), 0.1f );

    add_frog(rScene.m_frogs, glm::translate(glm::mat4x4{1.0f}, {0.0f, 1.0f, 0.0f}), 0.1f );

    add_frog(rScene.m_frogs, glm::translate(glm::mat4x4{1.0f}, {0.0f, 1.0f, 0.0f}), 0.1f );

    add_frog(rScene.m_frogs, glm::translate(glm::mat4x4{1.0f}, {0.0f, 1.0f, 0.0f}), 1.5f );

    // global constrained
    rScene.m_frogs.m_baits.emplace_back(
                FrogDyn::Bait{ {-1,  {0.0f, 7.0f, 0.0f}},
                               {0,   {0.0f, 1.0f, 0.0f}} });

    rScene.m_frogs.m_baits.emplace_back(
                FrogDyn::Bait{ {0,  {0.0f, -1.0f, 0.0f}},
                               {1,   {0.0f, 1.0f, 0.0f}} });

    rScene.m_frogs.m_baits.emplace_back(
                FrogDyn::Bait{ {1,   {0.0f, -1.0f, 0.0f}},
                               {2,   {0.0f,  1.0f, 0.0f}} });

    rScene.m_frogs.m_baits.emplace_back(
                FrogDyn::Bait{ {2,   {0.0f, -1.0f, 0.0f}},
                               {3,   {0.0f,  1.0f, 0.0f}} });

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

    rScene.m_cstSteps = 2;

    rScene.m_mat = LoadMaterialDefault();
    rScene.m_cube = GenMeshCube(0.5f, 0.5f, 0.5f);
    rScene.m_ui = LoadRenderTexture(GetScreenWidth(), GetScreenHeight());

    add_frog(rScene.m_frogs, glm::translate(glm::mat4x4{1.0f}, {0.0f, 1.0f, 0.0f}), 6.0f );

    add_frog(rScene.m_frogs, glm::translate(glm::mat4x4{1.0f}, {0.0f, 1.0f, 0.0f}), 0.2f );

    add_frog(rScene.m_frogs, glm::translate(glm::mat4x4{1.0f}, {0.0f, 1.0f, 0.0f}), 5.0f );

    add_frog(rScene.m_frogs, glm::translate(glm::mat4x4{1.0f}, {0.0f, 1.0f, 0.0f}), 4.0f );

    // add balls
    rScene.m_frogs.m_balls.resize_ids(20);
    rScene.m_frogs.m_balls.resize_data(20);
    rScene.m_frogs.m_ballPos.resize_ids(20);
    rScene.m_frogs.m_ballPos.resize_data(20);

    rScene.m_frogs.m_canCollide = {0, 2};

    rScene.m_frogs.m_balls.emplace(0, { {{0.0f, 0.0f, 0.0f}, 0.4f, 1, 1}, {{0.0f, 0.5f, 0.0f}, 0.4f, 1, 1}, {{-1.0f, 3.0f, 0.0f}, 0.4f, 1, 1} } );
    rScene.m_frogs.m_ballPos.emplace(0, 3);

    rScene.m_frogs.m_balls.emplace(2, { {{0.0f, 0.0f, 0.0f}, 0.7f, 1, 1}, {{0.0f, -0.5f, 0.0f}, 0.7f, 1, 1} });
    rScene.m_frogs.m_ballPos.emplace(2, 2);


    rScene.m_frogs.m_baits.emplace_back(
                FrogDyn::Bait{ {-1,  {1.0f, 5.0f, 0.0f}},
                               {0,   {0.0f, 2.0f, 0.0f}} });

    rScene.m_frogs.m_baits.emplace_back(
                FrogDyn::Bait{ {0,  {0.0f, -2.0f, 0.0f}},
                               {1,   {0.0f, 2.0f, 0.0f}} });

    rScene.m_frogs.m_baits.emplace_back(
                FrogDyn::Bait{ {-1,  {-1.0f, 5.0f, 0.0f}},
                               {2,   {0.0f,  2.0f, 0.0f}} });

    rScene.m_frogs.m_baits.emplace_back(
                FrogDyn::Bait{ {2,   {0.0f, -2.0f, 0.0f}},
                               {3,   {0.0f,  2.0f, 0.0f}} });

    rScene.m_camera.target = Vector3{ 0.0f, 2.0f, 0.0f };
    rScene.m_camera.up = Vector3{ 0.0f, 1.0f, 0.0f };
    rScene.m_camera.fovy = 50.0f;
    rScene.m_camera.projection = CAMERA_PERSPECTIVE;

    return [pScene = std::move(pScene)] (GameState &rGame) -> void
    {
        draw_scene(*pScene);
    };
}

