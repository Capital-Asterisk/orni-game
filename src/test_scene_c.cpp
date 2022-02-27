#include "scenes.hpp"
#include "salads.hpp"

#include <raylib.h>
#include <raymath.h>
#include <rlgl.h>

#include <glm/gtx/transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/matrix_interpolation.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtx/vector_angle.hpp>

#include <iostream>
#include <memory>
#include <unordered_map>

namespace orni
{



struct TestSceneC
{
    ~TestSceneC()
    {
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
            UnloadMaterial(rMat);
        }
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
    }

    Salads_t                m_salads;
    FrogDyn                 m_frogs;

    Characters_t            m_characters;

    Inputs                  m_inputs;
    ToolGrab                m_toolGrab;

    tinygltf::Model         m_gltf;
    std::vector<Image>      m_gltfRayImage;
    std::vector<Texture>    m_gltfRayTextures;
    std::vector<Material>   m_gltfRayMaterials;

    Camera3D                m_camera;

    Mesh                    m_cube;
    Material                m_mat;

    RenderTexture2D         m_ui;

    int                     m_cstSteps = 15;
    float                   m_extPercent = 0.5f;

    float                   m_camDist{3.0f};
    float                   m_camYaw{0.0f};
    float                   m_camPitch{0.0f};
    float                   m_time{0.0f};
    bool                    m_gravity{false};


};

static void update_apples(Apples &rApples, meshdeform::Joints const& rJoints)
{
    for (apple_id_t id = 0; id < rApples.m_ids.capacity(); id ++)
    {
        if ( ! rApples.m_ids.exists(id))
        {
            continue;
        }

        auto const &apl = rApples.m_data[id];

        rApples.m_dataOut[id] = rJoints.m_nodeTf[apl.m_jointParent] * apl.m_tf;
    }
}

// lol
float const eyeDepth = 0.03f;
float const eyeRadius = 0.09f;
float const eyeFov = glm::atan(eyeRadius / eyeDepth) * 2.0f;

static glm::vec2 calc_eye_pos(glm::mat4x4 const& eyeTf, glm::vec3 tgt)
{

    glm::mat4 eyeMatrix = glm::perspective(eyeFov, 1.0f, 0.001f, 1000.0f)
                        * glm::lookAt(glm::vec3(eyeTf[3]) - glm::vec3(eyeTf[2]) * eyeDepth,
                                      glm::vec3(eyeTf[3]) + glm::vec3(eyeTf[2]),
                                      glm::vec3(eyeTf[1]));

    glm::vec3 foo = eyeMatrix * glm::vec4(tgt, 1.0f);
    foo.x /= foo.z;
    foo.y /= foo.z;

    return glm::vec2(foo);
}

static bool eye_visible(glm::mat4x4 const& eyeTf, glm::vec3 tgt)
{
    float const ang = glm::angle(glm::vec3(eyeTf[2]), tgt - glm::vec3(eyeTf[3]));
    return ang < eyeFov;
}

static void draw_iris(Texture2D texture, int i, glm::vec2 pos)
{
    constexpr float const w = 128.0f, h = 128.0f, r = 64.0f, c = 32.0f;

    float const len = glm::min(glm::length(pos) * r, c);
    pos = glm::normalize(pos) * glm::vec2{-1.0f, 1.0f} * len;

    float const x = 128 * i, y = 128;
    float const xlapP = glm::clamp(pos.x, 0.0f, w);
    float const xlapN = glm::clamp(pos.x, -w, 0.0f);
    DrawTexturePro(texture,
                   Rectangle{-xlapN,            0,          w - xlapP + xlapN,  -h},
                   Rectangle{x + pos.x - xlapN, y + pos.y,  w - xlapP + xlapN,  h},
                   Vector2{0, 0}, 0.0f, WHITE);
}

static void update_inputs_rl(Camera const& cam, Inputs& rInputs)
{
    Vector2 mousePos = GetMousePosition();
    Ray ray = GetMouseRay(mousePos, cam);

    rInputs.m_mousePos      = {mousePos.x, mousePos.y};
    rInputs.m_mouseOrig     = {ray.position.x, ray.position.y, ray.position.z};
    rInputs.m_mouseDir      = {ray.direction.x, ray.direction.y, ray.direction.z};
}

static McRaySalad lazor_salads(glm::vec3 origin, glm::vec3 dir, Salads_t const& salads)
{
    McRaySalad mcRayOut;
    mcRayOut.m_mcray.m_dist = std::numeric_limits<float>::max();
    mcRayOut.m_salad = -1;
    for (int i = 0; i < salads.size(); i ++)
    {
        if ( ! bool(salads[i]))
        {
            continue;
        }

        SaladModel const& salad = *salads[i];

        McRay mcray = shoop_da_woop_salad(origin, dir, salad);
        if (mcRayOut.m_mcray.m_dist > mcray.m_dist)
        {
            mcRayOut.m_mcray = mcray;
            mcRayOut.m_salad = i;
        }
    }
    return mcRayOut;
}

static int paw_default_base_attribute(meshdeform::MeshJoints const& joints, unsigned short const *pInd)
{
    // Figure out which joint is associated with a triangle
    // 3 vertices each with 4 joints each with a weight
    // make a mini histogram!
    int                             count       = 0;
    int                             top         = 0;
    static constexpr int            fox         = 12;
    std::array<unsigned char, fox>  cntJoints;
    std::array<float, fox>          cntWeights;


    // loop through 3 vertices
    for (int i = 0; i < 3; ++i)
    {
        int const           index           = (*pInd) * 4;
        unsigned char const *pJointInd      = &joints.m_pJointsIn[index];
        float const         *pWeight        = &joints.m_pWeightsIn[index];
        //std::cout << "Vertex " << (*pInd) << "\n";
        for (int j = 0; j < 4; ++j)
        {
            if ((*pWeight) < 0.001)
            {
                continue;
            }
            //std::cout << "* Joint: " << int(*pJointInd) << " @" << (*pWeight) << "\n";

            auto const  pBegin  = cntJoints.begin();
            auto const  pEnd    = cntJoints.begin() + count;
            int         found   = std::distance(pBegin, std::find(pBegin, pEnd, *pJointInd));

            float &rTotalWeight = cntWeights[found];

            if (found == count)
            {
                // new entry
                cntJoints[found]  = *pJointInd;
                rTotalWeight = *pWeight;
                ++count;
            }
            else
            {
                // already exists
                rTotalWeight += *pWeight;
            }

            // track highest on the fly
            if (cntWeights[top] < rTotalWeight)
            {
                top = found;
            }

            ++pJointInd;
            ++pWeight;
        }
        ++pInd;
    }

    return cntJoints[top];
}

static float rand_dist(float dist)
{
    float woot = GetRandomValue(-65536, 65536) / 65536.0f;
    return woot * woot * glm::sign(woot) * dist;
}

static void update_expressions(Soul &rSoul, float delta)
{
    rSoul.m_blinkCdn -= delta;

    if (rSoul.m_blinkCdn <= 0.0f)
    {
        rSoul.m_blinkCdn = rand_dist(rSoul.m_blinkPeriodMargin) + rSoul.m_blinkPeriodAvg;
    }

    rSoul.m_breathCycle += delta * rSoul.m_breathSpeed;
    rSoul.m_breathCycle -= float(rSoul.m_breathCycle > 1.0f);
}

bool g_limits{true};

static void update_tool_grab(
        Salads_t const& salads,
        WetJoints const& wet,
        FrogDyn& rFrogs,
        Inputs& rInputs,
        ToolGrab& rToolGrab)
{
    bool const selected = rInputs.m_selected == rToolGrab.m_id;

    if (rToolGrab.m_active)
    {
        if (!IsMouseButtonDown(MOUSE_BUTTON_LEFT))
        {
            // release
            rToolGrab.m_active = false;            
            int id = rToolGrab.m_grabs.back().baitId;
            rToolGrab.m_grabs.pop_back();
            auto found = std::find_if(rFrogs.m_baits.begin(), rFrogs.m_baits.end(), [id] (FrogDyn::Bait const& bait) { return bait.m_id == id; });

            rFrogs.m_baits.erase(found);
        }
        else if (IsKeyPressed(KEY_SPACE))
        {
            rToolGrab.m_active = false;
        }
    }
    else if (selected && IsMouseButtonPressed(MOUSE_BUTTON_LEFT))
    {
        // do grab
        if (rInputs.m_lazor.m_salad == -1)
        {
            rInputs.m_lazor = lazor_salads(rInputs.m_mouseOrig, rInputs.m_mouseDir, salads);
        }

        if (rInputs.m_lazor.m_salad != -1)
        {
            SaladModel const&               salad       = *salads.at(rInputs.m_lazor.m_salad);

            //std::cout << "Connected joint: " << int(cntJoints[top]) << "\n";

            // get joint, then linear search for connected frog
            int const joint = paw_default_base_attribute(salad.m_spookM, &salad.m_rayMesh.indices[rInputs.m_lazor.m_mcray.m_index * 3]);

            auto foundIt = std::find_if(
                        wet.m_hoppers.begin(), wet.m_hoppers.end(),
                        [joint] (WetJoints::Hopper const& hopper)
            {
                return hopper.m_joint == joint;
            });

            if (foundIt != wet.m_hoppers.end())
            {
                ToolGrab::Grab &rGrab = rToolGrab.m_grabs.emplace_back();

                frog_id_t const frog = foundIt->m_frog;

                // add bait
                glm::vec3 const down{0.0f, -1.0f, 0.0f};
                glm::vec3 const fwd{0.0f, 0.0f, 1.0f};
                glm::vec3 const lazorHit = rInputs.m_mouseOrig + rInputs.m_mouseDir * rInputs.m_lazor.m_mcray.m_dist;
                glm::vec3 const frogInvHit = glm::inverse(rFrogs.m_tf[frog]) * glm::vec4(lazorHit, 1.0f);

                rToolGrab.m_active = true;
                rGrab.baitId = 420 + rToolGrab.m_grabs.size();

                FrogDyn::Bait &rBait = rFrogs.m_baits.emplace_back(
                            FrogDyn::Bait{ {-1,  lazorHit, down, fwd},
                                           {frog,   frogInvHit, down, fwd} });
                rBait.m_id = rGrab.baitId;
                if (g_limits)
                {
                    rBait.m_forceLim = 15.2f * 10.0f * 3.0f;
                }
                //rBait.m_strengthMul = 0.2f;

                //rFrogs.m_extImp[frog].m_lin += rInputs.m_mouseDir * 10.0f;
            }
        }
    }

     if (rToolGrab.m_active)
     {
         ToolGrab::Grab &rGrab = rToolGrab.m_grabs.back();

         for (FrogDyn::Bait &rBait : rFrogs.m_baits)
         {
             if (rBait.m_id != rGrab.baitId)
             {
                 continue;
             }

             float const dist = glm::length(rBait.m_a.m_offset - rInputs.m_mouseOrig);
             rBait.m_a.m_offset = rInputs.m_mouseOrig + rInputs.m_mouseDir * dist;

             rFrogs.m_vel[rBait.m_b.m_id].m_lin *= 0.0f;
             //rFrogs.m_vel[rBait.m_b.m_id].m_ang *= 0.0f;

             break;
         }
     }
}

static float breath_cycle(float t) noexcept
{
    float const tau = glm::pi<float>() * 2.0f;
    return glm::sin(tau*(t - 0.25f)) + 0.2f * sin(tau*2*t) + 1.0;
}

static void draw_scene(TestSceneC &rScene, GameState &rGame)
{
    auto &t = rScene.m_time;

    float delta = 1.0f / 60.0f;//GetFrameTime();


    rScene.m_camYaw += (float(IsKeyDown(KEY_RIGHT)) - float(IsKeyDown(KEY_LEFT))) * delta * 3.14159f;
    rScene.m_camPitch += (float(IsKeyDown(KEY_UP)) - float(IsKeyDown(KEY_DOWN))) * delta * 3.14159f;
    rScene.m_camDist += (float(IsKeyDown(KEY_Z)) - float(IsKeyDown(KEY_X))) * delta * 2.0f;

    reinterpret_cast<glm::vec3&>(rScene.m_camera.position) = glm::quat(glm::vec3{rScene.m_camPitch, rScene.m_camYaw, 0.0f}) * glm::vec3{0.0f, 0.0f, rScene.m_camDist} + reinterpret_cast<glm::vec3&>(rScene.m_camera.target);

    CharB &rChar = rScene.m_characters.begin()->second;

    update_apples(rChar.m_apples, rChar.m_joints);

    glm::vec2 eyePosL = calc_eye_pos(rChar.m_apples.m_dataOut[rChar.m_eyeL.m_apple], reinterpret_cast<glm::vec3&>(rScene.m_camera.position));
    glm::vec2 eyePosR = calc_eye_pos(rChar.m_apples.m_dataOut[rChar.m_eyeR.m_apple], reinterpret_cast<glm::vec3&>(rScene.m_camera.position));


    float pushMag = (float(IsKeyDown(KEY_D)) - float(IsKeyDown(KEY_A))) * delta * 5.0f;

    glm::vec3 const push = glm::vec3{
                glm::sin(rScene.m_camYaw + 1.5708f),
                0.0f,
                glm::cos(rScene.m_camYaw + 1.5708f)} * pushMag;

    update_inputs_rl(rScene.m_camera, rScene.m_inputs);
    rScene.m_inputs.m_selected = 0;
    rScene.m_inputs.m_lazor.m_salad = -1;
    update_tool_grab(rScene.m_salads, rChar.m_wetJoints, rScene.m_frogs, rScene.m_inputs, rScene.m_toolGrab);

    if (IsKeyPressed(KEY_G))
    {
        rScene.m_gravity = !rScene.m_gravity;
    }

    if (IsKeyPressed(KEY_E))
    {
        g_limits = !g_limits;
    }

    if (IsKeyDown(KEY_R))
    {
        for (int i = 0; i < rScene.m_frogs.m_ids.size(); i ++)
        {
            rScene.m_frogs.m_vel[i].m_ang.x += 0.5f;
        }
    }

    if (IsKeyDown(KEY_T))
    {
        for (int i = 0; i < rScene.m_frogs.m_ids.size(); i ++)
        {
            rScene.m_frogs.m_vel[i].m_lin.y += 0.1f;
        }
    }

    if (IsKeyDown(KEY_F))
    {
//        rScene.m_frogs.m_scale[14] += 0.05f;
//        rScene.m_frogs.m_mass[14] += 0.02f;
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

        apply_baits(rScene.m_frogs, {64.0f, 0.3f, 16.0f, 0.5f, 8.0f, 0.2f, 8.0f, 0.25f}, smldelta);

        // beak-head spring workaround
        glm::vec3 const axis = rScene.m_frogs.m_tf[rChar.m_frogBeak][0];
        rScene.m_frogs.m_cstImp[rChar.m_frogBeak].m_ang -= axis * smldelta * 2.0f;
        rScene.m_frogs.m_cstImp[rChar.m_frogHead].m_ang += axis * smldelta * 2.0f;

        apply_cst_forces(rScene.m_frogs, smldelta);

        calc_balls_pos(rScene.m_frogs);

        calc_frog_collisions(rScene.m_frogs);
    }

    float const angdrag = delta * 1.0f;
    float const lindrag = delta * 2.0f;

    // apply drag
    for (int i = 0; i < rScene.m_frogs.m_ids.size(); i ++)
    {
        float const anglength = glm::length(rScene.m_frogs.m_vel[i].m_ang);
        if (anglength > 0.001)
        {
            rScene.m_frogs.m_vel[i].m_ang = rScene.m_frogs.m_vel[i].m_ang / anglength * glm::max(0.0f, anglength - angdrag);
        }

        float const linlength = glm::length(rScene.m_frogs.m_vel[i].m_lin);
        if (linlength > 0.001)
        {
            rScene.m_frogs.m_vel[i].m_lin = rScene.m_frogs.m_vel[i].m_lin / linlength * glm::max(0.0f, linlength - lindrag);
        }

        //rScene.m_frogs.m_vel[i].m_lin *= 0.998f;
    }

    // Update hoppers
    for (WetJoints::Hopper const& hopper : rChar.m_wetJoints.m_hoppers)
    {
        float const scale = rScene.m_frogs.m_scale[hopper.m_frog];
        rChar.m_joints.m_nodeTf.at(hopper.m_joint) = glm::translate(rScene.m_frogs.m_tf.at(hopper.m_frog) * glm::scale(glm::vec3{scale, scale, scale}), glm::vec3{0, hopper.m_yoffset, 0});
    }

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

        rlUpdateVertexBuffer(rScene.m_salads[i]->m_rayMesh.vboId[0], rScene.m_salads[i]->m_rayMesh.vertices, rScene.m_salads[i]->m_rayMesh.vertexCount*3*sizeof(float), 0);    // Update vertex position
        rlUpdateVertexBuffer(rScene.m_salads[i]->m_rayMesh.vboId[2], rScene.m_salads[i]->m_rayMesh.normals, rScene.m_salads[i]->m_rayMesh.vertexCount*3*sizeof(float), 0);     // Update vertex normals

    }

    update_expressions(rChar.m_soul, delta);

    // breathing
    rScene.m_frogs.m_scale[rChar.m_frogBelly] = 1.0f + 0.2f * breath_cycle(rChar.m_soul.m_breathCycle);

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

    //reinterpret_cast<glm::vec3&>(rScene.m_camera.target) =  avgPos;
    //reinterpret_cast<glm::vec3&>(rScene.m_camera.position) = glm::quat(glm::vec3{rScene.m_camPitch, rScene.m_camYaw, 0.0f}) * glm::vec3{0.0f, 0.0f, rScene.m_camDist} + reinterpret_cast<glm::vec3&>(rScene.m_camera.target);

    BeginDrawing();


        ClearBackground(Color{ 10, 10, 10, 100 });

        // draw eyes


        BeginTextureMode(rChar.m_eyeTexture);

            ClearBackground(Color{ 0, 0, 0, 0 });

            int const g = 128;
            int gx;
            int gy;

            // blink for 5 frames
            if (rChar.m_soul.m_blinkCdn < 5.0f/60.0f)
            {
                gx = 0; gy = 2;
            }
            else
            {
                gx = 0; gy = 1;
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

        BeginMode3D(rScene.m_camera);
            BeginBlendMode(BLEND_ALPHA);
                DrawGrid(10, 1.0f);

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


            BeginMode3D(rScene.m_camera);

#if 0
                // draw frogs
                for (frog_id_t id = 0; id < rScene.m_frogs.m_ids.capacity(); id ++)
                {
                    if (!rScene.m_frogs.m_ids.exists(id))
                    {
                        continue;
                    }

                    auto tf = rScene.m_frogs.m_tf[id];
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

            DrawTextEx(*rGame.m_pFont, "What is the best way to immobilize caT?", Vector2{10.0, 100.0}, 20, 0, WHITE);

            DrawRectangle(10, 10, 50, 25, Color{255, 255, 255, 255});

            DrawRectangle(10, 35, 50, 25, Color{203, 194, 201, 255});

            //DrawTexture(rChar.m_eyeTexture.texture, 0, 0, WHITE);

            BeginBlendMode(BLEND_CUSTOM);
            rlSetBlendFactors(0, 770, 32774);

                DrawRectangle(30, 20, 50, 50, Color{255, 0, 0, 100});

            EndBlendMode();

        EndTextureMode();

        auto &rTex = rScene.m_ui.texture;
        DrawTextureRec(rTex, Rectangle{0.0f, 0.0f, float(rTex.width), -float(rTex.height) * 1.0f}, Vector2{0.0f, 0.0f}, Color{255, 255, 255, 255});
    EndDrawing();

    t += delta;
}


#if defined(PLATFORM_DESKTOP)
    #define GLSL_VERSION            330
#else   // PLATFORM_RPI, PLATFORM_ANDROID, PLATFORM_WEB
    #define GLSL_VERSION            100
#endif

SceneFunc_t gen_test_scene_c()
{
    std::shared_ptr<TestSceneC> pScene = std::make_shared<TestSceneC>();
    TestSceneC &rScene = *pScene;

    std::string err;
    std::string warn;

    tinygltf::TinyGLTF loader;

    rScene.m_toolGrab.m_id = 0;

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

    rScene.m_salads.reserve(10);
    metal_pipe(rScene.m_gltf, 0, rScene.m_characters, rScene.m_salads, rScene.m_gltfRayMaterials, rScene.m_frogs);

    rScene.m_cube = GenMeshCube(0.04545f, 0.1f, 0.01136f);
    rScene.m_mat = LoadMaterialDefault();


    rScene.m_ui = LoadRenderTexture(GetScreenWidth(), GetScreenHeight());


    rScene.m_camera.target = Vector3{ 0.0f, 1.5156f, -0.033483f };
    rScene.m_camera.up = Vector3{ 0.0f, 1.0f, 0.0f };
    rScene.m_camera.fovy = 40.0f;
    rScene.m_camera.projection = CAMERA_PERSPECTIVE;

    return [pScene = std::move(pScene)] (GameState &rGame) -> void
    {
        draw_scene(*pScene, rGame);
    };
}

} // namespace orni
