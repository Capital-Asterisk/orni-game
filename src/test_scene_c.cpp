#include "scenes.hpp"
#include "frogdyn.hpp"

#include "mesh_deform.hpp"

#include <tiny_gltf.h>

#include <raylib.h>
#include <raymath.h>
#include <rlgl.h>

#include <glm/gtx/transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/matrix_interpolation.hpp>
#include <glm/gtx/string_cast.hpp>


#include <iostream>
#include <memory>
#include <unordered_map>

using namespace frogdyn;

namespace orni
{

using lgrn::IdRegistry;

using salad_id_t = int;
using apple_id_t = int;

struct Apple
{
    glm::mat4x4 m_tf;
    int m_jointParent;
    // int charId
};

struct Apples
{
    IdRegistry<apple_id_t>  m_ids;
    std::vector<Apple>      m_data;
    std::vector<glm::mat4>  m_dataOut;

    apple_id_t create(Apple apl)
    {
        apple_id_t const aplId = m_ids.create();
        m_data.resize(m_ids.capacity());
        m_dataOut.resize(m_ids.capacity());
        m_data[aplId] = apl;
        return aplId;
    }
};


struct SaladModel
{
    glm::vec3 const         *m_pPosIn;
    glm::vec3 const         *m_pNrmIn;

    std::vector<glm::vec3>  m_Pos;
    std::vector<glm::vec3>  m_Nrm;

    meshdeform::Targets     m_tgt;
    meshdeform::MeshJoints  m_spookM;
    int                     m_spookId;
    salad_id_t              m_sockOnId;
    Mesh                    m_rayMesh;
    Model                   m_rayModel;
};

struct WetJoints : meshdeform::Joints
{
    // Controls a Bait and a Joint with a frog parent
    struct Scorpion
    {
        frog_id_t           m_frog;
        glm::mat4x4         m_tf;
        int                 m_joint;
    };

    // frog directly controls a Joint
    struct Hopper
    {
        frog_id_t           m_frog;
        int                 m_joint;
        float               m_yoffset;
    };

    std::vector<Scorpion>   m_scorpions; // order-dependent
    std::vector<Hopper>     m_hoppers;
};

struct CharB
{
    struct Eye
    {
        apple_id_t          m_apple;
        glm::ivec2          m_irisPos;
    };


    Material                m_eyeMaterial;
    Texture                 m_eyeSheet;
    RenderTexture           m_eyeTexture;
    Eye                     m_eyeL;
    Eye                     m_eyeR;

    meshdeform::Joints      m_joints;
    WetJoints               m_wetJoints;

};

using Characters_t = std::unordered_map<int, CharB>;

struct TestSceneC
{
    std::vector<SaladModel> m_salads;
    Apples                  m_apples;
    std::vector<WetJoints>  m_spooks;
    FrogDyn                 m_frogs;

    Characters_t            m_characters;

    tinygltf::Model         m_gltf;
    std::vector<Image>      m_gltfRayImage;
    std::vector<Texture>    m_gltfRayTextures;
    std::vector<Material>   m_gltfRayMaterials;

    Camera3D                m_camera;

    Mesh                    m_cube;
    Material                m_mat;

    RenderTexture2D         m_ui;

    int                     m_cstSteps = 12;
    float                   m_extPercent = 0.5f;

    float                   m_camDist{3.0f};
    float                   m_camYaw{0.0f};
    float                   m_camPitch{0.0f};
    float                   m_time{0.0f};
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

static glm::vec2 calc_eye_pos(glm::mat4x4 const& eyeTf, glm::vec3 tgt)
{
    float eyeDepth = 0.03f;
    float eyeRadius = 0.09f;
    glm::mat4 eyeMatrix = glm::perspective(glm::atan(eyeRadius / eyeDepth) * 2.0f, 1.0f, 0.001f, 100.0f)
                        * glm::lookAt(glm::vec3(eyeTf[3]) - glm::vec3(eyeTf[2]) * eyeDepth,
                                      glm::vec3(eyeTf[3]) + glm::vec3(eyeTf[2]),
                                      glm::vec3(eyeTf[1]));

    glm::vec3 foo = eyeMatrix * glm::vec4(tgt, 1.0f);
    foo.x /= foo.z;
    foo.y /= foo.z;

    return glm::vec2(foo);
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

static void draw_scene(TestSceneC &rScene)
{
    auto &t = rScene.m_time;

    float delta = GetFrameTime();


    rScene.m_camYaw += (float(IsKeyDown(KEY_RIGHT)) - float(IsKeyDown(KEY_LEFT))) * delta * 3.14159f;
    rScene.m_camPitch += (float(IsKeyDown(KEY_UP)) - float(IsKeyDown(KEY_DOWN))) * delta * 3.14159f;
    rScene.m_camDist += (float(IsKeyDown(KEY_Z)) - float(IsKeyDown(KEY_X))) * delta * 2.0f;

    reinterpret_cast<glm::vec3&>(rScene.m_camera.position) = glm::quat(glm::vec3{rScene.m_camPitch, rScene.m_camYaw, 0.0f}) * glm::vec3{0.0f, 0.0f, rScene.m_camDist} + reinterpret_cast<glm::vec3&>(rScene.m_camera.target);
    //rScene.m_camera.position = Vector3{ std::sin(rScene.m_camAngle) * rScene.m_camDist, 3.0f, std::cos(rScene.m_camAngle) * rScene.m_camDist };
    //rScene.m_camera.position = Vector3{ std::sin(t * 3.14159f * 0.1f) * 2.0f, 2.0f, std::cos(t * 3.14159f * 0.1f) * 2.0f };

    //rScene.m_characters[0].m_joints.m_nodeTf[3] *= glm::axisAngleMatrix(glm::vec3{1.0f, 0.0f, 0.0f}, std::sin(t * glm::pi<float>() * 2.0f * 2.066666f) * 0.1f);
    //rScene.m_characters[0].m_joints.m_nodeTf[3] *= glm::scale(glm::vec3{0.998f, 0.998f, 0.998f});

    CharB &rChar = rScene.m_characters.begin()->second;

    update_apples(rScene.m_apples, rChar.m_joints);

    glm::vec2 eyePosL = calc_eye_pos(rScene.m_apples.m_dataOut[rChar.m_eyeL.m_apple], reinterpret_cast<glm::vec3&>(rScene.m_camera.position));
    glm::vec2 eyePosR = calc_eye_pos(rScene.m_apples.m_dataOut[rChar.m_eyeR.m_apple], reinterpret_cast<glm::vec3&>(rScene.m_camera.position));


    float pushMag = (float(IsKeyDown(KEY_D)) - float(IsKeyDown(KEY_A))) * delta * 5.0f;

    glm::vec3 const push = glm::vec3{
                glm::sin(rScene.m_camYaw + 1.5708f),
                0.0f,
                glm::cos(rScene.m_camYaw + 1.5708f)} * pushMag;

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

    for (int i = 0; i < 1; i ++)
    {
        // apply gravity
        //rScene.m_frogs.m_extImp[i].m_lin.y -= 9.81f * delta * rScene.m_frogs.m_mass[i];

        rScene.m_frogs.m_extImp[i].m_lin += push;

        auto p = glm::vec3(rScene.m_frogs.m_tf[i][3]);
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

    for (WetJoints::Hopper const& hopper : rChar.m_wetJoints.m_hoppers)
    {
        rChar.m_joints.m_nodeTf.at(hopper.m_joint) = glm::translate(rScene.m_frogs.m_tf.at(hopper.m_frog), glm::vec3{0, hopper.m_yoffset, 0});
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
                rScene.m_salads[i].m_spookM,
                rScene.m_salads[i].m_tgt,
                rScene.m_salads[i].m_pPosIn,
                rScene.m_salads[i].m_pNrmIn,
                0,
                rScene.m_salads[i].m_rayMesh.vertexCount,
                rScene.m_salads[i].m_Pos.data(),
                rScene.m_salads[i].m_Nrm.data());

        rlUpdateVertexBuffer(rScene.m_salads[i].m_rayMesh.vboId[0], rScene.m_salads[i].m_rayMesh.vertices, rScene.m_salads[i].m_rayMesh.vertexCount*3*sizeof(float), 0);    // Update vertex position
        rlUpdateVertexBuffer(rScene.m_salads[i].m_rayMesh.vboId[2], rScene.m_salads[i].m_rayMesh.normals, rScene.m_salads[i].m_rayMesh.vertexCount*3*sizeof(float), 0);     // Update vertex normals

    }

    BeginDrawing();

        ClearBackground(Color{ 10, 10, 10, 100 });

        // draw eyes


        BeginTextureMode(rChar.m_eyeTexture);

            ClearBackground(Color{ 0, 0, 0, 0 });

            // white stuff 1
            DrawTexturePro(rChar.m_eyeSheet, Rectangle{128, 128, -128, -128}, Rectangle{0, 128, 128, 128}, Vector2{0, 0}, 0.0f, WHITE);
            DrawTexturePro(rChar.m_eyeSheet, Rectangle{128, 128, 128, -128}, Rectangle{128, 128, 128, 128}, Vector2{0, 0}, 0.0f, WHITE);

            // irises
            draw_iris(rChar.m_eyeSheet, 0, eyePosR);
            draw_iris(rChar.m_eyeSheet, 1, eyePosL);

            // erase overlap
            BeginBlendMode(BLEND_CUSTOM);
                rlSetBlendFactors(0, 770, 32774);
                DrawTexturePro(rChar.m_eyeSheet, Rectangle{128, 128, -128, -128}, Rectangle{0, 128, 128, 128}, Vector2{0, 0}, 0.0f, WHITE);
                DrawTexturePro(rChar.m_eyeSheet, Rectangle{128, 128, 128, -128}, Rectangle{128, 128, 128, 128}, Vector2{0, 0}, 0.0f, WHITE);
            EndBlendMode();

            // outline
            DrawTexturePro(rChar.m_eyeSheet, Rectangle{128, 0, -128, -128}, Rectangle{0, 128, 128, 128}, Vector2{0, 0}, 0.0f, WHITE);
            DrawTexturePro(rChar.m_eyeSheet, Rectangle{128, 0, 128, -128}, Rectangle{128, 128, 128, 128}, Vector2{0, 0}, 0.0f, WHITE);
        EndTextureMode();

        BeginMode3D(rScene.m_camera);
            BeginBlendMode(BLEND_ALPHA);
                DrawGrid(10, 1.0f);

                DrawModel(rScene.m_salads[0].m_rayModel, Vector3{0.0f, 0.0f, 0.0f}, 1.0f, WHITE);
                DrawModel(rScene.m_salads[2].m_rayModel, Vector3{0.0f, 0.0f, 0.0f}, 1.0f, WHITE);
                DrawModel(rScene.m_salads[3].m_rayModel, Vector3{0.0f, 0.0f, 0.0f}, 1.0f, WHITE);
                DrawModel(rScene.m_salads[1].m_rayModel, Vector3{0.0f, 0.0f, 0.0f}, 1.0f, WHITE);
            EndBlendMode();
        EndMode3D();

        BeginTextureMode(rScene.m_ui);
            ClearBackground(Color{ 0, 0, 0, 0 });


            BeginMode3D(rScene.m_camera);
/*
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

                DrawSphere(Vector3{0.0f, 0.0f, 0.0f}, 0.05f, Color{255, 0, 0, 255});
                if (glm::mod(t, 0.5f) < 0.5f/2.0f)
                {
                    DrawSphere(reinterpret_cast<Vector3&>(rScene.m_apples.m_dataOut[rChar.m_eyeL.m_apple][3]), 0.01f, GREEN);
                }

                */
            EndMode3D();

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

template<typename T>
struct Burger
{
    T const *m_data;
    std::size_t m_count;
};

template<typename T>
Burger<T> drivethrough(tinygltf::Model const& gltf, int accessorId)
{
    auto const &access  = gltf.accessors.at(accessorId);
    auto const &view    = gltf.bufferViews.at(access.bufferView);
    auto const &buffer  = gltf.buffers.at(view.buffer);
    assert(view.byteStride == 0);
    return {reinterpret_cast<T const*>(&buffer.data[view.byteOffset + access.byteOffset]), access.count};
}

static glm::mat4 node_transform(tinygltf::Node const& node)
{
    glm::vec3 const scale = (node.scale.size() == 3) ? glm::vec3{node.scale[0], node.scale[1], node.scale[2]} : glm::vec3{1.0f};
    glm::vec3 const translate = (node.translation.size() == 3) ? glm::vec3{node.translation[0], node.translation[1], node.translation[2]} : glm::vec3{0.0f};
    glm::quat const rot = (node.rotation.size() == 4) ? glm::quat(node.rotation[3], node.rotation[0], node.rotation[1], node.rotation[2]) : glm::quat{};
    return glm::translate(translate) * glm::mat4(rot) * glm::scale(scale);
}

static int zero_workaround_lol = 0;

static void metal_rod(
        std::vector<int> const& nodeToJoint,
        tinygltf::Model const& gltf,
        int nodeId,
        int parentId,
        int level,
        FrogDyn& rFrogs,
        Apples& rApples,
        CharB& rChar,
        glm::mat4x4 parentTfWorld)
{
    auto const &node = gltf.nodes[nodeId];

    glm::mat4x4 const nodeTf = node_transform(node);
    glm::mat4x4 const nodeTfWorld = parentTfWorld * nodeTf;

    // check for mass
    if (node.extras.IsObject()) // wtf?
    {
        float comwok = 0.12f;

        if (tinygltf::Value const& massVal = node.extras.Get("mass"); massVal.IsNumber())
        {
            float mass = massVal.GetNumberAsDouble();
            std::cout << "MASSSS! " << mass << "\n";

            // Make a frog
            glm::mat4x4 frogTfWorld = nodeTfWorld;
            frogTfWorld[3] += frogTfWorld[1] * comwok;
            frog_id_t const frogId = frogdyn::add_frog(rFrogs, frogTfWorld, mass);

            // Add hopper
            rChar.m_wetJoints.m_hoppers.emplace_back(WetJoints::Hopper{frogId, nodeToJoint[nodeId], -comwok});

            // connect to parent
            if (level != 0)
            {
                int const parentJoint = nodeToJoint[parentId];
                assert(parentJoint != -1);
                auto foundIt = std::find_if(rChar.m_wetJoints.m_hoppers.begin(), rChar.m_wetJoints.m_hoppers.end(), [parentJoint] (WetJoints::Hopper const &wet) -> bool { return wet.m_joint == parentJoint; });
                assert(foundIt != rChar.m_wetJoints.m_hoppers.end());
                frog_id_t const parentFrogId = foundIt->m_joint;
                glm::vec3 offset{0.0f, -comwok, 0.0f};
                rFrogs.m_baits.emplace_back(
                            FrogDyn::Bait{ {parentFrogId, glm::vec3{nodeTf[3]} + offset},
                                           {frogId,   offset} });
            }
        }
    }

    for (int childId : node.children)
    {
        auto const &child = gltf.nodes.at(childId);

        if (nodeToJoint[childId] == -1)
        {
            // Not bone

            if (child.name.rfind("Apl_EyeL") == 0)
            {
                rChar.m_eyeL.m_apple = rApples.create(Apple{ node_transform(child), nodeToJoint[nodeId] });
            }
            else if (child.name.rfind("Apl_EyeR") == 0)
            {
                rChar.m_eyeR.m_apple = rApples.create(Apple{ node_transform(child), nodeToJoint[nodeId] });
            }
        }
        else
        {
            // Bone

            std::cout << child.name << "\n";

            metal_rod(nodeToJoint, gltf, childId, nodeId, level + 1, rFrogs, rApples, rChar, nodeTfWorld);
        }
    }
}

static void metal_bar(
        tinygltf::Model const&      gltf,
        int                         nodeId,
        CharB&                      rChar,
        std::vector<SaladModel>&    rSalads,
        Apples&                     rApples,
        std::vector<WetJoints> &    rSpooks,
        std::vector<Material>&      rMaterials,
        FrogDyn&                    rFrogs)
{
    auto const &node = gltf.nodes.at(nodeId);

    int useSkin = -1;

    // init eyes
    rChar.m_eyeSheet = LoadTexture("eyessheet.png");
    rChar.m_eyeTexture = LoadRenderTexture(256, 256);
    rChar.m_eyeMaterial = LoadMaterialDefault();
    SetMaterialTexture(&rChar.m_eyeMaterial, MATERIAL_MAP_DIFFUSE, rChar.m_eyeTexture.texture);

    // load mesh children
    for (int childId : node.children)
    {
        auto const &child = gltf.nodes.at(childId);

        if (child.mesh != -1)
        {
            SaladModel &rSalad = rSalads.emplace_back();

            if (child.skin != -1)
            {
                if (useSkin == -1)
                {
                    useSkin = child.skin;
                }
                else
                {
                    assert(useSkin == child.skin);
                }


                // make a salad model for this node

                auto const &mesh = gltf.meshes.at(child.mesh);
                auto const &prim = mesh.primitives.at(0);

                // Load vertex data
                {
                    auto posData = drivethrough<glm::vec3>(gltf, prim.attributes.at("POSITION"));
                    rSalad.m_pPosIn = posData.m_data;
                    rSalad.m_Pos.resize(posData.m_count);
                    rSalad.m_rayMesh.vertices = const_cast<float*>(reinterpret_cast<float const*>(rSalad.m_Pos.data())); // LOL!!!
                    //rSalad.m_rayMesh.vertices = const_cast<float*>(reinterpret_cast<float const*>(posData.m_data)); // LOL!!!
                    rSalad.m_rayMesh.vertexCount = posData.m_count;
                }
                {
                    auto nrmData = drivethrough<glm::vec3>(gltf, prim.attributes.at("NORMAL"));
                    rSalad.m_pNrmIn = nrmData.m_data;
                    rSalad.m_Nrm.resize(nrmData.m_count);
                    rSalad.m_rayMesh.normals = const_cast<float*>(reinterpret_cast<float const*>(rSalad.m_Nrm.data())); // LOL!!!
                }
                {
                    auto txCrdData = drivethrough<float>(gltf, prim.attributes.at("TEXCOORD_0"));
                    rSalad.m_rayMesh.texcoords = const_cast<float*>(txCrdData.m_data);
                }
                {
                    auto jointData = drivethrough<unsigned char>(gltf, prim.attributes.at("JOINTS_0"));
                    rSalad.m_spookM.m_pJointsIn = jointData.m_data;
                }
                {
                    auto weightData = drivethrough<float>(gltf, prim.attributes.at("WEIGHTS_0"));
                    rSalad.m_spookM.m_pWeightsIn = weightData.m_data;
                }

                // Load index data
                {
                    auto idxData = drivethrough<unsigned short>(gltf, prim.indices);
                    rSalad.m_rayMesh.indices = const_cast<unsigned short*>(idxData.m_data); // LOL!!!!!
                    rSalad.m_rayMesh.triangleCount = idxData.m_count / 3;
                }


                UploadMesh(&rSalad.m_rayMesh, true);
                rSalad.m_rayModel = {};
                rSalad.m_rayModel.transform = MatrixIdentity();
                rSalad.m_rayModel.meshCount = 1;
                rSalad.m_rayModel.meshes = &rSalad.m_rayMesh;

                rSalad.m_rayModel.materialCount = 1;
                rSalad.m_rayModel.meshMaterial = &zero_workaround_lol;

                if (child.name.rfind("Eyes") == 0)
                {
                    rSalad.m_rayModel.materials = &rChar.m_eyeMaterial;
                }
                else
                {
                    rSalad.m_rayModel.materials = &rMaterials.at(mesh.primitives.at(0).material);
                }
            }
            else
            {
                std::cout << "no skin?\n";
            }
        }
    }

    // load bones

    auto const &skin = gltf.skins.at(useSkin);
    int const jointCount = skin.joints.size();
    auto invBindData = drivethrough<glm::mat4x4>(gltf, skin.inverseBindMatrices);

    rChar.m_joints.m_pInverseBindIn = invBindData.m_data;

    // initialize world transforms
    rChar.m_joints.m_nodeTf.resize(jointCount);
    rChar.m_joints.m_jointTf.resize(jointCount);
    for (int i = 0; i < jointCount; i ++)
    {
        glm::mat4x4 const &rInvMatrix = rChar.m_joints.m_pInverseBindIn[i];
        rChar.m_joints.m_nodeTf[i] = glm::inverse(rInvMatrix);
    }

    // generate node-to-joint vector
    std::vector<int> nodeToJoint(gltf.nodes.size(), -1);
    {
        int currBone = 0;
        for (int wat : skin.joints)
        {
            nodeToJoint[wat] = currBone;
            currBone ++;
        }
    }

    // Iterate bones
    metal_rod(nodeToJoint, gltf, nodeId, -1, -1, rFrogs, rApples, rChar, glm::mat4x4(1.0f));
}

static void metal_pipe(
        tinygltf::Model const&      gltf,
        int                         sceneId,
        Characters_t&               rChars,
        std::vector<SaladModel>&    rSalads,
        Apples&                     rApples,
        std::vector<WetJoints> &    rSpooks,
        std::vector<Material>&      rMaterials,
        FrogDyn&                    rFrogs)
{
    auto const &scene = gltf.scenes.at(sceneId);

    // Load characters
    for (int nodeId : scene.nodes)
    {
        if (gltf.nodes[nodeId].name.rfind("Ch_") == 0)
        {
            // yeah, this is a character to load
            // TODO: multiple character support?
            CharB &rChar = rChars.emplace(0, CharB{}).first->second;
            metal_bar(gltf, nodeId, rChar, rSalads, rApples, rSpooks, rMaterials, rFrogs);
        }
    }
}

SceneFunc_t gen_test_scene_c()
{
    std::shared_ptr<TestSceneC> pScene = std::make_shared<TestSceneC>();
    TestSceneC &rScene = *pScene;

    std::string err;
    std::string warn;

    tinygltf::TinyGLTF loader;

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

    loader.LoadASCIIFromFile(&rScene.m_gltf, &err, &warn, "salad0.gltf");

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
    metal_pipe(rScene.m_gltf, 0, rScene.m_characters, rScene.m_salads, rScene.m_apples, rScene.m_spooks, rScene.m_gltfRayMaterials, rScene.m_frogs);

    rScene.m_cube = GenMeshCube(0.04545f, 0.1f, 0.01136f);
    rScene.m_mat = LoadMaterialDefault();

    rScene.m_ui = LoadRenderTexture(GetScreenWidth(), GetScreenHeight());


    rScene.m_camera.target = Vector3{ 0.0f, 1.5156f, -0.033483f };
    rScene.m_camera.up = Vector3{ 0.0f, 1.0f, 0.0f };
    rScene.m_camera.fovy = 40.0f;
    rScene.m_camera.projection = CAMERA_PERSPECTIVE;

    return [pScene = std::move(pScene)] (GameState &rGame) -> void
    {
        draw_scene(*pScene);
    };
}

} // namespace orni
