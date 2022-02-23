#include "scenes.hpp"
#include "salads.hpp"

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

namespace orni
{



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
        float smldelta = delta * cstPercent / float(rScene.m_cstSteps);

        apply_baits(rScene.m_frogs, {16.0f, 0.1f, 16.0f, 0.1f}, smldelta);

        apply_cst_forces(rScene.m_frogs, smldelta);

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
