#include "scenes.hpp"

#include "mesh_deform.hpp"

#include <tiny_gltf.h>

#include <raylib.h>
#include <raymath.h>
#include <rlgl.h>

#include <memory>

using namespace orni;

struct LazyAnimMesh
{
    glm::vec3 const* m_pPosIn;
    glm::vec3 const* m_pNrmIn;

    std::vector<glm::vec3> m_Pos;
    std::vector<glm::vec3> m_Nrm;
};

struct TestSceneA
{
    Camera3D m_camera;

    meshdeform::MeshJoints m_meshSpooky;
    meshdeform::Joints m_spooky;
    meshdeform::Targets m_tgt;
    LazyAnimMesh m_animMesh;

    tinygltf::Model m_gltf;
    Mesh m_mesh;
    Material m_mat;

    RenderTexture2D m_ui;

    float m_time = 0.0f;
};


static void draw_scene(TestSceneA &rScene)
{
    auto &t = rScene.m_time;

    t += GetFrameTime();

    rScene.m_camera.position = Vector3{ std::sin(t * 3.14159f * 0.1f) * 8.0f, 5.0f, std::cos(t * 3.14159f * 0.1f) * 8.0f };

    rScene.m_tgt.m_targets[0].m_weight = std::sin(t * 6.0f) + 1.0f;

    rScene.m_spooky.m_nodeTf[0][3].x = std::sin(t * 4.0f) * 1.0f;
    rScene.m_spooky.m_nodeTf[1][3].x = std::sin(t * 4.0f + 1.0f) * 1.0f;
    rScene.m_spooky.m_nodeTf[2][3].x = std::sin(t * 4.0f + 2.0f) * 1.0f;
    rScene.m_spooky.m_nodeTf[3][3].x = std::sin(t * 4.0f + 3.0f) * 1.0f;
    //rScene.m_spooky.m_nodeTf[3] = glm::rotate(rScene.m_spooky.m_nodeTf[3], 0.1f, glm::vec3(1.0f, 0.0f, 0.0f));

    meshdeform::calculate_joint_transforms(
            (glm::mat4x4(1.0f)),
            rScene.m_spooky.m_pInverseBindIn,
            rScene.m_spooky.m_nodeTf.data(),
            0,
            rScene.m_spooky.m_jointTf.size(),
            rScene.m_spooky.m_jointTf.data());

    meshdeform::apply_vertex_transform(
            rScene.m_spooky,
            rScene.m_meshSpooky,
            rScene.m_tgt,
            rScene.m_animMesh.m_pPosIn,
            rScene.m_animMesh.m_pNrmIn,
            0,
            rScene.m_mesh.vertexCount,
            rScene.m_animMesh.m_Pos.data(),
            rScene.m_animMesh.m_Nrm.data());

    rlUpdateVertexBuffer(rScene.m_mesh.vboId[0], rScene.m_mesh.vertices, rScene.m_mesh.vertexCount*3*sizeof(float), 0);    // Update vertex position
    rlUpdateVertexBuffer(rScene.m_mesh.vboId[2], rScene.m_mesh.normals, rScene.m_mesh.vertexCount*3*sizeof(float), 0);     // Update vertex normals


    BeginDrawing();

        ClearBackground(Color{ 10, 10, 10, 240 });

        BeginMode3D(rScene.m_camera);

            DrawGrid(10, 1.0f);
            DrawMesh(rScene.m_mesh, rScene.m_mat, MatrixIdentity());




        EndMode3D();

        BeginTextureMode(rScene.m_ui);
            ClearBackground(Color{ 0, 0, 0, 0 });
            BeginMode3D(rScene.m_camera);
            DrawSphere(Vector3{0.0f, 0.0f, 0.0f}, 0.05f, Color{255, 0, 0, 255});

            // draw all the joints
            for (int i = 0; i < rScene.m_spooky.m_jointTf.size(); i ++)
            {
                glm::vec4 const &ni = rScene.m_spooky.m_nodeTf[i][3];
                DrawSphere(Vector3{ni.x, ni.y, ni.z}, 0.05f, Color{100, 255, 0, 255});

            }
            EndMode3D();
        EndTextureMode();

        auto &rTex = rScene.m_ui.texture;
        DrawTextureRec(rTex, Rectangle{0.0f, 0.0f, float(rTex.width), -float(rTex.height) * 1.0f}, Vector2{0.0f, 0.0f}, Color{255, 255, 255, 255});
    EndDrawing();
}


SceneFunc_t orni::gen_test_scene_a()
{
    std::shared_ptr<TestSceneA> pScene = std::make_shared<TestSceneA>();
    TestSceneA &rScene = *pScene;

    std::string err;
    std::string warn;

    tinygltf::TinyGLTF loader;
    loader.LoadBinaryFromFile(&rScene.m_gltf, &err, &warn, "stick5.glb");

    auto &rMeshes = rScene.m_gltf.meshes;
    auto &rAccessors = rScene.m_gltf.accessors;
    auto &rViews = rScene.m_gltf.bufferViews;
    auto &rBuffers = rScene.m_gltf.buffers;
    auto &rSkins = rScene.m_gltf.skins;

    rScene.m_mesh = Mesh();

    {
    auto &rVrtxAccessor = rAccessors[rMeshes[0].primitives[0].attributes["POSITION"]];
    rScene.m_mesh.vertexCount = rVrtxAccessor.count;
    auto &rVrtxView = rViews[rVrtxAccessor.bufferView];
    rScene.m_animMesh.m_pPosIn = reinterpret_cast<glm::vec3*>(&rBuffers[rVrtxView.buffer].data[rVrtxView.byteOffset + rVrtxAccessor.byteOffset]);
    rScene.m_animMesh.m_Pos.resize(rScene.m_mesh.vertexCount);
    rScene.m_mesh.vertices = reinterpret_cast<float*>(rScene.m_animMesh.m_Pos.data());
    }

    {
    auto &rNrmlAccessor = rAccessors[rMeshes[0].primitives[0].attributes["NORMAL"]];
    auto &rNrmlView = rViews[rNrmlAccessor.bufferView];
    rScene.m_animMesh.m_pNrmIn = reinterpret_cast<glm::vec3*>(&rBuffers[rNrmlView.buffer].data[rNrmlView.byteOffset + rNrmlAccessor.byteOffset]);
    rScene.m_animMesh.m_Nrm.resize(rScene.m_mesh.vertexCount);
    rScene.m_mesh.normals = reinterpret_cast<float*>(rScene.m_animMesh.m_Nrm.data());
    }

    {
    auto &rNrmlAccessor = rAccessors[rMeshes[0].primitives[0].attributes["JOINTS_0"]];
    auto &rNrmlView = rViews[rNrmlAccessor.bufferView];
    rScene.m_meshSpooky.m_pJointsIn = reinterpret_cast<unsigned char*>(&rBuffers[rNrmlView.buffer].data[rNrmlView.byteOffset + rNrmlAccessor.byteOffset]);
    }

    {
    auto &rNrmlAccessor = rAccessors[rMeshes[0].primitives[0].attributes["WEIGHTS_0"]];
    auto &rNrmlView = rViews[rNrmlAccessor.bufferView];
    rScene.m_meshSpooky.m_pWeightsIn = reinterpret_cast<float*>(&rBuffers[rNrmlView.buffer].data[rNrmlView.byteOffset + rNrmlAccessor.byteOffset]);
    }

    {
    int const jointCount = rSkins[0].joints.size();
    auto &rNrmlAccessor = rAccessors[rSkins[0].inverseBindMatrices];
    auto &rNrmlView = rViews[rNrmlAccessor.bufferView];

    rScene.m_spooky.m_pInverseBindIn = reinterpret_cast<glm::mat4x4*>(&rBuffers[rNrmlView.buffer].data[rNrmlView.byteOffset + rNrmlAccessor.byteOffset]);

    // initialize world transforms
    rScene.m_spooky.m_nodeTf.resize(jointCount);
    rScene.m_spooky.m_jointTf.resize(jointCount);
    for (int i = 0; i < jointCount; i ++)
    {
        glm::mat4x4 const &rInvMatrix = rScene.m_spooky.m_pInverseBindIn[i];
        rScene.m_spooky.m_nodeTf[i] = glm::mat4x4(1.0f);
        rScene.m_spooky.m_nodeTf[i][3].y = i * 0.5f;
    }

    }

    {
    auto &rIndxAccessor = rAccessors[rMeshes[0].primitives[0].indices];
    rScene.m_mesh.triangleCount = rIndxAccessor.count;
    auto &rIndxView = rViews[rIndxAccessor.bufferView];
    rScene.m_mesh.indices = reinterpret_cast<unsigned short*>(&rBuffers[rIndxView.buffer].data[rIndxView.byteOffset + rIndxAccessor.byteOffset]);
    }

    // add targets
    {
    auto const& targets = rMeshes[0].primitives[0].targets;
    int const tgtCount = targets.size();

    rScene.m_tgt.m_targets.resize(tgtCount);
    for (int i = 0; i < tgtCount; i ++)
    {
        meshdeform::Targets::Target &rTgt = rScene.m_tgt.m_targets[i];
        {
        auto &rVrtxAccessor = rAccessors[targets[i].at("POSITION")];
        rScene.m_mesh.vertexCount = rVrtxAccessor.count;
        auto &rVrtxView = rViews[rVrtxAccessor.bufferView];
        rTgt.m_pPosIn = reinterpret_cast<glm::vec3*>(&rBuffers[rVrtxView.buffer].data[rVrtxView.byteOffset + rVrtxAccessor.byteOffset]);
        }
        {
        auto &rVrtxAccessor = rAccessors[targets[i].at("NORMAL")];
        rScene.m_mesh.vertexCount = rVrtxAccessor.count;
        auto &rVrtxView = rViews[rVrtxAccessor.bufferView];
        rTgt.m_pNrmIn = reinterpret_cast<glm::vec3*>(&rBuffers[rVrtxView.buffer].data[rVrtxView.byteOffset + rVrtxAccessor.byteOffset]);
        }
    }
    }

    rScene.m_mat = LoadMaterialDefault();

    rScene.m_ui = LoadRenderTexture(GetScreenWidth(), GetScreenHeight());

    UploadMesh(&rScene.m_mesh, true);


    rScene.m_camera.target = Vector3{ 0.0f, 0.0f, 0.0f };
    rScene.m_camera.up = Vector3{ 0.0f, 1.0f, 0.0f };
    rScene.m_camera.fovy = 50.0f;
    rScene.m_camera.projection = CAMERA_PERSPECTIVE;

    return [pScene = std::move(pScene)] (GameState &rGame) -> void
    {
        draw_scene(*pScene);
    };
}
