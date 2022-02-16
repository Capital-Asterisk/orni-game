#include "scenes.hpp"

#include "mesh_deform.hpp"

#include <tiny_gltf.h>

#include <raylib.h>
#include <raymath.h>
#include <rlgl.h>

#include <memory>

namespace orni
{

struct AnimMesh
{
    glm::vec3 const* m_pPosIn;
    glm::vec3 const* m_pNrmIn;

    std::vector<glm::vec3> m_Pos;
    std::vector<glm::vec3> m_Nrm;
};

struct SaladModel
{
    AnimMesh                m_anim;
    meshdeform::Targets     m_tgt;
    int                     m_joints;
    Model                   m_rayModel;

};

struct TestSceneC
{
    Camera3D m_camera;

    //meshdeform::Joints m_spooky;
    //meshdeform::Targets m_tgt;

    tinygltf::Model m_gltf;
    Mesh m_mesh;
    Material m_mat;

    RenderTexture2D m_ui;

    float m_time = 0.0f;
};


static void draw_scene(TestSceneC &rScene)
{
    auto &t = rScene.m_time;

    t += GetFrameTime();

    rScene.m_camera.position = Vector3{ std::sin(t * 3.14159f * 0.1f) * 8.0f, 5.0f, std::cos(t * 3.14159f * 0.1f) * 8.0f };

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

            EndMode3D();
        EndTextureMode();

        auto &rTex = rScene.m_ui.texture;
        DrawTextureRec(rTex, Rectangle{0.0f, 0.0f, float(rTex.width), -float(rTex.height) * 1.0f}, Vector2{0.0f, 0.0f}, Color{255, 255, 255, 255});
    EndDrawing();
}


SceneFunc_t gen_test_scene_c()
{
    std::shared_ptr<TestSceneC> pScene = std::make_shared<TestSceneC>();
    TestSceneC &rScene = *pScene;

    std::string err;
    std::string warn;

    tinygltf::TinyGLTF loader;

    loader.SetImageLoader([] (
            tinygltf::Image *img, const int index,
            std::string *pErr, std::string *pWarn, int foo, int bar,
            const unsigned char *pData, int size, void *pUser) -> bool
    {
        Image rayImg = LoadImageFromMemory(".png", pData, size);
        return true;
    }, nullptr);

    //image, image_idx, err, warn, 0, 0, &img.at(0),
    //static_cast<int>(img.size()), load_image_user_data);

    loader.LoadASCIIFromFile(&rScene.m_gltf, &err, &warn, "salad0.gltf");

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
    rScene.m_mesh.vertices = reinterpret_cast<float*>(&rBuffers[rVrtxView.buffer].data[rVrtxView.byteOffset + rVrtxAccessor.byteOffset]);
    }

    {
    auto &rNrmlAccessor = rAccessors[rMeshes[0].primitives[0].attributes["NORMAL"]];
    auto &rNrmlView = rViews[rNrmlAccessor.bufferView];
    rScene.m_mesh.normals = reinterpret_cast<float*>(&rBuffers[rNrmlView.buffer].data[rNrmlView.byteOffset + rNrmlAccessor.byteOffset]);
    }


    {
    auto &rIndxAccessor = rAccessors[rMeshes[0].primitives[0].indices];
    rScene.m_mesh.triangleCount = rIndxAccessor.count;
    auto &rIndxView = rViews[rIndxAccessor.bufferView];
    rScene.m_mesh.indices = reinterpret_cast<unsigned short*>(&rBuffers[rIndxView.buffer].data[rIndxView.byteOffset + rIndxAccessor.byteOffset]);
    }

    // add targets

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

} // namespace orni
