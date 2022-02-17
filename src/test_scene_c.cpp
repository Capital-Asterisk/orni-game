#include "scenes.hpp"
#include "frogdyn.hpp"

#include "mesh_deform.hpp"

#include <tiny_gltf.h>

#include <raylib.h>
#include <raymath.h>
#include <rlgl.h>

#include <memory>
#include <iostream>

using namespace frogdyn;

namespace orni
{

using salad_id_t = int;

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
};

struct WetJoints : meshdeform::Joints
{
    struct Scorpion
    {
        glm::mat4x4         m_tf;
        int                 m_jointParent;
        int                 m_jointChild;
    };

    struct Hopper
    {
        int                 m_joint;
        frog_id_t           m_frog;
    };

    std::vector<Scorpion>   m_scorpions;
    std::vector<Hopper>     m_hoppers;
};

struct CharB
{
    RenderTexture           m_eyeTexture;
    meshdeform::Joints      m_joints;
};

using Characters_t = std::unordered_map<int, CharB>;

struct TestSceneC
{
    std::vector<SaladModel> m_salads;
    std::vector<WetJoints>  m_spooks;
    FrogDyn                 m_frogs;

    Characters_t            m_characters;

    tinygltf::Model         m_gltf;

    Camera3D                m_camera;

    Material                m_mat;

    RenderTexture2D         m_ui;

    float                   m_time{0.0f};
};


static void draw_scene(TestSceneC &rScene)
{
    auto &t = rScene.m_time;

    t += GetFrameTime();

    rScene.m_camera.position = Vector3{ std::sin(t * 3.14159f * 0.1f) * 2.0f, 2.0f, std::cos(t * 3.14159f * 0.1f) * 2.0f };

    meshdeform::calculate_joint_transforms(
            (glm::mat4x4(1.0f)),
            rScene.m_characters[0].m_joints.m_pInverseBindIn,
            rScene.m_characters[0].m_joints.m_nodeTf.data(),
            0,
            rScene.m_characters[0].m_joints.m_jointTf.size(),
            rScene.m_characters[0].m_joints.m_jointTf.data());

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

        ClearBackground(Color{ 10, 10, 10, 240 });

        BeginMode3D(rScene.m_camera);

            DrawGrid(10, 1.0f);
            DrawMesh(rScene.m_salads[0].m_rayMesh, rScene.m_mat, MatrixIdentity());
            DrawMesh(rScene.m_salads[1].m_rayMesh, rScene.m_mat, MatrixIdentity());
            DrawMesh(rScene.m_salads[2].m_rayMesh, rScene.m_mat, MatrixIdentity());
            DrawMesh(rScene.m_salads[3].m_rayMesh, rScene.m_mat, MatrixIdentity());

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

static void metal_bar(
        tinygltf::Model const&      gltf,
        int                         nodeId,
        CharB&                      rChar,
        std::vector<SaladModel>&    rSalads,
        std::vector<WetJoints> &    rSpooks,
        FrogDyn&                    rFrogs)
{
    auto const &root = gltf.nodes.at(nodeId);

    int useSkin = -1;

    // load mesh children
    for (int childId : root.children)
    {
        auto const &child = gltf.nodes.at(childId);

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
            SaladModel &rSalad = rSalads.emplace_back();
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
        }
        else
        {
            std::cout << "no skin?\n";
        }

        if (child.name.rfind("Eyes") == 0)
        {

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

}

static void metal_pipe(
        tinygltf::Model const&      gltf,
        int                         sceneId,
        Characters_t&               rChars,
        std::vector<SaladModel>&    rSalads,
        std::vector<WetJoints> &    rSpooks,
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
            metal_bar(gltf, nodeId, rChar, rSalads, rSpooks, rFrogs);
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

    rScene.m_salads.reserve(10);
    metal_pipe(rScene.m_gltf, 0, rScene.m_characters, rScene.m_salads, rScene.m_spooks, rScene.m_frogs);

    rScene.m_mat = LoadMaterialDefault();

    rScene.m_ui = LoadRenderTexture(GetScreenWidth(), GetScreenHeight());


    rScene.m_camera.target = Vector3{ 0.0f, 1.3f, 0.0f };
    rScene.m_camera.up = Vector3{ 0.0f, 1.0f, 0.0f };
    rScene.m_camera.fovy = 50.0f;
    rScene.m_camera.projection = CAMERA_PERSPECTIVE;

    return [pScene = std::move(pScene)] (GameState &rGame) -> void
    {
        draw_scene(*pScene);
    };
}

} // namespace orni
