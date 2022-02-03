

#include <raylib.h>
#include <raymath.h>
#include <rlgl.h>
#include <glm/mat4x4.hpp>
#include <glm/ext.hpp>


#include <tiny_gltf.h>

#if defined(PLATFORM_WEB)
    #include <emscripten/emscripten.h>
#endif

#include <iostream>

struct Spooky
{
    glm::mat4x4 const   *m_pInverseBindIn;
    float const         *m_pWeightsIn;
    unsigned char const *m_pJointsIn;


    std::vector<glm::mat4x4>    m_nodeTf;
    std::vector<glm::mat4x4>    m_jointTf;
};

struct Targets
{
    struct Target
    {
        glm::vec3 const *m_pPosIn;
        glm::vec3 const *m_pNrmIn;
        float m_weight;
    };

    std::vector<Target> m_targets;
};

struct LazyAnimMesh
{
    glm::vec3 const* m_pPosIn;
    glm::vec3 const* m_pNrmIn;

    std::vector<glm::vec3> m_Pos;
    std::vector<glm::vec3> m_Nrm;
};

struct TestScene
{
    Camera3D m_camera;

    Spooky m_spooky;
    Targets m_tgt;
    LazyAnimMesh m_animMesh;

    tinygltf::Model m_gltf;
    Mesh m_mesh;
    Material m_mat;

    RenderTexture2D m_ui;
};

void calculate_joint_transforms(
        int const           jointCount,
        glm::mat4x4 const   &invWorldTf,
        glm::mat4x4 const   *pInverseBind,
        glm::mat4x4 const   *pNodeTf,
        std::size_t const   first,
        std::size_t const   last,
        glm::mat4x4         *pJointTfOut)
{
    glm::mat4x4 const *pInverseBindCurr = &pInverseBind[first];
    glm::mat4x4 const *pNodeTfCurr      = &pNodeTf[first];

    glm::mat4x4 *pJointTfOutCurr = &pJointTfOut[first];

    for (std::size_t i = first; i < last; i ++)
    {
        *pJointTfOutCurr = invWorldTf * (*pNodeTfCurr) * (*pInverseBindCurr);

        pInverseBindCurr ++;
        pNodeTfCurr ++;
        pJointTfOutCurr ++;
    }
}

void apply_vertex_transform(
        Spooky const        &spooky,
        Targets const&      tgt,
        glm::vec3 const*    pPosIn,
        glm::vec3 const*    pNrmIn,
        std::size_t const   first,
        std::size_t const   last,
        glm::vec3           *pPosOut,
        glm::vec3           *pNrmOut)
{
    glm::vec3 *pPosOutCurr = &pPosOut[first];
    glm::vec3 *pNrmOutCurr = &pNrmOut[first];

    glm::vec3 const *pPosInCurr = &pPosIn[first];
    glm::vec3 const *pNrmInCurr = &pNrmIn[first];

    float const *pWeightCurr = &spooky.m_pWeightsIn[first * 4];
    unsigned char const *pJointCurr = &spooky.m_pJointsIn[first * 4];

    for (std::size_t i = first; i < last; i ++)
    {
        glm::vec3 pos(*pPosInCurr);
        glm::vec3 nrm(*pNrmOutCurr);

        pPosInCurr ++;
        pNrmInCurr ++;

        // Do targets
        for (Targets::Target const &tgt : tgt.m_targets)
        {
            if (tgt.m_weight == 0)
            {
                continue;
            }
            pos += tgt.m_pPosIn[i] * tgt.m_weight;
            nrm += tgt.m_pNrmIn[i] * tgt.m_weight;
        }

        glm::mat4x4 jointMatrix(0.0f);

        // Do joints
        for (int j = 0; j < 4; j ++)
        {
            float const weight = pWeightCurr[j];

            if (weight == 0.0f)
            {
                break;
            }

            unsigned char const joint = pJointCurr[j];

            jointMatrix += weight * spooky.m_jointTf[joint];

            //posJoint += glm::vec3(jointMatrix * glm::vec4(posTgt, 1.0f) * weight);
            //nrmJoint += glm::vec3(jointMatrix * glm::vec4(nrmTgt, 1.0f) * weight);
        }

        pos = glm::vec3(jointMatrix * glm::vec4(pos, 1.0f));
        nrm = glm::vec3(jointMatrix * glm::vec4(nrm, 0.0f));

        pWeightCurr += 4;
        pJointCurr += 4;

        *pPosOutCurr = pos;
        *pNrmOutCurr = glm::normalize(nrm);

        pPosOutCurr ++;
        pNrmOutCurr ++;
    }
}


TestScene g_scene;

float t = 0.0f;

void update_draw_frame()
{
    t += GetFrameTime();

    g_scene.m_camera.position = Vector3{ std::sin(t * 3.14159f * 0.1f) * 8.0f, 5.0f, std::cos(t * 3.14159f * 0.1f) * 8.0f };

    g_scene.m_tgt.m_targets[0].m_weight = std::sin(t * 6.0f) + 1.0f;

    g_scene.m_spooky.m_nodeTf[0][3].x = std::sin(t * 4.0f) * 1.0f;
    g_scene.m_spooky.m_nodeTf[1][3].x = std::sin(t * 4.0f + 1.0f) * 1.0f;
    g_scene.m_spooky.m_nodeTf[2][3].x = std::sin(t * 4.0f + 2.0f) * 1.0f;
    g_scene.m_spooky.m_nodeTf[3][3].x = std::sin(t * 4.0f + 3.0f) * 1.0f;
    //g_scene.m_spooky.m_nodeTf[3] = glm::rotate(g_scene.m_spooky.m_nodeTf[3], 0.1f, glm::vec3(1.0f, 0.0f, 0.0f));

    calculate_joint_transforms(
            g_scene.m_spooky.m_jointTf.size(),
            (glm::mat4x4(1.0f)),
            g_scene.m_spooky.m_pInverseBindIn,
            g_scene.m_spooky.m_nodeTf.data(),
            0,
            g_scene.m_spooky.m_jointTf.size(),
            g_scene.m_spooky.m_jointTf.data());

    apply_vertex_transform(
            g_scene.m_spooky,
            g_scene.m_tgt,
            g_scene.m_animMesh.m_pPosIn,
            g_scene.m_animMesh.m_pNrmIn,
            0,
            g_scene.m_mesh.vertexCount,
            g_scene.m_animMesh.m_Pos.data(),
            g_scene.m_animMesh.m_Nrm.data());

    rlUpdateVertexBuffer(g_scene.m_mesh.vboId[0], g_scene.m_mesh.vertices, g_scene.m_mesh.vertexCount*3*sizeof(float), 0);    // Update vertex position
    rlUpdateVertexBuffer(g_scene.m_mesh.vboId[2], g_scene.m_mesh.normals, g_scene.m_mesh.vertexCount*3*sizeof(float), 0);     // Update vertex normals


    BeginDrawing();

        ClearBackground(Color{ 10, 10, 10, 240 });

        BeginMode3D(g_scene.m_camera);

            DrawGrid(10, 1.0f);
            DrawMesh(g_scene.m_mesh, g_scene.m_mat, MatrixIdentity());




        EndMode3D();

        BeginTextureMode(g_scene.m_ui);
            ClearBackground(Color{ 0, 0, 0, 0 });
            BeginMode3D(g_scene.m_camera);
            // draw all the joints
            for (int i = 0; i < g_scene.m_spooky.m_jointTf.size(); i ++)
            {
                glm::vec4 const &ni = g_scene.m_spooky.m_nodeTf[i][3];
                DrawSphere(Vector3{ni.x, ni.y, ni.z}, 0.05f, Color{100, 255, 0, 255});

            }
            EndMode3D();
        EndTextureMode();

        auto &rTex = g_scene.m_ui.texture;
        DrawTextureRec(rTex, Rectangle{0.0f, rTex.height * 1.0f, rTex.width * 1.0f, -rTex.height * 1.0f}, Vector2{0.0f, 0.0f}, Color{255, 255, 255, 255});
    EndDrawing();



}

int main(int argc, char** argv)
{
    const int screenWidth = 800;
    const int screenHeight = 600;

    std::string err;
    std::string warn;

    tinygltf::TinyGLTF loader;
    loader.LoadASCIIFromFile(&g_scene.m_gltf, &err, &warn, "stick4.gltf");

    auto &fish = g_scene;
    auto &rMeshes = g_scene.m_gltf.meshes;
    auto &rAccessors = g_scene.m_gltf.accessors;
    auto &rViews = g_scene.m_gltf.bufferViews;
    auto &rBuffers = g_scene.m_gltf.buffers;
    auto &rSkins = g_scene.m_gltf.skins;

    g_scene.m_mesh = Mesh();

    {
    auto &rVrtxAccessor = rAccessors[rMeshes[0].primitives[0].attributes["POSITION"]];
    g_scene.m_mesh.vertexCount = rVrtxAccessor.count;
    auto &rVrtxView = rViews[rVrtxAccessor.bufferView];
    g_scene.m_animMesh.m_pPosIn = reinterpret_cast<glm::vec3*>(&rBuffers[rVrtxView.buffer].data[rVrtxView.byteOffset + rVrtxAccessor.byteOffset]);
    g_scene.m_animMesh.m_Pos.resize(g_scene.m_mesh.vertexCount);
    g_scene.m_mesh.vertices = reinterpret_cast<float*>(g_scene.m_animMesh.m_Pos.data());
    }

    {
    auto &rNrmlAccessor = rAccessors[rMeshes[0].primitives[0].attributes["NORMAL"]];
    auto &rNrmlView = rViews[rNrmlAccessor.bufferView];
    g_scene.m_animMesh.m_pNrmIn = reinterpret_cast<glm::vec3*>(&rBuffers[rNrmlView.buffer].data[rNrmlView.byteOffset + rNrmlAccessor.byteOffset]);
    g_scene.m_animMesh.m_Nrm.resize(g_scene.m_mesh.vertexCount);
    g_scene.m_mesh.normals = reinterpret_cast<float*>(g_scene.m_animMesh.m_Nrm.data());
    }

    {
    auto &rNrmlAccessor = rAccessors[rMeshes[0].primitives[0].attributes["JOINTS_0"]];
    auto &rNrmlView = rViews[rNrmlAccessor.bufferView];
    g_scene.m_spooky.m_pJointsIn = reinterpret_cast<unsigned char*>(&rBuffers[rNrmlView.buffer].data[rNrmlView.byteOffset + rNrmlAccessor.byteOffset]);
    }

    {
    auto &rNrmlAccessor = rAccessors[rMeshes[0].primitives[0].attributes["WEIGHTS_0"]];
    auto &rNrmlView = rViews[rNrmlAccessor.bufferView];
    g_scene.m_spooky.m_pWeightsIn = reinterpret_cast<float*>(&rBuffers[rNrmlView.buffer].data[rNrmlView.byteOffset + rNrmlAccessor.byteOffset]);
    }

    {
    int const jointCount = rSkins[0].joints.size();
    auto &rNrmlAccessor = rAccessors[rSkins[0].inverseBindMatrices];
    auto &rNrmlView = rViews[rNrmlAccessor.bufferView];

    g_scene.m_spooky.m_pInverseBindIn = reinterpret_cast<glm::mat4x4*>(&rBuffers[rNrmlView.buffer].data[rNrmlView.byteOffset + rNrmlAccessor.byteOffset]);

    // initialize world transforms
    g_scene.m_spooky.m_nodeTf.resize(jointCount);
    g_scene.m_spooky.m_jointTf.resize(jointCount);
    for (int i = 0; i < jointCount; i ++)
    {
        glm::mat4x4 const &rInvMatrix = g_scene.m_spooky.m_pInverseBindIn[i];
        g_scene.m_spooky.m_nodeTf[i] = glm::mat4x4(1.0f);
        g_scene.m_spooky.m_nodeTf[i][3].y = i * 0.5f;
    }

    std::cout << "lol\n";
    }

    {
    auto &rIndxAccessor = rAccessors[rMeshes[0].primitives[0].indices];
    g_scene.m_mesh.triangleCount = rIndxAccessor.count;
    auto &rIndxView = rViews[rIndxAccessor.bufferView];
    g_scene.m_mesh.indices = reinterpret_cast<unsigned short*>(&rBuffers[rIndxView.buffer].data[rIndxView.byteOffset + rIndxAccessor.byteOffset]);
    }

    // add targets
    {
    auto const& targets = rMeshes[0].primitives[0].targets;
    int const tgtCount = targets.size();

    g_scene.m_tgt.m_targets.resize(tgtCount);
    for (int i = 0; i < tgtCount; i ++)
    {
        Targets::Target &rTgt = g_scene.m_tgt.m_targets[i];
        {
        auto &rVrtxAccessor = rAccessors[targets[i].at("POSITION")];
        g_scene.m_mesh.vertexCount = rVrtxAccessor.count;
        auto &rVrtxView = rViews[rVrtxAccessor.bufferView];
        rTgt.m_pPosIn = reinterpret_cast<glm::vec3*>(&rBuffers[rVrtxView.buffer].data[rVrtxView.byteOffset + rVrtxAccessor.byteOffset]);
        }
        {
        auto &rVrtxAccessor = rAccessors[targets[i].at("NORMAL")];
        g_scene.m_mesh.vertexCount = rVrtxAccessor.count;
        auto &rVrtxView = rViews[rVrtxAccessor.bufferView];
        rTgt.m_pNrmIn = reinterpret_cast<glm::vec3*>(&rBuffers[rVrtxView.buffer].data[rVrtxView.byteOffset + rVrtxAccessor.byteOffset]);
        }
    }
    }

    SetAudioStreamBufferSizeDefault(14400 / 60 * 6);
    SetConfigFlags(FLAG_WINDOW_TRANSPARENT);
    InitWindow(screenWidth, screenHeight, "Nice");

    g_scene.m_mat = LoadMaterialDefault();

    g_scene.m_ui = LoadRenderTexture(screenWidth, screenHeight);

    UploadMesh(&g_scene.m_mesh, true);


    g_scene.m_camera.target = Vector3{ 0.0f, 0.0f, 0.0f };
    g_scene.m_camera.up = Vector3{ 0.0f, 1.0f, 0.0f };
    g_scene.m_camera.fovy = 50.0f;
    g_scene.m_camera.projection = CAMERA_PERSPECTIVE;

#if defined(PLATFORM_WEB)
    emscripten_set_main_loop(UpdateDrawFrame, 60, 1);
#else
    

    SetTargetFPS(60);

    while (!WindowShouldClose())
    {
        update_draw_frame();
    }

#endif

    return 0;
}
