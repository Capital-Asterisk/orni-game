#pragma once

#include "frogdyn.hpp"
#include "mesh_deform.hpp"

#include <raylib.h>

#include <tiny_gltf.h>

namespace orni
{

// LOL
using namespace frogdyn;
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

glm::mat4 node_transform(tinygltf::Node const& node);

void metal_rod(
        std::vector<int> const& nodeToJoint,
        tinygltf::Model const& gltf,
        int nodeId,
        int parentId,
        int level,
        FrogDyn& rFrogs,
        Apples& rApples,
        CharB& rChar,
        glm::mat4x4 parentTfWorld);

void metal_bar(
        tinygltf::Model const&      gltf,
        int                         nodeId,
        CharB&                      rChar,
        std::vector<SaladModel>&    rSalads,
        Apples&                     rApples,
        std::vector<WetJoints> &    rSpooks,
        std::vector<Material>&      rMaterials,
        FrogDyn&                    rFrogs);

void metal_pipe(
        tinygltf::Model const&      gltf,
        int                         sceneId,
        Characters_t&               rChars,
        std::vector<SaladModel>&    rSalads,
        Apples&                     rApples,
        std::vector<WetJoints> &    rSpooks,
        std::vector<Material>&      rMaterials,
        FrogDyn&                    rFrogs);

}
