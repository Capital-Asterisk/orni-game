#pragma once

#include "frogdyn.hpp"
#include "mesh_deform.hpp"

#include <raylib.h>

#include <tiny_gltf.h>

#include <unordered_map>

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
    Mesh                    m_rayMesh;
    Model                   m_rayModel;
};

using Salads_t = std::vector< std::unique_ptr<SaladModel> >;

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


struct Soul
{
    float m_blinkPeriodAvg{6.0f};
    float m_blinkPeriodMargin{2.0f};
    float m_blinkCdn{0.0f};

    float m_breathSpeed{0.35f};
    float m_breathCycle{0.0f};
};

struct CharB
{
    struct Eye
    {
        apple_id_t          m_apple;
        glm::ivec2          m_irisPos;
    };

    Soul                    m_soul;

    Material                m_eyeMaterial;
    Texture                 m_eyeSheet;
    RenderTexture           m_eyeTexture;
    Eye                     m_eyeL;
    Eye                     m_eyeR;

    std::vector<std::string> m_jointNames;
    meshdeform::Joints      m_joints;
    WetJoints               m_wetJoints;
    Apples                  m_apples;

    frog_id_t               m_frogBelly;
    frog_id_t               m_frogBeak;
    frog_id_t               m_frogHead;
};

using Characters_t = std::unordered_map<int, CharB>;

// gltf parsing

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
        CharB& rChar,
        glm::mat4x4 parentTfWorld);

void metal_bar(
        tinygltf::Model const&      gltf,
        int                         nodeId,
        CharB&                      rChar,
        Salads_t&                   rSalads,
        std::vector<WetJoints> &    rSpooks,
        std::vector<Material>&      rMaterials,
        FrogDyn&                    rFrogs);

void metal_pipe(
        tinygltf::Model const&      gltf,
        int                         sceneId,
        Characters_t&               rChars,
        Salads_t&                   rSalads,
        std::vector<Material>&      rMaterials,
        FrogDyn&                    rFrogs);

// Mesh lazoring

struct McRay
{
    glm::vec2   m_barypos;
    float       m_dist;
    int         m_index;
};

struct McRaySalad
{
    McRay           m_mcray;
    salad_id_t      m_salad;
};

McRay shoop_da_whoop(glm::vec3 origin, glm::vec3 dir, int triCount, glm::vec3 const* pVrt, unsigned short const* pInd);

McRay shoop_da_woop_salad(glm::vec3 origin, glm::vec3 dir, SaladModel const& salad);

// tools

using tool_id_t = int;

struct Inputs
{
    tool_id_t               m_selected;

    glm::vec2               m_mousePos;
    glm::vec3               m_mouseOrig;
    glm::vec3               m_mouseDir;
//    bool                    m_mouseLPrev;
//    bool                    m_mouseRPrev;
//    bool                    m_mouseL;
//    bool                    m_mouseR;

    McRaySalad              m_lazor;
};

struct ToolGrab
{
    struct Grab
    {
        //salad_id_t          m_salad;
        //int                 m_triangle;
        //glm::vec2           m_barypos;
        //glm::vec3           m_pullTo;
        int baitId;
    };

    tool_id_t               m_id;

    std::vector<Grab>       m_grabs;
    bool                    m_active;
};

}
