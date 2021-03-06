#pragma once

#include "frogdyn.hpp"
#include "mesh_deform.hpp"

#include <raylib.h>

#include <tiny_gltf.h>

#include <unordered_map>
#include <atomic>
#include <thread>
#include <condition_variable>
#include <mutex>

namespace orni
{

struct ThreadedLoop
{
    std::atomic<bool>       m_running;
    std::thread             m_updater;
    unsigned long           m_framesUpdated{0};
    unsigned long           m_framesRendered{0};
    std::mutex              m_syncFramesMtx;
    std::condition_variable m_syncFramesCv;
    std::mutex              m_updaterBusy;
};

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
    float m_emPissed{0.0f};
    float m_emWorry{0.0f};
    float m_emBadIdea{0.0f};
    float m_emShy{0.0f};

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
        glm::vec2           m_texturePos;
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
    frog_id_t               m_frogTailTip;
};

using Characters_t = std::unordered_map<int, CharB>;

// gltf parsing

template<typename T>
struct Burger
{
    T const *m_data;
    std::size_t m_count;
};

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

using grab_id_t = int;

struct ToolGrab
{
    struct Grab
    {
        //salad_id_t          m_salad;
        //int                 m_triangle;
        //glm::vec2           m_barypos;
        //glm::vec3           m_pullTo;
        glm::vec3           m_pos;
        glm::vec2           m_screenPos;
        float               m_cntUpLastTouched;
        float               m_cntUpSelected; // negative when not selected
        int                 baitId;
        bool                m_visible;
        //bool                m_selected;
    };

    static constexpr float smc_radius{12};

    enum class FunnyMoments { None, Selected, Deselected };

    tool_id_t               m_id;

    std::vector<Grab>       m_grabs;
    int                     m_selected{-1};
    bool                    m_active;
    bool                    m_removeOnRelease;
};

struct GrabDisplay
{
    glm::vec2           m_screenPos;
    float               m_cntUpLastTouched;
    float               m_cntUpSelected;
    bool                m_visible;
    bool                m_selected;
};

struct ToolGrabRemover
{
    tool_id_t               m_id;
};

struct ToolGrabRotater
{
    tool_id_t               m_id;
};

inline bool g_limits{true};

bool update_tool_grab(
        Salads_t const& salads,
        WetJoints const& wet,
        FrogDyn& rFrogs,
        Inputs& rInputs,
        ToolGrab& rToolGrab);

void update_tool_grab_pos(
        Camera const& cam,
        Inputs const& inputs,
        FrogDyn& rFrogs,
        ToolGrab& rToolGrab,
        float delta);

bool update_tool_grab_rotate(ToolGrabRotater rotate, ToolGrab& rToolGrab, Inputs& rInputs, FrogDyn& rFrogs, Camera const& cam, float delta);

void update_grab_displays(ToolGrab const& grabs, FrogDyn const &frogs, std::vector<GrabDisplay> &rDisplay);

bool offset_camera_lazor(Salads_t const& salads, Inputs& rInputs, bool &rMouseMoved, glm::vec3 &rOffset, glm::vec3 com);
int paw_default_base_attribute(meshdeform::MeshJoints const& joints, unsigned short const *pInd);

void update_expressions(Soul &rSoul, float delta);

inline float breath_cycle(float t) noexcept
{
    float const tau = glm::pi<float>() * 2.0f;
    return glm::sin(tau*(t - 0.25f)) + 0.2f * sin(tau*2*t) + 1.0;
}

inline float rand_dist(float dist) noexcept
{
    float woot = GetRandomValue(-65536, 65536) / 65536.0f;
    return woot * woot * glm::sign(woot) * dist;
}

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


McRay shoop_da_whoop(glm::vec3 origin, glm::vec3 dir, int triCount, glm::vec3 const* pVrt, unsigned short const* pInd);

McRay shoop_da_woop_salad(glm::vec3 origin, glm::vec3 dir, SaladModel const& salad);

McRaySalad lazor_salads(glm::vec3 origin, glm::vec3 dir, Salads_t const& salads);

void update_apples(Apples &rApples, meshdeform::Joints const& rJoints);

void update_hoppers(std::vector<WetJoints::Hopper> const& hoppers, FrogDyn const& frogs, glm::mat4x4 *pNodeTf);

glm::vec2 calc_eye_pos(glm::mat4x4 const& eyeTf, glm::vec3 tgt);

bool eye_visible(glm::mat4x4 const& eyeTf, glm::vec3 tgt) noexcept;

void draw_iris(Texture2D texture, int i, glm::vec2 pos);

void update_inputs_rl(Camera const& cam, Inputs& rInputs);



}
