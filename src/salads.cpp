#include "salads.hpp"

#include <glm/gtx/transform.hpp>
#include <glm/gtx/intersect.hpp>
#include <glm/gtx/vector_angle.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <raymath.h>

#include <iostream>

namespace orni
{

bool update_tool_grab(
        Salads_t const& salads,
        WetJoints const& wet,
        FrogDyn& rFrogs,
        Inputs& rInputs,
        ToolGrab& rToolGrab)
{
    bool const toolSelected = rInputs.m_selected == rToolGrab.m_id;
    bool clickedSomething = false;

    if (rToolGrab.m_active)
    {
        if (!IsMouseButtonDown(MOUSE_BUTTON_LEFT))
        {
            // release
            rToolGrab.m_active = false;
            if (rToolGrab.m_removeOnRelease)
            {
                int id = rToolGrab.m_grabs.back().baitId;
                rToolGrab.m_grabs.pop_back();
                auto found = std::find_if(rFrogs.m_baits.begin(), rFrogs.m_baits.end(), [id] (FrogDyn::Bait const& bait) { return bait.m_id == id; });
                rFrogs.m_baits.erase(found);
                clickedSomething = true;
            }
        }
        else if (IsKeyPressed(KEY_SPACE))
        {
            rToolGrab.m_active = false;
        }
    }
    else if (toolSelected && IsMouseButtonPressed(MOUSE_BUTTON_RIGHT) && (rToolGrab.m_selected != -1))
    {
        if (!g_limits || rToolGrab.m_grabs.size() != 1)
        {
            // copy paste other parts of code for remove!
            std::iter_swap(rToolGrab.m_grabs.begin() + rToolGrab.m_selected, rToolGrab.m_grabs.end() - 1);
            int id = rToolGrab.m_grabs.back().baitId;
            rToolGrab.m_grabs.pop_back();
            auto found = std::find_if(rFrogs.m_baits.begin(), rFrogs.m_baits.end(), [id] (FrogDyn::Bait const& bait) { return bait.m_id == id; });
            rFrogs.m_baits.erase(found);
            clickedSomething = true;
        }
    }
    else if (toolSelected && IsMouseButtonPressed(MOUSE_BUTTON_LEFT))
    {
        // do grab
        if (rToolGrab.m_selected != -1)
        {
            // something selected, drag it
            std::iter_swap(rToolGrab.m_grabs.begin() + rToolGrab.m_selected, rToolGrab.m_grabs.end() - 1);
            rToolGrab.m_active = true;
            rToolGrab.m_removeOnRelease = false;
            clickedSomething = true;
        }
        else
        {
            if (rInputs.m_lazor.m_salad == -1)
            {
                rInputs.m_lazor = lazor_salads(rInputs.m_mouseOrig, rInputs.m_mouseDir, salads);
            }

            if (rInputs.m_lazor.m_salad != -1)
            {
                SaladModel const& salad = *salads.at(rInputs.m_lazor.m_salad);

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
                    float const scale = rFrogs.m_scale[frog];
                    glm::vec3 const frogInvHit = glm::inverse(rFrogs.m_tf[frog]) * glm::vec4(lazorHit, 1.0f);

                    rToolGrab.m_active = true;
                    rToolGrab.m_removeOnRelease = true;
                    clickedSomething = true;
                    rGrab.baitId = rFrogs.m_baitNextId++;

                    FrogDyn::Bait &rBait = rFrogs.m_baits.emplace_back(
                                FrogDyn::Bait{ {-1,  lazorHit, down, fwd},
                                               {frog,   frogInvHit / scale, down, fwd} });
                    rBait.m_id = rGrab.baitId;
                    if (g_limits)
                    {
                        rBait.m_forceLim = 15.2f * 10.0f * 4.0f;
                    }
                    //rBait.m_strengthMul = 0.2f;

                    //rFrogs.m_extImp[frog].m_lin += rInputs.m_mouseDir * 10.0f;
                }
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

    return clickedSomething;
}


void update_tool_grab_pos(
        Camera const& cam,
        Inputs const& inputs,
        FrogDyn& rFrogs,
        ToolGrab& rToolGrab,
        float delta)
{
    glm::vec3 camPos = reinterpret_cast<glm::vec3 const&>(cam.position);
    glm::vec3 camDir = glm::normalize(reinterpret_cast<glm::vec3 const&>(cam.target) - camPos);

    int selected = -1;
    for (int i = 0; i < rToolGrab.m_grabs.size(); ++i)
    {
        ToolGrab::Grab &rGrab = rToolGrab.m_grabs[i];
        // yes this is O(n*m) with questionable cache locality
        // still faster than python lol
        for (FrogDyn::Bait const &bait : rFrogs.m_baits)
        {
            if (bait.m_id == rGrab.baitId)
            {
                rGrab.m_pos = bait.m_a.m_offset;
                break;
            }
        }

        rGrab.m_visible = glm::dot(rGrab.m_pos - camPos, camDir) > 0.0f;
        Vector2 pos = GetWorldToScreen(reinterpret_cast<Vector3&>(rGrab.m_pos), cam);
        rGrab.m_screenPos = reinterpret_cast<glm::vec2&>(pos);

        if (selected == -1 && rToolGrab.m_id == inputs.m_selected)
        {
            if (glm::distance(rGrab.m_screenPos, inputs.m_mousePos) < ToolGrab::smc_radius)
            {
                selected = i;
            }
        }

        rGrab.m_cntUpLastTouched += delta;
        rGrab.m_cntUpSelected += delta;
    }
    rToolGrab.m_selected = selected;
}

bool update_tool_grab_rotate(ToolGrabRotater rotate, ToolGrab& rToolGrab, Inputs& rInputs, FrogDyn& rFrogs, Camera const& cam, float delta)
{
    bool const toolSelected = rInputs.m_selected == rotate.m_id;

    if (toolSelected && IsMouseButtonDown(MOUSE_BUTTON_LEFT))
    {
        glm::vec3 const camPos = reinterpret_cast<glm::vec3 const&>(cam.position);
        glm::vec3 const camTgt = reinterpret_cast<glm::vec3 const&>(cam.target);
        glm::vec3 const camUp = reinterpret_cast<glm::vec3 const&>(cam.up);
        glm::vec3 const camDir = glm::normalize(camTgt - camPos);
        glm::vec3 const pitchAxis = glm::normalize(glm::cross(camDir, camUp));
        glm::vec3 const yawAxis = glm::normalize(glm::cross(camDir, pitchAxis));

        float const sensitivity = 0.25f;

        Vector2 const mouseDelta = GetMouseDelta();
        float pitch = mouseDelta.y * delta * sensitivity;
        float yaw = -mouseDelta.x * delta * sensitivity;

        glm::mat4x4 rotation = glm::rotate(pitch, pitchAxis) * glm::rotate(yaw, yawAxis);

        for (ToolGrab::Grab &rGrab : rToolGrab.m_grabs)
        {
            // same thing here
            for (FrogDyn::Bait &rBait : rFrogs.m_baits)
            {
                if (rBait.m_id == rGrab.baitId)
                {
                    rBait.m_a.m_offset = glm::vec3(rotation * glm::vec4(rBait.m_a.m_offset - camTgt, 1.0f)) + camTgt;
                    break;
                }
            }
        }
        return true;
    }

    return false;
}

void update_grab_displays(ToolGrab const& grabs, FrogDyn const &frogs, std::vector<GrabDisplay> &rDisplay)
{
    int size = grabs.m_grabs.size();
    rDisplay.resize(size);

    for (int i = 0; i < size; i ++)
    {
        ToolGrab::Grab const& grab = grabs.m_grabs[i];
        rDisplay[i] = {grab.m_screenPos, grab.m_cntUpLastTouched, grab.m_cntUpSelected, grab.m_visible, i == grabs.m_selected};
    }
}

bool offset_camera_lazor(Salads_t const& salads, Inputs& rInputs, bool &rMouseMoved, glm::vec3 &rOffset, glm::vec3 com)
{

    if (IsMouseButtonPressed(MOUSE_BUTTON_RIGHT))
    {
        rMouseMoved = false;
    }
    if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT))
    {
        Vector2 mouseDelta = GetMouseDelta();
        rMouseMoved |= (mouseDelta.x != 0.0f || mouseDelta.y != 0.0f);
    }

    if (IsMouseButtonReleased(MOUSE_BUTTON_RIGHT) && !rMouseMoved)
    {
        if (rInputs.m_lazor.m_salad == -1)
        {
            rInputs.m_lazor = lazor_salads(rInputs.m_mouseOrig, rInputs.m_mouseDir, salads);
        }

        if (rInputs.m_lazor.m_salad != -1)
        {
            glm::vec3 const lazorHit = rInputs.m_mouseOrig + rInputs.m_mouseDir * rInputs.m_lazor.m_mcray.m_dist;
            rOffset = lazorHit - com;
            return true;
        }
        else
        {
            rOffset = glm::vec3{0.0f};
        }
    }
    return false;
}

int paw_default_base_attribute(meshdeform::MeshJoints const& joints, unsigned short const *pInd)
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

glm::mat4 node_transform(tinygltf::Node const& node)
{
    glm::vec3 const scale = (node.scale.size() == 3) ? glm::vec3{node.scale[0], node.scale[1], node.scale[2]} : glm::vec3{1.0f};
    glm::vec3 const translate = (node.translation.size() == 3) ? glm::vec3{node.translation[0], node.translation[1], node.translation[2]} : glm::vec3{0.0f};
    glm::quat const rot = (node.rotation.size() == 4) ? glm::quat(node.rotation[3], node.rotation[0], node.rotation[1], node.rotation[2]) : glm::quat{};
    return glm::translate(translate) * glm::mat4(rot) * glm::scale(scale);
}

// write some templates and never use them
template<typename T>
static T crab(tinygltf::Value const& value, std::string const& name, T&& def)
{
    tinygltf::Type valType = tinygltf::Type::NULL_TYPE;
    if constexpr (std::is_same_v<T, double>)
    {
        valType = tinygltf::Type::REAL_TYPE;
    }
    else if constexpr (std::is_same_v<T, std::string>)
    {
        valType = tinygltf::Type::STRING_TYPE;
    }
    // add more lol

    if (tinygltf::Value const& patty = value.Get(name); patty.Type() == valType)
    {
        return patty.Get<T>();
    }
    else
    {
        return def;
    }
}

static int zero_workaround_lol = 0;

void metal_rod(
        std::vector<int> const& nodeToJoint,
        tinygltf::Model const& gltf,
        int nodeId,
        int parentId,
        int level,
        FrogDyn& rFrogs,
        CharB& rChar,
        glm::mat4x4 parentTfWorld)
{
    auto const &node = gltf.nodes[nodeId];

    // record name
    if (level != -1)
    {
        rChar.m_jointNames.at(nodeToJoint[nodeId]) = node.name;
    }

    glm::mat4x4 const nodeTf = node_transform(node);
    glm::mat4x4 const nodeTfWorld = parentTfWorld * nodeTf;

    frog_id_t frogId = lgrn::id_null<frog_id_t>();
    float comwok = 0.12f;
    glm::vec3 offset{0.0f, -comwok, 0.0f}; // bone length estimation lol

    // check for mass
    if (node.extras.IsObject()) // sure
    {


        if (tinygltf::Value const& massVal = node.extras.Get("mass"); massVal.IsNumber())
        {
            float mass = massVal.GetNumberAsDouble();
            //std::cout << "MASSSS! " << mass << "\n";

            // Make a frog
            glm::mat4x4 frogTfWorld     = nodeTfWorld;
            glm::mat4x4 const frogTfInv = glm::inverse(frogTfWorld);
            frogTfWorld[3] += nodeTfWorld[1] * comwok;
            frogId = frogdyn::add_frog(rFrogs, frogTfWorld, mass);

            std::cout << "Bone[" << node.name << "] is Frog[" << frogId << "] with mass: " << mass << "\n";

            // Add hopper
            rChar.m_wetJoints.m_hoppers.emplace_back(WetJoints::Hopper{frogId, nodeToJoint[nodeId], -comwok});

            // lazy record important frogs
            if (node.name == "belly")
            {
                rChar.m_frogBelly = frogId;
            }
            else if (node.name == "beak")
            {
                rChar.m_frogBeak = frogId;
            }
            else if (node.name == "head")
            {
                rChar.m_frogHead = frogId;
            }
            else if (node.name == "tail3")
            {
                rChar.m_frogTailTip = frogId;
            }

            // connect to parent
            if (level != 0)
            {
                int const parentJoint = nodeToJoint[parentId];
                assert(parentJoint != -1);
                auto foundIt = std::find_if(rChar.m_wetJoints.m_hoppers.begin(), rChar.m_wetJoints.m_hoppers.end(), [parentJoint] (WetJoints::Hopper const &wet) -> bool { return wet.m_joint == parentJoint; });
                assert(foundIt != rChar.m_wetJoints.m_hoppers.end());

                frog_id_t const     parentFrogId = foundIt->m_joint;

                // Create bait
                FrogDyn::Bait       &rBait      = rFrogs.m_baits.emplace_back();

                glm::vec4 const     side        {glm::vec3{nodeTfWorld[0]}, 0.0f};
                glm::vec4 const     dir         {glm::vec3{nodeTfWorld[1]}, 0.0f}; // Y is forward for blender

                glm::mat4x4 const   parentTfInv = glm::inverse(parentTfWorld);

                rBait.m_a = {parentFrogId, glm::vec3{nodeTf[3]} + offset, parentTfInv * dir, parentTfInv * side};
                rBait.m_b = {frogId, offset, frogTfInv * dir, frogTfInv * side};

                rBait.m_angDrag = 3.0f;

                if (tinygltf::Value const& patty = node.extras.Get("tlim"); patty.IsNumber())
                {
                    rBait.m_doTwistLim = true;
                    rBait.m_twistRange = glm::radians(patty.GetNumberAsDouble());
                }
                if (tinygltf::Value const& patty = node.extras.Get("tspr"); patty.IsNumber())
                {
                    rBait.m_doTwistSpring = true;
                    rBait.m_twistSpring = patty.GetNumberAsDouble();
                }
                if (tinygltf::Value const& patty = node.extras.Get("clim"); patty.IsNumber())
                {
                    rBait.m_doConeLim = true;
                    rBait.m_coneRange = glm::max(glm::radians(patty.GetNumberAsDouble()), 0.0125);
                }
                if (tinygltf::Value const& patty = node.extras.Get("cspr"); patty.IsNumber())
                {
                    rBait.m_doConeSpring = true;
                    rBait.m_coneSpring = patty.GetNumberAsDouble();
                }
                if (tinygltf::Value const& patty = node.extras.Get("align"); patty.IsString())
                {
                    if (patty.Get<std::string>() == "x")
                    {
                        rBait.m_doAlign = true;
                    }
                    else
                    {
                        assert(0); // not yet implemented!
                    }
                }
            }
        }
    }

    std::vector<FrogDyn::Ball> balls;

    for (int childId : node.children)
    {
        auto const &child = gltf.nodes.at(childId);

        if (nodeToJoint[childId] == -1)
        {
            // Not bone

            if (frogId != lgrn::id_null<frog_id_t>() && child.name.rfind("Ball") == 0)
            {
                std::cout << "BALL!\n";
                glm::mat4x4 const ballTf = node_transform(child);
                balls.emplace_back(FrogDyn::Ball{glm::vec3(ballTf[3]) + offset, glm::length(glm::vec3(ballTf[1])), 0, 0});
            }

            if (child.name.rfind("Apl_EyeL") == 0)
            {
                rChar.m_eyeL.m_apple = rChar.m_apples.create(Apple{ node_transform(child), nodeToJoint[nodeId] });
            }
            else if (child.name.rfind("Apl_EyeR") == 0)
            {
                rChar.m_eyeR.m_apple = rChar.m_apples.create(Apple{ node_transform(child), nodeToJoint[nodeId] });
            }
        }
        else
        {
            // Bone

            //std::cout << child.name << "\n";

            metal_rod(nodeToJoint, gltf, childId, nodeId, level + 1, rFrogs, rChar, nodeTfWorld);
        }
    }

    if (!balls.empty())
    {
        // balls contain pee
        rFrogs.m_balls.emplace(frogId, balls.begin(), balls.end());
        rFrogs.m_ballPos.emplace(frogId, balls.size());
        rFrogs.m_canCollide.push_back(frogId);
    }
}

void metal_bar(
        tinygltf::Model const&      gltf,
        int                         nodeId,
        CharB&                      rChar,
        Salads_t&                   rSalads,
        std::vector<Material>&      rMaterials,
        FrogDyn&                    rFrogs)
{
    auto const &node = gltf.nodes.at(nodeId);

    int useSkin = -1;

    // init eyes
    rChar.m_eyeSheet = LoadTexture("resources/textures/eyessheet.png");
    rChar.m_eyeTexture = LoadRenderTexture(256, 256);
    rChar.m_eyeMaterial = LoadMaterialDefault();
    SetMaterialTexture(&rChar.m_eyeMaterial, MATERIAL_MAP_DIFFUSE, rChar.m_eyeTexture.texture);

    // load mesh children
    for (int childId : node.children)
    {
        auto const &child = gltf.nodes.at(childId);

        if (child.mesh != -1)
        {
            SaladModel &rSalad = *rSalads.emplace_back(std::make_unique<SaladModel>());

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
    rChar.m_jointNames.resize(jointCount);
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

    rFrogs.m_balls.resize_ids(rFrogs.m_ids.capacity() + jointCount);
    rFrogs.m_ballPos.resize_ids(rFrogs.m_ids.capacity() + jointCount);
    rFrogs.m_balls.resize_data(gltf.nodes.size()); // ultra-safe assumption
    rFrogs.m_ballPos.resize_data(gltf.nodes.size());

    // Iterate bones
    metal_rod(nodeToJoint, gltf, nodeId, -1, -1, rFrogs, rChar, glm::mat4x4(1.0f));
}

void metal_pipe(
        tinygltf::Model const&      gltf,
        int                         sceneId,
        Characters_t&               rChars,
        Salads_t&                   rSalads,
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
            metal_bar(gltf, nodeId, rChar, rSalads, rMaterials, rFrogs);
        }
    }
}

McRay shoop_da_whoop(glm::vec3 origin, glm::vec3 dir, int triCount, glm::vec3 const* pVrt, unsigned short const* pInd)
{
    glm::vec2   outBarypos;
    float       outDist     = std::numeric_limits<float>::max();
    int         outIndex    = -1;

    unsigned short const *pIndCurr = pInd;

    for (int i = 0; i < triCount; i ++)
    {
        glm::vec2 barypos;
        float dist;
        if (glm::intersectRayTriangle<float, glm::defaultp>(origin, dir, pVrt[pIndCurr[0]], pVrt[pIndCurr[1]], pVrt[pIndCurr[2]], barypos, dist))
        {
            if (outDist > dist)
            {
                outDist = dist;
                outIndex = i;
            }
        }

        pIndCurr += 3;
    }

    return {outBarypos, outDist, outIndex};
}

McRay shoop_da_woop_salad(glm::vec3 origin, glm::vec3 dir, SaladModel const& salad)
{
    return shoop_da_whoop(origin, dir, salad.m_rayMesh.triangleCount, salad.m_Pos.data(), salad.m_rayMesh.indices);
}

McRaySalad lazor_salads(glm::vec3 origin, glm::vec3 dir, Salads_t const& salads)
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

void update_expressions(Soul &rSoul, float delta)
{
    rSoul.m_blinkCdn -= delta;

    if (rSoul.m_blinkCdn <= 0.0f)
    {
        rSoul.m_blinkCdn = rand_dist(rSoul.m_blinkPeriodMargin) + rSoul.m_blinkPeriodAvg;
    }

    rSoul.m_breathCycle += delta * rSoul.m_breathSpeed;
    rSoul.m_breathCycle -= float(rSoul.m_breathCycle > 1.0f);
}

void update_apples(Apples &rApples, meshdeform::Joints const& rJoints)
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

void update_hoppers(std::vector<WetJoints::Hopper> const& hoppers, FrogDyn const& frogs, glm::mat4x4 *pNodeTf)
{
    // Update hoppers
    for (WetJoints::Hopper const& hopper : hoppers)
    {
        float const scale = frogs.m_scale[hopper.m_frog];
        pNodeTf[hopper.m_joint] = glm::translate(frogs.m_tf.at(hopper.m_frog) * glm::scale(glm::vec3{scale, scale, scale}), glm::vec3{0, hopper.m_yoffset, 0});
    }

}

// lol
float const eyeDepth = 0.03f;
float const eyeRadius = 0.09f;
float const eyeFov = glm::atan(eyeRadius / eyeDepth) * 2.0f;

glm::vec2 calc_eye_pos(glm::mat4x4 const& eyeTf, glm::vec3 tgt)
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

bool eye_visible(glm::mat4x4 const& eyeTf, glm::vec3 tgt) noexcept
{
    float const ang = glm::angle(glm::vec3(eyeTf[2]), tgt - glm::vec3(eyeTf[3]));
    return ang < eyeFov;
}

void draw_iris(Texture2D texture, int i, glm::vec2 pos)
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

void update_inputs_rl(Camera const& cam, Inputs& rInputs)
{
    Vector2 mousePos = GetMousePosition();
    Ray ray = GetMouseRay(mousePos, cam);

    rInputs.m_mousePos      = {mousePos.x, mousePos.y};
    rInputs.m_mouseOrig     = {ray.position.x, ray.position.y, ray.position.z};
    rInputs.m_mouseDir      = {ray.direction.x, ray.direction.y, ray.direction.z};
}

}
