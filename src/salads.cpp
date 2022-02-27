#include "salads.hpp"

#include <glm/gtx/transform.hpp>
#include <glm/gtx/intersect.hpp>

#include <raymath.h>

#include <iostream>

namespace orni
{

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

}
