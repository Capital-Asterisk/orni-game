/**
 * SPDX-License-Identifier: MIT
 * SPDX-FileCopyrightText: 2022 Neal Nicdao <chrisnicdao0@gmail.com>
 */
#pragma once

#include <glm/mat4x4.hpp>

#include <vector>

namespace meshdeform
{

struct MeshJoints
{
    float const         *m_pWeightsIn;
    unsigned char const *m_pJointsIn;
};

struct Joints
{
    glm::mat4x4 const   *m_pInverseBindIn;

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

void calculate_joint_transforms(
        glm::mat4x4 const   &invWorldTf,
        glm::mat4x4 const   *pInverseBind,
        glm::mat4x4 const   *pNodeTf,
        std::size_t const   first,
        std::size_t const   last,
        glm::mat4x4         *pJointTfOut) noexcept;



void apply_vertex_transform(
        Joints const        &joints,
        MeshJoints const    &meshJoints,
        Targets const&      tgt,
        glm::vec3 const*    pPosIn,
        glm::vec3 const*    pNrmIn,
        std::size_t const   first,
        std::size_t const   last,
        glm::vec3           *pPosOut,
        glm::vec3           *pNrmOut) noexcept;



}
