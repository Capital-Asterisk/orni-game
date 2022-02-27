/**
 * SPDX-License-Identifier: MIT
 * SPDX-FileCopyrightText: 2022 Neal Nicdao <chrisnicdao0@gmail.com>
 */
#include "mesh_deform.hpp"

using namespace meshdeform;

void meshdeform::calculate_joint_transforms(
        glm::mat4x4 const   &invWorldTf,
        glm::mat4x4 const   *pInverseBind,
        glm::mat4x4 const   *pNodeTf,
        std::size_t const   first,
        std::size_t const   last,
        glm::mat4x4         *pJointTfOut) noexcept
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

void meshdeform::apply_vertex_transform(
        Joints const        &joints,
        MeshJoints const    &meshJoints,
        Targets const&      tgt,
        glm::vec3 const*    pPosIn,
        glm::vec3 const*    pNrmIn,
        std::size_t const   first,
        std::size_t const   last,
        glm::vec3           *pPosOut,
        glm::vec3           *pNrmOut) noexcept
{
    glm::vec3 *pPosOutCurr = &pPosOut[first];
    glm::vec3 *pNrmOutCurr = &pNrmOut[first];

    glm::vec3 const *pPosInCurr = &pPosIn[first];
    glm::vec3 const *pNrmInCurr = &pNrmIn[first];

    float const *pWeightCurr = &meshJoints.m_pWeightsIn[first * 4];
    unsigned char const *pJointCurr = &meshJoints.m_pJointsIn[first * 4];

    for (std::size_t i = first; i < last; i ++)
    {
        glm::vec3 pos(*pPosInCurr);
        glm::vec3 nrm(*pNrmInCurr);

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

            jointMatrix += weight * joints.m_jointTf[joint];
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
