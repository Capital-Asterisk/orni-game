/**
 * SPDX-License-Identifier: MIT
 * SPDX-FileCopyrightText: 2022 Neal Nicdao <chrisnicdao0@gmail.com>
 */
#include "frogdyn.hpp"

#include <glm/gtx/norm.hpp>
#include <glm/gtx/matrix_interpolation.hpp>

using namespace frogdyn;

void frogdyn::apply_baits(FrogDyn &rDyn)
{
    // calculate baits
    for (int i = 0; i < rDyn.m_baits.size(); i ++)
    {
        auto const &rBait = rDyn.m_baits[i];
        glm::vec3 posA, velA, posB, velB, offsetRotA, offsetRotB;
        float massA, massB;

        if (rBait.m_a.m_id == -1)
        {
            // world-anchored
            offsetRotA = rBait.m_a.m_offset;
            posA = rBait.m_a.m_offset;
            velA = glm::vec3{0.0f};
            massA = 9999999.0f;
        }
        else
        {
            offsetRotA = rDyn.m_tf[rBait.m_a.m_id] * glm::vec4(rBait.m_a.m_offset, 0.0f);
            posA = offsetRotA + glm::vec3(rDyn.m_tf[rBait.m_a.m_id][3]);

            glm::vec3 const velTan = glm::cross(rDyn.m_vel[rBait.m_a.m_id].m_ang, offsetRotA);
            velA = rDyn.m_vel[rBait.m_a.m_id].m_lin + velTan;
            massA = rDyn.m_mass[rBait.m_a.m_id];
        }

        offsetRotB = (rDyn.m_tf[rBait.m_b.m_id] * glm::vec4(rBait.m_b.m_offset, 0.0f));
        posB = offsetRotB + glm::vec3(rDyn.m_tf[rBait.m_b.m_id][3]);

        glm::vec3 const velTan = glm::cross(rDyn.m_vel[rBait.m_b.m_id].m_ang, offsetRotB);
        velB = rDyn.m_vel[rBait.m_b.m_id].m_lin + velTan;

        massB = rDyn.m_mass[rBait.m_b.m_id];



        glm::vec3 const     posRel  = posB - posA;
        glm::vec3 const     velRel  = velB - velA;

        // PD control lol what
        glm::vec3 const     p   = posRel * 32.0f;
        glm::vec3 const     d   = velRel * 0.5f;

        float const lenBSq = glm::length2(rBait.m_b.m_offset);
        float const lenASq = glm::length2(rBait.m_a.m_offset);

        float const minMass = glm::min(massA, massB);
        float const maxLen  = glm::max(lenASq, lenBSq);

        glm::vec3 const totalLinImp = (p + d) * minMass;
        glm::vec3 const totalAngImp = totalLinImp * 1.1f / (maxLen + 1.0f);

        if (rBait.m_a.m_id != -1)
        {
            rDyn.m_cstImp[rBait.m_a.m_id].m_lin += totalLinImp;
            rDyn.m_cstImp[rBait.m_a.m_id].m_ang += glm::cross(offsetRotA, totalAngImp);

        }

        rDyn.m_cstImp[rBait.m_b.m_id].m_lin -= totalLinImp;
        rDyn.m_cstImp[rBait.m_b.m_id].m_ang -= glm::cross(offsetRotB, totalAngImp);
    }
}

void frogdyn::apply_ext_forces(FrogDyn &rDyn, float delta)
{
    // apply external forces
    for (frog_id_t id = 0; id < rDyn.m_ids.capacity(); id ++)
    {
        if (!rDyn.m_ids.exists(id))
        {
            continue;
        }

        rDyn.m_vel[id].m_lin += (rDyn.m_extImp[id].m_lin) / rDyn.m_mass[id];
        rDyn.m_vel[id].m_ang += (rDyn.m_extImp[id].m_ang) / rDyn.m_mass[id];

        rDyn.m_tf[id][3] += glm::vec4(rDyn.m_vel[id].m_lin, 0.0f) * delta;
        auto const& w = rDyn.m_vel[id].m_ang;

        float const mag = glm::length(w);
        if (mag > 0.001f)
        {
            auto const saved = rDyn.m_tf[id][3];
            rDyn.m_tf[id][3] = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
            rDyn.m_tf[id] = glm::axisAngleMatrix(w / mag, mag * delta) * rDyn.m_tf[id];
            rDyn.m_tf[id][3] = saved;
        }

        rDyn.m_extImp[id].m_ang = glm::vec3(0.0f);
        rDyn.m_extImp[id].m_lin = glm::vec3(0.0f);
    }
}

void frogdyn::apply_cst_forces(FrogDyn &rDyn, float delta)
{
    // apply external forces
    for (frog_id_t id = 0; id < rDyn.m_ids.capacity(); id ++)
    {
        if (!rDyn.m_ids.exists(id))
        {
            continue;
        }

        rDyn.m_vel[id].m_lin += (rDyn.m_cstImp[id].m_lin) / rDyn.m_mass[id];
        rDyn.m_vel[id].m_ang += (rDyn.m_cstImp[id].m_ang) / rDyn.m_mass[id];

        rDyn.m_tf[id][3] += glm::vec4(rDyn.m_vel[id].m_lin, 0.0f) * delta;
        auto const& w = rDyn.m_vel[id].m_ang;

        float const mag = glm::length(w);
        if (mag > 0.001f)
        {
            auto const saved = rDyn.m_tf[id][3];
            rDyn.m_tf[id][3] = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
            rDyn.m_tf[id] = glm::axisAngleMatrix(w / mag, mag * delta) * rDyn.m_tf[id];
            rDyn.m_tf[id][3] = saved;
        }

        rDyn.m_cstImp[id].m_ang = glm::vec3(0.0f);
        rDyn.m_cstImp[id].m_lin = glm::vec3(0.0f);
    }
}

void frogdyn::calc_balls_pos(FrogDyn &rDyn)
{
    for (frog_id_t id : rDyn.m_canCollide)
    {
        auto const& tf      = rDyn.m_tf[id];
        auto &rAABB         = rDyn.m_aabb[id];
        auto const ballsIn  = rDyn.m_balls[id];
        auto posOut         = rDyn.m_ballPos[id];

        rAABB.m_min = glm::vec3{9999, 9999, 9999};
        rAABB.m_max = glm::vec3{-9999, -9999, -9999};

        for (unsigned int i = 0; i < ballsIn.size(); i ++)
        {
            auto const &ball = ballsIn[i];
            auto &rOut = posOut[i];
            rOut.m_pos = tf * glm::vec4(ball.m_pos, 1.0f);
            rOut.m_radius = ball.m_radius * rDyn.m_scale[id];

            rAABB.m_max = glm::max(rOut.m_pos + rOut.m_radius, rAABB.m_max);
            rAABB.m_min = glm::min(rOut.m_pos - rOut.m_radius, rAABB.m_min);
        }
    }
}

BallContact frogdyn::calc_ball_collisions(BallCollisions_t a, BallCollisions_t b)
{
    BallContact contact;
    // test each ball A against each ballB
    for (auto const& ballA : a)
    {
        for (auto const& ballB : b)
        {
            float const radC = ballA.m_radius + ballB.m_radius;
            float const distSq = glm::distance2(ballA.m_pos, ballB.m_pos);
            if (distSq < (radC * radC))
            {
                float const dist = glm::sqrt(distSq);

                contact.m_count ++;
                contact.m_depth += radC - dist;
                contact.m_pos += (ballA.m_pos * ballB.m_radius + ballB.m_pos * ballA.m_radius) / radC;
                contact.m_nrm += (ballB.m_pos - ballA.m_pos) / glm::sqrt(distSq);
            }
        }
    }
    return contact;
}

void frogdyn::calc_frog_collisions(FrogDyn &rDyn)
{
    //   0  1  2  3  4  5
    // 0
    // 1 x
    // 2 x  x
    // 3 x  x  x
    // 4 x  x  x  x
    // 5 x  x  x  x  x

    for (int i = 1; i < rDyn.m_canCollide.size(); i ++)
    {
        for (int j = 0; j < i; j ++)
        {
            frog_id_t const a = rDyn.m_canCollide[i];
            frog_id_t const b = rDyn.m_canCollide[j];

            if (aabb_intersect(rDyn.m_aabb[a], rDyn.m_aabb[b]))
            {
                // AABBs intersect, now start checking balls
                BallContact contact = calc_ball_collisions(rDyn.m_ballPos[a], rDyn.m_ballPos[b]);
                if (contact.m_count != 0)
                {
                    glm::vec3 const pos = contact.m_pos / float(contact.m_count);
                    glm::vec3 const nrm = contact.m_nrm / float(contact.m_count);

                    glm::vec3 const posA = pos - glm::vec3(rDyn.m_tf[a][3]);
                    glm::vec3 const posB = pos - glm::vec3(rDyn.m_tf[b][3]);

                    //glm::vec3 const velTanA = glm::cross(rDyn.m_vel[a].m_ang, posA);
                    //glm::vec3 const velTanB = glm::cross(rDyn.m_vel[b].m_ang, posB);

                    //rDyn.m_vel[b].m_lin + velTanB - rDyn.m_vel[a].m_lin - velTanA;
                    glm::vec3 const velRel = rDyn.m_vel[b].m_lin - rDyn.m_vel[a].m_lin;

                    float const hiteachotherness = glm::dot(velRel, nrm);
                    glm::vec3 const velTan = velRel - nrm * hiteachotherness;

                    float const minMass = glm::min(rDyn.m_mass[a], rDyn.m_mass[b]);


                    glm::vec3 const nrmImp = nrm * glm::min(0.0f, glm::dot(velRel, nrm) - contact.m_depth) * 0.5f;


                    //float const friction = 0.0f;
                    //glm::vec3 const totalImp = (nrmImp + velTan * friction) * minMass;
                    glm::vec3 const totalImp = nrmImp * minMass;

                    rDyn.m_cstImp[a].m_lin += totalImp;
                    rDyn.m_cstImp[a].m_ang += glm::cross(posA, totalImp);

                    rDyn.m_cstImp[b].m_lin -= totalImp;
                    rDyn.m_cstImp[b].m_ang -= glm::cross(posB, totalImp);


                    glm::vec3 const posnrm = pos + velTan * 1.0f;


                    //DrawSphereWires(reinterpret_cast<Vector3 const&>(pos), 0.1f, 5, 6, Color{255, 0, 255, 255});
                    //DrawLine3D(reinterpret_cast<Vector3 const&>(pos), reinterpret_cast<Vector3 const&>(posnrm), Color{0, 255, 0, 255});

                }
            }
        }
    }
}
