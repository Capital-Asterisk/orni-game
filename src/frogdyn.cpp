/**
 * SPDX-License-Identifier: MIT
 * SPDX-FileCopyrightText: 2022 Neal Nicdao <chrisnicdao0@gmail.com>
 */
#include "frogdyn.hpp"

#include <glm/gtx/norm.hpp>
#include <glm/gtx/matrix_interpolation.hpp>
#include <glm/gtx/vector_angle.hpp>
#include <glm/gtx/projection.hpp>

using namespace frogdyn;

void frogdyn::apply_baits(FrogDyn &rDyn, BaitOptions opt, float delta)
{
    static glm::mat4x4 const mcidentity{1.0f};
    static LinAng const linAngZero{};

    // calculate baits
    for (int i = 0; i < rDyn.m_baits.size(); i ++)
    {
        auto const &rBait = rDyn.m_baits[i];

        bool const worldA = rBait.m_a.m_id == -1;

        // get stuff for A
        glm::mat4x4 const&  tfA         = worldA ? mcidentity : rDyn.m_tf[rBait.m_a.m_id];
        float const         scaleA      = worldA ? 1.0f : rDyn.m_scale[rBait.m_a.m_id];
        glm::vec3 const     offsetRotA  = tfA * glm::vec4(rBait.m_a.m_offset * scaleA, 0.0f);
        glm::vec3 const     posA        = offsetRotA + glm::vec3(tfA[3]);
        LinAng const&       velA        = worldA ? linAngZero : rDyn.m_vel[rBait.m_a.m_id];
        glm::vec3 const     pointVelA   = worldA ? glm::vec3{0.0f} : (velA.m_lin + glm::cross(velA.m_ang, offsetRotA));
        float const         massA       = worldA ? 9999999.0f : rDyn.m_mass[rBait.m_a.m_id];

        // get stuff for B
        glm::mat4x4 const&  tfB         = rDyn.m_tf[rBait.m_b.m_id];
        float const         scaleB      = rDyn.m_scale[rBait.m_b.m_id];
        glm::vec3 const     offsetRotB  = (tfB * glm::vec4(rBait.m_b.m_offset * scaleB, 0.0f));
        glm::vec3 const     posB        = offsetRotB + glm::vec3(tfB[3]);
        LinAng const&       velB        = rDyn.m_vel[rBait.m_b.m_id];
        glm::vec3 const     pointVelB   = velB.m_lin + glm::cross(velB.m_ang, offsetRotB);
        float const         massB       = rDyn.m_mass[rBait.m_b.m_id];

        glm::vec3 const     posRel      = posB - posA;
        glm::vec3 const     velRel      = pointVelB - pointVelA;

        //float const lenBSq = glm::length2(rBait.m_b.m_offset);
        //float const lenASq = glm::length2(rBait.m_a.m_offset);
        //float const maxLen  = glm::max(lenASq, lenBSq);

        float const minMass = glm::min(massA, massB);

        glm::vec3 angImpA{0.0f}, angImpB{0.0f};

        glm::vec3 const     wdirA       = tfA * glm::vec4{rBait.m_a.m_dir, 0.0f};
        glm::vec3 const     wdirB       = tfB * glm::vec4{rBait.m_b.m_dir, 0.0f};

        glm::vec3 const     wsideA      = tfA * glm::vec4{rBait.m_a.m_side, 0.0f};

        if (rBait.m_doTwistLim || rBait.m_doTwistSpring)
        {
            // need to do twist calculations


            glm::vec3 const wsideB      = tfB * glm::vec4{rBait.m_b.m_side, 0.0f};

            glm::quat const dragon      = glm::rotation(wdirB, wdirA);
            float const     roll        = glm::orientedAngle(wsideA, dragon * wsideB, wdirA);

            float const     lim         = rBait.m_twistRange;

            if (rBait.m_doTwistLim && (-lim > roll || roll > lim))
            {
                float const rollSpd     = glm::dot(velA.m_ang, wdirA) - glm::dot(velB.m_ang, wdirB);

                float const phog = (roll > 0.0f) ? (roll - lim) : (roll + lim);
                angImpA += (phog * opt.m_twistP - rollSpd * opt.m_twistD) * wdirA * minMass;
                angImpB += (-phog * opt.m_twistP + rollSpd * opt.m_twistD) * wdirB * minMass;
            }

            if (rBait.m_doTwistSpring)
            {
                angImpA += (roll * rBait.m_twistSpring) * wdirA * delta;
                angImpB += (-roll * rBait.m_twistSpring) * wdirB * delta;
            }
        }

        if (rBait.m_doConeLim || rBait.m_doConeSpring)
        {
            // do cone calculations
            float const     angDiff     = glm::angle(wdirA, wdirB);
            glm::vec3 const dirCross    = glm::cross(wdirA, wdirB);
            glm::vec3 const axis        = glm::normalize(glm::cross(wdirA, wdirB));

            //std::cout << std::setprecision(5);
            //std::cout << "verify: " << (glm::asin(glm::length(dirCross)) - angDiff) << "(" << (glm::asin(glm::length(dirCross))) << " - " << angDiff << ")\n";

            if (angDiff > 0.004)
            {

            float const     lim         = rBait.m_coneRange;

            if (rBait.m_doConeLim && angDiff > lim)
            {
                float const swaySpd    = glm::dot(velA.m_ang, axis) - glm::dot(velB.m_ang, axis);

                float const phog = angDiff - lim;
                angImpA += (phog * opt.m_coneP - swaySpd * opt.m_coneD) * axis * minMass;
                angImpB += (-phog * opt.m_coneP + swaySpd * opt.m_coneD) * axis * minMass;
            }

            if (rBait.m_doConeSpring)
            {
                angImpA += (angDiff * rBait.m_coneSpring) * axis * delta;
                angImpB += (-angDiff * rBait.m_coneSpring) * axis * delta;
            }
            }
        }

        if (rBait.m_doAlign)
        {
            glm::vec3 const dirBFlat    = glm::normalize(wdirB - glm::proj(wdirB, wsideA));

            glm::quat const     gender  = glm::rotation(dirBFlat, wdirB);
            glm::vec3 const     axis    = glm::axis(gender);
            float const         angle   = glm::angle(gender);
            float const         swaySpd = glm::dot(velA.m_ang, axis) - glm::dot(velB.m_ang, axis);

            angImpA += (angle * opt.m_alignP - swaySpd * opt.m_alignD) * axis * minMass;
            angImpB += (-angle * opt.m_alignP + swaySpd * opt.m_alignD) * axis * minMass;
        }


        glm::vec3 const     normalLinImp = (posRel * opt.m_linP + velRel * opt.m_linD) * minMass;
        glm::vec3 const     normalAngImp = normalLinImp;

        if (! worldA)
        {
            rDyn.m_cstImp[rBait.m_a.m_id].m_lin += normalLinImp;
            rDyn.m_cstImp[rBait.m_a.m_id].m_ang += glm::cross(offsetRotA, normalAngImp) + angImpA;
        }

        rDyn.m_cstImp[rBait.m_b.m_id].m_lin += -normalLinImp;
        rDyn.m_cstImp[rBait.m_b.m_id].m_ang += -glm::cross(offsetRotB, normalAngImp) + angImpB;
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

frog_id_t frogdyn::add_frog(FrogDyn &rDyn, glm::mat4x4 const& tf, float mass)
{
    frog_id_t const id = rDyn.m_ids.create();
    std::size_t const capacity = rDyn.m_ids.capacity();
    rDyn.m_cstImp   .resize(capacity);
    rDyn.m_extImp   .resize(capacity);
    rDyn.m_vel      .resize(capacity);
    rDyn.m_mass     .resize(capacity);
    rDyn.m_tf       .resize(capacity);
    rDyn.m_scale    .resize(capacity);
    rDyn.m_aabb     .resize(capacity);
    rDyn.m_balls    .resize_ids(capacity);

    rDyn.m_scale[id]    = 1.0f;
    rDyn.m_mass[id]     = mass;
    rDyn.m_tf[id]       = tf;

    return id;
}
