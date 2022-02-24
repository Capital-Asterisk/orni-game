/**
 * SPDX-License-Identifier: MIT
 * SPDX-FileCopyrightText: 2022 Neal Nicdao <chrisnicdao0@gmail.com>
 */
#pragma once

#include <glm/mat4x4.hpp>
#include <glm/gtc/quaternion.hpp>

#include <longeron/containers/intarray_multimap.hpp>
#include <longeron/id_management/registry.hpp>

#include <vector>
#include <variant>

namespace frogdyn
{

using frog_id_t = int;

struct AABB
{
    glm::vec3 m_min, m_max;
};

struct LinAng
{
    glm::vec3 m_lin;
    glm::vec3 m_ang;
};

struct FrogDyn
{
    struct Ball
    {
        glm::vec3   m_pos;
        float       m_radius;
        int         m_bitType;
        int         m_bitCollide;
    };

    struct BallOut
    {
        glm::vec3   m_pos;
        float       m_radius;
    };

    struct Bait
    {
        struct Insect
        {
            frog_id_t   m_id;
            glm::vec3   m_offset;

            // A->B
            glm::vec3 m_dir;
            glm::vec3 m_side;
        };

        Bait() = default;
        Bait(Insect a, Insect b) : m_a{a}, m_b{b} { }

        Insect m_a, m_b;

        // rotation needed to rotate B into A's space
        glm::quat m_origin{};

        float m_twistRange;
        float m_twistSpring;

        float m_coneRange;
        float m_coneSpring;

        bool m_doTwistLim{false};
        bool m_doTwistSpring{false};
        bool m_doConeLim{false};
        bool m_doConeSpring{false};
        bool m_doAlign{false};

    };

    struct CollisionCheck
    {
        frog_id_t m_a, m_b;
    };



    lgrn::IdRegistry<frog_id_t> m_ids;
    std::vector<glm::mat4x4>    m_tf;
    std::vector<LinAng>         m_extImp;
    std::vector<LinAng>         m_cstImp;
    std::vector<LinAng>         m_vel;
    std::vector<float>          m_mass;
    std::vector<float>          m_scale;

    std::vector<AABB> m_aabb;

    std::vector<frog_id_t> m_canCollide;
    lgrn::IntArrayMultiMap<frog_id_t, Ball>     m_balls;
    lgrn::IntArrayMultiMap<frog_id_t, BallOut>  m_ballPos;

    //std::vector<Ball> m_ballData;

    std::vector<Bait> m_baits;
};


using BallCollisions_t = lgrn::IntArrayMultiMap<frog_id_t, FrogDyn::BallOut>::Span;

struct BallContact
{
    glm::vec3   m_pos{0.0f};
    glm::vec3   m_nrm{0.0f};
    float       m_depth{0.0f};
    int         m_count{0};
};

constexpr bool aabb_intersect(AABB const& a, AABB const& b) noexcept
{
    return     (a.m_min.x < b.m_max.x) && (a.m_max.x > b.m_min.x)
            && (a.m_min.y < b.m_max.y) && (a.m_max.y > b.m_min.y)
            && (a.m_min.z < b.m_max.z) && (a.m_max.z > b.m_min.z);
}

struct BaitOptions
{
    float m_linP, m_linD;
    float m_twistP, m_twistD;
    float m_coneP, m_coneD;
    float m_alignP, m_alignD;
};

void apply_baits(FrogDyn &rDyn, BaitOptions opt, float delta);
void apply_ext_forces(FrogDyn &rDyn, float delta);
void apply_cst_forces(FrogDyn &rDyn, float delta);
void calc_balls_pos(FrogDyn &rDyn);
BallContact calc_ball_collisions(BallCollisions_t a, BallCollisions_t b);
void calc_frog_collisions(FrogDyn &rDyn);
frog_id_t add_frog(FrogDyn &rDyn, glm::mat4x4 const& tf, float mass);


}
