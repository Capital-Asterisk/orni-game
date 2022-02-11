#pragma once

#include "game_state.hpp"

#include <functional>

namespace orni
{

using SceneFunc_t = std::function<void(orni::GameState&)>;

SceneFunc_t gen_test_scene_a();

SceneFunc_t gen_test_scene_b();

SceneFunc_t gen_test_scene_b_b();

}
