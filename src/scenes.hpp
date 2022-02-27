#pragma once

#include "game_state.hpp"

#include <functional>

namespace orni
{

using SceneFunc_t = std::function<void(orni::GameState&)>;

SceneFunc_t gen_test_scene_a();

SceneFunc_t gen_test_scene_b();

SceneFunc_t gen_test_scene_b_b();

SceneFunc_t gen_test_scene_c(GameState &rGame);

SceneFunc_t gen_the_game_scene(GameState &rGame);


}
