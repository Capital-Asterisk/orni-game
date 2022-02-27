#pragma once

struct Font;

namespace orni
{

using SceneId_t = int;

struct GameState
{
    // put persistent game data here
    Font *m_pFont;

    // Scene change request, literally just an integer lol
    SceneId_t m_sceneChange{0};

};



}
