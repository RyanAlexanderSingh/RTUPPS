////////////////////////////////////////////////////////////////////////////////
//
// (C) Juanmi, Ryan and Sam 2015
//

#include <vector>

namespace octet {
  /// Scene containing a box with octet.
  class RTUPPS : public app {

    struct Particle{
      vec3 position;
      vec3 velocity;
      float invmass;
      uint32_t phase;
    };

    // scene for drawing box
    ref<visual_scene> app_scene;

    std::vector<Particle> particles;

  public:
    /// this is called when we construct the class before everything is initialised.
    RTUPPS(int argc, char **argv) : app(argc, argv) {
    }

    /// this is called once OpenGL is initialized
    void app_init() {
      
      Particle p;
      p.position = vec3(0.0f);
      material *mat;
       
      app_scene = new visual_scene();
      app_scene->create_default_camera_and_lights();  
    }

    /// this is called to draw the world
    void draw_world(int x, int y, int w, int h) {
      int vx = 0, vy = 0;
      get_viewport_size(vx, vy);
      app_scene->begin_render(vx, vy);

      // update matrices. assume 30 fps.
      app_scene->update(1.0f / 30);

      // draw the scene
      app_scene->render((float)vx / vy);
    }
  };
}
