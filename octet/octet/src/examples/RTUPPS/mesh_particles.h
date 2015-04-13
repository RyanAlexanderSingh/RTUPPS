/// This file contains the class mesh_particle
///
/// Created by Sam Hayhurst, Juanmi Huertas and Ryan Singh
///
///

#ifndef _MESH_PARTICLES_INCLUDED_
#define _MESH_PARTICLES_INCLUDED_

namespace octet{
  enum {_NUM_PARTICLES_ = 1000};


  // this function converts three floats into a RGBA 8 bit color
  static uint32_t make_color(float r, float g, float b) {
    return 0xff000000 + ((int)(r*255.0f) << 0) + ((int)(g*255.0f) << 8) + ((int)(b*255.0f) << 16);
  }


  /// @brief This particle_basic contains all the info for each particle to be drawn
  struct particle_basic{
    vec3 pos;
    uint32_t phase;
  };

  /// @brief This particle_more contains more information needed to run the simulation
  struct particle_more{
    vec3 vel;
    float invmass;
  };

  namespace scene{

    /// @brief This is the class mesh_particles. This class contains a mesh with a set of particles with RTUPPS
    /// This class contains it's own data structure to define the vetices of the mesh
    /// This class will ignore the "index" structure in openGL, and will work only with vertices
    class mesh_particles : public mesh{

      dynarray<particle_basic> particles_basic;
      dynarray<particle_more> particles_more;
      size_t num_vertexes;
    public:
      mesh_particles(){}

      /// @brief This will initilize the mesh!
      void init(int type = 0){
        if (type == 0){          // Initializate the particles with fixed positions
          for (int i = 0; i < 10; ++i){
            for (int j = 0; j < 10; ++j){
              for (int k = 0; k < 10; ++k){
                particle_basic new_particle;
                new_particle.pos = vec3(i-5, k-5, j-5);
                new_particle.phase = 0;
                particles_basic.push_back(new_particle);
                particle_more more_particle;
                more_particle.vel = vec3(0, 0, 0);
                more_particle.invmass = 0.5f;
                particles_more.push_back(more_particle);
              }
            }
          }
        }
        // Add particles to the mesh
        num_vertexes = particles_basic.size();
        size_t num_indices = particles_basic.size();
        allocate(sizeof(particle_basic) * num_vertexes, sizeof(uint32_t)*num_indices);
        set_params(sizeof(particle_basic), num_indices, num_vertexes, GL_POINTS, GL_UNSIGNED_INT);
        //This will set up the attributes with position and color!
        clear_attributes();
        add_attribute(attribute_pos, 3, GL_FLOAT, 0);
        add_attribute(attribute_color, 4, GL_UNSIGNED_BYTE, 12, GL_TRUE);

        //This block below is just copying into the attribute all the particles generated before
        gl_resource::wolock vlock(get_vertices());
        particle_basic* vtx = (particle_basic*)vlock.u8();
        gl_resource::wolock ilock(get_indices());
        uint32_t* idx = ilock.u32();
        for (unsigned i = 0; i != num_vertexes; ++i){
          vtx->pos = particles_basic[i].pos;
          if (particles_basic[i].phase == 0)
            vtx->phase = make_color(1.f,0.0f,0.0f);
          ++vtx;
          *idx = i;
          ++idx;
        }
      }

      /// Serialise
      void visit(visitor &v){
        mesh::visit(v);
      }

      /// @brief This functions is where the "simulation" loop has to be written! 
      void update(){
        //This two lines take control on the points so we can modify them
        gl_resource::wolock vlock(get_vertices());
        particle_basic* vtx = (particle_basic*)vlock.u8();

        //We go through all the vertexes, changing their y position!
        for (unsigned i = 0; i != num_vertexes; ++i){
          vtx->pos.y() -= 0.01f;
          ++vtx;
        }
      }
    };
  }
}

#endif