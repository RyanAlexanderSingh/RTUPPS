/// This file contains the class mesh_particle
///
/// Created by Sam Hayhurst, Juanmi Huertas and Ryan Singh
///
///

#ifndef _MESH_PARTICLES_INCLUDED_
#define _MESH_PARTICLES_INCLUDED_
#include <utility>
#include <chrono>
#include <ctime>
#include <math.h>

namespace octet{
  // grid size determines the half extents of the simulation space (ie 50 -> cube with 100 cells in each dimension)
  enum { _NUM_PARTICLES_ = 1000, _PARTICLE_DIAM = 1, _GRID_SIZE = 5, _SMOOTHING_H_ = 1, _REST_DENSITY_ = 1000 };

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
    unsigned cell_id;
    float lambda;
    vec3 pos_prev;
    vec3 vel;
    int index;
  };

  namespace scene{

    /// @brief This is the class mesh_particles. This class contains a mesh with a set of particles with RTUPPS
    /// This class contains it's own data structure to define the vetices of the mesh
    /// This class will ignore the "index" structure in openGL, and will work only with vertices
    /// For the detection of collisions we use a "cell grid". 
    class mesh_particles : public mesh{
      dynarray<particle_basic> particles_basic;
      dynarray<particle_more> particles_more;
      float particle_mass;
      float particle_invmass;
      float time_inc;    //just making this local as 5 functions use it(can pass through as parameter though) feel free to change
      size_t num_particles;
      size_t stabilizationIterations;
      size_t solverIterations;
      float particle_radius;
      vec2 gravity = (0.0f, -9.8f);
      std::chrono::time_point<std::chrono::system_clock> before; //used to obtain the time increment
      std::array<std::vector<uint8_t>, _NUM_PARTICLES_> grid_particles_id; //This is an array of vectors (one per cell grid) -- extremelly space inefficient!

      void apply_external_forces(){
      }

      void apply_viscosity(){
      }

      void advance_particles(){
      }

      void update_neighbours(){
      }

      void double_density_relaxation(){
      }

      void resolve_collisions(){
      }

      void update_velocity(){
      }

      /// @brief This function updates a vector of particle indices within the simulation loop,
      // http://docs.nvidia.com/cuda/samples/5_Simulations/particles/doc/particles.pdf
      // Attempting the "Building the Grid using Atomic operations"
      void find_neighbouring_particles(){
        // for each particle in the particle list determine its cell index
        unsigned int particle_id = 0;
        for each (particle_basic pb in particles_basic){
          int cell_index = 0;
          //obtain index of the cell for that particle
          cell_index += std::floor((pb.pos.x() + _GRID_SIZE) / _PARTICLE_DIAM);
          cell_index += _GRID_SIZE * 2 * std::floor((pb.pos.y() + _GRID_SIZE) / _PARTICLE_DIAM);
          cell_index += std::pow(_GRID_SIZE * 2, 2) * std::floor((pb.pos.z() + _GRID_SIZE) / _PARTICLE_DIAM);
          //put that particle into the cell
          grid_particles_id[cell_index].push_back(particle_id); //telling the cell which particles does it have
          particles_more[particle_id].cell_id = cell_index; //telling the particle which one is his cell id
          ++particle_id;
        }
      }

      /// @brief This is the simulation loop for only fluid simulation 
      /// http://mmacklin.com/pbf_sig_preprint.pdf
      void simulation_fluids(){
        //Calculate increment of time
        std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
        std::chrono::duration<float> elapsed_seconds = now - before;
        before = now;
        time_inc = elapsed_seconds.count();

        /// NEW CODE GOES HERE <=
        apply_external_forces(); //uses time_inc

        apply_viscosity(); //uses time_inc

        advance_particles(); //uses time_inc

        update_neighbours();

        double_density_relaxation(); //uses time_inc

        resolve_collisions();

        update_velocity(); //uses time_inc
      }

    public:
      mesh_particles() : num_particles(0), stabilizationIterations(0), solverIterations(0){}

      /// @brief This will initilize the mesh!
      void init(int type = 1, int n_stabilization = 10, int n_solver = 10){
        stabilizationIterations = n_stabilization;
        solverIterations = n_solver;
        particle_radius = _PARTICLE_DIAM * 0.5f;
        particle_mass = 1.0f;
        int num_particles = 1;

        if (type == 0){          // Initializate the particles with fixed positions
          for (int i = 0; i < num_particles; ++i){
            for (int j = 0; j < num_particles; ++j){
              for (int k = 0; k < num_particles; ++k){
                particle_basic new_particle;
                new_particle.pos = vec3(i - 5, k - 5, j - 5);
                new_particle.phase = 0;
                particles_basic.push_back(new_particle);
                particle_more more_particle;
                more_particle.vel = vec3(0, 0, 0);
                particles_more.push_back(more_particle);
              }
            }
          }
        }
        else { // used for testing grid particle location and postions
          particle_basic new_particle;
          new_particle.pos = vec3(0.5f, 0.5f, 0.5f); // should be 555 
          new_particle.phase = 0;
          particles_basic.push_back(new_particle);
          new_particle.pos = vec3(-4.5f, -4.5f, -4.5f); // should be 0
          particles_basic.push_back(new_particle);
          new_particle.pos = vec3(4.5f, 4.5f, 4.5f);    // should be 999 
          particles_basic.push_back(new_particle);
          new_particle.pos = vec3(-3.5f, -4.5f, -4.5f); // should be 1
          particles_basic.push_back(new_particle);
          new_particle.pos = vec3(-4.5f, -3.5f, -4.5f); // should be 10
          particles_basic.push_back(new_particle);
          new_particle.pos = vec3(-4.5f, -4.5f, -3.5f); // should be 100
          particles_basic.push_back(new_particle);
        }
        // Add particles to the mesh
        num_particles = particles_basic.size();
        num_particles = particles_basic.size();
        allocate(sizeof(particle_basic) * num_particles, sizeof(uint32_t)*num_particles);
        set_params(sizeof(particle_basic), num_particles, num_particles, GL_POINTS, GL_UNSIGNED_INT);
        //This will set up the attributes with position and color!
        clear_attributes();
        add_attribute(attribute_pos, 3, GL_FLOAT, 0);
        add_attribute(attribute_color, 4, GL_UNSIGNED_BYTE, 12, GL_TRUE);
        glPointSize(5.0f);

        //This block below is just copying into the attribute all the particles generated before
        gl_resource::wolock vlock(get_vertices());
        particle_basic* vtx = (particle_basic*)vlock.u8();
        gl_resource::wolock ilock(get_indices());
        uint32_t* idx = ilock.u32();
        for (unsigned i = 0; i != num_particles; ++i){
          vtx->pos = particles_basic[i].pos;
          if (particles_basic[i].phase == 0)
            vtx->phase = make_color(1.f, 0.0f, 0.0f);
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
        //We call the simulation fluids, that will be the algorithm for fluid simulation
        simulation_fluids();

        //This two lines take control on the points so we can modify them
        gl_resource::wolock vlock(get_vertices());
        particle_basic* vtx = (particle_basic*)vlock.u8();

        //This updates the position with the new position calculated in the simulation before
        for (unsigned i = 0; i != num_particles; ++i){
          vtx->pos = particles_basic[i].pos;
        }
      }
    };
  }
}

#endif