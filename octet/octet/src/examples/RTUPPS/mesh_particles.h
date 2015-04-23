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
    vec3 pos_predicted;
    vec3 vel;
    unsigned cell_id;
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
      size_t num_particles;
      size_t stabilizationIterations;
      size_t solverIterations;
      float particle_radius;
      std::chrono::time_point<std::chrono::system_clock> before; //used to obtain the time increment
      std::array<std::vector<uint8_t>, _NUM_PARTICLES_> grid_particles_id; //This is an array of vectors (one per cell grid) -- extremelly space inefficient!

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

      /// @brief This is the kernel function to stimate the density
      /// This can be done following this paper
      /// http://www.physics.umu.se/digitalAssets/60/60425_constraint-fluids-ieee.pdf
      /// We are implementing the Poly6 kernel as seen in page 12 of:
      /// http://image.diku.dk/kenny/download/vriphys10_course/sph.pdf
      float kernel_function_poly6(float r){
        if (r > _SMOOTHING_H_ || r < 0.0f) return 0.0f;
        float h_2 = _SMOOTHING_H_*_SMOOTHING_H_;
        float r_2 = r*r;
        float left = 315.0f / (64.f * 3.141592f*h_2*h_2);
        float right = (h_2 - r_2);
        return left*right*right*right;
      }

      /// @brief This is the kernel function to stimate the gradient of the density
      /// This can be done following this paper
      /// http://www.physics.umu.se/digitalAssets/60/60425_constraint-fluids-ieee.pdf
      /// For the gradient kernel for spiky we used http://www8.cs.umu.se/kurser/TDBD24/VT06/lectures/sphsurvivalkit.pdf
      float gradient_kernel_function_spiky(float r){
        if (r > _SMOOTHING_H_ || r < 0.0f) return 0.0f;
        float h_2 = _SMOOTHING_H_ * _SMOOTHING_H_;
        float left = -45.0f / (3.141592f * h_2 * h_2 * h_2);
        float right = _SMOOTHING_H_ - r;
        return left*right*right;
      }

      /// @brief This function will obtain the scaling factor for the id_particle particle
      /// To solve this we have to implement equation (11) from the paper
      float obtain_scaling_factor(unsigned i){
        float epsilon = 1.0f; //Relaxation parameter
        //Obtain density constraint
        float density_estimator = 0.0f;
        float gradient_constraint_sum = 0.0f;
        std::vector<uint8_t> * grid_cell = &grid_particles_id[particles_more[i].cell_id];
        for each (uint8_t j in (*grid_cell)){
          vec3 distance = particles_basic[i].pos - particles_basic[j].pos;
          float distancef = distance.length();
          density_estimator += kernel_function_poly6(distancef); //W(pi-pj,h)
          float gradient = gradient_kernel_function_spiky(distancef);
          gradient_constraint_sum += (gradient*gradient);
        }
        density_estimator *= particle_mass;
        float density_constraint = density_estimator/_REST_DENSITY_ - 1.0f;

        return (-1.0f * density_constraint) / (gradient_constraint_sum + epsilon);
      }

      /// @brief This is the simulation loop for only fluid simulation 
      /// http://mmacklin.com/pbf_sig_preprint.pdf
      void simulation_fluids(){
        //Calculate increment of time
        std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
        std::chrono::duration<float> elapsed_seconds = now - before;
        before = now;
        float time_inc = elapsed_seconds.count();

        //Here starts the Algorithm 1 from the Siggraph paper
        //For all particles i do
        for (unsigned i = 0; i != num_particles; ++i){
          // Apply forces v[i] = v[i] + time_inc*fext(particle[i])
          float f_ext = 0;
          particles_more[i].vel += time_inc*f_ext;
          // Predict position particle[i]' = particle[i] + time_inc*v[i]
          particles_more[i].pos_predicted = particles_basic[i].pos + time_inc*particles_more[i].vel;
        }

        //Find neighbouring particles
        find_neighbouring_particles();

        // while iter < solverIterations do
        for (unsigned iter = 0; iter != stabilizationIterations; ++iter){
          std::array<float, _NUM_PARTICLES_> lambda;
          // for all particles i do
          for (unsigned i = 0; i != num_particles; ++i){
            // Calculate lambda (aka Scaling Factor)
            lambda[i] = obtain_scaling_factor(i);
          }
          // for all particles i do
          for (unsigned i = 0; i != num_particles; ++i){
            // Calculate increment position
            // perform collision detection and reponse
          }
          // for all particles i do
          for (unsigned i = 0; i != num_particles; ++i){
            // update estimated position with the increment
          }
        }

        // for all particles i do
        for (unsigned i = 0; i != num_particles; ++i){
          // update velocit v[i] = 1/temp_inc * (particle[i]' - particle[i])
          // apply velocity confinement and XSPH viscosity
          // update positions particle[i] = particle[i]' or apply sleeping
        }
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