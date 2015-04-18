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
  enum {_NUM_PARTICLES_ = 1000, _PARTICLE_DIAM = 1, _GRID_SIZE = 5};

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
    float invmass;
  };

  namespace scene{

    /// @brief This is the class mesh_particles. This class contains a mesh with a set of particles with RTUPPS
    /// This class contains it's own data structure to define the vetices of the mesh
    /// This class will ignore the "index" structure in openGL, and will work only with vertices
    /// For the detection of collisions we use a "cell grid". 
    class mesh_particles : public mesh{
      dynarray<particle_basic> particles_basic;
      dynarray<particle_more> particles_more;
      size_t num_particles;
      size_t stabilizationIterations;
      size_t solverIterations; 
      float particle_radius;
      std::chrono::time_point<std::chrono::system_clock> before; //used to obtain the time increment
      std::vector<std::pair<uint8_t, uint8_t>> grid_particle_hash; // used to store the cell grid location of each particle
      std::vector<int> grid_particle_count; // used to store the number of particles per grid cell (should be capped)


      /// @brief This function is required for the neighbouring part, to obtain the hash of a particle
      std::pair<uint8_t, uint8_t> calculate_hash_particle(int particle_id){
        std::pair<uint8_t, uint8_t> pair;
        int cell_id = 0; // TO DO!!! need to update this to obtain the cell_id
        // plan here is to break the 'world' or simulation environment into cells of the same size as the 
        // diameter of each particle, this could be altered

        pair.first = cell_id;
        pair.second = particle_id;
        return pair;
      }

      /// @brief This function updates the particle locations within the discrete grid,
      /// cell indexing starts at 0 when x = y = z = -gridsize * cell_diameter
      /// indexing increases in x, y then z respectively
      void update_particle_grid(){
        // reset the particle count for each cell
        grid_particle_count.clear();
        // for each particle in the particle list determine its cell index
        for each (particle_basic pb in particles_basic){
          int cell_index = 0;
          cell_index += std::floor((pb.pos.x() + _GRID_SIZE) / _PARTICLE_DIAM);
          cell_index += _GRID_SIZE * 2 * std::floor((pb.pos.y() + _GRID_SIZE) / _PARTICLE_DIAM);
          cell_index += std::pow(_GRID_SIZE * 2, 2) * std::floor((pb.pos.z() + _GRID_SIZE) / _PARTICLE_DIAM);
          printf("cell index: %i\n", cell_index);
        }
      }

      /// @brief This function will calculate all neighbouring particles for each particle (using particle grid)
      void find_neighbouring_particles(){
        for (unsigned i = 0; i != num_particles; ++i){
          grid_particle_hash.push_back(calculate_hash_particle(i));
        }
        // TO DO!!!!
        //http://docs.nvidia.com/cuda/samples/5_Simulations/particles/doc/particles.pdf
        // Attempting the "Building the Grid using Sorting"
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
          // for all particles i do
            // Calculate lambda
          // for all particles i do
            // Calculate increment position
            // perform collision detection and reponse
          // for all particles i do
            // update new position
        }

        // for all particles i do
        for (unsigned i = 0; i != num_particles; ++i){
          // update velocit v[i] = 1/temp_inc * (particle[i]' - particle[i])
          // apply velocity confinement and XSPH viscosity
          // update positions particle[i] = particle[i]' or apply sleeping
        }
      }

      /// @brief This is a candidate for the simulation of all type of particles
      void simulation_all(){
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
          // Apply mass scaling mass[i]' = mass[i]*e^(-k*h(particle[i]'))
          float e = 2.71828f; //base of natural logarithms
          int k = 3; //coefficient k (set values from 1 to 5)
          float h = particles_more[i].pos_predicted.y();
          particles_more[i].invmass = particles_more[i].invmass * pow(e, (-k * h));
        }

        //For all particles i do
        for (unsigned i = 0; i != num_particles; ++i){
          // Find neighobring particles set 
          // Find solid contacts using particle set found above
        }

        // While iter < stabilizationIterations do
        for (unsigned iter = 0; iter != stabilizationIterations; ++iter){
          // increment praticle = 0, n = 0
          // solve contact constraints for increment particle, n
          // update particle[i] = particle[i] + increment particle/n
          // update particle' = particle' + increment particle/n
        }

        // while iter < solverIterations do
        for (unsigned iter = 0; iter != stabilizationIterations; ++iter){
          // for each constraint group G do
          // increment particle = 0, n = 0
          // solve all contraints in G for increment particle, n
          // update particle' = particle' + increment particle/n
        }

        // for all particles i do
        for (unsigned i = 0; i != num_particles; ++i){
          // update velocity v[i] = 1/temp_inc * (particle[i]' - particle[i])
          // advect diffuse particles, I think this is for fluid particles only
          // apply internal forces fdrag, fvort
          // update positions particle[i] = particle[i]' or apply sleeping
        }
      }
    public:
      mesh_particles() : num_particles(0), stabilizationIterations(0), solverIterations(0){}

      /// @brief This will initilize the mesh!
      void init(int type = 0, int n_stabilization = 10, int n_solver = 10){
        stabilizationIterations = n_stabilization;
        solverIterations = n_solver;
        particle_radius = _PARTICLE_DIAM * 0.5f;
        int num_particles = 1;
        grid_particle_count.reserve(sizeof(int) * std::pow(num_particles, 3));

        if (type == 0){          // Initializate the particles with fixed positions
          for (int i = 0; i < num_particles; ++i){
            for (int j = 0; j < num_particles; ++j){
              for (int k = 0; k < num_particles; ++k){
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
        else { // used for testing grid particle location and postions
          particle_basic new_particle;
          new_particle.pos = vec3(0.5f, 0.5f, 0.5f);
          new_particle.phase = 0;
          particles_basic.push_back(new_particle);
          new_particle.pos = vec3(-4.5f, -4.5f, -4.5f);
          particles_basic.push_back(new_particle);
          new_particle.pos = vec3(4.5f, 4.5f, 4.5f);
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
        // update particle grid positions
        update_particle_grid();

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