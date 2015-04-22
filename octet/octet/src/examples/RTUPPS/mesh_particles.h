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
  enum { _NUM_PARTICLES_ = 1000, _PARTICLE_DIAM = 1, _GRID_SIZE = 5 };

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
    float invmass;
    float mass;
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

      /// @brief This function updates a vector of particle indices within the simulation loop,
      // This loop is what is used to run collision calcs on.
      // This fucntion returns true if the input particle has nearby particles
      // http://docs.nvidia.com/cuda/samples/5_Simulations/particles/doc/particles.pdf
      // Attempting the "Building the Grid using Sorting"
      void find_neighbouring_particles(std::vector<int> &neighbouring_particles, int particle_id){
        // find cell index of the given particle
        int particle_cell_index = grid_particle_hash[particle_id].second;
        //printf("cell_index: %i\n", particle_cell_index);
        //printf("part_id: %i\n", particle_id);

        // check whether this particular particle has any neighbouring particles
        for (int k = -1; k != 2; ++k){  // loops over cells in z _GRID_SIZE ^ 2 for each increment
          for (int j = -1; j != 2; ++j){  // loops over cells in y _GRID_SIZE for each increment
            for (int i = -1; i != 2; ++i){ // loops over cells in x 1 for each increment
              int cell_id = particle_cell_index + k * std::pow((int)_GRID_SIZE, 2) + j * _GRID_SIZE + i;
              //printf("cell_id: %i\n", cell_id);
              // if there are particles in the selcted cell add their ids to the neighbours vector
              if (cell_id < 0)
                break;
              int num_neighbours = grid_particle_count[cell_id];
              if (num_neighbours != 0){
                for (int index = 0; index != num_neighbours; ++index){ // find the number of particles
                  for (int part = 0; part != grid_particle_hash.size(); ++part){  // search through the vector of pairs searching whether their cells match, this could be done more efficiently
                    if (grid_particle_hash[part].second == cell_id && grid_particle_hash[part].first != particle_id){ // avoid adding input particle to neighbours
                      neighbouring_particles.push_back(grid_particle_hash[part].first);
                    }
                  }
                }
              }
            }
          }
        }
      }

      /// @brief This fucntion loops through each particle in the simulation and assigns it a cell
      /// within a 3D grid defined above. This function is to be run once per simulation step and updates
      /// the grid particle hash and count vectors which are then used by the find_neighbouring particles function
      /// cell indexing starts at 0 when x = y = z = -gridsize * cell_diameter
      /// indexing increases in x, y then z respectively
      void sort_particles_into_grid(){
        grid_particle_count.assign(std::pow(_GRID_SIZE * 2, 3), 0); // reset cell particle counts to zero
        grid_particle_hash.clear(); // this should be changed to be an access and change rather than a clear and create
        // for each particle in the particle list determine its cell index
        unsigned int particle_id = 0;
        for each (particle_basic pb in particles_basic){
          int cell_index = 0;
          cell_index += std::floor((pb.pos.x() + _GRID_SIZE) / _PARTICLE_DIAM);
          cell_index += _GRID_SIZE * 2 * std::floor((pb.pos.y() + _GRID_SIZE) / _PARTICLE_DIAM);
          cell_index += std::pow(_GRID_SIZE * 2, 2) * std::floor((pb.pos.z() + _GRID_SIZE) / _PARTICLE_DIAM);
          //printf("cell index: %i\n", cell_index);
          grid_particle_count[cell_index]++; // increment the number of particles in that particular cell;
          grid_particle_hash.push_back(std::pair<uint8_t, uint8_t>(particle_id, cell_index));
          ++particle_id;
        }
      }

      /// @brief This is the kernel function to stimate the density
      /// This can be done following this paper
      /// http://www.physics.umu.se/digitalAssets/60/60425_constraint-fluids-ieee.pdf
      /// and this wikipedia link
      /// http://en.wikipedia.org/wiki/Kernel_smoother
      float kernel_function(unsigned i, unsigned j){
        return 0.5f;
      }

      /// @brief This function will obtain the scaling factor for the id_particle particle
      /// To solve this we have to implement equation (11) from the paper
      float obtain_scaling_factor(unsigned i){
        float epsilon = 1.0f; //Relaxation parameter
        //Obtain density constraint
        float density_estimator = 0.0f;
        for (unsigned j = 0; j != num_particles; ++j){ // <= THIS LINE IS WRONG. Has to be done with the neighbour particles only! FIX!
          density_estimator += particles_more[j].mass*kernel_function(i,j);
        }
        float density_constraint = 0.0f;
        float gradient_constraint_sum = 0.0f;

        return 0.0f; // temporary return value to fix build. Sam
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

        //Sort particles using 3D grid
        sort_particles_into_grid();

        //Find neighbouring particles
        std::vector<int> vec; // stores the particle ids of all neighbouring particles to given particle id.
        for (int i = 0; i < particles_basic.size(); ++i){
          find_neighbouring_particles(vec, i); // something needs to be done with this vector, not sure what
        }

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
        int num_particles = 1;
        grid_particle_count.assign(std::pow(_GRID_SIZE * 2, 3), 0);  // assign all values inside vector 0;

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
                more_particle.invmass = 0.5f;
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

        sort_particles_into_grid();
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