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
  enum { _NUM_PARTICLES_ = 100, _PARTICLE_DIAM = 1, _GRID_SIZE = 20, _SMOOTHING_H_ = 1, _REST_DENSITY_ = 1000 };

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
    vec3 pos_prev;
    vec3 vel;
    int index;
    std::vector<int> neighbours; // neighbouring particles
  };

  struct DistanceField_Cell{
    float distance;
    vec3 normal;
    vec3 tangent;
  };

  /// @brief This helper fucntion returns the maximum value in vector 
  /// This is menat to be work around specifically for an indexing issue in
  /// the double density function; ie this should ideally be designed out.
  unsigned get_max_value(std::vector<int> vec){
    unsigned val = 0;
    for (size_t i = 0; i != vec.size(); i++)
    {
      if (vec[i] > val){
        val = vec[i];
      }
    }
    return val;
  }

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

      float k; // stiffness
      float k_near; // stiffness
      float p_0; // rest_density
      float sigma; // viscosity's linear dependence
      float beta; // viscosity's quadratic dependence
      float collision_radius;
      float friction;
      float collision_softness;
      float total_time;

      size_t num_particles;
      int grid_size;
      float particle_radius;
      std::chrono::time_point<std::chrono::system_clock> before; //used to obtain the time increment
      std::map<int, std::vector<int>> grid_particles_id; //This is an array of vectors (one per cell grid) -- extremelly space inefficient!
      std::vector<DistanceField_Cell> distance_field;

      /// @brief This function will apply the external forces to the particle
      /// Right now, the only external force is the gravity, although it could be interesting to add user interaction
      void apply_external_forces(float time_inc){
        for (unsigned p = 0; p != num_particles; ++p){
          // Apply gravity
          float inc_vel = -9.8f * time_inc;
          //if (p == 0) printf("\nStart Velocity: %f\n", particles_more[p].vel.y());
          particles_more[p].vel.y() += inc_vel;
          //if (p == 0) printf("Mid Velocity: %f", particles_more[p].vel.y());
          //potentially for the future use mouse x,y inputs
          //p.vec += forces_from_touch_input(p);
        }
      }

      /// @brief This function applies the viscosity to the particles.
      /// To work, it needs correct values for sigma and beta (linear and quadratic dependence of the particles)
      void apply_viscosity(float time_inc){
        for (unsigned p = 0; p != num_particles; ++p){
          unsigned cell_id = particles_more[p].cell_id;
          unsigned size_neighbours = particles_more[p].neighbours.size();
          for (unsigned i = 0; i != size_neighbours; ++i){
            unsigned n = particles_more[p].neighbours[i];
            if (i != n){
              vec3 distance = particles_basic[n].pos - particles_basic[p].pos;
              float vel_inward = (particles_more[p].vel - particles_more[n].vel).dot(distance);
              if (vel_inward > 0.0f){
                float length = distance.length();
                vel_inward /= length;
                float q = length / particle_radius;
                vec3 impulse = 0.5f*time_inc*(1.0f - q)*(sigma*vel_inward + beta*vel_inward*vel_inward)*distance;
                particles_more[p].vel -= impulse;
              }
            }
          }
        }
      }

      /// @brief This function will change the position of the particles
      /// This modification will be done using the velocity of the particle,
      /// it will not be visible until we update the info of the shader
      void advance_particles(float time_inc){
        for (unsigned p = 0; p != num_particles; ++p){
          particles_more[p].pos_prev = particles_basic[p].pos;
          //if (p == 0) printf("\nUsing velocity: %f", particles_more[p].vel.y());
          vec3 inc_pos = vec3(0.0f, particles_more[p].vel.y()*time_inc, 0.0f);
          particles_basic[p].pos += inc_pos;
        }
      }

      void print_neighbours(){
        printf("\nParticle 0 neighbours\n");
        unsigned neighoburs_size = particles_more[0].neighbours.size();
        for (unsigned p = 0; p != neighoburs_size; ++p){
          printf("p: %i ", particles_more[0].neighbours[p]);
        }
        printf("\nParticle 6 neighbours\n");
        neighoburs_size = particles_more[6].neighbours.size();
        for (unsigned p = 0; p != neighoburs_size; ++p){
          printf("p: %i ", particles_more[6].neighbours[p]);
        }
        printf("\n");
      }

      /// @brief This function updates the list of neighbouring particles that each
      /// particle holds, the function find neighbouring particles calculates a list of
      /// of possible neighbours by sorting the particles into a grid
      void update_neighbours(){
        // update the grid, sort the paricles and their ids into the grid
        update_particle_grid_positions();
        for each (particle_more &pm in particles_more){
          // clear the list of neighbouring particles
          pm.neighbours.clear();
          std::vector<int> possible_neighbours = get_possible_neighbours(pm.cell_id, pm.index);
          for each (int particle_id in possible_neighbours){
            // calculate the displacement given possible nighbour is from the particle in question
            float disp = (particles_basic[particle_id].pos - particles_basic[pm.index].pos).length();
            if (disp < _PARTICLE_DIAM * 0.5f){
              pm.neighbours.push_back(particle_id);
            }
          }
        }
      }

      /// @brief This function returns a vector of particle ids that are possible neighbours to a 
      /// particle in a given cell id.
      std::vector<int> get_possible_neighbours(unsigned int cell_id, unsigned int particle_id){
        std::vector<int> possible_neigbours;
        // detect edge cases and resize the for loops accordingly
        int left_index = (cell_id % _GRID_SIZE == 0) ? 0 : -1;            // is the cell in the left most column?
        int right_index = (cell_id % (_GRID_SIZE - 1) == 0) ? 1 : 2;       // is the cell in the right most column?
        int up_index = (cell_id >= _GRID_SIZE * (_GRID_SIZE - 1)) ? 1 : 2; // is the cell in the top row?
        int down_index = (cell_id <= _GRID_SIZE) ? 0 : -1;                // is the cell in the first row?

        for (int j = down_index; j != up_index; ++j){      // loop columns
          for (int i = left_index; i != right_index; i++){ // loop rows
            unsigned int id = cell_id + i + _GRID_SIZE * j;
            // detect the case where a particle is added to its own list of neighbours
            if (i == 0 && j == 0){ // run check only for particles own cell 
              for (size_t i = 0; i < grid_particles_id[id].size(); ++i)
              {
                if (grid_particles_id[id][i] != particle_id){
                  possible_neigbours.push_back(grid_particles_id[id][i]);
                }
              }
            }
            else{
              possible_neigbours.insert(possible_neigbours.end(), grid_particles_id[id].begin(), grid_particles_id[id].end());
            }
          }
        }
        return possible_neigbours;
      }

      /// @brief This function is in charge of the density relaxation
      /// This function is done following the paper (using a formulation of the SPH paradigm).
      /// We are avoiding to use a particle with itself, in case that it's considered its own 
      ///   neighbour
      void double_density_relaxation(float time_inc){
        for (unsigned i = 0; i != num_particles; ++i){
          float density = 0;
          float density_near = 0;
          std::array<float, _NUM_PARTICLES_> distances;
          unsigned cell_id = particles_more[i].cell_id;
          unsigned size_neighbours = particles_more[i].neighbours.size();
          for (unsigned j = 0; j != size_neighbours; ++j){
            unsigned n = particles_more[i].neighbours[j];
            if (n != i){ //to avoid being moved by itself!
              distances[n] = vec3(particles_basic[i].pos - particles_basic[n].pos).length();
              float q = 1.0f - (distances[n] / particle_radius);
              density += q*q;
              density_near += q*q*q;
            }
          }
          density = k * (density - p_0);
          density_near = k_near * density_near;
          vec3 delta = vec3(0);
          for (unsigned j = 0; j != size_neighbours; ++j){
            unsigned n = particles_more[i].neighbours[j];
            if (i != n){ //to avoid being moved by itself!
              float q = 1.0f - (distances[n] / particle_radius);
              vec3 direction = vec3(particles_basic[i].pos - particles_basic[n].pos) / distances[n];
              vec3 D = 0.5*time_inc*time_inc*(density*q + density_near*q*q)*direction;
              particles_basic[n].pos += D;
              delta -= D;
            }
          }
          particles_basic[i].pos += delta;
        }
      }

      /// @brief This function is in charge of resolving the collision between the particles and the world (bounding box!)
      /// This function requires the construction and updating of the distance field
      void resolve_collisions(float time_inc){
        for (unsigned p = 0; p != num_particles; ++p){
          int index = 2; //distanceField.GetIndex(p.pos)
          if (index != -1){
            float distance = particles_basic[p].pos.y() + 10.0f; //distanceField.GetDistance(index)
            if (distance < 0.0f){ //distance > -collisionRadius
              vec3 direction = (particles_basic[p].pos - particles_more[p].pos_prev).normalize();
              vec3 normal = vec3(0.0f, 1.0f, 0.0f); //distanceField.getNormal(index);
              vec3 tangent = vec3(1.0f, 0.0f, 0.0f); //PerpendicularCCW(normal);
              tangent = time_inc*friction*(direction.dot(tangent))*tangent;
              particles_basic[p].pos -= tangent;
              particles_basic[p].pos -= collision_softness*(distance)*normal;
            }
          }
        }
      }

      /// @brief This function updates the velocity of a particle by its previous and current position
      void update_velocity(float time_inc){
        for (unsigned p = 0; p != num_particles; ++p){
          vec3 new_velocity = (particles_basic[p].pos - particles_more[p].pos_prev) / time_inc;
          particles_more[p].vel = vec3(new_velocity);
        }
      }

      /// @brief This function updates the grid (an array of vectors) with a vector for each cell
      /// that contains the particle ids of the particles currently within that cell
      void update_particle_grid_positions(){
        // first clear the grid particle hash
        grid_particles_id.clear();
        // for each particle in the particle list determine its cell index
        unsigned int particle_id = 0;
        for each (particle_basic pb in particles_basic){
          int cell_index = 0;
          //obtain index of the cell for that particle
          cell_index += std::floor((pb.pos.x() + _GRID_SIZE) / _PARTICLE_DIAM);
          cell_index += _GRID_SIZE * 2 * std::floor((pb.pos.y() + _GRID_SIZE) / _PARTICLE_DIAM);
          //put that particle into the cell
          grid_particles_id[cell_index].push_back(particle_id); //telling the cell which particles does it have
          if (cell_index < 0 || cell_index > _GRID_SIZE *_GRID_SIZE * 4){
            // printf("WRONG!\n");
          }
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

        total_time += 0.01f;
        float time_inc = total_time; //Use this line while debugging
        //float time_inc = elapsed_seconds.count()*0.5f; //Use this line for the release
        printf(" \n\ntime_inc.. %f", time_inc);
        apply_external_forces(time_inc);

        //    printf("Applying viscosity.. ");
        //     apply_viscosity(time_inc);

        apply_external_forces(time_inc);

        //printf("Applying viscosity.. ");
        //apply_viscosity(time_inc);

        //printf("Advancing particles.. ");
        advance_particles(time_inc);

        //printf("Updating neighbours.. ");
        update_neighbours();
        print_neighbours();
        //printf("Relaxing double density.. ");
        //double_density_relaxation(time_inc);

        //printf("Colliding resolutions.. ");
        resolve_collisions(time_inc);

        //printf("Updating velocity.. ");
        update_velocity(time_inc);
      }

    public:
      mesh_particles() : num_particles(0){}

      /// @brief This will initilize the mesh!
      void init(int type = 0){
        particle_radius = _PARTICLE_DIAM * 0.5f;
        particle_mass = 1.0f;
        total_time = 0.0f;
        collision_radius = particle_radius + 0.1f;
        friction = 0.5f;
        collision_softness = 0.6f;
        k = 1.0f; // stiffness
        k_near = 1.0f; // stiffness
        p_0 = 1.0f; // rest_density
        num_particles = _NUM_PARTICLES_;
        grid_size = _GRID_SIZE;
        int sqrt_num_particles = sqrt(num_particles);
        int index = 0;

        if (type == 0){          // Initializate the particles with fixed positions
          for (int i = 0; i < sqrt_num_particles; ++i){
            for (int j = 0; j < sqrt_num_particles; ++j){
              particle_basic new_particle;
              new_particle.pos = vec3(i * 0.5f - 5, j * 0.5f - 5, -5);
              new_particle.phase = 0;
              particles_basic.push_back(new_particle);
              particle_more more_particle;
              more_particle.vel = vec3(0.0f, 10.0f, 0.0f);
              more_particle.index = index++;
              more_particle.cell_id = -1;
              more_particle.neighbours.reserve(36); //maximum amount of particle neighbours: 9 cells with 4 particles in each
              particles_more.push_back(more_particle);
            }
          }
        }

        // add the particles to the grid
        update_neighbours();

        // Add particles to the mesh
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

        before = std::chrono::system_clock::now();
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
          vtx++;
        }
      }
    };
  }
}

#endif