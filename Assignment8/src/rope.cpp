#include <iostream>
#include <vector>

#include "CGL/vector2D.h"

#include "mass.h"
#include "rope.h"
#include "spring.h"

namespace CGL {

    Rope::Rope(Vector2D start, Vector2D end, int num_nodes, float node_mass, float k, vector<int> pinned_nodes)
    {
        // TODO (Part 1): Create a rope starting at `start`, ending at `end`, and containing `num_nodes` nodes.
        float length = (end - start).norm();
        float interval = length / (num_nodes - 1);
        Vector2D direction = (end - start) / length;
 
        masses.push_back(new Mass(start, node_mass, false));
        for (int i = 1; i < num_nodes; i++) {
            Mass* prev_node = masses[i - 1];
            Mass* curr_node = new Mass(
                start + i * interval * direction, node_mass, false);
            
            masses.push_back(curr_node);
            springs.push_back(new Spring(prev_node, curr_node, k));
        }

       // Comment-in this part when you implement the constructor
        for (auto &i : pinned_nodes) {
           masses[i]->pinned = true;
        }
    }

    void Rope::simulateEuler(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {
            // TODO (Part 2): Use Hooke's law to calculate the force on a node
            float elastic_length = (s->m2->position - s->m1->position).norm();
            float magnitude = s->k * (elastic_length - s->rest_length);
            Vector2D direction = (s->m2->position - s->m1->position) / elastic_length;

            Vector2D force_m1to2 = magnitude * direction;

            s->m1->forces += force_m1to2;
            s->m2->forces -= force_m1to2;
        }

        float damping_factor = 0.005;
        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                // TODO (Part 2): Add the force due to gravity, then compute the new velocity and position
                m->forces += gravity;
                m->forces -= damping_factor * m->velocity;
                m->velocity += m->forces / m->mass * delta_t;
                m->position += m->velocity * delta_t;
                
                // TODO (Part 2): Add global damping
            }

            // Reset all forces on each mass
            m->forces = Vector2D(0, 0);
        }
    }

    void Rope::simulateVerlet(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {
            // TODO (Part 3): Simulate one timestep of the rope using explicit Verlet （solving constraints)
            float elastic_length = (s->m2->position - s->m1->position).norm();
            float magnitude = s->k * (elastic_length - s->rest_length);
            Vector2D direction = (s->m2->position - s->m1->position) / elastic_length;

            Vector2D force_m1to2 = magnitude * direction;

            s->m1->forces += force_m1to2;
            s->m2->forces -= force_m1to2;
        }

        float damping_factor = 0.0001;
        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                // TODO (Part 3.1): Set the new position of the rope mass
                m->forces += gravity;
                Vector2D acceleration = m->forces / m->mass;

                m->position = m->last_position 
                    + (1.0 - damping_factor) * (m->last_position - m->start_position)
                        + acceleration * std::pow(delta_t, 2);
                
                m->start_position = m->last_position;
                m->last_position = m->position;
                
                // TODO (Part 4): Add global Verlet damping
            }

            m->forces = Vector2D(0.0, 0.0);
        }
    }
}
