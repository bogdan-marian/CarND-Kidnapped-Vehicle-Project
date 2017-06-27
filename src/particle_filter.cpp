/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

  // Number of priticles
  num_particles = 100;

  // Random number generator
  std::default_random_engine gen;

  // Normal distributions for initial coordinates
  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[2]);

  // Initialize the weights and particles
  particles.resize(num_particles);
  weights.resize(num_particles);
  for (int i = 0; i < num_particles; i++)
  {
    double cur_x = dist_x(gen);
    double cur_y = dist_y(gen);
    double cur_theta = dist_theta(gen);
    double cur_weight = 1.;

    //id, x, y, theta, weight
    Particle cur_particle {i, cur_x, cur_y, cur_theta, cur_weight};
    particles[i] = cur_particle;
    weights[i] = cur_weight;
  }

  // Mark object as initialized
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

  // Random number generator
  std::default_random_engine gen;

  // Update position for each particle
  for (Particle& cur_particle : particles)
  {
    // Calculate the mean for the x, y and theta
    double new_x_mean = cur_particle.x + (std::abs(yaw_rate) < 0.0000001 ?
      velocity * delta_t * cos(cur_particle.theta) :
      velocity / yaw_rate * (sin(cur_particle.theta + yaw_rate * delta_t) - sin(cur_particle.theta)));
    double new_y_mean = cur_particle.y + (std::abs(yaw_rate) < 0.0000001 ?
      velocity * delta_t * sin(cur_particle.theta) :
      velocity / yaw_rate * (cos(cur_particle.theta) - cos(cur_particle.theta + yaw_rate * delta_t)));
    double new_theta_mean = cur_particle.theta + yaw_rate * delta_t;

    // Normal distributions for new coordinates
    std::normal_distribution<double> dist_x(new_x_mean, std_pos[0]);
    std::normal_distribution<double> dist_y(new_y_mean, std_pos[1]);
    std::normal_distribution<double> dist_theta(new_theta_mean, std_pos[2]);

    // Update particle with new estimate
    cur_particle.x = dist_x(gen);
    cur_particle.y = dist_y(gen);
    cur_particle.theta = dist_theta(gen);
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

  // Iterate through observations
  for (LandmarkObs& cur_obs : observations)
  {
    // Distance to nearest prediction
    double min_distance = INFINITY;
    // Nearest prediction id
    int min_id = 0;

    // For each predicted landmark
    for (const LandmarkObs& cur_pred : predicted)
    {
      // Distance from current observation to current prediction
      double cur_dist = dist(cur_pred.x, cur_pred.y, cur_obs.x, cur_obs.y);

      // If current distance is less then last known minimal distance,
      // update the minimum distance and nearest prediction id
      if (cur_dist < min_distance)
      {
        min_distance = cur_dist;
        min_id = cur_pred.id;
      }
    }

    // Assign nearest prediction id to currently known nearest id
    cur_obs.id = min_id;
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
