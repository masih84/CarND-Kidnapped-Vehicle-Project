/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   *  Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   *  Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 20;  //  Set the number of particles
  std::default_random_engine gen;

  //  creates a normal (Gaussian) distribution for x, y and theta
  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[2]);

  Particle sample;
  for (int i = 0; i < num_particles; ++i) {

	  // set initial particle based on normal distribution
	  sample.x = dist_x(gen);
	  sample.y = dist_y(gen);
	  sample.theta = dist_theta(gen);
	  sample.weight = 1.;
	  sample.id = i;
	  particles.push_back(sample);

	  // set initial weights to 1;
	  weights.push_back(1.);
	  
	  // std::cout << "samplex is " << sample.x << std::endl;
	  // std::cout << "sampley is " << sample.y << std::endl;
	  // std::cout << "sample is is " << sample.id << std::endl;
  }
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
   
	std::default_random_engine gen;

	// This line creates a normal (Gaussian) distribution for x
	double sample_x, sample_y, sample_theta;
	double updated_x, updated_y, updated_theta;

	// min yaw_rate limits  
	double yaw_rate_limit = 0.1;

	for (int i = 0; i < num_particles; ++i) {
		//  creates a normal (Gaussian) distribution for x, y and theta
		std::normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
		std::normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
		std::normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);

		sample_x = dist_x(gen);
		sample_y = dist_y(gen);
		sample_theta = dist_theta(gen);

		// if yaw_rate greater than its limits, update with yaw_rate  
		if (fabs(yaw_rate) > yaw_rate_limit) {

			updated_x = sample_x + velocity / yaw_rate * (sin(sample_theta + yaw_rate * delta_t) - sin(sample_theta));
			updated_y = sample_y + velocity / yaw_rate * (-1.*cos(sample_theta + yaw_rate * delta_t) + cos(sample_theta));
			updated_theta = sample_theta + yaw_rate * delta_t;
		}
		else {
			updated_x = sample_x + velocity  * delta_t * cos(sample_theta);
			updated_y = sample_y + velocity  * delta_t * sin(sample_theta);
			updated_theta = sample_theta + yaw_rate * delta_t;
		}

		particles[i].x = updated_x;
		particles[i].y = updated_y;
		particles[i].theta = updated_theta;
		/*
		std::cout << "updated_x is " << particles[i].x << std::endl;
		std::cout << "updated_y is " << particles[i].y << std::endl;
		std::cout << "updated_theta is " << particles[i].theta << std::endl;
		std::cout << "particle is is " << particles[i].id << std::endl;
		*/
	}
	std::string ready;
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted,
	vector<LandmarkObs>& observations) {
	/**
	 * Find the predicted measurement that is closest to each
	 *   observed measurement and assign the observed measurement to this
	 *   particular landmark.
	 * NOTE: this method will NOT be called by the grading code. But you will
	 *   probably find it useful to implement this method and use it as a helper
	 *   during the updateWeights phase.
	 */
	double d, min_d =10000;
	for (int i = 0; i < observations.size(); i++) {
		// std::cout << "initial observation id is " << observations[i].id <<std::endl;

		for (int j = 0; j < predicted.size(); j++) {
			// calculate distance between predicted and observations
			d = dist(predicted[j].x, predicted[j].y, observations[i].x, observations[i].y);
			if ((j == 0) || ((j > 0) && (d < min_d))) {
				// update obs index id the distance is smaller than min_d
				observations[i].id = predicted[j].id; //land mark id
				min_d = d;
			}
		}
		// std::cout << "final observation id is " << observations[i].id << std::endl;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   
   */
	// Step 1 define prediction vector
	vector<LandmarkObs> predicts; // THIS is in MAP coordinate systm (x,y);
	LandmarkObs predicted;
	std::vector<double> wi;
	vector<LandmarkObs> observations_map;
	wi.clear();

	// step 2 loop in all land_mark and find the ones in sensor range

	double sum_particle_weights = 0; // sum of all weights
	for (int j = 0; j < num_particles; ++j) {


		double p_X, p_Y, p_theta, d;

		p_X = particles[j].x;
		p_Y = particles[j].y;
		p_theta = particles[j].theta;
		predicts.clear();
		/*
		std::cout << " sensor range is  " << sensor_range << std::endl;
		std::cout << " particle " << j << " p_X is " << p_X << std::endl;
		std::cout << " particle " << j << " p_Y is " << p_Y << std::endl;
		std::cout << " particle " << j << " p_theta is " << p_theta << std::endl;
		*/
		int k = 0;
		for (int i = 0; i < map_landmarks.landmark_list.size(); i++) {
			double l_X = map_landmarks.landmark_list[i].x_f;
			double l_Y = map_landmarks.landmark_list[i].y_f;
			int l_id = i;
			d = dist(l_X, l_Y, p_X, p_Y);
			if (d <= sensor_range) {
				
				//Push back predicts with landmark in sensor range.
				predicted.x = l_X;
				predicted.y = l_Y;
				predicted.id = l_id;
				predicts.push_back(predicted);
				/*
				std::cout << " particle " << j << " predicted.x is " << predicts[k].x << std::endl;
				std::cout << " particle " << j << " predicted.y is " << predicts[k].x << std::endl;
				std::cout << " particle " << j << " landmark index is " << predicts[k].id << std::endl;
				*/
				k += 1;

			}

		}
		/*
		std::cout << "number of prediction for particle " << j << " is :" << predicts.size() << std::endl;
		std::cout << "number of observation is :" << observations.size() << std::endl;
		*/
		observations_map.clear();
		for (int i = 0; i < observations.size(); i++) {
			/*
			std::cout << " particle " << j << " observations[i].x is " << observations[i].x << std::endl;
			std::cout << " particle " << j << " observations[i].y is " << observations[i].y << std::endl;
			*/
			// Trnasfom obseravation into Map coordinate

			LandmarkObs obp_map = Transformations(observations[i], p_X, p_Y, p_theta);

			observations_map.push_back(obp_map);
			/*
			std::cout << " particle " << j << " observations_map[i].x after transform is " << observations_map[i].x << std::endl;
			std::cout << " particle " << j << " observations_map[i].y after transform is " << observations_map[i].y << std::endl;
			*/
		}
		std::string test;

		dataAssociation(predicts, observations_map);

		// define inputs
		double sig_x, sig_y, x_obs, y_obs, mu_x, mu_y, weight;
		int land_mark_id;

		// defined final weight
		double final_weight = 1.;
		for (int i = 0; i < observations_map.size(); i++) {
			// OBS1 values
			sig_x = std_landmark[0];
			sig_y = std_landmark[1];
			x_obs = observations_map[i].x;
			y_obs = observations_map[i].y;
			land_mark_id = observations_map[i].id;
			mu_x = map_landmarks.landmark_list[land_mark_id].x_f;
			mu_y = map_landmarks.landmark_list[land_mark_id].y_f;
			/*
			std::cout << " ma_obs " << i << " x_obs is " << x_obs << std::endl;
			std::cout << " ma_obs " << i << " y_obs is " << y_obs << std::endl;

			std::cout << " map_landmarks " << land_mark_id << " mu_x is " << mu_x << std::endl;
			std::cout << " map_landmarks " << land_mark_id << " mu_y is " << mu_y << std::endl;
			*/
			// Calculate OBS1 weight
			weight = multiv_prob(sig_x, sig_y, x_obs, y_obs, mu_x, mu_y);
			
			// std::cout << "weight of observartion "  << i << " particle " << j << " is :" << weight << std::endl;

			final_weight *= weight;
		}
		sum_particle_weights += final_weight;
		wi.push_back(final_weight);
		// std::cout << "Final weight of particle " << j << " is :" << final_weight << std::endl;
	}

	// normalize wights
	for (int j = 0; j < num_particles; ++j) {
		particles[j].weight = wi[j] / sum_particle_weights;
		weights[j] = particles[j].weight;
	}


}

void ParticleFilter::resample() {
  /**
   * Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   
	*/ 
	// define random number generator
	std::default_random_engine generator;

	// define discrete_distribution based on particle whights
	std::discrete_distribution<int> distribution(weights.begin(), weights.end());

	std::vector<Particle> Reasmpled_particles;
	Particle Reasmpled_particle;
	Reasmpled_particles.clear();
	// roll with the number of particles and resample based on discrete_distribution
	for (int i = 0; i < num_particles; ++i) {
		// generate a particle index
		int partic_id = distribution(generator);
		// std::cout << "resample id is " << partic_id << " with weigh of " << weights[partic_id] << std::endl;
		//Vector now has 1 element @ index 0, so modify it.
		Reasmpled_particle = particles[partic_id];
		Reasmpled_particle.id = i;

		Reasmpled_particles.push_back(Reasmpled_particle);
	}
	particles.clear();
	// replace particles with resampled particles
	particles = Reasmpled_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}