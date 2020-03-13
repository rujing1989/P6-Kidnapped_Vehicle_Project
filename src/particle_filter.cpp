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
static std::default_random_engine gen;
void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 501;  // TODO: Set the number of 
  double std_x=std[0];
  double std_y=std[1];
  double std_theta=std[2];
  static std::normal_distribution<double> dist_x(0,std_x);
  static std::normal_distribution<double> dist_y(0,std_y);
  static std::normal_distribution<double> dist_theta(0,std_theta);
  for(int i=0;i<num_particles;i++){
    Particle p;
    p.id=i;
    p.x=x;
    p.y=y;
    p.theta=theta;
    p.weight=1.0;
    p.x+=dist_x(gen);
    p.y+=dist_y(gen);
    p.theta+=dist_theta(gen);
    
    particles.push_back(p);// why no declare before 
  }
  is_initialized= true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  double std_x=std_pos[0];
  double std_y=std_pos[1];
  double std_theta=std_pos[2];
  static std::normal_distribution<double> dist_x(0,std_x);
  static std::normal_distribution<double> dist_y(0,std_y);
  static std::normal_distribution<double> dist_theta(0,std_theta);
  for(int i=0;i<num_particles;i++){
    if(fabs(yaw_rate)<0.00001){
      particles[i].x=particles[i].x+velocity*delta_t*cos(particles[i].theta);
      particles[i].y=particles[i].y+velocity*delta_t*sin(particles[i].theta);
    }
    else{
      particles[i].x=particles[i].x+velocity/yaw_rate*(sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta));
    	particles[i].y=particles[i].y+velocity/yaw_rate*(cos(particles[i].theta)-cos(particles[i].theta+yaw_rate*delta_t));
    	particles[i].theta=particles[i].theta+yaw_rate*delta_t;
    }
    particles[i].x+=dist_x(gen);
    particles[i].y+=dist_y(gen);
    particles[i].theta+=dist_theta(gen);
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  for (unsigned int i = 0;i<observations.size();i++){
    LandmarkObs o=observations[i];
    double min_dist=std::numeric_limits<double>::max();//init munimum distance to maximun possible,and iterativ find smaller
    int mapID=-1;

    for(unsigned int j=0;j<predicted.size();j++){
      LandmarkObs p=predicted[j];
      double Distance=dist(o.x,o.y,p.x,p.y);
      //if the distance is less than min, stored the id and update min.
      if(Distance<min_dist){
        min_dist=Distance;
        mapID=p.id;
      }  
    }
    //update the observation
    observations[i].id=mapID;
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
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
  for (int i=0;i<num_particles;i++){
    double p_x=particles[i].x;
    double p_y=particles[i].y;
    double p_theta=particles[i].theta;
    //creata a vector to hold the map landmark locations predicted to be within sensor range of particle; vector predictions contain the Landmark,near by particles(map coordinate)
    vector<LandmarkObs> predictions;
    for(unsigned int j=0;j<map_landmarks.landmark_list.size();j++){
      float lm_x=map_landmarks.landmark_list[j].x_f;
      float lm_y=map_landmarks.landmark_list[j].y_f;
      int lm_id=map_landmarks.landmark_list[j].id_i;
      //only consider landmarks within sensor range of the particle(rather than using the 'dist' considering a circular)
      if((fabs(lm_x-p_x)<=sensor_range)&&(fabs(lm_y-p_y)<=sensor_range)){
        predictions.push_back(LandmarkObs{lm_id,lm_x,lm_y});
      }
    }
    //create and populate a copy of the list of observations transformed from vehicle coordinates to map coordinates
    vector<LandmarkObs> transformed_os;
    for (unsigned int j=0;j<observations.size();j++){
      double t_x=cos(p_theta)*observations[j].x-sin(p_theta)*observations[j].y+p_x;
      double t_y=sin(p_theta)*observations[j].x+cos(p_theta)*observations[j].y+p_y;
      transformed_os.push_back(LandmarkObs{observations[j].id,t_x,t_y});
    }
    dataAssociation(predictions,transformed_os);
    particles[i].weight=1.0;
    for(unsigned int j=0;j<transformed_os.size();j++){
      double detected_x,detected_y,pr_x,pr_y;
      detected_x=transformed_os[j].x;
      detected_y=transformed_os[j].y;
      int associated_prediction=transformed_os[j].id;
      //get the x,y coordinates of the vector <predictions>
      for (unsigned int k=0;k<predictions.size();k++){
        //the  real landmark near the particleand ?==detected landmark
        if(predictions[k].id==associated_prediction){
          pr_x=predictions[k].x;
          pr_y=predictions[k].y;
        }
      }
      double s_x=std_landmark[0];
      double s_y=std_landmark[1];
      double obs_w = ( 1/(2*M_PI*s_x*s_y)) * exp( -( pow(pr_x-detected_x,2)/(2*pow(s_x, 2)) + (pow(pr_y-detected_y,2)/(2*pow(s_y, 2))) ) );
      particles[i].weight*=obs_w;
    }
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  vector<Particle> new_particles;
  vector<double> weights;
  for(int i=0;i<num_particles;i++){
    weights.push_back(particles[i].weight);
  }
  //generate random starting index for resampling wheel
  std::uniform_int_distribution<int> uniinitdist(0,num_particles-1);
  auto index=uniinitdist(gen);
  double max_weight=*std::max_element(weights.begin(),weights.end());
  //unifor radom distribution[0.0,max_weight)
  std::uniform_real_distribution<double> uniweight(0.0,max_weight);
  double beta=0.0;
  //spin the resampling wheel
  for(int i=0;i<num_particles;i++){
    beta+=uniweight(gen)*2.0;
    while (beta>weights[index]){
      beta-=weights[index];
      index=(index+1)%num_particles;
    }
    new_particles.push_back(particles[index]);
  }
  particles=new_particles;
  

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