#ifndef HELPERS_H
#define HELPERS_H

#include <math.h>
#include <string>
#include <vector>

// for convenience
using std::string;
using std::vector;

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
//   else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

//
// Helper functions related to waypoints and converting from XY to Frenet
//   or vice versa
//

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Calculate distance between two points
double distance(double x1, double y1, double x2, double y2) {
  return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

// Calculate closest waypoint to current x, y position
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, 
                    const vector<double> &maps_y) {
  double closestLen = 100000; //large number
  int closestWaypoint = 0;

  for (int i = 0; i < maps_x.size(); ++i) {
    double map_x = maps_x[i];
    double map_y = maps_y[i];
    double dist = distance(x,y,map_x,map_y);
    if (dist < closestLen) {
      closestLen = dist;
      closestWaypoint = i;
    }
  }

  return closestWaypoint;
}

// Returns next waypoint of the closest waypoint
int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, 
                 const vector<double> &maps_y) {
  int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

  double map_x = maps_x[closestWaypoint];
  double map_y = maps_y[closestWaypoint];

  double heading = atan2((map_y-y),(map_x-x));

  double angle = fabs(theta-heading);
  angle = std::min(2*pi() - angle, angle);

  if (angle > pi()/2) {
    ++closestWaypoint;
    if (closestWaypoint == maps_x.size()) {
      closestWaypoint = 0;
    }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, 
                         const vector<double> &maps_x, 
                         const vector<double> &maps_y) {
  int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

  int prev_wp;
  prev_wp = next_wp-1;
  if (next_wp == 0) {
    prev_wp  = maps_x.size()-1;
  }

  double n_x = maps_x[next_wp]-maps_x[prev_wp];
  double n_y = maps_y[next_wp]-maps_y[prev_wp];
  double x_x = x - maps_x[prev_wp];
  double x_y = y - maps_y[prev_wp];

  // find the projection of x onto n
  double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
  double proj_x = proj_norm*n_x;
  double proj_y = proj_norm*n_y;

  double frenet_d = distance(x_x,x_y,proj_x,proj_y);

  //see if d value is positive or negative by comparing it to a center point
  double center_x = 1000-maps_x[prev_wp];
  double center_y = 2000-maps_y[prev_wp];
  double centerToPos = distance(center_x,center_y,x_x,x_y);
  double centerToRef = distance(center_x,center_y,proj_x,proj_y);

  if (centerToPos <= centerToRef) {
    frenet_d *= -1;
  }

  // calculate s value
  double frenet_s = 0;
  for (int i = 0; i < prev_wp; ++i) {
    frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
  }

  frenet_s += distance(0,0,proj_x,proj_y);

  return {frenet_s,frenet_d};
}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, 
                     const vector<double> &maps_x, 
                     const vector<double> &maps_y) {
  int prev_wp = -1;

  while (s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1))) {
    ++prev_wp;
  }

  int wp2 = (prev_wp+1)%maps_x.size();

  double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),
                         (maps_x[wp2]-maps_x[prev_wp]));
  // the x,y,s along the segment
  double seg_s = (s-maps_s[prev_wp]);

  double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
  double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

  double perp_heading = heading-pi()/2;

  double x = seg_x + d*cos(perp_heading);
  double y = seg_y + d*sin(perp_heading);

  return {x,y};
}
vector<vector<double>> check_cars(int lane,double car_s, int &previous_size, const vector<vector<double>> &sensor_fusion){
  vector<vector<double>> adjacent_cars;
  for(int i = 0; i<sensor_fusion.size();i++){
    float d = sensor_fusion[i][6];
    if((d<(2+4*lane+2)) && (d>(2+4*lane-2))){
      double check_car_s = sensor_fusion[i][5];
      double vx = sensor_fusion[i][3];
      double vy = sensor_fusion[i][4];
      double check_speed = sqrt(vx*vx+vy*vy);
      check_car_s += ((double)previous_size*0.02*check_speed);
      
      if(fabs(check_car_s - car_s)<60){
        adjacent_cars.push_back(sensor_fusion[i]);
      } 
    }
  }
  return adjacent_cars;
}

// calculate the inefficiency cost according to the lane speed
float inefficiency_cost(int next_lane, double &ref_v,double car_s, int &previous_size, const vector<vector<double>> &sensor_fusion){
  int i_cost = 0;
  vector<vector<double>> adjacent_cars = check_cars(next_lane,car_s,previous_size, sensor_fusion);
  
  //if the lane is empty, cost = 0
  if(adjacent_cars.size() ==0){
    return i_cost;
  }
  //else check if there's any car in the front and if closest car in the front is faster
  else{
    double min_dis = 100000000;
    vector<double> target_car = {};
    double target_speed;
    for(int j=0; j<adjacent_cars.size(); j++){
      double check_car_s = adjacent_cars[j][5];
      double vx = adjacent_cars[j][3];
      double vy = adjacent_cars[j][4];
      double check_speed = sqrt(vx*vx+vy*vy);
      check_car_s += ((double)previous_size*0.02*check_speed);
      if((check_car_s > car_s) && (check_car_s - car_s)<min_dis){
        min_dis = check_car_s - car_s;
        target_car =  adjacent_cars[j];
        target_speed = check_speed;
      }
    }
    if(target_car.size() ==0){
      return i_cost;
    }
    else{
      //if target lane is faster, add cost according to the speed
      if(target_speed>ref_v){
        return i_cost += ref_v/target_speed;
      }
      //if the target lane is slower, add large cost
      else{
        return i_cost +=300;
      }
    }
  }
}

//calculate the safety cost 
float safety_cost(int &lane, double ref_v,int target_lane, double car_s, int &previous_size, const vector<vector<double>> &sensor_fusion){
  int s_cost = 0;
  //if the target lane is out of the 3 right lane, then add very big cost to avoid accident
  if(target_lane<0 || target_lane>2){
    return s_cost += 1000000;
  }
  // found adjacent cars
  vector<vector<double>> adjacent_cars = check_cars(target_lane,car_s,previous_size, sensor_fusion);
  
  //check if there's any car is too close to ego car
  for(int j=0; j<adjacent_cars.size(); j++){
    double check_car_s = adjacent_cars[j][5];
    double vx = adjacent_cars[j][3];
    double vy = adjacent_cars[j][4];
    double check_speed = sqrt(vx*vx+vy*vy);
    check_car_s += ((double)previous_size*0.02*check_speed);
    //check if it is safe to change lane
    if(fabs(check_car_s-car_s)<15){
      // make keep lane has better safety cost
      if(target_lane == lane){
        return s_cost += 100;
      }
      else if(check_speed<ref_v || check_car_s-car_s > -10){
        return s_cost += 50;
      }
      return s_cost += 10000;
    }
  }
  return s_cost;
}

float choose_next_lane(int &lane,double ref_v, double car_s, int &previous_size, const vector<vector<double>> &sensor_fusion){
  int next_lane;
  double min_cost=100000000000.0;
  for(int i=-1; i<2; i++){
    int target_lane = lane + i;
    double cost = inefficiency_cost(target_lane, ref_v,car_s,previous_size, sensor_fusion) + safety_cost(lane, ref_v, target_lane, car_s, previous_size, sensor_fusion);
    if(cost<min_cost){
      min_cost = cost;
      next_lane = target_lane;
    }
  }
  return next_lane;
  
}

#endif  // HELPERS_H