#include "fe_lens.hpp"



double vis_scale = 1.0;

NewtonRaphson::NewtonRaphson(double tol, unsigned int max_iter) : tol_(tol), max_iter_(max_iter) {}


double NewtonRaphson::solve(double x0, std::function<double (double)> (f), std::function<double (double)> (f_prime)) {
  double x = x0;

  //std::cout << "Newton-Raphson: " << x << std::endl;
    
  for (unsigned int i=0; i<max_iter_; i++) {
    double fx = f(x);
    double dfx = f_prime(x);
    double x1 = x - fx / dfx;
    //std::cout << "NR " << i << ":\t" << x << " " << fx << " " << dfx << " " << x1 << std::endl;
    if (std::abs(x1 - x) < tol_) {
      return x1;
    }
    x = x1;
  }
  return std::numeric_limits<double>::quiet_NaN(); // return NaN if the method doesn't converge
}



FisheyeLens::FisheyeLens(double fx, double fy, double cx, double cy,
              double k1, double k2, double k3, double k4): 
              fx_(fx), fy_(fy), cx_(cx), cy_(cy), k1_(k1), k2_(k2), k3_(k3), k4_(k4) {}


double FisheyeLens::RTheta(double theta) {
  return k1_*theta + k2_*pow(theta,3) + k3_*pow(theta,5) + k4_*pow(theta,7);
}


double FisheyeLens::Rd(double theta) {
  return theta * ( 1 + k1_*pow(theta,2) + k2_*pow(theta,4) + k3_*pow(theta,6) + k4_*pow(theta,8) );
}


double FisheyeLens::RThetaInv(double r_theta, double x0) {
  //std::cout << "RThetaInv: " << r_theta << std::endl;
//  auto func = [this, r_theta](double th) {
//    return k1_*pow(th,1) + k2_*pow(th,3) + k3_*pow(th,5) + k4_*pow(th,7) - r_theta;
//  };
//
//  auto func_prime = [this](double th) {
//    return k1_*1*pow(th,0) - k2_*3*pow(th,2) - k3_*5*pow(th,4) - k4_*7*pow(th,6);
//  };

  double r_d = r_theta;

  auto func = [this, r_d](double th) {
    return  th * ( 1.0 + k1_*pow(th,2) + k2_*pow(th,4) + k3_*pow(th,6) + k4_*pow(th,8)) - r_d;
  };

  auto func_prime = [this](double th) {
    return  1.0 + 3*k1_*pow(th,2) + 5*k2_*pow(th,4) + 7*k3_*pow(th,6) + 9*k4_*pow(th,8);
  };

  NewtonRaphson solver(1e-6, 1000);
  double theta_solver = solver.solve(x0, func, func_prime);

  return theta_solver;
}



cv::Point2f FisheyeLens::Compute2D(double theta, double phi, bool world_coord) {
  double r_d = Rd(theta);
  double x_d = r_d * cos(phi);
  double y_d = r_d * sin(phi);
  //std::cout << "r_d: " << r_d << " " << x_d << " " << y_d << std::endl;

  if (world_coord) {
    return cv::Point2f(x_d, y_d);
  } else {
    double x = fx_ * x_d + cx_;
    double y = fy_ * y_d + cy_;
    return cv::Point2f(x, y);
  }

  //  double r_theta = RTheta(theta);
  //  double x = r_theta * cos(phi);
  //  double y = r_theta * sin(phi);
  //  return cv::Point2f(x, y);
}

/*
std::vector<double> FisheyeLens::Compute3D(double x, double y, bool sim, double x0) {

  if (x > 0.99) {
    x = 0.99;
  } else if (x < -0.99) {
    x = -0.99;
  }

  if (y > 0.99) {
    y = 0.99;
  } else if (y < -0.99) {
    y = -0.99;
  }
  
  double phi = 0.0;
  double theta = 0.0;

  if (x == 0.0 && y == 0.0) {
    return std::vector<double>({theta, phi});
  }

  if (x == 0.0) {
    if (y >= 0.0) {
      phi = M_PI / 2.0;
    } else {
      phi = -M_PI / 2.0;
    }
  } else {
    phi = atan(y / x);
  }

  double r_theta = std::abs(x / cos(phi));
  theta = RThetaInv(r_theta, x0);
  std::cout << "( " << x << ", " << y << " ) -> [ " << theta << ", " << phi << "] " << r_theta << std::endl;
  return std::vector<double>({theta, phi});
}
*/



std::vector<double> FisheyeLens::Compute3D(double x, double y, bool world_coord, double x0) {  
  double phi = 0.0;
  double theta = 0.0;

  double x_d = x;
  double y_d = y;

  // If we're in camera coordinate system, convert to world coordinate system
  if (!world_coord) {
    x_d = (x - cx_) / fx_;
    y_d = (y - cy_) / fy_;
  }

  // If the point is at the center of the image, return 0,0
  if (x_d == 0.0 && y_d == 0.0) {
    return std::vector<double>({theta, phi});
  }

  // If the point is on the y-axis, set phi to pi/2 or -pi/2
  if (x_d == 0.0) {
    phi = (y_d >= 0.0) ? (M_PI / 2.0) : (-M_PI / 2.0);
  } else {
    phi = atan(y_d / x_d);
  }

  // If the point is on the negative x-axis, sum 180 degrees to phi
  if (x_d < 0.0) {
    phi += M_PI;
  }

  // Compute the distance from the center of the image to the distorted point
  double r_d = sqrt(pow(x_d, 2) + pow(y_d, 2));

  theta = RThetaInv(r_d, x0);

  //std::cout << "( " << x << ", " << y << " ) -> " 
  //          << "( " << x_d << ", " << y_d << " ) -> " 
  //          << "[ " << theta << ", " << phi << "] " 
  //          << r_d << std::endl;
            
  return std::vector<double>({theta, phi});
}


double FisheyeLens::ComputeError(std::vector<std::vector<double>> v1, std::vector<std::vector<double>> v2) {
  double error_theta = 0.0;
  double error_phi = 0.0;

  for (unsigned int i = 0; i < v1.size(); i++) {
    double error_theta_i = pow(v1[i][0] - v2[i][0], 2);
    double error_phi_i = pow(v1[i][1] - v2[i][1], 2);
    error_theta += error_theta_i;
    error_phi += error_phi_i;
  }
  error_theta = sqrt(error_theta / v1.size());
  error_phi = sqrt(error_phi / v1.size());
  std::cout << "Error theta: " << error_theta << std::endl;
  std::cout << "Error phi: " << error_phi << std::endl;

  return sqrt(pow(error_theta, 2) + pow(error_phi, 2));
}



double FisheyeLens::f() {
  double lx = 2 * cx_;
  double ly = 2 * cy_;

  double f = (lx / (lx + ly)) * fx_ + (ly / (lx + ly)) * fy_;

  return f;
}



Visualizer::Visualizer(std::string window_name, double scale):
  window_name_(window_name), scale_(scale) {
  //cv::viz::Viz3d window_(window_name_);
  //cv::namedWindow(window_name_, cv::WINDOW_AUTOSIZE);
}


void Visualizer::AddCloud(std::vector<cv::Point3f> cloud, std::vector<cv::Vec3b> color) {
  point_clouds_.push_back(cloud);
  colors_.push_back(color);
}


void Visualizer::AddCloud(std::vector<cv::Point3f> cloud, cv::Vec3b color) {
  std::vector<cv::Vec3b> colors;
  for (unsigned int i = 0; i < cloud.size(); i++) {
    colors.push_back(color);
  }
  AddCloud(cloud, colors);
}


void Visualizer::Render() {
  window_ = cv::viz::Viz3d(window_name_);

  for (unsigned int i=0; i<point_clouds_.size(); i++) {
    std::string name = "cloud" + std::to_string(i);
    cv::viz::WCloud cloud(point_clouds_[i], colors_[i]);
    window_.showWidget(name, cloud);
  }

  // Add a coordinate system
  cv::viz::WCoordinateSystem world_coor(vis_scale);
  window_.showWidget("World", world_coor);

  window_.spin();
}



