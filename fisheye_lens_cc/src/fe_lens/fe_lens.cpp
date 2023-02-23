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


double FisheyeLens::FocalLength() {
  double lx = 2 * cx_;
  double ly = 2 * cy_;

  double f = (lx / (lx + ly)) * fx_ + (ly / (lx + ly)) * fy_;

  return f;
}


//cv::Point2f FisheyeLens::Compute2D(double theta, double phi) {
//  double r_theta = RTheta(theta);
//  double x = r_theta * cos(phi);
//  double y = r_theta * sin(phi);
//  return cv::Point2f(x, y);
//}


cv::Point2f FisheyeLens::Compute2D(double theta, double phi, bool sim) {
  double r_d = Rd(theta);
  double x_d = r_d * cos(phi);
  double y_d = r_d * sin(phi);
  std::cout << "r_d: " << r_d << " " << x_d << " " << y_d << std::endl;
  double x = fx_ * x_d + cx_;
  double y = fy_ * y_d + cy_;

  // We want to compare in world coordinates
  x = x_d;
  y = y_d;

  return cv::Point2f(x, y);
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



std::vector<double> FisheyeLens::Compute3D(double x, double y, bool sim, double x0) {  
  double phi = 0.0;
  double theta = 0.0;

  double x_d = (x - cx_) / fx_;
  double y_d = (y - cy_) / fy_;

  if (x_d == 0.0) {
    if (y_d >= 0.0) {
      phi = M_PI / 2.0;
    } else {
      phi = -M_PI / 2.0;
    }
  } else {
    phi = atan(y_d / x_d);
  }


  // We want to compare in world coordinates
  if (sim) {
    x_d = x;
    y_d = y;
  }

  if (x_d == 0.0 && y_d == 0.0) {
    return std::vector<double>({theta, phi});
  }

  if (x_d == 0.0) {
    if (y_d >= 0.0) {
      phi = M_PI / 2.0;
    } else {
      phi = -M_PI / 2.0;
    }
  } else {
    phi = atan(y_d / x_d);
  }

  if (x_d < 0.0) {
    phi += M_PI;
  }

  double r_d = sqrt(pow(x_d, 2) + pow(y_d, 2));

  theta = RThetaInv(r_d, x0);

  std::cout << "( " << x << ", " << y << " ) -> " 
            << "( " << x_d << ", " << y_d << " ) -> " 
            << "[ " << theta << ", " << phi << "] " 
            << r_d << std::endl;
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




double FocalLength(const cv::Mat &K) {
  double fx = K.at<double>(0, 0);
  double fy = K.at<double>(1, 1);
  double cx = K.at<double>(0, 2);
  double cy = K.at<double>(1, 2);

  double lx = 2 * cx;
  double ly = 2 * cy;

  double f = (lx / (lx + ly)) * fx + (ly / (lx + ly)) * fy;

  return f;
}


cv::Point3f compute3dOverLens(cv::Point2f p, float f, const cv::Mat &D) {
  float r = sqrt(pow(p.x, 2) + pow(p.y, 2));
  float theta_d = atan2(r, f);
  float theta = theta_d / (1 + D.at<double>(0, 0)*pow(theta_d,2) + D.at<double>(1, 0)*pow(theta_d,4) + D.at<double>(2, 0)*pow(theta_d,6) + D.at<double>(3, 0)*pow(theta_d,8));


  // Compute phi, the angle between the x-axis and the point
  float phi = atan2(p.y, p.x);

  float x = f * sin(theta) * cos(phi);
  float y = f * sin(theta) * sin(phi);
  float z = f * cos(theta);

  return cv::Point3f(f * sin(theta) * cos(phi), f * sin(theta) * sin(phi), f * cos(theta));
  return cv::Point3f(-x, -y, z);
}




cv::Point3f compute3dOverLens2(cv::Point2f p, float f, const cv::Mat &D) {

  double d1 = D.at<double>(0, 0);
  double d2 = D.at<double>(1, 0);
  double d3 = D.at<double>(2, 0);
  double d4 = D.at<double>(3, 0);

  // Compute phi, the angle between the x-axis and the point
  float phi = atan2(p.y, p.x);

  float r_theta = p.x / (f * cos(phi));

  // Use NewtonRaphson to solve for theta
  float theta = r_theta;
  float theta_new = theta - (theta - d1*pow(theta,1) - d2*pow(theta,3) - d3*pow(theta,5) - d4*pow(theta,7) - r_theta) / (1 - d1*1*pow(theta,0) - d2*3*pow(theta,2) - d3*5*pow(theta,4) - d4*7*pow(theta,6));

  while (abs(theta_new - theta) > 0.0001) {
    theta = theta_new;
    theta_new = theta - (theta - d1*pow(theta,1) - d2*pow(theta,3) - d3*pow(theta,5) - d4*pow(theta,7) - r_theta) / (1 - d1*1*pow(theta,0) - d2*3*pow(theta,2) - d3*5*pow(theta,4) - d4*7*pow(theta,6));
  }




/*
  auto func = [d1, d2, d3, d4](double th) {
    std::cout << d1 << "\t" << d2 << "\t" << d3 << "\t" << d4 << "\t" << th << std::endl;
    return d1*pow(th,1) + d2*pow(th,3) + d3*pow(th,5) + d4*pow(th,7);
  };

  auto func_prime = [d1, d2, d3, d4](double th) {
    return d1*1*pow(th,0) - d2*3*pow(th,2) - d3*5*pow(th,4) - d4*7*pow(th,6);
  };


  NewtonRaphson solver(0.1, 1e-6, 10000);
  double theta_solver = solver.solve(func, func_prime);


  std::cout << "theta: " << theta << "\t" << theta_solver << std::endl;
*/








  float x = (f/2) * sin(theta) * cos(phi);
  float y = (f/2) * sin(theta) * sin(phi);
  float z = (f/2) * cos(theta);






/*

  float r = sqrt(pow(p.x, 2) + pow(p.y, 2));
  float theta_d = atan2(r, f);
  float theta = theta_d / (1 + D.at<double>(0, 0)*pow(theta_d,2) + D.at<double>(1, 0)*pow(theta_d,4) + D.at<double>(2, 0)*pow(theta_d,6) + D.at<double>(3, 0)*pow(theta_d,8));


  float x = f * sin(theta) * cos(phi);
  float y = f * sin(theta) * sin(phi);
  float z = f * cos(theta);

  return cv::Point3f(f * sin(theta) * cos(phi), f * sin(theta) * sin(phi), f * cos(theta));
  return cv::Point3f(-x, -y, z);
*/
  return cv::Point3f(x,y,z);
}




cv::Point2f compute2dFromLens(cv::Point3f p, float f, float theta, float phi, const cv::Mat &D) {
  float d1 = D.at<double>(0, 0);
  float d2 = D.at<double>(1, 0);
  float d3 = D.at<double>(2, 0);
  float d4 = D.at<double>(3, 0);

  float r = f * (d1*pow(theta,1) + d2*pow(theta,3) + d3*pow(theta,5) + d4*pow(theta,7));

  return cv::Point2f(r*cos(phi), r*sin(phi));
}


// Compute the 3D surface of a Kannala-Brandt fisheye lens
std::vector<cv::Point3f> computeFisheyeSurface(std::vector<cv::Point2f> points, const cv::Mat &K, const cv::Mat &D, int width, int height){
  float f = FocalLength(K);

  // Compute the surface of the fisheye lens from the undistorted points
  std::vector<cv::Point3f> points3d;
  for (unsigned int i = 0; i < points.size(); i++) {
    points3d.push_back(compute3dOverLens2(points[i], f, D));
  }

  return points3d;
}



std::vector<cv::Point3f> computeImageFromFisheyeSurface(std::vector<cv::Point3f> points, const cv::Mat &K, const cv::Mat &D, int width, int height) {
  float f = FocalLength(K);


  // Compute the surface of the fisheye lens from the undistorted points
  std::vector<cv::Point3f> points3d;
  for (unsigned int i = 0; i < points.size(); i++) {

    // Phi from any of the other components
    //float phi = acos(points[i].x / (f * sin(theta)));

    float d1 = D.at<double>(0, 0);
    float d2 = D.at<double>(1, 0);
    float d3 = D.at<double>(2, 0);
    float d4 = D.at<double>(3, 0);
    float fx = K.at<double>(0, 0);
    float fy = K.at<double>(1, 1);
    float cx = K.at<double>(0, 2);
    float cy = K.at<double>(1, 2);

    // Theta from third component
    float theta = acos(points[i].z / f);

    float xc = points[i].x - cx;
    float yc = points[i].y - cy;

    float theta_d = theta * (1+ d1*pow(theta,2) + d2*pow(theta,4) + d3*pow(theta,6) + d4*pow(theta,8));
    float r = f * theta;
    float x_d = theta_d * xc / r;
    float y_d = theta_d * yc / r;
    float u = fx * x_d + cx;
    float v = fy * y_d + cy;

    points3d.push_back(cv::Point3f(u, v, -2*f));
  }

  return points3d;
}



std::vector<cv::Point3f> computeImageFromFisheyeSurface2(std::vector<cv::Point3f> points, const cv::Mat &K, const cv::Mat &D, int width, int height) {
  float f = FocalLength(K);


  // Compute the surface of the fisheye lens from the undistorted points
  std::vector<cv::Point3f> points3d;
  for (unsigned int i = 0; i < points.size(); i++) {

    // Phi from any of the other components
    //float phi = acos(points[i].x / (f * sin(theta)));

    float d1 = D.at<double>(0, 0);
    float d2 = D.at<double>(1, 0);
    float d3 = D.at<double>(2, 0);
    float d4 = D.at<double>(3, 0);
    //float fx = K.at<double>(0, 0);
    //float fy = K.at<double>(1, 1);
    //float cx = K.at<double>(0, 2);
    //float cy = K.at<double>(1, 2);


    // Applying the distortion model for computing r(theta) and (x,y) from (theta, phi)
    float r_3d = sqrt(pow(points[i].x, 2) + pow(points[i].y, 2));
    float theta = atan2(r_3d, points[i].z);

    float phi = atan2(points[i].y, points[i].x);

    float r_theta = f * (d1*pow(theta,1) + d2*pow(theta,3) + d3*pow(theta,5) + d4*pow(theta,7));
    
    float x = r_theta * cos(phi);
    float y = r_theta * sin(phi);
    
    points3d.push_back(cv::Point3f(x, y, -f));
  }

  return points3d;
}



std::vector<cv::Point3f> compute3DFrom2D(std::vector<cv::Point2f> points, const cv::Mat &K, const cv::Mat &D, int width, int height) {
  float f = FocalLength(K);


  // Compute the surface of the fisheye lens from the undistorted points
  std::vector<cv::Point3f> points3d;
  for (unsigned int i = 0; i < points.size(); i++) {

    float d1 = D.at<double>(0, 0);
    float d2 = D.at<double>(1, 0);
    float d3 = D.at<double>(2, 0);
    float d4 = D.at<double>(3, 0);
    float fx = K.at<double>(0, 0);
    float fy = K.at<double>(1, 1);
    float cx = K.at<double>(0, 2);
    float cy = K.at<double>(1, 2);
    
    float u = points[i].x;
    float v = points[i].y;

    float x_d = (u - cx) / fx;
    float y_d = (v - cy) / fy;

    float r = sqrt(pow(x_d, 2) + pow(y_d, 2));
    float theta_d = r / f;
    float theta = theta_d / (1+ d1*pow(theta_d,2) + d2*pow(theta_d,4) + d3*pow(theta_d,6) + d4*pow(theta_d,8));

    float phi = atan2(y_d, x_d);

    float x = f * sin(theta) * cos(phi);
    float y = f * sin(theta) * sin(phi);
    float z = f * cos(theta);

    points3d.push_back(cv::Point3f(x, y, 3*z));
  }

  return points3d;
}





void display3dSurface(std::vector<cv::Point3f> points, std::vector<cv::Vec3b> colors) {
  cv::viz::Viz3d window("3D Surface");

  if (colors.size() > 0) {
    cv::viz::WCloud cloud(points, colors);
    window.showWidget("cloud", cloud);
  } else {
    cv::viz::WCloud cloud(points, cv::viz::Color::white());
    window.showWidget("cloud", cloud);
  }
  
  // Add a coordinate system
  cv::viz::WCoordinateSystem world_coor(vis_scale);
  window.showWidget("World", world_coor);

  window.spin();
}



void display3dSurfaceAndImage(std::vector<std::vector<cv::Point3f>> points, std::vector<cv::Vec3b> colors) {
  cv::viz::Viz3d window("3D Surface");
  
  for (unsigned int i=0; i<points.size(); i++) {
    std::string name = "cloud" + std::to_string(i);
    cv::viz::WCloud cloud(points[i], colors);
    window.showWidget(name, cloud);
  }
  
  // Add a coordinate system
  cv::viz::WCoordinateSystem world_coor(vis_scale);
  window.showWidget("World", world_coor);
  
  window.spin();
}




void display3dSurfaceAndImage(std::vector<std::vector<cv::Point3f>> points, std::vector<std::vector<cv::Vec3b>> colors) {
  cv::viz::Viz3d window("3D Surfaces");
  
  for (unsigned int i=0; i<points.size(); i++) {
    std::string name = "cloud" + std::to_string(i);
    cv::viz::WCloud cloud(points[i], colors[i]);
    window.showWidget(name, cloud);
  }

  // Add a coordinate system
  cv::viz::WCoordinateSystem world_coor(vis_scale);
  window.showWidget("World", world_coor);

  window.spin();
}
