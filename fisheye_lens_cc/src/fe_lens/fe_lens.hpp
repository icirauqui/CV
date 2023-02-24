#ifndef FE_LENS_HPP
#define FE_LENS_HPP

#include <iostream>
#include <vector>
#include <cmath>

#include <opencv2/opencv.hpp>
#include <opencv2/viz.hpp> // include the viz module

#include <boost/math/tools/polynomial.hpp>
#include <boost/math/tools/roots.hpp>


class NewtonRaphson {

// Instructions for Newton-Raphson solver
//   double f(double x) { return x * x - 2; }
//   double f_prime(double x) { return 2 * x; }
//   NewtonRaphson solver(1, 1e-6, 100);
//   double x = solver.solve(f, f_prime);

public:
  NewtonRaphson(double tol, unsigned int max_iter);

  double solve(double x0, std::function<double (double)> (f), std::function<double (double)> (f_prime));

private:
  double tol_;
  unsigned int max_iter_;
};




class FisheyeLens {
public:
  FisheyeLens(double fx, double fy, double cx, double cy,
              double k1, double k2, double k3, double k4);

  double RTheta(double theta);

  double Rd(double theta);

  // Inverse of RTheta with Newton-Raphson
  double RThetaInv(double r_theta, double x0 = 0.1);

  cv::Point2f Compute2D(double theta, double phi, bool world_coord = false);

  std::vector<double> Compute3D(double x, double y, bool world_coord = false, double x0 = 0.1);

  double ComputeError(std::vector<std::vector<double>> v1, std::vector<std::vector<double>> v2);

  double f();

  inline double fx() { return fx_; }
  inline double fy() { return fy_; }
  inline double cx() { return cx_; }
  inline double cy() { return cy_; }
  inline double k1() { return k1_; }
  inline double k2() { return k2_; }
  inline double k3() { return k3_; }
  inline double k4() { return k4_; }
  inline double k5() { return k5_; }

private:
  double fx_ = 0.0;
  double fy_ = 0.0;
  double cx_ = 0.0;
  double cy_ = 0.0;
  double k1_ = 0.0;
  double k2_ = 0.0;
  double k3_ = 0.0;
  double k4_ = 0.0;
  double k5_ = 0.0;
};



class Visualizer {

public: 
  Visualizer(std::string window_name = "3D", double scale = 1.0);

  void AddCloud(std::vector<cv::Point3f> cloud, std::vector<cv::Vec3b> color);

  void AddCloud(std::vector<cv::Point3f> cloud, cv::Vec3b color = cv::Vec3b(100,100,100));

  void Render();

private:

  std::vector<std::vector<cv::Point3f>> point_clouds_;
  std::vector<std::vector<cv::Vec3b>> colors_;

  cv::viz::Viz3d window_;
  std::string window_name_ = "3D";
  double scale_ = 1.0;

};

#endif

