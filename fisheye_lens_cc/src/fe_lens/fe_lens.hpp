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
  //double f(double x) { return x * x - 2; }
  //double f_prime(double x) { return 2 * x; }
  //NewtonRaphson solver(1, 1e-6, 100);
  //double x = solver.solve(f, f_prime);

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

  double FocalLength();

  cv::Point2f Compute2D(double theta, double phi, bool sim = false);

  std::vector<double> Compute3D(double x, double y, bool sim = false, double x0 = 0.1);

  double ComputeError(std::vector<std::vector<double>> v1, std::vector<std::vector<double>> v2);

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






// Compute focal lenght from camera intrinsics
double FocalLength(const cv::Mat &K);

// Compute the 3D position of a 2D feature, i.e., it's position over the lens surface in camera rf (pixel coords)
cv::Point3f compute3dOverLens(cv::Point2f p, float f, const cv::Mat &D);

// Compute the 3D position of a 2D feature, i.e., it's position over the lens surface in camera rf (pixel coords)
cv::Point3f compute3dOverLens2(cv::Point2f p, float f, const cv::Mat &D);

// Opposite from the previous one, computes the 2D x,y coordinates from the feature position over the lens surface
cv::Point2f compute2dFromLens(cv::Point3f p, float f, float theta, float phi, const cv::Mat &D);

// Compute the 3D surface of a Kannala-Brandt fisheye lens
std::vector<cv::Point3f> computeFisheyeSurface(std::vector<cv::Point2f> points, const cv::Mat &K, const cv::Mat &D, int width, int height);

// Opposite from the prevoous one, computes the 2D image projection from the feature position over the lens surface. Follows OpenCV docs
std::vector<cv::Point3f> computeImageFromFisheyeSurface(std::vector<cv::Point3f> points, const cv::Mat &K, const cv::Mat &D, int width, int height);

// Opposite from the prevoous one, computes the 2D image projection from the feature position over the lens surface. Follows KB equations
std::vector<cv::Point3f> computeImageFromFisheyeSurface2(std::vector<cv::Point3f> points, const cv::Mat &K, const cv::Mat &D, int width, int height);

// Same as before, different approach
std::vector<cv::Point3f> compute3DFrom2D(std::vector<cv::Point2f> points, const cv::Mat &K, const cv::Mat &D, int width, int height);

// Display the 3D surface of a fisheye lens (or any other 3D surface)
void display3dSurface(std::vector<cv::Point3f> points, std::vector<cv::Vec3b> colors = std::vector<cv::Vec3b>());

// Display the 3D surface of a fisheye lens (or any other 3D surface) and the image (also in 3D space)
void display3dSurfaceAndImage(std::vector<std::vector<cv::Point3f>> points, std::vector<cv::Vec3b> colors);

void display3dSurfaceAndImage(std::vector<std::vector<cv::Point3f>> points, std::vector<std::vector<cv::Vec3b>> colors);


#endif

