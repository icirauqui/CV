#include "fe_lens.hpp"

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
  for (int i = 0; i < points.size(); i++) {
    points3d.push_back(compute3dOverLens(points[i], f, D));
  }

  return points3d;
}



std::vector<cv::Point3f> computeImageFromFisheyeSurface(std::vector<cv::Point3f> points, const cv::Mat &K, const cv::Mat &D, int width, int height) {
  float f = FocalLength(K);


  // Compute the surface of the fisheye lens from the undistorted points
  std::vector<cv::Point3f> points3d;
  for (int i = 0; i < points.size(); i++) {
    // Theta from third component
    float theta = acos(points[i].z / f);

    // Phi from any of the other components
    float phi = acos(points[i].x / (f * sin(theta)));

    float d1 = D.at<double>(0, 0);
    float d2 = D.at<double>(1, 0);
    float d3 = D.at<double>(2, 0);
    float d4 = D.at<double>(3, 0);
    float fx = K.at<double>(0, 0);
    float fy = K.at<double>(1, 1);
    float cx = K.at<double>(0, 2);
    float cy = K.at<double>(1, 2);
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



std::vector<cv::Point3f> compute3DFrom2D(std::vector<cv::Point2f> points, const cv::Mat &K, const cv::Mat &D, int width, int height) {
  float f = FocalLength(K);


  // Compute the surface of the fisheye lens from the undistorted points
  std::vector<cv::Point3f> points3d;
  for (int i = 0; i < points.size(); i++) {

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

  window.spin();
}



void display3dSurfaceAndImage(std::vector<std::vector<cv::Point3f>> points, std::vector<cv::Vec3b> colors) {
  cv::viz::Viz3d window("3D Surface");
  
  for (unsigned int i=0; i<points.size(); i++) {
    std::string name = "cloud" + std::to_string(i);
    cv::viz::WCloud cloud(points[i], colors);
    window.showWidget(name, cloud);
  }
  
  window.spin();
}
