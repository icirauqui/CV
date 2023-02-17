#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/viz.hpp> // include the viz module


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


void display3dSurface(std::vector<cv::Point3f> points, std::vector<cv::Vec3b> colors = std::vector<cv::Vec3b>()) {
  cv::viz::Viz3d window("3D Surface");

  if (colors.size() > 0) {
    cv::viz::WCloud cloud(points, colors);
    window.showWidget("cloud", cloud);
  } else {
    cv::viz::WCloud cloud(points, cv::viz::Color::white());
    window.showWidget("cloud", cloud);
  }

  //cv::viz::WCloud cloud(points, cv::viz::Color::white());
  //window.showWidget("cloud", cloud);
  window.spin();
}


void display3dSurfaceAndImage(std::vector<std::vector<cv::Point3f>> points, std::vector<cv::Vec3b> colors) {
  cv::viz::Viz3d window("3D Surface");
  
  for (unsigned int i=0; i<points.size(); i++) {
    std::string name = "cloud" + std::to_string(i);
    cv::viz::WCloud cloud(points[i], colors);
    window.showWidget(name, cloud);
  }
  //cv::viz::WCloud cloud_surface(points_surface, colors);
  //cv::viz::WCloud cloud_image(points_image, colors);

  //window.showWidget("points_surface", cloud_surface);
  //window.showWidget("points_image", cloud_image);
  
  window.spin();
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


    //cv::Point2f p = compute2dFromLens(points[i], f, theta, phi, D);

    //points3d.push_back(cv::Point3f(-p.x, -p.y, -2*f));  // Clear the 2*
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




int main()
{
  // Set the camera intrinsics
  cv::Mat K = cv::Mat::eye(3, 3, CV_64F);
  K.at<double>(0, 0) = 717.2104;
  K.at<double>(1, 1) = 717.4816;
  K.at<double>(0, 2) = 735.3566;
  K.at<double>(1, 2) = 552.7982;
  cv::Mat D = cv::Mat::zeros(4, 1, CV_64F);
  D.at<double>(0, 0) = -0.1389272;
  D.at<double>(1, 0) = -0.001239606;
  D.at<double>(2, 0) = 0.0009125824;
  D.at<double>(3, 0) = -0.00004071615;
  float f = FocalLength(K);

  int im_width = 1440;
  int im_height = 1080;

  float cx = K.at<double>(0, 2);
  float cy = K.at<double>(1, 2);

  cv::Mat im1 = imread("images/1.png", cv::IMREAD_COLOR);
  std::vector<cv::Point2f> points;
  std::vector<cv::Vec3b> colors;
  std::vector<cv::Point3f> image3d;
  for (int x = 0; x < im1.cols; x++) {
    for (int y = 0; y < im1.rows; y++) {
      if (im1.at<cv::Vec3b>(y, x)[0] != 0 || im1.at<cv::Vec3b>(y, x)[1] != 0 || im1.at<cv::Vec3b>(y, x)[2] != 0) {
        points.push_back(cv::Point2f(x - cx, y - cy));
        colors.push_back(im1.at<cv::Vec3b>(y, x));
        image3d.push_back(cv::Point3f(-(x - cx), -(y - cy), -f));
      }
    }
  }


  std::vector<cv::Point3f> points3d = computeFisheyeSurface(points, K, D, im_width, im_height);

  std::vector<cv::Point3f> image_reconstr = computeImageFromFisheyeSurface(points3d, K, D, im_width, im_height);

  std::vector<cv::Point3f> points3d_surface = compute3DFrom2D(points, K, D, im_width, im_height);

  std::vector<std::vector<cv::Point3f>> point_clouds;
  point_clouds.push_back(points3d);
  point_clouds.push_back(image3d);
  //point_clouds.push_back(image_reconstr);
  //point_clouds.push_back(points3d_surface);
  display3dSurfaceAndImage(point_clouds, colors);

  //display3dSurface(points3d, colors);

  return 0;
}
