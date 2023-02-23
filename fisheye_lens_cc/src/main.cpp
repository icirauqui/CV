#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/viz.hpp> // include the viz module

#include <fstream>


#include "fe_lens/fe_lens.hpp"


void BuildLens() {
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

  std::vector<cv::Point3f> image_reconstr_2 = computeImageFromFisheyeSurface2(points3d, K, D, im_width, im_height);

  std::vector<cv::Point3f> points3d_surface = compute3DFrom2D(points, K, D, im_width, im_height);

  std::vector<std::vector<cv::Point3f>> point_clouds;
  point_clouds.push_back(points3d);
  //point_clouds.push_back(image3d);
  point_clouds.push_back(image_reconstr_2);
  display3dSurfaceAndImage(point_clouds, colors);
}

void FullLens() {
    // Set the camera intrinsics
  FisheyeLens lens(717.2104, 717.4816, 735.3566, 552.7982, -0.1389272, -0.001239606, 0.0009125824, -0.00004071615);
  cv::Mat im = imread("images/1.png", cv::IMREAD_COLOR);

  // Define max angles for semi-sphere
  double theta_max = 60*M_PI/180;  // 60 degrees
  double phi_max = 2 * M_PI;  // 360 degrees
  std::cout << "Theta max: " << theta_max << std::endl;
  std::cout << "Phi max: " << phi_max << std::endl;

  // Generate 3D points over the semi-sphere
  std::vector<std::vector<double>> coords3d;
  coords3d.push_back(std::vector<double>({0.0, 0.0}));  // Center
  double resolution = 0.1;
  for (double theta = resolution; theta < theta_max; theta += resolution) {
    for (double phi = 0.0; phi < phi_max; phi += resolution) {
      coords3d.push_back(std::vector<double>({theta, phi}));
    }
  }

  // Project the semi-sphere onto the image plane
  std::vector<cv::Point2f> coords2d;
  std::vector<std::vector<double>> coords3d_reconstr;
  //for (auto coord : coords3d) {
  for (unsigned int i=0; i<coords3d.size(); i++) {
    double theta = coords3d[i][0];
    double phi = coords3d[i][1];
    cv::Point2f coord = lens.Compute2D(theta, phi);
    //std::cout << "theta/phi: " << theta << " " << phi << " -> " << coord.x << " " << coord.y << std::endl;
    coords2d.push_back(coord);
  }

  // Project back the image points onto the semi-sphere
  int it = 0;
  for (auto coord : coords2d) {
    it++;
    double x = coord.x;
    double y = coord.y;
    coords3d_reconstr.push_back(lens.Compute3D(x, y));
    //if (it == 3)
    //  break;
  }

  // Compute the error 
  double error = lens.ComputeError(coords3d, coords3d_reconstr);
  std::cout << "Global error = " << error << std::endl;

  float scale = 1;

  // - - - - Visualization - - - - - - - - - - - - - - - - - -

  // Original Lens
  std::vector<cv::Point3f> points_lens;
  std::vector<cv::Vec3b> colors_lens;
  for (auto coord: coords3d) {
    double x = scale * sin(coord[0]) * cos(coord[1]);
    double y = scale * sin(coord[0]) * sin(coord[1]);
    double z = scale * cos(coord[0]);
    points_lens.push_back(cv::Point3f(x, y, z));
    colors_lens.push_back(cv::Vec3b(255, 150, 150));
  }
  
  // Projected image
  std::vector<cv::Point3f> points_image;
  std::vector<cv::Vec3b> colors_image;
  for (auto coord: coords2d) {
    double x = scale * coord.x;
    double y = scale * coord.y;
    double z = 0;
    points_image.push_back(cv::Point3f(x, y, z));
    colors_image.push_back(cv::Vec3b(50, 50, 150));
  }  

  // Reconstructed lens
  std::vector<cv::Point3f> points_lens_reconstr;
  std::vector<cv::Vec3b> colors_lens_reconstr;
  for (auto coord: coords3d_reconstr) {
    double x = scale * sin(coord[0]) * cos(coord[1]);
    double y = scale * sin(coord[0]) * sin(coord[1]);
    double z = scale * cos(coord[0]);
    points_lens_reconstr.push_back(cv::Point3f(x, y, z));
    colors_lens_reconstr.push_back(cv::Vec3b(100, 250, 100));
  }

  // Aggregate results and visualize
  std::vector<std::vector<cv::Point3f>> point_clouds;
  std::vector<std::vector<cv::Vec3b>> colors;
  point_clouds.push_back(points_lens);
  point_clouds.push_back(points_image);
  point_clouds.push_back(points_lens_reconstr);
  colors.push_back(colors_lens);
  colors.push_back(colors_image);
  colors.push_back(colors_lens_reconstr);
  display3dSurfaceAndImage(point_clouds, colors);
}



void ImgMethod2() {
    // Set the camera intrinsics
  FisheyeLens lens(717.2104, 717.4816, 735.3566, 552.7982, -0.1389272, -0.001239606, 0.0009125824, -0.00004071615);
  
  cv::Mat im1 = imread("images/1.png", cv::IMREAD_COLOR);
  std::vector<cv::Point2f> points;
  std::vector<cv::Vec3b> colors_original;
  std::vector<cv::Point3f> image3d;
  std::cout << "Image size: " << im1.cols << "x" << im1.rows << std::endl;
  int resolution = 1;
  for (int x = 0; x < im1.cols; x+=resolution) {
    for (int y = 0; y < im1.rows; y+=resolution) {
      //if (im1.at<cv::Vec3b>(y, x)[0] != 0 || im1.at<cv::Vec3b>(y, x)[1] != 0 || im1.at<cv::Vec3b>(y, x)[2] != 0) {
        double xd = (x - lens.cx()) / lens.fx();
        double yd = (y - lens.cy()) / lens.fy();
        //double xd = x;
        //double yd = y;
        points.push_back(cv::Point2f(x, y));
        colors_original.push_back(im1.at<cv::Vec3b>(y, x));
        image3d.push_back(cv::Point3f(-xd, -yd, -lens.FocalLength()));
      //}
    }
  }

  // Project back the image points onto the semi-sphere
  std::vector<std::vector<double>> coords3d_reconstr;
  double f = 1.0; //lens.FocalLength();
  int it = 0;
  for (auto coord : points) {
    //std::cout << it++ << "/" << points.size() << std::endl;
    it++;
    double x = coord.x;
    double y = coord.y;
    //std::cout << x << " " << y << std::endl;
    std::vector<double> coord3d = lens.Compute3D(x, y, false, 1.0);
    //for (auto &c : coord3d) { c*=f; }
    coords3d_reconstr.push_back(coord3d);
    //if (it == 3) 
    //  break;
  }

  // Save points to csv
  std::ofstream myfile;
  myfile.open("/home/icirauqui/w0rkspace/CV/fisheye_lens_cc/points.csv");
  for (auto coord : coords3d_reconstr) {
    myfile << coord[0] << "," << coord[1] << "," << coord[2] << std::endl;
  }
  myfile.close();



  float scale = 1;
  float offset = 0.5;

  // - - - - Visualization - - - - - - - - - - - - - - - - - -
 
  // Projected image
  std::vector<cv::Point3f> points_image;
  std::vector<cv::Vec3b> colors_image;
  for (auto coord: image3d) {
    double x = scale * coord.x;
    double y = scale * coord.y;
    double z = - offset;
    points_image.push_back(cv::Point3f(x, y, z));
    colors_image.push_back(cv::Vec3b(50, 50, 150));
  }  

  // Reconstructed lens
  std::vector<cv::Point3f> points_lens_reconstr;
  std::vector<cv::Vec3b> colors_lens_reconstr;
  for (auto coord: coords3d_reconstr) {
    double x = scale * sin(coord[0]) * cos(coord[1]);
    double y = scale * sin(coord[0]) * sin(coord[1]);
    double z = offset + scale * cos(coord[0]);
    if (z - offset < 0) {
      x = 0.0;
      y = 0.0;
      z = 0.0;
    }
    points_lens_reconstr.push_back(cv::Point3f(x, y, z));
    colors_lens_reconstr.push_back(cv::Vec3b(100, 250, 100));
  }

  // Aggregate results and visualize
  std::vector<std::vector<cv::Point3f>> point_clouds;
  std::vector<std::vector<cv::Vec3b>> colors;
  point_clouds.push_back(points_image);
  point_clouds.push_back(points_lens_reconstr);
  colors.push_back(colors_original);
  colors.push_back(colors_original);
  display3dSurfaceAndImage(point_clouds, colors);
}

int main() {

  //FullLens();

  ImgMethod2();

  //display3dSurface(points3d, colors);

  return 0;
}
