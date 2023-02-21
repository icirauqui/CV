#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/viz.hpp> // include the viz module


#include "fe_lens/fe_lens.hpp"





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
  display3dSurfaceAndImage(point_clouds, colors);

  //display3dSurface(points3d, colors);

  return 0;
}
