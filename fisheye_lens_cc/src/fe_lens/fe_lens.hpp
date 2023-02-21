#ifndef FE_LENS_HPP
#define FE_LENS_HPP

#include <iostream>
#include <vector>

#include <opencv2/opencv.hpp>
#include <opencv2/viz.hpp> // include the viz module


// Compute focal lenght from camera intrinsics
double FocalLength(const cv::Mat &K);

// Compute the 3D position of a 2D feature, i.e., it's position over the lens surface in camera rf (pixel coords)
cv::Point3f compute3dOverLens(cv::Point2f p, float f, const cv::Mat &D);

// Opposite from the previous one, computes the 2D x,y coordinates from the feature position over the lens surface
cv::Point2f compute2dFromLens(cv::Point3f p, float f, float theta, float phi, const cv::Mat &D);

// Compute the 3D surface of a Kannala-Brandt fisheye lens
std::vector<cv::Point3f> computeFisheyeSurface(std::vector<cv::Point2f> points, const cv::Mat &K, const cv::Mat &D, int width, int height);

// Opposite from the prevoous one, computes the 2D image projection from the feature position over the lens surface
std::vector<cv::Point3f> computeImageFromFisheyeSurface(std::vector<cv::Point3f> points, const cv::Mat &K, const cv::Mat &D, int width, int height);

// Same as before, different approach
std::vector<cv::Point3f> compute3DFrom2D(std::vector<cv::Point2f> points, const cv::Mat &K, const cv::Mat &D, int width, int height);

// Display the 3D surface of a fisheye lens (or any other 3D surface)
void display3dSurface(std::vector<cv::Point3f> points, std::vector<cv::Vec3b> colors = std::vector<cv::Vec3b>());

// Display the 3D surface of a fisheye lens (or any other 3D surface) and the image (also in 3D space)
void display3dSurfaceAndImage(std::vector<std::vector<cv::Point3f>> points, std::vector<cv::Vec3b> colors);


#endif

