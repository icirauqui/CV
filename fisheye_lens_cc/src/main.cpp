#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/viz.hpp> // include the viz module

#include <fstream>


#include "fe_lens/fe_lens.hpp"


#include "src/ang_matcher/ang_matcher.h"


void FullLens() {
    // Set the camera intrinsics
  FisheyeLens lens(717.2104, 717.4816, 735.3566, 552.7982, 
                   -0.1389272, -0.001239606, 0.0009125824, -0.00004071615);

  // Define max angles for semi-sphere
  double theta_max = 60 * M_PI / 180;  // 60 degrees
  double phi_max = 2 * M_PI;  // 360 degrees
  std::cout << "Theta max: " << theta_max << std::endl;
  std::cout << "Phi max: " << phi_max << std::endl;

  // Generate 3D points over the semi-sphere
  std::vector<std::vector<double>> coords3d;
  coords3d.push_back(std::vector<double>({0.0, 0.0}));  // Center
  double resolution = 0.05;
  for (double theta = resolution; theta < theta_max; theta += resolution) {
    for (double phi = 0.0; phi < phi_max; phi += resolution) {
      coords3d.push_back(std::vector<double>({theta, phi}));
    }
  }

  // Project the semi-sphere onto the image plane
  std::vector<cv::Point2f> coords2d;
  //for (auto coord : coords3d) {
  for (unsigned int i=0; i<coords3d.size(); i++) {
    double theta = coords3d[i][0];
    double phi = coords3d[i][1];
    cv::Point2f coord = lens.Compute2D(theta, phi, true);
    coords2d.push_back(coord);
  }

  // Project back the image points onto the semi-sphere
  std::vector<std::vector<double>> coords3d_reconstr;
  for (auto coord : coords2d) {
    double x = coord.x;
    double y = coord.y;
    coords3d_reconstr.push_back(lens.Compute3D(x, y, true));
    //if (it == 3)
    //  break;
  }

  // Compute the error 
  double error = lens.ComputeError(coords3d, coords3d_reconstr);
  std::cout << "Global error = " << error << std::endl;


  // - - - - Visualization - - - - - - - - - - - - - - - - - -

  float scale = 1;
  float offset = 0.75;

  // Original Lens
  std::vector<cv::Point3f> points_lens;
  std::vector<cv::Vec3b> colors_lens;
  for (auto coord: coords3d) {
    double x = scale * sin(coord[0]) * cos(coord[1]);
    double y = scale * sin(coord[0]) * sin(coord[1]);
    double z = offset + scale * cos(coord[0]);
    points_lens.push_back(cv::Point3f(x, y, z));
    colors_lens.push_back(cv::Vec3b(255, 255, 255));
  }
  
  // Projected image
  std::vector<cv::Point3f> points_image;
  std::vector<cv::Vec3b> colors_image;
  for (auto coord: coords2d) {
    double x = scale * coord.x;
    double y = scale * coord.y;
    double z = -offset;
    points_image.push_back(cv::Point3f(x, y, z));
    colors_image.push_back(cv::Vec3b(0, 0, 255));
  }  

  // Reconstructed lens
  std::vector<cv::Point3f> points_lens_reconstr;
  std::vector<cv::Vec3b> colors_lens_reconstr;
  for (auto coord: coords3d_reconstr) {
    double x = scale * sin(coord[0]) * cos(coord[1]);
    double y = scale * sin(coord[0]) * sin(coord[1]);
    double z = 2*offset + scale * cos(coord[0]);
    points_lens_reconstr.push_back(cv::Point3f(x, y, z));
    colors_lens_reconstr.push_back(cv::Vec3b(0, 255, 0));
  }


  Visualizer vis;
  vis.AddCloud(points_lens, colors_lens);
  vis.AddCloud(points_image, colors_image);
  vis.AddCloud(points_lens_reconstr, colors_lens_reconstr);
  vis.Render();
}



void ImgMethod() {
    // Set the camera intrinsics
  FisheyeLens lens(717.2104, 717.4816, 735.3566, 552.7982, 
                   -0.1389272, -0.001239606, 0.0009125824, -0.00004071615);
  
  cv::Mat im1 = imread("images/1.png", cv::IMREAD_COLOR);
  std::vector<cv::Point2f> points;
  std::vector<cv::Vec3b> colors_original;
  std::vector<cv::Point3f> image3d;
  std::cout << "Image size: " << im1.cols << "x" << im1.rows << std::endl;
  int resolution = 1;
  for (int x = 0; x < im1.cols; x+=resolution) {
    for (int y = 0; y < im1.rows; y+=resolution) {
      double xd = (x - lens.cx()) / lens.fx();
      double yd = (y - lens.cy()) / lens.fy();
      points.push_back(cv::Point2f(x, y));
      colors_original.push_back(im1.at<cv::Vec3b>(y, x));
      image3d.push_back(cv::Point3f(-xd, -yd, -lens.f()));
    }
  }

  // Project back the image points onto the semi-sphere
  std::vector<std::vector<double>> coords3d_reconstr;
  //double f = 1.0; //lens.f();
  for (auto coord : points) {
    double x = coord.x;
    double y = coord.y;
    std::vector<double> coord3d = lens.Compute3D(x, y, false, 1.0);
    coords3d_reconstr.push_back(coord3d);
  }

  // Save points to csv
  std::ofstream myfile;
  myfile.open("/home/icirauqui/w0rkspace/CV/fisheye_lens_cc/points.csv");
  for (auto coord : coords3d_reconstr) {
    myfile << coord[0] << "," << coord[1] << "," << coord[2] << std::endl;
  }
  myfile.close();




  // - - - - Visualization - - - - - - - - - - - - - - - - - -
 
  float scale = 1;
  float offset = 0.75;

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
  }

  Visualizer vis;
  vis.AddCloud(points_image, colors_original);
  vis.AddCloud(points_lens_reconstr, colors_original);
  vis.Render();
}







void ImgMatching() {
  FisheyeLens lens(717.2104, 717.4816, 735.3566, 552.7982, 
                   -0.1389272, -0.001239606, 0.0009125824, -0.00004071615);

  // Load images  
  std::cout << " 1. Loading images" << std::endl;
  std::vector<Image> imgs;
  imgs.push_back(Image(imread("images/1.png", cv::IMREAD_COLOR), &lens));
  imgs.push_back(Image(imread("images/2.png", cv::IMREAD_COLOR), &lens));



  std::cout << " 2. Detecting features" << std::endl;

  // Detect features (parameters from COLMAP)
  int max_features = 1000; //8192;
  int num_octaves = 4;
  int octave_resolution = 3;
  float peak_threshold = 0.02 / octave_resolution;  // 0.04
  float edge_threshold = 10;
  float sigma = 1.6;

  cv::Ptr<cv::SIFT> f2d = cv::SIFT::create(max_features, num_octaves, peak_threshold, edge_threshold, sigma);
  for (auto &img : imgs) {
    f2d->detect(img.image_, img.kps_, cv::noArray());
    img.contours_ = am::GetContours(img.image_, 20, 3, false);
    img.kps_ = am::KeypointsInContour(img.contours_[0], img.kps_);
    f2d->compute(img.image_, img.kps_, img.desc_);
    std::cout << "      Keypoints: " << img.kps_.size() 
              << "      Descriptors: " << img.desc_.rows << std::endl;
  }
  



  std::cout << " 3. Matching features" << std::endl;

  //float max_ratio = 0.8f;           //COLMAP
  //float max_distance = 0.7f;        //COLMAP
  //int max_num_matches = 32768;      //COLMAP
  //float max_error = 4.0f;           //COLMAP
  //float confidence = 0.999f;        //COLMAP
  //int max_num_trials = 10000;       //COLMAP
  //float min_inliner_ratio = 0.25f;  //COLMAP
  //int min_num_inliers = 15;         //COLMAP

  std::vector<cv::DMatch> matches = am::MatchFLANN(imgs[0].desc_, imgs[1].desc_, 0.7f);
  std::cout << " 3.1. Matches: " << matches.size() << std::endl;



  std::cout << " 4. Compute F and epilines" << std::endl;

  // Compute F and epilines
  std::vector<cv::Point2f> points1, points2;
  for (unsigned int i = 0; i < matches.size(); i++) {
    points1.push_back(imgs[0].kps_[matches[i].queryIdx].pt);
    points2.push_back(imgs[1].kps_[matches[i].trainIdx].pt);
  }
  cv::Mat F12 = cv::findFundamentalMat(points1, points2);

  std::cout << " 4.1 Decompose E" << std::endl;
  cv::Mat Kp = lens.K();
  Kp.convertTo(Kp, CV_64F);
  cv::Mat E = am::EfromF(F12, Kp);

  cv::Mat R1, R2, t;
  cv::decomposeEssentialMat(E, R1, R2, t);

  std::cout << " 4.2. R1" << std::endl;
  std::cout << R1 << std::endl;
  std::cout << " 4.3. R2" << std::endl;
  std::cout << R2 << std::endl;
  std::cout << " 4.4. t" << std::endl;
  std::cout << t << std::endl;


  double tx = t.at<double>(0, 0);
  double ty = t.at<double>(0, 1);
  double tz = t.at<double>(0, 2);

  
  







  std::cout << " 5. Compute matches by distance and angle" << std::endl;
    
  bool cross_check = true;
  bool draw_inline = false;
  bool draw_global = false;
  
  //float th_epiline = 4.0;
  float th_sampson = 4.0;
  //float th_angle2d = DegToRad(1.0);
  float th_angle3d = DegToRad(4.0);
  double th_sift = 100.0;

  cv::Point3f c1g(0.0, 0.0, 0.0);
  cv::Point3f c2g = c1g + cv::Point3f(tx, ty, tz);

  cv::Point3f cc = lens.c();


  // Match by distance threshold
  am::AngMatcher angmatcher(imgs[0].kps_, imgs[1].kps_, imgs[0].desc_, imgs[1].desc_,
                            F12, imgs[0].image_, imgs[1].image_, 2*lens.cx(), 2*lens.cy(),
                            lens.f(), cc, cc, c1g, c2g, R1, R2, t, Kp, lens.D());

  //angmatcher.Match("sampson", th_sampson, th_sift, cross_check, draw_inline, draw_global);
  angmatcher.Match("angle3d", th_angle3d, th_sift, cross_check, draw_inline, draw_global);

  //am.ViewMatches("epiline", "epiline desc matches", 0.5);



  std::cout << " 6. Compare matches" << std::endl;
  int report_level = 1;
  //am.CompareMatches("epiline", "sampson", report_level);
  //am.CompareMatches("epiline", "angle2d", report_level);
  //am.CompareMatches("epiline", "angle3d", report_level);
  //am.CompareMatches("sampson", "angle2d", report_level);
  angmatcher.CompareMatches("sampson", "angle3d", report_level);
  //am.CompareMatches("angle2d", "angle3d", report_level);


  std::cout << " 7. Compare matches for specific query keypoint" << std::endl;
  //angmatcher.ViewCandidatesCompare("sampson", "angle3d",  32);
  //angmatcher.ViewCandidatesCompare("sampson", "angle3d", 395);
  //angmatcher.ViewCandidatesCompare("sampson", "angle3d", 409);
  //angmatcher.ViewCandidatesCompare("sampson", "angle3d", 430);
  //angmatcher.ViewCandidatesCompare("sampson", "angle3d", 473);
  //angmatcher.ViewCandidatesCompare("sampson", "angle3d", 642);
  //angmatcher.ViewCandidatesCompare("sampson", "angle3d", 644);



  // For kp = 430 in img 0, draw epipolar plane
  cv::Point2f kp = imgs[0].kps_[430].pt;
  std::vector<double> kp3d_v = lens.Compute3D(kp.x, kp.y, true);



  



  //cv::waitKey(0);






  // - - - - Visualization - - - - - - - - - - - - - - - - - -
 
  float scale = 1;
  cv::Vec3d offset(0, 0, 0);

  std::vector<std::vector<cv::Point3f>> points_images;
  for (unsigned int i=0; i<imgs.size(); i++) {
    // Projected image
    std::vector<cv::Point3f> points_image;
    for (auto coord: imgs[i].image3d_) {
      double x = (i*tx) + (i*offset(0)) + scale * coord.x;
      double y = (i*ty) + (i*offset(1)) + scale * coord.y;
      double z = (i*tz) + (i*offset(2)) + 0.0;
      points_image.push_back(cv::Point3f(x, y, z));
    }  
    points_images.push_back(points_image);
  }



  // Reconstructed lens
  std::vector<std::vector<cv::Point3f>> points_lens;

  for (unsigned int i=0; i<imgs.size(); i++) {

    std::vector<cv::Point3f> points_lens_reconstr;
    for (auto coord: imgs[i].lens_->Lens3dReconstr()) {
      double x = (i*tx) + (i*offset(0)) + scale * sin(coord[0]) * cos(coord[1]);
      double y = (i*ty) + (i*offset(1)) + scale * sin(coord[0]) * sin(coord[1]);
      double z = (i*tz) + (i*offset(2)) + scale * cos(coord[0]);
      if (z < 0) {
        x = 0.0;
        y = 0.0;
        z = 0.0;
      }
      points_lens_reconstr.push_back(cv::Point3f(x, y, z));
    }
    points_lens.push_back(points_lens_reconstr);
  }



  // Draw a plane in 3D that goes through c1g, c2g and kp3d
  //double kp3d_theta = kp3d_v[0];
  //double kp3d_phi = kp3d_v[1];
  //double kp3d_x = sin(kp3d_theta) * cos(kp3d_phi);
  //double kp3d_y = sin(kp3d_theta) * sin(kp3d_phi);
  //double kp3d_z = cos(kp3d_theta);
  //cv::Point3f kp3d(kp3d_x, kp3d_y, kp3d_z);
  //cv::Point3f n = (c2g - c1g).cross(kp3d - c1g);
  //n = n / cv::norm(n);

  

  




  Visualizer vis;
  for (unsigned int i=0; i<points_images.size(); i++) {
    vis.AddCloud(points_images[i], imgs[i].colors_);
    vis.AddCloud(points_lens[i], imgs[i].colors_);
  }

  vis.AddPlane(kp3d_v, c1g, c2g, 5);

  // Draw epipolar plane
  //cv::viz::WPlane plane_widget(cv::Point3d(0, 0, 0), 
  //                             cv::Vec3d(n.x, n.y, n.z), 
  //                             cv::Vec3d(0, 0, 1),
  //                             cv::Size2d(10, 10));
  //plane_widget.setRenderingProperty(cv::viz::LINE_WIDTH, 2.0);
  //plane_widget.setRenderingProperty(cv::viz::OPACITY, 0.5);
  //vis.AddWidget(plane_widget);



  vis.Render();




/*


  // - - - - Visualization - - - - - - - - - - - - - - - - - -
 
  float scale = 1;
  cv::Vec3d offset(0, 0, 0);

  std::vector<std::vector<cv::Point3f>> points_images;
  for (unsigned int i=0; i<images.size(); i++) {
    // Projected image
    std::vector<cv::Point3f> points_image;
    for (auto coord: image3d) {
      double x = (i*tx) + (i*offset(0)) + scale * coord.x;
      double y = (i*ty) + (i*offset(1)) + scale * coord.y;
      double z = (i*tz) + (i*offset(2)) + 0.0;
      points_image.push_back(cv::Point3f(x, y, z));
    }  
    points_images.push_back(points_image);
  }



  // Reconstructed lens
  std::vector<cv::Point3f> points_lens_reconstr;
  for (auto coord: images_3d[0]) {
    double x = scale * sin(coord[0]) * cos(coord[1]);
    double y = scale * sin(coord[0]) * sin(coord[1]);
    double z = scale * cos(coord[0]);
    if (z < 0) {
      x = 0.0;
      y = 0.0;
      z = 0.0;
    }
    points_lens_reconstr.push_back(cv::Point3f(x, y, z));
  }

  Visualizer vis;
  vis.AddCloud(points_images[0], colors_original);
  vis.AddCloud(points_images[1], colors_original);
  vis.AddCloud(points_lens_reconstr, colors_original);
  vis.Render();

  */
}



int main() {

  //FullLens();

  //ImgMethod();

  ImgMatching();

  return 0;
}
