# opencv.BUILD file

load("@rules_foreign_cc//foreign_cc:defs.bzl", "cmake")

filegroup(
    name = "srcs",
    srcs = glob(["**"]),
    visibility = ["//visibility:public"],
)
cmake(
    name = "opencv",
    generate_args = ["-GNinja"],
    #additional_inputs=["@opencv_contrib//:modules"],
    build_data=["@opencv_contrib//:modules"],
    cache_entries = {
        "CMAKE_BUILD_TYPE": "RELEASE",
        "INSTALL_C_EXAMPLES": "OFF",
        "INSTALL_PYTHON_EXAMPLES": "OFF",
        "OPENCV_GENERATE_PKGCONFIG": "ON",
        "OPENCV_ENALBLE_NONFREE": "ON",
        "BUILD_EXAMPLES": "OFF",
        #"WITH_VTK": "ON",
        "BUILD_SHARED_LIBS": "ON",
        "BUILD_opencv_world": "ON",
        "OPENCV_EXTRA_MODULES_PATH": "$$EXT_BUILD_ROOT$$/external/opencv_contrib/modules",
        "WITH_QT": "OFF",          
        #"BUILD_LIST": "core,highgui,imgcodecs,imgproc",
    },
    lib_source = ":srcs",
    #lib_source = "@opencv//:all",
    #out_static_libs = ["libopencv_world.a"],
    out_shared_libs = [
        #"libopencv_world.so",
        "libopencv_world.so.407",
        "libopencv_world.so.4.7.0",
        ],
    #out_static_libs = [
    #    "libopencv_alphamat.a",
    #    "libopencv_aruco.a",
    #    "libopencv_barcode.a",
    #    "libopencv_bgsegm.a",
    #    "libopencv_bioinspired.a",
    #    "libopencv_calib3d.a",
    #    "libopencv_ccalib.a",
    #    "libopencv_core.a",
    #    #"libopencv_cvv.a",
    #    "libopencv_datasets.a",
    #    "libopencv_dnn_objdetect.a",
    #    "libopencv_dnn.a",
    #    "libopencv_dnn_superres.a",
    #    "libopencv_dpm.a",
    #    "libopencv_face.a",
    #    "libopencv_features2d.a",
    #    "libopencv_flann.a",
    #    "libopencv_freetype.a",
    #    "libopencv_fuzzy.a",
    #    #"libopencv_gapi.a",
    #    "libopencv_hdf.a",
    #    "libopencv_hfs.a",
    #    "libopencv_highgui.a",
    #    "libopencv_imgcodecs.a",
    #    "libopencv_img_hash.a",
    #    "libopencv_imgproc.a",
    #    "libopencv_intensity_transform.a",
    #    "libopencv_line_descriptor.a",
    #    "libopencv_mcc.a",
    #    "libopencv_ml.a",
    #    "libopencv_objdetect.a",
    #    "libopencv_optflow.a",
    #    "libopencv_phase_unwrapping.a",
    #    "libopencv_photo.a",
    #    "libopencv_plot.a",
    #    "libopencv_quality.a",
    #    "libopencv_rapid.a",
    #    "libopencv_reg.a",
    #    "libopencv_rgbd.a",
    #    "libopencv_saliency.a",
    #    #"libopencv_sfm.a",
    #    "libopencv_shape.a",
    #    "libopencv_stereo.a",
    #    "libopencv_stitching.a",
    #    "libopencv_structured_light.a",
    #    "libopencv_superres.a",
    #    "libopencv_surface_matching.a",
    #    "libopencv_text.a",
    #    "libopencv_tracking.a",
    #    "libopencv_videoio.a",
    #    "libopencv_video.a",
    #    "libopencv_videostab.a",
    #    "libopencv_viz.a",
    #    "libopencv_wechat_qrcode.a",
    #    "libopencv_xfeatures2d.a",
    #    "libopencv_ximgproc.a",
    #    "libopencv_xobjdetect.a",
    #    "libopencv_xphoto.a"
    #],
    #out_shared_libs = [
    #    "libopencv_alphamat.so.4.7.0",
    #    "libopencv_aruco.so.4.7.0",
    #    "libopencv_barcode.so.4.7.0",
    #    "libopencv_bgsegm.so.4.7.0",
    #    "libopencv_bioinspired.so.4.7.0",
    #    "libopencv_calib3d.so.4.7.0",
    #    "libopencv_ccalib.so.4.7.0",
    #    "libopencv_core.so.4.7.0",
    #    #"libopencv_cvv.so.4.7.0",
    #    "libopencv_datasets.so.4.7.0",
    #    "libopencv_dnn_objdetect.so.4.7.0",
    #    "libopencv_dnn.so.4.7.0",
    #    "libopencv_dnn_superres.so.4.7.0",
    #    "libopencv_dpm.so.4.7.0",
    #    "libopencv_face.so.4.7.0",
    #    "libopencv_features2d.so.4.7.0",
    #    "libopencv_flann.so.4.7.0",
    #    "libopencv_freetype.so.4.7.0",
    #    "libopencv_fuzzy.so.4.7.0",
    #    #"libopencv_gapi.so.4.7.0",
    #    "libopencv_hdf.so.4.7.0",
    #    "libopencv_hfs.so.4.7.0",
    #    "libopencv_highgui.so.4.7.0",
    #    "libopencv_imgcodecs.so.4.7.0",
    #    "libopencv_img_hash.so.4.7.0",
    #    "libopencv_imgproc.so.4.7.0",
    #    "libopencv_intensity_transform.so.4.7.0",
    #    "libopencv_line_descriptor.so.4.7.0",
    #    "libopencv_mcc.so.4.7.0",
    #    "libopencv_ml.so.4.7.0",
    #    "libopencv_objdetect.so.4.7.0",
    #    "libopencv_optflow.so.4.7.0",
    #    "libopencv_phase_unwrapping.so.4.7.0",
    #    "libopencv_photo.so.4.7.0",
    #    "libopencv_plot.so.4.7.0",
    #    "libopencv_quality.so.4.7.0",
    #    "libopencv_rapid.so.4.7.0",
    #    "libopencv_reg.so.4.7.0",
    #    "libopencv_rgbd.so.4.7.0",
    #    "libopencv_saliency.so.4.7.0",
    #    #"libopencv_sfm.so.4.7.0",
    #    "libopencv_shape.so.4.7.0",
    #    "libopencv_stereo.so.4.7.0",
    #    "libopencv_stitching.so.4.7.0",
    #    "libopencv_structured_light.so.4.7.0",
    #    "libopencv_superres.so.4.7.0",
    #    "libopencv_surface_matching.so.4.7.0",
    #    "libopencv_text.so.4.7.0",
    #    "libopencv_tracking.so.4.7.0",
    #    "libopencv_videoio.so.4.7.0",
    #    "libopencv_video.so.4.7.0",
    #    "libopencv_videostab.so.4.7.0",
    #    "libopencv_viz.so.4.7.0",
    #    "libopencv_wechat_qrcode.so.4.7.0",
    #    "libopencv_xfeatures2d.so.4.7.0",
    #    "libopencv_ximgproc.so.4.7.0",
    #    "libopencv_xobjdetect.so.4.7.0",
    #    "libopencv_xphoto.so.4.7.0"
    #],
    out_include_dir = "include/opencv4",
    visibility = ["//visibility:public"],
)