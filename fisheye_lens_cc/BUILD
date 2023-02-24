load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_library")



cc_binary(
    name = "main",
    srcs = [
        #"src/analysis.hpp", 
        "src/main.cpp",
        ],
    deps = [
        "//src/fe_lens:fe_lens",
        "//src/ang_matcher:ang_matcher",
        "@opencv//:opencv",
    ],
    includes = ["src"],
    data = glob(["images/**"]),
    visibility = ["//visibility:public"],
)