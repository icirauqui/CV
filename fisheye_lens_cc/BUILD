load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_library")



cc_binary(
    name = "main",
    srcs = ["src/main.cpp"],
    deps = [
        "//src/fe_lens:fe_lens",
        "@opencv//:opencv",
    ],
    includes = ["src"],
    data = glob(["images/**"]),
    visibility = ["//visibility:public"],
)