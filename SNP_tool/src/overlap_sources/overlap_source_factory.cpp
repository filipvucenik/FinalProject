#include "overlap_source_factory.hpp"
#include <dlfcn.h>
#include <string>
#include <iostream>
#include <filesystem>


OverlapSource * create_overlap_source(const std::string &library, const std::string &args) {
    std::string path = "../lib/lib"+library+"_overlap_source.so";
    void *handle = dlopen(path.c_str(), RTLD_LAZY);
    if (!handle) {
        std::cerr << "Cannot open library: " << dlerror() << std::endl;
        exit(1);
    }
    typedef OverlapSource * (*create_t)(std::string);
    create_t create = (create_t) dlsym(handle, "create");
    if (!create) {
        std::cerr << "Cannot load symbol create: " << dlerror() << std::endl;
        exit(1);
    }
    return create(args);
}