// Generic header for inclusion and library design
#include "SymFun/SymFun.hpp"
#include <map>
#include <string>

class Descriptor {
    // Base class for all descriptors
    // Figure out how to make it work with enzyme
    // leaving blank for now
public:
    std::string descriptor_kind;
    std::string descriptor_param_file;
    std::map<std::string, void *> descriptor_map;

    Descriptor();
    Descriptor(std::string& descriptor_name);
    Descriptor(std::string& descriptor_name, std::string& descriptor_params);
//    void getDescriptor(std::string,double *);
    void initDescriptor();
    ~Descriptor();

private:
    std::map<std::string, void *> descriptor_data_map;
    std::map<std::string, void *> descriptor_grad_map;
    double bhor2ang = 0.529177;
};

