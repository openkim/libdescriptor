#include "../../Descriptors.hpp"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

int main(){
    int n_atoms = 259;
    int width = 51;
    int n_contributing = 64;
    double coords[259 * 3];
    double d_coords[259 * 3];
    double grad[64][3];
    double dzeta[64 * 51];
    double kliffzeta[64][51];
    int * neighbors[64];
    int n_neigh[64];
    int * species;
    auto desc = new double[51 * 64];
    int image[259];

    species = new int[259];
    for (int i = 0; i < 259; i++){
        species[i] = 0;
    }

    fstream image_file("image.txt");
    for (int i = 0; i < n_atoms; i++){
        image_file >> image[i];
    }
    image_file.close();

    fstream coords_file("coords.txt");
    for (int i = 0; i < n_atoms; i++){
        for(int j = 0; j < 3; j++){
             coords_file >> coords[3 * i + j];
             d_coords[3 * i + j] = 0.0;
        }
    }
    coords_file.close();

    fstream grad_file("grad.txt");
    for (int i = 0; i < n_contributing; i++){
        for(int j = 0; j < 3; j++){
             grad_file >> grad[i][j];
        }
    }
    grad_file.close();

    fstream dzeta_file("dedzeta.txt");
    for (int i = 0; i < n_contributing; i++){
        for(int j = 0; j < width; j++){
             dzeta_file >> dzeta[i * 51 + j];
        }
    }
    dzeta_file.close();

    fstream kliff_zeta_file("kliff_desc.dat");
    double tmp;
    for (int i = 0; i < n_contributing; i++){
        for(int j = 0; j < width; j++){
             kliff_zeta_file >> kliffzeta[i][j];
        }
    }
    kliff_zeta_file.close();

    fstream neighbor_file("neighs.txt");
    for (int i = 0; i < n_contributing; i++){
        neighbor_file >> n_neigh[i];
        neighbors[i] = new int [n_neigh[i]];
        for(int j = 0; j < n_neigh[i]; j++){
             neighbor_file >> neighbors[i][j];
        }
    }
    neighbor_file.close();
    
    std::cout << "Loaded data" << std::endl;

    string file_name = "descriptor.dat";
    auto dbs = Descriptor::DescriptorKind::initDescriptor(file_name,Descriptor::KindSymmetryFunctions);
    std::cout << "Initialized descriptor" << std::endl;
    cout << dbs->descriptor_param_file << "\n";
    cout <<dbs->width <<"\n";

    // FWD check
    double aggregated_error = 0;
    //for(int i = 0; i < 64; i++) {
    int i = 0;
        Descriptor::compute_single_atom(i,
                                        259,
                                        species,
                                        neighbors[i],
                                        n_neigh[i],
                                        coords,
                                        desc + i * width,
                                        dbs);
        for (int j = 0; j < 51; j++){
            aggregated_error += abs(kliffzeta[i][j] - *(desc + i * 51 + j));
        }
    //}

    // print first 10 elements of both
    for (int i = 0; i < 10; i++){
        cout << kliffzeta[0][i] << " " << *(desc + i) << "\n";
    }


    if (aggregated_error == 0.){
        cout << "Descriptor agrees exactly.\n";
    } else if(aggregated_error < 1e-6){
        cout << "Descriptor agrees, tol < 1e-6\n";
    } else {
        cout << "Descriptor error > 1e-6 | " << aggregated_error << " |\n";
    }

    std::cout << "computed fwd\n";
//    aggregated_error = 0;
//    //for(int i = 0; i < 64; i++) {
    Descriptor::gradient_single_atom(i,
                                     259,
                                     species,
                                     neighbors[i],
                                     n_neigh[i],
                                     coords,
                                     d_coords,
                                     desc + i * width,
                                     dzeta + i * width,
                                     dbs);
//    //}
//
//       // Sum forces
    double forces[64][3];
    for(int i = 0; i < n_atoms; i++){
        for(int j = 0; j < 3; j++){
            forces[image[i]][j] += *(d_coords + i * 3 + j);
        }
    }
    // print 64 x 3 forces
    for (int i = 0; i < 64; i++){
        for (int j = 0; j < 3; j++){
            cout << forces[i][j] << " " << grad[i][j] ;
        }
        cout << "\n";
    }
//
//
//    for(int i = 0; i < n_contributing; i++){
//        for(int j = 0; j < 3; j++){
//            aggregated_error += abs(grad[i][j] - forces[i][j]);
//        }
//    }
//    cout << "Forces aggregated errors:" << aggregated_error << "\n";

    delete[] desc;
    delete[] species;
    for(int i = 0; i < n_contributing; i++){
        delete[] neighbors[i];
    }
}
