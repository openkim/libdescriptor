// Generic header for inclusion and library design

#ifndef __DESCRIPTOR_HPP_
#define __DESCRIPTOR_HPP_

#include <map>
#include <string>
#include <vector>
#include <stdexcept>

/*!
 * \brief Namespace containing all required classes and methods for generating descriptors, as well as their derivatives.
 *
 * The namespace Descriptor contains wrappers to all of the implemented descriptors.
 * The wrapper functions are used, instead of more direct inheritance, as Enzyme as of now (v0.0.41)
 * cannot accept pointer to class member function. This design might get more streamlined, as more mature version of Enzyme
 * is released.
 */
namespace Descriptor {
    /// This enumerator lists all the kind of descriptors that are available at run time.
    enum AvailableDescriptor {
        KindSymmetryFunctions, //!< For selecting Behler Symmetry Functions (SymmetryFunctions)
        KindBispectrum //!< For selecting Bispectrum descriptor
    };

    class DescriptorKind;

    /*!
     * <b> "TO BE IMPLEMENTED"</b>
     * This method will use Forward differentiation method (Jacobian-vector product) for computing derivatives.
     * It is useful when descriptor are large.
     * @param n_atoms number of total contributing atoms
     * @param species species code of <i>all</i> atoms
     * @param neighbor_list list of neighbor lists
     * @param number_of_neighs list of number of neighbors
     * @param coordinates array containing coordinated of all atoms
     * @param desc pointer for storing computed descriptor
     * @param descriptor_kind Initialized descriptor kind class
     */
    inline void fwd_diff(int n_atoms, int *species, int *neighbor_list, int *number_of_neighs,
                         double *coordinates, double *desc,
                         DescriptorKind *descriptor_kind) { throw std::logic_error("Function not implemented yet."); }

    /*!
     *  <b> "TO BE IMPLEMENTED"</b>
     *  This method will compute gradient for all contributing atoms using reverse diff (vector-Jacobian product)
     * @param n_atoms
     * @param species
     * @param neighbor_list
     * @param number_of_neighs
     * @param coordinates
     * @param d_coordinates
     * @param desc
     * @param dE_dzeta
     * @param descriptor_to_diff
     */
    void gradient(int n_atoms, int *species, int *neighbor_list, int *number_of_neighs,
                  double *coordinates, double *d_coordinates, double *desc,
                  double *dE_dzeta, DescriptorKind *descriptor_to_diff);

    /*!
     * This method computes the gradient of single atoms. It is somewhat inefficient as it needs to create Enzyme
     * virtual object for every atom. In future it will be used for more selective problems, and gradient() would
     * be the main workhorse.
     * @param index
     * @param n_atoms
     * @param species
     * @param neighbor_list
     * @param number_of_neighs
     * @param coordinates
     * @param d_coordinates
     * @param desc
     * @param dE_dzeta
     * @param descriptor_to_diff
     */
    void gradient_single_atom(int index, int n_atoms, int *species, int *neighbor_list, int number_of_neighs,
                              double *coordinates, double *d_coordinates, double *desc,
                              double *dE_dzeta, DescriptorKind *descriptor_to_diff);

    void num_gradient_single_atom(int index, int n_atoms, int *species, int *neighbor_list, int number_of_neighs,
                              double *coordinates, double *d_coordinates, double *dE_dzeta,
                              DescriptorKind *descriptor_to_diff);

    /*!
     *  <b> "TO BE IMPLEMENTED"</b>
     *  Compute full Jacobian using reverse diff. Might be useful in some cases
     * @param n_atoms
     * @param species
     * @param neighbor_list
     * @param number_of_neighs
     * @param coordinates
     * @param d_coordinates
     * @param desc
     * @param dzeta_dr
     * @param descriptor_to_diff
     */
    inline void jacobian(int n_atoms, int *species, int *neighbor_list, int *number_of_neighs,
                         double *coordinates, double *d_coordinates, double *desc,
                         double *dzeta_dr, DescriptorKind *descriptor_to_diff) {
        throw std::logic_error("Function not implemented yet");
    }

    /*!
     * Computes the descriptor for all contributing atoms, and stores them in pointer desc
     * @param n_atoms
     * @param species
     * @param neighbor_list
     * @param number_of_neighs
     * @param coordinates
     * @param desc
     * @param descriptor_kind
     */
    void compute(int n_atoms, int *species, int *neighbor_list, int *number_of_neighs,
                 double *coordinates, double *desc,
                 DescriptorKind *descriptor_kind);

    /*!
     * Computes the descriptors for single atom (at position @param index) and stores its descriptor in @param desc
     * @param index
     * @param n_atoms
     * @param species
     * @param number_of_neighs
     * @param number_of_neigh_list
     * @param coordinates
     * @param desc
     * @param descriptor_kind
     */
    void compute_single_atom(int index, int n_atoms, int *species, int *number_of_neighs, int number_of_neigh_list,
                             double *coordinates, double *desc,
                             DescriptorKind *descriptor_kind);
}


/*!
 * \brief Base class for all descriptors.
 *
 * This will be the parent class for all descriptors. To ensure compatibility with
 * Enzyme, class members and datastructures will be kept simple. Structures like lists
 * of lists are excruciatingly slow to compile through using enzyme. Enzyme plays well with
 * more C-like code. This is the main base class that will be passed on to the TorchMLModel driver. New descriptors can
 * be added by extending the class and implementing the DescriptorKind::compute() method.
 */
class Descriptor::DescriptorKind {

public:
    AvailableDescriptor descriptor_kind; //!< Kind of instantiated descriptor, will be used in creating clone for AD
    std::string descriptor_param_file; //!< Full path to descriptor parameter file.
    int width; //!< Dimension of the descriptor

    DescriptorKind() = default;

    static DescriptorKind *
    initDescriptor(AvailableDescriptor availableDescriptorKind); //!< Initialize an empty descriptor of a kind.
    //!< Used in creating virtual objects for AD.
    /*!
     * This function will initialize a required descriptor kind class from hyperparameter file generated by KIM-API.
     * @param file_name Fully qualified path of ascii file containing information required by the class. Usually it is
     * left for the class to figure out how to use it. In current context, all classes have a constructor,
     * which takes in the file name, and reads the file for instantiation.
     * @param availableDescriptorKind Instantiate required class of kind availableDescriptorKind  and read
     * parameters form file file_name
     * @return Pointer to instantiated class.
     * \see SymmetryFunctions::initFromFile() Bispectrum::initFromFile()
     */
    static DescriptorKind *initDescriptor(std::string &file_name,
                                          AvailableDescriptor availableDescriptorKind);

    /*!
     * compute() method contains the logic for calculating the descriptor for a single molecular environment.
     * Expected inputs include superset of all required parameters. In case some other information is needed in future,
     * we can try and expand it. For making it compatible with Enzyme few recommendations are, i) make it more C-like,
     * ii) avoid calling functions which dynamically allocate memory, and iii) avoid initializing any class members
     * inside the compute() function. A more "functional" programming principle is recommended. It is advisable
     * to avoid too much of memory allocation-deallocation, and rather use the constructor to allocate maximum
     * required memory and use it accordingly.
     *
     * @param index index of atom for which descriptor is to be computed. (i)
     * @param n_atoms total number of contributing atoms (This parameter is usually not used, as neighbor list contain
     * enough information)
     * @param species pointer to integer array containing species code for all atoms, note that this is usually numbered
     * 0...n_species - 1, and have no relation with atomic number. So if atomic numbers are needed, you need to
     * accordingly assign them inside compute() function.
     * @param neighbor_lists pointer to integer array containing indices of all the neighbor atoms
     * @param number_of_neighbors integer containing length of neighbor list
     * @param coordinates array containing all of the required coordinates, usually it is much simple to just pass
     * entire coordinate list. Then during differentiation you can avoid additional step of manually adding individual
     * derivatives.
     * @param desc pointer pointing to location where descriptor will be saved. It should have at least ::width size
     * memory available, else you will run into segfaults.
     */
    virtual void compute(int index,
                         int n_atoms,
                         int *species,
                         int *neighbor_lists,
                         int number_of_neighbors,
                         double *coordinates,
                         double *desc) = 0;

    // TODO: Use const pointers
    // virtual void clone_empty(DescriptorKind * descriptorKind){};
    // TODO: Cant make it virtual, enzyme segfaults. But every class must have
    // empty constructor to differentiate against
    // Due to enzyme issue I cannot yet make destructor as virtual. Therefore to prevent memory leak, the model driver
    // needs to reinterpret the pointer before destroying. It will be bit ugly hack but necessary at this moment
     ~DescriptorKind();
//    virtual ~DescriptorKind();
    /*!
     */
};


#endif // __DESCRIPTOR_HPP_