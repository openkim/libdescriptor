// Generic header for inclusion and library design

#ifndef __DESCRIPTOR_DR_HPP_
#define __DESCRIPTOR_DR_HPP_

#include <map>
#include <string>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <numeric>

#include <cmath>
#include <cstring>
#include <map>
#include <set>
#include <memory>
#include <algorithm>
#include <fstream>

/*!
 * \brief Namespace containing all required classes and methods for generating descriptors, as well as their derivatives.
 *
 * The namespace Descriptor_dr contains wrappers to all of the implemented descriptors.
 * The wrapper functions are used, instead of more direct inheritance, as Enzyme as of now (v0.0.41)
 * cannot accept pointer to class member function. This design might get more streamlined, as more mature version of Enzyme
 * is released.
 */
namespace Descriptor_dr {
    /// This enumerator lists all the kind of descriptors that are available at run time.
    enum AvailableDescriptor {
        KindSymmetryFunctions, //!< For selecting Behler Symmetry Functions (SymmetryFunctions)
//        KindBispectrum, //!< For selecting Bispectrum descriptor
//        KindSOAP, //!< For selecting Smooth Overlap of Atomic Position (SOAP) descriptor
//        KindXi //!< For selecting Xi descriptor
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

//    void num_gradient_single_atom(int index, int n_atoms, int *species, int *neighbor_list, int number_of_neighs,
//                              double *coordinates, double *d_coordinates, double *dE_dzeta,
//                              DescriptorKind *descriptor_to_diff);

    /*!
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
    void jacobian(int n_atoms, int n_total_atoms, int *species, int *neighbor_list, int *number_of_neighs,
                         double *coordinates, double *J_coordinates, DescriptorKind *descriptor_to_diff);

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
class Descriptor_dr::DescriptorKind {

public:
    AvailableDescriptor descriptor_kind; //!< Kind of instantiated descriptor, will be used in creating clone for AD
    std::string descriptor_param_file; //!< Full path to descriptor parameter file.
    int width=-1; //!< Dimension of the descriptor

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
                         double *distances, // natoms x natoms
                         double *desc) = 0;

    // TODO: Use const pointers
    // virtual void clone_empty(DescriptorKind * descriptorKind){};
    // TODO: Cant make it virtual, enzyme segfaults. But every class must have
    // empty constructor to differentiate against
    // Due to enzyme issue I cannot yet make destructor as virtual. Therefore to prevent memory leak, the model driver
    // needs to reinterpret the pointer before destroying. It will be bit ugly hack but necessary at this moment
//     ~DescriptorKind();
    virtual ~DescriptorKind();
    /*!
     */
     // ***************************************************************************************************************
     // Specialized descriptor initializing overloads
     // **************************************************************************************************************

    // ********Symmetry Functions********
    static DescriptorKind *
    initDescriptor(AvailableDescriptor availableDescriptorKind, std::vector<std::string> *species,
                   std::string *cutoff_function, double *cutoff_matrix,
                   std::vector<std::string> *symmetry_function_types, std::vector<int> *symmetry_function_sizes,
                   std::vector<double> *symmetry_function_parameters);

//    // ********Bispectrum********
//    static DescriptorKind *
//    initDescriptor(AvailableDescriptor availableDescriptorKind, double rfac0_in, int twojmax_in, int diagonalstyle_in,
//                   int use_shared_arrays_in, double rmin0_in, int switch_flag_in, int bzero_flag_in,
//                   double * cutoff_array, std::vector<std::string> * species, std::vector<double> * weights);
//
//    // ********SOAP********
//    static DescriptorKind *
//    initDescriptor(AvailableDescriptor availableDescriptorKind, int n_max, int l_max, double cutoff,
//                               std::vector<std::string> &species, std::string radial_basis, double eta);
//
//    // ********Xi********
//    DescriptorKind *
//    initDescriptor(AvailableDescriptor availableDescriptorKind, int q, double cutoff, std::vector<std::string> &species,
//                   std::string& radial_basis);
};

/*!
 * \brief This function formats messages, filename, line number and function
 * name into an std::ostringstream object
 *
 * \param message1      Starting message
 * \param fileName      File name
 * \param lineNumber    Line number
 * \param functionName  Function name
 * \param message2      Ending message
 *
 * \returns The combined std::ostringstream object as a string
 */
std::string
FormatMessageFileLineFunctionMessage(std::string const & message1,
                                     std::string const & fileName,
                                     long lineNumber,
                                     std::string const & functionName,
                                     std::string const & message2);

#ifdef LOG_ERROR
#undef LOG_ERROR
#endif

/*!
 * \brief Helper macro for printing error message
 *
 */
#define LOG_ERROR(msg)                                           \
  {                                                              \
    std::ostringstream ss;                                       \
    ss << msg;                                                   \
    std::string _Messagef_(FormatMessageFileLineFunctionMessage( \
        "Error ", __FILE__, __LINE__, __FUNCTION__, ss.str()));  \
    std::cerr << _Messagef_;                                     \
  }

/*! \class _Array_Basic The basic STL like container similar to `std::vector` to
 * handle multi-dimensional arrays.
 *
 * \brief An STL like container similar to <a
 * href="https://en.cppreference.com/w/cpp/container/vector">std::vector</a>
 * that encapsulates dynamic size arrays in a sequence container
 *
 * \tparam DataType The type of the elements. Default (double)
 */
template<class DataType = double>
class _Array_Basic
{
 public:
  /*!
   * \brief Construct a new _Array_Basic object
   *
   */
  _Array_Basic();

  /*!
   * \brief Construct a new _Array_Basic object
   *
   * \param count The size of the container
   */
  _Array_Basic(std::size_t const count);

  /*!
   * \brief Construct a new _Array_Basic object
   *
   * \param count The size of the container
   * \param value The value to initialize elements of the container with
   */
  _Array_Basic(std::size_t const count, DataType const value);

  /*!
   * \brief Construct a new _Array_Basic object
   *
   * \param count The size of the container
   * \param array The array of data to initialize elements of the container
   * with.
   */
  _Array_Basic(std::size_t const count, DataType const * array);

  /*!
   * \brief Construct a new _Array_Basic object
   * Copy constructor. Constructs the container with the copy of the contents of
   * other.
   *
   * \param other Another container to be used as source to initialize the
   * elements of the container with
   */
  _Array_Basic(_Array_Basic<DataType> const & other);

  /*!
   * \brief Construct a new _Array_Basic object
   * Move constructor. Constructs the container with the contents of other using
   * move semantics.
   *
   * \param other Another container to be used as source to initialize the
   * elements of the container with
   */
  _Array_Basic(_Array_Basic<DataType> && other);

  /*!
   * \brief Destroy the _Array_Basic object
   *
   */
  ~_Array_Basic();

  /*!
   * \brief Copy assignment operator. Replaces the contents with a copy of the
   * contents of other
   *
   * \param other Another container to use as data source
   *
   * \return _Array_Basic<DataType>&
   */
  _Array_Basic<DataType> & operator=(_Array_Basic<DataType> const & other);

  /*!
   * \brief Move assignment operator. Replaces the contents with those of other
   * using move semantics
   *
   * \param other Another container to use as data source
   *
   * \return _Array_Basic<DataType>&
   */
  _Array_Basic<DataType> & operator=(_Array_Basic<DataType> && other);

  /*!
   * \brief Returns pointer to the underlying array serving as element storage.
   *
   * \return DataType*
   */
  inline DataType const * data() const noexcept;
  inline DataType * data() noexcept;

  /*!
   * \brief Returns the number of elements in the container.
   *
   * \return std::size_t
   */
  inline std::size_t size() const;

  /*!
   * \brief Erases all elements from the container.
   *
   */
  inline void clear() noexcept;

  /*!
   * \brief Requests the removal of unused capacity.
   *
   */
  inline void shrink_to_fit();

  /*!
   * \brief Returns the number of elements that the container has currently
   * allocated space for.
   *
   * \return std::size_t
   */
  inline std::size_t capacity() const noexcept;

  /*!
   * \brief Appends the given element value to the end of the container.
   *
   * \param value
   */
  inline void push_back(DataType const & value);
  inline void push_back(DataType && value);

 protected:
  /*!
   * \brief Check the index range based on container size
   *
   * \param _n Index
   */
  inline void _range_check(int _n) const;

  /*!
   * \brief Check the index range based on input size
   *
   * \param _n Index
   * \param tsize Input size
   */
  inline void _range_check(int _n, std::size_t tsize) const;

 protected:
  /*! Dynamic contiguous array */
  std::vector<DataType> m;
};

/*! \class Array1DView A 1-dimensional STL like container.
 *
 * \brief A 1-dimensional STL like container that encapsulates dynamic size
 * arrays in a sequence container
 *
 * \tparam DataType The type of the elements. Default (double)
 */
template<class DataType = double>
class Array1DView
{
 public:
  /*!
   * \brief Construct a new Array1DView object
   *
   * \param count The size of the container
   * \param array The array of data to initialize elements of the container
   * with.
   */
  Array1DView(std::size_t const count, DataType * array);
  Array1DView(std::size_t const count, DataType const * array);

  /*!
   * \brief Construct a new Array1DView object
   * Copy constructor. Constructs the container with the copy of the contents of
   * other.
   *
   * \param other Another container to be used as source to initialize the
   * elements of the container with
   */
  Array1DView(Array1DView<DataType> const & other);

  /*!
   * \brief Destroy the Array1DView object
   *
   */
  ~Array1DView();

  /*!
   * \brief Returns pointer to the underlying array serving as element storage.
   *
   * \return DataType*
   */
  inline DataType const * data() const noexcept;
  inline DataType * data() noexcept;

  /*!
   * \brief Returns the element at specified location \b i.
   * No bounds checking is performed.
   *
   * \param i Position of the element to return
   *
   * \return const DataType The requested element.
   */
  inline const DataType operator()(int i) const;
  inline DataType & operator()(int i);

  /*!
   * \brief Returns the element at specified location \b i , with bounds
   * checking.
   *
   * \param i Position of the element to return
   *
   * \return const DataType The requested element.
   */
  inline DataType const at(int i) const;
  inline DataType & at(int i);

  /*!
   * \brief Returns the element at specified location \b i.
   * No bounds checking is performed.
   *
   * \param i Position of the element to return
   *
   * \return const DataType The requested element.
   */
  const DataType operator[](int i) const;
  DataType & operator[](int i);

 private:
  Array1DView() = delete;

  Array1DView<DataType> & operator=(Array1DView<DataType> const & other)
      = delete;

  Array1DView<DataType> & operator=(Array1DView<DataType> && other) = delete;

 protected:
  /*!
   * \brief Check the index range based on input size
   *
   * \param _n Index
   * \param tsize Input size
   */
  inline void _range_check(int _n, std::size_t tsize) const;

 protected:
  /*! The extent of the container in the 1st mode */
  std::size_t _extentZero;

  /*! Data pointer */
  DataType * const m;
};

template<class DataType = double>
class Array2DView
{
 public:
  /*!
   * \brief Construct a new Array2DView object
   *
   * \param extentZero The extent of the container in the 1st mode
   * \param extentOne The extent of the container in the 2nd mode
   * \param array The array of data to set the pointer of the container to it.
   */
  Array2DView(std::size_t const extentZero,
              std::size_t const extentOne,
              DataType * array);

  Array2DView(std::size_t const extentZero,
              std::size_t const extentOne,
              DataType const * array);

  /*!
   * \brief Construct a new Array2DView object
   * Copy constructor. Constructs the container with the copy of the contents of
   * other.
   *
   * \param other Another container to be used as source to initialize the
   * elements of the container with
   */
  Array2DView(Array2DView<DataType> const & other);

  /*!
   * \brief Destroy the Array2D object
   *
   */
  ~Array2DView();

  /*!
   * \brief Returns pointer to the underlying array serving as element storage.
   *
   * \return DataType*
   */
  inline DataType const * data() const noexcept;
  inline DataType * data() noexcept;

  inline Array1DView<DataType> data_1D(int i);

  /*!
   * \brief Returns the element at specified location \b (i, j).
   * No bounds checking is performed.
   *
   * \param i Position of the element in the 1st mode
   * \param j Position of the element in the 2nd mode
   *
   * \return const DataType The requested element.
   */
  inline const DataType operator()(int i, int j) const;
  inline DataType & operator()(int i, int j);

  /*!
   * \brief Returns the element at specified location \b (i, j) , with bounds
   * checking.
   *
   * \param i Position of the element in the 1st mode
   * \param j Position of the element in the 2nd mode
   *
   * \return const DataType The requested element.
   */
  inline DataType const at(int i, int j) const;
  inline DataType & at(int i, int j);

  /*! \class j_operator A helper class to provide multidimensional array access
   * semantics.
   *
   * \brief To provide 2-dimensional array access semantics, operator[] has to
   * return a reference to a 1D vector, which has to have its own operator[]
   * which returns a reference to the element.
   */
  class j_operator
  {
   public:
    /*!
     * \brief Construct a new j_operator object
     *
     * \param _array Refernce to Array2D class
     * \param i Position of the element in the 1st mode
     */
    j_operator(Array2DView<DataType> & _array, int i);

    /*!
     * \brief Provide array-like access and returns the element at specified
     * location \b [i][j]. No bounds checking is performed.
     *
     * \param j Position of the element in the 2nd mode
     *
     * \return const DataType The requested element.
     */
    const DataType operator[](int j) const;
    DataType & operator[](int j);

   private:
    /*! Refernce to Array2D class */
    Array2DView<DataType> & j_array;

    std::size_t _i;
  };

  /*!
   * \brief Provide array-like access and returns the element at specified
   * location \b [i][j]. No bounds checking is performed.
   *
   * \param i Position of the element in the 1st mode
   * \param j Position of the element in the 2nd mode
   *
   * \return const DataType The requested element.
   *
   * \note
   * To provide multidimensional array access semantics, we are using multiple
   * overloads for \code operator[] \endcode . For speed one should avoid this
   * complexity, uses \code operator() \endcode as \code (i, j) \endcode
   * directly.
   */
  const j_operator operator[](int i) const;
  j_operator operator[](int i);

 private:
  Array2DView() = delete;

  Array2DView<DataType> & operator=(Array2DView<DataType> const & other)
      = delete;

  Array2DView<DataType> & operator=(Array2DView<DataType> && other) = delete;

 protected:
  /*!
   * \brief Check the index range based on input size
   *
   * \param _n Index
   * \param tsize Input size
   */
  inline void _range_check(int _n, std::size_t tsize) const;

 protected:
  /*! The extent of the container in the 1st mode */
  std::size_t _extentZero;

  /*! The extent of the container in the 2nd mode */
  std::size_t _extentOne;

  /*! Data pointer */
  DataType * const m;
};

template<class DataType = double>
class Array3DView
{
 public:
  Array3DView(std::size_t const extentZero,
              std::size_t const extentOne,
              std::size_t const extentTwo,
              DataType * array);

  Array3DView(std::size_t const extentZero,
              std::size_t const extentOne,
              std::size_t const extentTwo,
              DataType const * array);

  Array3DView(Array3DView<DataType> const & other);

  ~Array3DView();

  inline DataType const * data() const noexcept;
  inline DataType * data() noexcept;

  inline Array2DView<DataType> data_2D(int i);

  inline Array1DView<DataType> data_1D(int i, int j);

  inline const DataType operator()(int i, int j, int k) const;
  inline DataType & operator()(int i, int j, int k);

  inline DataType const at(int i, int j, int k) const;
  inline DataType & at(int i, int j, int k);

  class j_operator
  {
   public:
    j_operator(Array3DView<DataType> & _array, int i);

    class k_operator
    {
     public:
      k_operator(Array3DView<DataType> & _array, int i, int j);

      const DataType operator[](int k) const;
      DataType & operator[](int k);

     private:
      Array3DView<DataType> & k_array;

      std::size_t _i;
      std::size_t _j;
    };

    const k_operator operator[](int j) const;

    k_operator operator[](int j);

   private:
    Array3DView<DataType> & j_array;

    std::size_t _i;
  };

  const j_operator operator[](int i) const;
  j_operator operator[](int i);

 private:
  Array3DView() = delete;

  Array3DView<DataType> & operator=(Array3DView<DataType> const & other)
      = delete;

  Array3DView<DataType> & operator=(Array3DView<DataType> && other) = delete;

 protected:
  /*!
   * \brief Check the index range based on input size
   *
   * \param _n Index
   * \param tsize Input size
   */
  inline void _range_check(int _n, std::size_t tsize) const;

 protected:
  std::size_t _extentZero;
  std::size_t _extentOne;
  std::size_t _extentTwo;

  /*! Data pointer */
  DataType * const m;
};

template<class DataType = double>
class Array4DView
{
 public:
  Array4DView(std::size_t const extentZero,
              std::size_t const extentOne,
              std::size_t const extentTwo,
              std::size_t const extentThree,
              DataType * array);

  Array4DView(std::size_t const extentZero,
              std::size_t const extentOne,
              std::size_t const extentTwo,
              std::size_t const extentThree,
              DataType const * array);

  Array4DView(Array4DView<DataType> const & other);

  ~Array4DView();

  inline DataType const * data() const noexcept;
  inline DataType * data() noexcept;

  inline Array3DView<DataType> data_3D(int i);

  inline Array2DView<DataType> data_2D(int i, int j);

  inline Array1DView<DataType> data_1D(int i, int j, int k);

  inline const DataType operator()(int i, int j, int k, int l) const;
  inline DataType & operator()(int i, int j, int k, int l);

  inline DataType const at(int i, int j, int k, int l) const;
  inline DataType & at(int i, int j, int k, int l);

  class j_operator
  {
   public:
    j_operator(Array4DView<DataType> & _array, int i);

    class k_operator
    {
     public:
      k_operator(Array4DView<DataType> & _array, int i, int j);

      class l_operator
      {
       public:
        l_operator(Array4DView<DataType> & _array, int i, int j, int k);

        const DataType operator[](int l) const;
        DataType & operator[](int l);

       private:
        Array4DView<DataType> & l_array;

        std::size_t _i;
        std::size_t _j;
        std::size_t _k;
      };

      const l_operator operator[](int k) const;
      l_operator operator[](int k);

     private:
      Array4DView<DataType> & k_array;

      std::size_t _i;
      std::size_t _j;
    };

    const k_operator operator[](int j) const;
    k_operator operator[](int j);

   private:
    Array4DView<DataType> & j_array;

    std::size_t _i;
  };

  const j_operator operator[](int i) const;
  j_operator operator[](int i);

 private:
  Array4DView() = delete;

  Array4DView<DataType> & operator=(Array4DView<DataType> const & other)
      = delete;

  Array4DView<DataType> & operator=(Array4DView<DataType> && other) = delete;

 protected:
  /*!
   * \brief Check the index range based on input size
   *
   * \param _n Index
   * \param tsize Input size
   */
  inline void _range_check(int _n, std::size_t tsize) const;

 protected:
  std::size_t _extentZero;
  std::size_t _extentOne;
  std::size_t _extentTwo;
  std::size_t _extentThree;

  DataType * const m;
};

/*! \class Array2D A 2-dimensional STL like container.
 *
 * \brief A 2-dimensional STL like container that encapsulates dynamic size
 * arrays with a 2-dimensional shape in a sequence container
 *
 * \tparam DataType The type of the elements. Default (double)
 */
template<class DataType = double>
class Array2D : public _Array_Basic<DataType>
{
 public:
  /*!
   * \brief Construct a new Array2D object
   *
   */
  Array2D();

  /*!
   * \brief Construct a new Array2D object
   *
   * \param extentZero The extent of the container in the 1st mode
   * \param extentOne The extent of the container in the 2nd mode
   */
  Array2D(std::size_t const extentZero, std::size_t const extentOne);

  /*!
   * \brief Construct a new Array2D object
   *
   * \param extentZero The extent of the container in the 1st mode
   * \param extentOne The extent of the container in the 2nd mode
   * \param value The value to initialize elements of the container with
   */
  Array2D(std::size_t const extentZero,
          std::size_t const extentOne,
          DataType const value);

  /*!
   * \brief Construct a new Array2D object
   *
   * \param extentZero The extent of the container in the 1st mode
   * \param extentOne The extent of the container in the 2nd mode
   * \param array The array of data to initialize elements of the container in a
   * row-major format.
   */
  Array2D(std::size_t const extentZero,
          std::size_t const extentOne,
          DataType const * array);

  /*!
   * \brief Construct a new Array2D object
   * Copy constructor. Constructs the container with the copy of the contents of
   * other.
   *
   * \param other Another container to be used as source to initialize the
   * elements of the container with
   */
  Array2D(Array2D<DataType> const & other);

  /*!
   * \brief Construct a new Array2D object
   * Move constructor. Constructs the container with the contents of other using
   * move semantics.
   *
   * \param other Another container to be used as source to initialize the
   * elements of the container with
   */
  Array2D(Array2D<DataType> && other);

  /*!
   * \brief Destroy the Array2D object
   *
   */
  ~Array2D();

  /*!
   * \brief Copy assignment operator. Replaces the contents with a copy of the
   * contents of other
   *
   * \param other Another container to use as data source
   *
   * \return Array2D<DataType>&
   */
  Array2D<DataType> & operator=(Array2D<DataType> const & other);

  /*!
   * \brief Move assignment operator. Replaces the contents with those of other
   * using move semantics
   *
   * \param other Another container to use as data source
   *
   * \return Array2D<DataType>&
   */
  Array2D<DataType> & operator=(Array2D<DataType> && other);

  /*!
   * \brief Returns Array1DView to the underlying starting element at row \c i.
   * \sa Array1DView
   *
   * \param i Row index
   * \return Array1DView<DataType>
   */
  inline Array1DView<DataType> data_1D(int i);

  /*!
   * \brief Resizes the container to contain \c extentZero times \c extentOne
   * elements.
   *
   * \param extentZero
   * \param extentOne
   */
  inline void resize(int const extentZero, int const extentOne);

  /*!
   * \brief Resizes the container to contain \c extentZero times \c extentOne
   * elements.
   *
   * \param extentZero
   * \param extentOne
   * \param new_value The new value to initialize the new elements with
   */
  inline void
  resize(int const extentZero, int const extentOne, DataType const new_value);

  /*!
   * \brief Resizes the container to contain \c extentZero times \c extentOne
   * elements.
   *
   * \param extentZero
   * \param extentOne
   * \param new_array The new array of data to initialize elements of the
   * container with.
   */
  inline void
  resize(int const extentZero, int const extentOne, DataType const * new_array);

  /*!
   * \brief Returns the element at specified location \b (i, j).
   * No bounds checking is performed.
   *
   * \param i Position of the element in the 1st mode
   * \param j Position of the element in the 2nd mode
   *
   * \return const DataType The requested element.
   */
  inline const DataType operator()(int i, int j) const;
  inline DataType & operator()(int i, int j);

  /*!
   * \brief Returns the element at specified location \b (i, j) , with bounds
   * checking.
   *
   * \param i Position of the element in the 1st mode
   * \param j Position of the element in the 2nd mode
   *
   * \return const DataType The requested element.
   */
  inline DataType const at(int i, int j) const;
  inline DataType & at(int i, int j);

  /*! \class j_operator A helper class to provide multidimensional array access
   * semantics.
   *
   * \brief To provide 2-dimensional array access semantics, operator[] has to
   * return a reference to a 1D vector, which has to have its own operator[]
   * which returns a reference to the element.
   */
  class j_operator
  {
   public:
    /*!
     * \brief Construct a new j_operator object
     *
     * \param _array Refernce to Array2D class
     * \param i Position of the element in the 1st mode
     */
    j_operator(Array2D<DataType> & _array, int i);

    /*!
     * \brief Provide array-like access and returns the element at specified
     * location \b [i][j]. No bounds checking is performed.
     *
     * \param j Position of the element in the 2nd mode
     *
     * \return const DataType The requested element.
     */
    const DataType operator[](int j) const;
    DataType & operator[](int j);

   private:
    /*! Refernce to Array2D class */
    Array2D<DataType> & j_array;

    std::size_t _i;
  };

  /*!
   * \brief Provide array-like access and returns the element at specified
   * location \b [i][j]. No bounds checking is performed.
   *
   * \param i Position of the element in the 1st mode
   * \param j Position of the element in the 2nd mode
   *
   * \return const DataType The requested element.
   *
   * \note
   * To provide multidimensional array access semantics, we are using multiple
   * overloads for \code operator[] \endcode . For speed one should avoid this
   * complexity, uses \code operator() \endcode as \code (i, j) \endcode
   * directly.
   */
  const j_operator operator[](int i) const;
  j_operator operator[](int i);

 protected:
  /*! The extent of the container in the 1st mode */
  std::size_t _extentZero;

  /*! The extent of the container in the 2nd mode */
  std::size_t _extentOne;
};

template<class DataType = double>
class Array3D : public _Array_Basic<DataType>
{
 public:
  Array3D();

  Array3D(std::size_t const extentZero,
          std::size_t const extentOne,
          std::size_t const extentTwo);

  Array3D(std::size_t const extentZero,
          std::size_t const extentOne,
          std::size_t const extentTwo,
          DataType const value);

  Array3D(std::size_t const extentZero,
          std::size_t const extentOne,
          std::size_t const extentTwo,
          DataType const * array);

  Array3D(Array3D<DataType> const & other);

  Array3D(Array3D<DataType> && other);

  ~Array3D();

  Array3D<DataType> & operator=(Array3D<DataType> const & other);

  Array3D<DataType> & operator=(Array3D<DataType> && other);

  inline Array2DView<DataType> data_2D(int i);

  inline Array1DView<DataType> data_1D(int i, int j);

  inline void
  resize(int const extentZero, int const extentOne, int const extentTwo);

  inline void resize(int const extentZero,
                     int const extentOne,
                     int const extentTwo,
                     DataType const new_value);

  inline void resize(int const extentZero,
                     int const extentOne,
                     int const extentTwo,
                     DataType const * new_array);

  inline const DataType operator()(int i, int j, int k) const;
  inline DataType & operator()(int i, int j, int k);

  inline DataType const at(int i, int j, int k) const;
  inline DataType & at(int i, int j, int k);

  class j_operator
  {
   public:
    j_operator(Array3D<DataType> & _array, int i);

    class k_operator
    {
     public:
      k_operator(Array3D<DataType> & _array, int i, int j);

      const DataType operator[](int k) const;
      DataType & operator[](int k);

     private:
      Array3D<DataType> & k_array;

      std::size_t _i;
      std::size_t _j;
    };

    const k_operator operator[](int j) const;
    k_operator operator[](int j);

   private:
    Array3D<DataType> & j_array;

    std::size_t _i;
  };

  const j_operator operator[](int i) const;
  j_operator operator[](int i);

 protected:
  std::size_t _extentZero;
  std::size_t _extentOne;
  std::size_t _extentTwo;
};

template<class DataType = double>
class Array4D : public _Array_Basic<DataType>
{
 public:
  Array4D();

  Array4D(std::size_t const extentZero,
          std::size_t const extentOne,
          std::size_t const extentTwo,
          std::size_t const extentThree);

  Array4D(std::size_t const extentZero,
          std::size_t const extentOne,
          std::size_t const extentTwo,
          std::size_t const extentThree,
          DataType const value);

  Array4D(std::size_t const extentZero,
          std::size_t const extentOne,
          std::size_t const extentTwo,
          std::size_t const extentThree,
          DataType const * array);

  Array4D(Array4D<DataType> const & other);

  Array4D(Array4D<DataType> && other);

  ~Array4D();

  Array4D<DataType> & operator=(Array4D<DataType> const & other);

  Array4D<DataType> & operator=(Array4D<DataType> && other);

  inline Array3DView<DataType> data_3D(int i);

  inline Array2DView<DataType> data_2D(int i, int j);

  inline Array1DView<DataType> data_1D(int i, int j, int k);

  inline void resize(int const extentZero,
                     int const extentOne,
                     int const extentTwo,
                     int const extentThree);

  inline void resize(int const extentZero,
                     int const extentOne,
                     int const extentTwo,
                     int const extentThree,
                     DataType const new_value);

  inline void resize(int const extentZero,
                     int const extentOne,
                     int const extentTwo,
                     int const extentThree,
                     DataType const * new_array);

  inline const DataType operator()(int i, int j, int k, int l) const;
  inline DataType & operator()(int i, int j, int k, int l);

  inline DataType const at(int i, int j, int k, int l) const;
  inline DataType & at(int i, int j, int k, int l);

  class j_operator
  {
   public:
    j_operator(Array4D<DataType> & _array, int i);

    class k_operator
    {
     public:
      k_operator(Array4D<DataType> & _array, int i, int j);

      class l_operator
      {
       public:
        l_operator(Array4D<DataType> & _array, int i, int j, int k);

        const DataType operator[](int l) const;
        DataType & operator[](int l);

       private:
        Array4D<DataType> & l_array;

        std::size_t _i;
        std::size_t _j;
        std::size_t _k;
      };

      const l_operator operator[](int k) const;
      l_operator operator[](int k);

     private:
      Array4D<DataType> & k_array;

      std::size_t _i;
      std::size_t _j;
    };

    const k_operator operator[](int j) const;
    k_operator operator[](int j);

   private:
    Array4D<DataType> & j_array;

    std::size_t _i;
  };

  const j_operator operator[](int i) const;
  j_operator operator[](int i);

 protected:
  std::size_t _extentZero;
  std::size_t _extentOne;
  std::size_t _extentTwo;
  std::size_t _extentThree;
};

template<class DataType = double>
class Array5D : public _Array_Basic<DataType>
{
 public:
  Array5D();

  Array5D(std::size_t const extentZero,
          std::size_t const extentOne,
          std::size_t const extentTwo,
          std::size_t const extentThree,
          std::size_t const extentFour);

  Array5D(std::size_t const extentZero,
          std::size_t const extentOne,
          std::size_t const extentTwo,
          std::size_t const extentThree,
          std::size_t const extentFour,
          DataType const value);

  Array5D(std::size_t const extentZero,
          std::size_t const extentOne,
          std::size_t const extentTwo,
          std::size_t const extentThree,
          std::size_t const extentFour,
          DataType const * array);

  Array5D(Array5D<DataType> const & other);

  Array5D(Array5D<DataType> && other);

  ~Array5D();

  Array5D<DataType> & operator=(Array5D<DataType> const & other);

  Array5D<DataType> & operator=(Array5D<DataType> && other);

  inline Array4DView<DataType> data_4D(int i);

  inline Array3DView<DataType> data_3D(int i, int j);

  inline Array2DView<DataType> data_2D(int i, int j, int k);

  inline Array1DView<DataType> data_1D(int i, int j, int k, int l);

  inline void resize(int const extentZero,
                     int const extentOne,
                     int const extentTwo,
                     int const extentThree,
                     int const extentFour);

  inline void resize(int const extentZero,
                     int const extentOne,
                     int const extentTwo,
                     int const extentThree,
                     int const extentFour,
                     DataType const new_value);

  inline void resize(int const extentZero,
                     int const extentOne,
                     int const extentTwo,
                     int const extentThree,
                     int const extentFour,
                     DataType const * new_array);

  inline const DataType operator()(int i, int j, int k, int l, int n) const;
  inline DataType & operator()(int i, int j, int k, int l, int n);

  inline DataType const at(int i, int j, int k, int l, int n) const;
  inline DataType & at(int i, int j, int k, int l, int n);

  class j_operator
  {
   public:
    j_operator(Array5D<DataType> & _array, int i);

    class k_operator
    {
     public:
      k_operator(Array5D<DataType> & _array, int i, int j);

      class l_operator
      {
       public:
        l_operator(Array5D<DataType> & _array, int i, int j, int k);

        class n_operator
        {
         public:
          n_operator(Array5D<DataType> & _array, int i, int j, int k, int l);

          const DataType operator[](int n) const;
          DataType & operator[](int n);

         private:
          Array5D<DataType> & n_array;

          std::size_t _i;
          std::size_t _j;
          std::size_t _k;
          std::size_t _l;
        };

        const n_operator operator[](int l) const;
        n_operator operator[](int l);

       private:
        Array5D<DataType> & l_array;

        std::size_t _i;
        std::size_t _j;
        std::size_t _k;
      };

      const l_operator operator[](int k) const;
      l_operator operator[](int k);

     private:
      Array5D<DataType> & k_array;

      std::size_t _i;
      std::size_t _j;
    };

    const k_operator operator[](int j) const;
    k_operator operator[](int j);

   private:
    Array5D<DataType> & j_array;

    std::size_t _i;
  };

  const j_operator operator[](int i) const;
  j_operator operator[](int i);

 protected:
  std::size_t _extentZero;
  std::size_t _extentOne;
  std::size_t _extentTwo;
  std::size_t _extentThree;
  std::size_t _extentFour;
};

// --------------------------- Implementation --------------------------- //

template<class DataType>
_Array_Basic<DataType>::_Array_Basic()
{
}

template<class DataType>
_Array_Basic<DataType>::_Array_Basic(std::size_t const count) :
    m(count, static_cast<DataType>(0))
{
}

template<class DataType>
_Array_Basic<DataType>::_Array_Basic(std::size_t const count,
                                     DataType const value) :
    m(count, value)
{
}

template<class DataType>
_Array_Basic<DataType>::_Array_Basic(std::size_t const count,
                                     DataType const * array) :
    m(array, array + count)
{
}

template<class DataType>
_Array_Basic<DataType>::_Array_Basic(_Array_Basic<DataType> const & other) :
    m(other.m)
{
}

template<class DataType>
_Array_Basic<DataType>::_Array_Basic(_Array_Basic<DataType> && other) :
    m(std::move(other.m))
{
}

template<class DataType>
_Array_Basic<DataType>::~_Array_Basic()
{
}

template<class DataType>
_Array_Basic<DataType> &
_Array_Basic<DataType>::operator=(_Array_Basic<DataType> const & other)
{
  m.resize(other.size());
  std::copy(other.m.begin(), other.m.end(), m.begin());
  return *this;
}

template<class DataType>
_Array_Basic<DataType> &
_Array_Basic<DataType>::operator=(_Array_Basic<DataType> && other)
{
  m = std::move(other.m);
  return *this;
}

template<class DataType>
inline DataType const * _Array_Basic<DataType>::data() const noexcept
{
  return m.data();
}

template<class DataType>
inline DataType * _Array_Basic<DataType>::data() noexcept
{
  return m.data();
}

template<class DataType>
inline std::size_t _Array_Basic<DataType>::size() const
{
  return m.size();
}

template<class DataType>
inline void _Array_Basic<DataType>::clear() noexcept
{
  m.clear();
}

template<class DataType>
inline void _Array_Basic<DataType>::shrink_to_fit()
{
  m.shrink_to_fit();
}

template<class DataType>
inline std::size_t _Array_Basic<DataType>::capacity() const noexcept
{
  return m.capacity();
}

template<class DataType>
inline void _Array_Basic<DataType>::push_back(DataType const & value)
{
  m.push_back(value);
}

template<class DataType>
inline void _Array_Basic<DataType>::push_back(DataType && value)
{
  m.push_back(value);
}

template<class DataType>
inline void _Array_Basic<DataType>::_range_check(int _n) const
{
  if (_n >= size())
  {
    LOG_ERROR("The input index is out of range! " + std::to_string(_n)
              + " >= " + std::to_string(size()));
    std::abort();
  }
}

template<class DataType>
inline void _Array_Basic<DataType>::_range_check(int _n,
                                                 std::size_t tsize) const
{
  if (_n >= tsize)
  {
    LOG_ERROR("The input index is out of range! " + std::to_string(_n)
              + " >= " + std::to_string(tsize));
    std::abort();
  }
}

template<class DataType>
Array1DView<DataType>::Array1DView(std::size_t const count, DataType * array) :
    _extentZero(count), m(array)
{
}

template<class DataType>
Array1DView<DataType>::Array1DView(Array1DView<DataType> const & other) :
    _extentZero(other._extentZero), m(other.m)
{
}

template<class DataType>
Array1DView<DataType>::~Array1DView()
{
}

template<class DataType>
inline DataType const * Array1DView<DataType>::data() const noexcept
{
  return m;
}

template<class DataType>
inline DataType * Array1DView<DataType>::data() noexcept
{
  return m;
}

template<class DataType>
inline const DataType Array1DView<DataType>::operator()(int i) const
{
  return m[i];
}

template<class DataType>
inline DataType & Array1DView<DataType>::operator()(int i)
{
  return m[i];
}

template<class DataType>
inline DataType & Array1DView<DataType>::at(int i)
{
  _range_check(i, _extentZero);
  return m[i];
}

template<class DataType>
inline DataType const Array1DView<DataType>::at(int i) const
{
  _range_check(i, _extentZero);
  return m[i];
}

template<class DataType>
const DataType Array1DView<DataType>::operator[](int i) const
{
  return m[i];
}

template<class DataType>
DataType & Array1DView<DataType>::operator[](int i)
{
  return m[i];
}

template<class DataType>
inline void Array1DView<DataType>::_range_check(int _n, std::size_t tsize) const
{
  if (_n >= tsize)
  {
    LOG_ERROR("The input index is out of range! " + std::to_string(_n)
              + " >= " + std::to_string(tsize));
    std::abort();
  }
}

template<class DataType>
Array2DView<DataType>::Array2DView(std::size_t const extentZero,
                                   std::size_t const extentOne,
                                   DataType * array) :
    _extentZero(extentZero), _extentOne(extentOne), m(array)
{
}

template<class DataType>
Array2DView<DataType>::Array2DView(std::size_t const extentZero,
                                   std::size_t const extentOne,
                                   DataType const * array) :
    _extentZero(extentZero),
    _extentOne(extentOne),
    m(const_cast<DataType *>(array))
{
}

template<class DataType>
Array2DView<DataType>::Array2DView(Array2DView<DataType> const & other) :
    _extentZero(other._extentZero), _extentOne(other._extentOne), m(other.m)
{
}

template<class DataType>
Array2DView<DataType>::~Array2DView()
{
}

template<class DataType>
inline DataType const * Array2DView<DataType>::data() const noexcept
{
  return m;
}

template<class DataType>
inline DataType * Array2DView<DataType>::data() noexcept
{
  return m;
}

template<class DataType>
inline Array1DView<DataType> Array2DView<DataType>::data_1D(int i)
{
  return Array1DView<DataType>(_extentOne, m + i * _extentOne);
}

template<class DataType>
inline const DataType Array2DView<DataType>::operator()(int i, int j) const
{
  std::size_t const _n = i * _extentOne + j;
  return m[_n];
}

template<class DataType>
inline DataType & Array2DView<DataType>::operator()(int i, int j)
{
  std::size_t const _n = i * _extentOne + j;
  return m[_n];
}

template<class DataType>
inline DataType & Array2DView<DataType>::at(int i, int j)
{
  _range_check(i, _extentZero);
  _range_check(j, _extentOne);
  std::size_t const _n = i * _extentOne + j;
  return m[_n];
}

template<class DataType>
inline DataType const Array2DView<DataType>::at(int i, int j) const
{
  _range_check(i, _extentZero);
  _range_check(j, _extentOne);
  std::size_t const _n = i * _extentOne + j;
  return m[_n];
}

template<class DataType>
Array2DView<DataType>::j_operator::j_operator(Array2DView<DataType> & _array,
                                              int i) :
    j_array(_array), _i(i)
{
}

template<class DataType>
const DataType Array2DView<DataType>::j_operator::operator[](int j) const
{
  std::size_t const _n = _i * j_array._extentOne + j;
  return j_array.m[_n];
}

template<class DataType>
DataType & Array2DView<DataType>::j_operator::operator[](int j)
{
  std::size_t const _n = _i * j_array._extentOne + j;
  return j_array.m[_n];
}

template<class DataType>
const typename Array2DView<DataType>::j_operator
Array2DView<DataType>::operator[](int i) const
{
  return j_operator(*this, i);
}

template<class DataType>
typename Array2DView<DataType>::j_operator
Array2DView<DataType>::operator[](int i)
{
  return j_operator(*this, i);
}

template<class DataType>
inline void Array2DView<DataType>::_range_check(int _n, std::size_t tsize) const
{
  if (_n >= tsize)
  {
    LOG_ERROR("The input index is out of range! " + std::to_string(_n)
              + " >= " + std::to_string(tsize));
    std::abort();
  }
}

template<class DataType>
Array3DView<DataType>::Array3DView(std::size_t const extentZero,
                                   std::size_t const extentOne,
                                   std::size_t const extentTwo,
                                   DataType * array) :
    _extentZero(extentZero),
    _extentOne(extentOne),
    _extentTwo(extentTwo),
    m(array)
{
}

template<class DataType>
Array3DView<DataType>::Array3DView(std::size_t const extentZero,
                                   std::size_t const extentOne,
                                   std::size_t const extentTwo,
                                   DataType const * array) :
    _extentZero(extentZero),
    _extentOne(extentOne),
    _extentTwo(extentTwo),
    m(const_cast<DataType *>(array))
{
}

template<class DataType>
Array3DView<DataType>::Array3DView(Array3DView<DataType> const & other) :
    _extentZero(other._extentZero),
    _extentOne(other._extentOne),
    _extentTwo(other._extentTwo),
    m(other.m)
{
}

template<class DataType>
Array3DView<DataType>::~Array3DView()
{
}

template<class DataType>
inline DataType const * Array3DView<DataType>::data() const noexcept
{
  return m;
}

template<class DataType>
inline DataType * Array3DView<DataType>::data() noexcept
{
  return m;
}

template<class DataType>
inline Array2DView<DataType> Array3DView<DataType>::data_2D(int i)
{
  return Array2DView<DataType>(
      _extentOne, _extentTwo, m + i * _extentOne * _extentTwo);
}

template<class DataType>
inline Array1DView<DataType> Array3DView<DataType>::data_1D(int i, int j)
{
  return Array1DView<DataType>(_extentTwo,
                               m + (i * _extentOne + j) * _extentTwo);
}

template<class DataType>
inline const DataType
Array3DView<DataType>::operator()(int i, int j, int k) const
{
  std::size_t const _n = (i * _extentOne + j) * _extentTwo + k;
  return m[_n];
}

template<class DataType>
inline DataType & Array3DView<DataType>::operator()(int i, int j, int k)
{
  std::size_t const _n = (i * _extentOne + j) * _extentTwo + k;
  return m[_n];
}

template<class DataType>
inline DataType const Array3DView<DataType>::at(int i, int j, int k) const
{
  _range_check(i, _extentZero);
  _range_check(j, _extentOne);
  _range_check(k, _extentTwo);
  std::size_t const _n = (i * _extentOne + j) * _extentTwo + k;
  return m[_n];
}

template<class DataType>
inline DataType & Array3DView<DataType>::at(int i, int j, int k)
{
  _range_check(i, _extentZero);
  _range_check(j, _extentOne);
  _range_check(k, _extentTwo);
  std::size_t const _n = (i * _extentOne + j) * _extentTwo + k;
  return m[_n];
}

template<class DataType>
Array3DView<DataType>::j_operator::j_operator(Array3DView<DataType> & _array,
                                              int i) :
    j_array(_array), _i(i)
{
}

template<class DataType>
Array3DView<DataType>::j_operator::k_operator::k_operator(
    Array3DView<DataType> & _array, int i, int j) :
    k_array(_array), _i(i), _j(j)
{
}

template<class DataType>
const DataType
Array3DView<DataType>::j_operator::k_operator::operator[](int k) const
{
  std::size_t const _n
      = (_i * k_array._extentOne + _j) * k_array._extentTwo + k;
  return k_array.m[_n];
}

template<class DataType>
DataType & Array3DView<DataType>::j_operator::k_operator::operator[](int k)
{
  std::size_t const _n
      = (_i * k_array._extentOne + _j) * k_array._extentTwo + k;
  return k_array.m[_n];
}

template<class DataType>
const typename Array3DView<DataType>::j_operator::k_operator
Array3DView<DataType>::j_operator::operator[](int j) const
{
  return k_operator(j_array, _i, j);
}

template<class DataType>
typename Array3DView<DataType>::j_operator::k_operator
Array3DView<DataType>::j_operator::operator[](int j)
{
  return k_operator(j_array, _i, j);
}

template<class DataType>
const typename Array3DView<DataType>::j_operator
Array3DView<DataType>::operator[](int i) const
{
  return j_operator(*this, i);
}

template<class DataType>
typename Array3DView<DataType>::j_operator
Array3DView<DataType>::operator[](int i)
{
  return j_operator(*this, i);
}

template<class DataType>
inline void Array3DView<DataType>::_range_check(int _n, std::size_t tsize) const
{
  if (_n >= tsize)
  {
    LOG_ERROR("The input index is out of range! " + std::to_string(_n)
              + " >= " + std::to_string(tsize));
    std::abort();
  }
}

template<class DataType>
Array4DView<DataType>::Array4DView(std::size_t const extentZero,
                                   std::size_t const extentOne,
                                   std::size_t const extentTwo,
                                   std::size_t const extentThree,
                                   DataType * array) :
    _extentZero(extentZero),
    _extentOne(extentOne),
    _extentTwo(extentTwo),
    _extentThree(extentThree),
    m(array)
{
}

template<class DataType>
Array4DView<DataType>::Array4DView(std::size_t const extentZero,
                                   std::size_t const extentOne,
                                   std::size_t const extentTwo,
                                   std::size_t const extentThree,
                                   DataType const * array) :
    _extentZero(extentZero),
    _extentOne(extentOne),
    _extentTwo(extentTwo),
    _extentThree(extentThree),
    m(const_cast<DataType *>(array))
{
}

template<class DataType>
Array4DView<DataType>::Array4DView(Array4DView<DataType> const & other) :
    _extentZero(other._extentZero),
    _extentOne(other._extentOne),
    _extentTwo(other._extentTwo),
    _extentThree(other._extentThree),
    m(other.m)
{
}

template<class DataType>
Array4DView<DataType>::~Array4DView()
{
}

template<class DataType>
inline DataType const * Array4DView<DataType>::data() const noexcept
{
  return m;
}

template<class DataType>
inline DataType * Array4DView<DataType>::data() noexcept
{
  return m;
}

template<class DataType>
inline Array3DView<DataType> Array4DView<DataType>::data_3D(int i)
{
  std::size_t const _n = i * _extentOne * _extentTwo * _extentThree;
  return Array3DView<DataType>(_extentOne, _extentTwo, _extentThree, m + _n);
}

template<class DataType>
inline Array2DView<DataType> Array4DView<DataType>::data_2D(int i, int j)
{
  std::size_t const _n = (i * _extentOne + j) * _extentTwo * _extentThree;
  return Array2DView<DataType>(_extentTwo, _extentThree, m + _n);
}

template<class DataType>
inline Array1DView<DataType> Array4DView<DataType>::data_1D(int i, int j, int k)
{
  std::size_t const _n = ((i * _extentOne + j) * _extentTwo + k) * _extentThree;
  return Array1DView<DataType>(_extentThree, m + _n);
}

template<class DataType>
inline const DataType
Array4DView<DataType>::operator()(int i, int j, int k, int l) const
{
  std::size_t const _n
      = ((i * _extentOne + j) * _extentTwo + k) * _extentThree + l;
  return m[_n];
}

template<class DataType>
inline DataType & Array4DView<DataType>::operator()(int i, int j, int k, int l)
{
  std::size_t const _n
      = ((i * _extentOne + j) * _extentTwo + k) * _extentThree + l;
  return m[_n];
}

template<class DataType>
inline DataType const
Array4DView<DataType>::at(int i, int j, int k, int l) const
{
  _range_check(i, _extentZero);
  _range_check(j, _extentOne);
  _range_check(k, _extentTwo);
  _range_check(l, _extentThree);
  std::size_t const _n
      = ((i * _extentOne + j) * _extentTwo + k) * _extentThree + l;
  return m[_n];
}

template<class DataType>
inline DataType & Array4DView<DataType>::at(int i, int j, int k, int l)
{
  _range_check(i, _extentZero);
  _range_check(j, _extentOne);
  _range_check(k, _extentTwo);
  _range_check(l, _extentThree);
  std::size_t const _n
      = ((i * _extentOne + j) * _extentTwo + k) * _extentThree + l;
  return m[_n];
}

template<class DataType>
Array4DView<DataType>::j_operator::j_operator(Array4DView<DataType> & _array,
                                              int i) :
    j_array(_array), _i(i)
{
}

template<class DataType>
Array4DView<DataType>::j_operator::k_operator::k_operator(
    Array4DView<DataType> & _array, int i, int j) :
    k_array(_array), _i(i), _j(j)
{
}

template<class DataType>
Array4DView<DataType>::j_operator::k_operator::l_operator::l_operator(
    Array4DView<DataType> & _array, int i, int j, int k) :
    l_array(_array), _i(i), _j(j), _k(k)
{
}

template<class DataType>
const DataType
Array4DView<DataType>::j_operator::k_operator::l_operator::operator[](
    int l) const
{
  std::size_t const _n
      = ((_i * l_array._extentOne + _j) * l_array._extentTwo + _k)
            * l_array._extentThree
        + l;
  return l_array.m[_n];
}

template<class DataType>
DataType &
Array4DView<DataType>::j_operator::k_operator::l_operator::operator[](int l)
{
  std::size_t const _n
      = ((_i * l_array._extentOne + _j) * l_array._extentTwo + _k)
            * l_array._extentThree
        + l;
  return l_array.m[_n];
}

template<class DataType>
const typename Array4DView<DataType>::j_operator::k_operator::l_operator
Array4DView<DataType>::j_operator::k_operator::operator[](int k) const
{
  return l_operator(k_array, _i, _j, k);
}

template<class DataType>
typename Array4DView<DataType>::j_operator::k_operator::l_operator
Array4DView<DataType>::j_operator::k_operator::operator[](int k)
{
  return l_operator(k_array, _i, _j, k);
}

template<class DataType>
const typename Array4DView<DataType>::j_operator::k_operator
Array4DView<DataType>::j_operator::operator[](int j) const
{
  return k_operator(j_array, _i, j);
}

template<class DataType>
typename Array4DView<DataType>::j_operator::k_operator
Array4DView<DataType>::j_operator::operator[](int j)
{
  return k_operator(j_array, _i, j);
}

template<class DataType>
const typename Array4DView<DataType>::j_operator
Array4DView<DataType>::operator[](int i) const
{
  return j_operator(*this, i);
}

template<class DataType>
typename Array4DView<DataType>::j_operator
Array4DView<DataType>::operator[](int i)
{
  return j_operator(*this, i);
}

template<class DataType>
inline void Array4DView<DataType>::_range_check(int _n, std::size_t tsize) const
{
  if (_n >= tsize)
  {
    LOG_ERROR("The input index is out of range! " + std::to_string(_n)
              + " >= " + std::to_string(tsize));
    std::abort();
  }
}

template<class DataType>
Array2D<DataType>::Array2D() :
    _Array_Basic<DataType>(), _extentZero(0), _extentOne(0)
{
}

template<class DataType>
Array2D<DataType>::Array2D(std::size_t const extentZero,
                           std::size_t const extentOne) :
    _Array_Basic<DataType>(extentZero * extentOne),
    _extentZero(extentZero),
    _extentOne(extentOne)
{
}

template<class DataType>
Array2D<DataType>::Array2D(std::size_t const extentZero,
                           std::size_t const extentOne,
                           DataType const value) :
    _Array_Basic<DataType>(extentZero * extentOne, value),
    _extentZero(extentZero),
    _extentOne(extentOne)
{
}

template<class DataType>
Array2D<DataType>::Array2D(std::size_t const extentZero,
                           std::size_t const extentOne,
                           DataType const * array) :
    _Array_Basic<DataType>(extentZero * extentOne, array),
    _extentZero(extentZero),
    _extentOne(extentOne)
{
}

template<class DataType>
Array2D<DataType>::Array2D(Array2D<DataType> const & other) :
    _Array_Basic<DataType>(other),
    _extentZero(other._extentZero),
    _extentOne(other._extentOne)
{
}

template<class DataType>
Array2D<DataType>::Array2D(Array2D<DataType> && other) :
    _Array_Basic<DataType>(std::move(other)),
    _extentZero(other._extentZero),
    _extentOne(other._extentOne)
{
}

template<class DataType>
Array2D<DataType>::~Array2D()
{
}

template<class DataType>
Array2D<DataType> &
Array2D<DataType>::operator=(Array2D<DataType> const & other)
{
  _Array_Basic<DataType>::operator=(other);
  _extentZero = other._extentZero;
  _extentOne = other._extentOne;
  return *this;
}

template<class DataType>
Array2D<DataType> & Array2D<DataType>::operator=(Array2D<DataType> && other)
{
  _Array_Basic<DataType>::operator=(std::move(other));
  _extentZero = other._extentZero;
  _extentOne = other._extentOne;
  return *this;
}

template<class DataType>
inline Array1DView<DataType> Array2D<DataType>::data_1D(int i)
{
  return Array1DView<DataType>(_extentOne, this->m.data() + i * _extentOne);
}

template<class DataType>
inline void Array2D<DataType>::resize(int const extentZero, int const extentOne)
{
  _extentZero = extentZero;
  _extentOne = extentOne;
  std::size_t const _n = _extentZero * _extentOne;
  this->m.resize(_n, static_cast<DataType>(0));
}

template<class DataType>
inline void Array2D<DataType>::resize(int const extentZero,
                                      int const extentOne,
                                      DataType const new_value)
{
  _extentZero = extentZero;
  _extentOne = extentOne;
  std::size_t const _n = _extentZero * _extentOne;
  this->m.resize(_n, new_value);
}

template<class DataType>
inline void Array2D<DataType>::resize(int const extentZero,
                                      int const extentOne,
                                      DataType const * new_array)
{
  _extentZero = extentZero;
  _extentOne = extentOne;
  std::size_t const _n = _extentZero * _extentOne;
  this->m.resize(_n);
  std::copy(new_array, new_array + _n, this->m.data());
}

template<class DataType>
inline const DataType Array2D<DataType>::operator()(int i, int j) const
{
  std::size_t const _n = i * _extentOne + j;
  return this->m[_n];
}

template<class DataType>
inline DataType & Array2D<DataType>::operator()(int i, int j)
{
  std::size_t const _n = i * _extentOne + j;
  return this->m[_n];
}

template<class DataType>
inline DataType & Array2D<DataType>::at(int i, int j)
{
  this->_range_check(i, _extentZero);
  this->_range_check(j, _extentOne);
  std::size_t const _n = i * _extentOne + j;
  return this->m[_n];
}

template<class DataType>
inline DataType const Array2D<DataType>::at(int i, int j) const
{
  this->_range_check(i, _extentZero);
  this->_range_check(j, _extentOne);
  std::size_t const _n = i * _extentOne + j;
  return this->m[_n];
}

template<class DataType>
Array2D<DataType>::j_operator::j_operator(Array2D<DataType> & _array, int i) :
    j_array(_array), _i(i)
{
}

template<class DataType>
const DataType Array2D<DataType>::j_operator::operator[](int j) const
{
  std::size_t const _n = _i * j_array._extentOne + j;
  return j_array.m[_n];
}

template<class DataType>
DataType & Array2D<DataType>::j_operator::operator[](int j)
{
  std::size_t const _n = _i * j_array._extentOne + j;
  return j_array.m[_n];
}

template<class DataType>
const typename Array2D<DataType>::j_operator
Array2D<DataType>::operator[](int i) const
{
  return j_operator(*this, i);
}

template<class DataType>
typename Array2D<DataType>::j_operator Array2D<DataType>::operator[](int i)
{
  return j_operator(*this, i);
}

template<class DataType>
Array3D<DataType>::Array3D() :
    _Array_Basic<DataType>(), _extentZero(0), _extentOne(0), _extentTwo(0)
{
}

template<class DataType>
Array3D<DataType>::Array3D(std::size_t const extentZero,
                           std::size_t const extentOne,
                           std::size_t const extentTwo) :
    _Array_Basic<DataType>(extentZero * extentOne * extentTwo),
    _extentZero(extentZero),
    _extentOne(extentOne),
    _extentTwo(extentTwo)
{
}

template<class DataType>
Array3D<DataType>::Array3D(std::size_t const extentZero,
                           std::size_t const extentOne,
                           std::size_t const extentTwo,
                           DataType const value) :
    _Array_Basic<DataType>(extentZero * extentOne * extentTwo, value),
    _extentZero(extentZero),
    _extentOne(extentOne),
    _extentTwo(extentTwo)
{
}

template<class DataType>
Array3D<DataType>::Array3D(std::size_t const extentZero,
                           std::size_t const extentOne,
                           std::size_t const extentTwo,
                           DataType const * array) :
    _Array_Basic<DataType>(extentZero * extentOne * extentTwo, array),
    _extentZero(extentZero),
    _extentOne(extentOne),
    _extentTwo(extentTwo)
{
}

template<class DataType>
Array3D<DataType>::Array3D(Array3D<DataType> const & other) :
    _Array_Basic<DataType>(other),
    _extentZero(other._extentZero),
    _extentOne(other._extentOne),
    _extentTwo(other._extentTwo)
{
}

template<class DataType>
Array3D<DataType>::Array3D(Array3D<DataType> && other) :
    _Array_Basic<DataType>(std::move(other)),
    _extentZero(other._extentZero),
    _extentOne(other._extentOne),
    _extentTwo(other._extentTwo)
{
}

template<class DataType>
Array3D<DataType>::~Array3D()
{
}

template<class DataType>
Array3D<DataType> &
Array3D<DataType>::operator=(Array3D<DataType> const & other)
{
  _Array_Basic<DataType>::operator=(other);
  _extentZero = other._extentZero;
  _extentOne = other._extentOne;
  _extentTwo = other._extentTwo;
  return *this;
}

template<class DataType>
Array3D<DataType> & Array3D<DataType>::operator=(Array3D<DataType> && other)
{
  _Array_Basic<DataType>::operator=(std::move(other));
  _extentZero = other._extentZero;
  _extentOne = other._extentOne;
  _extentTwo = other._extentTwo;
  return *this;
}

template<class DataType>
inline Array2DView<DataType> Array3D<DataType>::data_2D(int i)
{
  return Array2DView<DataType>(
      _extentOne, _extentTwo, this->m.data() + i * _extentOne * _extentTwo);
}

template<class DataType>
inline Array1DView<DataType> Array3D<DataType>::data_1D(int i, int j)
{
  return Array1DView<DataType>(
      _extentTwo, this->m.data() + (i * _extentOne + j) * _extentTwo);
}

template<class DataType>
inline void Array3D<DataType>::resize(int const extentZero,
                                      int const extentOne,
                                      int const extentTwo)
{
  _extentZero = extentZero;
  _extentOne = extentOne;
  _extentTwo = extentTwo;
  std::size_t const _n = _extentZero * _extentOne * _extentTwo;
  this->m.resize(_n, static_cast<DataType>(0));
}

template<class DataType>
inline void Array3D<DataType>::resize(int const extentZero,
                                      int const extentOne,
                                      int const extentTwo,
                                      DataType const new_value)
{
  _extentZero = extentZero;
  _extentOne = extentOne;
  _extentTwo = extentTwo;
  std::size_t const _n = _extentZero * _extentOne * _extentTwo;
  this->m.resize(_n, new_value);
}

template<class DataType>
inline void Array3D<DataType>::resize(int const extentZero,
                                      int const extentOne,
                                      int const extentTwo,
                                      DataType const * new_array)
{
  _extentZero = extentZero;
  _extentOne = extentOne;
  _extentTwo = extentTwo;
  std::size_t const _n = _extentZero * _extentOne * _extentTwo;
  this->m.resize(_n);
  std::copy(new_array, new_array + _n, this->m.data());
}

template<class DataType>
inline const DataType Array3D<DataType>::operator()(int i, int j, int k) const
{
  std::size_t const _n = (i * _extentOne + j) * _extentTwo + k;
  return this->m[_n];
}

template<class DataType>
inline DataType & Array3D<DataType>::operator()(int i, int j, int k)
{
  std::size_t const _n = (i * _extentOne + j) * _extentTwo + k;
  return this->m[_n];
}

template<class DataType>
inline DataType const Array3D<DataType>::at(int i, int j, int k) const
{
  this->_range_check(i, _extentZero);
  this->_range_check(j, _extentOne);
  this->_range_check(k, _extentTwo);
  std::size_t const _n = (i * _extentOne + j) * _extentTwo + k;
  return this->m[_n];
}

template<class DataType>
inline DataType & Array3D<DataType>::at(int i, int j, int k)
{
  this->_range_check(i, _extentZero);
  this->_range_check(j, _extentOne);
  this->_range_check(k, _extentTwo);
  std::size_t const _n = (i * _extentOne + j) * _extentTwo + k;
  return this->m[_n];
}

template<class DataType>
Array3D<DataType>::j_operator::j_operator(Array3D<DataType> & _array, int i) :
    j_array(_array), _i(i)
{
}

template<class DataType>
Array3D<DataType>::j_operator::k_operator::k_operator(
    Array3D<DataType> & _array, int i, int j) :
    k_array(_array), _i(i), _j(j)
{
}

template<class DataType>
const DataType
Array3D<DataType>::j_operator::k_operator::operator[](int k) const
{
  std::size_t const _n
      = (_i * k_array._extentOne + _j) * k_array._extentTwo + k;
  return k_array.m[_n];
}

template<class DataType>
DataType & Array3D<DataType>::j_operator::k_operator::operator[](int k)
{
  std::size_t const _n
      = (_i * k_array._extentOne + _j) * k_array._extentTwo + k;
  return k_array.m[_n];
}

template<class DataType>
const typename Array3D<DataType>::j_operator::k_operator
Array3D<DataType>::j_operator::operator[](int j) const
{
  return k_operator(j_array, _i, j);
}

template<class DataType>
typename Array3D<DataType>::j_operator::k_operator
Array3D<DataType>::j_operator::operator[](int j)
{
  return k_operator(j_array, _i, j);
}

template<class DataType>
const typename Array3D<DataType>::j_operator
Array3D<DataType>::operator[](int i) const
{
  return j_operator(*this, i);
}

template<class DataType>
typename Array3D<DataType>::j_operator Array3D<DataType>::operator[](int i)
{
  return j_operator(*this, i);
}

template<class DataType>
Array4D<DataType>::Array4D() :
    _Array_Basic<DataType>(),
    _extentZero(0),
    _extentOne(0),
    _extentTwo(0),
    _extentThree(0)
{
}

template<class DataType>
Array4D<DataType>::Array4D(std::size_t const extentZero,
                           std::size_t const extentOne,
                           std::size_t const extentTwo,
                           std::size_t const extentThree) :
    _Array_Basic<DataType>(extentZero * extentOne * extentTwo * extentThree),
    _extentZero(extentZero),
    _extentOne(extentOne),
    _extentTwo(extentTwo),
    _extentThree(extentThree)
{
}

template<class DataType>
Array4D<DataType>::Array4D(std::size_t const extentZero,
                           std::size_t const extentOne,
                           std::size_t const extentTwo,
                           std::size_t const extentThree,
                           DataType const value) :
    _Array_Basic<DataType>(extentZero * extentOne * extentTwo * extentThree,
                           value),
    _extentZero(extentZero),
    _extentOne(extentOne),
    _extentTwo(extentTwo),
    _extentThree(extentThree)
{
}

template<class DataType>
Array4D<DataType>::Array4D(std::size_t const extentZero,
                           std::size_t const extentOne,
                           std::size_t const extentTwo,
                           std::size_t const extentThree,
                           DataType const * array) :
    _Array_Basic<DataType>(extentZero * extentOne * extentTwo * extentThree,
                           array),
    _extentZero(extentZero),
    _extentOne(extentOne),
    _extentTwo(extentTwo),
    _extentThree(extentThree)
{
}

template<class DataType>
Array4D<DataType>::Array4D(Array4D<DataType> const & other) :
    _Array_Basic<DataType>(other),
    _extentZero(other._extentZero),
    _extentOne(other._extentOne),
    _extentTwo(other._extentTwo),
    _extentThree(other._extentThree)
{
}

template<class DataType>
Array4D<DataType>::Array4D(Array4D<DataType> && other) :
    _Array_Basic<DataType>(std::move(other)),
    _extentZero(other._extentZero),
    _extentOne(other._extentOne),
    _extentTwo(other._extentTwo),
    _extentThree(other._extentThree)
{
}

template<class DataType>
Array4D<DataType>::~Array4D()
{
}

template<class DataType>
Array4D<DataType> &
Array4D<DataType>::operator=(Array4D<DataType> const & other)
{
  _Array_Basic<DataType>::operator=(other);
  _extentZero = other._extentZero;
  _extentOne = other._extentOne;
  _extentTwo = other._extentTwo;
  _extentThree = other._extentThree;
  return *this;
}

template<class DataType>
Array4D<DataType> & Array4D<DataType>::operator=(Array4D<DataType> && other)
{
  _Array_Basic<DataType>::operator=(std::move(other));
  _extentZero = other._extentZero;
  _extentOne = other._extentOne;
  _extentTwo = other._extentTwo;
  _extentThree = other._extentThree;
  return *this;
}

template<class DataType>
inline Array3DView<DataType> Array4D<DataType>::data_3D(int i)
{
  std::size_t const _n = i * _extentOne * _extentTwo * _extentThree;
  return Array3DView<DataType>(
      _extentOne, _extentTwo, _extentThree, this->m.data() + _n);
}

template<class DataType>
inline Array2DView<DataType> Array4D<DataType>::data_2D(int i, int j)
{
  std::size_t const _n = (i * _extentOne + j) * _extentTwo * _extentThree;
  return Array2DView<DataType>(_extentTwo, _extentThree, this->m.data() + _n);
}

template<class DataType>
inline Array1DView<DataType> Array4D<DataType>::data_1D(int i, int j, int k)
{
  std::size_t const _n = ((i * _extentOne + j) * _extentTwo + k) * _extentThree;
  return Array1DView<DataType>(_extentThree, this->m.data() + _n);
}

template<class DataType>
inline void Array4D<DataType>::resize(int const extentZero,
                                      int const extentOne,
                                      int const extentTwo,
                                      int const extentThree)
{
  _extentZero = extentZero;
  _extentOne = extentOne;
  _extentTwo = extentTwo;
  _extentThree = extentThree;
  std::size_t const _n = _extentZero * _extentOne * _extentTwo * _extentThree;
  this->m.resize(_n, static_cast<DataType>(0));
}

template<class DataType>
inline void Array4D<DataType>::resize(int const extentZero,
                                      int const extentOne,
                                      int const extentTwo,
                                      int const extentThree,
                                      DataType const new_value)
{
  _extentZero = extentZero;
  _extentOne = extentOne;
  _extentTwo = extentTwo;
  _extentThree = extentThree;
  std::size_t const _n = _extentZero * _extentOne * _extentTwo * _extentThree;
  this->m.resize(_n, new_value);
}

template<class DataType>
inline void Array4D<DataType>::resize(int const extentZero,
                                      int const extentOne,
                                      int const extentTwo,
                                      int const extentThree,
                                      DataType const * new_array)
{
  _extentZero = extentZero;
  _extentOne = extentOne;
  _extentTwo = extentTwo;
  _extentThree = extentThree;
  std::size_t const _n = _extentZero * _extentOne * _extentTwo * _extentThree;
  this->m.resize(_n);
  std::copy(new_array, new_array + _n, this->m.data());
}

template<class DataType>
inline const DataType
Array4D<DataType>::operator()(int i, int j, int k, int l) const
{
  std::size_t const _n
      = ((i * _extentOne + j) * _extentTwo + k) * _extentThree + l;
  return this->m[_n];
}

template<class DataType>
inline DataType & Array4D<DataType>::operator()(int i, int j, int k, int l)
{
  std::size_t const _n
      = ((i * _extentOne + j) * _extentTwo + k) * _extentThree + l;
  return this->m[_n];
}

template<class DataType>
inline DataType const Array4D<DataType>::at(int i, int j, int k, int l) const
{
  this->_range_check(i, _extentZero);
  this->_range_check(j, _extentOne);
  this->_range_check(k, _extentTwo);
  this->_range_check(l, _extentThree);
  std::size_t const _n
      = ((i * _extentOne + j) * _extentTwo + k) * _extentThree + l;
  return this->m[_n];
}

template<class DataType>
inline DataType & Array4D<DataType>::at(int i, int j, int k, int l)
{
  this->_range_check(i, _extentZero);
  this->_range_check(j, _extentOne);
  this->_range_check(k, _extentTwo);
  this->_range_check(l, _extentThree);
  std::size_t const _n
      = ((i * _extentOne + j) * _extentTwo + k) * _extentThree + l;
  return this->m[_n];
}

template<class DataType>
Array4D<DataType>::j_operator::j_operator(Array4D<DataType> & _array, int i) :
    j_array(_array), _i(i)
{
}

template<class DataType>
Array4D<DataType>::j_operator::k_operator::k_operator(
    Array4D<DataType> & _array, int i, int j) :
    k_array(_array), _i(i), _j(j)
{
}

template<class DataType>
Array4D<DataType>::j_operator::k_operator::l_operator::l_operator(
    Array4D<DataType> & _array, int i, int j, int k) :
    l_array(_array), _i(i), _j(j), _k(k)
{
}

template<class DataType>
const DataType
Array4D<DataType>::j_operator::k_operator::l_operator::operator[](int l) const
{
  std::size_t const _n
      = ((_i * l_array._extentOne + _j) * l_array._extentTwo + _k)
            * l_array._extentThree
        + l;
  return l_array.m[_n];
}

template<class DataType>
DataType &
Array4D<DataType>::j_operator::k_operator::l_operator::operator[](int l)
{
  std::size_t const _n
      = ((_i * l_array._extentOne + _j) * l_array._extentTwo + _k)
            * l_array._extentThree
        + l;
  return l_array.m[_n];
}

template<class DataType>
const typename Array4D<DataType>::j_operator::k_operator::l_operator
Array4D<DataType>::j_operator::k_operator::operator[](int k) const
{
  return l_operator(k_array, _i, _j, k);
}

template<class DataType>
typename Array4D<DataType>::j_operator::k_operator::l_operator
Array4D<DataType>::j_operator::k_operator::operator[](int k)
{
  return l_operator(k_array, _i, _j, k);
}

template<class DataType>
const typename Array4D<DataType>::j_operator::k_operator
Array4D<DataType>::j_operator::operator[](int j) const
{
  return k_operator(j_array, _i, j);
}

template<class DataType>
typename Array4D<DataType>::j_operator::k_operator
Array4D<DataType>::j_operator::operator[](int j)
{
  return k_operator(j_array, _i, j);
}

template<class DataType>
const typename Array4D<DataType>::j_operator
Array4D<DataType>::operator[](int i) const
{
  return j_operator(*this, i);
}

template<class DataType>
typename Array4D<DataType>::j_operator Array4D<DataType>::operator[](int i)
{
  return j_operator(*this, i);
}

template<class DataType>
Array5D<DataType>::Array5D() :
    _Array_Basic<DataType>(),
    _extentZero(0),
    _extentOne(0),
    _extentTwo(0),
    _extentThree(0),
    _extentFour(0)
{
}

template<class DataType>
Array5D<DataType>::Array5D(std::size_t const extentZero,
                           std::size_t const extentOne,
                           std::size_t const extentTwo,
                           std::size_t const extentThree,
                           std::size_t const extentFour) :
    _Array_Basic<DataType>(extentZero * extentOne * extentTwo * extentThree
                           * extentFour),
    _extentZero(extentZero),
    _extentOne(extentOne),
    _extentTwo(extentTwo),
    _extentThree(extentThree),
    _extentFour(extentFour)
{
}

template<class DataType>
Array5D<DataType>::Array5D(std::size_t const extentZero,
                           std::size_t const extentOne,
                           std::size_t const extentTwo,
                           std::size_t const extentThree,
                           std::size_t const extentFour,
                           DataType const value) :
    _Array_Basic<DataType>(
        extentZero * extentOne * extentTwo * extentThree * extentFour, value),
    _extentZero(extentZero),
    _extentOne(extentOne),
    _extentTwo(extentTwo),
    _extentThree(extentThree),
    _extentFour(extentFour)
{
}

template<class DataType>
Array5D<DataType>::Array5D(std::size_t const extentZero,
                           std::size_t const extentOne,
                           std::size_t const extentTwo,
                           std::size_t const extentThree,
                           std::size_t const extentFour,
                           DataType const * array) :
    _Array_Basic<DataType>(
        extentZero * extentOne * extentTwo * extentThree * extentFour, array),
    _extentZero(extentZero),
    _extentOne(extentOne),
    _extentTwo(extentTwo),
    _extentThree(extentThree),
    _extentFour(extentFour)
{
}

template<class DataType>
Array5D<DataType>::Array5D(Array5D<DataType> const & other) :
    _Array_Basic<DataType>(other),
    _extentZero(other._extentZero),
    _extentOne(other._extentOne),
    _extentTwo(other._extentTwo),
    _extentThree(other._extentThree),
    _extentFour(other._extentFour)
{
}

template<class DataType>
Array5D<DataType>::Array5D(Array5D<DataType> && other) :
    _Array_Basic<DataType>(std::move(other)),
    _extentZero(other._extentZero),
    _extentOne(other._extentOne),
    _extentTwo(other._extentTwo),
    _extentThree(other._extentThree),
    _extentFour(other._extentFour)
{
}

template<class DataType>
Array5D<DataType>::~Array5D()
{
}

template<class DataType>
Array5D<DataType> &
Array5D<DataType>::operator=(Array5D<DataType> const & other)
{
  _Array_Basic<DataType>::operator=(other);
  _extentZero = other._extentZero;
  _extentOne = other._extentOne;
  _extentTwo = other._extentTwo;
  _extentThree = other._extentThree;
  _extentFour = other._extentFour;
  return *this;
}

template<class DataType>
Array5D<DataType> & Array5D<DataType>::operator=(Array5D<DataType> && other)
{
  _Array_Basic<DataType>::operator=(std::move(other));
  _extentZero = other._extentZero;
  _extentOne = other._extentOne;
  _extentTwo = other._extentTwo;
  _extentThree = other._extentThree;
  _extentFour = other._extentFour;
  return *this;
}

template<class DataType>
inline Array4DView<DataType> Array5D<DataType>::data_4D(int i)
{
  std::size_t const _n
      = i * _extentOne * _extentTwo * _extentThree * _extentFour;
  return Array4DView<DataType>(
      _extentOne, _extentTwo, _extentThree, _extentFour, this->m.data() + _n);
}

template<class DataType>
inline Array3DView<DataType> Array5D<DataType>::data_3D(int i, int j)
{
  std::size_t const _n
      = (i * _extentOne + j) * _extentTwo * _extentThree * _extentFour;
  return Array3DView<DataType>(
      _extentTwo, _extentThree, _extentFour, this->m.data() + _n);
}

template<class DataType>
inline Array2DView<DataType> Array5D<DataType>::data_2D(int i, int j, int k)
{
  std::size_t const _n
      = ((i * _extentOne + j) * _extentTwo + k) * _extentThree * _extentFour;
  return Array2DView<DataType>(_extentThree, _extentFour, this->m.data() + _n);
}

template<class DataType>
inline Array1DView<DataType>
Array5D<DataType>::data_1D(int i, int j, int k, int l)
{
  std::size_t const _n
      = (((i * _extentOne + j) * _extentTwo + k) * _extentThree + l)
        * _extentFour;
  return Array1DView<DataType>(_extentFour, this->m.data() + _n);
}

template<class DataType>
inline void Array5D<DataType>::resize(int const extentZero,
                                      int const extentOne,
                                      int const extentTwo,
                                      int const extentThree,
                                      int const extentFour)
{
  _extentZero = extentZero;
  _extentOne = extentOne;
  _extentTwo = extentTwo;
  _extentThree = extentThree;
  _extentFour = extentFour;
  std::size_t const _n
      = _extentZero * _extentOne * _extentTwo * _extentThree * _extentFour;
  this->m.resize(_n, static_cast<DataType>(0));
}

template<class DataType>
inline void Array5D<DataType>::resize(int const extentZero,
                                      int const extentOne,
                                      int const extentTwo,
                                      int const extentThree,
                                      int const extentFour,
                                      DataType const new_value)
{
  _extentZero = extentZero;
  _extentOne = extentOne;
  _extentTwo = extentTwo;
  _extentThree = extentThree;
  _extentFour = extentFour;
  std::size_t const _n
      = _extentZero * _extentOne * _extentTwo * _extentThree * _extentFour;
  this->m.resize(_n, new_value);
}

template<class DataType>
inline void Array5D<DataType>::resize(int const extentZero,
                                      int const extentOne,
                                      int const extentTwo,
                                      int const extentThree,
                                      int const extentFour,
                                      DataType const * new_array)
{
  _extentZero = extentZero;
  _extentOne = extentOne;
  _extentTwo = extentTwo;
  _extentThree = extentThree;
  _extentFour = extentFour;
  std::size_t const _n
      = _extentZero * _extentOne * _extentTwo * _extentThree * _extentFour;
  this->m.resize(_n);
  std::copy(new_array, new_array + _n, this->m.data());
}

template<class DataType>
inline const DataType
Array5D<DataType>::operator()(int i, int j, int k, int l, int n) const
{
  std::size_t const _n
      = (((i * _extentOne + j) * _extentTwo + k) * _extentThree + l)
            * _extentFour
        + n;
  return this->m[_n];
}

template<class DataType>
inline DataType &
Array5D<DataType>::operator()(int i, int j, int k, int l, int n)
{
  std::size_t const _n
      = (((i * _extentOne + j) * _extentTwo + k) * _extentThree + l)
            * _extentFour
        + n;
  return this->m[_n];
}

template<class DataType>
inline DataType const
Array5D<DataType>::at(int i, int j, int k, int l, int n) const
{
  this->_range_check(i, _extentZero);
  this->_range_check(j, _extentOne);
  this->_range_check(k, _extentTwo);
  this->_range_check(l, _extentThree);
  std::size_t const _n
      = (((i * _extentOne + j) * _extentTwo + k) * _extentThree + l)
            * _extentFour
        + n;
  return this->m[_n];
}

template<class DataType>
inline DataType & Array5D<DataType>::at(int i, int j, int k, int l, int n)
{
  this->_range_check(i, _extentZero);
  this->_range_check(j, _extentOne);
  this->_range_check(k, _extentTwo);
  this->_range_check(l, _extentThree);
  this->_range_check(n, _extentFour);
  std::size_t const _n
      = (((i * _extentOne + j) * _extentTwo + k) * _extentThree + l)
            * _extentFour
        + n;
  return this->m[_n];
}

template<class DataType>
Array5D<DataType>::j_operator::j_operator(Array5D<DataType> & _array, int i) :
    j_array(_array), _i(i)
{
}

template<class DataType>
Array5D<DataType>::j_operator::k_operator::k_operator(
    Array5D<DataType> & _array, int i, int j) :
    k_array(_array), _i(i), _j(j)
{
}

template<class DataType>
Array5D<DataType>::j_operator::k_operator::l_operator::l_operator(
    Array5D<DataType> & _array, int i, int j, int k) :
    l_array(_array), _i(i), _j(j), _k(k)
{
}

template<class DataType>
Array5D<DataType>::j_operator::k_operator::l_operator::n_operator::n_operator(
    Array5D<DataType> & _array, int i, int j, int k, int l) :
    n_array(_array), _i(i), _j(j), _k(k), _l(l)
{
}

template<class DataType>
const DataType
Array5D<DataType>::j_operator::k_operator::l_operator::n_operator::operator[](
    int n) const
{
  std::size_t const _n
      = (((_i * n_array._extentOne + _j) * n_array._extentTwo + _k)
             * n_array._extentThree
         + _l)
            * n_array._extentFour
        + n;
  return n_array.m[_n];
}

template<class DataType>
DataType &
Array5D<DataType>::j_operator::k_operator::l_operator::n_operator::operator[](
    int n)
{
  std::size_t const _n
      = (((_i * n_array._extentOne + _j) * n_array._extentTwo + _k)
             * n_array._extentThree
         + _l)
            * n_array._extentFour
        + n;
  return n_array.m[_n];
}

template<class DataType>
const typename Array5D<DataType>::j_operator::k_operator::l_operator::n_operator
Array5D<DataType>::j_operator::k_operator::l_operator::operator[](int l) const
{
  return n_operator(l_array, _i, _j, _k, l);
}

template<class DataType>
typename Array5D<DataType>::j_operator::k_operator::l_operator::n_operator
Array5D<DataType>::j_operator::k_operator::l_operator::operator[](int l)
{
  return n_operator(l_array, _i, _j, _k, l);
}

template<class DataType>
const typename Array5D<DataType>::j_operator::k_operator::l_operator
Array5D<DataType>::j_operator::k_operator::operator[](int k) const
{
  return l_operator(k_array, _i, _j, k);
}

template<class DataType>
typename Array5D<DataType>::j_operator::k_operator::l_operator
Array5D<DataType>::j_operator::k_operator::operator[](int k)
{
  return l_operator(k_array, _i, _j, k);
}

template<class DataType>
const typename Array5D<DataType>::j_operator::k_operator
Array5D<DataType>::j_operator::operator[](int j) const
{
  return k_operator(j_array, _i, j);
}

template<class DataType>
typename Array5D<DataType>::j_operator::k_operator
Array5D<DataType>::j_operator::operator[](int j)
{
  return k_operator(j_array, _i, j);
}

template<class DataType>
const typename Array5D<DataType>::j_operator
Array5D<DataType>::operator[](int i) const
{
  return j_operator(*this, i);
}

template<class DataType>
typename Array5D<DataType>::j_operator Array5D<DataType>::operator[](int i)
{
  return j_operator(*this, i);
}

#define DIM 3

typedef double VectorOfSizeDIM[DIM];

using namespace Descriptor_dr;

class SymmetryFunctions_dr final: public DescriptorKind {
public:
    // final because otherwise I would need "virtual" destructor on DescriptorKind, and that makes enzyme spit out
    // warnings
    // In KLIFF a utility will create in memory file for initiaization?
    // TODO create appropriate constructor
    explicit SymmetryFunctions_dr(std::string &file_name);

    void initFromFile(std::string &file_name);

    SymmetryFunctions_dr() {};
    SymmetryFunctions_dr(std::vector<std::string>* species,
                      std::string* cutoff_function,
                      double * cutoff_matrix,
                      std::vector<std::string>* symmetry_function_types,
                      std::vector<int>* symmetry_function_sizes,
                      std::vector<double>* symmetry_function_parameters);

    void compute(int index,
                 int n_atoms,
                 int *species,
                 int *neigh_list,
                 int number_of_neigh,
                 double *distances,
                 double *zeta) override;

    void clone_empty(DescriptorKind *descriptorKind);

    inline void set_species(std::vector<std::string> &species);

    inline void get_species(std::vector<std::string> &species);

    inline int get_num_species();

    void set_cutoff(char const *name,
                    std::size_t Nspecies,
                    double const *rcut_2D);

    inline double get_cutoff(int iCode, int jCode);

    void add_descriptor(char const *name,
                        double const *values,
                        int row,
                        int col);

    int get_num_descriptors();


private:
    int n_species = -1;

    bool has_three_body_;

    double bhor2ang = 0.529177;
    std::vector<std::string> species_;
    std::vector<int> name_;
    std::vector<int> starting_index_;
    Array2D<double> rcut_2D_;
    std::vector<Array2D<double> > params_;
    std::vector<int> num_param_sets_;
    std::vector<int> num_params_;
};

inline double cut_cos(double r, double rcut);

// Symmetry functions: Jorg Behler, J. Chem. Phys. 134, 074106, 2011.

/*!
 * \brief Radial symmetry function \c g1 suitable for describing the
 * radial environment of atom \c i
 *
 * \c g1 is the sum of the cutoff functions with respect to
 * all neighboring atoms
 *
 * \param r The distance between atoms \c i and \c j
 * \param rcut The cutoff radius
 * \param phi Radial symmetry function \c g1 value
 */
void sym_g1(double r, double rcut, double &phi);

/*!
 * \brief Radial symmetry function \c g2 suitable for describing the
 * radial environment of atom \c i
 *
 * \c g2 is the sum of Gaussians multiplied by cutoff functions. The shifted
 * \g2 is suitable to describe a spherical shell around the reference atom.
 *
 * \param eta Defines the width of the Gaussians
 * \param Rs Shifts the center of the Gaussians to a certain radial distance
 * \param r The distance between atoms \c i and \c j
 * \param rcut The cutoff radius
 * \param phi Radial symmetry function \c g2 value
 */
void sym_g2(double eta, double Rs, double r, double rcut, double &phi);

/*!
 * \brief Radial symmetry function \c g3 suitable for describing the
 * radial environment of atom \c i
 *
 * \c g3 is a damped cosine function with a period length adjusted
 * by parameter \c kappa
 *
 * \param kappa Adjusting the period length of the damped cosine function
 * \param r The distance between atoms \c i and \c j
 * \param rcut The cutoff radius
 * \param phi Radial symmetry function \c g3 value
 *
 * \note
 * Care must be taken when using the \c g3 function. Neighboring atoms can
 * cancel each other’s due to the existence of positive and negative function
 * values. It is recommended to use \c g3 in combination with other symmetry
 * functions.
 */
void sym_g3(double kappa, double r, double rcut, double &phi);

/*!
 * \brief Summations of cosine functions of the angles centered at atom \c i.
 *
 * The angular part must be symmetric with respect to \c 180 angle
 *
 * \param zeta Provides the angular resolution. High values yield a narrower
 * range of nonzero symmetry function values \param lambda Shifting the maxima
 * of the cosine function (can have the values +1 or -1) \param eta
 * Controlling the radial resolution of the radial function. \param r The
 * distance between atoms \c i and \c j \param rcut The cutoff radius \param
 * phi Function \c g4 value
 */
void sym_g4(double zeta, double lambda, double eta, double const *r, double const *rcut, double &phi);

/*!
 * \brief Summations of cosine functions of the angles centered at atom \c i.
 *
 * The angular part must be symmetric with respect to \c 180 angle.
 * In function \c g5 there is no constraint on the distance between atoms
 * resulting in a larger number of terms in the summation.
 *
 * \param zeta Provides the angular resolution. High values yield a narrower
 * range of nonzero symmetry function values \param lambda Shifting the maxima
 * of the cosine function (can have the values +1 or -1) \param eta
 * Controlling the radial resolution of the radial function. \param r The
 * distance between atoms \c i and \c j \param rcut The cutoff radius \param
 * phi Function \c g4 value
 */
void sym_g5(double zeta, double lambda, double eta, double const *r, double const *rcut, double &phi);
#undef LOG_ERROR


/*! \file file_io_utils.hpp
 * \brief Simple file I/O utilities. For reading driver parameter files.
 *
 */
namespace FileIOUtils
{
    std::ifstream open_file(const std::string& file_name);

    void get_next_data_line(std::ifstream& file, std::string& line);

    void parse_int_params(std::string& line, std::vector<int>& params, int num_params);
    void parse_double_params(std::string& line, std::vector<double>& params, int num_params);
    void parse_string_params(std::string& line, std::vector<std::string>& params, int num_params);
    void parse_bool_params(std::string& line, std::vector<bool>& params, int num_params);

}

/*! \fn std::ifstream open_file(const std::string& file_name)
 * \brief Open a file for reading.
 *
 * \param file_name Name of the file to open.
 * \return std::ifstream object for the file.
 */
inline std::ifstream FileIOUtils::open_file(const std::string &file_name)
{
    std::ifstream file(file_name);
    if (!file.is_open())
    {
        throw std::runtime_error("Could not open file: " + file_name);
    }
    return file;
}

/*! \fn void get_next_data_line(std::ifstream& file, std::string& line)
 * \brief Get next line of data from a file. It ignores comments and blank lines.
 *
 * \param file File to read from.
 * \param line String to store the line in.
 */
inline void FileIOUtils::get_next_data_line(std::ifstream &file, std::string &line)
{
    while (line[0] == '#' || line.empty())
    {
        std::getline(file, line);
    }
}

/*! \fn void parse_int_params(std::string& line, std::vector<int>& params, int num_params)
 * \brief Parse a line of data to extract num_params number of integers from it.
 * It stops the moment it has extracted num_params number of integers and throws
 * runtime error if it does not find num_params number of integers in the line.
 *
 * \param line Line of data to parse.
 * \param params Vector to store the parsed data in.
 * \param num_params Number of integers to parse.
 */
inline void FileIOUtils::parse_int_params(std::string &line, std::vector<int>& params, int num_params)
{
    std::string param;
    std::stringstream ss(line);
    int read_params = 0;

    while(!ss.eof() && read_params < num_params)
    {
        ss >> param;
        try {
            params.push_back(std::stoi(param));
            read_params++;
        }
        catch (std::invalid_argument& e)
        {
            // Move on to next parameter
        }
    }
    if (read_params != num_params)
    {
        throw std::runtime_error("Could not read all int parameters");
    }
}

/*! \fn void parse_double_params(std::string& line, std::vector<double>& params, int num_params)
 * \brief Parse a line of data to extract num_params number of doubles from it.
 * It stops the moment it has extracted num_params number of doubles and throws
 * runtime error if it does not find num_params number of doubles in the line.
 *
 * \param line Line of data to parse.
 * \param params Vector to store the parsed data in.
 * \param num_params Number of doubles to parse.
 */
inline void FileIOUtils::parse_double_params(std::string &line, std::vector<double>& params, int num_params)
{
    std::string param;
    std::stringstream ss(line);
    int read_params = 0;

    while(!ss.eof() && read_params < num_params)
    {
        ss >> param;
        try {
            params.push_back(std::stod(param));
            read_params++;
        }
        catch (std::invalid_argument& e)
        {
            // Move on to next parameter
        }
    }
    if (read_params != num_params)
    {
        throw std::runtime_error("Could not read all double parameters");
    }
}

/*! \fn void parse_string_params(std::string& line, std::vector<std::string>& params, int num_params)
 * \brief Parse a line of data to extract num_params number of strings from it.
 * As every data is a string, data, it just returns the first num_params number of string fragments,
 * including numbers etc. so it does not check if the data is actually a string.
 *
 * \param line Line of data to parse.
 * \param params Vector to store the parsed data in.
 * \param num_params Number of strings to parse.
 */
inline void FileIOUtils::parse_string_params(std::string &line, std::vector<std::string>& params, int num_params)
{
    std::string param;
    std::stringstream ss(line);
    int read_params = 0;

    while(!ss.eof() && read_params < num_params)
    {
        ss >> param;
        params.push_back(param);
        read_params++;
    }
    if (read_params != num_params)
    {
        throw std::runtime_error("Could not read all string parameters");
    }
}

/*! \fn void parse_bool_params(std::string& line, std::vector<bool>& params, int num_params)
 * \brief Parse a line of data to extract num_params number of bools from it.
 * every string "true" is considered a boolean manifestation.
 *
 * \param line Line of data to parse.
 * \param params Vector to store the parsed data in.
 * \param num_params Number of bools to parse.
 */
inline void FileIOUtils::parse_bool_params(std::string &line, std::vector<bool>& params, int num_params)
{
    std::string param;
    std::stringstream ss(line);
    int read_params = 0;

    while(!ss.eof() && read_params < num_params)
    {
        ss >> param;
        if (param == "true" || param == "True" || param == "TRUE")
        {
            params.push_back(true);
            read_params++;
        }
        else if (param == "false" || param == "False" || param == "FALSE")
        {
            params.push_back(false);
            read_params++;
        }
    }
    if (read_params != num_params)
    {
        throw std::runtime_error("Could not read all bool parameters");
    }
}

#endif // __DESCRIPTOR_DR_HPP_