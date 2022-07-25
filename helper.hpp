#ifndef KLIFF_HELPER_HPP_
#define KLIFF_HELPER_HPP_

#include <algorithm>
#include <iostream>
#include <sstream>
#include <vector>

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
FormatMessageFileLineFunctionMessage(std::string const &message1,
                                     std::string const &fileName,
                                     long lineNumber,
                                     std::string const &functionName,
                                     std::string const &message2);

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


template<class DataType = double>
class _Array_Basic {
public:
    _Array_Basic();

    _Array_Basic(std::size_t const count);

    _Array_Basic(std::size_t const count, DataType const value);

    _Array_Basic(std::size_t const count, DataType const *array);

    _Array_Basic(_Array_Basic<DataType> const &other);

    _Array_Basic(_Array_Basic<DataType> &&other);

    ~_Array_Basic();

    _Array_Basic<DataType> &operator=(_Array_Basic<DataType> const &other);

    _Array_Basic<DataType> &operator=(_Array_Basic<DataType> &&other);

    inline DataType const *data() const noexcept;

    inline DataType *data() noexcept;

    inline std::size_t size() const;

    inline void clear() noexcept;

    inline void shrink_to_fit();

    inline std::size_t capacity() const noexcept;

    inline void push_back(DataType const &value);

    inline void push_back(DataType &&value);

    std::vector<DataType> m;

protected:
    inline void _range_check(int _n) const;

    inline void _range_check(int _n, std::size_t tsize) const;

};

template<class DataType = double>
class Array1DView {
public:
    Array1DView(std::size_t const count, DataType *array);

    Array1DView(std::size_t const count, DataType const *array);

    Array1DView(Array1DView<DataType> const &other);

    ~Array1DView();

    inline DataType const *data() const noexcept;

    inline DataType *data() noexcept;

    inline const DataType operator()(int i) const;

    inline DataType &operator()(int i);

    inline DataType const at(int i) const;

    inline DataType &at(int i);

    const DataType operator[](int i) const;

    DataType &operator[](int i);

private:
    Array1DView() = delete;

    Array1DView<DataType> &operator=(Array1DView<DataType> const &other)
    = delete;

    Array1DView<DataType> &operator=(Array1DView<DataType> &&other) = delete;

protected:
    inline void _range_check(int _n, std::size_t tsize) const;

protected:
    /*! The extent of the container in the 1st mode */
    std::size_t _extentZero;

    /*! Data pointer */
    DataType *const m;
};

template<class DataType = double>
class Array2DView {
public:
    Array2DView(std::size_t const extentZero,
                std::size_t const extentOne,
                DataType *array);

    Array2DView(std::size_t const extentZero,
                std::size_t const extentOne,
                DataType const *array);

    Array2DView(Array2DView<DataType> const &other);

    ~Array2DView();

    inline DataType const *data() const noexcept;

    inline DataType *data() noexcept;

    inline Array1DView<DataType> data_1D(int i);

    inline const DataType operator()(int i, int j) const;

    inline DataType &operator()(int i, int j);

    inline DataType const at(int i, int j) const;

    inline DataType &at(int i, int j);

    class j_operator {
    public:
        j_operator(Array2DView<DataType> &_array, int i);

        const DataType operator[](int j) const;

        DataType &operator[](int j);

    private:
        /*! Refernce to Array2D class */
        Array2DView<DataType> &j_array;

        std::size_t _i;
    };

    const j_operator operator[](int i) const;

    j_operator operator[](int i);

private:
    Array2DView() = delete;

    Array2DView<DataType> &operator=(Array2DView<DataType> const &other)
    = delete;

    Array2DView<DataType> &operator=(Array2DView<DataType> &&other) = delete;

protected:
    inline void _range_check(int _n, std::size_t tsize) const;

protected:
    /*! The extent of the container in the 1st mode */
    std::size_t _extentZero;

    /*! The extent of the container in the 2nd mode */
    std::size_t _extentOne;

    /*! Data pointer */
    DataType *const m;
};

template<class DataType = double>
class Array2D : public _Array_Basic<DataType> {
public:
    Array2D();

    Array2D(std::size_t const extentZero, std::size_t const extentOne);

    Array2D(std::size_t const extentZero,
            std::size_t const extentOne,
            DataType const value);

    Array2D(std::size_t const extentZero,
            std::size_t const extentOne,
            DataType const *array);

    Array2D(Array2D<DataType> const &other);

    Array2D(Array2D<DataType> &&other);

    ~Array2D();

    Array2D<DataType> &operator=(Array2D<DataType> const &other);

    Array2D<DataType> &operator=(Array2D<DataType> &&other);

    inline Array1DView<DataType> data_1D(int i);

    inline void resize(int const extentZero, int const extentOne);

    inline void
    resize(int const extentZero, int const extentOne, DataType const new_value);

    inline void
    resize(int const extentZero, int const extentOne, DataType const *new_array);

    inline const DataType operator()(int i, int j) const;

    inline DataType &operator()(int i, int j);

    inline DataType const at(int i, int j) const;

    inline DataType &at(int i, int j);

    class j_operator {
    public:
        j_operator(Array2D<DataType> &_array, int i);

        const DataType operator[](int j) const;

        DataType &operator[](int j);

    private:
        /*! Refernce to Array2D class */
        Array2D<DataType> &j_array;

        std::size_t _i;
    };

    const j_operator operator[](int i) const;

    j_operator operator[](int i);

protected:
    /*! The extent of the container in the 1st mode */
    std::size_t _extentZero;

    /*! The extent of the container in the 2nd mode */
    std::size_t _extentOne;
};

void getNextDataLine(FILE *const filePtr,
                     char *const nextLine,
                     int const maxSize,
                     int *endOfFileFlag);

int getXdouble(char *linePtr, const int N, double *list);

int getXint(char *linePtr, const int N, int *list);

void lowerCase(char *linePtr);

void lowerCase(std::string &InputLineArg);

// --------------------------- Implementation --------------------------- //

template<class DataType>
_Array_Basic<DataType>::_Array_Basic() {
}

template<class DataType>
_Array_Basic<DataType>::_Array_Basic(std::size_t const count) :
        m(count, static_cast<DataType>(0)) {
}

template<class DataType>
_Array_Basic<DataType>::_Array_Basic(std::size_t const count,
                                     DataType const value) :
        m(count, value) {
}

template<class DataType>
_Array_Basic<DataType>::_Array_Basic(std::size_t const count,
                                     DataType const *array) :
        m(array, array + count) {
}

template<class DataType>
_Array_Basic<DataType>::_Array_Basic(_Array_Basic<DataType> const &other) :
        m(other.m) {
}

template<class DataType>
_Array_Basic<DataType>::_Array_Basic(_Array_Basic<DataType> &&other) :
        m(std::move(other.m)) {
}

template<class DataType>
_Array_Basic<DataType>::~_Array_Basic() {
}

template<class DataType>
_Array_Basic<DataType> &
_Array_Basic<DataType>::operator=(_Array_Basic<DataType> const &other) {
    m.resize(other.size());
    std::copy(other.m.begin(), other.m.end(), m.begin());
    return *this;
}

template<class DataType>
_Array_Basic<DataType> &
_Array_Basic<DataType>::operator=(_Array_Basic<DataType> &&other) {
    m = std::move(other.m);
    return *this;
}

template<class DataType>
inline DataType const *_Array_Basic<DataType>::data() const noexcept {
    return m.data();
}

template<class DataType>
inline DataType *_Array_Basic<DataType>::data() noexcept {
    return m.data();
}

template<class DataType>
inline std::size_t _Array_Basic<DataType>::size() const {
    return m.size();
}

template<class DataType>
inline void _Array_Basic<DataType>::clear() noexcept {
    m.clear();
}

template<class DataType>
inline void _Array_Basic<DataType>::shrink_to_fit() {
    m.shrink_to_fit();
}

template<class DataType>
inline std::size_t _Array_Basic<DataType>::capacity() const noexcept {
    return m.capacity();
}

template<class DataType>
inline void _Array_Basic<DataType>::push_back(DataType const &value) {
    m.push_back(value);
}

template<class DataType>
inline void _Array_Basic<DataType>::push_back(DataType &&value) {
    m.push_back(value);
}

template<class DataType>
inline void _Array_Basic<DataType>::_range_check(int _n) const {
    if (_n >= size()) {
        LOG_ERROR("The input index is out of range! " + std::to_string(_n)
                  + " >= " + std::to_string(size()));
        std::abort();
    }
}

template<class DataType>
inline void _Array_Basic<DataType>::_range_check(int _n,
                                                 std::size_t tsize) const {
    if (_n >= tsize) {
        LOG_ERROR("The input index is out of range! " + std::to_string(_n)
                  + " >= " + std::to_string(tsize));
        std::abort();
    }
}

template<class DataType>
Array1DView<DataType>::Array1DView(std::size_t const count, DataType *array) :
        _extentZero(count), m(array) {
}

template<class DataType>
Array1DView<DataType>::Array1DView(Array1DView<DataType> const &other) :
        _extentZero(other._extentZero), m(other.m) {
}

template<class DataType>
Array1DView<DataType>::~Array1DView() {
}

template<class DataType>
inline DataType const *Array1DView<DataType>::data() const noexcept {
    return m;
}

template<class DataType>
inline DataType *Array1DView<DataType>::data() noexcept {
    return m;
}

template<class DataType>
inline const DataType Array1DView<DataType>::operator()(int i) const {
    return m[i];
}

template<class DataType>
inline DataType &Array1DView<DataType>::operator()(int i) {
    return m[i];
}

template<class DataType>
inline DataType &Array1DView<DataType>::at(int i) {
    _range_check(i, _extentZero);
    return m[i];
}

template<class DataType>
inline DataType const Array1DView<DataType>::at(int i) const {
    _range_check(i, _extentZero);
    return m[i];
}

template<class DataType>
const DataType Array1DView<DataType>::operator[](int i) const {
    return m[i];
}

template<class DataType>
DataType &Array1DView<DataType>::operator[](int i) {
    return m[i];
}

template<class DataType>
inline void Array1DView<DataType>::_range_check(int _n, std::size_t tsize) const {
    if (_n >= tsize) {
        LOG_ERROR("The input index is out of range! " + std::to_string(_n)
                  + " >= " + std::to_string(tsize));
        std::abort();
    }
}

template<class DataType>
Array2DView<DataType>::Array2DView(std::size_t const extentZero,
                                   std::size_t const extentOne,
                                   DataType *array) :
        _extentZero(extentZero), _extentOne(extentOne), m(array) {
}

template<class DataType>
Array2DView<DataType>::Array2DView(std::size_t const extentZero,
                                   std::size_t const extentOne,
                                   DataType const *array) :
        _extentZero(extentZero),
        _extentOne(extentOne),
        m(const_cast<DataType *>(array)) {
}

template<class DataType>
Array2DView<DataType>::Array2DView(Array2DView<DataType> const &other) :
        _extentZero(other._extentZero), _extentOne(other._extentOne), m(other.m) {
}

template<class DataType>
Array2DView<DataType>::~Array2DView() {
}

template<class DataType>
inline DataType const *Array2DView<DataType>::data() const noexcept {
    return m;
}

template<class DataType>
inline DataType *Array2DView<DataType>::data() noexcept {
    return m;
}

template<class DataType>
inline Array1DView<DataType> Array2DView<DataType>::data_1D(int i) {
    return Array1DView<DataType>(_extentOne, m + i * _extentOne);
}

template<class DataType>
inline const DataType Array2DView<DataType>::operator()(int i, int j) const {
    std::size_t const _n = i * _extentOne + j;
    return m[_n];
}

template<class DataType>
inline DataType &Array2DView<DataType>::operator()(int i, int j) {
    std::size_t const _n = i * _extentOne + j;
    return m[_n];
}

template<class DataType>
inline DataType &Array2DView<DataType>::at(int i, int j) {
    _range_check(i, _extentZero);
    _range_check(j, _extentOne);
    std::size_t const _n = i * _extentOne + j;
    return m[_n];
}

template<class DataType>
inline DataType const Array2DView<DataType>::at(int i, int j) const {
    _range_check(i, _extentZero);
    _range_check(j, _extentOne);
    std::size_t const _n = i * _extentOne + j;
    return m[_n];
}

template<class DataType>
Array2DView<DataType>::j_operator::j_operator(Array2DView<DataType> &_array,
                                              int i) :
        j_array(_array), _i(i) {
}

template<class DataType>
const DataType Array2DView<DataType>::j_operator::operator[](int j) const {
    std::size_t const _n = _i * j_array._extentOne + j;
    return j_array.m[_n];
}

template<class DataType>
DataType &Array2DView<DataType>::j_operator::operator[](int j) {
    std::size_t const _n = _i * j_array._extentOne + j;
    return j_array.m[_n];
}

template<class DataType>
const typename Array2DView<DataType>::j_operator
Array2DView<DataType>::operator[](int i) const {
    return j_operator(*this, i);
}

template<class DataType>
typename Array2DView<DataType>::j_operator
Array2DView<DataType>::operator[](int i) {
    return j_operator(*this, i);
}

template<class DataType>
inline void Array2DView<DataType>::_range_check(int _n, std::size_t tsize) const {
    if (_n >= tsize) {
        LOG_ERROR("The input index is out of range! " + std::to_string(_n)
                  + " >= " + std::to_string(tsize));
        std::abort();
    }
}

template<class DataType>
Array2D<DataType>::Array2D() :
        _Array_Basic<DataType>(), _extentZero(0), _extentOne(0) {
}

template<class DataType>
Array2D<DataType>::Array2D(std::size_t const extentZero,
                           std::size_t const extentOne) :
        _Array_Basic<DataType>(extentZero * extentOne),
        _extentZero(extentZero),
        _extentOne(extentOne) {
}

template<class DataType>
Array2D<DataType>::Array2D(std::size_t const extentZero,
                           std::size_t const extentOne,
                           DataType const value) :
        _Array_Basic<DataType>(extentZero * extentOne, value),
        _extentZero(extentZero),
        _extentOne(extentOne) {
}

template<class DataType>
Array2D<DataType>::Array2D(std::size_t const extentZero,
                           std::size_t const extentOne,
                           DataType const *array) :
        _Array_Basic<DataType>(extentZero * extentOne, array),
        _extentZero(extentZero),
        _extentOne(extentOne) {
}

template<class DataType>
Array2D<DataType>::Array2D(Array2D<DataType> const &other) :
        _Array_Basic<DataType>(other),
        _extentZero(other._extentZero),
        _extentOne(other._extentOne) {
}

template<class DataType>
Array2D<DataType>::Array2D(Array2D<DataType> &&other) :
        _Array_Basic<DataType>(std::move(other)),
        _extentZero(other._extentZero),
        _extentOne(other._extentOne) {
}

template<class DataType>
Array2D<DataType>::~Array2D() {
}

template<class DataType>
Array2D<DataType> &
Array2D<DataType>::operator=(Array2D<DataType> const &other) {
    _Array_Basic<DataType>::operator=(other);
    _extentZero = other._extentZero;
    _extentOne = other._extentOne;
    return *this;
}

template<class DataType>
Array2D<DataType> &Array2D<DataType>::operator=(Array2D<DataType> &&other) {
    _Array_Basic<DataType>::operator=(std::move(other));
    _extentZero = other._extentZero;
    _extentOne = other._extentOne;
    return *this;
}

template<class DataType>
inline Array1DView<DataType> Array2D<DataType>::data_1D(int i) {
    return Array1DView<DataType>(_extentOne, this->m.data() + i * _extentOne);
}

template<class DataType>
inline void Array2D<DataType>::resize(int const extentZero, int const extentOne) {
    _extentZero = extentZero;
    _extentOne = extentOne;
    std::size_t const _n = _extentZero * _extentOne;
    this->m.resize(_n, static_cast<DataType>(0));
}

template<class DataType>
inline void Array2D<DataType>::resize(int const extentZero,
                                      int const extentOne,
                                      DataType const new_value) {
    _extentZero = extentZero;
    _extentOne = extentOne;
    std::size_t const _n = _extentZero * _extentOne;
    this->m.resize(_n, new_value);
}

template<class DataType>
inline void Array2D<DataType>::resize(int const extentZero,
                                      int const extentOne,
                                      DataType const *new_array) {
    _extentZero = extentZero;
    _extentOne = extentOne;
    std::size_t const _n = _extentZero * _extentOne;
    this->m.resize(_n);
    std::copy(new_array, new_array + _n, this->m.data());
}

template<class DataType>
inline const DataType Array2D<DataType>::operator()(int i, int j) const {
    std::size_t const _n = i * _extentOne + j;
    return this->m[_n];
}

template<class DataType>
inline DataType &Array2D<DataType>::operator()(int i, int j) {
    std::size_t const _n = i * _extentOne + j;
    return this->m[_n];
}

template<class DataType>
inline DataType &Array2D<DataType>::at(int i, int j) {
    this->_range_check(i, _extentZero);
    this->_range_check(j, _extentOne);
    std::size_t const _n = i * _extentOne + j;
    return this->m[_n];
}

template<class DataType>
inline DataType const Array2D<DataType>::at(int i, int j) const {
    this->_range_check(i, _extentZero);
    this->_range_check(j, _extentOne);
    std::size_t const _n = i * _extentOne + j;
    return this->m[_n];
}

template<class DataType>
Array2D<DataType>::j_operator::j_operator(Array2D<DataType> &_array, int i) :
        j_array(_array), _i(i) {
}

template<class DataType>
const DataType Array2D<DataType>::j_operator::operator[](int j) const {
    std::size_t const _n = _i * j_array._extentOne + j;
    return j_array.m[_n];
}

template<class DataType>
DataType &Array2D<DataType>::j_operator::operator[](int j) {
    std::size_t const _n = _i * j_array._extentOne + j;
    return j_array.m[_n];
}

template<class DataType>
const typename Array2D<DataType>::j_operator
Array2D<DataType>::operator[](int i) const {
    return j_operator(*this, i);
}

template<class DataType>
typename Array2D<DataType>::j_operator Array2D<DataType>::operator[](int i) {
    return j_operator(*this, i);
}

#undef LOG_ERROR

#endif  // KLIFF_HELPER_HPP_
