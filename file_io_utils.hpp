//
// Simple file I/O utilities. For reading driver parameter files.
//
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <iostream>


#ifndef LIBDESCRIPTOR_FILE_IO_UTILS_HPP
#define LIBDESCRIPTOR_FILE_IO_UTILS_HPP

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
#endif //LIBDESCRIPTOR_FILE_IO_UTILS_HPP
