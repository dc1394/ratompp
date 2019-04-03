#ifndef __RATOM_PARAMDB_H__
#define __RATOM_PARAMDB_H__

#include <fstream>  // for std::ifstream
#include <memory>   // for std::unique_ptr
#include <string>   // for std::string
#include <tuple>    // for std::tuple

/** \brief Database of input parameters.
*
* \note It is represented as a hash table for pai (parameter, value).
*   Both, parameter and value are strings. The conversion from string into
*   required value is done after parsing and reading from file
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/


class ParamDb final : private std::map<std::string, std::string>
{
public:
	ParamDb(const char* path);
    ~ParamDb() = default;

	void ReadParams();
	void WriteParams() const;

	std::string Get(std::string const & param) const;
    std::size_t GetSize_t(std::string const & param) const;
    double GetDouble(std::string const & param) const;
    long int GetLong(std::string const & param) const;
    bool GetBool(std::string const & param) const;
    
    std::unique_ptr<FILE, decltype(&std::fclose)> OpenFile(const char* ext, const char* mode) const;
    std::ifstream OpenFile() const;

	void GetPath(std::string& path) const { path = m_path; }

private:
    std::tuple<int, std::string, std::string> ReadOneParam(std::ifstream & ifs) const;

    // Path to input file
	std::string m_path;
};

#endif

