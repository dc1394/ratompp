#ifndef __RATOM_PARAMDB_H__
#define __RATOM_PARAMDB_H__


/** \brief Database of input parameters.
*
* \note It is represented as a hash table for pai (parameter, value).
*   Both, parameter and value are strings. The conversion from string into
*   required value is done after parsing and reading from file
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/


class ParamDb : private std::map<std::string, std::string>
{
public:
	ParamDb(const char* path);
	~ParamDb(void);

	void ReadParams();
	void WriteParams() const;

	const char* Get(const char* param) const;
	size_t GetSize_t(const char* param) const;
	double GetDouble(const char* param) const;
	long int GetLong(const char* param) const;
	bool GetBool(const char* param) const;


	FILE* OpenFile(const char* ext, const char* mode) const;

	void GetPath(std::string& path) const { path = m_path; }

private:
	int ReadOneParam(FILE* in, char* param, char* val) const;

private:
        // Path to input file
	std::string m_path;
};

#endif

