#include "stdafx.h"
#include "paramdb.h"
#include <array>                                        // for std::array
#include <vector>                                       // for std::vector
#include <boost/algorithm/string/classification.hpp>    // for boost::is_any_of
#include <boost/algorithm/string/split.hpp>             // for boost::algorithm::split
#include <boost/format.hpp>                             // for boost::format

extern void Intro(FILE* out);

//
// Constructor
//
ParamDb::ParamDb(const char* path) : m_path(path)
{
}

//
// Destructor
//
ParamDb::~ParamDb(void)
{
}

//
// Returns value of parameter "param"
//
std::string ParamDb::Get(std::string const & param) const
{
    auto const iter = find(param);

    if (iter != end()) {
        return iter->second;
    }
    else {
        throw std::invalid_argument((boost::format("Parametr '%s' could not be found.") % param).str());
	}
}

std::size_t ParamDb::GetSize_t(std::string const & param) const
{
	return static_cast<std::size_t>(GetLong(param));
}

double ParamDb::GetDouble(std::string const & param) const
{
	return std::stod(Get(param));
}

long int ParamDb::GetLong(std::string const & param) const
{
	return std::stol(Get(param));
}

bool ParamDb::GetBool(std::string const & param) const
{
	return "Yes" == Get(param);
}


//
// Reads input parameters form file
//
void ParamDb::ReadParams(void)
{
    auto ifs = OpenFile();

	while(true)
	{
        auto const ret = ReadOneParam(ifs);
		if(std::get<0>(ret) != 0)
			break;
        insert(std::pair <std::string, std::string>(std::get<1>(ret), std::get<2>(ret)));
	}
}

//
// Reads one parameter from file.
// Returns "0", if parameter is read.
// Returns "1", if end of file is found.
//
std::tuple<int, std::string, std::string> ParamDb::ReadOneParam(std::ifstream & ifs) const
{
    std::string line;

    while (true)
	{
        std::getline(ifs, line);
		
        if (ifs.eof()) {
            return std::forward_as_tuple(1, std::string(), std::string());
        }

        if (line.empty()) {
            continue;
        }
        // Skip comments and empty lines
        else if (!(line[0] == '#')) {
            break;
        }
	}
    
    std::vector<std::string> tokens;
    boost::algorithm::split(tokens, line, boost::is_any_of(" "));

	if (tokens.size() != 2)
	{
        throw std::invalid_argument((boost::format("Error during reading file. Line '%s'.") % line).str());
	}

    return std::forward_as_tuple(0, tokens[0], tokens[1]);
}

//
// Returns pointer FILE* for opend file.
// "ext" - extenstion of file.
// "mode" - mode of opening.
//
FILE* ParamDb::OpenFile(const char* ext, const char* mode) const
{
std::string tmp(m_path);
FILE* file;

	if(ext)
	{
		tmp += std::string(".");
		tmp += std::string(ext);
	}

	file = fopen(tmp.c_str(), mode);

	if(!file)
	{
        throw std::invalid_argument((boost::format("Cannot open file '%s'.") % tmp).str());
	}
	return file;
}

//
// Returns pointer FILE* for opend file.
// "ext" - extenstion of file.
// "mode" - mode of opening.
//
std::ifstream ParamDb::OpenFile() const
{
    auto ifs = std::ifstream(m_path, std::ios::in);

    if (!ifs)
    {
        throw std::invalid_argument((boost::format("Cannot open file '%s'.") % m_path).str());
    }

    return ifs;
}

//
// Log of read input files
//
void ParamDb::WriteParams(void) const
{
FILE* out = OpenFile("out", "wt");

	::Intro(out);

	time_t tt;
	time(&tt);
	fprintf(out, "  C A L C U L A T I O N   D A T E   %s\n", ctime(&tt));

	fprintf(out, "\n\n"
		"*************************************\n"
		"*                                   *\n"
		"*      I N P U T   P A R A M S      *\n"
		"*                                   *\n"
		"*************************************\n");

	for (auto && i = begin(); i != end(); ++i)
		fprintf(out, "%-30s   %s\n\n", i->first.c_str(), i->second.c_str());

	fclose(out);
}
