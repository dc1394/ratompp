#include "stdafx.h"
#include "paramdb.h"
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
const char* ParamDb::Get(const char* param) const
{
std::map <std::string, std::string>::const_iterator iter = find(param);

	if(iter != end())
		return iter->second.c_str();
	else
	{
        throw std::invalid_argument((boost::format("Parametr '%s' could not be found.") % param).str());
	}
}

size_t ParamDb::GetSize_t(const char* param) const
{
	return static_cast<size_t>(GetLong(param));
}

double ParamDb::GetDouble(const char* param) const
{
	return ::atof(Get(param));
}

long int ParamDb::GetLong(const char* param) const
{
	return ::atol(Get(param));
}

bool ParamDb::GetBool(const char* param) const
{
	return (::strcmp("Yes", Get(param)) == 0);
}


//
// Reads input parameters form file
//
void ParamDb::ReadParams(void)
{
    FILE* in = OpenFile(NULL, "rt");

	while(true)
	{
        auto const ret = ReadOneParam(in);
		if(std::get<0>(ret) != 0)
			break;
        insert(std::pair <std::string, std::string>(std::get<1>(ret), std::get<2>(ret)));
	}

	fclose(in);
}

//
// Reads one parameter from file.
// Returns "0", if parameter is read.
// Returns "1", if end of file is found.
//
std::tuple<int, std::string, std::string> ParamDb::ReadOneParam(FILE* in) const
{
    char line[500 + 1];

	while (true)
	{
		auto const s = std::fgets(line, 500, in);
		if (s == nullptr || feof(in))
            return std::make_tuple(1, std::string(), std::string());

		// Skip comments and empty lines
		if(!(line[0] == '#' || line[0] == '\n'))
			break;
	}
    
    std::string s(line);
    std::vector<std::string> tokens;
    boost::algorithm::split(tokens, s, boost::is_any_of(" \n"));

	if (tokens.size() != 3)
	{
        throw std::invalid_argument((boost::format("Error during reading file. Line '%s'.") % line).str());
	}

    return std::make_tuple(0, tokens[0], tokens[1]);
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

	std::map <std::string, std::string>::const_iterator i;

	for(i = begin(); i != end(); ++i)
		fprintf(out, "%-30s   %s\n\n", i->first.c_str(), i->second.c_str());

	fclose(out);
}
