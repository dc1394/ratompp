#include "stdafx.h"
#include "paramdb.h"

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
		char str[100];
		sprintf(str, "Parametr '%s' could not be found.", param);
		throw std::invalid_argument(str);
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
char param[50], val[50];

	while(true)
	{
		if(ReadOneParam(in, param, val) != 0)
			break;
		 insert( std::pair <std::string, std::string>(param, val) );
	}

	fclose(in);
}

//
// Reads one parameter from file.
// Returns "0", if parameter is read.
// Returns "1", if end of file is found.
//
int ParamDb::ReadOneParam(FILE* in, char* param, char* val) const
{
char line[500 + 1];
char *s;

	while(1)
	{
		s = fgets(line, 500, in);
		if(s == NULL || feof(in))
			return 1;

		// Skip comments and empty lines
		if(!(line[0] == '#' || line[0] == '\n'))
			break;
	}

	if(sscanf(line, "%s %s", param, val) != 2)
	{
		char str[1000];
		sprintf(str, "Error during reading file. Line '%s'.", line);
		throw std::invalid_argument(str);
	}
	return 0;
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
		char str[1000];
		sprintf(str, "Cannot open file '%s'.", tmp.c_str());
		throw std::invalid_argument(str);
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
