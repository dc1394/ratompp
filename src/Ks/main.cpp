/** \brief Main file for RAtom application
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/


#include "stdafx.h"
#include "nonlinks.h"

void Intro(FILE* out);

int main(int argc, char* argv[])
{
	Intro(stdout);

	if(argc != 2)
	{
		printf("Usage: ratom name\n\n");
		return 1;
	}

	try
	{
		ks::NonLinKs ks(argv[1]);
		return ks.Run();
	}
	catch(std::exception& e)
	{
		printf("ERROR! %s\n", e.what());
		return 1;
	}
}

void Intro(FILE* out)
{
	fprintf(out,
		"===============================================================================\n"
		" RRRR   AAA  TTTTT  OOO   M   M      Zbigniew Romanowski          \n"
		" R   R A   A  TTT  O   O  MM MM                                   \n"         
		" RRRR  AAAAA   T   O   O  M M M      VERSION   1.4.1              \n"
		" R  R  A   A   T   O   O  M   M                                   \n"    
		" R   R A   A   T    OOO   M   M      compilation date: %s         \n"
		"===============================================================================\n\n\n", __DATE__);
}

