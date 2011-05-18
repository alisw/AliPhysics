/******************************************************************************
 *                      T H E R M I N A T O R                                 *
 *                   THERMal heavy-IoN generATOR                              *
 *                           version 1.0                                      *
 *                                                                            *
 * Authors of the model: Wojciech Broniowski, Wojciech.Broniowski@ifj.edu.pl, *
 *                       Wojciech Florkowski, Wojciech.Florkowski@ifj.edu.pl  *
 * Authors of the code:  Adam Kisiel, kisiel@if.pw.edu.pl                     *
 *                       Tomasz Taluc, ttaluc@if.pw.edu.pl                    *
 * Code designers: Adam Kisiel, Tomasz Taluc, Wojciech Broniowski,            *
 *                 Wojciech Florkowski                                        *
 *                                                                            *
 * For the detailed description of the program and furhter references         * 
 * to the description of the model plesase refer to: nucl-th/0504047,         *
 * accessibile at: http://www.arxiv.org/nucl-th/0504047                       *
 *                                                                            *
 * Homepage: http://hirg.if.pw.edu.pl/en/therminator/                         *
 *                                                                            *
 * This code can be freely used and redistributed. However if you decide to   *
 * make modifications to the code, please contact the authors, especially     *
 * if you plan to publish the results obtained with such modified code.       *
 * Any publication of results obtained using this code must include the       *
 * reference to nucl-th/0504047 and the published version of it, when         *
 * available.                                                                 *
 *                                                                            *
 *****************************************************************************/
#include "THGlobal.h"
#include "ReadPar.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iosfwd>
#include <stdlib.h>

ReadPar::ReadPar()
{
  fname = 0;
}

ReadPar::ReadPar(const char *aFName)
{
  fname = strdup(aFName);
  readFile(aFName);
}

ReadPar::~ReadPar()
{
  if (fname)
    free(fname);
}

int ReadPar::readFile(const char *aFName) throw (int)
{
  option read_opt;
  STR buf_str;
  char buff[200];
  char dummy[200];

  std::ifstream       infile(aFName);
  std::istringstream  *instream;
  
  if (!infile.is_open())
    throw RP_Exception_NoParFile;

  while (!infile.eof())
    {
      infile.getline(buff, 200);
      instream = new std::istringstream(buff);
      memset(dummy,0,200);
      *instream >> dummy;
      
      PRINT_DEBUG_3("Read " << dummy);;
      read_opt.keyword = dummy;
      memset(dummy,0,200);
      *instream >> dummy;

      PRINT_DEBUG_3("Read " << dummy);
      if (strstr(dummy,"="))
	{
	  dummy[0]='\0';
	  
	  memset(dummy,0,200);
	  *instream >> dummy;
	  PRINT_DEBUG_3("Read " << dummy);
	  
	  read_opt.value = dummy;
	  options.push_back(read_opt);
	}
    }
  infile.close();

  return 0; 
}

int ReadPar::printOptions()
{
  VOPT::iterator c;

  for (c=options.begin(); c != options.end(); c++)
    PRINT_DEBUG_3("Keyword: " << c->keyword << " Value: " << c->value);

  return 0;
}

STR ReadPar::getPar(const char *name) throw(STR)
{
  VOPT::iterator c;
  STR pname(name);

  for (c=options.begin(); c != options.end(); c++)
    if (c->keyword == pname)
      {
	PRINT_DEBUG_2("Returning value " << c->value << " for keyword " << c->keyword);
	return c->value;
      }
  throw *(new STR(name));

  //  return TString("");
}



