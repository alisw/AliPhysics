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
#ifndef _THERMINATOR_READPAR_
#define _THERMINATOR_READPAR_
#include <iostream>
#include <string>
#include <vector>
#include <exception>
#include "TString.h"

typedef TString STR;

// Ecxeption values
#define RP_Exception_UnknownException 0
#define RP_Exception_NoSuchParamter   1
#define RP_Exception_NoParFile        2

struct struct_option 
{
  STR keyword;
  STR value;
};

typedef struct struct_option option;
typedef std::vector<option> VOPT;

class ReadPar 
{
 private:
  char *fname;
  VOPT options;
  
 public:
  ReadPar(); // Default constructor
  ReadPar(const char *aFName);
  
  int readFile(const char *aFName) throw(int); 
  int printOptions();
  STR getPar(const char *name) throw(STR);
  
};

#endif
