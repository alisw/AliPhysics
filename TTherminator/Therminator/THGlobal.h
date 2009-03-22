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
#ifndef _CF_GLOBAL_H_
#define _CF_GLOBAL_H_

// Define global types

using namespace std;

// Define compilation specific variables

#define PRINT_MESSAGE(_mes) cout << _mes << endl;

#if _DEBUG_LEVEL_==0
#define PRINT_DEBUG_3(_mes) {}
#define PRINT_DEBUG_2(_mes) {}
#define PRINT_DEBUG_1(_mes) {}
#elif _DEBUG_LEVEL_==1
#define PRINT_DEBUG_3(_mes) {}
#define PRINT_DEBUG_2(_mes) {}
#define PRINT_DEBUG_1(_mes) cerr << _mes << endl;
#elif _DEBUG_LEVEL_==2
#define PRINT_DEBUG_3(_mes) {}
#define PRINT_DEBUG_2(_mes) cerr << _mes << endl;
#define PRINT_DEBUG_1(_mes) cerr << _mes << endl;
#elif _DEBUG_LEVEL_==3
#define PRINT_DEBUG_3(_mes) cerr << _mes << endl;
#define PRINT_DEBUG_2(_mes) cerr << _mes << endl;
#define PRINT_DEBUG_1(_mes) cerr << _mes << endl;
#endif

#ifdef _GCC2_
#define STDIOS ios
#endif

#ifdef _GCC3_
#define STDIOS ios_base
#endif

#endif
