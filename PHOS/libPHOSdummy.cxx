/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
$Log$
*/

#include <stdio.h>
/*

    dummy entry points for phos

    */

#ifdef WIN32
#define shlist SHLIST
#define shinit SHINIT
#define shevnt SHEVNT
#define reconsfirst RECONSFIRST
#define type_of_call  _stdcall
#else
#define shlist shlist_
#define shinit shinit_
#define shevnt shevnt_
#define reconsfirst reconsfirst_
#define type_of_call
#endif

#define DUMMY(name) \
extern "C" type_of_call void name() {\
     printf("Dummy version of \"" #name "\" reached \n"); \
     }

DUMMY(shlist)
DUMMY(shinit)
DUMMY(shevnt)
DUMMY(reconsfirst)


