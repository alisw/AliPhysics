//-*- Mode: C++ -*-
// $Id$

#ifndef COMMONDEFS_H
#define COMMONDEFS_H

/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2006                                       *
 *                                                                        * 
 * Author: Per Thomas Hille perthi@fys.uio.no for the ALICE HLT Project.  *
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                * 
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//#include "PhosConst.h"

#define PHOSCRYSTALS	(PHOS_MODS*PHOS_ROWS*PHOS_COLS)  // Total number of PHOS crystals

//#define unsigned long int PHOS_CHANNELS	(PHOS_GAINS*PHOS_CRYSTALS) // Total number of PHOS channels
//#define unsigned long int MP_MAP_FILE_NAME	"phosmp.map" // Shared memory map file name
//#define unsigned long int MP_MAP_FILE_SIZE	(PHOS_CHANNELS*1024*8) // Shared memory map file size
//#define unsigned long int MP_RESULT_DIR	"mp_result" // Directory to store result to

#define  PHOSCHANNELS	(PHOS_GAINS*PHOS_CRYSTALS) // Total number of PHOS channels
#define  MP_MAP_FILE_NAME	"phosmp.map" // Shared memory map file name
#define  MP_MAP_FILE_SIZE	(PHOS_CHANNELS*1024*8) // Shared memory map file size
#define  MP_RESULT_DIR	"mp_result" // Directory to store result to


////#define unsigned long int TRUE	1 // General purpose definition
////#define FALSE	0 // General purpose definition


#define MIN(x,y) ((x)>(y)?(y):(x))
#define MAX(x,y) ((x)>(y)?(x):(y))
#define ADJUST_TO_HIGHER(size,step) ((((size)+(step)-1)/(step))*(step))
#define ADJUST_TO_LOWER(size,step) (((size)/(step))*(step))

/////////////////////////////////////////////
// Find index of maximum in array of any type
/////////////////////////////////////////////

#define __IMPLEMENT__findIndexOfMax(type)			\
	inline int findIndexOfMax(type *data, int count)	\
	{ int m=0; for(int i=1; i<count; i++) if(data[i]>data[m]) m=i; return m;}
__IMPLEMENT__findIndexOfMax(char)
__IMPLEMENT__findIndexOfMax(unsigned char)
__IMPLEMENT__findIndexOfMax(int)
__IMPLEMENT__findIndexOfMax(unsigned int)
__IMPLEMENT__findIndexOfMax(short int)
__IMPLEMENT__findIndexOfMax(unsigned short int)
__IMPLEMENT__findIndexOfMax(long long)
__IMPLEMENT__findIndexOfMax(unsigned long long)
__IMPLEMENT__findIndexOfMax(float)
__IMPLEMENT__findIndexOfMax(double)
	
/////////////////////////////////////////////
// Find index of minimum in array of any type
/////////////////////////////////////////////
#define __IMPLEMENT__findIndexOfMin(type)			\
	inline int findIndexOfMin(type *data, int count)	\
	{ int m=0; for(int i=1; i<count; i++) if(data[i]<data[m]) m=i; return m;}
__IMPLEMENT__findIndexOfMin(char)
__IMPLEMENT__findIndexOfMin(unsigned char)
__IMPLEMENT__findIndexOfMin(int)
__IMPLEMENT__findIndexOfMin(unsigned int)
__IMPLEMENT__findIndexOfMin(short int)
__IMPLEMENT__findIndexOfMin(unsigned short int)
__IMPLEMENT__findIndexOfMin(long long)
__IMPLEMENT__findIndexOfMin(unsigned long long)
__IMPLEMENT__findIndexOfMin(float)
__IMPLEMENT__findIndexOfMin(double)
	
//////////////////////////////////////
// Fill struct x of any type with zero
//////////////////////////////////////
#define FILL_ZERO(x)	memset(&(x),0,sizeof(x))

//////////////
// END OF FILE
//////////////

#endif


