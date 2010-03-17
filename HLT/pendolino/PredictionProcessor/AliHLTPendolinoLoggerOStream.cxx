// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Sebastian Bablok <Sebastian.Bablok@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTPendolinoLoggerOStream.cxx
    @author Sebastian Bablok
    @date   
    @brief  
*/

#include "AliHLTPendolinoLoggerOStream.h"

#include <iostream>

using namespace std;


ClassImp(AliHLTPendolinoLoggerOStream)


AliHLTPendolinoLoggerOStream::AliHLTPendolinoLoggerOStream() {
	// C-tor of AliHLTPendolinoLoggerOStream

}

AliHLTPendolinoLoggerOStream::~AliHLTPendolinoLoggerOStream() {
	// D-tor of AliHLTPendolinoLoggerOStream

}

void AliHLTPendolinoLoggerOStream::log(const char* detector, const char* msg) {
	// logging function, which prints message to command line
	cout << "   === " << detector << ": " << msg << endl;

}

