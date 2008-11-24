//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSCHANNELDATASTRUCT_H
#define ALIHLTPHOSCHANNELDATASTRUCT_H

/**************************************************************************
 * Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved.      *
 *                                                                        *
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*
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
//#include "AliHLTPHOSCommonDefs.h"

#include "Rtypes.h"

//#include "AliHLTPHOSConstants.h"
//using namespace PhosHLTConst;

struct AliHLTPHOSChannelDataStruct
{
  Float_t fEnergy;
  Float_t fTime;
  UShort_t fChannelID;
  Short_t fCrazyness;
};

#endif

