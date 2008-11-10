//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSVALIDCELLDATASTRUCT_H
#define ALIHLTPHOSVALIDCELLDATASTRUCT_H 

/***************************************************************************
 * Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved.       *
 *                                                                         *
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project. *
 * Contributors are mentioned in the code where appropriate.               *
 *                                                                         *
 * Permission to use, copy, modify and distribute this software and its    *
 * documentation strictly for non-commercial purposes is hereby granted    *
 * without fee, provided that the above copyright notice appears in all    *
 * copies and that both the copyright notice and this permission notice    *
 * appear in the supporting documentation. The authors make no claims      *
 * about the suitability of this software for any purpose. It is           *
 * provided "as is" without express or implied warranty.                   *
 **************************************************************************/

#include "AliHLTDataTypes.h"
#include "Rtypes.h"



//#define N_DATA_SAMPLES 70
//#define N_DATA_SAMPLES 510
#define N_DATA_SAMPLES 300

struct AliHLTPHOSValidCellDataStruct
{
  Int_t fX;
  Int_t fZ;
  Int_t fGain;
  Float_t fEnergy;
  Float_t fTime;
  Int_t fCrazyness;

  //  Int_t fNSamples;

  // Int_t *fData;

  //  unsigned int  fDataOffset;
  //  Int_t *fData2;
  
  //  Int_t fData[N_DATA_SAMPLES];

};



#endif

