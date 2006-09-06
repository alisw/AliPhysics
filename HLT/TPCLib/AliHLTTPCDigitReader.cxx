// $Id$

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Matthias Richter <Matthias.Richter@ift.uib.no>                *
 *          Timm Steinbeck <timm@kip.uni-heidelberg.de>                   *
 *          Jochen Thaeder <thaeder@kip.uni-heidelberg.de>                *
 *          for The ALICE Off-line Project.                               *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTTPCDigitReader.cxx
    @author Timm Steinbeck, Jochen Thaeder, Matthias Richter
    @date   
    @brief  An abstract reader class for TPC data.
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTTPCDigitReader.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTStdIncludes.h"

ClassImp(AliHLTTPCDigitReader)

AliHLTTPCDigitReader::AliHLTTPCDigitReader(){
}

AliHLTTPCDigitReader::~AliHLTTPCDigitReader(){
}

int AliHLTTPCDigitReader::InitBlock(void* ptr,unsigned long size,Int_t firstrow,Int_t lastrow, Int_t patch, Int_t slice){
  if (patch<0 || patch>=AliHLTTPCTransform::GetNumberOfPatches()) {
    HLTError("invalid readout partition number %d", patch);
    return -EINVAL;
  }
  if (firstrow!=AliHLTTPCTransform::GetFirstRow(patch)) {
    HLTWarning("The firstrow parameter does not match the layout of the readout partition %d "
	       "(firstrow=%d). Parameter is ignored", patch, AliHLTTPCTransform::GetFirstRow(patch));
  }
  if (lastrow!=AliHLTTPCTransform::GetLastRow(patch)) {
    HLTWarning("The lastrow parameter does not match the layout of the readout partition %d "
	       "(lastrow=%d). Parameter is ignored", patch, AliHLTTPCTransform::GetLastRow(patch));
  }
  return InitBlock(ptr, size, patch, slice);
}

