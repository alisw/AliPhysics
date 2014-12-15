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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Projection Chamber clusters objects                                //
//
//  Origin: Marian Ivanov , GSI Darmstadt
//                                                                           //
//                                                                           //
//                                                                          //
///////////////////////////////////////////////////////////////////////////////
#include "AliTPCParam.h" 
#include "AliTPCPRF2D.h"

#include "TObjArray.h"
#include "AliSegmentID.h" 
#include "AliSegmentArray.h" 

#include "AliDigits.h"
#include "AliSimDigits.h"
#include "AliDigitsArray.h" 
#include "AliTPCDigitsArray.h"
#include <TDirectory.h>



//_____________________________________________________________________________

ClassImp(AliTPCDigitsArray) 

AliTPCDigitsArray::AliTPCDigitsArray(Bool_t sim)
                  :AliDigitsArray(),
		   fBSim(kFALSE),
		   fCompression(0),
		   fTrackLevel(0)
{
  //
  //default constructor
  fParam = 0;
  fBSim = sim;
  if ( sim == kTRUE) SetClass("AliSimDigits");
  else
    SetClass("AliDigits");
  fParam = 0;
  //  fPRF   = 0;
  //fRF    = 0;  
  fCompression = 1;
  fTrackLevel = 3;
}

AliTPCDigitsArray::~AliTPCDigitsArray()
{
  //
  
  //
}

AliDigits *  AliTPCDigitsArray::CreateRow(Int_t sector, Int_t row)
{
  //
  //create digits row  
  //
  //if row just exist - delete it
  AliTPCParam * param = (AliTPCParam*)fParam;
  Int_t index = param->GetIndex(sector,row);  
  AliDigits * dig = (AliDigits *)(*this)[index];
  if (dig !=0) delete dig;

  dig = (AliDigits *) AddSegment(index);
  if (dig == 0) return 0;
  dig->Allocate(param->GetMaxTBin(),param->GetNPads(sector,row));  
  if (fBSim == kTRUE) ((AliSimDigits*) dig)->AllocateTrack(fTrackLevel);
  return dig;
}


AliDigits * AliTPCDigitsArray::GetRow(Int_t sector,Int_t row)
{
  //
  //return clusters ((AliTPCDigitsRow *) per given sector and padrow
  //
  if (fParam==0) return 0;
  Int_t index = ((AliTPCParam*)fParam)->GetIndex(sector,row);  
  return (AliDigits *)(*this)[index];
}

AliDigits * AliTPCDigitsArray::LoadRow(Int_t sector,Int_t row)
{
  //
  //return clusters ((AliTPCDigitsRow *) per given sector and padrow
  //
  if (fParam==0) return 0;
  Int_t index = ((AliTPCParam*)fParam)->GetIndex(sector,row);  
  return (AliDigits *)LoadSegment(index);
}

Bool_t  AliTPCDigitsArray::StoreRow(Int_t sector,Int_t row)
{
  //
  //return clusters ((AliTPCDigitsRow *) per given sector and padrow
  //
  AliTPCParam * param = (AliTPCParam*)fParam;
  if (fParam==0) return 0;
  Int_t index = param->GetIndex(sector,row);  
  ( (AliDigits *)At(index))->CompresBuffer(fCompression,param->GetZeroSup());
  if (fBSim == kTRUE) ( (AliSimDigits *)At(index))->CompresTrackBuffer(1);
  StoreSegment(index);
  return kTRUE;
}

Bool_t  AliTPCDigitsArray::ClearRow(Int_t sector,Int_t row)
{
  //
  //return clusters ((AliTPCDigitsRow *) per given sector and padrow
  //
  if (fParam==0) return 0;
  Int_t index = ((AliTPCParam*)fParam)->GetIndex(sector,row);  
  ClearSegment(index);
  return kTRUE;
}



Bool_t AliTPCDigitsArray::Setup(AliDetectorParam *param)
{
  //
  //setup  function to adjust array parameters
  //
  if (param==0) return kFALSE;
  if (fParam !=0) delete fParam;
  //  fParam = new AliTPCParam((AliTPCParam&)(*param));
  fParam = param;
  return MakeArray(((AliTPCParam*)fParam)->GetNRowsTotal());
}

Bool_t AliTPCDigitsArray::Update()
{
  //
  //setup  function to adjust array parameters
  //
  if (fParam ==0 ) return kFALSE;
  if (fTree!=0) return MakeDictionary( ((AliTPCParam*)fParam)->GetNRowsTotal()) ;
  ((AliTPCParam*)fParam)->Update();
  return MakeArray(((AliTPCParam*)fParam)->GetNRowsTotal());
}
