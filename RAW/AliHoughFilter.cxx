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
// high level filter algorithm for TPC using a hough transformation          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <TStopwatch.h>

#include "AliL3StandardIncludes.h"
#include "AliL3Logging.h"
#include "AliL3Transform.h"
#include "AliL3Hough.h"
#include "AliLog.h"

#include "AliHoughFilter.h"


ClassImp(AliHoughFilter)

//_____________________________________________________________________________
AliHoughFilter::AliHoughFilter()
{
// default constructor

  AliL3Log::fgLevel = AliL3Log::kError;
  if (AliDebugLevel() > 0) AliL3Log::fgLevel = AliL3Log::kWarning;
  if (AliDebugLevel() > 1) AliL3Log::fgLevel = AliL3Log::kInformational;
  if (AliDebugLevel() > 2) AliL3Log::fgLevel = AliL3Log::kDebug;
  if (!AliL3Transform::Init("./", kFALSE)) {
    AliError("HLT initialization failed!");
  }
}

//_____________________________________________________________________________
Bool_t AliHoughFilter::Filter(AliRawEvent* event, AliESD* esd)
{
// TPC hough transformation

  TStopwatch timer;
  timer.Start();

  Float_t ptmin = 0.1*AliL3Transform::GetSolenoidField();

  AliL3Hough *hough1 = new AliL3Hough();
    
  hough1->SetThreshold(4);
  hough1->CalcTransformerParams(ptmin);
  hough1->SetPeakThreshold(70,-1);
  // Attention Z of the vertex to be taken from the event head!
  // So far for debug purposes it is fixed by hand...
  hough1->Init(100,4,event,3.82147);
  hough1->SetAddHistograms();

  AliL3Hough *hough2 = new AliL3Hough();

  hough2->SetThreshold(4);
  hough2->CalcTransformerParams(ptmin);
  hough2->SetPeakThreshold(70,-1);
  hough2->Init(100,4,event,3.82147);
  hough2->SetAddHistograms();

  Int_t nglobaltracks = 0;
  /* In case we run HLT code in 2 threads */
  hough1->StartProcessInThread(0,17);
  hough2->StartProcessInThread(18,35);

  if(hough1->WaitForThreadFinish())
    AliFatal(" Can not join the required thread! ");
  if(hough2->WaitForThreadFinish())
    AliFatal(" Can not join the required thread! ");
  
  /* In case we run HLT code in the main thread
  for(Int_t slice=0; slice<=17; slice++) 
    {
      hough1->ReadData(slice,0);
      hough1->Transform();
      hough1->AddAllHistogramsRows();
      hough1->FindTrackCandidatesRow();
      hough1->AddTracks();
    }
  for(Int_t slice=18; slice<=35; slice++)
    {
      hough2->ReadData(slice,0);
      hough2->Transform();
      hough2->AddAllHistogramsRows();
      hough2->FindTrackCandidatesRow();
      hough2->AddTracks();
    }
  */

  nglobaltracks += hough1->FillESD(esd);
  nglobaltracks += hough2->FillESD(esd);

  /* In case we want to debug the ESD
  gSystem->MakeDirectory("hough1");
  hough1->WriteTracks("./hough1");
  gSystem->MakeDirectory("hough2");
  hough2->WriteTracks("./hough2");
  */

  delete hough1;
  delete hough2;

  printf("Filter has found %d TPC tracks in %f seconds\n", nglobaltracks,timer.RealTime());

  return kFALSE;
}
