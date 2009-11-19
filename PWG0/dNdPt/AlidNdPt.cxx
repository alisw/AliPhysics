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

#include "AliESDtrackCuts.h"  
#include "AlidNdPtEventCuts.h"
#include "AlidNdPtAcceptanceCuts.h"
#include "AlidNdPt.h"

using namespace std;

ClassImp(AlidNdPt)

//_____________________________________________________________________________
AlidNdPt::AlidNdPt(): TNamed()
, fdNdPtEventCuts(0)
, fdNdPtAcceptanceCuts(0)
, fEsdTrackCuts(0)
, fUseMCInfo(kFALSE)
, fAnalysisMode(AlidNdPtHelper::kTPC) 
, fTrigger(AliPWG0Helper::kMB1) 
{
  // default constructor
}

//_____________________________________________________________________________
AlidNdPt::AlidNdPt(Char_t* name, Char_t* title): TNamed(name,title)
, fdNdPtEventCuts(0)
, fdNdPtAcceptanceCuts(0)
, fEsdTrackCuts(0)
, fUseMCInfo(kFALSE)
, fAnalysisMode(AlidNdPtHelper::kTPC) 
, fTrigger(AliPWG0Helper::kMB1) 
{
  // constructor
}

//_____________________________________________________________________________
AlidNdPt::~AlidNdPt(){
  // destructor
  if(fdNdPtEventCuts) delete fdNdPtEventCuts; fdNdPtEventCuts=NULL; 
  if(fdNdPtAcceptanceCuts) delete fdNdPtAcceptanceCuts; fdNdPtAcceptanceCuts=NULL;
  if(fEsdTrackCuts) delete fEsdTrackCuts; fEsdTrackCuts=NULL;
}
