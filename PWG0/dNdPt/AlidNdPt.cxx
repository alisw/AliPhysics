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
 // last change: 2011-04-04 by M.Knichel

#include "AliESDtrackCuts.h"  
#include "AlidNdPtEventCuts.h"
#include "AlidNdPtAcceptanceCuts.h"
#include "AliPhysicsSelection.h"
#include "AlidNdPtBackgroundCuts.h"
#include "AlidNdPt.h"

using namespace std;

ClassImp(AlidNdPt)

//_____________________________________________________________________________
AlidNdPt::AlidNdPt(): TNamed()
, fdNdPtEventCuts(0)
, fdNdPtAcceptanceCuts(0)
, fdNdPtRecAcceptanceCuts(0)
, fEsdTrackCuts(0)
, fUseMCInfo(kFALSE)
, fAnalysisMode(AlidNdPtHelper::kTPC) 
, fTrigger(AliTriggerAnalysis::kMB1) 
, fTriggerClass(0) 
, fParticleMode(AlidNdPtHelper::kAllPart) 
, fPhysicsSelection(0)
, fdNdPtBackgroundCuts(0)
, fAnalyseOutput(kFALSE)
, fMergeTHnSparse(kTRUE)
{
  // default constructor
}

//_____________________________________________________________________________
AlidNdPt::AlidNdPt(Char_t* name, Char_t* title): TNamed(name,title)
, fdNdPtEventCuts(0)
, fdNdPtAcceptanceCuts(0)
, fdNdPtRecAcceptanceCuts(0)
, fEsdTrackCuts(0)
, fUseMCInfo(kFALSE)
, fAnalysisMode(AlidNdPtHelper::kTPC) 
, fTrigger(AliTriggerAnalysis::kMB1) 
, fTriggerClass(0) 
, fParticleMode(AlidNdPtHelper::kAllPart) 
, fPhysicsSelection(0)
, fdNdPtBackgroundCuts(0)
, fAnalyseOutput(kFALSE)
, fMergeTHnSparse(kTRUE)
{
  // constructor
}

AlidNdPt::AlidNdPt(const AlidNdPt&): TNamed()
, fdNdPtEventCuts(0)
, fdNdPtAcceptanceCuts(0)
, fdNdPtRecAcceptanceCuts(0)
, fEsdTrackCuts(0)
, fUseMCInfo(kFALSE)
, fAnalysisMode(AlidNdPtHelper::kTPC) 
, fTrigger(AliTriggerAnalysis::kMB1) 
, fTriggerClass(0) 
, fParticleMode(AlidNdPtHelper::kAllPart) 
, fPhysicsSelection(0)
, fdNdPtBackgroundCuts(0)
, fAnalyseOutput(kFALSE)
, fMergeTHnSparse(kTRUE)
{
  // not implemented
}

AlidNdPt& AlidNdPt::operator=(const AlidNdPt&)
{
  // not implemented
  return *this;
}

//_____________________________________________________________________________
AlidNdPt::~AlidNdPt() {
  // destructor
  if(fdNdPtEventCuts) delete fdNdPtEventCuts; fdNdPtEventCuts=NULL; 
  if(fdNdPtAcceptanceCuts) delete fdNdPtAcceptanceCuts; fdNdPtAcceptanceCuts=NULL;
  if(fdNdPtRecAcceptanceCuts) delete fdNdPtRecAcceptanceCuts; fdNdPtRecAcceptanceCuts=NULL;  
  if(fEsdTrackCuts) delete fEsdTrackCuts; fEsdTrackCuts=NULL;
  if(fPhysicsSelection) delete fPhysicsSelection; fPhysicsSelection=NULL;
  if(fdNdPtBackgroundCuts) delete fdNdPtBackgroundCuts; fdNdPtBackgroundCuts=NULL;
}

//_____________________________________________________________________________
Double_t * AlidNdPt::CreateLogAxis(Int_t nbins, Double_t xmin, Double_t xmax) {
  // retun pointer to the array with log axis
  // it is user responsibility to delete the array
 
  Double_t logxmin = TMath::Log10(xmin);
  Double_t logxmax = TMath::Log10(xmax);
  Double_t binwidth = (logxmax-logxmin)/nbins;
  
  Double_t *xbins =  new Double_t[nbins+1];

  xbins[0] = xmin;
  for (Int_t i=1;i<=nbins;i++) {
    xbins[i] = xmin + TMath::Power(10,logxmin+i*binwidth);
  }

return xbins;
}
//_____________________________________________________________________________

Double_t* AlidNdPt::CloneArray(Int_t n, Double_t* source)
{
    if (!source || n==0) return 0;
    Double_t* dest = new Double_t[n];
    for (Int_t i=0; i<n ; i++) { dest[i] = source[i]; }
    return dest;
}
