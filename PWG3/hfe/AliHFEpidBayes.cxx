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
//
// Class for Bayesian PID
// Implements the abstract base class AliHFEpidBase
// 
// Class contains Bayesian specific cuts and QA histograms
// several detectors added for combined PID
// 
//
// Authors: 
//
//   Yvonne Pachmayer <pachmay@physi.uni-heidelberg.de>
//  
#include <TF1.h>
#include <TMath.h>

#include "AliAODPid.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"
#include "AliLog.h"
#include "AliMCParticle.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"

#include "AliHFEpidTPC.h"
#include "AliHFEpidTOF.h"
#include "AliHFEpidTRD.h"
#include "AliHFEpidBayes.h"
#include "AliHFEpidQAmanager.h"

ClassImp(AliHFEpidBayes)

//___________________________________________________________________
AliHFEpidBayes::AliHFEpidBayes() :
  // add a list here
    AliHFEpidBase()
    , fPIDCombined(NULL)
    , fDetMask(10)
    , fDetMaskDefault(10)
    , fDetMaskDefaultandTRD(14)
    , fpidthres(0.9)
{
  //
  // default  constructor
  // 


}

//___________________________________________________________________
AliHFEpidBayes::AliHFEpidBayes(const char* name) :
  // add a list here
    AliHFEpidBase(name)
    , fPIDCombined(NULL)
    , fDetMask(10)
    , fDetMaskDefault(10)
    , fDetMaskDefaultandTRD(14)
    , fpidthres(0.9)
{
  //
  // default  constructor
  // 
  //

}

//____________________________________________________________
AliHFEpidBayes::AliHFEpidBayes(const AliHFEpidBayes &ref):
    AliHFEpidBase("")
    , fPIDCombined(NULL)
    , fDetMask(10)
    , fDetMaskDefault(10)
    , fDetMaskDefaultandTRD(14)
    , fpidthres(0.9)
{
  //
  // Copy constructor
  //
  ref.Copy(*this);
}

//____________________________________________________________
AliHFEpidBayes &AliHFEpidBayes::operator=(const AliHFEpidBayes &ref){
  //
  // Make assignment
  //
  if(this != &ref){
    ref.Copy(*this);
  } 
  return *this;
}

//___________________________________________________________________
AliHFEpidBayes::~AliHFEpidBayes(){
  //
  // Destructor
  //
}

//___________________________________________________________________
void AliHFEpidBayes::Copy(TObject &ref) const {
  //
  // Performs the copying of the object
  //
  AliHFEpidBayes &target = dynamic_cast<AliHFEpidBayes &>(ref);

  target.fPIDCombined           = fPIDCombined;
  target.fDetMask               = fDetMask;
  target.fDetMaskDefault        = fDetMaskDefault;
  target.fDetMaskDefaultandTRD  = fDetMaskDefaultandTRD;
  target.fpidthres              = fpidthres;

  AliHFEpidBase::Copy(ref);

}

//___________________________________________________________________
Bool_t AliHFEpidBayes::InitializePID(Int_t /*run*/){
  //
  // Initialize Bayes PID
  //

    fPIDCombined=new AliPIDCombined;
    fDetMaskDefault=AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF;
    fDetMaskDefaultandTRD=AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF+AliPIDResponse::kDetTRD;
    return kTRUE;
}

//___________________________________________________________________
Int_t AliHFEpidBayes::IsSelected(const AliHFEpidObject *track, AliHFEpidQAmanager *pidqa) const
{
  //
  // identification of electrons via Bayes pid
  //

    if(pidqa) pidqa->ProcessTrack(track, AliHFEpid::kBAYESpid, AliHFEdetPIDqa::kBeforePID);
    AliDebug(1, "Doing Bayes PID based");

    if (!fkPIDResponse)
    {
	AliFatal("This Task needs the PID response attached to the inputHandler");
	return 0;
    }

    Bool_t used=kFALSE;
    Bool_t usedTOF=kFALSE;
    Bool_t usedTPC=kFALSE;
    Bool_t usedTRD=kFALSE;



    if((fDetMask==fDetMaskDefault)||(fDetMask==fDetMaskDefaultandTRD))
    {
	Double_t probTPC[AliPID::kSPECIES]={0.};
	usedTPC = fkPIDResponse->ComputePIDProbability(AliPIDResponse::kDetTPC,(AliVTrack*)track->GetRecTrack(),AliPID::kSPECIES,probTPC);
	Double_t probTOF[AliPID::kSPECIES]={0.};
	usedTOF = fkPIDResponse->ComputePIDProbability(AliPIDResponse::kDetTOF,(AliVTrack*)track->GetRecTrack(),AliPID::kSPECIES,probTOF);
	if((usedTOF==1)&&(usedTPC==1)) used=1;
        else used=0;

	if(fDetMask==fDetMaskDefaultandTRD)
	{
	    Double_t probTRD[AliPID::kSPECIES]={0.};
	    usedTRD = fkPIDResponse->ComputePIDProbability(AliPIDResponse::kDetTRD,(AliVTrack*)track->GetRecTrack(),AliPID::kSPECIES,probTRD);
	    if((usedTOF==1)&&(usedTPC==1)&&(usedTRD==1)) used=1;
	    else used=0;
	} 
    }




    Double_t fprobComb[AliPID::kSPECIES]={0.};
    CalcCombProb(track,fkPIDResponse, fprobComb);

    Int_t pdg=0;
    if(fprobComb[AliPID::kElectron]>fpidthres) pdg=11;

    if((pidqa) && (pdg==11) && (used==1)) pidqa->ProcessTrack(track, AliHFEpid::kBAYESpid, AliHFEdetPIDqa::kAfterPID);

    return pdg;


}

void AliHFEpidBayes::CalcCombProb(const AliHFEpidObject *track,const AliPIDResponse *fPIDResponse, Double_t* fprobTPCTOF) const
{
  //
  // The name says it all
  //

    for(Int_t i=0;i<AliPID::kSPECIES;i++)
    {
	fprobTPCTOF[i]=0.;
    }
    fPIDCombined->SetEnablePriors(kFALSE);
    fPIDCombined->SetDetectorMask(fDetMask);

    fPIDCombined->ComputeProbabilities((AliVTrack*)track->GetRecTrack(), fPIDResponse, fprobTPCTOF);

}

