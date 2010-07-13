/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercialf purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

/////////////////////////////////////////////////////
//
//  Monitor V0 for TRD
//
//  Authors:                                          
//  Markus Heide <mheide@uni-muenster.de> 
//////////////////////////////////////////////////////
     
#include "TPDGCode.h"
#include "TFile.h"
#include "TTree.h"
#include "TEventList.h"
#include "TObjArray.h"
#include "TH2.h"
#include "TH2F.h"
#include "TH1I.h"
#include "TList.h"

#include "AliLog.h"
#include "AliESDtrack.h"
#include "AliPID.h"
#include "AliTRDv0Monitor.h"
#include "info/AliTRDv0Info.h"


ClassImp(AliTRDv0Monitor)

//________________________________________________________________________
AliTRDv0Monitor::AliTRDv0Monitor() 
  :AliTRDrecoTask()
   ,fOutput(NULL)
   ,hQualityReductions(NULL)
   ,fV0s(NULL)
   ,fData(NULL)
   ,fInfo(NULL)
   ,fP(-1.) 
{
  //
  // Default constructor
  //
  SetNameTitle("v0Monitor", "Reference V0 Monitor");
}

//________________________________________________________________________
AliTRDv0Monitor::AliTRDv0Monitor(const char *name, const char *title) 
  :AliTRDrecoTask(name, title)
   ,fOutput(NULL)
   ,hQualityReductions(NULL)
   ,fV0s(NULL)
   ,fData(NULL)
   ,fInfo(NULL)
   ,fP(-1.)
{
  //
  // Default constructor
  //
  DefineInput(2, TObjArray::Class()); // v0 list
  DefineInput(3, TObjArray::Class()); // pid info list 
  DefineOutput(1, TList::Class());
  
}


//________________________________________________________________________
AliTRDv0Monitor::~AliTRDv0Monitor() 
{
  if(fOutput) delete fOutput;
}
//________________________________________________________________________
void AliTRDv0Monitor::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

 
 fOutput = new TList;
 fOutput->SetName("V0Monitoring");
 

 const char *samplename[AliPID::kSPECIES] = {"electrons","muons","pions","kaons","protons"};
 const char *decayname[AliTRDv0Info::kNDecays] = {"gamma","K0s","Lambda","AntiLambda"};
 const char *detectorname[kNDets] = {"ITS","TPC","TOF"};

 hQualityReductions = new TH1I(Form("hQualityReductions"),Form("Number of tracks cut out by different quality cut steps"),11,-9,2);
 fOutput->Add(hQualityReductions);

 for(Int_t ipart = 0;ipart < AliPID::kSPECIES; ipart++){
   hCutReductions[ipart] = new TH1I(Form("hCutReductions_%s",samplename[ipart]),Form("Number of tracks cut out by different cut steps for %s",samplename[ipart]),19,-17,2);
   fOutput->Add(hCutReductions[ipart]);
   for(Int_t idetector = 0; idetector < kNDets; idetector++){
     hDetPID[idetector][ipart] = new TH2F(Form("hDetector_%s_%s",detectorname[idetector],samplename[ipart]),Form("%s Likelihood for %s vs. momentum",detectorname[idetector], samplename[ipart]),100,0.,10.,100, 0., 1.);

     fOutput->Add(hDetPID[idetector][ipart]);
   }
   hComPID[ipart] = new TH2F(Form("hComPID_%s",samplename[ipart]),Form("Combined TPC/TOF PID: Likelihood for %s",samplename[ipart]),100,0.,10.,100,0.,1.);

     fOutput->Add(hComPID[ipart]);

     for(Int_t cutstep = 0; cutstep < kNCutSteps; cutstep++){
       hTPCdEdx[ipart][cutstep] = new TH2F(Form("hTPCdEdx_%s_[%d]",samplename[ipart],cutstep),Form("TPC dE/dx for %s [%d]",samplename[ipart],cutstep),100,0.,10.,300,0.,300.);
       
       fOutput->Add(hTPCdEdx[ipart][cutstep]);
     }
 }

 for(Int_t iDecay = 0; iDecay < AliTRDv0Info::kNDecays; iDecay++){   
   for(Int_t cutstep =0; cutstep < kNCutSteps; cutstep++){
     hV0Chi2ndf[iDecay][cutstep] =  new TH2F(Form("hV0Chi2ndf_%s_[%d]",decayname[iDecay],cutstep),Form("Chi2/NDF vs. momentum"),100,0.,10.,500, 0., 500.);
     
     fOutput->Add(hV0Chi2ndf[iDecay][cutstep]);
     
     hPsiPair[iDecay][cutstep] =  new TH2F(Form("hV0PsiPair_%s_[%d]",decayname[iDecay],cutstep),Form("Psi_pair vs. momentum"),100,0.,10.,200, 0., 1.6);
     
     fOutput->Add(hPsiPair[iDecay][cutstep]);

     hPointAngle[iDecay][cutstep] =  new TH2F(Form("hPointAngle_%s_[%d]",decayname[iDecay],cutstep),Form("Pointing Angle vs. momentum"),100,0.,10.,500, 0., 1.6);     
     fOutput->Add(hPointAngle[iDecay][cutstep]);

     hDCA[iDecay][cutstep] =  new TH2F(Form("hDCA_%s_[%d]",decayname[iDecay],cutstep),Form("V0 Daughter DCA vs. momentum"),100,0.,10.,500, 0., 1.);
     
     fOutput->Add(hDCA[iDecay][cutstep]);

     hOpenAngle[iDecay][cutstep] =  new TH2F(Form("hOpenAngle_%s_[%d]",decayname[iDecay],cutstep),Form("Opening Angle vs. momentum"),100,0.,10.,500, 0., 1.6);
     
     fOutput->Add(hOpenAngle[iDecay][cutstep]);

     hRadius[iDecay][cutstep] =  new TH2F(Form("hRadius_%s_[%d]",decayname[iDecay],cutstep),Form("V0 Generation Radius vs. momentum"),100,0.,10.,500, 0., 150.);
     
     fOutput->Add(hRadius[iDecay][cutstep]);

    

   
   }

   hInvMass[iDecay] =  new TH2F(Form("hInvMass_%s",decayname[iDecay]),Form("Invariant Mass vs. momentum"),100,0.,10.,500, 0., 2.);
     
   fOutput->Add(hInvMass[iDecay]); 
 
 } 
 
 /*TH1F *hV0mcPID[AliPID::kSPECIES][AliPID::kSPECIES];
   Int_t nPBins = 200;
   
   
   
   for(Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++){
   for(Int_t iSample = 0; iSample < AliPID::kSPECIES; iSample++){
   
   hV0mcPID[iSample][iSpecies] = new TH1F(Form("hV0mcPID_%s_is_%s",name[iSample],name[iSpecies]),Form("%s contained in %s sample",name[iSpecies],name[iSample]), nPBins, 0.2, 13.);
   }
   }*/
}
//________________________________________________________________________
void AliTRDv0Monitor::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  
  if(!(fTracks = dynamic_cast<TObjArray*>(GetInputData(1)))) return;
  if(!(fV0s    = dynamic_cast<TObjArray*>(GetInputData(2)))) return;
  if(!(fInfo   = dynamic_cast<TObjArray*>(GetInputData(3)))) return;
  
  
  AliTRDtrackInfo     *track = NULL;
  AliTRDv0Info *v0(NULL);

 

  for(Int_t itrk=0; itrk<fTracks->GetEntriesFast(); itrk++){
    track = (AliTRDtrackInfo*)fTracks->UncheckedAt(itrk);
    for(Int_t iv(0); iv<fV0s->GetEntriesFast(); iv++){
      if(!(v0 = (AliTRDv0Info*)fV0s->At(iv))) continue;
      if(!v0->HasTrack(track)) continue;
      ULong_t status = track->GetStatus();
      if(!(status&AliESDtrack::kTRDpid)) continue;
      
      hQualityReductions->Fill(v0->fQuality);//fills integer codes for tracks cut out by track/V0 quality cuts
      
      if(!(v0->fQuality == 1))continue;

      for(Int_t part = 0; part < AliPID::kSPECIES; part++){
	hCutReductions[part]->Fill(v0->GetPID(part,track));//fill in numbers of tracks eliminated by different PID cuts
      }

     
      for(Int_t idecay(0), part(-1); idecay <  AliTRDv0Info::kNDecays; idecay++){//loop over decay types considered for reference data
	
	if(idecay ==  AliTRDv0Info::kLambda){ //protons and pions from Lambda
	  part = AliPID::kProton;
	} else if(idecay == AliTRDv0Info::kAntiLambda) { //antiprotons and pions from Anti-Lambda     
	  part = AliPID::kProton;
	} else if(idecay ==   AliTRDv0Info::kK0s) {//pions from K0s
	  part = AliPID::kPion;
	} else if(idecay ==  AliTRDv0Info::kGamma) {//electrons from conversions
	  part = AliPID::kElectron;
	} 
	
	//fill histograms with track/V0 quality cuts only
	hPsiPair[idecay][0]->Fill(v0->fV0Momentum,v0->fPsiPair);//Angle between daughter momentum plane and plane perpendicular to magnetic field
	hInvMass[idecay]->Fill(v0->fV0Momentum,v0->fInvMass[idecay]);//Invariant mass
	hPointAngle[idecay][0]->Fill(v0->fV0Momentum,v0->fPointingAngle);// = TMath::ACos(esdv0->GetV0CosineOfPointingAngle()); // Cosine of pointing angle
	hOpenAngle[idecay][0]->Fill(v0->fV0Momentum,v0->fOpenAngle);// opening angle between daughters
	hDCA[idecay][0]->Fill(v0->fV0Momentum,v0->fDCA);// Distance of closest approach of daughter tracks	
	hV0Chi2ndf[idecay][0]->Fill(v0->fV0Momentum,v0->fChi2ndf[idecay]);//Kalman Filter Chi2/NDF
	hRadius[idecay][0]->Fill(v0->fV0Momentum,v0->fRadius);//distance of decay/conversion from primary vertex in x-y plane
	
	if(v0->HasTrack(track) == -1){	  
	  hTPCdEdx[part][0]->Fill(v0->fTrackN->P(),v0->fTPCdEdx[AliTRDv0Info::kNeg]);//TPC dE/dx for negative track	 
	}
	else if(v0->HasTrack(track) == 1){
	  hTPCdEdx[part][0]->Fill(v0->fTrackP->P(),v0->fTPCdEdx[AliTRDv0Info::kPos]);//TPC dE/dx for positive track
	}
	//fill histograms after invariant mass cuts
	if((v0->fInvMass[idecay] < v0->fUpInvMass[idecay][0])&&(v0->fInvMass[idecay]> v0->fDownInvMass[idecay])){
	  hV0Chi2ndf[idecay][1]->Fill(v0->fV0Momentum,v0->fChi2ndf[idecay]);
	  hPsiPair[idecay][1]->Fill(v0->fV0Momentum,v0->fPsiPair);
	  hPointAngle[idecay][1]->Fill(v0->fV0Momentum,v0->fPointingAngle);
	  hOpenAngle[idecay][1]->Fill(v0->fV0Momentum,v0->fOpenAngle);
	  hDCA[idecay][1]->Fill(v0->fV0Momentum,v0->fDCA);
	  hRadius[idecay][1]->Fill(v0->fV0Momentum,v0->fRadius);
	  if(v0->HasTrack(track) == -1)
	    hTPCdEdx[part][1]->Fill(v0->fTrackN->P(),v0->fTPCdEdx[AliTRDv0Info::kNeg]);
	  else if(v0->HasTrack(track) == 1)
	    hTPCdEdx[part][1]->Fill(v0->fTrackP->P(),v0->fTPCdEdx[AliTRDv0Info::kPos]);

	}

	//fill histograms after all reference selection cuts
	if(v0->GetPID(part,track)==1){
	  hV0Chi2ndf[idecay][2]->Fill(v0->fV0Momentum,v0->fChi2ndf[idecay]);
	  hPsiPair[idecay][2]->Fill(v0->fV0Momentum,v0->fPsiPair);
	  hPointAngle[idecay][2]->Fill(v0->fV0Momentum,v0->fPointingAngle);
	  hOpenAngle[idecay][2]->Fill(v0->fV0Momentum,v0->fOpenAngle);
	  hDCA[idecay][2]->Fill(v0->fV0Momentum,v0->fDCA);
	  hRadius[idecay][2]->Fill(v0->fV0Momentum,v0->fRadius);
	  if(v0->HasTrack(track) == -1)
	    hTPCdEdx[part][2]->Fill(v0->fTrackN->P(),v0->fTPCdEdx[AliTRDv0Info::kNeg]);
	  else if(v0->HasTrack(track) == 1)
	    hTPCdEdx[part][2]->Fill(v0->fTrackP->P(),v0->fTPCdEdx[AliTRDv0Info::kPos]);

	}
      }
    }
  }
  PostData(1, fOutput);
}
//_______________________________________________________________________

/*Int_t AliTRDv0Monitor::GetPDG(Int_t index)
{//yet to come



}*/
//_______________________________________________________________________

