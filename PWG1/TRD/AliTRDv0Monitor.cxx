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

#include "TObjArray.h"
#include "TH2.h"
#include "TH2F.h"
#include "TH1I.h"
#include "TCanvas.h"

#include "AliLog.h"
#include "AliESDtrack.h"
#include "AliPID.h"
#include "AliTRDv0Monitor.h"
#include "info/AliTRDv0Info.h"
#include "info/AliTRDeventInfo.h"


ClassImp(AliTRDv0Monitor)

//________________________________________________________________________
AliTRDv0Monitor::AliTRDv0Monitor() 
  :AliTRDrecoTask()
  ,fhQualityReductions(NULL)
  ,fV0s(NULL)
  ,fData(NULL)
  ,fInfo(NULL)
  ,fP(-1.) 
{
  //
  // Default constructor
  //
  SetNameTitle("TRDv0Monitor", "V0 Monitor for TRD PID");
}

//________________________________________________________________________
AliTRDv0Monitor::AliTRDv0Monitor(const char *name) 
  :AliTRDrecoTask(name, "V0 Monitor for TRD PID")
  ,fhQualityReductions(NULL)
  ,fV0s(NULL)
  ,fData(NULL)
  ,fInfo(NULL)
  ,fP(-1.)
{
  //
  // Default constructor
  //
  DefineInput(3, TObjArray::Class()); // v0 list
  DefineInput(4, TObjArray::Class()); // pid info list
}



// //____________________________________________________________________
// void AliTRDv0Monitor::MakeSummary(){//makes a summary with potentially nice reference figures
//   //TCanvas *cOut = new TCanvas("v0MonitorSummary1", "Summary 1 for task V0Monitor", 1024, 768);
//   //cOut->cd();
//   //GetRefFigure(4);
//   //cOut->SaveAs("V0MonitorSummary.gif");
// 
//   //cOut = new TCanvas("v0MonitorSummary2","Summary 2 for task V0Monitor", 1024, 768);
//   //cOut->cd();
//   //GetRefFigure(5);
//   //cOut->SaveAs("V0MonitorSummary2.gif");
// }

//________________________________________________________________________
Bool_t AliTRDv0Monitor::GetRefFigure(Int_t /*ifig*/)
{
  //creating reference figures

  AliInfo("Implementation on going ...");
  return kTRUE;
}

//________________________________________________________________________
TObjArray* AliTRDv0Monitor::Histos()
{
  // Create histograms
  // Called once
  if(fContainer) return fContainer;

  fContainer = new TObjArray(kNPlots); fContainer->SetOwner();
  fContainer->SetName("V0Monitoring");
  
  const char *samplename[AliPID::kSPECIES] = {"electrons","muons","pions","kaons","protons"};
  const char *decayname[AliTRDv0Info::kNDecays] = {"gamma","K0s","Lambda","AntiLambda"};
  const char *detectorname[kNDets] = {"ITS","TPC","TOF"};
  
  fhQualityReductions = new TH1I(Form("fhQualityReductions"),Form("Number of tracks cut out by different quality cut steps"),11,-9,2);
  fContainer->Add(fhQualityReductions);
  
  for(Int_t ipart = 0;ipart < AliPID::kSPECIES; ipart++){
    fhCutReductions[ipart] = new TH1I(Form("fhCutReductions_%s",samplename[ipart]),Form("Number of tracks cut out by different cut steps for %s",samplename[ipart]),19,-17,2);
    fContainer->Add(fhCutReductions[ipart]);
    for(Int_t idetector = 0; idetector < kNDets; idetector++){
      fhDetPID[idetector][ipart] = new TH2F(Form("fhDetector_%s_%s",detectorname[idetector],samplename[ipart]),Form("%s Likelihood for %s vs. momentum",detectorname[idetector], samplename[ipart]),100,0.,10.,100, 0., 1.);
  
      fContainer->Add(fhDetPID[idetector][ipart]);
    }
    fhComPID[ipart] = new TH2F(Form("fhComPID_%s",samplename[ipart]),Form("Combined TPC/TOF PID: Likelihood for %s",samplename[ipart]),100,0.,10.,100,0.,1.);
  
    fContainer->Add(fhComPID[ipart]);

    for(Int_t cutstep = 0; cutstep < kNCutSteps; cutstep++){
      fhTPCdEdx[ipart][cutstep] = new TH2F(Form("fhTPCdEdx_%s_[%d]",samplename[ipart],cutstep),Form("TPC dE/dx for %s [%d]",samplename[ipart],cutstep),100,0.,10.,300,0.,300.);
      
      fContainer->Add(fhTPCdEdx[ipart][cutstep]);
    }
  }
  
  for(Int_t iDecay = 0; iDecay < AliTRDv0Info::kNDecays; iDecay++){   
    for(Int_t cutstep =0; cutstep < kNCutSteps; cutstep++){
      fhV0Chi2ndf[iDecay][cutstep] =  new TH2F(Form("fhV0Chi2ndf_%s_[%d]",decayname[iDecay],cutstep),Form("Chi2/NDF vs. momentum"),100,0.,10.,500, 0., 500.);
      
      fContainer->Add(fhV0Chi2ndf[iDecay][cutstep]);
      
      fhPsiPair[iDecay][cutstep] =  new TH2F(Form("fhV0PsiPair_%s_[%d]",decayname[iDecay],cutstep),Form("Psi_pair vs. momentum"),100,0.,10.,200, 0., 1.6);
      
      fContainer->Add(fhPsiPair[iDecay][cutstep]);
  
      fhPointAngle[iDecay][cutstep] =  new TH2F(Form("fhPointAngle_%s_[%d]",decayname[iDecay],cutstep),Form("Pointing Angle vs. momentum"),100,0.,10.,500, 0., 1.6);     
      fContainer->Add(fhPointAngle[iDecay][cutstep]);
  
      fhDCA[iDecay][cutstep] =  new TH2F(Form("fhDCA_%s_[%d]",decayname[iDecay],cutstep),Form("V0 Daughter DCA vs. momentum"),100,0.,10.,500, 0., 1.);
      
      fContainer->Add(fhDCA[iDecay][cutstep]);
  
      fhOpenAngle[iDecay][cutstep] =  new TH2F(Form("fhOpenAngle_%s_[%d]",decayname[iDecay],cutstep),Form("Opening Angle vs. momentum"),100,0.,10.,500, 0., 1.6);
      
      fContainer->Add(fhOpenAngle[iDecay][cutstep]);
  
      fhRadius[iDecay][cutstep] =  new TH2F(Form("fhRadius_%s_[%d]",decayname[iDecay],cutstep),Form("V0 Generation Radius vs. momentum"),100,0.,10.,500, 0., 150.);
      
      fContainer->Add(fhRadius[iDecay][cutstep]);
    }
  
    fhInvMass[iDecay] =  new TH2F(Form("fhInvMass_%s",decayname[iDecay]),Form("Invariant Mass vs. momentum"),100,0.,10.,500, 0., 2.);
      
    fContainer->Add(fhInvMass[iDecay]); 
  } 

/*TH1F *hV0mcPID[AliPID::kSPECIES][AliPID::kSPECIES];
  Int_t nPBins = 200;
  
  
  
  for(Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++){
  for(Int_t iSample = 0; iSample < AliPID::kSPECIES; iSample++){
  
  fhV0mcPID[iSample][iSpecies] = new TH1F(Form("fhV0mcPID_%s_is_%s",name[iSample],name[iSpecies]),Form("%s contained in %s sample",name[iSpecies],name[iSample]), nPBins, 0.2, 13.);
  }
  }*/

  return fContainer;
}


//________________________________________________________________________
void AliTRDv0Monitor::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  if(!(fTracks = dynamic_cast<TObjArray*>(GetInputData(1)))) return;
  if(!(fEvent  = dynamic_cast<AliTRDeventInfo*>(GetInputData(2)))) return;
  if(!(fV0s    = dynamic_cast<TObjArray*>(GetInputData(3)))) return;
  if(!(fInfo   = dynamic_cast<TObjArray*>(GetInputData(4)))) return;
  
  
  AliTRDtrackInfo     *track = NULL;
  AliTRDv0Info *v0(NULL);

  for(Int_t itrk=0; itrk<fTracks->GetEntriesFast(); itrk++){
    track = (AliTRDtrackInfo*)fTracks->UncheckedAt(itrk);
    for(Int_t iv(0); iv<fV0s->GetEntriesFast(); iv++){
      if(!(v0 = (AliTRDv0Info*)fV0s->At(iv))) continue;
      if(!v0->HasTrack(track)) continue;
      ULong_t status = track->GetStatus();
      if(!(status&AliESDtrack::kTRDpid)) continue;
      
      fhQualityReductions->Fill(v0->GetQuality());//fills integer codes for tracks cut out by track/V0 quality cuts
      
      if(!(v0->GetQuality() == 1)) continue;

      for(Int_t part = 0; part < AliPID::kSPECIES; part++){
        fhCutReductions[part]->Fill(v0->GetPID(part,track));//fill in numbers of tracks eliminated by different PID cuts
      }

    
      for(Int_t idecay(0), part(-1); idecay <  Int_t(AliTRDv0Info::kNDecays); idecay++){//loop over decay types considered for reference data
        switch(idecay){
        case AliTRDv0Info::kLambda: //protons and pions from Lambda
        case AliTRDv0Info::kAntiLambda: //antiprotons and pions from Anti-Lambda     
          part = AliPID::kProton;
          break;
        case AliTRDv0Info::kK0s: //pions from K0s
          part = AliPID::kPion;
          break;
        case  AliTRDv0Info::kGamma: //electrons from conversions
          part = AliPID::kElectron;
          break;
        }
        
        //fill histograms with track/V0 quality cuts only
        fhPsiPair[idecay][0]->Fill(v0->GetV0Momentum(),v0->GetPsiPair());//Angle between daughter momentum plane and plane perpendicular to magnetic field
        fhInvMass[idecay]->Fill(v0->GetV0Momentum(),v0->GetInvMass(idecay));//Invariant mass
        fhPointAngle[idecay][0]->Fill(v0->GetV0Momentum(),v0->GetPointingAngle());// = TMath::ACos(esdv0->GetV0CosineOfPointingAngle()); // Cosine of pointing angle
        fhOpenAngle[idecay][0]->Fill(v0->GetV0Momentum(),v0->GetOpenAngle());// opening angle between daughters
        fhDCA[idecay][0]->Fill(v0->GetV0Momentum(),v0->GetDCA());// Distance of closest approach of daughter tracks	
        fhV0Chi2ndf[idecay][0]->Fill(v0->GetV0Momentum(),v0->GetChi2ndf(idecay));//Kalman Filter Chi2/NDF
        fhRadius[idecay][0]->Fill(v0->GetV0Momentum(),v0->GetRadius());//distance of decay/conversion from primary vertex in x-y plane
      
        if(v0->HasTrack(track) == -1){	  
          fhTPCdEdx[part][0]->Fill(v0->GetV0Daughter(-1)->P(),v0->GetTPCdEdx(AliTRDv0Info::kNeg));//TPC dE/dx for negative track
        } else if(v0->HasTrack(track) == 1){
          fhTPCdEdx[part][0]->Fill(v0->GetV0Daughter(1)->P(),v0->GetTPCdEdx(AliTRDv0Info::kPos));//TPC dE/dx for positive track
        }
      
        //fill histograms after invariant mass cuts
        if((v0->GetInvMass(idecay) < v0->GetUpInvMass(idecay,0))&&(v0->GetInvMass(idecay)> v0->GetDownInvMass(idecay))){
          fhV0Chi2ndf[idecay][1]->Fill(v0->GetV0Momentum(),v0->GetChi2ndf(idecay));
          fhPsiPair[idecay][1]->Fill(v0->GetV0Momentum(),v0->GetPsiPair());
          fhPointAngle[idecay][1]->Fill(v0->GetV0Momentum(),v0->GetPointingAngle());
          fhOpenAngle[idecay][1]->Fill(v0->GetV0Momentum(),v0->GetOpenAngle());
          fhDCA[idecay][1]->Fill(v0->GetV0Momentum(),v0->GetDCA());
          fhRadius[idecay][1]->Fill(v0->GetV0Momentum(),v0->GetRadius());
          if(v0->HasTrack(track) == -1)
            fhTPCdEdx[part][1]->Fill(v0->GetV0Daughter(-1)->P(),v0->GetTPCdEdx(AliTRDv0Info::kNeg));
          else if(v0->HasTrack(track) == 1)
            fhTPCdEdx[part][1]->Fill(v0->GetV0Daughter(1)->P(),v0->GetTPCdEdx(AliTRDv0Info::kPos));
      
        }

        //fill histograms after all reference selection cuts
        if(v0->GetPID(part,track)==1){
          fhV0Chi2ndf[idecay][2]->Fill(v0->GetV0Momentum(),v0->GetChi2ndf(idecay));
          fhPsiPair[idecay][2]->Fill(v0->GetV0Momentum(),v0->GetPsiPair());
          fhPointAngle[idecay][2]->Fill(v0->GetV0Momentum(),v0->GetPointingAngle());
          fhOpenAngle[idecay][2]->Fill(v0->GetV0Momentum(),v0->GetOpenAngle());
          fhDCA[idecay][2]->Fill(v0->GetV0Momentum(),v0->GetDCA());
          fhRadius[idecay][2]->Fill(v0->GetV0Momentum(),v0->GetRadius());
          if(v0->HasTrack(track) == -1)
            fhTPCdEdx[part][2]->Fill(v0->GetV0Daughter(-1)->P(),v0->GetTPCdEdx(AliTRDv0Info::kNeg));
          else if(v0->HasTrack(track) == 1)
            fhTPCdEdx[part][2]->Fill(v0->GetV0Daughter(1)->P(),v0->GetTPCdEdx(AliTRDv0Info::kPos));
        }
      }
    }
  }
}
