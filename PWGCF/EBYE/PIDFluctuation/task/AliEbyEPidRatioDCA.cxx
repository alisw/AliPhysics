/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: ALICE Offline.                                                 *
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


//=========================================================================//
//             AliEbyE Analysis for Particle Ratio Fluctuation             //
//                   Deepika Rathee  | Satyajit Jena                       //
//                   drathee@cern.ch | sjena@cern.ch                       //
//                  Date: Wed Jul  9 18:38:30 CEST 2014                    //
//          New approch to find particle ratio to reduce memory            //
//                             (Test Only)                                 //
//        Copied from NetParticle Classes
//        Origin: Authors: Jochen Thaeder <jochen@thaeder.de>
//                         Michael Weber <m.weber@cern.ch>
//=========================================================================//

#include "TMath.h"
#include "TAxis.h"

#include "AliESDEvent.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliESDtrackCuts.h"

#include "AliAODEvent.h"
#include "AliAODMCParticle.h"

#include "AliEbyEPidRatioDCA.h"

using namespace std;

ClassImp(AliEbyEPidRatioDCA)

//________________________________________________________________________
AliEbyEPidRatioDCA::AliEbyEPidRatioDCA() :
  AliEbyEPidRatioBase("DCA", "DCA"),
  fESDTrackCutsBkg(NULL),
  fHnDCA(NULL) {
  // Constructor   


  AliLog::SetClassDebugLevel("AliEbyEPidRatioDCA",10);
}

//________________________________________________________________________
AliEbyEPidRatioDCA::~AliEbyEPidRatioDCA() {
  // Destructor
}

//________________________________________________________________________
void AliEbyEPidRatioDCA::Process() {
   Float_t etaRange[2];
  fESDTrackCutsBkg->GetEtaRange(etaRange[0],etaRange[1]);

  Float_t ptRange[2];
  fESDTrackCutsBkg->GetPtRange(ptRange[0],ptRange[1]); 

  // -- Track Loop
  for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {
    
    AliVTrack *track = (fESD) ? static_cast<AliVTrack*>(fESD->GetTrack(idxTrack)) : static_cast<AliVTrack*>(fAOD->GetTrack(idxTrack)); 

    // -- Check if track is accepted for basic parameters
    if (!fHelper->IsTrackAcceptedBasicCharged(track))
      continue;
    
    // -- Check if accepted - ESD
    if (fESD && !fESDTrackCutsBkg->AcceptTrack(dynamic_cast<AliESDtrack*>(track)))
      continue;
    
    // -- Check if accepted - AOD
    if (fAOD){
      AliAODTrack * trackAOD = dynamic_cast<AliAODTrack*>(track);
      
      if (!trackAOD) {
	AliError("Pointer to dynamic_cast<AliAODTrack*>(track) = ZERO");
	continue;
      }
      if (!trackAOD->TestFilterBit(fAODtrackCutBit))
	continue;

      // -- Check if in pT and eta range (is done in ESDTrackCuts for ESDs)
      if(!(track->Pt() > ptRange[0] && track->Pt() <= ptRange[1] && TMath::Abs(track->Eta()) <= etaRange[1]))
	continue;
    }

    Int_t gPdgCode = 0;
    Int_t iPid = 0;
    Double_t pid[3];
    if      (fHelper->IsTrackAcceptedPID(track, pid, (AliPID::kPion)))  {  iPid = 1; gPdgCode = 211;}
    else if (fHelper->IsTrackAcceptedPID(track, pid, (AliPID::kKaon)))  {  iPid = 2; gPdgCode = 321;}
    else if (fHelper->IsTrackAcceptedPID(track, pid, (AliPID::kProton))){  iPid = 3; gPdgCode = 2212;}
    else { iPid = 0; gPdgCode = 0;}

    //  cout << " --- DCA ---- " << iPid << "  " << gPdgCode << endl;
    Double_t yP;
    if (!fHelper->IsTrackAcceptedRapidity(track, yP, iPid))
      continue;
  
    Bool_t isDCArAccepted = fHelper->IsTrackAcceptedDCA(track);

    // -- Check if accepted with thighter DCA cuts
    // ?!?!? How to mimic this in AODs?
    if (fESD && !fESDTrackCuts->AcceptTrack(dynamic_cast<AliESDtrack*>(track)))
      isDCArAccepted = kFALSE;
    
    // -- Check for contamination 
    Int_t contIdx = (fIsMC) ? GetContIdxTrack(TMath::Abs(track->GetLabel()), track->Charge(), gPdgCode) : 1;
    
    // -- Get DCAs (dca_r, dca_z, sigma_xy, sigma_xy_z, sigma_z)
    Float_t dca[2], cov[3]; // 
    if (fESD)
      (dynamic_cast<AliESDtrack*>(track))->GetImpactParameters(dca, cov);
    else 
      dca[0] = 1.; 

    // -- Fill THnSparse 
    
    if(iPid != 0) {   
      Double_t hnDCA[10] = {fCentralityBin,0.,  
			    static_cast<Double_t>(track->Charge()), 
			    track->Eta(), 
			    yP, 
			    track->Phi(), 
			    track->Pt(), 
			    static_cast<Double_t>(contIdx),
			    static_cast<Double_t>(isDCArAccepted), 
			    dca[0]};
      fHnDCA->Fill(hnDCA);
    }      
    
    Double_t hnDCA[10] = {fCentralityBin, 
			  static_cast<Double_t>(iPid), 
			  static_cast<Double_t>(track->Charge()), 
			  track->Eta(), 
			  yP, 
			  track->Phi(), 
			  track->Pt(), 
			  static_cast<Double_t>(contIdx),
			  static_cast<Double_t>(isDCArAccepted), 
			  dca[0]};
      fHnDCA->Fill(hnDCA);


    } // for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {

  return;
}      

//________________________________________________________________________
void AliEbyEPidRatioDCA::CreateHistograms() {
  Int_t    binHnDCA[10] = {AliEbyEPidRatioHelper::fgkfHistNBinsCent,4,
			   AliEbyEPidRatioHelper::fgkfHistNBinsSign,      
			   AliEbyEPidRatioHelper::fgkfHistNBinsEta,       
			   AliEbyEPidRatioHelper::fgkfHistNBinsRap,  
			   AliEbyEPidRatioHelper::fgkfHistNBinsPhi,        
			   AliEbyEPidRatioHelper::fgkfHistNBinsPt,   
			   4,  2,  77};      
  
  Double_t minHnDCA[10] = {AliEbyEPidRatioHelper::fgkfHistRangeCent[0],-0.5, 
			   AliEbyEPidRatioHelper::fgkfHistRangeSign[0], 
			   AliEbyEPidRatioHelper::fgkfHistRangeEta[0], 
			   AliEbyEPidRatioHelper::fgkfHistRangeRap[0],  
			   AliEbyEPidRatioHelper::fgkfHistRangePhi[0], 
			   AliEbyEPidRatioHelper::fgkfHistRangePt[0],   
			   0.5, -0.5, -3.};
  
  Double_t maxHnDCA[10] = {AliEbyEPidRatioHelper::fgkfHistRangeCent[1],3.5,
			   AliEbyEPidRatioHelper::fgkfHistRangeSign[1], 
			   AliEbyEPidRatioHelper::fgkfHistRangeEta[1], 
			   AliEbyEPidRatioHelper::fgkfHistRangeRap[1],  
			   AliEbyEPidRatioHelper::fgkfHistRangePhi[1], 
			   AliEbyEPidRatioHelper::fgkfHistRangePt[1],   
			   4.5, 1.5, 3.};


  fHnDCA = new THnSparseD("hnDCA", "cent:pid:etaRec:yRec:phiRec:ptRec:sign:contPart:contSign:DCArAccepted:DCAr", 10, binHnDCA, minHnDCA, maxHnDCA);
  
  fHnDCA->Sumw2();

  fHnDCA->GetAxis(0)->SetTitle("centrality");                   //  0-5|5-10|10-20|20-30|30-40|40-50|50-60|60-70|70-80|80-90 --> 10 bins
  fHnDCA->GetAxis(1)->SetTitle("N_{ch}|N_{#pi}|N_{K}|N_{p}");   //  0 | 1 | 2 | 3 
  fHnDCA->GetAxis(2)->SetTitle("sign");                         //  -1 | 0 | +1 
  fHnDCA->GetAxis(3)->SetTitle("#eta_{Rec}");                   //  eta  [-0.9, 0.9]
  fHnDCA->GetAxis(4)->SetTitle("#it{y}_{Rec}");                 //  rapidity  [-0.5, 0.5]
  fHnDCA->GetAxis(5)->SetTitle("#varphi_{Rec} (rad)");          //  phi  [ 0. , 2Pi]
  fHnDCA->GetAxis(6)->SetTitle("#it{p}_{T,Rec} (GeV/#it{c})");  //  pT   [ 0.2, 2.6]
  fHnDCA->GetAxis(7)->SetTitle("contPart");                     //  1  primary | 2 missId | 3 from WeakDecay | 4 p from Material
  fHnDCA->GetAxis(8)->SetTitle("DCArAccepted");                 //  0 not accepted | 1 accepted 
  fHnDCA->GetAxis(9)->SetTitle("DCAr");                         //  DCAr [-3, 3]


  fHelper->BinLogAxis(fHnDCA,  6, fESDTrackCuts);
  fHelper->BinLogAxis(fHnDCA,  6, fESDTrackCutsBkg);

  // -- Set binning for DCAr
  Double_t binsDCAr[77] = {-3.,-2.85,-2.7,-2.55,-2.4,-2.25,-2.1,-1.95,-1.8,-1.65,-1.5,-1.35,-1.2,-1.05,-0.9,-0.75,-0.6,-0.45,-0.3,-0.285,-0.27,-0.255,-0.24,-0.225,-0.21,-0.195,-0.18,-0.165,-0.15,-0.135,-0.12,-0.105,-0.09,-0.075,-0.06,-0.045,-0.03,-0.015,0.,0.015,0.03,0.045,0.06,0.075,0.09,0.105,0.12,0.135,0.15,0.165,0.18,0.195,0.21,0.225,0.24,0.255,0.27,0.285,0.3,0.45,0.6,0.75,0.9,1.05,1.2,1.35,1.5,1.65,1.8,1.95,2.1,2.25,2.4,2.55,2.7,2.85,3.};
  fHnDCA->GetAxis(9)->Set(76, binsDCAr);

  // ------------------------------------------------------------------
  
  return;
}

//________________________________________________________________________
Int_t AliEbyEPidRatioDCA::GetContIdxTrack(Int_t label, Int_t sign, Int_t gPdgCode) {
  Int_t contIdx = -1;

  AliVParticle* particle = (fESD) ? fMCEvent->GetTrack(label) : static_cast<AliVParticle*>(fArrayMC->At(label));
  if (!particle)
    return contIdx;

  Bool_t isPhysicalPrimary        = (fESD) ? fStack->IsPhysicalPrimary(label): (static_cast<AliAODMCParticle*>(particle))->IsPhysicalPrimary();
  Bool_t isSecondaryFromWeakDecay = (fESD) ? fStack->IsSecondaryFromWeakDecay(label) : (static_cast<AliAODMCParticle*>(particle))->IsSecondaryFromWeakDecay();
  Bool_t isSecondaryFromMaterial  = (fESD) ? fStack->IsSecondaryFromMaterial(label)  : (static_cast<AliAODMCParticle*>(particle))->IsSecondaryFromMaterial();
  
  if (isPhysicalPrimary) {
    if (gPdgCode == 0) {
      // -- Check if correctly identified 
      if (particle->PdgCode() == (sign*gPdgCode))
	contIdx = 1;   
      // -- MissIdentification
      else 
	contIdx = 2;
    }
    else
      contIdx = 1;   
  }

  // -- Check if secondaries from material or weak decay
  else if(isSecondaryFromWeakDecay)
    contIdx = 3;
  else if (isSecondaryFromMaterial)
    contIdx = 4;
  else
    contIdx = -1;

  return contIdx;
}



