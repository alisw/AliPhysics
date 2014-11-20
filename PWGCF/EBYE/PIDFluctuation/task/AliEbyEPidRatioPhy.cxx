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
#include "TProfile.h" 
#include "TProfile2D.h" 
#include "TH2D.h" 
#include "TH3D.h" 
#include "TRandom3.h"

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"

#include "AliCentrality.h"
#include "AliTracker.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"

#include "AliEbyEPidRatioPhy.h"

using namespace std;

ClassImp(AliEbyEPidRatioPhy)
//________________________________________________________________________
AliEbyEPidRatioPhy::AliEbyEPidRatioPhy() : 
  AliEbyEPidRatioBase("Dist", "Dist"),
  fOutList(NULL),
  fOrder(8),
  fNNp(7),
  fNp(NULL),
  fNpPt(NULL),
  fNMCNp(7),
  fMCNp(NULL),
  fMCNpPt(NULL),
  fRedFactp(NULL),
  fPtBinHist(NULL), 
  fIsQA(kFALSE), 
  fIsSub(kFALSE), 
  fIsBS(kFALSE), 
  fIsPer(kFALSE), fRan(0),
  fHnTrackUnCorrRec(NULL), fHnTrackUnCorrMc(NULL) {
  AliLog::SetClassDebugLevel("AliEbyEPidRatioPhy",10);
}

//________________________________________________________________________
AliEbyEPidRatioPhy::~AliEbyEPidRatioPhy() {
  // Destructor

 
  for (Int_t ii = 0; ii < fNNp; ++ii) {
    for (Int_t kk = 0; kk < 2; ++kk)
      if (fNp[ii][kk]) delete[] fNp[ii][kk];
    if (fNp[ii]) delete[] fNp[ii];
  }
  if (fNp) delete[] fNp;


 for (Int_t ii = 0; ii < fNMCNp; ++ii) {
    for (Int_t kk = 0; kk < 2; ++kk)
      if (fMCNp[ii][kk]) delete[] fMCNp[ii][kk];
    if (fMCNp[ii]) delete[] fMCNp[ii];
  }
  if (fMCNp) delete[] fMCNp;

  for (Int_t ii = 0; ii <= fOrder; ++ii) 
    if (fRedFactp[ii]) delete[] fRedFactp[ii];
  if (fRedFactp) delete[] fRedFactp;

  if(fPtBinHist) delete fPtBinHist;

 for (Int_t ii = 0; ii < fNNp; ++ii) {
    for (Int_t kk = 0; kk < 2; ++kk)
      for (Int_t jj = 0; jj < 2; ++jj) {
	if (fNpPt[ii][kk][jj]) delete[] fNpPt[ii][kk][jj];
	if(fNpPt[ii][kk]) delete[] fNpPt[ii][kk];
      }
    if (fNpPt[ii]) delete[] fNpPt[ii];
  }
  if (fNpPt) delete[] fNpPt;

 for (Int_t ii = 0; ii < fNNp; ++ii) {
    for (Int_t kk = 0; kk < 2; ++kk)
      for (Int_t jj = 0; jj < 2; ++jj) {
	if (fMCNpPt[ii][kk][jj]) delete[] fMCNpPt[ii][kk][jj];
	if(fMCNpPt[ii][kk]) delete[] fMCNpPt[ii][kk];
      }
    if (fMCNpPt[ii]) delete[] fMCNpPt[ii];
  }
  if (fMCNpPt) delete[] fMCNpPt;

  return;
}

//________________________________________________________________________
void AliEbyEPidRatioPhy::Process() {
  ProcessTracks();
  if (fIsMC)
    ProcessParticles();
  return;
}

//________________________________________________________________________
void AliEbyEPidRatioPhy::Init() {
 
  fNp = new Int_t**[fNNp];
  for (Int_t ii = 0 ; ii < fNNp; ++ii) {
    fNp[ii] = new Int_t*[4];
    for (Int_t kk = 0 ; kk < 4; ++kk)
      fNp[ii][kk] = new Int_t[2];
  }

 fNpPt = new Int_t***[fNNp];
  for (Int_t ii = 0 ; ii < fNNp; ++ii) {
    fNpPt[ii] = new Int_t**[4];
    for (Int_t kk = 0 ; kk < 4; ++kk) {
      fNpPt[ii][kk] = new Int_t*[2];
      for (Int_t ll = 0 ; ll < 2; ++ll)
	fNpPt[ii][kk][ll] = new Int_t[AliEbyEPidRatioHelper::fgkfHistNBinsPt];
    }
  }

  fMCNp = new Int_t**[fNMCNp];
  for (Int_t ii = 0 ; ii < fNMCNp; ++ii) {
    fMCNp[ii] = new Int_t*[4];
    for (Int_t kk = 0 ; kk < 4; ++kk)
      fMCNp[ii][kk] = new Int_t[2];
  }
  
 fMCNpPt = new Int_t***[fNMCNp];
  for (Int_t ii = 0 ; ii < fNMCNp; ++ii) {
    fMCNpPt[ii] = new Int_t**[4];
    for (Int_t kk = 0 ; kk < 4; ++kk) {
      fMCNpPt[ii][kk] = new Int_t*[2];
      for (Int_t ll = 0 ; ll < 2; ++ll)
	fMCNpPt[ii][kk][ll] = new Int_t[AliEbyEPidRatioHelper::fgkfHistNBinsPt];
    }
  }

  fRedFactp = new Double_t*[fOrder+1];
  for (Int_t ii = 0 ; ii <= fOrder; ++ii)
    fRedFactp[ii] = new Double_t[2];
  Printf(" >>>> AliEbyEPidRatioEffContExtra - inside");
}

//________________________________________________________________________
void AliEbyEPidRatioPhy::Reset() {
  for (Int_t ii = 0; ii < fNNp; ++ii) 
    for (Int_t kk = 0; kk < 4; ++kk)
      for (Int_t jj = 0; jj < 2; ++jj)
	fNp[ii][kk][jj] = 0;
  for (Int_t ii = 0; ii < fNMCNp; ++ii) 
    for (Int_t kk = 0; kk < 4; ++kk)
      for (Int_t jj = 0; jj < 2; ++jj)
	fMCNp[ii][kk][jj] = 0;
  for (Int_t ii = 0; ii <= fOrder; ++ii) 
    for (Int_t jj = 0; jj < 2; ++jj)
      fRedFactp[ii][jj] = 1.;

  for (Int_t ii = 0; ii < fNNp; ++ii) 
    for (Int_t kk = 0; kk < 4; ++kk)
      for (Int_t jj = 0; jj < 2; ++jj)
	for (Int_t ll = 0; ll < AliEbyEPidRatioHelper::fgkfHistNBinsPt; ++ll)
	  fNpPt[ii][kk][jj][ll] = 0;

 for (Int_t ii = 0; ii < fNMCNp; ++ii) 
    for (Int_t kk = 0; kk < 4; ++kk)
      for (Int_t jj = 0; jj < 2; ++jj)
	for (Int_t ll = 0; ll < AliEbyEPidRatioHelper::fgkfHistNBinsPt; ++ll)
	  fMCNpPt[ii][kk][jj][ll] = 0;

}

//________________________________________________________________________
void AliEbyEPidRatioPhy::CreateHistograms() {

  Float_t etaRange[2];
  fESDTrackCuts->GetEtaRange(etaRange[0],etaRange[1]);
  Float_t ptRange[2];
  fESDTrackCuts->GetPtRange(ptRange[0],ptRange[1]);

  
  
  
  Int_t    binHnUnCorr[7] = {AliEbyEPidRatioHelper::fgkfHistNBinsCent, 4,
			     AliEbyEPidRatioHelper::fgkfHistNBinsSign,
			     AliEbyEPidRatioHelper::fgkfHistNBinsEta, 
			     AliEbyEPidRatioHelper::fgkfHistNBinsRap,  
			     AliEbyEPidRatioHelper::fgkfHistNBinsPhi,   
			     AliEbyEPidRatioHelper::fgkfHistNBinsPt };      
  
  Double_t minHnUnCorr[7] = {AliEbyEPidRatioHelper::fgkfHistRangeCent[0],-0.5, 
			     AliEbyEPidRatioHelper::fgkfHistRangeSign[0],
			     AliEbyEPidRatioHelper::fgkfHistRangeEta[0], 
			     AliEbyEPidRatioHelper::fgkfHistRangeRap[0],  
			     AliEbyEPidRatioHelper::fgkfHistRangePhi[0], 
			     AliEbyEPidRatioHelper::fgkfHistRangePt[0]};  
			     
  
  Double_t maxHnUnCorr[7] = {AliEbyEPidRatioHelper::fgkfHistRangeCent[1],3.5, 
			     AliEbyEPidRatioHelper::fgkfHistRangeSign[1],
			     AliEbyEPidRatioHelper::fgkfHistRangeEta[1], 
			     AliEbyEPidRatioHelper::fgkfHistRangeRap[1],  
			     AliEbyEPidRatioHelper::fgkfHistRangePhi[1], 
			     AliEbyEPidRatioHelper::fgkfHistRangePt[1]};  
			     
  
  // -- UnCorrected

  if (fIsQA) {
    fOutList->Add(new THnSparseD("hnTrackUnCorrRec", "cent:pid:sign:eta:y:phi:pt", 7, binHnUnCorr, minHnUnCorr, maxHnUnCorr));  
    fHnTrackUnCorrRec = static_cast<THnSparseD*>(fOutList->Last());
    fHnTrackUnCorrRec->Sumw2(); 
    fHnTrackUnCorrRec->GetAxis(0)->SetTitle("centrality");
    fHnTrackUnCorrRec->GetAxis(1)->SetTitle("N #pi K  P");
    fHnTrackUnCorrRec->GetAxis(2)->SetTitle("sign");
    fHnTrackUnCorrRec->GetAxis(3)->SetTitle("#eta");
    fHnTrackUnCorrRec->GetAxis(4)->SetTitle("#it{y}");
    fHnTrackUnCorrRec->GetAxis(5)->SetTitle("#varphi");
    fHnTrackUnCorrRec->GetAxis(6)->SetTitle("#it{p}_{T} (GeV/#it{c})");

    if (fIsMC)  {
      fOutList->Add(new THnSparseD("hnTrackUnCorrMc", "cent:pid:sign:eta:y:phi:pt", 7, binHnUnCorr, minHnUnCorr, maxHnUnCorr));  
      fHnTrackUnCorrMc = static_cast<THnSparseD*>(fOutList->Last());
      fHnTrackUnCorrMc->Sumw2(); 
      fHnTrackUnCorrMc->GetAxis(0)->SetTitle("centrality");
      fHnTrackUnCorrMc->GetAxis(1)->SetTitle("N #pi K  P");
      fHnTrackUnCorrMc->GetAxis(2)->SetTitle("sign");
      fHnTrackUnCorrMc->GetAxis(3)->SetTitle("#eta");
      fHnTrackUnCorrMc->GetAxis(4)->SetTitle("#it{y}");
      fHnTrackUnCorrMc->GetAxis(5)->SetTitle("#varphi");
      fHnTrackUnCorrMc->GetAxis(6)->SetTitle("#it{p}_{T} (GeV/#it{c})");
      
    }

  }





  // fHelper->BinLogAxis(fHnTrackUnCorr, 4, fESDTrackCuts);

  // for (Int_t idx = 1 ; idx <= fHnTrackUnCorr->GetAxis(4)->GetNbins(); ++idx)
  //  printf("%02d |  %f > %f < %f\n", idx, fHnTrackUnCorr->GetAxis(4)->GetBinLowEdge(idx), fHnTrackUnCorr->GetAxis(4)->GetBinCenter(idx), fHnTrackUnCorr->GetAxis(4)->GetBinUpEdge(idx));


  fRan = new TRandom3();
  fRan->SetSeed();

  TString sTitle("");
  fPtBinHist = new TH1F("hPtBinHist","Make the pT Bins",AliEbyEPidRatioHelper::fgkfHistNBinsPt, AliEbyEPidRatioHelper::fgkfHistRangePt[0], AliEbyEPidRatioHelper::fgkfHistRangePt[1]);
  if (fIsRatio) AddHistSetRatio("Ratio",       Form("%s, #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), ptRange[0], ptRange[1]));
  
  AddHistSetCent("Phy",    Form("%s, #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), ptRange[0], ptRange[1]));
  if (fIsDetectorWise) AddHistSetCent("PhyTPC", Form("%s, #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), ptRange[0], fHelper->GetMinPtForTOFRequired()));
  if (fIsDetectorWise) AddHistSetCent("PhyTOF", Form("%s, #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), fHelper->GetMinPtForTOFRequired(), ptRange[1]));
  
  if (fIsPtBin) AddHistSetCentPt("PhyBin",    Form("%s, #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), ptRange[0], ptRange[1]));  

  if (fIsSub) AddHistInGroup("PhySS",    Form("%s, #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), ptRange[0], ptRange[1]),fHelper->GetNSubSamples());  
  if (fIsBS)  AddHistInGroup("PhyBS",    Form("%s, #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), ptRange[0], ptRange[1]),fHelper->GetNSubSamples());  


#if USE_PHI
  AddHistSetCent("Phyphi", Form("%s,#it{p}_{T} [%.1f,%.1f], #varphi [%.2f,%.2f]", sTitle.Data(), ptRange[0], ptRange[1], fHelper->GetPhiMin(), fHelper->GetPhiMax()));
  if (fIsDetectorWise) AddHistSetCent("PhyTPCphi",Form("%s, #it{p}_{T} [%.1f,%.1f], #varphi [%.2f,%.2f]",
						       sTitle.Data(), ptRange[0], fHelper->GetMinPtForTOFRequired(), fHelper->GetPhiMin(), fHelper->GetPhiMax()));
  if (fIsDetectorWise) AddHistSetCent("PhyTOFphi",Form("%s, #it{p}_{T} [%.1f,%.1f], #varphi [%.2f,%.2f]",
						       sTitle.Data(), fHelper->GetMinPtForTOFRequired(), ptRange[1], fHelper->GetPhiMin(), fHelper->GetPhiMax()));
  if (fIsPtBin) AddHistSetCentPt("PhyBinPhi", Form("%s, #it{p}_{T} [%.1f,%.1f], #varphi [%.2f,%.2f]", 
				     sTitle.Data(), ptRange[0], ptRange[1], fHelper->GetPhiMin(), fHelper->GetPhiMax()));
#endif

  if (fIsMC) {
    TString sMCTitle("");
  
    if (fIsRatio) AddHistSetRatio("MCRatio",       Form("%s, #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), ptRange[0], ptRange[1]));

    AddHistSetCent("MC",      Form("%s", sTitle.Data()));
    AddHistSetCent("MCpt",    Form("%s, #it{p}_{T} [%.1f,%.1f]", sMCTitle.Data(), ptRange[0], ptRange[1]));
    if (fIsDetectorWise) AddHistSetCent("MCTPC",   Form("%s, #it{p}_{T} [%.1f,%.1f]", sMCTitle.Data(), ptRange[0], fHelper->GetMinPtForTOFRequired()));
    if (fIsDetectorWise) AddHistSetCent("MCTOF",   Form("%s, #it{p}_{T} [%.1f,%.1f]", sMCTitle.Data(), fHelper->GetMinPtForTOFRequired(), ptRange[1]));
    if (fIsPtBin) AddHistSetCentPt("MCBin",        Form("%s, #it{p}_{T} [%.1f,%.1f]", sMCTitle.Data(), ptRange[0], ptRange[1])); 
    if (fIsSub) AddHistInGroup("McSS",      Form("%s, #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), ptRange[0], ptRange[1]),fHelper->GetNSubSamples());  
    if (fIsBS)  AddHistInGroup("McBS",      Form("%s, #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), ptRange[0], ptRange[1]),fHelper->GetNSubSamples());  

#if USE_PHI
    AddHistSetCent("MCPhi",   Form("%s, #it{p}_{T} [%.1f,%.1f], #varphi [%.2f,%.2f]",sMCTitle.Data(), ptRange[0], ptRange[1], fHelper->GetPhiMin(), fHelper->GetPhiMax()));
    if (fIsDetectorWise) AddHistSetCent("MCTPCphi",Form("%s, #it{p}_{T} [%.1f,%.1f], #varphi [%.2f,%.2f]",
							sMCTitle.Data(), ptRange[0], fHelper->GetMinPtForTOFRequired(), fHelper->GetPhiMin(), fHelper->GetPhiMax()));
    if (fIsDetectorWise) AddHistSetCent("MCTOFphi",Form("%s, #it{p}_{T} [%.1f,%.1f], #varphi [%.2f,%.2f]",
							sMCTitle.Data(), fHelper->GetMinPtForTOFRequired(), ptRange[1], fHelper->GetPhiMin(), fHelper->GetPhiMax()));
 if (fIsPtBin) AddHistSetCentPt("MCBinPhi", Form("%s, #it{p}_{T} [%.1f,%.1f], #varphi [%.2f,%.2f]", 
					sMCTitle.Data(), ptRange[0], ptRange[1], fHelper->GetPhiMin(), fHelper->GetPhiMax()));

  
#endif
    
  }
  
  return;
}

//________________________________________________________________________
Int_t AliEbyEPidRatioPhy::ProcessTracks() {
  Float_t etaRange[2];
  fESDTrackCuts->GetEtaRange(etaRange[0],etaRange[1]);
  Float_t ptRange[2];
  fESDTrackCuts->GetPtRange(ptRange[0],ptRange[1]);
  
  // -- Track Loop
  for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {
    AliVTrack *track = (fESD) ? static_cast<AliVTrack*>(fESD->GetTrack(idxTrack)) : static_cast<AliVTrack*>(fAOD->GetTrack(idxTrack)); 
    // -- Check if track is accepted for basic parameters
    if (!fHelper->IsTrackAcceptedBasicCharged(track))
      continue;
    // -- Check if accepted - ESD
    if (fESD && !fESDTrackCuts->AcceptTrack(dynamic_cast<AliESDtrack*>(track)))
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
      if(!(track->Pt() > ptRange[0] && track->Pt() <= ptRange[1] && TMath::Abs(track->Eta()) <= etaRange[1]))
	continue;
    }

    if (!fHelper->IsTrackAcceptedDCA(track))
      continue;
    
    Int_t iPid = 0;
    Double_t pid[3];
    if      (fHelper->IsTrackAcceptedPID(track, pid, (AliPID::kPion)))   iPid = 1;
    else if (fHelper->IsTrackAcceptedPID(track, pid, (AliPID::kKaon)))   iPid = 2;
    else if (fHelper->IsTrackAcceptedPID(track, pid, (AliPID::kProton))) iPid = 3;
    else iPid = 0;
   

    Double_t yP;
    if (iPid != 0 && !fHelper->IsTrackAcceptedRapidity(track, yP, iPid))
      continue;
    
    
    if (fIsQA) {
      if (iPid != 0) {
	Double_t aTrack[7] = {static_cast<Double_t>(fCentralityBin), 0, static_cast<Double_t>(track->Charge()), track->Eta(),yP,track->Phi(),track->Pt()};
	fHnTrackUnCorrRec->Fill(aTrack);
      }
      Double_t aTrack[7] = {static_cast<Double_t>(fCentralityBin), static_cast<Double_t>(iPid), static_cast<Double_t>(track->Charge()), track->Eta(),yP,track->Phi(),track->Pt()};
      fHnTrackUnCorrRec->Fill(aTrack);
    }

   
    
    Int_t idxPart = (track->Charge() < 0) ? 0 : 1;
    // -- in pt Range
    fNp[0][0][idxPart] += 1;
    if(iPid != 0) fNp[0][iPid][idxPart] += 1;

    Int_t idxPt = fPtBinHist->FindBin(track->Pt()) -1;

    if (idxPt>=0 && idxPt < AliEbyEPidRatioHelper::fgkfHistNBinsPt) {
      fNpPt[0][0][idxPart][idxPt] += 1;
      if(iPid != 0) fNpPt[0][iPid][idxPart][idxPt] += 1;
    }

    // -- in TPC pt Range
    if (track->Pt() <= fHelper->GetMinPtForTOFRequired()) {
      fNp[1][0][idxPart] += 1;
      if(iPid != 0) fNp[1][iPid][idxPart] += 1;
    }
    // -- in TPC+TOF pt Range
    if (track->Pt() > fHelper->GetMinPtForTOFRequired()){
      fNp[2][0][idxPart] += 1;
      if(iPid != 0) fNp[2][iPid][idxPart] += 1;
    }

#if USE_PHI
    if(!fHelper->IsTrackAcceptedPhi(track))
      continue;
    
    // -- in pt Range
    fNp[3][iPid][idxPart] += 1;
    if(iPid != 0) fNp[3][0][idxPart] += 1;

    if (idxPt>=0 && idxPt < AliEbyEPidRatioHelper::fgkfHistNBinsPt) {
      fNpPt[1][0][idxPart][idxPt] += 1;
      if(iPid != 0) fNpPt[1][iPid][idxPart][idxPt] += 1;
    }

    // -- in TPC pt Range
    if (track->Pt() <= fHelper->GetMinPtForTOFRequired()) {
      fNp[4][0][idxPart] += 1;
      if(iPid != 0)fNp[4][iPid][idxPart] += 1;
    }
    // -- in TPC+TOF pt Range
    if (track->Pt() > fHelper->GetMinPtForTOFRequired()) {
      fNp[5][0][idxPart] += 1;
      if(iPid != 0) fNp[5][iPid][idxPart] += 1;
    }
#endif
  } // for (Int_t idxTrack = 0; idxTrack < fESD->GetNumberOfTracks(); ++idxTrack) {
 
  FillHistSetCent("Phy",    0, kFALSE);
  if (fIsDetectorWise) FillHistSetCent("PhyTPC", 1, kFALSE);
  if (fIsDetectorWise) FillHistSetCent("PhyTOF", 2, kFALSE);
  
 if (fIsRatio) FillHistSetRatio("Ratio",   0, kFALSE);
 if (fIsPtBin) FillHistSetCentPt("PhyBin", 0, kFALSE);

 if (fIsSub) FillHistInGroup("PhySS", 0,fHelper->GetSubSampleIdx(),kFALSE);

 

 if (fIsBS)  {
   Int_t sap = fHelper->GetNSubSamples();
   for (Int_t i = 0; i < sap; ++i)  
     FillHistInGroup("PhyBS",0,fRan->Integer(sap),kFALSE);
 }

#if USE_PHI
  FillHistSetCent("Phyphi",    3, kFALSE);
  if (fIsDetectorWise) FillHistSetCent("PhyTPCphi", 4, kFALSE);
  if (fIsDetectorWise) FillHistSetCent("PhyTOFphi", 5, kFALSE);
  if (fIsPtBin) FillHistSetCentPt("PhyBinPhi", 1, kFALSE);
 #endif

  return 0;
}

//________________________________________________________________________
Int_t AliEbyEPidRatioPhy::ProcessParticles() {
  // -- Process primary particles from the stack and fill histograms
  
  Float_t etaRange[2];
  fESDTrackCuts->GetEtaRange(etaRange[0],etaRange[1]);
  Float_t ptRange[2];
  fESDTrackCuts->GetPtRange(ptRange[0],ptRange[1]);


  for (Int_t idxMC = 0; idxMC < fStack->GetNprimary(); ++idxMC) {
    AliVParticle* particle = (fESD) ? fMCEvent->GetTrack(idxMC) : NULL;

    if (!particle) 
      continue;
    if (!fHelper->IsParticleAcceptedBasicCharged(particle, idxMC))
      continue;
    
    Int_t iPid = 0;  
    if      (TMath::Abs(particle->PdgCode()) ==  211) iPid = 1; // pion
    else if (TMath::Abs(particle->PdgCode()) ==  321) iPid = 2; // kaon
    else if (TMath::Abs(particle->PdgCode()) == 2212) iPid = 3; // proton
    else    iPid = 0;
    
    Double_t yMC;
    if ((iPid != 0) && !fHelper->IsParticleAcceptedRapidity(particle, yMC, iPid))
      continue;

    // -- Check eta window -- for charged particles
    if ((iPid == 0) && TMath::Abs(particle->Eta()) > etaRange[1])
      continue;
    
    Int_t idxPart = (particle->PdgCode() < 0) ? 0 : 1;
   

    if (fIsQA) {
      Float_t signMC    = (particle->PdgCode() < 0) ? -1. : 1.;
      if (iPid != 0) {
	Double_t aTrack[7] = {static_cast<Double_t>(fCentralityBin), 0, static_cast<Double_t>(signMC), particle->Eta(),yMC,particle->Phi(),particle->Pt()};
	fHnTrackUnCorrMc->Fill(aTrack);
      }
      Double_t aTrack[7] = {static_cast<Double_t>(fCentralityBin), static_cast<Double_t>(iPid), static_cast<Double_t>(signMC), particle->Eta(),yMC,particle->Phi(),particle->Pt()};
      fHnTrackUnCorrMc->Fill(aTrack);
    }


    // idx 0
    fMCNp[0][0][idxPart]    += 1.;        
    if(iPid != 0) fMCNp[0][iPid][idxPart] += 1.;        
    
    // -- Check main pt window
    if (!(particle->Pt() > ptRange[0] && particle->Pt() <= ptRange[1]))
      continue;
    
    // -- in pt Range idx 1
    fMCNp[1][0][idxPart]    += 1.;        
    if(iPid != 0)fMCNp[1][iPid][idxPart] += 1.;        
    

    

    Int_t idxPt = fPtBinHist->FindBin(particle->Pt()) -1;
    
    if (idxPt>=0 && idxPt < AliEbyEPidRatioHelper::fgkfHistNBinsPt) {
      fMCNpPt[0][0][idxPart][idxPt] += 1;
      if(iPid != 0) fMCNpPt[0][iPid][idxPart][idxPt] += 1;
    }


      
#if USE_PHI
    if(!fHelper->IsParticleAcceptedPhi(particle))
      continue;
    // idxPhi = 1;

    if (idxPt>=0 && idxPt < AliEbyEPidRatioHelper::fgkfHistNBinsPt) {
      fMCNpPt[1][0][idxPart][idxPt] += 1;
      if(iPid != 0) fMCNpPt[1][iPid][idxPart][idxPt] += 1;
    }

    // -- in pt Range
    fMCNp[2][0][idxPart]    += 1;
    if(iPid != 0)fMCNp[2][iPid][idxPart] += 1;
    
  
#endif
  } // for (Int_t idxMC = 0; idxMC < nPart; ++idxMC) {
  
  FillHistSetCent("MC",   0, kTRUE);
  FillHistSetCent("MCpt", 1, kTRUE);

  if (fIsRatio) FillHistSetRatio("MCRatio", 0, kTRUE);
  if (fIsPtBin) FillHistSetCentPt("MCBin",  0, kTRUE);

  if (fIsSub) FillHistInGroup("McSS", 0,fHelper->GetSubSampleIdx(),kTRUE);
  if (fIsBS)  {
    Int_t sap = fHelper->GetNSubSamples();
    for (Int_t i = 0; i < sap; ++i)  
      FillHistInGroup("McBS",0,fRan->Integer(sap),kTRUE);
  }
  
#if USE_PHI
  FillHistSetCent("MCphi",  2, kTRUE);
  if (fIsPtBin)  FillHistSetCentPt("MCBinPhi", 1, kTRUE);
#endif

  return 0;
}

//________________________________________________________________________
void  AliEbyEPidRatioPhy::AddHistSetCent(const Char_t *name, const Char_t *title)  {
  TString sName(name);
  TString sTitle(title);
 
  Float_t etaRange[2];

  fESDTrackCuts->GetEtaRange(etaRange[0],etaRange[1]);

  //TList *list[4];
  fOutList->Add(new TList);
  TList *list =  static_cast<TList*>(fOutList->Last());
  list->SetName(Form("f%s", name));
  list->SetOwner(kTRUE);
  

  Int_t nBinsCent         =  AliEbyEPidRatioHelper::fgkfHistNBinsCent;
  Double_t centBinRange[] = {AliEbyEPidRatioHelper::fgkfHistRangeCent[0], AliEbyEPidRatioHelper::fgkfHistRangeCent[1]};

 
  for (Int_t iPid = 0; iPid < 4; ++iPid) {
    TString sNetTitle(Form("%s - %s", AliEbyEPidRatioHelper::fgkPidLatex[iPid][1], AliEbyEPidRatioHelper::fgkPidLatex[iPid][0]));

    sTitle = (iPid != 0 ) ? Form("|y| < %.1f", fHelper->GetRapidityMax()) : Form(" |#eta|<%.1f", etaRange[1]);

    list->Add(new TProfile(Form("fProfTot%sPlus%s", AliEbyEPidRatioHelper::fgkPidName[iPid],name), 
			   Form("(%s) : %s;Centrality(100);(%s)",AliEbyEPidRatioHelper::fgkPidName[iPid], sTitle.Data(), sNetTitle.Data()),
			   100,-0.5,99.5));

    list->Add(new TProfile(Form("fProfTot%sMinus%s", AliEbyEPidRatioHelper::fgkPidName[iPid],name), 
			   Form("(%s) : %s;Centrality(100);(%s)",AliEbyEPidRatioHelper::fgkPidName[iPid], sTitle.Data(), sNetTitle.Data()),
			   100,-0.5,99.5));


    

    for (Int_t idx = 1; idx <= fOrder; ++idx) {
      list->Add(new TProfile(Form("fProf%s%sNet%dM", AliEbyEPidRatioHelper::fgkPidName[iPid],name, idx), 
			     Form("(%s)^{%d} : %s;Centrality(100);(%s)^{%d}",sNetTitle.Data(), idx, sTitle.Data(), sNetTitle.Data(), idx),
			     100,-0.5,99.5));
    }
    
    for (Int_t ii = 0; ii <= fOrder; ++ii) {
      for (Int_t kk = 0; kk <= fOrder; ++kk) {
	list->Add(new TProfile(Form("fProf%s%sNetF%02d%02d", AliEbyEPidRatioHelper::fgkPidName[iPid], name, ii, kk),
			       Form("f_{%02d%02d} : %s;Centrality(100);f_{%02d%02d}", ii, kk, sTitle.Data(), ii, kk),
			       100,-0.5,99.5));
      }
    }
  
  }  
   
  for (Int_t iPid = 0; iPid < 4; ++iPid) {
    TString sNetTitle(Form("%s - %s", AliEbyEPidRatioHelper::fgkPidLatex[iPid][1], AliEbyEPidRatioHelper::fgkPidLatex[iPid][0]));
    sTitle = (iPid != 0 ) ? Form(" |y|<%.1f", fHelper->GetRapidityMax()) : Form(" |#eta| < %.1f", etaRange[1]);

    list->Add(new TProfile(Form("fProfBinTot%sPlus%s", AliEbyEPidRatioHelper::fgkPidName[iPid],name), 
			   Form("(%s) : %s;Centrality(11);(%s)",AliEbyEPidRatioHelper::fgkPidName[iPid], sTitle.Data(), sNetTitle.Data()),
			   nBinsCent, centBinRange[0], centBinRange[1]));

    list->Add(new TProfile(Form("fProfBinTot%sMinus%s", AliEbyEPidRatioHelper::fgkPidName[iPid],name), 
			   Form("(%s) : %s;Centrality(11);(%s)",AliEbyEPidRatioHelper::fgkPidName[iPid], sTitle.Data(), sNetTitle.Data()),
			   nBinsCent, centBinRange[0], centBinRange[1]));



    for (Int_t idx = 1; idx <= fOrder; ++idx) {
      list->Add(new TProfile(Form("fProfBin%s%sNet%dM", AliEbyEPidRatioHelper::fgkPidName[iPid],name, idx), 
			     Form("(%s)^{%d} : %s;Centrality(11);(%s)^{%d}", sNetTitle.Data(), idx, sTitle.Data(), sNetTitle.Data(), idx),
			     nBinsCent, centBinRange[0], centBinRange[1]));
    }
    
    for (Int_t ii = 0; ii <= fOrder; ++ii) {
      for (Int_t kk = 0; kk <= fOrder; ++kk) {
	list->Add(new TProfile(Form("fProfBin%s%sNetF%02d%02d", AliEbyEPidRatioHelper::fgkPidName[iPid], name, ii, kk),
			       Form("f_{%02d%02d} : %s;Centrality(11);f_{%02d%02d}", ii, kk, sTitle.Data(), ii, kk),
			       nBinsCent, centBinRange[0], centBinRange[1]));
      }
    }
  
  }  
  
  for (Int_t iPhy = 0; iPhy < 46; ++iPhy) { 
    list->Add(new TProfile(Form("fProf%sNu%02d",name,iPhy),Form("Physics Variable for index %d | %s ; Centrality;",iPhy,name),100,-0.5,99.5));
  }
  for (Int_t iPhy = 0; iPhy < 46; ++iPhy) { 
    list->Add(new TProfile(Form("fProfBin%sNu%02d",name,iPhy),Form("Physics Variable for index %d | %s ; Centrality;",iPhy,name),nBinsCent, centBinRange[0], centBinRange[1]));
  }
  
  return;
}



//________________________________________________________________________
void  AliEbyEPidRatioPhy::AddHistSetRatio(const Char_t *name, const Char_t *title)  {
  TString sName(name);
  TString sTitle(title);
 
  Float_t etaRange[2];
  fESDTrackCuts->GetEtaRange(etaRange[0],etaRange[1]);

  //TList *list[4];
  fOutList->Add(new TList);
  TList *list =  static_cast<TList*>(fOutList->Last());
  list->SetName(Form("f%s", name));
  list->SetOwner(kTRUE);
  
  Int_t    nRbin  = 15000;
  Double_t mRat[] = {0,1.5};

  Int_t nBinsCent         =  (fIsPer) ? 100 : AliEbyEPidRatioHelper::fgkfHistNBinsCent;
  Double_t centBinRange[2];  
  centBinRange[0]  =  (fIsPer) ?  0   : AliEbyEPidRatioHelper::fgkfHistRangeCent[0];
  centBinRange[1]  =  (fIsPer) ?  100 : AliEbyEPidRatioHelper::fgkfHistRangeCent[1];

  TString xyz = Form("|y| < %.1f",fHelper->GetRapidityMax()); 

  list->Add(new TH2F(Form("fHistRatioKPi%s",name), 
		     Form("(%s %s) : K/#pi;Centrality(11);K/#pi", xyz.Data(), sTitle.Data()),
		     nBinsCent, centBinRange[0], centBinRange[1], nRbin,mRat[0],mRat[1]));
  
  list->Add(new TH2F(Form("fHistRatioKpPip%s",name), 
		     Form("(%s %s) : K^{+}/#pi^{+};Centrality(11);K^{+}/#pi^{+}", xyz.Data(), sTitle.Data()),
		     nBinsCent, centBinRange[0], centBinRange[1], nRbin,mRat[0],mRat[1]));
  
  list->Add(new TH2F(Form("fHistRatioKmPip%s",name), 
		     Form("(%s %s) : K^{-}/#pi^{+};Centrality(11);K^{-}/#pi^{+}", xyz.Data(), sTitle.Data()),
		     nBinsCent, centBinRange[0], centBinRange[1], nRbin,mRat[0],mRat[1]));
  
  list->Add(new TH2F(Form("fHistRatioKmPim%s",name), 
		     Form("(%s %s) : K^{-}/#pi^{-};Centrality(11);K^{-}/#pi^{-}", xyz.Data(), sTitle.Data()),
		     nBinsCent, centBinRange[0], centBinRange[1], nRbin,mRat[0],mRat[1]));



 list->Add(new TH2F(Form("fHistRatioPK%s",name), 
		     Form("(%s %s) : P/K;Centrality(11);P/K", xyz.Data(), sTitle.Data()),
		     nBinsCent, centBinRange[0], centBinRange[1], nRbin,mRat[0],mRat[1]));
  
  list->Add(new TH2F(Form("fHistRatioPpKp%s",name), 
		     Form("(%s %s) : P/K^{+};Centrality(11);P/K^{+}", xyz.Data(), sTitle.Data()),
		     nBinsCent, centBinRange[0], centBinRange[1], nRbin,mRat[0],mRat[1]));
  
  list->Add(new TH2F(Form("fHistRatioPmKp%s",name), 
		     Form("(%s %s) : #bar{P}/K^{+};Centrality(11);#bar{P}/K^{+}", xyz.Data(), sTitle.Data()),
		     nBinsCent, centBinRange[0], centBinRange[1], nRbin,mRat[0],mRat[1]));
  
  list->Add(new TH2F(Form("fHistRatioPmKm%s",name), 
		     Form("(%s %s) : #bar{P}/K^{-};Centrality(11);#bar{P}/K^{-}", xyz.Data(), sTitle.Data()),
		     nBinsCent, centBinRange[0], centBinRange[1], nRbin,mRat[0],mRat[1]));



 list->Add(new TH2F(Form("fHistRatioPPi%s",name), 
		     Form("(%s %s) : P/#pi;Centrality(11);K/#pi", xyz.Data(), sTitle.Data()),
		     nBinsCent, centBinRange[0], centBinRange[1], nRbin,mRat[0],mRat[1]));
  
  list->Add(new TH2F(Form("fHistRatioPpPip%s",name), 
		     Form("(%s %s) : P/#pi^{+};Centrality(11);P/#pi^{+}", xyz.Data(), sTitle.Data()),
		     nBinsCent, centBinRange[0], centBinRange[1], nRbin,mRat[0],mRat[1]));
  
  list->Add(new TH2F(Form("fHistRatioPmPip%s",name), 
		     Form("(%s %s) : #bar{P}/#pi^{+};Centrality(11);#bar{P}/#pi^{+}", xyz.Data(), sTitle.Data()),
		     nBinsCent, centBinRange[0], centBinRange[1], nRbin,mRat[0],mRat[1]));
  
  list->Add(new TH2F(Form("fHistRatioPmPim%s",name), 
		     Form("(%s %s) : #bar{P}/#pi^{-};Centrality(11);#bar{P}/#pi^{-}", xyz.Data(), sTitle.Data()),
		     nBinsCent, centBinRange[0], centBinRange[1], nRbin,mRat[0],mRat[1]));
  

  //------- - - -  -  -   -   - - -   - --- - - --- - - - - - -- --------
  Int_t bin[4] = {2800,2200,1200,600}; 
  Int_t bd[] = {1,2,2,2};
  for (Int_t iPid = 0; iPid < 4; ++iPid) {
    Int_t bb = bin[iPid];

    for (Int_t iNet = 0; iNet < 4; ++iNet) {
      Int_t bn     = (iNet == 3) ?  500   : bb/bd[iNet]; 
      Float_t blow = (iNet == 3) ? -250.5 : -0.5;
      Float_t bup  = (iNet == 3) ?  249.5 : bn-0.5;


      list->Add(new TH2F(Form("fHistDist%s%s%s",name, AliEbyEPidRatioHelper::fgkPidName[iPid], AliEbyEPidRatioHelper::fgkNetHistName[iNet]), 
			 Form("(%s %s) : %s Distribution;Centrality(11);%s_{(%s)}", xyz.Data(), sTitle.Data(), 
			      AliEbyEPidRatioHelper::fgkPidShLatex[iPid],AliEbyEPidRatioHelper::fgkPidShLatex[iPid],AliEbyEPidRatioHelper::fgkNetHistLatex[iNet]),
			 nBinsCent, centBinRange[0], centBinRange[1], bn, blow,bup));    
    }
  }
  



return;
}

//________________________________________________________________________
void AliEbyEPidRatioPhy::FillHistSetCent(const Char_t *name, Int_t idx, Bool_t isMC)  {
  Int_t ***np = (isMC) ? fMCNp : fNp;
  
  Float_t centralityBin = fHelper->GetCentralityBin();
  Float_t centralityPer = fHelper->GetCentralityPercentile();
  
  TList *list = static_cast<TList*>(fOutList->FindObject(Form("f%s",name)));
  
  for (Int_t iPid = 0; iPid < 4; ++iPid) {
    Int_t deltaNp = np[idx][iPid][1]-np[idx][iPid][0];  
    Double_t delta = 1.;

    (static_cast<TProfile*>(list->FindObject(Form("fProfBinTot%sPlus%s", AliEbyEPidRatioHelper::fgkPidName[iPid], name))))->Fill(centralityBin, np[idx][iPid][1]);
    (static_cast<TProfile*>(list->FindObject(Form("fProfTot%sPlus%s", AliEbyEPidRatioHelper::fgkPidName[iPid], name))))->Fill(centralityPer, np[idx][iPid][1]);

    (static_cast<TProfile*>(list->FindObject(Form("fProfBinTot%sMinus%s", AliEbyEPidRatioHelper::fgkPidName[iPid], name))))->Fill(centralityBin, np[idx][iPid][0]);
    (static_cast<TProfile*>(list->FindObject(Form("fProfTot%sMinus%s", AliEbyEPidRatioHelper::fgkPidName[iPid], name))))->Fill(centralityPer, np[idx][iPid][0]);


    for (Int_t idxOrder = 1; idxOrder <= fOrder; ++idxOrder) {
      delta *= deltaNp;

      (static_cast<TProfile*>(list->FindObject(Form("fProfBin%s%sNet%dM", AliEbyEPidRatioHelper::fgkPidName[iPid], name, idxOrder))))->Fill(centralityBin, delta);
      (static_cast<TProfile*>(list->FindObject(Form("fProf%s%sNet%dM", AliEbyEPidRatioHelper::fgkPidName[iPid], name, idxOrder))))->Fill(centralityPer, delta);
    }
    
    for (Int_t idxOrder = 0; idxOrder <= fOrder; ++ idxOrder) {
      fRedFactp[idxOrder][0]  = 1.;
      fRedFactp[idxOrder][1]  = 1.;
    }
    
    for (Int_t idxOrder = 1; idxOrder <= fOrder; ++ idxOrder) {
      fRedFactp[idxOrder][0]  = fRedFactp[idxOrder-1][0]  * Double_t(np[idx][iPid][0]-(idxOrder-1));
      fRedFactp[idxOrder][1]  = fRedFactp[idxOrder-1][1]  * Double_t(np[idx][iPid][1]-(idxOrder-1));
    }
    
    for (Int_t ii = 0; ii <= fOrder; ++ii) {  
      for (Int_t kk = 0; kk <= fOrder; ++kk) { 
	Double_t fik = fRedFactp[ii][1] * fRedFactp[kk][0];   
	(static_cast<TProfile*>(list->FindObject(Form("fProfBin%s%sNetF%02d%02d", AliEbyEPidRatioHelper::fgkPidName[iPid], name, ii, kk))))->Fill(centralityBin, fik);
	(static_cast<TProfile*>(list->FindObject(Form("fProf%s%sNetF%02d%02d", AliEbyEPidRatioHelper::fgkPidName[iPid], name, ii, kk))))->Fill(centralityPer, fik);
      }
    }
  }
 
  //Printf("%6d %20s %6.2f %6d %6d %6d %6d  %6d %6d %6d %6d", idx, name, centralityBin,
  //	 np[idx][0][1],  np[idx][0][0], 
  //	 np[idx][1][1],  np[idx][1][0], 
  ///	 np[idx][2][1],  np[idx][2][0], 
  //	 np[idx][3][1],  np[idx][3][0]);
  //

   Int_t a[6][4]; Int_t b[22];
   for (Int_t iPid = 0; iPid < 4; ++iPid) {
     a[0][iPid] = np[idx][iPid][1]+np[idx][iPid][0];       // 0  n+ + n-
     a[1][iPid] = np[idx][iPid][1];                        // 1  n+
     a[2][iPid] = np[idx][iPid][0];                        // 2  n-
     a[3][iPid] = np[idx][iPid][1]*np[idx][iPid][0];       // 3  n+ . n-
     a[4][iPid] = np[idx][iPid][1]*(np[idx][iPid][1]-1);   // 4  n+ (n+ - 1)
     a[5][iPid] = np[idx][iPid][0]*(np[idx][iPid][0]-1);   // 5  n- (n- - 1)
     
     // Printf("%6d %20s %6.2f %6d %6d %6d ", idx, name, centralityBin,
     //	   a[0][iPid], a[1][iPid], a[2][iPid]);

  }
  
  b[0]  = a[0][0]*a[0][2];       // 24 N   K
  b[1]  = a[0][1]*a[0][2];       // 25 Pi  K
  b[2]  = a[1][1]*a[1][2];       // 26 pi+ k+
  b[3]  = a[1][1]*a[2][2];       // 27 pi+ k-
  b[4]  = a[2][1]*a[1][2];       // 28 pi- k+  
  b[5]  = a[2][1]*a[2][2];       // 29 pi- k-
  
  b[6]  = a[0][0]*a[0][3];       // 30 N   P
  b[7]  = a[0][2]*a[0][3];       // 31 K   P
  b[8]  = a[1][2]*a[1][3];       // 32 k+  p+
  b[9]  = a[1][2]*a[2][3];       // 33 k+  p-
  b[10] = a[2][2]*a[1][3];       // 34 k-  p+
  b[11] = a[2][2]*a[2][3];       // 35 k-  p-
  
  b[12] = a[0][0]*a[0][1];       // 36 N  Pi
  b[13] = a[0][3]*a[0][1];       // 37 P  Pi
  b[14] = a[1][3]*a[1][1];       // 38 p+ pi+
  b[15] = a[1][3]*a[2][1];       // 39 p+ pi-
  b[16] = a[2][3]*a[1][1];       // 40 p- pi+
  b[17] = a[2][3]*a[2][1];       // 41 p- pi-
  
  b[18] = a[0][0]*(a[0][0] - 1); // 42 N ( N - 1 )
  b[19] = a[0][1]*(a[0][1] - 1); // 43 Pi( Pi- 1 )
  b[20] = a[0][2]*(a[0][1] - 1); // 44 K ( K - 1 )
  b[21] = a[0][3]*(a[0][3] - 1); // 45 P ( P - 1 )
  // TList *list_nu = static_cast<TList*>(fOutList->FindObject(Form("f%s_nu",name)));
  Int_t k = 0;
  for (Int_t j = 0; j < 4; j++) {
    for (Int_t i = 0; i < 6; i++) {
      (static_cast<TProfile*>(list->FindObject(Form("fProfBin%sNu%02d", name,k))))->Fill(centralityBin,a[i][j]); 
      (static_cast<TProfile*>(list->FindObject(Form("fProf%sNu%02d", name,k))))->Fill(centralityPer,a[i][j]); 
      k++;
    }
  }

  for (Int_t j = 0; j < 22; j++) {
    (static_cast<TProfile*>(list->FindObject(Form("fProfBin%sNu%02d", name,j+24))))->Fill(centralityBin,b[j]); 
    (static_cast<TProfile*>(list->FindObject(Form("fProf%sNu%02d", name,j+24))))->Fill(centralityPer,b[j]); 
  }
  
  return;
}


//________________________________________________________________________
void AliEbyEPidRatioPhy::FillHistSetRatio(const Char_t *name, Int_t idx, Bool_t isMC)  {
   
  Int_t ***np = (isMC) ? fMCNp : fNp;

  Float_t centralityBin = (fIsPer) ? fHelper->GetCentralityPercentile() : fHelper->GetCentralityBin();

  TList *list = static_cast<TList*>(fOutList->FindObject(Form("f%s",name)));
    
  
  if((Double_t)np[idx][1][1]+(Double_t)np[idx][1][0] != 0 ) {
    Double_t KPi = ((Double_t)np[idx][2][1]+(Double_t)np[idx][2][0])/((Double_t)np[idx][1][1]+(Double_t)np[idx][1][0]);
    Double_t PPi = ((Double_t)np[idx][3][1]+(Double_t)np[idx][3][0])/((Double_t)np[idx][1][1]+(Double_t)np[idx][1][0]);
    
    (static_cast<TH2F*>(list->FindObject(Form("fHistRatioKPi%s",name))))->Fill(centralityBin, KPi);
    (static_cast<TH2F*>(list->FindObject(Form("fHistRatioPPi%s",name))))->Fill(centralityBin,   PPi);
  }

  if((Double_t)np[idx][2][1]+(Double_t)np[idx][2][0] != 0 ){
    Double_t PK = ((Double_t)np[idx][3][1]+(Double_t)np[idx][3][0])/((Double_t)np[idx][2][1]+(Double_t)np[idx][2][0]);
    (static_cast<TH2F*>(list->FindObject(Form("fHistRatioPK%s",name))))->Fill(centralityBin, PK);
  }

  if ((Double_t)np[idx][1][1] != 0 ) {
    Double_t KpPip  = ((Double_t)np[idx][2][1])/((Double_t)np[idx][1][1]); 
    Double_t KmPip  = ((Double_t)np[idx][2][0])/((Double_t)np[idx][1][1]); 
    Double_t PpPip  = ((Double_t)np[idx][3][1])/((Double_t)np[idx][1][1]); 
    Double_t PmPip =  ((Double_t)np[idx][3][0])/((Double_t)np[idx][1][1]);

    (static_cast<TH2F*>(list->FindObject(Form("fHistRatioKpPip%s",name))))->Fill(centralityBin, KpPip);
    (static_cast<TH2F*>(list->FindObject(Form("fHistRatioKmPip%s",name))))->Fill(centralityBin, KmPip);
    (static_cast<TH2F*>(list->FindObject(Form("fHistRatioPpPip%s",name))))->Fill(centralityBin, PpPip);
    (static_cast<TH2F*>(list->FindObject(Form("fHistRatioPmPip%s",name))))->Fill(centralityBin, PmPip);
  }  

  if ((Double_t)np[idx][1][0] != 0) {
    Double_t KmPim = ((Double_t)np[idx][2][0])/((Double_t)np[idx][1][0]);
    Double_t PmPim = ((Double_t)np[idx][3][0])/((Double_t)np[idx][1][0]);
    (static_cast<TH2F*>(list->FindObject(Form("fHistRatioKmPim%s",name))))->Fill(centralityBin, KmPim);
    (static_cast<TH2F*>(list->FindObject(Form("fHistRatioPmPim%s",name))))->Fill(centralityBin, PmPim);
  }
  
  if ((Double_t)np[idx][2][1] != 0 ) { 
    Double_t PpKp  = ((Double_t)np[idx][3][1])/((Double_t)np[idx][2][1]); 
    Double_t PmKp =  ((Double_t)np[idx][3][0])/((Double_t)np[idx][2][1]); 
   (static_cast<TH2F*>(list->FindObject(Form("fHistRatioPpKp%s",name))))->Fill(centralityBin, PpKp);
   (static_cast<TH2F*>(list->FindObject(Form("fHistRatioPmKp%s",name))))->Fill(centralityBin, PmKp);
  }
  
  if ((Double_t)np[idx][2][0] != 0) {
   Double_t PmKm = ((Double_t)np[idx][3][0])/((Double_t)np[idx][2][0]);
   (static_cast<TH2F*>(list->FindObject(Form("fHistRatioPmKm%s",name))))->Fill(centralityBin, PmKm);
  }
  
  for (Int_t iPid = 0; iPid < 4; ++iPid) {
    (static_cast<TH2F*>(list->FindObject(Form("fHistDist%s%s%s",name, AliEbyEPidRatioHelper::fgkPidName[iPid], AliEbyEPidRatioHelper::fgkNetHistName[0]))))->Fill(centralityBin, np[idx][iPid][1]+np[idx][iPid][0]); 
    (static_cast<TH2F*>(list->FindObject(Form("fHistDist%s%s%s",name, AliEbyEPidRatioHelper::fgkPidName[iPid], AliEbyEPidRatioHelper::fgkNetHistName[1]))))->Fill(centralityBin, np[idx][iPid][1]); 
    (static_cast<TH2F*>(list->FindObject(Form("fHistDist%s%s%s",name, AliEbyEPidRatioHelper::fgkPidName[iPid], AliEbyEPidRatioHelper::fgkNetHistName[2]))))->Fill(centralityBin,                  np[idx][iPid][0]); 
    (static_cast<TH2F*>(list->FindObject(Form("fHistDist%s%s%s",name, AliEbyEPidRatioHelper::fgkPidName[iPid], AliEbyEPidRatioHelper::fgkNetHistName[3]))))->Fill(centralityBin, np[idx][iPid][1]-np[idx][iPid][0]); 
  }
  

  return;
}


//________________________________________________________________________
void  AliEbyEPidRatioPhy::AddHistSetCentPt(const Char_t *name, const Char_t *title)  {
  TString sName(name);
  TString sTitle(title);
  
  Float_t etaRange[2];
  fESDTrackCuts->GetEtaRange(etaRange[0],etaRange[1]);

  //TList *list[4];
  fOutList->Add(new TList);
  TList *list =  static_cast<TList*>(fOutList->Last());
  list->SetName(Form("f%s", name));
  list->SetOwner(kTRUE);
  

  Int_t nBinsCent         =  AliEbyEPidRatioHelper::fgkfHistNBinsCent;
  Double_t centBinRange[] = {AliEbyEPidRatioHelper::fgkfHistRangeCent[0], AliEbyEPidRatioHelper::fgkfHistRangeCent[1]};
  Int_t nBinsPt           =  AliEbyEPidRatioHelper::fgkfHistNBinsPt;
  Double_t ptBinRange[]   = {-0.5, AliEbyEPidRatioHelper::fgkfHistNBinsPt - 0.5};
  
  for (Int_t iPid = 0; iPid < 4; ++iPid) {

    TString sNetTitle(Form("%s - %s", AliEbyEPidRatioHelper::fgkPidLatex[iPid][1], AliEbyEPidRatioHelper::fgkPidLatex[iPid][0]));
    sTitle = (iPid != 0 ) ? Form("|y| < %.1f", fHelper->GetRapidityMax()) : Form(" |#eta|<%.1f", etaRange[1]);

    for (Int_t idx = 1; idx <= fOrder; ++idx) {
      list->Add(new TProfile2D(Form("fProf%s%sNetPt%dM", AliEbyEPidRatioHelper::fgkPidName[iPid],name, idx), 
			       Form("(%s)^{%d} Pt: %s;Centrality(100);(%s)^{%d}",sNetTitle.Data(), idx, sTitle.Data(), sNetTitle.Data(), idx),
			       100,-0.5,99.5, nBinsPt, ptBinRange[0], ptBinRange[1]));
    }
    
    for (Int_t ii = 0; ii <= fOrder; ++ii) {
      for (Int_t kk = 0; kk <= fOrder; ++kk) {
	list->Add(new TProfile2D(Form("fProf%s%sNetPtF%02d%02d", AliEbyEPidRatioHelper::fgkPidName[iPid], name, ii, kk),
				 Form("f_{%02d%02d} Pt: %s;Centrality(100);f_{%02d%02d}", ii, kk, sTitle.Data(), ii, kk),
				 100,-0.5,99.5, nBinsPt, ptBinRange[0], ptBinRange[1]));
      }
    }
    
    
    for (Int_t idx = 1; idx <= fOrder; ++idx) {
      list->Add(new TProfile2D(Form("fProfBin%s%sNetPt%dM", AliEbyEPidRatioHelper::fgkPidName[iPid],name, idx), 
			       Form("(%s)^{%d} Pt: %s;Centrality(11);(%s)^{%d}", sNetTitle.Data(), idx, sTitle.Data(), sNetTitle.Data(), idx),
			       nBinsCent, centBinRange[0], centBinRange[1],nBinsPt, ptBinRange[0], ptBinRange[1]));
    }
    
    for (Int_t ii = 0; ii <= fOrder; ++ii) {
      for (Int_t kk = 0; kk <= fOrder; ++kk) {
	list->Add(new TProfile2D(Form("fProfBin%s%sNetPtF%02d%02d", AliEbyEPidRatioHelper::fgkPidName[iPid], name, ii, kk),
				 Form("f_{%02d%02d} Pt : %s;Centrality(11);f_{%02d%02d}", ii, kk, sTitle.Data(), ii, kk),
				 nBinsCent, centBinRange[0], centBinRange[1],nBinsPt, ptBinRange[0], ptBinRange[1]));
      }
    }
  } // iPid  
  

  for (Int_t iPhy = 0; iPhy < 46; ++iPhy) { 
    list->Add(new TProfile2D(Form("fProf%sNuPt%02d",name,iPhy),
			     Form("Physics Variable for index %d | %s in PtBin ; Centrality;",iPhy,name),
			     100,-0.5,99.5, nBinsPt, ptBinRange[0], ptBinRange[1]));
  }
  for (Int_t iPhy = 0; iPhy < 46; ++iPhy) { 
    list->Add(new TProfile2D(Form("fProfBin%sNuPt%02d",name,iPhy),
			     Form("Physics Variable for index %d | %s in Pt Bin; Centrality;",iPhy,name),
			     nBinsCent, centBinRange[0], centBinRange[1],nBinsPt, ptBinRange[0], ptBinRange[1]));
  }
  

  return;
}


//________________________________________________________________________
void AliEbyEPidRatioPhy::FillHistSetCentPt(const Char_t *name, Int_t idx, Bool_t isMC)  {
 
  Int_t ****np = (isMC) ? fMCNpPt : fNpPt;
  
  Float_t centralityBin = fHelper->GetCentralityBin();
  Float_t centralityPer = fHelper->GetCentralityPercentile();
  
  TList *list = static_cast<TList*>(fOutList->FindObject(Form("f%s",name)));
  
  
  for (Int_t idxPt  = 0; idxPt < AliEbyEPidRatioHelper::fgkfHistNBinsPt; ++idxPt) {
    for (Int_t iPid = 0; iPid < 4; ++iPid) {
      Int_t deltaNp = np[idx][iPid][1][idxPt]-np[idx][iPid][0][idxPt];  
      Double_t delta = 1.;
      for (Int_t idxOrder = 1; idxOrder <= fOrder; ++idxOrder) {
	delta *= deltaNp;
	(static_cast<TProfile2D*>(list->FindObject(Form("fProfBin%s%sNetPt%dM", AliEbyEPidRatioHelper::fgkPidName[iPid], name, idxOrder))))->Fill(centralityBin,idxPt, delta);
	(static_cast<TProfile2D*>(list->FindObject(Form("fProf%s%sNetPt%dM", AliEbyEPidRatioHelper::fgkPidName[iPid], name, idxOrder))))->Fill(centralityPer,idxPt, delta);
      }
    
      for (Int_t idxOrder = 0; idxOrder <= fOrder; ++ idxOrder) {
	fRedFactp[idxOrder][0]  = 1.;
	fRedFactp[idxOrder][1]  = 1.;
      }
      
      for (Int_t idxOrder = 1; idxOrder <= fOrder; ++ idxOrder) {
	fRedFactp[idxOrder][0]  = fRedFactp[idxOrder-1][0] * Double_t(np[idx][iPid][0][idxPt]-(idxOrder-1));
	fRedFactp[idxOrder][1]  = fRedFactp[idxOrder-1][1] * Double_t(np[idx][iPid][1][idxPt]-(idxOrder-1));
      }
      
      for (Int_t ii = 0; ii <= fOrder; ++ii) {   
	for (Int_t kk = 0; kk <= fOrder; ++kk) { 
	  Double_t fik = fRedFactp[ii][1] * fRedFactp[kk][0];   // n1 *n2 -> p * pbar
	  (static_cast<TProfile2D*>(list->FindObject(Form("fProfBin%s%sNetPtF%02d%02d", AliEbyEPidRatioHelper::fgkPidName[iPid], name, ii, kk))))->Fill(centralityBin,idxPt, fik);
	  (static_cast<TProfile2D*>(list->FindObject(Form("fProf%s%sNetPtF%02d%02d", AliEbyEPidRatioHelper::fgkPidName[iPid], name, ii, kk))))->Fill(centralityPer,idxPt, fik);
	}
      }
    }


    Int_t a[6][4]; Int_t b[22];
    for (Int_t iPid = 0; iPid < 4; ++iPid) {
      a[0][iPid] = np[idx][iPid][1][idxPt]+np[idx][iPid][0][idxPt];       // 0  n+ + n-
      a[1][iPid] = np[idx][iPid][1][idxPt];                        // 1  n+
      a[2][iPid] = np[idx][iPid][0][idxPt];                        // 2  n-
      a[3][iPid] = np[idx][iPid][1][idxPt]*np[idx][iPid][0][idxPt];       // 3  n+ . n-
      a[4][iPid] = np[idx][iPid][1][idxPt]*(np[idx][iPid][1][idxPt]-1);   // 4  n+ (n+ - 1)
      a[5][iPid] = np[idx][iPid][0][idxPt]*(np[idx][iPid][0][idxPt]-1);   // 5  n- (n- - 1)
      
      // Printf("%6d %20s %6.2f %6d %6d %6d ", idx, name, centralityBin,
      //	   a[0][iPid], a[1][iPid], a[2][iPid]);
      
    }
    
    b[0]  = a[0][0]*a[0][2];       // 24 N   K
    b[1]  = a[0][1]*a[0][2];       // 25 Pi  K
    b[2]  = a[1][1]*a[1][2];       // 26 pi+ k+
    b[3]  = a[1][1]*a[2][2];       // 27 pi+ k-
    b[4]  = a[2][1]*a[1][2];       // 28 pi- k+  
    b[5]  = a[2][1]*a[2][2];       // 29 pi- k-
    
    b[6]  = a[0][0]*a[0][3];       // 30 N   P
    b[7]  = a[0][2]*a[0][3];       // 31 K   P
    b[8]  = a[1][2]*a[1][3];       // 32 k+  p+
    b[9]  = a[1][2]*a[2][3];       // 33 k+  p-
    b[10] = a[2][2]*a[1][3];       // 34 k-  p+
    b[11] = a[2][2]*a[2][3];       // 35 k-  p-
    
    b[12] = a[0][0]*a[0][1];       // 36 N  Pi
    b[13] = a[0][3]*a[0][1];       // 37 P  Pi
    b[14] = a[1][3]*a[1][1];       // 38 p+ pi+
    b[15] = a[1][3]*a[2][1];       // 39 p+ pi-
    b[16] = a[2][3]*a[1][1];       // 40 p- pi+
    b[17] = a[2][3]*a[2][1];       // 41 p- pi-
    
    b[18] = a[0][0]*(a[0][0] - 1); // 42 N ( N - 1 )
    b[19] = a[0][1]*(a[0][1] - 1); // 43 Pi( Pi- 1 )
    b[20] = a[0][2]*(a[0][1] - 1); // 44 K ( K - 1 )
    b[21] = a[0][3]*(a[0][3] - 1); // 45 P ( P - 1 )
    
    Int_t k = 0;
    for (Int_t j = 0; j < 4; j++) {
      for (Int_t i = 0; i < 6; i++) {
	(static_cast<TProfile2D*>(list->FindObject(Form("fProfBin%sNuPt%02d", name,k))))->Fill(centralityBin,idxPt,a[i][j]); 
	(static_cast<TProfile2D*>(list->FindObject(Form("fProf%sNuPt%02d", name,k))))->Fill(centralityPer,idxPt,a[i][j]); 
	k++;
      }
    }
    
    for (Int_t j = 0; j < 22; j++) {
      (static_cast<TProfile2D*>(list->FindObject(Form("fProfBin%sNuPt%02d", name,j+24))))->Fill(centralityBin,idxPt, b[j]); 
      (static_cast<TProfile2D*>(list->FindObject(Form("fProf%sNuPt%02d", name,j+24))))->Fill(centralityPer,idxPt,b[j]); 
    }
    
  }//idxPt
  
  return;
}


//________________________________________________________________________
void  AliEbyEPidRatioPhy::AddHistInGroup(const Char_t *name, const Char_t *title, Int_t nSample)  {

  

  TString sName(name);
  TString sTitle(title);
  
 

  Float_t etaRange[2];
  fESDTrackCuts->GetEtaRange(etaRange[0],etaRange[1]);

  //TList *list[4];
  fOutList->Add(new TList);
  TList *list =  static_cast<TList*>(fOutList->Last());
  list->SetName(Form("f%s", name));
  list->SetOwner(kTRUE);
  

  TString tname = (fIsPer) ? Form("%sPer", name) : Form("%sBin", name);
  Int_t nBinsCent         =  (fIsPer) ? 100 : AliEbyEPidRatioHelper::fgkfHistNBinsCent;
  Double_t centBinRange[2];  
  centBinRange[0]  =  (fIsPer) ?  0   : AliEbyEPidRatioHelper::fgkfHistRangeCent[0];
  centBinRange[1]  =  (fIsPer) ?  100 : AliEbyEPidRatioHelper::fgkfHistRangeCent[1];
  Int_t nBinsPt           =  AliEbyEPidRatioHelper::fgkfHistNBinsPt;
  Double_t ptBinRange[]   = {-0.5, AliEbyEPidRatioHelper::fgkfHistNBinsPt - 0.5};

  for (Int_t iSub = 0; iSub <= nSample; ++iSub) {
    
    list->Add(new TList);
    TList *listSub = static_cast<TList*>(list->Last());
    listSub->SetName(Form("%s%02d",name, iSub));
    listSub->SetOwner(kTRUE);
    
    for (Int_t iPid = 0; iPid < 4; ++iPid) {
    
      TString sNetTitle(Form("%s - %s", AliEbyEPidRatioHelper::fgkPidLatex[iPid][1], AliEbyEPidRatioHelper::fgkPidLatex[iPid][0]));
      sTitle = (iPid != 0 ) ? Form("|y| < %.1f", fHelper->GetRapidityMax()) : Form(" |#eta|<%.1f", etaRange[1]);
      
      listSub->Add(new TProfile(Form("fProfS%02dTot%sPlus%s", iSub, AliEbyEPidRatioHelper::fgkPidName[iPid],tname.Data()), 
				Form("(%s) : %s;Centrality(11);(%s)",AliEbyEPidRatioHelper::fgkPidName[iPid], sTitle.Data(), sNetTitle.Data()),
				nBinsCent, centBinRange[0], centBinRange[1]));
      
      listSub->Add(new TProfile(Form("fProfS%02dTot%sMinus%s",iSub, AliEbyEPidRatioHelper::fgkPidName[iPid],tname.Data()), 
				Form("(%s) : %s;Centrality(11);(%s)",AliEbyEPidRatioHelper::fgkPidName[iPid], sTitle.Data(), sNetTitle.Data()),
				nBinsCent, centBinRange[0], centBinRange[1]));
      
      
      
      for (Int_t idx = 1; idx <= fOrder; ++idx) {
	listSub->Add(new TProfile(Form("fProfS%02d%s%sNet%dM",iSub, AliEbyEPidRatioHelper::fgkPidName[iPid],tname.Data(), idx), 
				  Form("(%s)^{%d} : %s;Centrality(11);(%s)^{%d}", sNetTitle.Data(), idx, sTitle.Data(), sNetTitle.Data(), idx),
				  nBinsCent, centBinRange[0], centBinRange[1]));
      }
      
      for (Int_t ii = 0; ii <= fOrder; ++ii) {
	for (Int_t kk = 0; kk <= fOrder; ++kk) {
	  listSub->Add(new TProfile(Form("fProfS%02d%s%sNetF%02d%02d",iSub, AliEbyEPidRatioHelper::fgkPidName[iPid], tname.Data(), ii, kk),
				    Form("f_{%02d%02d} : %s;Centrality(11);f_{%02d%02d}", ii, kk, sTitle.Data(), ii, kk),
				    nBinsCent, centBinRange[0], centBinRange[1]));
	}
      }
      
      //-------------------------------------------
      if (fIsPtBin) {
	
	for (Int_t idx = 1; idx <= fOrder; ++idx) {
	  listSub->Add(new TProfile2D(Form("fProfS%02d%s%sNetPt%dM",iSub, AliEbyEPidRatioHelper::fgkPidName[iPid],tname.Data(), idx), 
				      Form("(%s)^{%d} Pt: %s;Centrality(11);(%s)^{%d}", sNetTitle.Data(), idx, sTitle.Data(), sNetTitle.Data(), idx),
				      nBinsCent, centBinRange[0], centBinRange[1],nBinsPt, ptBinRange[0], ptBinRange[1]));
	}
	
	for (Int_t ii = 0; ii <= fOrder; ++ii) {
	  for (Int_t kk = 0; kk <= fOrder; ++kk) {
	    listSub->Add(new TProfile2D(Form("fProfS%02d%s%sNetPtF%02d%02d",iSub, AliEbyEPidRatioHelper::fgkPidName[iPid], tname.Data(), ii, kk),
					Form("f_{%02d%02d} Pt : %s;Centrality(11);f_{%02d%02d}", ii, kk, sTitle.Data(), ii, kk),
					nBinsCent, centBinRange[0], centBinRange[1],nBinsPt, ptBinRange[0], ptBinRange[1]));
	  }
	}
      } // ptbin
      
    }  // iPid
    //-------------------------------------------
       
    for (Int_t iPhy = 0; iPhy < 46; ++iPhy) { 
      listSub->Add(new TProfile(Form("fProfS%02d%sNu%02d",iSub,tname.Data(),iPhy),
				Form("Physics Variable for index %d | %s | Sub S%02d; Centrality;",iPhy,tname.Data(), iSub),
				nBinsCent, centBinRange[0], centBinRange[1]));
    }
    
    
    if (fIsPtBin) {
      for (Int_t iPhy = 0; iPhy < 46; ++iPhy) { 
	listSub->Add(new TProfile2D(Form("fProfS%02d%sNuPt%02d",iSub,tname.Data(),iPhy),
				    Form("Physics Variable for index %d | %s in Pt Bin; Centrality;",iPhy,tname.Data()),
				    nBinsCent, centBinRange[0], centBinRange[1],nBinsPt, ptBinRange[0], ptBinRange[1]));
      }
    }
    
  }// isub
  
  return;
}

//________________________________________________________________________
void AliEbyEPidRatioPhy::FillHistInGroup(const Char_t *name, Int_t idx, Int_t iSub, Bool_t isMC)  {

  Float_t centralityBin = (fIsPer) ? fHelper->GetCentralityPercentile() : fHelper->GetCentralityBin();
  TString tname = (fIsPer) ? Form("%sPer", name) : Form("%sBin", name);
  
  TList *list    = static_cast<TList*>(fOutList->FindObject(Form("f%s",name)));
  TList *listSub = static_cast<TList*>(list->FindObject(Form("%s%02d",name, iSub)));

  
  if (fIsPtBin) {

    Int_t ****nppt = (isMC) ? fMCNpPt : fNpPt;
    
    for (Int_t idxPt  = 0; idxPt < AliEbyEPidRatioHelper::fgkfHistNBinsPt; ++idxPt) {
      for (Int_t iPid = 0; iPid < 4; ++iPid) {
	Int_t deltaNp = nppt[idx][iPid][1][idxPt]-nppt[idx][iPid][0][idxPt];  
	Double_t delta = 1.;
	for (Int_t idxOrder = 1; idxOrder <= fOrder; ++idxOrder) {
	  delta *= deltaNp;
	  (static_cast<TProfile2D*>(listSub->FindObject(Form("fProfS%02d%s%sNetPt%dM",iSub, AliEbyEPidRatioHelper::fgkPidName[iPid], tname.Data(), idxOrder))))->Fill(centralityBin,idxPt, delta);
	}
	
	for (Int_t idxOrder = 0; idxOrder <= fOrder; ++ idxOrder) {
	  fRedFactp[idxOrder][0]  = 1.;
	  fRedFactp[idxOrder][1]  = 1.;
	}
	

	for (Int_t idxOrder = 1; idxOrder <= fOrder; ++ idxOrder) {
	  fRedFactp[idxOrder][0]  = fRedFactp[idxOrder-1][0] * Double_t(nppt[idx][iPid][0][idxPt]-(idxOrder-1));
	  fRedFactp[idxOrder][1]  = fRedFactp[idxOrder-1][1] * Double_t(nppt[idx][iPid][1][idxPt]-(idxOrder-1));
	}
	
	for (Int_t ii = 0; ii <= fOrder; ++ii) {   
	  for (Int_t kk = 0; kk <= fOrder; ++kk) { 
	    Double_t fik = fRedFactp[ii][1] * fRedFactp[kk][0];   // n1 *n2 -> p * pbar
	    (static_cast<TProfile2D*>(listSub->FindObject(Form("fProfS%02d%s%sNetPtF%02d%02d",iSub, AliEbyEPidRatioHelper::fgkPidName[iPid], tname.Data(), ii, kk))))->Fill(centralityBin,idxPt, fik);
	  }
	}
      }
      
      
      Int_t a[6][4]; Int_t b[22];
      for (Int_t iPid = 0; iPid < 4; ++iPid) {
	a[0][iPid] = nppt[idx][iPid][1][idxPt]+nppt[idx][iPid][0][idxPt];       // 0  n+ + n-
	a[1][iPid] = nppt[idx][iPid][1][idxPt];                        // 1  n+
	a[2][iPid] = nppt[idx][iPid][0][idxPt];                        // 2  n-
	a[3][iPid] = nppt[idx][iPid][1][idxPt]*nppt[idx][iPid][0][idxPt];       // 3  n+ . n-
	a[4][iPid] = nppt[idx][iPid][1][idxPt]*(nppt[idx][iPid][1][idxPt]-1);   // 4  n+ (n+ - 1)
	a[5][iPid] = nppt[idx][iPid][0][idxPt]*(nppt[idx][iPid][0][idxPt]-1);   // 5  n- (n- - 1)
	
	// Printf("%6d %20s %6.2f %6d %6d %6d ",iSub, idx, name, centralityBin,
	//	   a[0][iPid], a[1][iPid], a[2][iPid]);
	
      }
      
      b[0]  = a[0][0]*a[0][2];       // 24 N   K
      b[1]  = a[0][1]*a[0][2];       // 25 Pi  K
      b[2]  = a[1][1]*a[1][2];       // 26 pi+ k+
      b[3]  = a[1][1]*a[2][2];       // 27 pi+ k-
      b[4]  = a[2][1]*a[1][2];       // 28 pi- k+  
      b[5]  = a[2][1]*a[2][2];       // 29 pi- k-
      
      b[6]  = a[0][0]*a[0][3];       // 30 N   P
      b[7]  = a[0][2]*a[0][3];       // 31 K   P
      b[8]  = a[1][2]*a[1][3];       // 32 k+  p+
      b[9]  = a[1][2]*a[2][3];       // 33 k+  p-
      b[10] = a[2][2]*a[1][3];       // 34 k-  p+
      b[11] = a[2][2]*a[2][3];       // 35 k-  p-
    
      b[12] = a[0][0]*a[0][1];       // 36 N  Pi
      b[13] = a[0][3]*a[0][1];       // 37 P  Pi
      b[14] = a[1][3]*a[1][1];       // 38 p+ pi+
      b[15] = a[1][3]*a[2][1];       // 39 p+ pi-
      b[16] = a[2][3]*a[1][1];       // 40 p- pi+
      b[17] = a[2][3]*a[2][1];       // 41 p- pi-
      
      b[18] = a[0][0]*(a[0][0] - 1); // 42 N ( N - 1 )
      b[19] = a[0][1]*(a[0][1] - 1); // 43 Pi( Pi- 1 )
      b[20] = a[0][2]*(a[0][1] - 1); // 44 K ( K - 1 )
      b[21] = a[0][3]*(a[0][3] - 1); // 45 P ( P - 1 )
    
      Int_t k = 0;
      for (Int_t j = 0; j < 4; j++) {
	for (Int_t i = 0; i < 6; i++) {
	  (static_cast<TProfile2D*>(listSub->FindObject(Form("fProfS%02d%sNuPt%02d",iSub, tname.Data(),k))))->Fill(centralityBin,idxPt,a[i][j]); 
	  k++;
	}
      }
      
      for (Int_t j = 0; j < 22; j++) {
	(static_cast<TProfile2D*>(listSub->FindObject(Form("fProfS%02d%sNuPt%02d",iSub, tname.Data(),j+24))))->Fill(centralityBin,idxPt, b[j]); 
      }
      
    }//idxPt
  } // fIsPtBin

  //------------------------------------------------------------

  Int_t ***np = (isMC) ? fMCNp : fNp;
      
  for (Int_t iPid = 0; iPid < 4; ++iPid) {
    Int_t deltaNp = np[idx][iPid][1]-np[idx][iPid][0];  
    Double_t delta = 1.;
    
    (static_cast<TProfile*>(listSub->FindObject(Form("fProfS%02dTot%sPlus%s",iSub, AliEbyEPidRatioHelper::fgkPidName[iPid], tname.Data()))))->Fill(centralityBin, np[idx][iPid][1]);
    (static_cast<TProfile*>(listSub->FindObject(Form("fProfS%02dTot%sMinus%s",iSub, AliEbyEPidRatioHelper::fgkPidName[iPid], tname.Data()))))->Fill(centralityBin, np[idx][iPid][0]);
    
    for (Int_t idxOrder = 1; idxOrder <= fOrder; ++idxOrder) {
      delta *= deltaNp;
      (static_cast<TProfile*>(listSub->FindObject(Form("fProfS%02d%s%sNet%dM",iSub, AliEbyEPidRatioHelper::fgkPidName[iPid], tname.Data(), idxOrder))))->Fill(centralityBin, delta);
    }
    
    for (Int_t idxOrder = 0; idxOrder <= fOrder; ++ idxOrder) {
      fRedFactp[idxOrder][0]  = 1.;
      fRedFactp[idxOrder][1]  = 1.;
    }
    
    for (Int_t idxOrder = 1; idxOrder <= fOrder; ++ idxOrder) {
      fRedFactp[idxOrder][0]  = fRedFactp[idxOrder-1][0]  * Double_t(np[idx][iPid][0]-(idxOrder-1));
      fRedFactp[idxOrder][1]  = fRedFactp[idxOrder-1][1]  * Double_t(np[idx][iPid][1]-(idxOrder-1));
    }
    
    for (Int_t ii = 0; ii <= fOrder; ++ii) {  
      for (Int_t kk = 0; kk <= fOrder; ++kk) { 
	Double_t fik = fRedFactp[ii][1] * fRedFactp[kk][0];   
       	(static_cast<TProfile*>(listSub->FindObject(Form("fProfS%02d%s%sNetF%02d%02d",iSub, AliEbyEPidRatioHelper::fgkPidName[iPid], tname.Data(), ii, kk))))->Fill(centralityBin, fik);
      }
    }
  }
 
  //Printf("%6d %20s %6.2f %6d %6d %6d %6d  %6d %6d %6d %6d", idx, name, centralityBin,
  //	 np[idx][0][1],  np[idx][0][0], 
  //	 np[idx][1][1],  np[idx][1][0], 
  ///	 np[idx][2][1],  np[idx][2][0], 
  //	 np[idx][3][1],  np[idx][3][0]);
  //

   Int_t a[6][4]; Int_t b[22];
   for (Int_t iPid = 0; iPid < 4; ++iPid) {
     a[0][iPid] = np[idx][iPid][1]+np[idx][iPid][0];       // 0  n+ + n-
     a[1][iPid] = np[idx][iPid][1];                        // 1  n+
     a[2][iPid] = np[idx][iPid][0];                        // 2  n-
     a[3][iPid] = np[idx][iPid][1]*np[idx][iPid][0];       // 3  n+ . n-
     a[4][iPid] = np[idx][iPid][1]*(np[idx][iPid][1]-1);   // 4  n+ (n+ - 1)
     a[5][iPid] = np[idx][iPid][0]*(np[idx][iPid][0]-1);   // 5  n- (n- - 1)
     
     // Printf("%6d %20s %6.2f %6d %6d %6d ", idx, name, centralityBin,
     //	   a[0][iPid], a[1][iPid], a[2][iPid]);

  }
  
  b[0]  = a[0][0]*a[0][2];       // 24 N   K
  b[1]  = a[0][1]*a[0][2];       // 25 Pi  K
  b[2]  = a[1][1]*a[1][2];       // 26 pi+ k+
  b[3]  = a[1][1]*a[2][2];       // 27 pi+ k-
  b[4]  = a[2][1]*a[1][2];       // 28 pi- k+  
  b[5]  = a[2][1]*a[2][2];       // 29 pi- k-
  
  b[6]  = a[0][0]*a[0][3];       // 30 N   P
  b[7]  = a[0][2]*a[0][3];       // 31 K   P
  b[8]  = a[1][2]*a[1][3];       // 32 k+  p+
  b[9]  = a[1][2]*a[2][3];       // 33 k+  p-
  b[10] = a[2][2]*a[1][3];       // 34 k-  p+
  b[11] = a[2][2]*a[2][3];       // 35 k-  p-
  
  b[12] = a[0][0]*a[0][1];       // 36 N  Pi
  b[13] = a[0][3]*a[0][1];       // 37 P  Pi
  b[14] = a[1][3]*a[1][1];       // 38 p+ pi+
  b[15] = a[1][3]*a[2][1];       // 39 p+ pi-
  b[16] = a[2][3]*a[1][1];       // 40 p- pi+
  b[17] = a[2][3]*a[2][1];       // 41 p- pi-
  
  b[18] = a[0][0]*(a[0][0] - 1); // 42 N ( N - 1 )
  b[19] = a[0][1]*(a[0][1] - 1); // 43 Pi( Pi- 1 )
  b[20] = a[0][2]*(a[0][1] - 1); // 44 K ( K - 1 )
  b[21] = a[0][3]*(a[0][3] - 1); // 45 P ( P - 1 )
  // TList *list_nu = static_cast<TList*>(fOutlistSub->FindObject(Form("f%s_nu",name)));
  Int_t k = 0;
  for (Int_t j = 0; j < 4; j++) {
    for (Int_t i = 0; i < 6; i++) {
      (static_cast<TProfile*>(listSub->FindObject(Form("fProfS%02d%sNu%02d",iSub, tname.Data(),k))))->Fill(centralityBin,a[i][j]); 
       k++;
    }
  }

  for (Int_t j = 0; j < 22; j++) {
    (static_cast<TProfile*>(listSub->FindObject(Form("fProfS%02d%sNu%02d",iSub, tname.Data(),j+24))))->Fill(centralityBin,b[j]); 
  }


  return;
}



