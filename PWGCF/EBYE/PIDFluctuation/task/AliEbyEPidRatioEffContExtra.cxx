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
//        Origin: Jochen Thaeder <jochen@thaeder.de>
//                Michael Weber <m.weber@cern.ch>
//=========================================================================//

#include "TMath.h"
#include "TAxis.h"

#include "AliESDEvent.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliESDtrackCuts.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"

#include "AliEbyEPidRatioEffContExtra.h"

using namespace std;

ClassImp(AliEbyEPidRatioEffContExtra)

//________________________________________________________________________
AliEbyEPidRatioEffContExtra::AliEbyEPidRatioEffContExtra() :
  AliEbyEPidRatioBase("EffCont", "EffCont"),
  fLabelsRec(NULL),
 
  fHnNchEMc(NULL),
  fHnNchERec(NULL),
  fHnNpiEMc(NULL),
  fHnNpiERec(NULL),
  fHnNkaEMc(NULL),
  fHnNkaERec(NULL),
  fHnNprEMc(NULL),
  fHnNprERec(NULL),

  fHnNchCMc(NULL),
  fHnNchCRec(NULL),
  fHnNpiCMc(NULL),
  fHnNpiCRec(NULL),
  fHnNkaCMc(NULL),
  fHnNkaCRec(NULL),
  fHnNprCMc(NULL),
  fHnNprCRec(NULL) {
  AliLog::SetClassDebugLevel("AliEbyEPidRatioEffContExtra",10);
}

//________________________________________________________________________
AliEbyEPidRatioEffContExtra::~AliEbyEPidRatioEffContExtra() {
  // Destructor


  for (Int_t ii = 0; ii < 2; ++ii) {
    for (Int_t kk = 0; kk < 4; ++kk)
      if (fLabelsRec[ii][kk]) delete[] fLabelsRec[ii][kk];
    if (fLabelsRec[ii]) delete[] fLabelsRec[ii];
  }
  if (fLabelsRec) delete[] fLabelsRec;

}


//________________________________________________________________________
void AliEbyEPidRatioEffContExtra::Process() {
  for(Int_t i = 0; i < 4; i++)  {
    FillMCLabels(i);
    FillMCEffHist(i);
  }
  return;
}      

//________________________________________________________________________
void AliEbyEPidRatioEffContExtra::Init() {
  fLabelsRec = new Int_t**[2];
  for (Int_t ii = 0 ; ii < 2; ++ii) {
    fLabelsRec[ii] = new Int_t*[4];
    for (Int_t kk = 0 ; kk < 4; ++kk)
      fLabelsRec[ii][kk] = NULL;
  }
  Printf(" >>>> AliEbyEPidRatioEffContExtra - inside");
}

//________________________________________________________________________
void AliEbyEPidRatioEffContExtra::CreateHistograms() {
   Int_t    binHnEffPIDMC[8] = {AliEbyEPidRatioHelper::fgkfHistNBinsCent, 
			       AliEbyEPidRatioHelper::fgkfHistNBinsSign, 
			       2, 2, 2,
			       AliEbyEPidRatioHelper::fgkfHistNBinsRap,     
			       AliEbyEPidRatioHelper::fgkfHistNBinsPhi, 
			       35};
  
  Double_t minHnEffPIDMC[8] = {AliEbyEPidRatioHelper::fgkfHistRangeCent[0], 
			       AliEbyEPidRatioHelper::fgkfHistRangeSign[0], 
			       -0.5,     -0.5,     -0.5,
			       AliEbyEPidRatioHelper::fgkfHistRangeRap[0],  
			       AliEbyEPidRatioHelper::fgkfHistRangePhi[0], 
			       0.2};  
  
  Double_t maxHnEffPIDMC[8] = {AliEbyEPidRatioHelper::fgkfHistRangeCent[1], 
			       AliEbyEPidRatioHelper::fgkfHistRangeSign[1], 
			       1.5,      1.5,      1.5,
			       AliEbyEPidRatioHelper::fgkfHistRangeRap[1],  
			       AliEbyEPidRatioHelper::fgkfHistRangePhi[1], 
			       3.0};  
  
  

  TString titilemc       = "cent:signMC:findable:recStatus:pidStatus:yMC:phiMC:ptMC";
 
  TString tiltlelaxmc[8]  = {"Centrality", 
			     "sign", 
			     "findable",
			     "recStatus",
			     "recPid",
			     "#it{y}_{MC}", 
			     "#varphi_{MC} (rad)", 
			     "#it{p}_{T,MC} (GeV/#it{c})"
  };

 
  //eff
  fHnNpiEMc  = new THnSparseF("hmNpiEffMc",titilemc.Data(),8,binHnEffPIDMC,minHnEffPIDMC, maxHnEffPIDMC);
  fHnNkaEMc  = new THnSparseF("hmNkaEffMc",titilemc.Data(),8,binHnEffPIDMC,minHnEffPIDMC, maxHnEffPIDMC);
  fHnNprEMc  = new THnSparseF("hmNprEffMc",titilemc.Data(),8,binHnEffPIDMC,minHnEffPIDMC, maxHnEffPIDMC);
  fHnNchEMc  = new THnSparseF("hmNchEffMc",titilemc.Data(),8,binHnEffPIDMC,minHnEffPIDMC, maxHnEffPIDMC);

  for (Int_t i = 0; i < 8; i++) {  
    fHnNchEMc->GetAxis(i)->SetTitle(tiltlelaxmc[i].Data());
    fHnNpiEMc->GetAxis(i)->SetTitle(tiltlelaxmc[i].Data());
    fHnNkaEMc->GetAxis(i)->SetTitle(tiltlelaxmc[i].Data());
    fHnNprEMc->GetAxis(i)->SetTitle(tiltlelaxmc[i].Data());
  }
  
  


  Int_t    binHnEffPID[5] = {AliEbyEPidRatioHelper::fgkfHistNBinsCent, 
			     AliEbyEPidRatioHelper::fgkfHistNBinsSign, 
			     AliEbyEPidRatioHelper::fgkfHistNBinsRap, 
			     AliEbyEPidRatioHelper::fgkfHistNBinsPhi,     
			     35};

  Double_t minHnEffPID[5] = {AliEbyEPidRatioHelper::fgkfHistRangeCent[0], 
			     AliEbyEPidRatioHelper::fgkfHistRangeSign[0],			  
			     AliEbyEPidRatioHelper::fgkfHistRangeRap[0], 
			     AliEbyEPidRatioHelper::fgkfHistRangePhi[0],  
			     0.2};

  
  Double_t maxHnEffPID[5] = {AliEbyEPidRatioHelper::fgkfHistRangeCent[1], 
			     AliEbyEPidRatioHelper::fgkfHistRangeSign[1], 
			     AliEbyEPidRatioHelper::fgkfHistRangeRap[1], 
			     AliEbyEPidRatioHelper::fgkfHistRangePhi[1],  
			     3.0};


  TString titilerec      = "cent:signRec:yRec:phiRec:ptRec";
  TString tiltlelaxrec[5]  = {"Centrality", 
			      "sign", 
			      "#it{y}_{Rec}", 
			      "#varphi_{Rec} (rad)", 
			      "#it{p}_{T,Rec} (GeV/#it{c})"};

 
  fHnNpiERec = new THnSparseF("hmNpiEffRec",titilerec.Data(),5,binHnEffPID, minHnEffPID, maxHnEffPID);
  fHnNkaERec = new THnSparseF("hmNkaEffRec",titilerec.Data(),5,binHnEffPID, minHnEffPID, maxHnEffPID);
  fHnNprERec = new THnSparseF("hmNprEffRec",titilerec.Data(),5,binHnEffPID, minHnEffPID, maxHnEffPID);
  fHnNchERec = new THnSparseF("hmNchEffRec",titilerec.Data(),5,binHnEffPID,minHnEffPID, maxHnEffPID);

  

  for (Int_t i = 0; i < 5; i++) { 
    fHnNchERec->GetAxis(i)->SetTitle(tiltlelaxrec[i].Data());
    fHnNpiERec->GetAxis(i)->SetTitle(tiltlelaxrec[i].Data());
    fHnNkaERec->GetAxis(i)->SetTitle(tiltlelaxrec[i].Data());
    fHnNprERec->GetAxis(i)->SetTitle(tiltlelaxrec[i].Data());
  }  

  //----- cont

  //----
  Int_t    binHnCont[6] = {AliEbyEPidRatioHelper::fgkfHistNBinsCent, 
			   AliEbyEPidRatioHelper::fgkfHistNBinsSign, 
			   8,                                        
			   AliEbyEPidRatioHelper::fgkfHistNBinsRap,  
			   AliEbyEPidRatioHelper::fgkfHistNBinsPhi, 
			   35};  
  
  Double_t minHnCont[6] = {AliEbyEPidRatioHelper::fgkfHistRangeCent[0], 
			   AliEbyEPidRatioHelper::fgkfHistRangeSign[0], 
			   0.5,                                            
			   AliEbyEPidRatioHelper::fgkfHistRangeRap[0],  
			   AliEbyEPidRatioHelper::fgkfHistRangePhi[0], 
			   0.2};
  
  
  
  Double_t maxHnCont[6] = {AliEbyEPidRatioHelper::fgkfHistRangeCent[1], 
			   AliEbyEPidRatioHelper::fgkfHistRangeSign[1],
			   8.5,                        
			   AliEbyEPidRatioHelper::fgkfHistRangeRap[1],  
			   AliEbyEPidRatioHelper::fgkfHistRangePhi[1], 
			   3.0};   




  TString titilecont     = "cent:signMC:contStatus:yMC:phiMC:ptMC";
  TString tiltlelaxcont[6] 
    = {"Centrality","sign","contPart","#it{y}_{MC}","#varphi_{MC} (rad)","#it{p}_{T,MC} (GeV/#it{c})"};
  
  fHnNpiCMc  = new THnSparseF("hmNpiContMc",titilecont.Data(),6,binHnCont,minHnCont, maxHnCont);
  fHnNkaCMc  = new THnSparseF("hmNkaContMc",titilecont.Data(),6,binHnCont,minHnCont, maxHnCont);
  fHnNprCMc  = new THnSparseF("hmNprContMc",titilecont.Data(),6,binHnCont,minHnCont, maxHnCont);
  fHnNchCMc  = new THnSparseF("hmNchContMc",titilecont.Data(),6,binHnCont,minHnCont, maxHnCont);

 for (Int_t i = 0; i < 6; i++) {  
    fHnNchCMc->GetAxis(i)->SetTitle(tiltlelaxcont[i].Data());
    fHnNpiCMc->GetAxis(i)->SetTitle(tiltlelaxcont[i].Data());
    fHnNkaCMc->GetAxis(i)->SetTitle(tiltlelaxcont[i].Data());
    fHnNprCMc->GetAxis(i)->SetTitle(tiltlelaxcont[i].Data());
  }

 

 fHnNpiCRec = new THnSparseF("hmNpiContRec",titilerec.Data(),5,binHnEffPID,minHnEffPID, maxHnEffPID);
 fHnNkaCRec = new THnSparseF("hmNkaContRec",titilerec.Data(),5,binHnEffPID,minHnEffPID, maxHnEffPID);
 fHnNprCRec = new THnSparseF("hmNprContRec",titilerec.Data(),5,binHnEffPID,minHnEffPID, maxHnEffPID);
 fHnNchCRec = new THnSparseF("hmNchContRec",titilerec.Data(),5,binHnEffPID,minHnEffPID, maxHnEffPID);

 for (Int_t i = 0; i < 5; i++) {  
   fHnNchCRec->GetAxis(i)->SetTitle(tiltlelaxrec[i].Data());
   fHnNpiCRec->GetAxis(i)->SetTitle(tiltlelaxrec[i].Data());
   fHnNkaCRec->GetAxis(i)->SetTitle(tiltlelaxrec[i].Data());
   fHnNprCRec->GetAxis(i)->SetTitle(tiltlelaxrec[i].Data());
 }  

 
  fHelper->BinLogAxis(fHnNchERec,4);
  fHelper->BinLogAxis(fHnNpiERec,4);
  fHelper->BinLogAxis(fHnNkaERec,4);
  fHelper->BinLogAxis(fHnNprERec,4);
  
  fHelper->BinLogAxis(fHnNchEMc,7);
  fHelper->BinLogAxis(fHnNpiEMc,7);
  fHelper->BinLogAxis(fHnNkaEMc,7);
  fHelper->BinLogAxis(fHnNprEMc,7);

  fHelper->BinLogAxis(fHnNchCMc,5);
  fHelper->BinLogAxis(fHnNpiCMc,5);
  fHelper->BinLogAxis(fHnNkaCMc,5);
  fHelper->BinLogAxis(fHnNprCMc,5);

  fHelper->BinLogAxis(fHnNchCRec,4);
  fHelper->BinLogAxis(fHnNpiCRec,4);
  fHelper->BinLogAxis(fHnNkaCRec,4);
  fHelper->BinLogAxis(fHnNprCRec,4);


  return;
}

//________________________________________________________________________
Int_t AliEbyEPidRatioEffContExtra::Setup() {
  for(Int_t i = 0; i < 4; i++) {
    fLabelsRec[0][i] = new Int_t[fNTracks];
    if(!fLabelsRec[0][i]) {
      AliError("Cannot create fLabelsRec[0]");
      return -1;
    }
    
    fLabelsRec[1][i] = new Int_t[fNTracks];
    if(!fLabelsRec[1][i]) {
      AliError("Cannot create fLabelsRec[1] for PID");
      return -1;
    }
    
    for(Int_t ii = 0; ii < fNTracks; ++ii) {
      fLabelsRec[0][i][ii] = 0;
      fLabelsRec[1][i][ii] = 0;
    }
  }
 
  return 0;
}

//________________________________________________________________________
void AliEbyEPidRatioEffContExtra::Reset() {
  // -- Reset eventwise

  for(Int_t i = 0; i < 4; i++) {
    for (Int_t ii = 0; ii < 2 ; ++ii) {
      if (fLabelsRec[ii][i])
	delete[] fLabelsRec[ii][i];
      fLabelsRec[ii][i] = NULL;
    }
  }
}

//________________________________________________________________________
void AliEbyEPidRatioEffContExtra::FillMCLabels(Int_t ipid) {
  // Fill MC labels
  // Loop over ESD tracks and fill arrays with MC lables
  //  fLabelsRec[0] : all Tracks
  //  fLabelsRec[1] : all Tracks accepted by PID of TPC
  // Check every accepted track if correctly identified
  //  otherwise check for contamination

  // -- Get ranges for AOD particles
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

      // -- Check if in pT and eta range (is done in ESDTrackCuts for ESDs)
      if(!(track->Pt() > ptRange[0] && track->Pt() <= ptRange[1] && TMath::Abs(track->Eta()) <= etaRange[1]))
	continue;
    }
   
    // -- Check if accepted in rapidity window -- for identified particles
    //    for charged particles -- eta check is done in the trackcuts
    Double_t yP;
    if (fHelper->GetUsePID(ipid) && !fHelper->IsTrackAcceptedRapidity(track, yP, ipid))
      continue;

 
    if (!fHelper->IsTrackAcceptedDCA(track))
      continue;

    Int_t label  = TMath::Abs(track->GetLabel()); 
    
    // continue;
    // -- Fill Label of all reconstructed
    fLabelsRec[0][ipid][idxTrack] = label;

    // -- Check if accepted by PID from TPC or TPC+TOF
    Double_t pid[3];
    if (!fHelper->IsTrackAcceptedPID(track, pid, fHelper->GetParticleSpecies(ipid))) // check it
      continue;

   

    // -- Fill Label of all reconstructed && recPid_TPC+TOF    
    fLabelsRec[1][ipid][idxTrack] = label;    
    
    // -- Check for contamination and fill contamination THnSparse
    CheckContTrack(track, ipid);

  } // for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {

  return;
}

//________________________________________________________________________
void AliEbyEPidRatioEffContExtra::CheckContTrack(AliVTrack *track, Int_t ipid) {
  // Check if particle is contamination or correctly identified for ESDs and AODs
  // Check for missidentified primaries and secondaries
  // Fill contamination THnSparse

  Int_t label     = TMath::Abs(track->GetLabel()); 
  Float_t signRec = track->Charge();

  
  AliVParticle* particle = (fESD) ? fMCEvent->GetTrack(label) : static_cast<AliVParticle*>(fArrayMC->At(label));
  if (!particle)
    return;

  Bool_t isPhysicalPrimary = (fESD) ? fStack->IsPhysicalPrimary(label): (static_cast<AliAODMCParticle*>(particle))->IsPhysicalPrimary();

  // -- Check if correctly identified 
  //    > return if correctly identified -> all ok, no action neededin this method
  //    > if PID required check -> for the correct (signed pdgcode) particle
  //    > no PID just check for primary 
  if (fHelper->GetUsePID(ipid)) {
    if (particle->PdgCode() == (signRec*fHelper->GetPdg(ipid)))
      if (isPhysicalPrimary)
	return;
  }
  else {
    if (isPhysicalPrimary)
      return;
  }

  // -- Check if secondaries from material or weak decay
  Bool_t isSecondaryFromWeakDecay = (fESD) ? fStack->IsSecondaryFromWeakDecay(label) : (static_cast<AliAODMCParticle*>(particle))->IsSecondaryFromWeakDecay();
  Bool_t isSecondaryFromMaterial  = (fESD) ? fStack->IsSecondaryFromMaterial(label)  : (static_cast<AliAODMCParticle*>(particle))->IsSecondaryFromMaterial();

  // -- Get PDG Charge of contaminating particle
  Float_t signMC = 0.;
  if      (particle->Charge() == 0.) signMC =  0.;
  else if (particle->Charge() <  0.) signMC = -1.;	
  else if (particle->Charge() >  0.) signMC =  1.;	

  // -- Get contaminating particle
  Double_t contPart = 0;
  if        (isSecondaryFromWeakDecay)                contPart = 7; // probeParticle from WeakDecay
  else if   (isSecondaryFromMaterial)                 contPart = 8; // probeParticle from Material
  else {
    if      (TMath::Abs(particle->PdgCode()) ==  211) contPart = 1; // pion
    else if (TMath::Abs(particle->PdgCode()) ==  321) contPart = 2; // kaon
    else if (TMath::Abs(particle->PdgCode()) == 2212) contPart = 3; // proton
    else if (TMath::Abs(particle->PdgCode()) ==   11) contPart = 4; // electron
    else if (TMath::Abs(particle->PdgCode()) ==   13) contPart = 5; // muon
    else                                              contPart = 6; // other
  }
  
  // -- Get Reconstructed y
  //    yRec = y for identified particles | yRec = eta for charged particles
  Double_t yRec  = 0.;
  fHelper->IsTrackAcceptedRapidity(track, yRec, ipid); 

 
  Double_t yetapid = (ipid == 0 ) ? particle->Eta() : particle->Y();
  Double_t yeta    = (ipid == 0 ) ? track->Eta() : yRec;

  Double_t hnContMc[6]  = {fCentralityBin,signMC,contPart,yetapid,particle->Phi(),particle->Pt()};
  Double_t hnContRec[5] = {fCentralityBin,signRec, yeta,track->Phi(),track->Pt()};
  
  if (ipid == 0) { fHnNchCRec->Fill(hnContRec); fHnNchCMc->Fill(hnContMc); }
  else if (ipid == 1) { fHnNpiCRec->Fill(hnContRec); fHnNpiCMc->Fill(hnContMc); }
  else if (ipid == 2) { fHnNkaCRec->Fill(hnContRec); fHnNkaCMc->Fill(hnContMc); }
  else if (ipid == 3) { fHnNprCRec->Fill(hnContRec); fHnNprCMc->Fill(hnContMc); }
  
  
}

//________________________________________________________________________
void AliEbyEPidRatioEffContExtra::FillMCEffHist(Int_t ipid) {
  // Fill efficiency THnSparse for ESDs

  Float_t etaRange[2];
  fESDTrackCuts->GetEtaRange(etaRange[0],etaRange[1]);

  Int_t nPart  = (fESD) ? fStack->GetNprimary() : fArrayMC->GetEntriesFast();

  for (Int_t idxMC = 0; idxMC < nPart; ++idxMC) {
    AliVParticle* particle = (fESD) ? fMCEvent->GetTrack(idxMC) : static_cast<AliVParticle*>(fArrayMC->At(idxMC));

    // -- Check basic MC properties -> charged physical primary
    if (!fHelper->IsParticleAcceptedBasicCharged(particle, idxMC))
      continue;

    // -- Check if accepted in rapidity window -- for identified particles
    Double_t yMC;
    if (fHelper->GetUsePID(ipid) && !fHelper->IsParticleAcceptedRapidity(particle, yMC, ipid))
      continue;

    // -- Check if accepted in eta window -- for charged particles
    if (!fHelper->GetUsePID(ipid) && TMath::Abs(particle->Eta()) > etaRange[1])
      continue;

    // -- Check if probeParticle / anti-probeParticle 
    //    > skip check if PID is not required
    if (fHelper->GetUsePID(ipid) && TMath::Abs(particle->PdgCode()) != fHelper->GetPdg(ipid))
      continue;
    
    // -- Get sign of particle
    Float_t signMC    = (particle->PdgCode() < 0) ? -1. : 1.;

    // -- Get if particle is findable --- not availible for AODs yet
    Float_t findable  = (fESD) ? Float_t(fHelper->IsParticleFindable(idxMC)) : 1.;

    // -- Get recStatus and pidStatus
    Float_t recStatus = 0.;
    Float_t recPid    = 0.;

    // -- Get Reconstructed values 
    Float_t etaRec  = 0.;
    Float_t phiRec  = 0.;
    Float_t ptRec   = 0.;
    Double_t yRec   = 0.;
    Float_t signRec = 0.;

    // -- Loop over all labels
    for (Int_t idxRec=0; idxRec < fNTracks; ++idxRec) {
      if (idxMC == fLabelsRec[0][ipid][idxRec]) {
	recStatus = 1.;
	
	if (idxMC == fLabelsRec[1][ipid][idxRec]) recPid = 1.;
	
        AliVTrack *track = NULL;
        if(fESD)
          track = fESD->GetTrack(idxRec);
        else if(fAOD)
          track = fAOD->GetTrack(idxRec);
	
        if (track) {
	  etaRec  = track->Eta();
          phiRec  = track->Phi();         
          ptRec   = track->Pt();
	  signRec = track->Charge();
          fHelper->IsTrackAcceptedRapidity(track, yRec, ipid); 
	  Double_t yeta    = (ipid == 0 ) ? etaRec : yRec;
	  Double_t hneffRec[5] = {fCentralityBin,signRec, yeta,phiRec,ptRec};
	  if (ipid == 0) { fHnNchERec->Fill(hneffRec); }
	  else if (ipid == 1) { fHnNpiERec->Fill(hneffRec); }
	  else if (ipid == 2) { fHnNkaERec->Fill(hneffRec); }
	  else if (ipid == 3) { fHnNprERec->Fill(hneffRec); }

        }     
        break;
      }
    } // for (Int_t idxRec=0; idxRec < fNTracks; ++idxRec) {  
    /*
    Double_t deltaPhi = particle->Phi()-phiRec;
    if (TMath::Abs(deltaPhi) > TMath::TwoPi()) {
      if (deltaPhi < 0)
	deltaPhi += TMath::TwoPi();
      else
    	deltaPhi -= TMath::TwoPi();
    }
    */
  
    Double_t yetapid = (ipid == 0 ) ? particle->Eta() : particle->Y();
   

    // Printf("%2d  %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f", ipid, yetapid, yeta, particle->Eta(), particle->Y(),  etaRec, yRec);    

    Double_t hneffMc[8]  = {fCentralityBin,signMC,findable, recStatus, recPid,yetapid,particle->Phi(),particle->Pt()};
    

    if (ipid == 0) { fHnNchEMc->Fill(hneffMc); }
    else if (ipid == 1) { fHnNpiEMc->Fill(hneffMc); }
    else if (ipid == 2) { fHnNkaEMc->Fill(hneffMc); }
    else if (ipid == 3) { fHnNprEMc->Fill(hneffMc); }
   

  } // for (Int_t idxMC = 0; idxMC < nPart; ++idxMC) {
  
  return;
}
