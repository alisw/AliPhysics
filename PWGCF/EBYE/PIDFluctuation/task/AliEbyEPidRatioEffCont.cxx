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

#include "AliEbyEPidRatioEffCont.h"

using namespace std;


ClassImp(AliEbyEPidRatioEffCont)

//________________________________________________________________________
AliEbyEPidRatioEffCont::AliEbyEPidRatioEffCont() :
  AliEbyEPidRatioBase("EffCont", "EffCont"),
  fLabelsRec(NULL),
  fHnEffMc(NULL),
  fHnContMc(NULL),
  fHnEffRec(NULL),
  fHnContRec(NULL) {
  // Constructor   

  AliLog::SetClassDebugLevel("AliEbyEPidRatioEffCont",10);
}

//________________________________________________________________________
AliEbyEPidRatioEffCont::~AliEbyEPidRatioEffCont() {
  // Destructor

  for (Int_t ii = 0; ii < 2; ++ii) {
    for (Int_t kk = 0; kk < 4; ++kk)
      if (fLabelsRec[ii][kk]) delete[] fLabelsRec[ii][kk];
    if (fLabelsRec[ii]) delete[] fLabelsRec[ii];
  }
  if (fLabelsRec) delete[] fLabelsRec;


}

//________________________________________________________________________
void AliEbyEPidRatioEffCont::Process() {
  // -- Process event

  // -- Setup (clean, create and fill) MC labels
  FillMCLabels();
 
  // -- Fill  MC histograms for efficiency studies
  FillMCEffHist();

  return;
}      

//________________________________________________________________________
void AliEbyEPidRatioEffCont::Init() {
  // -- Init eventwise

  fLabelsRec = new Int_t**[2];
  for (Int_t ii = 0 ; ii < 2; ++ii) {
    fLabelsRec[ii] = new Int_t*[4];
    for (Int_t kk = 0 ; kk < 4; ++kk)
      fLabelsRec[ii][kk] = NULL;
  }
}

//________________________________________________________________________
void AliEbyEPidRatioEffCont::CreateHistograms() {
  // Copied from NetParticle class
  Int_t    binHnEff[10] = { AliEbyEPidRatioHelper::fgkfHistNBinsCent, 4,   
			    AliEbyEPidRatioHelper::fgkfHistNBinsSign, 2,      2  ,      2  ,                       
			    AliEbyEPidRatioHelper::fgkfHistNBinsEta,     
			    AliEbyEPidRatioHelper::fgkfHistNBinsRap,     
			    AliEbyEPidRatioHelper::fgkfHistNBinsPhi,     
			    AliEbyEPidRatioHelper::fgkfHistNBinsPt};
  
  Double_t minHnEff[10] = {AliEbyEPidRatioHelper::fgkfHistRangeCent[0], -0.5, 
			   AliEbyEPidRatioHelper::fgkfHistRangeSign[0], -0.5, -0.5, -0.5,  
			   AliEbyEPidRatioHelper::fgkfHistRangeEta[0], 
			   AliEbyEPidRatioHelper::fgkfHistRangeRap[0],  
			   AliEbyEPidRatioHelper::fgkfHistRangePhi[0], 
			   AliEbyEPidRatioHelper::fgkfHistRangePt[0]};


  Double_t maxHnEff[10] = {AliEbyEPidRatioHelper::fgkfHistRangeCent[1], 3.5,
			   AliEbyEPidRatioHelper::fgkfHistRangeSign[1], 1.5,      1.5,      1.5,
			   AliEbyEPidRatioHelper::fgkfHistRangeEta[1],  
			   AliEbyEPidRatioHelper::fgkfHistRangeRap[1],  
			   AliEbyEPidRatioHelper::fgkfHistRangePhi[1], 
			   AliEbyEPidRatioHelper::fgkfHistRangePt[1]};
  
  fHnEffMc    = new THnSparseF("hnEffMc", "cent:pid:SignMC:findable:recStatus:pidStatus:etaMC:yMC:phiMC:ptMC", 10, binHnEff, minHnEff, maxHnEff);
  fHnEffRec   = new THnSparseF("hnEffRec", "cent:pid:SignMC:findable:recStatus:pidStatus:etaRec:yRec:phiRec:ptRec", 10, binHnEff, minHnEff, maxHnEff);
  //fHnEffRecMc = new THnSparseF("hnEffRecMMc", "cent:pid:SignMC:findable:recStatus:pidStatus:deltaEta:deltaY:deltaPhi:deltaPt" 10, binHnEff, minHnEff, maxHnEff);

  fHnEffMc->Sumw2();    
  fHnEffRec->Sumw2();    
 
  

  fHnEffMc->GetAxis(0)->SetTitle("centrality");                   //  0-5|5-10|10-20|20-30|30-40|40-50|50-60|60-70|70-80|80-90 --> 10 bins
  fHnEffMc->GetAxis(1)->SetTitle("N_{ch}|N_{#pi}|N_{K}|N_{p}");                //  0 | 1 | 2 | 3
  fHnEffMc->GetAxis(2)->SetTitle("sign");                         //  -1 | 0 | +1 
  fHnEffMc->GetAxis(3)->SetTitle("findable");                     //  0 not findable      |  1 findable
  fHnEffMc->GetAxis(4)->SetTitle("recStatus");                    //  0 not reconstructed |  1 reconstructed
  fHnEffMc->GetAxis(5)->SetTitle("recPid");                       //  0 not accepted      |  1 accepted
  fHnEffMc->GetAxis(6)->SetTitle("#eta_{MC}");                    //  eta  [-0.9, 0.9]
  fHnEffMc->GetAxis(7)->SetTitle("#it{y}_{MC}");                  //  rapidity  [-0.5, 0.5]
  fHnEffMc->GetAxis(8)->SetTitle("#varphi_{MC} (rad)");           //  phi  [ 0. , 2Pi]
  fHnEffMc->GetAxis(9)->SetTitle("#it{p}_{T,MC} (GeV/#it{c})");   //  pT   [ 0.2, 2.3]
  

  fHnEffRec->GetAxis(0)->SetTitle("centrality");                   //  0-5|5-10|10-20|20-30|30-40|40-50|50-60|60-70|70-80|80-90 --> 10 bins
  fHnEffRec->GetAxis(1)->SetTitle("N_{ch}|N_{#pi}|N_{K}|N_{p}");                //  0 | 1 | 2 | 3
  fHnEffRec->GetAxis(2)->SetTitle("sign");                         //  -1 | 0 | +1 
  fHnEffRec->GetAxis(3)->SetTitle("findable");                     //  0 not findable      |  1 findable
  fHnEffRec->GetAxis(4)->SetTitle("recStatus");                    //  0 not reconstructed |  1 reconstructed
  fHnEffRec->GetAxis(5)->SetTitle("recPid");                       //  0 not accepted      |  1 accepted
  fHnEffRec->GetAxis(6)->SetTitle("#eta_{Rec}");                   //  eta  [-0.9, 0.9]
  fHnEffRec->GetAxis(7)->SetTitle("#it{y}_{Rec}");                 //  rapidity  [-0.5, 0.5]
  fHnEffRec->GetAxis(8)->SetTitle("#varphi_{Rec} (rad)");          //  phi  [ 0. , 2Pi]
  fHnEffRec->GetAxis(9)->SetTitle("#it{p}_{T,Rec} (GeV/#it{c})");  //  pt   [ 0.2, 2.3]
 
  fHelper->BinLogAxis(fHnEffMc, 9);
  fHelper->BinLogAxis(fHnEffRec, 9);

  /* 
     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Copied from NetParticle
     
     creation -> findable -> reconstructed -> pid_TPC+TOF
        (1)         (2)          (3)            (4)      
                                                                      ||   findable | recStatus | recPid 
     1) all primary probeParticles_MC                                 ||      -           -         -
     2) all findable primary probeParticles_MC                        ||      x           -         -
     3) all reconstructed primary probeParticles_MC                   ||      x           x         -
     4) all reconstructed primary probeParticles_MC & recPid_TPC+TOF  ||      x           x         x
     
     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  */ 

  // ------------------------------------------------------------------
  // -- Create THnSparse - Cont
  // ------------------------------------------------------------------

  Int_t    binHnCont[8] = {AliEbyEPidRatioHelper::fgkfHistNBinsCent, 4, 
			   AliEbyEPidRatioHelper::fgkfHistNBinsSign, 8,   
			   AliEbyEPidRatioHelper::fgkfHistNBinsEta,     
			   AliEbyEPidRatioHelper::fgkfHistNBinsRap,  
			   AliEbyEPidRatioHelper::fgkfHistNBinsPhi,     
			   AliEbyEPidRatioHelper::fgkfHistNBinsPt};   
  
  Double_t minHnCont[8] = {AliEbyEPidRatioHelper::fgkfHistRangeCent[0], -0.5,
			   AliEbyEPidRatioHelper::fgkfHistRangeSign[0],  0.5,
			   AliEbyEPidRatioHelper::fgkfHistRangeEta[0], 
			   AliEbyEPidRatioHelper::fgkfHistRangeRap[0],  
			   AliEbyEPidRatioHelper::fgkfHistRangePhi[0], 
			   AliEbyEPidRatioHelper::fgkfHistRangePt[0]};   
  
  Double_t maxHnCont[8] = {AliEbyEPidRatioHelper::fgkfHistRangeCent[1], 3.5,
			   AliEbyEPidRatioHelper::fgkfHistRangeSign[1], 8.5,
			   AliEbyEPidRatioHelper::fgkfHistRangeEta[1], 
			   AliEbyEPidRatioHelper::fgkfHistRangeRap[1],  
			   AliEbyEPidRatioHelper::fgkfHistRangePhi[1], 
			   AliEbyEPidRatioHelper::fgkfHistRangePt[1]};   
  
  fHnContMc  = new THnSparseF("hnContMc", "cent:pid:SignMC:contPart:etaMC:yMC:phiMC:ptMC",8, binHnCont, minHnCont, maxHnCont);
  fHnContRec = new THnSparseF("hnContRec", "cent:pid:SignRec:contPart:etaRec:yRec:phiRec:ptRec",8, binHnCont, minHnCont, maxHnCont);

  fHnContMc->Sumw2();    
  fHnContRec->Sumw2();    

  fHnContMc->GetAxis(0)->SetTitle("centrality");                   //  0-5|5-10|10-20|20-30|30-40|40-50|50-60|60-70|70-80|80-90 --> 10 bins
  fHnContMc->GetAxis(1)->SetTitle("N_{ch}|N_{#pi}|N_{K}|N_{p}");                //  0 | 1 | 2 | 3
  fHnContMc->GetAxis(2)->SetTitle("sign");                         //  -1 | 0 | +1  
  fHnContMc->GetAxis(3)->SetTitle("contPart");                     //  1 pi | 2 K | 3 p | 4 e | 5 mu | 6 other | 7 p from WeakDecay | 8 p from Material
  fHnContMc->GetAxis(4)->SetTitle("#eta_{MC}");                    //  eta  [-0.9,0.9]
  fHnContMc->GetAxis(5)->SetTitle("#it{y}_{MC}");                  //  rapidity  [-0.5, 0.5]
  fHnContMc->GetAxis(6)->SetTitle("#varphi_{MC} (rad)");           //  phi  [ 0. ,2Pi]
  fHnContMc->GetAxis(7)->SetTitle("#it{p}_{T,MC} (GeV/#it{c})");   //  pT   [ 0.2,2.3]
  
  
  fHnContRec->GetAxis(0)->SetTitle("centrality");                   //  0-5|5-10|10-20|20-30|30-40|40-50|50-60|60-70|70-80|80-90 --> 10 bins
  fHnContRec->GetAxis(1)->SetTitle("N_{ch}|N_{#pi}|N_{K}|N_{p}");                //  0 | 1 | 2 | 3
  fHnContRec->GetAxis(2)->SetTitle("sign");                         //  -1 | 0 | +1  
  fHnContRec->GetAxis(3)->SetTitle("contPart");                     //  1 pi | 2 K | 3 p | 4 e | 5 mu | 6 other | 7 p from WeakDecay | 8 p from Material
  fHnContRec->GetAxis(4)->SetTitle("#eta_{Rec}");                   //  eta  [-0.9, 0.9]
  fHnContRec->GetAxis(5)->SetTitle("#it{y}_{Rec}");                 //  rapidity  [-0.5, 0.5]
  fHnContRec->GetAxis(6)->SetTitle("#varphi_{Rec} (rad)");          //  phi  [ 0. , 2Pi]
  fHnContRec->GetAxis(7)->SetTitle("#it{p}_{T,Rec} (GeV/#it{c})"); //  pt   [ 0.2, 2.3]
 

  //  fHnCont->GetAxis(12)->SetTitle("#eta_{MC}-#eta_{Rec}");                      //  eta  [-0.9, 0.9]
  // fHnCont->GetAxis(13)->SetTitle("#it{y}_{MC}-#it{y}_{Rec}");                  //  rapidity  [-0.5, 0.5]
  // fHnCont->GetAxis(14)->SetTitle("#varphi_{MC}-#varphi_{Rec} (rad)");          //  phi  [ -2Pi , 2Pi]
  // fHnCont->GetAxis(15)->SetTitle("#it{p}_{T,MC}-#it{p}_{T,Rec} (GeV/#it{c})"); //  pt   [ -2.3, 2.3]
  // fHnCont->GetAxis(16)->SetTitle("sign_{MC}-sign_{Rec}");                      //  -2 | 0 | +2 
  // fHnCont->GetAxis(17)->SetTitle("N_{ch}|N_{#pi}|N_{K}|N_{p}");                //  0 | 1 | 2 | 3

  fHelper->BinLogAxis(fHnContMc,  7);
  fHelper->BinLogAxis(fHnContRec, 7);

  return;
}

//________________________________________________________________________
Int_t AliEbyEPidRatioEffCont::Setup() {
  // -- Setup eventwise

  // -- Create label arrays
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
void AliEbyEPidRatioEffCont::Reset() {
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
void AliEbyEPidRatioEffCont::FillMCLabels() {
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

    Int_t gPdgCode = 0;
    
    Int_t iPid = 0;
    Double_t pid[3];
    if      (fHelper->IsTrackAcceptedPID(track, pid, (AliPID::kPion)))  {  iPid = 1; gPdgCode = 211;}
    else if (fHelper->IsTrackAcceptedPID(track, pid, (AliPID::kKaon)))  {  iPid = 2; gPdgCode = 321;}
    else if (fHelper->IsTrackAcceptedPID(track, pid, (AliPID::kProton))){  iPid = 3; gPdgCode = 2212;}
    else iPid = 0;

    //  cout << " --- EFF ---- " << iPid << "  " << gPdgCode << endl;
    
    Double_t yP;
    if (!fHelper->IsTrackAcceptedRapidity(track, yP, iPid))
      continue;

    if (!fHelper->IsTrackAcceptedDCA(track))
      continue;

    Int_t label  = TMath::Abs(track->GetLabel()); 
    
    // -- Fill Label of all reconstructed
    if(iPid != 0) fLabelsRec[0][0][idxTrack]    = label;
    fLabelsRec[0][iPid][idxTrack] = label;

    // -- Fill Label of all reconstructed && recPid_TPC+TOF    
    if(iPid != 0) fLabelsRec[1][0][idxTrack]    = label;    
    fLabelsRec[1][iPid][idxTrack] = label;    
    
    // -- Check for contamination and fill contamination THnSparse
    CheckContTrack(track, iPid, gPdgCode);
    if(iPid != 0) CheckContTrack(track, 0, 0);

  } // for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {

  return;
}

//________________________________________________________________________
void AliEbyEPidRatioEffCont::CheckContTrack(AliVTrack *track, Int_t iPid, Int_t gPdgCode) {
  Int_t label     = TMath::Abs(track->GetLabel()); 
  Float_t signRec = track->Charge();

  
  AliVParticle* particle = (fESD) ? fMCEvent->GetTrack(label) : static_cast<AliVParticle*>(fArrayMC->At(label));
  if (!particle)
    return;

  Bool_t isPhysicalPrimary = (fESD) ? fStack->IsPhysicalPrimary(label): (static_cast<AliAODMCParticle*>(particle))->IsPhysicalPrimary();
  if (iPid == 0) {
    if (particle->PdgCode() == (signRec*gPdgCode))
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
  Float_t contPart = 0;
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
  
  // cout << " --- CONT ---- " << iPid << "  " << gPdgCode << endl;

  // -- Get Reconstructed y
  //    yRec = y for identified particles | yRec = eta for charged particles
  Double_t yRec  = 0.;
  fHelper->IsTrackAcceptedRapidity(track, yRec, iPid); 

  Double_t deltaPhi = particle->Phi()-track->Phi();
  if (TMath::Abs(deltaPhi) > TMath::TwoPi()) {
    if (deltaPhi < 0)
      deltaPhi += TMath::TwoPi();
    else
      deltaPhi -= TMath::TwoPi();
  }

  Double_t hnContMc[8]  = {fCentralityBin,static_cast<Double_t>(iPid),signMC,static_cast<Double_t>(contPart),particle->Eta(),particle->Y(),particle->Phi(),particle->Pt()};
  Double_t hnContRec[8] = {fCentralityBin,static_cast<Double_t>(iPid),signRec,static_cast<Double_t>(contPart), track->Eta(),yRec,track->Phi(),track->Pt()};
  fHnContMc->Fill(hnContMc);
  fHnContRec->Fill(hnContRec);
   
}

//________________________________________________________________________
void AliEbyEPidRatioEffCont::FillMCEffHist() {
  // Fill efficiency THnSparse for ESDs

  Float_t etaRange[2];
  fESDTrackCuts->GetEtaRange(etaRange[0],etaRange[1]);

  Int_t nPart  = (fESD) ? fStack->GetNprimary() : fArrayMC->GetEntriesFast();

  for (Int_t idxMC = 0; idxMC < nPart; ++idxMC) {
    AliVParticle* particle = (fESD) ? fMCEvent->GetTrack(idxMC) : static_cast<AliVParticle*>(fArrayMC->At(idxMC));

    // -- Check basic MC properties -> charged physical primary
    if (!fHelper->IsParticleAcceptedBasicCharged(particle, idxMC))
      continue;

    Int_t iPid = 0;
    Int_t gPdgCode = 0;
    if ( TMath::Abs(particle->PdgCode())      ==  211 ) {  iPid = 1; gPdgCode = 211;}
    else if ( TMath::Abs(particle->PdgCode()) ==  321 ) {  iPid = 2; gPdgCode = 321;}
    else if ( TMath::Abs(particle->PdgCode()) == 2212 ) {  iPid = 3; gPdgCode = 2212;}
    else {iPid = 0; gPdgCode = 0;}

    // -- Check if accepted in rapidity window -- for identified particles
    Double_t yMC;
    if (iPid != 0) {
      if(!fHelper->IsParticleAcceptedRapidity(particle, yMC, iPid))
	continue;
    } else {
      // -- Check if accepted in eta window -- for charged particles
      if (TMath::Abs(particle->Eta()) > etaRange[1])
	continue;
    }
  
    // cout << particle->PdgCode() << "  " <<iPid << endl;
    
    // -- Check if probeParticle / anti-probeParticle 
    //    > skip check if PID is not required
    if (iPid != 0) { 
      if (TMath::Abs(particle->PdgCode()) != gPdgCode)
      continue;
    }
    
    // -- Get sign of particle
    Float_t signMC    = (particle->PdgCode() < 0) ? -1. : 1.;

    // -- Get if particle is findable --- not availible for AODs yet
    Float_t findable  = (fESD) ? Float_t(fHelper->IsParticleFindable(idxMC)) : 1.;

    // cout << findable << "  " << fHelper->IsParticleFindable(idxMC)<< endl;

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
      if (idxMC == fLabelsRec[0][iPid][idxRec]) {
	recStatus = 1.;
	
	if (idxMC == fLabelsRec[1][iPid][idxRec])
	  recPid = 1.;
	
        AliVTrack *track = NULL;
        if(fESD)
          track = fESD->GetTrack(idxRec);
        else if(fAOD)
          track = fAOD->GetTrack(idxRec);
	
        if (track) {
          // if no track present (which should not happen)
          // -> pt = 0. , which is not in the looked at range
	  
          // -- Get Reconstructed values
          etaRec  = track->Eta();
          phiRec  = track->Phi();         
          ptRec   = track->Pt();
	  signRec = track->Charge();
          fHelper->IsTrackAcceptedRapidity(track, yRec, iPid); // yRec = y for identified particles | yRec = eta for charged particles
        }      
        break;
      }
    } // for (Int_t idxRec=0; idxRec < fNTracks; ++idxRec) {  

    Double_t deltaPhi = particle->Phi()-phiRec;
    if (TMath::Abs(deltaPhi) > TMath::TwoPi()) {
      if (deltaPhi < 0)
	deltaPhi += TMath::TwoPi();
      else
    	deltaPhi -= TMath::TwoPi();
    }
    
    // if (signRec == 0) continue;

    if(iPid != 0) {
      Double_t hnEffMc[10]  = {fCentralityBin,0,
			       static_cast<Double_t>(signMC),
			       static_cast<Double_t>(findable), 
			       static_cast<Double_t>(recStatus),
			       static_cast<Double_t>(recPid),
			       particle->Eta(), particle->Y(), particle->Phi(),particle->Pt()};
      Double_t hnEffRec[10] = {fCentralityBin,0,
			       static_cast<Double_t>(signRec),
			       static_cast<Double_t>(findable), 
			       static_cast<Double_t>(recStatus),
			       static_cast<Double_t>(recPid),etaRec, yRec, phiRec, ptRec};
      fHnEffMc->Fill(hnEffMc);
      fHnEffRec->Fill(hnEffRec);
    }
    Double_t hnEffMc[10]  = {fCentralityBin,static_cast<Double_t>(iPid),
			     static_cast<Double_t>(signMC),
			     static_cast<Double_t>(findable), 
			     static_cast<Double_t>(recStatus),
			     static_cast<Double_t>(recPid),
			     particle->Eta(), particle->Y(), particle->Phi(),particle->Pt()};
    Double_t hnEffRec[10] = {fCentralityBin,
			     static_cast<Double_t>(iPid),
			     static_cast<Double_t>(signRec),
			     static_cast<Double_t>(findable), 
			     static_cast<Double_t>(recStatus),
			     static_cast<Double_t>(recPid),etaRec, yRec, phiRec, ptRec};
    fHnEffMc->Fill(hnEffMc);
    fHnEffRec->Fill(hnEffRec);
    
   
    //  cout << signMC << "  " << signRec << "  " << iPid << "  " << gPdgCode << endl;


  } // for (Int_t idxMC = 0; idxMC < nPart; ++idxMC) {
  
  return;
}
