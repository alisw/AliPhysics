/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------
//         AliAnalysisTaskSpectraAOD class
//-----------------------------------------------------------------

#include "TChain.h"
#include "TTree.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliVParticle.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskSpectraAOD.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisDataContainer.h"
#include "AliSpectraAODHistoManager.h"
#include "AliSpectraAODTrackCuts.h"
#include "AliSpectraAODEventCuts.h"
#include "AliCentrality.h"
#include "TProof.h"
#include "AliPID.h"
#include "AliVEvent.h"
#include "AliPIDResponse.h"
#include "AliStack.h"
#include <iostream>




using namespace AliSpectraNameSpace;
using namespace std;

ClassImp(AliAnalysisTaskSpectraAOD) // EX1 This stuff tells root to implement the streamer, inspector methods etc (we discussed about it today)

//________________________________________________________________________
AliAnalysisTaskSpectraAOD::AliAnalysisTaskSpectraAOD(const char *name) : AliAnalysisTaskSE(name), fAOD(0), fHistMan(0), fTrackCuts(0), fEventCuts(0), fIsMC(0), fPIDResponse(0), fNSigmaPID(0), fYCut(0)
{
  // Default constructor
  fNSigmaPID = 3.;
  fYCut = .5;
  DefineInput(0, TChain::Class());
  DefineOutput(1, AliSpectraAODHistoManager::Class());
  DefineOutput(2, AliSpectraAODEventCuts::Class());
  DefineOutput(3, AliSpectraAODTrackCuts::Class());

}
//________________________________________________________________________
Bool_t AliAnalysisTaskSpectraAOD::CheckYCut(AODParticleSpecies_t species, AliAODTrack* track) const
{
  // check if the rapidity is within the set range
  // note: masses are hardcoded for now. we could look them up in the pdg database, but that would mean accecing it 100k+ times per run ...
  Double_t y;
  if (species == kSpProton) { y = track->Y(9.38271999999999995e-01); }
  if ( species == kSpKaon ) { y = track->Y(4.93676999999999977e-01); }
  if ( species == kSpPion)  { y = track->Y(1.39570000000000000e-01); }
  if (TMath::Abs(y) > fYCut || y < -998.) return kFALSE;
  return kTRUE;
}
//____________________________________________________________________________
Bool_t AliAnalysisTaskSpectraAOD::CheckYCut(AliAODMCParticle* particle) const
{
  // check if the rapidity is within the set range
  Double_t y = particle->Y();
  if (TMath::Abs(y) > fYCut || y < -998.) return kFALSE;
  return kTRUE;
}
//________________________________________________________________________
void AliAnalysisTaskSpectraAOD::UserCreateOutputObjects()
{
  // create output objects
  fHistMan = new AliSpectraAODHistoManager("SpectraHistos");

  if (!fTrackCuts) AliFatal("Track Cuts should be set in the steering macro");
  if (!fEventCuts) AliFatal("Event Cuts should be set in the steering macro");

  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  PostData(1, fHistMan);
  PostData(2, fEventCuts);
  PostData(3, fTrackCuts);

}
//________________________________________________________________________
void AliAnalysisTaskSpectraAOD::UserExec(Option_t *)
{
  // main event loop
  fAOD = dynamic_cast<AliAODEvent*>(fInputEvent);
  if (strcmp(fAOD->ClassName(), "AliAODEvent"))
    {
      AliFatal("Not processing AODs");
    }
  
  if(!fPIDResponse) {
    AliError("Cannot get pid response");
    return;
  }
  //Selection on QVector, before ANY other selection on the event
  //Spectra MUST be normalized wrt events AFTER the selection on Qvector
  // Can we include this in fEventCuts
  Double_t Qx2EtaPos = 0, Qy2EtaPos = 0;
  Double_t Qx2EtaNeg = 0, Qy2EtaNeg = 0;
  Int_t multPos = 0;
  Int_t multNeg = 0;
  for(Int_t iT = 0; iT < fAOD->GetNumberOfTracks(); iT++) {
    AliAODTrack* aodTrack = fAOD->GetTrack(iT);
    if (!fTrackCuts->IsSelected(aodTrack)) continue;
    if (aodTrack->Eta() >= 0){
      multPos++;
      Qx2EtaPos += TMath::Cos(2*aodTrack->Phi()); 
      Qy2EtaPos += TMath::Sin(2*aodTrack->Phi());
    } else {
      multNeg++;
      Qx2EtaNeg += TMath::Cos(2*aodTrack->Phi()); 
      Qy2EtaNeg += TMath::Sin(2*aodTrack->Phi());
    }
  } 
  Double_t qPos=-999;
  if(multPos!=0)qPos= TMath::Sqrt((Qx2EtaPos*Qx2EtaPos + Qy2EtaPos*Qy2EtaPos)/multPos);
  Double_t qNeg=-999;
  if(multNeg!=0)qNeg= TMath::Sqrt((Qx2EtaNeg*Qx2EtaNeg + Qy2EtaNeg*Qy2EtaNeg)/multNeg);
  
  if((qPos>fTrackCuts->GetQvecMin() && qPos<fTrackCuts->GetQvecMax()) || (qNeg>fTrackCuts->GetQvecMin() && qNeg<fTrackCuts->GetQvecMax())){
    
    //check on centrality distribution
    fHistMan->GetPtHistogram("CentCheck")->Fill(fAOD->GetCentrality()->GetCentralityPercentile("V0M"),fAOD->GetHeader()->GetCentralityP()->GetCentralityPercentileUnchecked("V0M"));
    
    if(!fEventCuts->IsSelected(fAOD))return;//event selection
      
    //fill q distributions vs centrality, after all event selection
    fHistMan->GetqVecHistogram(kHistqVecPos)->Fill(qPos,fAOD->GetCentrality()->GetCentralityPercentile("V0M"));  // qVector distribution
    fHistMan->GetqVecHistogram(kHistqVecNeg)->Fill(qNeg,fAOD->GetCentrality()->GetCentralityPercentile("V0M"));  // qVector distribution
    
    //AliCentrality fAliCentral*;
    //   if ((fAOD->GetCentrality())->GetCentralityPercentile("V0M") > 5.) return;
    
    // First do MC to fill up the MC particle array, such that we can use it later
    TClonesArray *arrayMC = 0;
    if (fIsMC)
      {
	arrayMC = (TClonesArray*) fAOD->GetList()->FindObject(AliAODMCParticle::StdBranchName());
	if (!arrayMC) {
	  AliFatal("Error: MC particles branch not found!\n");
	}
	Int_t nMC = arrayMC->GetEntries();
	for (Int_t iMC = 0; iMC < nMC; iMC++)
	  {
	    AliAODMCParticle *partMC = (AliAODMCParticle*) arrayMC->At(iMC);
	    if(!partMC->Charge()) continue;//Skip neutrals
	    if(TMath::Abs(partMC->Eta()) > fTrackCuts->GetEta()) continue;
	    
	    fHistMan->GetPtHistogram(kHistPtGen)->Fill(partMC->Pt(),partMC->IsPhysicalPrimary());
	    // check for true PID + and fill P_t histos 
	    if (CheckYCut(partMC) ){// only primary vertices and y cut satisfied
	      if ( partMC->PdgCode() == 2212) { fHistMan->GetPtHistogram(kHistPtGenTruePrimaryProtonPlus)->Fill(partMC->Pt(),partMC->IsPhysicalPrimary()); } 
	      if ( partMC->PdgCode() == -2212) { fHistMan->GetPtHistogram(kHistPtGenTruePrimaryProtonMinus)->Fill(partMC->Pt(),partMC->IsPhysicalPrimary()); } 
	      if ( partMC->PdgCode() == 321)  { fHistMan->GetPtHistogram(kHistPtGenTruePrimaryKaonPlus)->Fill(partMC->Pt(),partMC->IsPhysicalPrimary()); } 
	      if ( partMC->PdgCode() == -321) { fHistMan->GetPtHistogram(kHistPtGenTruePrimaryKaonMinus)->Fill(partMC->Pt(),partMC->IsPhysicalPrimary()); } 
	      if ( partMC->PdgCode() == 211)  { fHistMan->GetPtHistogram(kHistPtGenTruePrimaryPionPlus)->Fill(partMC->Pt(),partMC->IsPhysicalPrimary()); } 
	      if ( partMC->PdgCode() == -211) { fHistMan->GetPtHistogram(kHistPtGenTruePrimaryPionMinus)->Fill(partMC->Pt(),partMC->IsPhysicalPrimary()); } 
	    }
	  }
      }
    
    //main loop on tracks
    for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++)
      {
	AliAODTrack* track = fAOD->GetTrack(iTracks);
	if (!fTrackCuts->IsSelected(track)) continue;
	
	//cut on q vectors track-by-track
	if ((qPos>fTrackCuts->GetQvecMin() && qPos<fTrackCuts->GetQvecMax() && track->Eta()<0) || (qNeg>fTrackCuts->GetQvecMin() && qNeg<fTrackCuts->GetQvecMax() && track->Eta()>=0)){
	  
	  //calculate DCA for AOD track
	  Double_t d[2], covd[3];
	  AliAODTrack* track_clone=(AliAODTrack*)track->Clone("track_clone"); // need to clone because PropagateToDCA updates the track parameters
	  Bool_t isDCA = track_clone->PropagateToDCA(fAOD->GetPrimaryVertex(),fAOD->GetMagneticField(),9999.,d,covd);
	  delete track_clone;
	  if(!isDCA)d[0]=-999;
	  
	  fHistMan->GetPtHistogram(kHistPtRec)->Fill(track->Pt(),d[0]);  // PT histo
	  //Response
	  fHistMan->GetPIDHistogram(kHistPIDTPC)->Fill(track->GetTPCmomentum(), track->GetTPCsignal()*track->Charge()); // PID histo
	  
	  AliVParticle *inEvHMain = dynamic_cast<AliVParticle *>(track);
	  Double_t nsigmaTPCkProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kProton));
	  Double_t nsigmaTPCkKaon = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kKaon)); 
	  Double_t nsigmaTPCkPion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kPion)); 
	  Double_t nsigmaTOFkProton=0,nsigmaTOFkKaon=0,nsigmaTOFkPion=0;
	  if(track->Pt()>fTrackCuts->GetPtTOFMatching()){
	    fHistMan->GetPIDHistogram(kHistPIDTOF)->Fill(track->P(),(track->GetTOFsignal()/100)*track->Charge()); // PID histo
	    nsigmaTOFkProton = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kProton));
	    nsigmaTOFkKaon = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kKaon)); 
	    nsigmaTOFkPion = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kPion)); 
	    
	    //TOF
	    fHistMan->GetPtHistogram(kHistNSigProtonTOF)->Fill(track->P(),fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kProton));
	    fHistMan->GetPtHistogram(kHistNSigKaonTOF)->Fill(track->P(),fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kKaon));
	    fHistMan->GetPtHistogram(kHistNSigPionTOF)->Fill(track->P(),fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kPion));
	    fHistMan->GetPtHistogram(kHistNSigProtonPtTOF)->Fill(track->Pt(),fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kProton));
	    fHistMan->GetPtHistogram(kHistNSigKaonPtTOF)->Fill(track->Pt(),fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kKaon));
	    fHistMan->GetPtHistogram(kHistNSigPionPtTOF)->Fill(track->Pt(),fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kPion));
	  }
	  
	  Double_t nsigmaTPCTOFkProton = TMath::Sqrt(nsigmaTPCkProton*nsigmaTPCkProton+nsigmaTOFkProton*nsigmaTOFkProton);
	  Double_t nsigmaTPCTOFkKaon = TMath::Sqrt(nsigmaTPCkKaon*nsigmaTPCkKaon+nsigmaTOFkKaon*nsigmaTOFkKaon);
	  Double_t nsigmaTPCTOFkPion = TMath::Sqrt(nsigmaTPCkPion*nsigmaTPCkPion+nsigmaTOFkPion*nsigmaTOFkPion);
	  
	  //TPC
	  fHistMan->GetPtHistogram(kHistNSigProtonTPC)->Fill(track->P(),fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kProton));
	  fHistMan->GetPtHistogram(kHistNSigKaonTPC)->Fill(track->P(),fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kKaon));
	  fHistMan->GetPtHistogram(kHistNSigPionTPC)->Fill(track->P(),fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kPion));
	  fHistMan->GetPtHistogram(kHistNSigProtonPtTPC)->Fill(track->Pt(),fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kProton));
	  fHistMan->GetPtHistogram(kHistNSigKaonPtTPC)->Fill(track->Pt(),fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kKaon));
	  fHistMan->GetPtHistogram(kHistNSigPionPtTPC)->Fill(track->Pt(),fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kPion));
	  //TPCTOF
	  fHistMan->GetPtHistogram(kHistNSigProtonTPCTOF)->Fill(track->P(),nsigmaTPCTOFkProton);
	  fHistMan->GetPtHistogram(kHistNSigKaonTPCTOF)->Fill(track->P(),nsigmaTPCTOFkKaon);
	  fHistMan->GetPtHistogram(kHistNSigPionTPCTOF)->Fill(track->P(),nsigmaTPCTOFkPion);
	  fHistMan->GetPtHistogram(kHistNSigProtonPtTPCTOF)->Fill(track->Pt(),nsigmaTPCTOFkProton);
	  fHistMan->GetPtHistogram(kHistNSigKaonPtTPCTOF)->Fill(track->Pt(),nsigmaTPCTOFkKaon);
	  fHistMan->GetPtHistogram(kHistNSigPionPtTPCTOF)->Fill(track->Pt(),nsigmaTPCTOFkPion);
	  
	  
	  if( ( nsigmaTPCTOFkKaon < nsigmaTPCTOFkPion ) && ( nsigmaTPCTOFkKaon < nsigmaTPCTOFkProton )) { 
	    if ((nsigmaTPCTOFkKaon > fNSigmaPID) || (!CheckYCut(kSpKaon, track) ) ) continue;
	    if ( track->Charge() > 0 ) { fHistMan->GetPtHistogram(kHistPtRecSigmaKaonPlus)->Fill(track->Pt(),d[0]); } 
	    else { fHistMan->GetPtHistogram(kHistPtRecSigmaKaonMinus)->Fill(track->Pt(),d[0]); } 
	  }
	  if( ( nsigmaTPCTOFkProton < nsigmaTPCTOFkKaon ) && ( nsigmaTPCTOFkProton < nsigmaTPCTOFkPion ) ) {
	    if ( nsigmaTPCTOFkProton > fNSigmaPID || (!CheckYCut(kSpProton, track) ) )  continue;
	    if ( track->Charge() > 0 ) { fHistMan->GetPtHistogram(kHistPtRecSigmaProtonPlus)->Fill(track->Pt(),d[0]); }
	    else { fHistMan->GetPtHistogram(kHistPtRecSigmaProtonMinus)->Fill(track->Pt(),d[0]); }
	  }
	  if( (nsigmaTPCTOFkPion < nsigmaTPCTOFkProton ) && ( nsigmaTPCTOFkPion < nsigmaTPCTOFkKaon ) ) {
	    if (nsigmaTPCTOFkPion > fNSigmaPID || (!CheckYCut(kSpPion, track) ) ) continue;
	    if ( track->Charge() > 0 )  { fHistMan->GetPtHistogram(kHistPtRecSigmaPionPlus)->Fill(track->Pt(),d[0]); }
	    else  { fHistMan->GetPtHistogram(kHistPtRecSigmaPionMinus)->Fill(track->Pt(),d[0]); }
	  }
	  /* MC Part */
	  // MF 22/02/2012
	  // fill mc DCA vs pt
	  if (arrayMC) {
	    AliAODMCParticle *partMC = (AliAODMCParticle*) arrayMC->At(TMath::Abs(track->GetLabel()));
	    if (!partMC) { 
	      AliError("Cannot get MC particle");
	      continue; }
	    if (partMC->IsPhysicalPrimary())fHistMan->GetPtHistogram(kHistPtRecPrimary)->Fill(track->Pt(),d[0]);  // PT histo
	    if (CheckYCut(partMC)){
	      // primaries, true pid
	      //25th Apr - nsigma cut added in addition to the PDG code
	      if ( partMC->PdgCode() == 2212 && nsigmaTPCTOFkProton < fNSigmaPID) { fHistMan->GetPtHistogram(kHistPtRecTrueProtonPlus)->Fill(track->Pt(),d[0]); 
		if (partMC->IsPhysicalPrimary()) {fHistMan->GetPtHistogram(kHistPtRecTruePrimaryProtonPlus)->Fill(track->Pt(),d[0]); }}
	      if ( partMC->PdgCode() == -2212 && nsigmaTPCTOFkProton < fNSigmaPID) { fHistMan->GetPtHistogram(kHistPtRecTrueProtonMinus)->Fill(track->Pt(),d[0]); 
		if (partMC->IsPhysicalPrimary()) {fHistMan->GetPtHistogram(kHistPtRecTruePrimaryProtonMinus)->Fill(track->Pt(),d[0]); }}
	      if ( partMC->PdgCode() == 321 && nsigmaTPCTOFkKaon < fNSigmaPID) { fHistMan->GetPtHistogram(kHistPtRecTrueKaonPlus)->Fill(track->Pt(),d[0]); 
		if (partMC->IsPhysicalPrimary()) {fHistMan->GetPtHistogram(kHistPtRecTruePrimaryKaonPlus)->Fill(track->Pt(),d[0]); }}
	      if ( partMC->PdgCode() == -321 && nsigmaTPCTOFkKaon < fNSigmaPID) { fHistMan->GetPtHistogram(kHistPtRecTrueKaonMinus)->Fill(track->Pt(),d[0]); 
		if (partMC->IsPhysicalPrimary()) {fHistMan->GetPtHistogram(kHistPtRecTruePrimaryKaonMinus)->Fill(track->Pt(),d[0]); }}
	      if ( partMC->PdgCode() == 211 && nsigmaTPCTOFkPion < fNSigmaPID) { fHistMan->GetPtHistogram(kHistPtRecTruePionPlus)->Fill(track->Pt(),d[0]); 
		if (partMC->IsPhysicalPrimary()) {fHistMan->GetPtHistogram(kHistPtRecTruePrimaryPionPlus)->Fill(track->Pt(),d[0]); }}
	      if ( partMC->PdgCode() == -211 && nsigmaTPCTOFkPion < fNSigmaPID) { fHistMan->GetPtHistogram(kHistPtRecTruePionMinus)->Fill(track->Pt(),d[0]);  
		if (partMC->IsPhysicalPrimary()) {fHistMan->GetPtHistogram(kHistPtRecTruePrimaryPionMinus)->Fill(track->Pt(),d[0]); }}
	      //25th Apr - Muons are added to Pions
	      if ( partMC->PdgCode() == 13 && nsigmaTPCTOFkPion < fNSigmaPID) { 
		fHistMan->GetPtHistogram(kHistPtRecTrueMuonPlus)->Fill(track->Pt(),d[0]); 
		if (partMC->IsPhysicalPrimary()) {
		  fHistMan->GetPtHistogram(kHistPtRecTruePrimaryMuonPlus)->Fill(track->Pt(),d[0]); 
		}}
	      if ( partMC->PdgCode() == -13 && nsigmaTPCTOFkPion < fNSigmaPID) { 
		fHistMan->GetPtHistogram(kHistPtRecTrueMuonMinus)->Fill(track->Pt(),d[0]); 
		if (partMC->IsPhysicalPrimary()) {
		  fHistMan->GetPtHistogram(kHistPtRecTruePrimaryMuonMinus)->Fill(track->Pt(),d[0]); 
		}}
	      
	      // primaries, sigma pid 
	      if (partMC->IsPhysicalPrimary()) { 
		if( ( nsigmaTPCTOFkKaon < nsigmaTPCTOFkPion ) && ( nsigmaTPCTOFkKaon < nsigmaTPCTOFkProton ) ) {
		  if ( (nsigmaTPCTOFkKaon > fNSigmaPID ) || (!CheckYCut(kSpKaon, track) ) ) continue; 
		  if ( track->Charge() > 0 ) { fHistMan->GetPtHistogram(kHistPtRecSigmaPrimaryKaonPlus)->Fill(track->Pt(),d[0]); } 
		  else { fHistMan->GetPtHistogram(kHistPtRecSigmaPrimaryKaonMinus)->Fill(track->Pt(),d[0]); } 
		}
		if( ( nsigmaTPCTOFkProton < nsigmaTPCTOFkKaon ) && ( nsigmaTPCTOFkProton < nsigmaTPCTOFkPion ) ) {
		  if ( (nsigmaTPCTOFkProton > fNSigmaPID ) || (!CheckYCut(kSpProton, track) ) ) continue;
		  if ( track->Charge() > 0 ) { fHistMan->GetPtHistogram(kHistPtRecSigmaPrimaryProtonPlus)->Fill(track->Pt(),d[0]); }
		  else { fHistMan->GetPtHistogram(kHistPtRecSigmaPrimaryProtonMinus)->Fill(track->Pt(),d[0]); }
		}
		if( (nsigmaTPCTOFkPion < nsigmaTPCTOFkProton ) && ( nsigmaTPCTOFkPion < nsigmaTPCTOFkKaon ) ) {
		  if ( ( nsigmaTPCTOFkPion > fNSigmaPID )  || (!CheckYCut(kSpPion, track) ) ) continue;
		  if ( track->Charge() > 0 )  { fHistMan->GetPtHistogram(kHistPtRecSigmaPrimaryPionPlus)->Fill(track->Pt(),d[0]); }
		  else  { fHistMan->GetPtHistogram(kHistPtRecSigmaPrimaryPionMinus)->Fill(track->Pt(),d[0]); }
		}
	      }//end if primaries
	      // MF 22/02/2012
	      // Distinguish weak decay and material
	      //done quickly, code can be improved
	      else {//to be added in a separate class
		Int_t mfl=-999,codemoth=-999;
		Int_t indexMoth=partMC->GetMother(); // FIXME ignore fakes? TO BE CHECKED, on ESD is GetFirstMother()
		if(indexMoth>=0){//is not fake
		  AliAODMCParticle* moth = (AliAODMCParticle*) arrayMC->At(indexMoth);
		  codemoth = TMath::Abs(moth->GetPdgCode());
		  mfl = Int_t (codemoth/ TMath::Power(10, Int_t(TMath::Log10(codemoth))));
		}
		//UInt_t flag; 
		//flag=partMC->GetStatus();
		//Printf("FLAG: %d",flag);
		//if(mfl==3 && (flag&AliAODMCParticle::kPDecay)!=0){//strangeness KPDecay not working on AOD, to be understood
		//Printf("\n\n\n STRANGENESS!!!!!");
		//if(codemoth!=-999)Printf("mfl:%d    codemoth%d",mfl,codemoth);
		//if(codemoth==310 || codemoth==130 || codemoth==321 || codemoth==3122 || codemoth==3112 ||
		//codemoth==3222 || codemoth==3312 || codemoth==3322 || codemoth==3334){//K0_S, K0_L, K^+-,lambda, sigma0,sigma+,xi-,xi0, omega
		if(mfl==3){//strangeness
		  if( ( nsigmaTPCkKaon < nsigmaTPCkPion ) && ( nsigmaTPCkKaon < nsigmaTPCkProton ) ) { 
		    if ( (nsigmaTPCkKaon > fNSigmaPID )  || (!CheckYCut(kSpKaon, track) ) ) continue;
		    if ( track->Charge() > 0 ) { fHistMan->GetPtHistogram(kHistPtRecSigmaSecondaryWeakDecayKaonPlus)->Fill(track->Pt(),d[0]); } 
		    else { fHistMan->GetPtHistogram(kHistPtRecSigmaSecondaryWeakDecayKaonMinus)->Fill(track->Pt(),d[0]); } 
		  }
		  if( ( nsigmaTPCkProton < nsigmaTPCkKaon ) && ( nsigmaTPCkProton < nsigmaTPCkPion ) ) {
		    if ( (nsigmaTPCkProton > fNSigmaPID )  || (!CheckYCut(kSpProton, track) ) ) continue;
		    if ( track->Charge() > 0 ) { fHistMan->GetPtHistogram(kHistPtRecSigmaSecondaryWeakDecayProtonPlus)->Fill(track->Pt(),d[0]); }
		    else { fHistMan->GetPtHistogram(kHistPtRecSigmaSecondaryWeakDecayProtonMinus)->Fill(track->Pt(),d[0]); }
		  }
		  if( (nsigmaTPCkPion < nsigmaTPCkProton ) && ( nsigmaTPCkPion < nsigmaTPCkKaon ) ) {
		    if ( ( nsigmaTPCkPion > fNSigmaPID )  || (!CheckYCut(kSpPion, track) ) ) continue;
		    if ( track->Charge() > 0 )  { fHistMan->GetPtHistogram(kHistPtRecSigmaSecondaryWeakDecayPionPlus)->Fill(track->Pt(),d[0]); }
		    else  { fHistMan->GetPtHistogram(kHistPtRecSigmaSecondaryWeakDecayPionMinus)->Fill(track->Pt(),d[0]); }
		  }
		}//end if strangeness
		else{//material
		  if( ( nsigmaTPCkKaon < nsigmaTPCkPion ) && ( nsigmaTPCkKaon < nsigmaTPCkProton ) ) { 
		    if ( (nsigmaTPCkKaon > fNSigmaPID )  || (!CheckYCut(kSpKaon, track) ) ) continue;
		    if ( track->Charge() > 0 ) { fHistMan->GetPtHistogram(kHistPtRecSigmaSecondaryMaterialKaonPlus)->Fill(track->Pt(),d[0]); } 
		    else { fHistMan->GetPtHistogram(kHistPtRecSigmaSecondaryMaterialKaonMinus)->Fill(track->Pt(),d[0]); } 
		  }
		  if( ( nsigmaTPCkProton < nsigmaTPCkKaon ) && ( nsigmaTPCkProton < nsigmaTPCkPion ) ) {
		    if ( (nsigmaTPCkProton > fNSigmaPID )  || (!CheckYCut(kSpProton, track) ) ) continue;
		    if ( track->Charge() > 0 ) { fHistMan->GetPtHistogram(kHistPtRecSigmaSecondaryMaterialProtonPlus)->Fill(track->Pt(),d[0]); }
		    else { fHistMan->GetPtHistogram(kHistPtRecSigmaSecondaryMaterialProtonMinus)->Fill(track->Pt(),d[0]); }
		  }
		  if( (nsigmaTPCkPion < nsigmaTPCkProton ) && ( nsigmaTPCkPion < nsigmaTPCkKaon ) ) {
		    if ( ( nsigmaTPCkPion > fNSigmaPID )  || (!CheckYCut(kSpPion, track) ) ) continue;
		    if ( track->Charge() > 0 )  { fHistMan->GetPtHistogram(kHistPtRecSigmaSecondaryMaterialPionPlus)->Fill(track->Pt(),d[0]); }
		    else  { fHistMan->GetPtHistogram(kHistPtRecSigmaSecondaryMaterialPionMinus)->Fill(track->Pt(),d[0]); }
		  }
		}//end if material
	      }//end else (IsPhysPrim)
	    }//end if(CheckY)
	  }//end if(arrayMC)
	}//end if on q vector (track)
      } // end loop on tracks
  }//end if on q vector (event)
  PostData(1, fHistMan);
  PostData(2, fEventCuts);
  PostData(3, fTrackCuts);
}

//_________________________________________________________________
void   AliAnalysisTaskSpectraAOD::Terminate(Option_t *)
{
  // Terminate
}
