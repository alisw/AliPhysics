//
// Hadron Trigger & Jet Substructure.
//
// Summer Student: C.Durante - Supervisor: L.Cunqueiro

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TTree.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TChain.h>
#include <TSystem.h>
#include <TFile.h>
#include <TKey.h>
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVector3.h"
#include "TVector2.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliEmcalParticle.h"
#include "AliMCEvent.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliPythiaInfo.h"
#include "TRandom3.h"
#include "AliPicoTrack.h"
#include "AliEmcalJetFinder.h"
#include "AliAODEvent.h"
#include "AliAnalysisTaskEmcalMissingEnergy.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskEmcalMissingEnergy)

//________________________________________________________________________
AliAnalysisTaskEmcalMissingEnergy::AliAnalysisTaskEmcalMissingEnergy() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalMissingEnergy", kTRUE),
  fContainer(0),
  fMinFractionShared(0),
  fJetShapeType(kData),
  fJetShapeSub(kNoSub),
  fJetSelection(kInclusive),
  fSubstructureVar(0),
  fPtThreshold(-9999.),
  fRMatching(0.2),
  fJetRadius(0.4),
  fSubJetRadius(0.2),
  fminpTTrig(20.),
  fmaxpTTrig(50.),
  fangWindowRecoil(0.6),
  fSemigoodCorrect(0),
  fHolePos(0),
  fHoleWidth(0), 
  fCentSelectOn(kTRUE),
  fCentMin(0),
  fCentMax(10),
  fh2ResponseUW(0x0),
  fh2ResponseW(0x0), 
  fPhiJetCorr6(0x0), 
  fPhiJetCorr7(0x0),
  fEtaJetCorr6(0x0),
  fEtaJetCorr7(0x0),
  fPtJetCorr(0x0),
  fPtJet(0x0),
  fhpTjetpT(0x0),
  fhPt(0x0),
  fhJetPt(0x0),
  fhTrackPt(0x0),
  fhTriggerPt(0x0),
  fhJetTriggeredPt(0x0),
  fhTriggerVSjetPt(0x0),
  fhdPhiTrigVSjetPt(0x0),
  fhdPhiTrigVSjetPtNOCUT(0x0),
  fhdPhiPartVSjetPt(0x0),
  fhdPhiPartVSjetPtVSPartPt(0x0),
  fhPhi(0x0),
  fSubstructure(0),
  fHadronTrigger(0),
  fhTau1(0x0),
  fhTau2(0x0),
  fhTau3(0x0),
  fhTau1vsTau2(0x0),
  fhTau1vsTau3(0x0),
  fhTau2vsTau3(0x0),
  fhTau2OverTau1(0x0),
  fhTau3OverTau2(0x0),
  fhTau2vsJetPt(0x0),
  fhTau2OverTau1vsJetPt(0x0),
  fhTau1vsTau2vsTau3(0x0),
  fhTrackPt_JT(0x0)

{
  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalMissingEnergy::AliAnalysisTaskEmcalMissingEnergy(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fContainer(0),
  fMinFractionShared(0),
  fJetShapeType(kData),
  fJetShapeSub(kNoSub),
  fJetSelection(kInclusive),
  fSubstructureVar(0),
  fPtThreshold(-9999.),
  fRMatching(0.2),
  fJetRadius(0.4),
  fSubJetRadius(0.2),
  fminpTTrig(20.),
  fmaxpTTrig(50.),
  fangWindowRecoil(0.6),
  fSemigoodCorrect(0),
  fHolePos(0),
  fHoleWidth(0), 
  fCentSelectOn(kTRUE),
  fCentMin(0),
  fCentMax(10),
  fh2ResponseUW(0x0),
  fh2ResponseW(0x0),
  fPhiJetCorr6(0x0), 
  fPhiJetCorr7(0x0),
  fEtaJetCorr6(0x0),
  fEtaJetCorr7(0x0),
  fPtJetCorr(0x0),
  fPtJet(0x0),
  fhpTjetpT(0x0),
  fhPt(0x0),
  fhJetPt(0x0),
  fhTrackPt(0x0),
  fhTriggerPt(0x0),
  fhJetTriggeredPt(0x0),
  fhTriggerVSjetPt(0x0),
  fhdPhiTrigVSjetPt(0x0),
  fhdPhiTrigVSjetPtNOCUT(0x0),
  fhdPhiPartVSjetPt(0x0),
  fhdPhiPartVSjetPtVSPartPt(0x0),
  fhPhi(0x0),
  fSubstructure(0),
  fHadronTrigger(0),
  fhTau1(0x0),
  fhTau2(0x0),
  fhTau3(0x0),
  fhTau1vsTau2(0x0),
  fhTau1vsTau3(0x0),
  fhTau2vsTau3(0x0),
  fhTau2OverTau1(0x0),
  fhTau3OverTau2(0x0),
  fhTau2vsJetPt(0x0),
  fhTau2OverTau1vsJetPt(0x0),
  fhTau1vsTau2vsTau3(0x0),
  fhTrackPt_JT(0x0)
  
{
  // Standard constructor.
  
  SetMakeGeneralHistograms(kTRUE);
  
  DefineOutput(1, TTree::Class());
  
}

//________________________________________________________________________
AliAnalysisTaskEmcalMissingEnergy::~AliAnalysisTaskEmcalMissingEnergy()
{
  // Destructor.
}

//________________________________________________________________________
 void AliAnalysisTaskEmcalMissingEnergy::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  // trees

  fSubstructure = new TTree("fSubstructure","fSubstructure");
  fHadronTrigger = new TTree("fHadronTrigger","fHadronTrigger");

  const Int_t nVarS = 5;  
  const Int_t nVarH = 5; 
  
  fSubstructureVar = new Float_t [nVarS];
  fHadronTriggerVar = new Float_t [nVarH];

  TString *fSubstructureVarNames = new TString [nVarS];
  TString *fHadronTriggerVarNames = new TString [nVarH];
  
  fSubstructureVarNames[0] = "JetPt";
  fSubstructureVarNames[1] = "Tau1";
  fSubstructureVarNames[2] = "Tau2";
  fSubstructureVarNames[3] = "Tau3";
  fSubstructureVarNames[4] = "DeltaR2hardest";

  fHadronTriggerVarNames[0] = "TriggerPt";
  fHadronTriggerVarNames[1] = "JetPt";
  fHadronTriggerVarNames[2] = "DeltaPhiJetTrigger";
  fHadronTriggerVarNames[3] = "PartPt";
  fHadronTriggerVarNames[4] = "DeltaPhiJetParticles";

  for(Int_t ivar = 0; ivar<nVarS; ivar++) {
    cout<<"looping over variables to create the branches of my Substructure Tree"<<endl;
    fSubstructure->Branch(fSubstructureVarNames[ivar].Data(), &fSubstructureVar[ivar], Form("%s/F", fSubstructureVarNames[ivar].Data()));
  }

  for(Int_t ivar = 0; ivar<nVarH; ivar++) {
    cout<<"looping over variables to create the branches of my HadronTrigger Tree"<<endl;
    fHadronTrigger->Branch(fHadronTriggerVarNames[ivar].Data(), &fHadronTriggerVar[ivar], Form("%s/F", fHadronTriggerVarNames[ivar].Data()));
  }

  // hadron trigger output

  fhJetPt= new TH1F("fhJetPt", "fhJetPt;Jet p_{T} (GeV)", 100, 0, 200);
  fOutput->Add(fhJetPt);

  fhTrackPt= new TH1F("fhTrackPt", "fhTrackPt;Track p_{T} (GeV)", 100, 0, 200);
  fOutput->Add(fhTrackPt);

  fhTrackPt_JT= new TH1F("fhTrackPt_JT", "fhTrackPt_JT;Track p_{T} (GeV)", 100, 0, 200);
  fOutput->Add(fhTrackPt_JT);

  fhTriggerPt= new TH1F("fhTriggerPt", "fhTriggerPt;Trigger p_{T} (GeV)", 100, 0, 200);
  fOutput->Add(fhTriggerPt);

  fhJetTriggeredPt= new TH1F("fhJetTriggeredPt", "fhJetTriggeredPt;Triggered Jet p_{T} (GeV)", 100, 0, 200);
  fOutput->Add(fhJetTriggeredPt);

  fhTriggerVSjetPt= new TH2F("fhTriggerVSjetPt", "fhTriggerVSjetPt;Trigger p_{T} (GeV);Jet p_{T} (GeV)", 100, 0, 200,  100, 0, 200);
  fOutput->Add(fhTriggerVSjetPt);

  fhdPhiTrigVSjetPtNOCUT= new TH2F("fhdPhiTrigVSjetPtNOCUT", "fhdPhiTrigVSjetPtNOCUT;#Delta#phi (Trig-Jet);Jet p_{T} (GeV)", 100, -0.5*TMath::Pi(), 1.5*TMath::Pi(),  100, 0, 200);
  fOutput->Add(fhdPhiTrigVSjetPtNOCUT);

  fhdPhiTrigVSjetPt= new TH2F("fhdPhiTrigVSjetPt", "fhdPhiTrigVSjetPt;#Delta#phi (Trig-Jet);Jet p_{T} (GeV)", 100, TMath::Pi() - 0.6, TMath::Pi(),  100, 0, 200);
  fOutput->Add(fhdPhiTrigVSjetPt);


  fhdPhiPartVSjetPt= new TH2F("fhdPhiPartVSjetPt", "fhdPhiPartVSjetPt;#Delta#phi (Part-Jet);Jet p_{T} (GeV)", 100, -0.5*TMath::Pi(), 1.5*TMath::Pi(),  100, 0, 200);
  fOutput->Add(fhdPhiPartVSjetPt);

  fhdPhiPartVSjetPtVSPartPt= new TH3F("fhdPhiPartVSjetPtVSPartPt", "fhdPhiPartVSjetPtVSPartPt;#Delta#phi (Part-Jet);Jet p_{T} (GeV);Part p_{T} (GeV)", 100, -0.5*TMath::Pi(), 1.5*TMath::Pi(),  100, 0, 200, 100, 0, 200);
  /*fhdPhiPartVSjetPtVSPartPt->GetXaxis()->SetTitle("Delta phi");
  fhdPhiPartVSjetPtVSPartPt->GetYaxis()->SetTitle("Jet Pt");
  fhdPhiPartVSjetPtVSPartPt->GetZaxis()->SetTitle("Particle Pt");*/
  fOutput->Add(fhdPhiPartVSjetPtVSPartPt);

  // substructure output

  fhTau1 = new TH1F("fhTau1", "fhTau1;#tau_{1}", 100, 0, 1);
  fOutput->Add(fhTau1);

  fhTau2 = new TH1F("fhTau2", "fhTau2;#tau_{2}", 100, 0, 1);
  fOutput->Add(fhTau2);

  fhTau3 = new TH1F("fhTau3", "fhTau3;#tau_{3}", 100, 0, 1);
  fOutput->Add(fhTau3);



  fhTau1vsTau2 = new TH2F("fhTau1vsTau2", "fhTau1vsTau2;#tau_{1};#tau_{2}", 100, 0, 1,  100, 0, 1);
  fOutput->Add(fhTau1vsTau2);
  
  fhTau1vsTau3 = new TH2F("fhTau1vsTau3", "fhTau1vsTau3;#tau_{1};#tau_{3}", 100, 0, 1,  100, 0, 1);
  fOutput->Add(fhTau1vsTau3);

  fhTau2vsTau3 = new TH2F("fhTau2vsTau3", "fhTau2vsTau3;#tau_{2};#tau_{3}", 100, 0, 1,  100, 0, 1);
  fOutput->Add(fhTau2vsTau3);



  fhTau1vsTau2vsTau3 = new TH3F("fhTau1vsTau2vsTau3", "fhTau1vsTau2vsTau3;#tau_{1};#tau_{2};#tau_{3}", 100, 0, 1,  100, 0, 1, 100, 0, 1);
  /* fhTau1vsTau2vsTau3->GetXaxis()->SetTitle("\tau _1");
  fhTau1vsTau2vsTau3->GetYaxis()->SetTitle("\tau _2");
  fhTau1vsTau2vsTau3->GetZaxis()->SetTitle("\tau _3");*/
  fOutput->Add(fhTau1vsTau2vsTau3);



  fhTau2OverTau1 = new TH1F("fhTau2OverTau1", "fhTau2OverTau1;#tau_{2} / #tau_{1}", 100, 0, 1);
  fOutput->Add(fhTau2OverTau1);

  fhTau3OverTau2 = new TH1F("fhTau3OverTau2", "fhTau3OverTau2;#tau_{3} / #tau_{2}", 100, 0, 1);
  fOutput->Add(fhTau3OverTau2);
  

  fhTau2vsJetPt = new TH2F("fhTau2vsJetPt", "fhTau2vsJetPt;#tau_{2}; Jet p_{T} (GeV)", 100, 0, 1,  100, 0, 200);
  fOutput->Add(fhTau2vsJetPt);

  fhTau2OverTau1vsJetPt = new TH2F("fhTau2OverTau1vsJetPt", "fhTau2OverTau1vsJetPt;#tau_{2} / #tau_{1}; Jet p_{T} (GeV)", 100, 0, 1,  100, 0, 200);
  fOutput->Add(fhTau2OverTau1vsJetPt);

  fOutput->Add(fSubstructure);
  fOutput->Add(fHadronTrigger);
  PostData(1, fOutput);  //  Post data for ALL output slots > 0 here.

}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalMissingEnergy::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalMissingEnergy::FillHistograms()
{
  /*  // Fill histograms.
  //cout<<"base container"<<endl;
  AliEmcalJet* jet1 = NULL;
  AliJetContainer *jetCont = GetJetContainer(0);
  Float_t kWeight=1;
  if (fCentSelectOn)
    if ((fCent>fCentMax) || (fCent<fCentMin)) return 0;
  
  AliAODTrack *triggerHadron = 0x0;
  */
  /////////////////////////////////////////// start hadron trigger
  
  const Double_t eta_max = 0.9;
  const Double_t track_pt_min = 0.15;
  const Double_t jet_pt_min = 10;
  const Double_t phi_h = 0.6;
  
  // const Double_t epsilon = 0.00001;

  AliParticleContainer *partCont = GetParticleContainer(0);
  TClonesArray *tracksArray = partCont->GetArray();
  
  if(!partCont || !tracksArray) return 0;
  //  AliAODTrack *track = 0x0;
  //  AliEmcalParticle *emcPart = 0x0;
  AliPicoTrack *trigger = 0x0;
  AliPicoTrack *track = 0x0;

  AliEmcalJet* jetTEMP = NULL;
  AliJetContainer *JetCont = GetJetContainer(0);
  
  /*      OLD WAY
  for(Int_t iTrack=0; iTrack <= tracksArray->GetEntriesFast(); iTrack++){
    
      picoTrack = static_cast<AliPicoTrack*>(tracksArray->At(iTrack));
      if (!picoTrack) continue;
      
      
      if(TMath::Abs(picoTrack->Eta())>eta_max) continue;
      if(picoTrack->Pt()<track_pt_min) continue;
      if(picoTrack->GetTrackType() == 2) continue;
  
      fhTrackPt->Fill(picoTrack->Pt());
      
      if(picoTrack->Pt() < fminpTTrig || picoTrack->Pt() > fmaxpTTrig) continue; // selecting trigger hadron ... !!!
      fhTriggerPt->Fill(picoTrack->Pt());
      
      if(JetCont) {
	JetCont->ResetCurrentID();
	while((jetTEMP = JetCont->GetNextAcceptJet())) {
	  if (!jetTEMP) continue;
	  if (jetTEMP->Pt() < jet_pt_min) continue;
	  if (TMath::Abs(jetTEMP->Eta()) > eta_max - fJetRadius) continue; // TO BE UPDATED
	 
	  fhJetPt->Fill(jetTEMP->Pt());

	  Float_t dphi = RelativePhi (picoTrack->Phi(), jetTEMP->Phi()); // -pi to pi
	            
	  fhdPhiTrigVSjetPtNOCUT->Fill(TMath::Abs(dphi),jetTEMP->Pt());

	  if (TMath::Abs(dphi) < TMath::Pi() - phi_h) continue; // trigger selection
	  
	  fhJetTriggeredPt->Fill(jetTEMP->Pt());
	  fhTriggerVSjetPt->Fill(picoTrack->Pt(),jetTEMP->Pt());
	  fhdPhiTrigVSjetPt->Fill(TMath::Abs(dphi),jetTEMP->Pt());
	}
      }
  }
  */

  for(Int_t iTrack=0; iTrack <= tracksArray->GetEntriesFast(); iTrack++){ 
    track = static_cast<AliPicoTrack*>(tracksArray->At(iTrack));
    if (!track || track->Pt()<1) continue;
    fhTrackPt->Fill(track->Pt());
  }
    
  if(JetCont) {
    JetCont->ResetCurrentID();
    while((jetTEMP = JetCont->GetNextAcceptJet())) {
      if (!jetTEMP) continue;
      if (jetTEMP->Pt() < jet_pt_min) continue;
      if (TMath::Abs(jetTEMP->Eta()) > eta_max - fJetRadius) continue;
      fhJetPt->Fill(jetTEMP->Pt());
    }
  }

  Int_t trigger_index = SelectTrigger(fminpTTrig,fmaxpTTrig);
      
  if(trigger_index >= 0) {
    trigger = static_cast<AliPicoTrack*>(tracksArray->At(trigger_index));
    if (trigger) {
      fhTriggerPt->Fill(trigger->Pt());
      
      if(JetCont) {
	JetCont->ResetCurrentID();
	while((jetTEMP = JetCont->GetNextAcceptJet())) {
	  if (!jetTEMP) continue;
	  if (jetTEMP->Pt() < jet_pt_min) continue;
	  if (TMath::Abs(jetTEMP->Eta()) > eta_max - fJetRadius) continue;
	  
	  Float_t dphi = RelativePhi(trigger->Phi(), jetTEMP->Phi()); // -pi to pi
    Float_t dphiFancy = RelativePhiFancy(trigger->Phi(), jetTEMP->Phi()); // -0.5pi to 1.5pi
	  
	  fhdPhiTrigVSjetPtNOCUT->Fill(dphiFancy,jetTEMP->Pt());
	  
	  if (TMath::Abs(dphi) < TMath::Pi() - phi_h) continue; // trigger selection
	  
	  fhJetTriggeredPt->Fill(jetTEMP->Pt());
	  fhTriggerVSjetPt->Fill(trigger->Pt(),jetTEMP->Pt());
	  fhdPhiTrigVSjetPt->Fill(TMath::Abs(dphi),jetTEMP->Pt());
	  
	  //now I have a trigger and a recoil jet. Let's loop over all particles!
	  for(Int_t iTrack=0; iTrack <= tracksArray->GetEntriesFast(); iTrack++){ 
	    track = static_cast<AliPicoTrack*>(tracksArray->At(iTrack));
	    if (!track) continue;

	    if(TMath::Abs(track->Eta())>0.9) continue;
	    if(track->Pt()<0.15) continue;
	    //if(track->GetTrackType() == 2) continue;

	    Float_t DPHI = RelativePhiFancy(jetTEMP->Phi(), track->Phi()); // -0.5 pi to 1.5 pi
	    fhdPhiPartVSjetPt->Fill(DPHI,jetTEMP->Pt());
	    fhdPhiPartVSjetPtVSPartPt->Fill(DPHI,jetTEMP->Pt(),track->Pt());
	    fhTrackPt_JT->Fill(track->Pt());

	    fHadronTriggerVar[0] = trigger->Pt();
	    fHadronTriggerVar[1] = jetTEMP->Pt();
	    fHadronTriggerVar[2] = TMath::Abs(dphi);
	    fHadronTriggerVar[3] = track->Pt();
	    fHadronTriggerVar[4] = DPHI;
	    
	    fHadronTrigger->Fill();
	  }
	}
      }
    }
  }
  ////////////////////////////////////////// end hadron trigger




  ////////////////////////////////////////// start jet substructure
  
  AliEmcalJet *mainJet = NULL; 
  AliEmcalJet *subJet = NULL;
  AliEmcalJet *subJet1hardest = NULL;
  AliEmcalJet *subJet2hardest = NULL;
  AliEmcalJet *subJet3hardest = NULL;
  
  //AliVParticle *mainJetParticle = 0x0; //Individual constituent tracks (particles) in a mainJet 
  //AliVParticle *subJetParticle = 0x0; //Individual constituent tracks (particles) in a subJet   
  
  AliPicoTrack *mainJetParticle = 0x0;
  AliPicoTrack *subJetPartiicle = 0x0;

  AliEmcalJetFinder *Reclusterer = new AliEmcalJetFinder("SubJetFinder"); //JetFinder Object for reclustered jets from MC particle Level data
                 
  const AliVVertex *vert = InputEvent()->GetPrimaryVertex();
  Double_t dVtx[3]={vert->GetX(),vert->GetY(),vert->GetZ()};

  Reclusterer->SetRadius(fSubJetRadius);
  Reclusterer->SetJetMinPt(track_pt_min);
  Reclusterer->SetJetAlgorithm(0); //  0  for anti-kt     1 for kt

  Int_t ReclustererOk = 0;
  Int_t SubJetCounter = 0;
  
  Double_t tau1_num = 0.;
  Double_t tau2_num = 0.;
  Double_t tau3_num = 0.;
  Double_t tau_den = 0.;
  Double_t dR1 = 0.;
  Double_t dR2 = 0.;
  Double_t dR3 = 0.;
  Double_t dRmin = 0.;
  Double_t dR2hardest = 0.;
  
  if(JetCont){
    JetCont->ResetCurrentID();
    
    while((mainJet = JetCont->GetNextAcceptJet())) {      // loop over main Jets found in the event
      if(!mainJet) continue;
      if (mainJet->Pt() < jet_pt_min) continue;
      if (TMath::Abs(mainJet->Eta()) > eta_max - fJetRadius) continue; // TO BE UPDATED
      
      if(Reclusterer->AliEmcalJetFinder::Filter(mainJet, JetCont, dVtx)){  //reclustering mainJets using the jetfinderobject Reclusterer                
	
	ReclustererOk = 1;
	SubJetCounter = Reclusterer->GetNumberOfJets(); // Number of reclustered SubJets in each original jet                                            
	for (Int_t i=0; i<SubJetCounter; i++){ //Loops through each SubJet in a reclustered jet 
	  subJet = Reclusterer->GetJet(i);
	  cout<<"SubJet: "<<i+1<<"/"<<SubJetCounter<<" | subJetPt/JetPt: "<<subJet->Pt()<<"/"<<mainJet->Pt()<<endl;
	}
	
	tau_den = 0.;
	for (Int_t i = 0; i < mainJet->GetNumberOfTracks(); i++) { //compute tau denominator once and for all
	  // mainJetParticle =  static_cast<AliVParticle*>(mainJet->TrackAt(i, JetCont->GetParticleContainer()->GetArray()));
	  mainJetParticle =  static_cast<AliPicoTrack*>(mainJet->TrackAt(i, JetCont->GetParticleContainer()->GetArray()));
	  tau_den += mainJetParticle->Pt() * fJetRadius;
	}
	
	if (tau_den <= 0. || 0 == SubJetCounter) continue;
	
	fSubstructureVar[0] = mainJet->Pt();
	
	if (1 == SubJetCounter) { // I have just 1 subjet, it's gonna be the hardest
	  subJet1hardest = Reclusterer->GetJet(0);
	  tau1_num = 0.;
	  for (Int_t i = 0; i < mainJet->GetNumberOfTracks(); i++) { // tau1 computation
	    //mainJetParticle =  static_cast<AliVParticle*>(mainJet->TrackAt(i, JetCont->GetParticleContainer()->GetArray()));
	    mainJetParticle =  static_cast<AliPicoTrack*>(mainJet->TrackAt(i, JetCont->GetParticleContainer()->GetArray()));
	    dR1 = R_distance(subJet1hardest,mainJetParticle);
	    //  dEta1 = mainJetParticle->Eta() - subJet1hardest->Eta();
	    // dPhi1 = RelativePhi(mainJetParticle->Phi(), subJet1hardest->Phi());
	    //  dR1 = TMath::Sqrt(dEta1*dEta1 + dPhi1*dPhi1);
	    tau1_num += mainJetParticle->Pt() * dR1;
	  }
	  
	  if (tau1_num == 0. ) continue;
	  fhTau1->Fill(tau1_num/tau_den);
	  fSubstructureVar[1] = tau1_num/tau_den;
	  fSubstructureVar[2] = -1;
	  fSubstructureVar[3] = -1;
	  fSubstructureVar[4] = -1;
	}
	
	else if (2 == SubJetCounter) { // I have just 2 subjets, the second is the hardest
	  subJet1hardest = Reclusterer->GetJet(1);
	  subJet2hardest = Reclusterer->GetJet(0);
	  tau1_num = 0.;
	  tau2_num = 0.;
	  dR2hardest = R_distance(subJet1hardest,subJet2hardest);
	  for (Int_t i = 0; i < mainJet->GetNumberOfTracks(); i++) { //tau1 & tau2 computation
	    //mainJetParticle =  static_cast<AliVParticle*>(mainJet->TrackAt(i, JetCont->GetParticleContainer()->GetArray()));
	    mainJetParticle =  static_cast<AliPicoTrack*>(mainJet->TrackAt(i, JetCont->GetParticleContainer()->GetArray()));
	    dR1 = R_distance(subJet1hardest,mainJetParticle);
	    dR2 = R_distance(subJet2hardest,mainJetParticle);
	    // dEta1 = mainJetParticle->Eta() - subJet1hardest->Eta();
	    // dPhi1 = RelativePhi(mainJetParticle->Phi(), subJet1hardest->Phi());
	    // dR1 = TMath::Sqrt(dEta1*dEta1 + dPhi1*dPhi1);
	    // dEta2 = mainJetParticle->Eta() - subJet2hardest->Eta();
	    // dPhi2 = (mainJetParticle->Phi(), subJet2hardest->Phi());
	    // dR2 = TMath::Sqrt(dEta2*dEta2 + dPhi2*dPhi2);
	       dRmin = Minimum(dR1, dR2);
	       tau1_num += mainJetParticle->Pt() * dR1;
	       tau2_num += mainJetParticle->Pt() * dRmin;
	  }
	  
	  if (tau1_num == 0. || tau2_num == 0.) continue;
	  fhTau1->Fill(tau1_num/tau_den);
	  fhTau2->Fill(tau2_num/tau_den);
	  fhTau1vsTau2->Fill(tau1_num/tau_den,tau2_num/tau_den);
	  fhTau2OverTau1->Fill(tau2_num/tau1_num);
	  fhTau2vsJetPt->Fill(tau2_num/tau_den,mainJet->Pt());
	  fhTau2OverTau1vsJetPt->Fill(tau2_num/tau1_num,mainJet->Pt());
	  
	  fSubstructureVar[1] = tau1_num/tau_den;
	  fSubstructureVar[2] = tau2_num/tau_den;
	  fSubstructureVar[3] = -1;
	  fSubstructureVar[4] = dR2hardest;
	}
	
	else if (SubJetCounter > 2) { // I have more than 2 subjets, the last is the hardest
	  subJet1hardest = Reclusterer->GetJet(SubJetCounter - 1);
	  subJet2hardest = Reclusterer->GetJet(SubJetCounter - 2);
	  subJet3hardest = Reclusterer->GetJet(SubJetCounter - 3);
	  tau1_num = 0.;
	  tau2_num = 0.;
	  tau3_num = 0.;
	  dR2hardest = R_distance(subJet1hardest,subJet2hardest);
	  for (Int_t i = 0; i < mainJet->GetNumberOfTracks(); i++) { //tau1 & tau2 computation
	    //mainJetParticle =  static_cast<AliVParticle*>(mainJet->TrackAt(i, JetCont->GetParticleContainer()->GetArray()));
	    mainJetParticle =  static_cast<AliPicoTrack*>(mainJet->TrackAt(i, JetCont->GetParticleContainer()->GetArray()));
	    dR1 = R_distance(subJet1hardest,mainJetParticle);
	    dR2 = R_distance(subJet2hardest,mainJetParticle);
	    dR3 = R_distance(subJet3hardest,mainJetParticle);
	    //   dEta1 = mainJetParticle->Eta() - subJet1hardest->Eta();
	    // dPhi1 = RelativePhi(mainJetParticle->Phi(), subJet1hardest->Phi());
	    //dR1 = TMath::Sqrt(dEta1*dEta1 + dPhi1*dPhi1);
	    //dEta2 = mainJetParticle->Eta() - subJet2hardest->Eta();
	    //dPhi2 = RelativePhi(mainJetParticle->Phi(), subJet2hardest->Phi());
	    // dR2 = TMath::Sqrt(dEta2*dEta2 + dPhi2*dPhi2);
	    // dEta3 = mainJetParticle->Eta() - subJet3hardest->Eta();
	    // dPhi3 = RelativePhi(mainJetParticle->Phi(), subJet3hardest->Phi());
	    //dR3 = TMath::Sqrt(dEta3*dEta3 + dPhi3*dPhi3);
	    
	    tau1_num += mainJetParticle->Pt() * dR1;
	    
	    dRmin = Minimum(dR1,dR2);
	    tau2_num += mainJetParticle->Pt() * dRmin;
	    
	    dRmin = Minimum(dR1,dR2,dR3);
	    tau3_num += mainJetParticle->Pt() * dRmin;
	  }
	  
	  if (tau1_num == 0. || tau2_num == 0. || tau3_num == 0.) continue;
	  fhTau1->Fill(tau1_num/tau_den);
	  fhTau2->Fill(tau2_num/tau_den);
	  fhTau3->Fill(tau3_num/tau_den);
	  fhTau1vsTau2->Fill(tau1_num/tau_den,tau2_num/tau_den);
	  fhTau1vsTau3->Fill(tau1_num/tau_den,tau3_num/tau_den);
	  fhTau2vsTau3->Fill(tau2_num/tau_den,tau3_num/tau_den);
	  fhTau1vsTau2vsTau3->Fill(tau1_num/tau_den,tau2_num/tau_den,tau3_num/tau_den); 
	  fhTau2OverTau1->Fill(tau2_num/tau1_num);
	  fhTau3OverTau2->Fill(tau3_num/tau2_num);
	  fhTau2vsJetPt->Fill(tau2_num/tau_den,mainJet->Pt());
	  fhTau2OverTau1vsJetPt->Fill(tau2_num/tau1_num,mainJet->Pt());
	  
	  fSubstructureVar[1] = tau1_num/tau_den;
	  fSubstructureVar[2] = tau2_num/tau_den;
	  fSubstructureVar[3] = tau3_num/tau_den;
	  fSubstructureVar[4] = dR2hardest;
	  
	}
	
	fSubstructure->Fill();
	
      }
    }
  }
  
  

   ////////////////////////////////////////// end jet substructure


   // OLD STUFF


  /*
  if (fJetSelection == kRecoil) {
    //Printf("Recoil jets!!!, fminpTTrig = %f, fmaxpTTrig = %f", fminpTTrig, fmaxpTTrig);
    Int_t triggerHadronLabel = SelectTrigger(fminpTTrig, fmaxpTTrig);
     
    
    if (triggerHadronLabel==-99999) {
      //Printf ("Trigger Hadron not found, return");
      return 0;}

   
    AliParticleContainer *partContAn = GetParticleContainer(0);
    TClonesArray *trackArrayAn = partContAn->GetArray();
    triggerHadron = static_cast<AliAODTrack*>(trackArrayAn->At(triggerHadronLabel));
  
    if (!triggerHadron) {
      //Printf("No Trigger hadron with the found label!!");
      return 0;
    }

    if(fSemigoodCorrect){
      Double_t disthole=RelativePhi(triggerHadron->Phi(),fHolePos);
      if(TMath::Abs(disthole)+fHoleWidth>TMath::Pi()-fangWindowRecoil){
        return 0;}
    }
   
    fhPt->Fill(triggerHadron->Pt());

  }
  
  if(jetCont) {
    jetCont->ResetCurrentID();
    while((jet1 = jetCont->GetNextAcceptJet())) {
      if (!jet1) continue;
      AliEmcalJet* jet2 = 0x0;
      AliEmcalJet* jet3=0x0;
      fPtJet->Fill(jet1->Pt());
      AliEmcalJet *jetUS = NULL;
      Int_t ifound=0;
      Int_t ilab=-1;
      
      if(fSemigoodCorrect && (fJetSelection != kRecoil)){
      Double_t disthole=RelativePhi(jet1->Phi(),fHolePos);
      if(TMath::Abs(disthole)<fHoleWidth){
      continue;}
    } 
 
      if (!(fJetShapeType == kData)) {
        AliPythiaInfo *partonsInfo = 0x0;
        if((fJetShapeType == kTrueDet) || (fJetShapeType == kDetEmb) || (fJetShapeType == kPythiaDef) || (fJetShapeType == kDetEmbPart) ){
          AliJetContainer *jetContTrue = GetJetContainer(1);
          AliJetContainer *jetContUS = GetJetContainer(2);
       
 
          if(fJetShapeSub==kConstSub){
            for(Int_t i = 0; i<jetContUS->GetNJets(); i++) {
              jetUS = jetContUS->GetJet(i);
              if(jetUS->GetLabel()==jet1->GetLabel()) {
                ifound++;
                if(ifound==1) ilab = i;
              }
            }
            if(ilab==-1) continue;
            jetUS=jetContUS->GetJet(ilab);
            
            jet2=jetUS->ClosestJet();
          }
          if(!(fJetShapeSub==kConstSub)) jet2 = jet1->ClosestJet();
          if (!jet2) {
            Printf("jet2 does not exist, returning");
            continue;
          }
          
          if(fJetShapeType==kDetEmbPart){
            AliJetContainer *jetContPart=GetJetContainer(3);
            jet3=jet2->ClosestJet();
           
	    if(!jet3){
	      Printf("jet3 does not exist, returning");
              continue;}
	    cout<<"jet 3 exists"<<jet3->Pt()<<endl;
}
            
          
          fh2ResponseUW->Fill(jet1->Pt(),jet2->Pt());
          
          Double_t fraction=0;
          if(!(fJetShapeSub==kConstSub))  fraction = jetCont->GetFractionSharedPt(jet1);
          if(fJetShapeSub==kConstSub) fraction = jetContUS->GetFractionSharedPt(jetUS);
          //if (fraction > 0.1) cout<<"***** hey a jet matched with fraction"<<fraction<<"  "<<jet1->Pt()<<" "<<jet2->Pt()<<" "<<fCent<<endl;
          
          if(fraction<fMinFractionShared) continue;
          //InputEvent()->Print();
          if (!(fJetShapeType == kPythiaDef)) {
            partonsInfo = (AliPythiaInfo*) jetContTrue->GetPythiaInfo();
            if(!partonsInfo) return 0;
          }
        }
        else {
          partonsInfo = (AliPythiaInfo*) jetCont->GetPythiaInfo();
          jet2=jet1;
          if(!partonsInfo) return 0;
        }
        
        if (!(fJetShapeType == kPythiaDef)){
          Double_t jp1=RelativePhi(jet2->Phi(),partonsInfo->GetPartonPhi6());
          Double_t detap1=(jet2->Eta())-(partonsInfo->GetPartonEta6());
          kWeight=partonsInfo->GetPythiaEventWeight();
          fh2ResponseW->Fill(jet1->Pt(),jet2->Pt(),kWeight);
          
          Float_t dRp1 = TMath::Sqrt(jp1 * jp1 + detap1 * detap1);
          fEtaJetCorr6->Fill(jet2->Eta(), partonsInfo->GetPartonEta6());
          fPhiJetCorr6->Fill(jet2->Phi(), partonsInfo->GetPartonPhi6());
          if(dRp1 < fRMatching) {
            fShapesVar[0] = partonsInfo->GetPartonFlag6();
            fPtJetCorr ->Fill(partonsInfo->GetPartonPt6(), jet2->Pt());
          }
          else {
            jp1=RelativePhi(jet2->Phi(),partonsInfo->GetPartonPhi7());
            detap1=(jet2->Eta())-(partonsInfo->GetPartonEta7());
            dRp1 = TMath::Sqrt(jp1 * jp1 + detap1 * detap1);
            fEtaJetCorr7->Fill(jet2->Eta(), partonsInfo->GetPartonEta7());
            fPhiJetCorr7->Fill(jet2->Phi(), partonsInfo->GetPartonPhi7());
            if(dRp1 < fRMatching) {
              fShapesVar[0] = partonsInfo->GetPartonFlag7();
              fPtJetCorr->Fill(partonsInfo->GetPartonPt7(), jet2->Pt());
            }
            else fShapesVar[0]=0;
          }
        }
        else
          fShapesVar[0] = -1.;
      }
      else
        fShapesVar[0] = 0.;
    
      Double_t ptSubtracted = 0;
      
     
      
      Float_t dphiRecoil = 0.;
      if (fJetSelection == kRecoil){
        dphiRecoil = RelativePhi(triggerHadron->Phi(), jet1->Phi());
        if (TMath::Abs(dphiRecoil) < (TMath::Pi() - fangWindowRecoil)) {
         // Printf("Recoil jets back to back not found! continuing");
          continue;
        }
        
        fhpTjetpT->Fill(triggerHadron->Pt(), jet1->Pt());
        //Printf(" ************ FILLING HISTOS****** shapeSub = %d, triggerHadron = %f, jet1 = %f", fJetShapeSub, triggerHadron->Pt(), jet1->Pt());
        fhPhi->Fill(RelativePhi(triggerHadron->Phi(), jet1->Phi()));
        
      }
      
      if ((((fJetShapeType == kData) || (fJetShapeType == kDetEmb)) && (fJetShapeSub == kConstSub))|| (fJetShapeType==kPythiaDef))
        ptSubtracted = jet1->Pt();
      
      else ptSubtracted  = jet1->Pt() - GetRhoVal(0)*jet1->Area();
      
      if ((fJetShapeType == kData) || (fJetShapeType== kDetEmb)||(fJetShapeType==kTrueDet) || (fJetShapeType==kPythiaDef))
        if (ptSubtracted < fPtThreshold) continue;
     
     
      if ((fCentSelectOn==kFALSE) && (jet1->GetNumberOfTracks() <= 1)) continue;
      
      fShapesVar[1] = ptSubtracted;
      fShapesVar[2] = GetJetpTD(jet1,0);
      //fShapesVar[3] = GetJetMass(jet1,0);
      // fShapesVar[4] = 1.*GetJetNumberOfConstituents(jet1,0);
      fShapesVar[3] = GetJetAngularity(jet1,0);
      fShapesVar[4] = GetJetCircularity(jet1,0);
      fShapesVar[5] = GetJetLeSub(jet1,0);
      fShapesVar[6] = GetSigma2(jet1,0);

      Float_t ptMatch=0., ptDMatch=0., massMatch=0., constMatch=0.,angulMatch=0.,circMatch=0., lesubMatch=0., sigma2Match=0.;
      Int_t kMatched = 0;
      if (fJetShapeType == kTrueDet || fJetShapeType == kDetEmb || fJetShapeType==kPythiaDef) {
        kMatched = 1;
        ptMatch=jet2->Pt();
        ptDMatch=GetJetpTD(jet2, kMatched);
        //massMatch=GetJetMass(jet2,kMatched);
	// constMatch=1.*GetJetNumberOfConstituents(jet2,kMatched);
        angulMatch=GetJetAngularity(jet2, kMatched);
        circMatch=GetJetCircularity(jet2, kMatched);
        lesubMatch=GetJetLeSub(jet2, kMatched);
        sigma2Match = GetSigma2(jet2, kMatched);
        
      }
      
       if (fJetShapeType == kDetEmbPart) {
	 if(fJetShapeSub==kConstSub) kMatched = 3;
         if(fJetShapeSub==kDerivSub) kMatched = 2;
        ptMatch=jet3->Pt();
        ptDMatch=GetJetpTD(jet3, kMatched);
	// massMatch=GetJetMass(jet3,kMatched);
	// constMatch=1.*GetJetNumberOfConstituents(jet3,kMatched);
        angulMatch=GetJetAngularity(jet3, kMatched);
        circMatch=GetJetCircularity(jet3, kMatched);
        lesubMatch=GetJetLeSub(jet3, kMatched);
        sigma2Match = GetSigma2(jet3, kMatched);
        
      }



      if (fJetShapeType == kTrue || fJetShapeType == kData) {
        kMatched = 0;
        ptMatch=0.;
        ptDMatch=0.;
        //massMatch=0.;
	// constMatch=0.;
        angulMatch=0.;
        circMatch=0.;
        lesubMatch=0.;
        sigma2Match =0.;
        
      }
      
      fShapesVar[7] = ptMatch;
      fShapesVar[8] = ptDMatch;
      // fShapesVar[11] = massMatch;
      // fShapesVar[12] = constMatch;
      fShapesVar[9] = angulMatch;
      fShapesVar[10] = circMatch;
      fShapesVar[11] = lesubMatch;
      fShapesVar[12] = sigma2Match;
      fShapesVar[13] = kWeight;
      fShapesVar[14] =  0.;  // updated
      fShapesVar[15] =  0.;  // updated
      fTreeObservableTagging->Fill();
      
    }
    
  }
  */
  return kTRUE;
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalMissingEnergy::GetJetMass(AliEmcalJet *jet,Int_t jetContNb=0) {
  //calc subtracted jet mass
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    return jet->GetSecondOrderSubtracted();
  else 
    return jet->M();
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalMissingEnergy::Angularity(AliEmcalJet *jet, Int_t jetContNb = 0){

  AliJetContainer *jetCont = GetJetContainer(jetContNb);
  if (!jet->GetNumberOfTracks())
    return 0; 
    Double_t den=0.;
    Double_t num = 0.;
    AliVParticle *vp1 = 0x0;
    for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
      vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
      
      if (!vp1){
        Printf("AliVParticle associated to constituent not found");
        continue;
      }
      
      Double_t dphi = RelativePhi(vp1->Phi(),jet->Phi());
      Double_t dr2 = (vp1->Eta()-jet->Eta())*(vp1->Eta()-jet->Eta()) + dphi*dphi;
      Double_t dr = TMath::Sqrt(dr2);
      num=num+vp1->Pt()*dr;
      den=den+vp1->Pt();
    }
    return num/den;
} 

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalMissingEnergy::GetJetAngularity(AliEmcalJet *jet, Int_t jetContNb = 0) {

  if((fJetShapeSub==kDerivSub) && (jetContNb==0))
    return jet->GetSecondOrderSubtractedAngularity();
  else
    return Angularity(jet, jetContNb);
 
}


//________________________________________________________________________
Float_t AliAnalysisTaskEmcalMissingEnergy::PTD(AliEmcalJet *jet, Int_t jetContNb = 0){

  AliJetContainer *jetCont = GetJetContainer(jetContNb);
  if (!jet->GetNumberOfTracks())
      return 0; 
    Double_t den = 0.;
    Double_t num = 0.;
    AliVParticle *vp1 = 0x0;
    for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
      vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
      
      if (!vp1){
        Printf("AliVParticle associated to constituent not found");
        continue;
      }
      
      num=num+vp1->Pt()*vp1->Pt();
      den=den+vp1->Pt();
    }
    return TMath::Sqrt(num)/den;
} 

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalMissingEnergy::GetJetpTD(AliEmcalJet *jet, Int_t jetContNb = 0) {
  //calc subtracted jet mass
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    return jet->GetSecondOrderSubtractedpTD();
  else
    return PTD(jet, jetContNb);
 
}

//_____________________________________________________________________________
Float_t AliAnalysisTaskEmcalMissingEnergy::Circularity(AliEmcalJet *jet, Int_t jetContNb = 0){

  AliJetContainer *jetCont = GetJetContainer(jetContNb);
  if (!jet->GetNumberOfTracks())
    return 0; 
  Double_t mxx = 0.;
  Double_t myy = 0.;
  Double_t mxy = 0.;
  int nc = 0;
  Double_t sump2 = 0.;
  Double_t pxjet=jet->Px();
  Double_t pyjet=jet->Py();
  Double_t pzjet=jet->Pz();
  
  
  //2 general normalized vectors perpendicular to the jet
  TVector3  ppJ1(pxjet, pyjet, pzjet);
  TVector3  ppJ3(- pxjet* pzjet, - pyjet * pzjet, pxjet * pxjet + pyjet * pyjet);
  ppJ3.SetMag(1.);
  TVector3  ppJ2(-pyjet, pxjet, 0);
  ppJ2.SetMag(1.);
  AliVParticle *vp1 = 0x0;
  for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
    vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));  
    
    if (!vp1){
      Printf("AliVParticle associated to constituent not found");
      continue;
    }
    
    TVector3 pp(vp1->Px(), vp1->Py(), vp1->Pz());
   
    //local frame
    TVector3 pLong = pp.Dot(ppJ1) / ppJ1.Mag2() * ppJ1;
    TVector3 pPerp = pp - pLong;
    //projection onto the two perpendicular vectors defined above
    
    Float_t ppjX = pPerp.Dot(ppJ2);
    Float_t ppjY = pPerp.Dot(ppJ3);
    Float_t ppjT = TMath::Sqrt(ppjX * ppjX + ppjY * ppjY);
    if(ppjT<=0) return 0;
    
    mxx += (ppjX * ppjX / ppjT);
    myy += (ppjY * ppjY / ppjT);
    mxy += (ppjX * ppjY / ppjT);
    nc++;
    sump2 += ppjT;}
  
  if(nc<2) return 0;
  if(sump2==0) return 0;
  // Sphericity Matrix
  Double_t ele[4] = {mxx / sump2, mxy / sump2, mxy / sump2, myy / sump2};
  TMatrixDSym m0(2,ele);
  
  // Find eigenvectors
  TMatrixDSymEigen m(m0);
  TVectorD eval(2);
  TMatrixD evecm = m.GetEigenVectors();
  eval  = m.GetEigenValues();
  // Largest eigenvector
  int jev = 0;
  //  cout<<eval[0]<<" "<<eval[1]<<endl;
  if (eval[0] < eval[1]) jev = 1;
  TVectorD evec0(2);
  // Principle axis
  evec0 = TMatrixDColumn(evecm, jev);
  Double_t compx=evec0[0];
  Double_t compy=evec0[1];
  TVector2 evec(compx, compy);
  Double_t circ=0;
  if(jev==1) circ=2*eval[0];
  if(jev==0) circ=2*eval[1];
  
  return circ;
}




//________________________________________________________________________
Float_t AliAnalysisTaskEmcalMissingEnergy::GetJetCircularity(AliEmcalJet *jet, Int_t jetContNb =0 ) {
  //calc subtracted jet mass
 
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    return jet->GetSecondOrderSubtractedCircularity();
  else
    return Circularity(jet, jetContNb);
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalMissingEnergy::LeSub(AliEmcalJet *jet, Int_t jetContNb =0 ){

  AliJetContainer *jetCont = GetJetContainer(jetContNb);
  if (!jet->GetNumberOfTracks())
    return 0;
  Double_t den=0.;
  Double_t num = 0.;
  AliVParticle *vp1 = 0x0;
  AliVParticle *vp2 = 0x0;
  std::vector<int> ordindex;
  ordindex=jet->SortConstituentsPt(jetCont->GetParticleContainer()->GetArray());
  //Printf("Nbof const = %d", jet->GetNumberOfTracks());
  //Printf("ordindex[0] = %d, ordindex[1] = %d", ordindex[0], ordindex[1]);
  
  if(ordindex.size()<2) return -1;
  
  vp1 = static_cast<AliVParticle*>(jet->TrackAt(ordindex[0], jetCont->GetParticleContainer()->GetArray()));
  if (!vp1){
    Printf("AliVParticle associated to Leading constituent not found");
    return -1;
  }
  
  vp2 = static_cast<AliVParticle*>(jet->TrackAt(ordindex[1], jetCont->GetParticleContainer()->GetArray()));
  if (!vp2){
    Printf("AliVParticle associated to Subleading constituent not found");
    return -1;
  }
  
  
  num=vp1->Pt();
  den=vp2->Pt();
  //Printf("vp1->Pt() =%f, vp2->Pt() =%f", vp1->Pt(), vp2->Pt());
  
return num-den;
} 

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalMissingEnergy::GetJetLeSub(AliEmcalJet *jet, Int_t jetContNb =0) {
  //calc subtracted jet mass
 
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    return jet->GetSecondOrderSubtractedLeSub();
  else
    return LeSub(jet, jetContNb);
 
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalMissingEnergy::GetJetNumberOfConstituents(AliEmcalJet *jet,Int_t jetContNb=0) {
  //calc subtracted jet mass
  
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    return jet->GetSecondOrderSubtractedConstituent();
  else
    return jet->GetNumberOfTracks();
 
}
   

//______________________________________________________________________________
Float_t AliAnalysisTaskEmcalMissingEnergy::Sigma2(AliEmcalJet *jet, Int_t jetContNb=0){

  AliJetContainer *jetCont = GetJetContainer(jetContNb);
  if (!jet->GetNumberOfTracks())
      return 0; 
      Double_t mxx    = 0.;
      Double_t myy    = 0.;
      Double_t mxy    = 0.;
      int  nc     = 0;
      Double_t sump2  = 0.;
       
     AliVParticle *vp1 = 0x0;
     for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
       vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
       
       if (!vp1){
         Printf("AliVParticle associated to constituent not found");
         continue;
       }
       
       Double_t ppt=vp1->Pt();
       Double_t dphi = RelativePhi(vp1->Phi(),jet->Phi());
     
       Double_t deta = vp1->Eta()-jet->Eta();
       mxx += ppt*ppt*deta*deta;
       myy += ppt*ppt*dphi*dphi;
       mxy -= ppt*ppt*deta*TMath::Abs(dphi);
       nc++;
       sump2 += ppt*ppt;
       
     }  
     if(nc<2) return 0;
     if(sump2==0) return 0;
     // Sphericity Matrix
     Double_t ele[4] = {mxx , mxy , mxy , myy };
     TMatrixDSym m0(2,ele);
     
     // Find eigenvectors
     TMatrixDSymEigen m(m0);
     TVectorD eval(2);
     TMatrixD evecm = m.GetEigenVectors();
     eval  = m.GetEigenValues();
     // Largest eigenvector
     int jev = 0;
     //  cout<<eval[0]<<" "<<eval[1]<<endl;
     if (eval[0] < eval[1]) jev = 1;
     TVectorD evec0(2);
     // Principle axis
     evec0 = TMatrixDColumn(evecm, jev);
     Double_t compx=evec0[0];
     Double_t compy=evec0[1];
     TVector2 evec(compx, compy);
     Double_t sig=0;
     if(jev==1) sig=TMath::Sqrt(TMath::Abs(eval[0])/sump2);
     if(jev==0) sig=TMath::Sqrt(TMath::Abs(eval[1])/sump2);
     
     return sig;
     
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalMissingEnergy::GetSigma2(AliEmcalJet *jet, Int_t jetContNb=0) {
  //calc subtracted jet mass
 
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    return jet->GetSecondOrderSubtractedSigma2();
  else
    return Sigma2(jet, jetContNb);
 
}

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalMissingEnergy::SelectTrigger(Float_t minpT, Float_t maxpT){

  AliParticleContainer *partCont = GetParticleContainer(0);
  TClonesArray *tracksArray = partCont->GetArray();
  
  if(!partCont || !tracksArray) return -99999;
  AliAODTrack *track = 0x0;
  AliEmcalParticle *emcPart = 0x0;
  AliPicoTrack *picoTrack = 0x0;
  
  TList *trackList = new TList();
  Int_t triggers[100];
  for (Int_t iTrigger=0; iTrigger<100; iTrigger++) triggers[iTrigger] = 0;
  Int_t iTT = 0;
  
  for(Int_t iTrack=0; iTrack <= tracksArray->GetEntriesFast(); iTrack++){
    
    if ((fJetShapeSub == kNoSub) || (fJetShapeSub == kDerivSub)) {
      picoTrack = static_cast<AliPicoTrack*>(tracksArray->At(iTrack));
      if (!picoTrack) continue;
      
      
      if(TMath::Abs(picoTrack->Eta())>0.9) continue;
      if(picoTrack->Pt()<0.15) continue;
      if(picoTrack->GetTrackType() == 2) continue;
    
      //if ((picoTrack->Pt()>8) && (picoTrack->Pt()<9)) Printf("picoTrackLabel = %d", picoTrack->GetTrackType());
      
      if ((picoTrack->Pt() >= minpT) && (picoTrack->Pt()< maxpT)) {
        trackList->Add(picoTrack);
        triggers[iTT] = iTrack;
        iTT++;
      }
    }
    else if (fJetShapeSub == kConstSub){
      emcPart = static_cast<AliEmcalParticle*>(tracksArray->At(iTrack));
      if (!emcPart) continue;
      if(TMath::Abs(emcPart->Eta())>0.9) continue;
      if (emcPart->Pt()<0.15) continue;
      
      if ((emcPart->Pt() >= minpT) && (emcPart->Pt()< maxpT)) {
        trackList->Add(emcPart);
        triggers[iTT] = iTrack;
        iTT++;
      }
    }
    else{
      track = static_cast<AliAODTrack*>(tracksArray->At(iTrack));
      if (!track) continue;
      if(TMath::Abs(track->Eta())>0.9) continue;
      if (track->Pt()<0.15) continue;
      if (!(track->TestFilterBit(768))) continue;
      
      if ((track->Pt() >= minpT) && (track->Pt()< maxpT)) {
        trackList->Add(track);
        triggers[iTT] = iTrack;
        iTT++;
        
      }
    }
  }

  if (iTT == 0) return -99999;
  Int_t nbRn = 0, index = 0 ; 
  TRandom3* random = new TRandom3(0); 
  nbRn = random->Integer(iTT);

  index = triggers[nbRn];
  //Printf("iTT Total= %d, nbRn = %d, Index = %d",iTT, nbRn, index );
  return index; 
  
}

//__________________________________________________________________________________
Double_t AliAnalysisTaskEmcalMissingEnergy::RelativePhi(Double_t mphi,Double_t vphi){

  if (vphi < -1*TMath::Pi()) vphi += (2*TMath::Pi());
  else if (vphi > TMath::Pi()) vphi -= (2*TMath::Pi());
  if (mphi < -1*TMath::Pi()) mphi += (2*TMath::Pi());
  else if (mphi > TMath::Pi()) mphi -= (2*TMath::Pi());
  double dphi = mphi-vphi;
  if (dphi < -1*TMath::Pi()) dphi += (2*TMath::Pi());
  else if (dphi > TMath::Pi()) dphi -= (2*TMath::Pi());
  return dphi;//dphi in [-Pi, Pi]
}


//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalMissingEnergy::RetrieveEventObjects() {
  //
  // retrieve event objects
  //
  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  return kTRUE;
}

//_______________________________________________________________________
void AliAnalysisTaskEmcalMissingEnergy::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.

  // fTreeObservableTagging = dynamic_cast<TTree*>(GetOutputData(1));
  // if (!fTreeObservableTagging){
  //   Printf("ERROR: fTreeObservable not available"); 
  //   return;
  // }

}

//_______________________________________________________________________
Double_t AliAnalysisTaskEmcalMissingEnergy::Minimum(Double_t x, Double_t y) 
{
  if (x <= y) return x;
  else return y;
}

//_______________________________________________________________________
Double_t AliAnalysisTaskEmcalMissingEnergy::Minimum(Double_t x, Double_t y, Double_t z) 
{
  Double_t minimum = x;
  if (y < minimum) minimum = y;
  if (z < minimum) minimum  = z;
  return minimum;
}

//_______________________________________________________________________
Double_t AliAnalysisTaskEmcalMissingEnergy::R_distance(AliEmcalJet *jet1, AliEmcalJet *jet2) 
{
  Double_t dPhi, dEta;
  dPhi = RelativePhi(jet1->Phi(),jet2->Phi());
  dEta = jet1->Eta() - jet2->Eta();
  return TMath::Sqrt(dEta*dEta + dPhi*dPhi);
}

//_______________________________________________________________________
Double_t AliAnalysisTaskEmcalMissingEnergy::R_distance(Double_t phi1, Double_t eta1, Double_t phi2, Double_t eta2) 
{
  Double_t dPhi, dEta;
  dPhi = RelativePhi(phi1,phi2);
  dEta = eta1 - eta2;
  return TMath::Sqrt(dEta*dEta + dPhi*dPhi);
}

//_______________________________________________________________________
Double_t AliAnalysisTaskEmcalMissingEnergy::R_distance(AliEmcalJet *jet,  AliPicoTrack *part) 
{
  Double_t dPhi, dEta;
  dPhi = RelativePhi(jet->Phi(),part->Phi());
  dEta = jet->Eta() - part->Eta();
  return TMath::Sqrt(dEta*dEta + dPhi*dPhi);
}

//__________________________________________________________________________________
Double_t AliAnalysisTaskEmcalMissingEnergy::RelativePhiFancy(Double_t mphi,Double_t vphi){

  Double_t dPhi = mphi - vphi;
  
  if(dPhi>2*TMath::Pi()) dPhi -= 2*TMath::Pi();
  if(dPhi<-2*TMath::Pi()) dPhi += 2*TMath::Pi();
  if(dPhi<-0.5*TMath::Pi()) dPhi += 2*TMath::Pi();
  if(dPhi>1.5*TMath::Pi()) dPhi -= 2*TMath::Pi();
  
  return dPhi; // [-0.5 pi, 1.5 pi]
}
