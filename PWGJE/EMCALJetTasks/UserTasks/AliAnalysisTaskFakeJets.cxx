//
// Trying to see if fake jets can be identified by their shape, collection of selected shapes
//
// Author: L. Cunqueiro
// n-subjetiness code copied from Nima Zhardosti's task AliAnalysisTaskSubetFraction

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
#include "TRandom3.h"
#include "AliEmcalJetFinder.h"
#include "AliAODEvent.h"
#include "AliAnalysisTaskFakeJets.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskFakeJets)

//________________________________________________________________________
AliAnalysisTaskFakeJets::AliAnalysisTaskFakeJets() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskFakeJets", kTRUE),
  fContainer(0),
  fMinFractionShared(0),
  fJetShapeType(kData),
  fJetShapeSub(kNoSub),
  fJetSelection(kInclusive),
  fShapesVar(0),
  fPtThreshold(-9999.),
  fRMatching(0.2),
  fminpTTrig(20.),
  fmaxpTTrig(50.),
  fangWindowRecoil(0.6),
  fSemigoodCorrect(0),
  fHolePos(0),
  fHoleWidth(0), 
  fCentSelectOn(kTRUE),
  fCentMin(0),
  fCentMax(10),
  fOneConstSelectOn(kFALSE),
  fDerivSubtrOrder(0),
    fSubjetRadius(0.2),
  fJetRadius(0.4),
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
  fhPhi(0x0),
  fNbOfConstvspT(0x0),
  fTreeFakeJets(0)

{
  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskFakeJets::AliAnalysisTaskFakeJets(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fContainer(0),
  fMinFractionShared(0),
  fJetShapeType(kData),
  fJetShapeSub(kNoSub),
  fJetSelection(kInclusive),
  fShapesVar(0),
  fPtThreshold(-9999.),
  fRMatching(0.2),
  fminpTTrig(20.),
  fmaxpTTrig(50.),
  fangWindowRecoil(0.6),
  fSemigoodCorrect(0),
  fHolePos(0),
  fHoleWidth(0), 
  fCentSelectOn(kTRUE),
  fCentMin(0),
  fCentMax(10),
  fOneConstSelectOn(kFALSE),
   fDerivSubtrOrder(0),
  fSubjetRadius(0.2),
  fJetRadius(0.4),
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
  fhPhi(0x0),
  fNbOfConstvspT(0x0),
  fTreeFakeJets(0)
  
{
  // Standard constructor.
  
  SetMakeGeneralHistograms(kTRUE);

  DefineOutput(1, TTree::Class());

}

//________________________________________________________________________
AliAnalysisTaskFakeJets::~AliAnalysisTaskFakeJets()
{
  // Destructor.
}

//________________________________________________________________________
 void AliAnalysisTaskFakeJets::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fTreeFakeJets = new TTree("fTreeFakeJets", "fTreeFakeJets");
  const Int_t nVar = 11;
  fShapesVar = new Float_t [nVar]; 
  TString *fShapesVarNames = new TString [nVar];


  fShapesVarNames[0] = "ptJet"; 
  fShapesVarNames[1] = "ptDJet"; 
  fShapesVarNames[2] = "angularity";
   fShapesVarNames[3] = "angularitysquared";
  fShapesVarNames[4]= "hardtrack";
  fShapesVarNames[5]="hard2track";
  fShapesVarNames[6]="corefrac";
  fShapesVarNames[7]="nsubjet1";
  fShapesVarNames[8]="nsubjet2";
  fShapesVarNames[9]="subjetfrac";
 fShapesVarNames[10]="mass";
  for(Int_t ivar=0; ivar < nVar; ivar++){
    cout<<"looping over variables"<<endl;
    fTreeFakeJets->Branch(fShapesVarNames[ivar].Data(), &fShapesVar[ivar], Form("%s/F", fShapesVarNames[ivar].Data()));

    //if( ivar == 4 )  fTreeObservableTagging->Branch(fShapesVarNames[ivar].Data(), &fShapesVar[ivar], Form("%s/I", fShapesVarNames[ivar].Data()));

  }

  fh2ResponseUW= new TH2F("fh2ResponseUW", "fh2ResponseUW", 100, 0, 200,  100, 0, 200); 
  fOutput->Add(fh2ResponseUW);
  fh2ResponseW= new TH2F("fh2ResponseW", "fh2ResponseW", 100, 0, 200,  100, 0, 200);
  fOutput->Add(fh2ResponseW);
  fPhiJetCorr6= new TH2F("fPhiJetCorr6", "fPhiJetCorr6", 50, 0, 2*TMath::Pi(), 50, 0, 2*TMath::Pi());
  fOutput->Add(fPhiJetCorr6);
  fEtaJetCorr6= new TH2F("fEtaJetCorr6", "fEtaJetCorr6", 50, -1.5, 1.5, 50, -1.5, 1.5);
  fOutput->Add(fEtaJetCorr6);
  
  fPhiJetCorr7= new TH2F("fPhiJetCorr7", "fPhiJetCorr7", 50, 0, 2*TMath::Pi(), 50, 0, 2*TMath::Pi());
  fOutput->Add(fPhiJetCorr7);
  fEtaJetCorr7= new TH2F("fEtaJetCorr7", "fEtaJetCorr7", 50, -1.5, 1.5, 50, -1.5, 1.5);
  fOutput->Add(fEtaJetCorr7);
  
  fPtJetCorr= new TH2F("fPtJetCorr", "fPtJetCorr", 100, 0, 200,  100, 0, 200);
  fOutput->Add(fPtJetCorr);
  fPtJet= new TH1F("fPtJet", "fPtJet", 100, 0, 200);
  fOutput->Add(fPtJet);
  
  fhpTjetpT= new TH2F("fhpTjetpT", "fhpTjetpT", 200, 0, 200,  200, 0, 200);
  fOutput->Add(fhpTjetpT);
  fhPt= new TH1F("fhPt", "fhPt", 200, 0, 200);
  fOutput->Add(fhPt);
  fhPhi= new TH1F("fhPhi", "fhPhi", 100, -TMath::Pi(), TMath::Pi());
  fOutput->Add(fhPhi);
  
  fNbOfConstvspT=new TH2F("fNbOfConstvspT", "fNbOfConstvspT", 100, 0, 100, 200, 0, 200);
  fOutput->Add(fNbOfConstvspT);
  
  fOutput->Add(fTreeFakeJets);
  TH1::AddDirectory(oldStatus);
  PostData(1, fOutput); // Post data for ALL output slots > 0 here.

}

//________________________________________________________________________
Bool_t AliAnalysisTaskFakeJets::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskFakeJets::FillHistograms()
{
  // Fill histograms.
  //cout<<"base container"<<endl;
  AliEmcalJet* jet1 = NULL;
  AliJetContainer *jetCont = GetJetContainer(0);
  Float_t kWeight=1;
  if (fCentSelectOn)
    if ((fCent>fCentMax) || (fCent<fCentMin)) return 0;
  
  AliAODTrack *triggerHadron = 0x0;
  
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
      
      fPtJet->Fill(jet1->Pt());
      AliEmcalJet *jetUS = NULL;
      Int_t ifound=0;
      Int_t ilab=-1;
      
      if(fSemigoodCorrect && (fJetSelection != kRecoil)){
      Double_t disthole=RelativePhi(jet1->Phi(),fHolePos);
      if(TMath::Abs(disthole)<fHoleWidth){
      continue;}
    }
 
     
       
      
       Double_t ptSubtracted = jet1->Pt();
       if (ptSubtracted < fPtThreshold) continue;
     
     if (fOneConstSelectOn == kTRUE) fNbOfConstvspT->Fill(GetJetNumberOfConstituents(jet1,0), ptSubtracted);
      
      if ((fCentSelectOn == kFALSE) && (jet1->GetNumberOfTracks() <= 1)) continue;

        AliEmcalJetFinder *Reclusterer1; //Object containg Subjets from Subtracted Hybrid Jets
        Reclusterer1 = Recluster(jet1, 0, fSubjetRadius, 0, 0, "SubJetFinder_1");

    
      fShapesVar[0] = ptSubtracted;
      fShapesVar[1] = GetJetpTD(jet1,0);
      fShapesVar[2] = GetJetAngularity(jet1,0);
      fShapesVar[3] = GetJetAngularitySquared(jet1,0);
      fShapesVar[4]= GetJetHardTrack(jet1,0);
      fShapesVar[5]=GetJetSecHardTrack(jet1,0);
      fShapesVar[6]=GetJetCoreFrac(jet1,0);
      fShapesVar[7]=NSubJettiness(jet1, 0, fJetRadius, Reclusterer1, 1, 0, 1);
      fShapesVar[8]=NSubJettiness(jet1, 0, fJetRadius, Reclusterer1, 2, 0, 1);
      fShapesVar[9]=GetSubjetFraction(jet1,0,fJetRadius,Reclusterer1);
     
      fShapesVar[10]=GetJetMass(jet1,0);
      fTreeFakeJets->Fill();
 
       }}
    
  
  
  return kTRUE;
  }

//________________________________________________________________________
Float_t AliAnalysisTaskFakeJets::GetJetMass(AliEmcalJet *jet,Int_t jetContNb=0) {
  //calc subtracted jet mass
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtracted();
      else return jet->GetShapeProperties()->GetSecondOrderSubtracted();
  else 
    return jet->M();
}

//________________________________________________________________________
Float_t AliAnalysisTaskFakeJets::Angularity(AliEmcalJet *jet, Int_t jetContNb = 0){

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
Float_t AliAnalysisTaskFakeJets::GetJetAngularity(AliEmcalJet *jet, Int_t jetContNb = 0) {

  if((fJetShapeSub==kDerivSub) && (jetContNb==0))
    if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedAngularity();
      else return jet->GetShapeProperties()->GetSecondOrderSubtractedAngularity();
  else
    return Angularity(jet, jetContNb);
 
}
//________________________________________________________________________
Float_t AliAnalysisTaskFakeJets::AngularitySquared(AliEmcalJet *jet, Int_t jetContNb = 0){

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
      num=num+vp1->Pt()*TMath::Sqrt(dr);
      den=den+vp1->Pt();
    }
    return num/den;
} 

//________________________________________________________________________
Float_t AliAnalysisTaskFakeJets::GetJetAngularitySquared(AliEmcalJet *jet, Int_t jetContNb = 0) {

  if((fJetShapeSub==kDerivSub) && (jetContNb==0))
    if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedAngularity();
      else return jet->GetShapeProperties()->GetSecondOrderSubtractedAngularity();
  else
    return AngularitySquared(jet, jetContNb);
 
}


//________________________________________________________________________
Float_t AliAnalysisTaskFakeJets::PTD(AliEmcalJet *jet, Int_t jetContNb = 0){

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
      
      num=num+vp1->Pt()*vp1->Pt();
      den=den+vp1->Pt();
    }
    return TMath::Sqrt(num)/den;
} 

//________________________________________________________________________
Float_t AliAnalysisTaskFakeJets::GetJetpTD(AliEmcalJet *jet, Int_t jetContNb = 0) {
  //calc subtracted jet mass
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedpTD();
      else return jet->GetShapeProperties()->GetSecondOrderSubtractedpTD();
  else
    return PTD(jet, jetContNb);
 
}




Float_t AliAnalysisTaskFakeJets::HardTrack(AliEmcalJet *jet, Int_t jetContNb =0 ){

  AliJetContainer *jetCont = GetJetContainer(jetContNb);
  if (!jet->GetNumberOfTracks())
    return 0;
  Double_t den=0.;
  Double_t num = 0.;
  AliVParticle *vp1 = 0x0;

  std::vector<int> ordindex;
  ordindex=jet->GetPtSortedTrackConstituentIndexes(jetCont->GetParticleContainer()->GetArray());
   if(ordindex.size()<1) return -1;
  
  vp1 = static_cast<AliVParticle*>(jet->TrackAt(ordindex[0], jetCont->GetParticleContainer()->GetArray()));
  if (!vp1){
    Printf("AliVParticle associated to Leading constituent not found");
    return -1;
  }
  
   
  num=vp1->Pt();


  
return num;
} 

//________________________________________________________________________
Float_t AliAnalysisTaskFakeJets::GetJetHardTrack(AliEmcalJet *jet, Int_t jetContNb =0) {
  
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedLeSub();
      else return jet->GetShapeProperties()->GetSecondOrderSubtractedLeSub();
  else
    return HardTrack(jet, jetContNb);
 
}


//________________________________________________________________________
Float_t AliAnalysisTaskFakeJets::SecHardTrack(AliEmcalJet *jet, Int_t jetContNb =0 ){

  AliJetContainer *jetCont = GetJetContainer(jetContNb);
  if (!jet->GetNumberOfTracks())
    return 0;
  Double_t den=0.;
  Double_t num = 0.;
  AliVParticle *vp1 = 0x0;
  AliVParticle *vp2 = 0x0;
  std::vector<int> ordindex;
  ordindex=jet->GetPtSortedTrackConstituentIndexes(jetCont->GetParticleContainer()->GetArray());
 
  
  if(ordindex.size()<2) return -1;
  
   
  vp2 = static_cast<AliVParticle*>(jet->TrackAt(ordindex[1], jetCont->GetParticleContainer()->GetArray()));
  if (!vp2){
    Printf("AliVParticle associated to Subleading constituent not found");
    return -1;
  }
  
  
    den=vp2->Pt();
  //Printf("vp1->Pt() =%f, vp2->Pt() =%f", vp1->Pt(), vp2->Pt());
  
return den;
} 

//________________________________________________________________________
Float_t AliAnalysisTaskFakeJets::GetJetSecHardTrack(AliEmcalJet *jet, Int_t jetContNb =0) {
  //calc subtracted jet mass
 
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedLeSub();
      else return jet->GetShapeProperties()->GetSecondOrderSubtractedLeSub();
  else
    return SecHardTrack(jet, jetContNb);
 
}

//----------------------------------------------------------------------
Double_t AliAnalysisTaskFakeJets::NSubJettiness(AliEmcalJet *Jet, Int_t JetContNb, Double_t JetRadius,  AliEmcalJetFinder *Reclusterer, Int_t N, Int_t A, Int_t B){
  AliJetContainer *JetCont = GetJetContainer(JetContNb);
  AliEmcalJet *SubJet=NULL;
  Double_t DeltaR1=0;
  Double_t DeltaR2=0;
  AliVParticle *JetParticle=0x0;
  Double_t SubJetiness_Numerator = 0;
  Double_t SubJetiness_Denominator = 0;
  Double_t Index=-2;
  Bool_t Error=kFALSE;
  if (!Jet->GetNumberOfTracks())    return -2;                                                                     
  if (Reclusterer->GetNumberOfJets() < N) return -2;
  for (Int_t i=0; i< Jet->GetNumberOfTracks(); i++){  //loops through all tracks (particles in the jet
    JetParticle = static_cast<AliVParticle*>(Jet->TrackAt(i, JetCont->GetParticleContainer()->GetArray()));       
    for (Int_t j=1; j<=N; j++){
      Index=SubJetOrdering(Jet,Reclusterer,j,0,kTRUE);
      if(Index==-999){
	Error = kTRUE;
	i=Jet->GetNumberOfTracks();
	break;
      }
      if(j==1){
	DeltaR1=TMath::Power((Reclusterer->GetJet(Index)->Pt()),A)*TMath::Power((TMath::Sqrt((((JetParticle->Eta())-(Reclusterer->GetJet(Index)->Eta()))*((JetParticle->Eta())- (Reclusterer->GetJet(Index)->Eta())))+((RelativePhi((Reclusterer->GetJet(Index)->Phi()),JetParticle->Phi()))*(RelativePhi((Reclusterer->GetJet(Index)->Phi()),JetParticle->Phi()))))),B);
      }
      else{
	DeltaR2=TMath::Power((Reclusterer->GetJet(Index)->Pt()),A)*TMath::Power((TMath::Sqrt((((JetParticle->Eta())-(Reclusterer->GetJet(Index)->Eta()))*((JetParticle->Eta())- (Reclusterer->GetJet(Index)->Eta())))+((RelativePhi((Reclusterer->GetJet(Index)->Phi()),JetParticle->Phi()))*(RelativePhi((Reclusterer->GetJet(Index)->Phi()),JetParticle->Phi()))))),B);
	if (DeltaR2<DeltaR1) DeltaR1=DeltaR2;
      }
    }
    SubJetiness_Numerator=SubJetiness_Numerator+(JetParticle->Pt()*DeltaR1);    
    if (A>=0) SubJetiness_Denominator=SubJetiness_Denominator+(TMath::Power((Reclusterer->GetJet(SubJetOrdering(Jet,Reclusterer,1,0,kTRUE))->Pt()),A)*JetParticle->Pt()*TMath::Power(JetRadius,B));       
    else  SubJetiness_Denominator=SubJetiness_Denominator+(TMath::Power((Reclusterer->GetJet(SubJetOrdering(Jet,Reclusterer,N,0,kTRUE))->Pt()),A)*JetParticle->Pt()*TMath::Power(JetRadius,B));  
  }
  if (SubJetiness_Denominator!=0 && !Error) return SubJetiness_Numerator/SubJetiness_Denominator;                                                                                  
  else return -2;
}
//____________________________________________________________________________

//----------------------------------------------------------------------
Double_t AliAnalysisTaskFakeJets::GetSubjetFraction(AliEmcalJet *Jet, Int_t JetContNb, Double_t JetRadius,  AliEmcalJetFinder *Reclusterer){
  AliJetContainer *JetCont = GetJetContainer(JetContNb);
  AliEmcalJet *SubJet=NULL;
  Double_t DeltaR1=0;
  Double_t DeltaR2=0;
  AliVParticle *JetParticle=0x0;
  Double_t SubJetiness_Numerator = 0;
  Double_t SubJetiness_Denominator = 0;
  Double_t Index=-2;
  if (Reclusterer->GetNumberOfJets() < 1) return -2;
  Index=SubJetOrdering(Jet,Reclusterer,1,0,kTRUE);
  if(Index==-999) return -2;
   SubJetiness_Numerator=(Reclusterer->GetJet(Index)->Pt());
  SubJetiness_Denominator=Jet->Pt();  
   return SubJetiness_Numerator/SubJetiness_Denominator; 

  
}
//__________________________________________________________________________________





Float_t AliAnalysisTaskFakeJets::CoreFrac(AliEmcalJet *jet, Int_t jetContNb = 0){

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
      if(dr<=fSubjetRadius) num=num+vp1->Pt();
      
    }
    return num/jet->Pt();
} 




//________________________________________________________________________
Float_t AliAnalysisTaskFakeJets::GetJetCoreFrac(AliEmcalJet *jet, Int_t jetContNb =0) {
  //calc subtracted jet mass
 
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedLeSub();
      else return jet->GetShapeProperties()->GetSecondOrderSubtractedLeSub();
  else
    return CoreFrac(jet, jetContNb);
 
}




//----------------------------------------------------------------------
Double_t AliAnalysisTaskFakeJets::SubJetOrdering(AliEmcalJet *Jet, AliEmcalJetFinder *Reclusterer, Int_t N, Int_t Type, Bool_t Index){
  AliEmcalJet *SubJet=NULL;
  Double_t SortingVariable;
  Int_t ArraySize =N+1;
  TArrayD *JetSorter = new TArrayD(ArraySize);
  TArrayD *JetIndexSorter = new TArrayD(ArraySize);
  for (Int_t i=0; i<ArraySize; i++){
    JetSorter->SetAt(0,i);
  }
  for (Int_t i=0; i<ArraySize; i++){
    JetIndexSorter->SetAt(0,i);
  }
  if(Reclusterer->GetNumberOfJets()<N) return -999;
  for (Int_t i=0; i<Reclusterer->GetNumberOfJets(); i++){
    SubJet=Reclusterer->GetJet(i);
    if (Type==0) SortingVariable=SubJet->Pt();
    else if (Type==1) SortingVariable=SubJet->E();
    else if (Type==2) SortingVariable=SubJet->M();
    for (Int_t j=0; j<N; j++){
      if (SortingVariable>JetSorter->GetAt(j)){
	for (Int_t k=N-1; k>=j; k--){
	  JetSorter->SetAt(JetSorter->GetAt(k),k+1);
	  JetIndexSorter->SetAt(JetIndexSorter->GetAt(k),k+1);
	}
	JetSorter->SetAt(SortingVariable,j);
	JetIndexSorter->SetAt(i,j);
	break;
      }
    }
  }
  if (!Index) return JetSorter->GetAt(N-1);
  else return JetIndexSorter->GetAt(N-1);
}



//returns -1 if the Nth hardest jet is requested where N>number of available jets
//type:  0=Pt  1=E  2=M 
//Index TRUE=returns index   FALSE=returns value of quantatiy in question












//________________________________________________________________________
Int_t AliAnalysisTaskFakeJets::SelectTrigger(Float_t minpT, Float_t maxpT){

  AliParticleContainer *partCont = GetParticleContainer(0);
  TClonesArray *tracksArray = partCont->GetArray();
  
  if(!partCont || !tracksArray) return -99999;
  AliAODTrack *track = 0x0;
  AliEmcalParticle *emcPart = 0x0;

  
  TList *trackList = new TList();
  Int_t triggers[100];
  for (Int_t iTrigger=0; iTrigger<100; iTrigger++) triggers[iTrigger] = 0;
  Int_t iTT = 0;
  
  for(Int_t iTrack=0; iTrack <= tracksArray->GetEntriesFast(); iTrack++){
    
   
    if (fJetShapeSub == kConstSub){
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

//----------------------------------------------------------------------
///////////////returns jet finder object containg subjets
AliEmcalJetFinder *AliAnalysisTaskFakeJets::Recluster(AliEmcalJet *Jet, Int_t JetContNb, Double_t SubJetRadius, Double_t SubJetMinPt, Int_t Algorithm, const char* Name){

  AliJetContainer *JetCont = GetJetContainer(JetContNb);
  AliEmcalJetFinder *Reclusterer = new AliEmcalJetFinder(Name); //JetFinder Object for reclustered jets                                                                 
  Reclusterer->SetRadius(SubJetRadius); 
  Reclusterer->SetJetMinPt(SubJetMinPt);
  Reclusterer->SetJetAlgorithm(Algorithm); //0 for anti-kt     1 for kt   
  const AliVVertex *vert = InputEvent()->GetPrimaryVertex();
  Double_t dVtx[3]={vert->GetX(),vert->GetY(),vert->GetZ()};
  if(Reclusterer->AliEmcalJetFinder::Filter(Jet, JetCont, dVtx)){;}  //reclustering jet1 using the jetfinderobject Reclusterer                            
  return Reclusterer;
}




//__________________________________________________________________________________
Double_t AliAnalysisTaskFakeJets::RelativePhi(Double_t mphi,Double_t vphi){

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
Bool_t AliAnalysisTaskFakeJets::RetrieveEventObjects() {
  //
  // retrieve event objects
  //
  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  return kTRUE;
}

//_______________________________________________________________________
void AliAnalysisTaskFakeJets::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.

  // fTreeObservableTagging = dynamic_cast<TTree*>(GetOutputData(1));
  // if (!fTreeObservableTagging){
  //   Printf("ERROR: fTreeObservableTagging not available"); 
  //   return;
  // }

}

