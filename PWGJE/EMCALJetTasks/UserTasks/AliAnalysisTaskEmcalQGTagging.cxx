//
// Jet QG tagging analysis task.
//
// Author: D. Caffarri, L. Cunqueiro 

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


#include "AliAODEvent.h"

#include "AliAnalysisTaskEmcalQGTagging.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskEmcalQGTagging)

//________________________________________________________________________
AliAnalysisTaskEmcalQGTagging::AliAnalysisTaskEmcalQGTagging() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalQGTagging", kTRUE),
  fContainer(0),
  fMinFractionShared(0),
  fJetShapeType(kData),
  fJetShapeSub(kNoSub),
  fJetSelection(kInclusive),
  fShapesVar(0),
  fPtThreshold(-9999.),
  fRMatching(0.3),
  fminpTTrig(20.),
  fmaxpTTrig(50.),
  fangWindowRecoil(0.6),
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
  fTreeObservableTagging(0)

{
  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalQGTagging::AliAnalysisTaskEmcalQGTagging(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fContainer(0),
  fMinFractionShared(0),
  fJetShapeType(kData),
  fJetShapeSub(kNoSub),
  fJetSelection(kInclusive),
  fShapesVar(0),
  fPtThreshold(-9999.),
  fRMatching(0.3),
  fminpTTrig(20.),
  fmaxpTTrig(50.),
  fangWindowRecoil(0.6),
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
  fTreeObservableTagging(0)
  
{
  // Standard constructor.
  
  SetMakeGeneralHistograms(kTRUE);

  DefineOutput(1, TTree::Class());

}

//________________________________________________________________________
AliAnalysisTaskEmcalQGTagging::~AliAnalysisTaskEmcalQGTagging()
{
  // Destructor.
}

//________________________________________________________________________
 void AliAnalysisTaskEmcalQGTagging::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fTreeObservableTagging = new TTree("fTreeJetShape", "fTreeJetShape");
  Int_t nVar = 18; 
  fShapesVar = new Float_t [nVar]; 
  TString *fShapesVarNames = new TString [nVar];

  fShapesVarNames[0] = "partonCode"; 
  fShapesVarNames[1] = "ptJet"; 
  fShapesVarNames[2] = "ptDJet"; 
  fShapesVarNames[3] = "mJet";
  fShapesVarNames[4] = "nbOfConst";
  fShapesVarNames[5] = "angularity";
  fShapesVarNames[6] = "circularity";
  fShapesVarNames[7] = "lesub";
  fShapesVarNames[8] = "sigma2";

  fShapesVarNames[9] = "ptJetMatch"; 
  fShapesVarNames[10] = "ptDJetMatch"; 
  fShapesVarNames[11] = "mJetMatch";
  fShapesVarNames[12] = "nbOfConstMatch";
  fShapesVarNames[13] = "angularityMatch";
  fShapesVarNames[14] = "circularityMatch";
  fShapesVarNames[15] = "lesubMatch";
  fShapesVarNames[16] = "sigma2Match";
  fShapesVarNames[17]="weightPythia";

  for(Int_t ivar=0; ivar < nVar; ivar++){
    cout<<"looping over variables"<<endl;
    fTreeObservableTagging->Branch(fShapesVarNames[ivar].Data(), &fShapesVar[ivar], Form("%s/F", fShapesVarNames[ivar].Data()));

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
  
  fhpTjetpT= new TH2F("fhpTjetpT", "fhpTjetpT", 100, 0, 200,  100, 0, 200);
  fOutput->Add(fhpTjetpT);
  fhPt= new TH1F("fhPt", "fhPt", 100, 0, 200);
  fOutput->Add(fhPt);
  fhPhi= new TH1F("fhPhi", "fhPhi", 100, -TMath::Pi(), TMath::Pi());
  fOutput->Add(fhPhi);
  
  fOutput->Add(fTreeObservableTagging);
  TH1::AddDirectory(oldStatus);
  PostData(1, fOutput); // Post data for ALL output slots > 0 here.

}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalQGTagging::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalQGTagging::FillHistograms()
{
  // Fill histograms.
  //cout<<"base container"<<endl;
  AliEmcalJet* jet1 = NULL;
  AliJetContainer *jetCont = GetJetContainer(0);
  Float_t kWeight=1;
  if(fCent>10) return 0;

  AliAODTrack *triggerHadron = 0x0;
  
  if (fJetSelection == kRecoil) {
    //Printf("Recoil jets!!!, fminpTTrig = %f, fmaxpTTrig = %f", fminpTTrig, fmaxpTTrig);
    Int_t triggerHadronLabel = SelectTrigger(fminpTTrig, fmaxpTTrig);
    if (triggerHadronLabel==-9999) {
      Printf ("Trigger Hadron not found, return");
      return 0;
    }
    AliParticleContainer *partContAn = GetParticleContainer(0);
    TClonesArray *trackArrayAn = partContAn->GetArray();
    triggerHadron = static_cast<AliAODTrack*>(trackArrayAn->At(triggerHadronLabel));
    if (!triggerHadron) {
      Printf("No Trigger hadron with the found label!!");
      return 0;
    }
  }
  
  if(jetCont) {
    jetCont->ResetCurrentID();
    while((jet1 = jetCont->GetNextAcceptJet())) {
      if (!jet1) continue;
      AliEmcalJet* jet2 = 0x0;
      fPtJet->Fill(jet1->Pt());
      AliEmcalJet *jetUS = NULL;
      Int_t ifound=0;
      Int_t ilab=-1;
      
      if (!(fJetShapeType == kData)) {
        AliPythiaInfo *partonsInfo = 0x0;
        
        if((fJetShapeType == kTrueDet) || (fJetShapeType == kDetEmb)){
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
            Printf("jet2 not exists, returning");
            continue;
          }
          
          fh2ResponseUW->Fill(jet1->Pt(),jet2->Pt());
          
          Double_t fraction=0;
          if(!(fJetShapeSub==kConstSub))  fraction = jetCont->GetFractionSharedPt(jet1);
          if(fJetShapeSub==kConstSub) fraction = jetContUS->GetFractionSharedPt(jetUS);
          cout<<"hey a jet"<<fraction<<" "<<jet1->Pt()<<" "<<jet2->Pt()<<" "<<fCent<<endl;
          
          if(fraction<fMinFractionShared) continue;
          partonsInfo = (AliPythiaInfo*) jetContTrue->GetPythiaInfo();
          if(!partonsInfo) return 0;
          
        }
        else {
          partonsInfo = (AliPythiaInfo*) jetCont->GetPythiaInfo();
          jet2=jet1;
          if(!partonsInfo) return 0;
        }
        
        Double_t jp1=(jet2->Phi())-(partonsInfo->GetPartonPhi6());
        Double_t detap1=(jet2->Eta())-(partonsInfo->GetPartonEta6());
        kWeight=partonsInfo->GetPythiaEventWeight();
        fh2ResponseW->Fill(jet1->Pt(),jet2->Pt(),kWeight);
        if (jp1< -1*TMath::Pi()) jp1 = (-2*TMath::Pi())-jp1;
        else if (jp1 > TMath::Pi()) jp1 = (2*TMath::Pi())-jp1;
        Float_t dRp1 = TMath::Sqrt(jp1 * jp1 + detap1 * detap1);
        fEtaJetCorr6->Fill(jet2->Eta(), partonsInfo->GetPartonEta6());
        fPhiJetCorr6->Fill(jet2->Phi(), partonsInfo->GetPartonPhi6());
        if(dRp1 < fRMatching) {
          fShapesVar[0] = partonsInfo->GetPartonFlag6();
          fPtJetCorr ->Fill(partonsInfo->GetPartonPt6(), jet2->Pt());
        }
        else {
          jp1=(jet2->Phi())-(partonsInfo->GetPartonPhi7());
          detap1=(jet2->Eta())-(partonsInfo->GetPartonEta7());
          if (jp1< -1*TMath::Pi()) jp1= (-2*TMath::Pi())-jp1;
          else if (jp1 > TMath::Pi()) jp1 = (2*TMath::Pi())-jp1;
          dRp1 = TMath::Sqrt(jp1 * jp1 + detap1 * detap1);
          fEtaJetCorr7->Fill(jet2->Eta(), partonsInfo->GetPartonEta7());
          fPhiJetCorr7->Fill(jet2->Phi(), partonsInfo->GetPartonPhi7());
          if(dRp1 < fRMatching) {
            fShapesVar[0] = partonsInfo->GetPartonFlag7();
            fPtJetCorr->Fill(partonsInfo->GetPartonPt7(), jet2->Pt());
          }
          else continue;
        }
      }
      else
        fShapesVar[0] = 0.;
      
      Double_t ptSubtracted = 0;
      
      Float_t dphiRecoil = 0.;
      if (fJetSelection == kRecoil){
        dphiRecoil = RelativePhi(triggerHadron->Phi(), jet1->Phi());
        if (TMath::Abs(dphiRecoil) < TMath::Pi() - fangWindowRecoil) {
         // Printf("Recoil jets back to back not found! continuing");
          continue;
        }
        //Printf("Recoil jets back to back found! filling histos!");
        fhpTjetpT->Fill(triggerHadron->Pt(), jet1->Pt());
        //Printf("triggerHadron =%f, jet1 = %f", triggerHadron->Pt(), jet1->Pt());
        fhPt->Fill(triggerHadron->Pt());
        fhPhi->Fill(RelativePhi(triggerHadron->Phi(), jet1->Phi()));
      }
      
      if (((fJetShapeType == kData) || (fJetShapeType == kDetEmb)) && (fJetShapeSub == kConstSub))
        ptSubtracted = jet1->Pt();
      else ptSubtracted  = jet1->Pt() - GetRhoVal(0)*jet1->Area();
      
      if ((fJetShapeType == kData || fJetShapeType== kDetEmb))
        if (ptSubtracted < fPtThreshold) continue;
      
      fShapesVar[1] = ptSubtracted;
      fShapesVar[2] = GetJetpTD(jet1,0);
      fShapesVar[3] = GetJetMass(jet1,0);
      fShapesVar[4] = 1.*GetJetNumberOfConstituents(jet1,0);
      fShapesVar[5] = GetJetAngularity(jet1,0);
      fShapesVar[6] = GetJetCircularity(jet1,0);
      fShapesVar[7] = GetJetLeSub(jet1,0);
      fShapesVar[8] = GetSigma2(jet1,0);
      
      Float_t ptMatch=0., ptDMatch=0., massMatch=0., constMatch=0.,angulMatch=0.,circMatch=0., lesubMatch=0., sigma2Match=0.;
      Int_t kMatched = 0;
      if (fJetShapeType == kTrueDet || fJetShapeType == kDetEmb) {
        kMatched = 1;
        ptMatch=jet2->Pt();
        ptDMatch=GetJetpTD(jet2, kMatched);
        massMatch=GetJetMass(jet2,kMatched);
        constMatch=1.*GetJetNumberOfConstituents(jet2,kMatched);
        angulMatch=GetJetAngularity(jet2, kMatched);
        circMatch=GetJetCircularity(jet2, kMatched);
        lesubMatch=GetJetLeSub(jet2, kMatched);
        sigma2Match = GetSigma2(jet2, kMatched);
        
      }
      
      if (fJetShapeType == kTrue || fJetShapeType == kData) {
        kMatched = 0;
        ptMatch=0.;
        ptDMatch=0.;
        massMatch=0.;
        constMatch=0.;
        angulMatch=0.;
        circMatch=0.;
        lesubMatch=0.;
        sigma2Match =0.;
        
      }
      
      fShapesVar[9] = ptMatch;
      fShapesVar[10] = ptDMatch;
      fShapesVar[11] = massMatch;
      fShapesVar[12] = constMatch;
      fShapesVar[13] = angulMatch;
      fShapesVar[14] = circMatch;
      fShapesVar[15] = lesubMatch;
      fShapesVar[16] = sigma2Match;
      fShapesVar[17] = kWeight;
      fTreeObservableTagging->Fill();
      
    }
    
  }
  
  return kTRUE;
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalQGTagging::GetJetMass(AliEmcalJet *jet,Int_t jetContNb=0) {
  //calc subtracted jet mass
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    return jet->GetSecondOrderSubtracted();
  else 
    return jet->M();
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalQGTagging::Angularity(AliEmcalJet *jet, Int_t jetContNb = 0){

  AliJetContainer *jetCont = GetJetContainer(jetContNb);
  if (!jet->GetNumberOfTracks())
      return 0; 
    Double_t den=0.;
    Double_t num = 0.;
    AliVParticle *vp1 = 0x0;
    for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
      vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));  
      Double_t dphi = vp1->Phi()-jet->Phi();
      if(dphi<-1.*TMath::Pi()) dphi+=TMath::TwoPi();
      if(dphi>TMath::Pi()) dphi-=TMath::TwoPi();
      Double_t dr2 = (vp1->Eta()-jet->Eta())*(vp1->Eta()-jet->Eta()) + dphi*dphi;
      Double_t dr = TMath::Sqrt(dr2);
      num=num+vp1->Pt()*dr;
      den=den+vp1->Pt();
    }
    return num/den;
} 

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalQGTagging::GetJetAngularity(AliEmcalJet *jet, Int_t jetContNb = 0) {

  if((fJetShapeSub==kDerivSub) && (jetContNb==0))
    return jet->GetSecondOrderSubtractedAngularity();
  else
    return Angularity(jet, jetContNb);
 
}


//________________________________________________________________________
Float_t AliAnalysisTaskEmcalQGTagging::PTD(AliEmcalJet *jet, Int_t jetContNb = 0){

  AliJetContainer *jetCont = GetJetContainer(jetContNb);
  if (!jet->GetNumberOfTracks())
      return 0; 
    Double_t den=0.;
    Double_t num = 0.;
    AliVParticle *vp1 = 0x0;
    for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
      vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));  
      num=num+vp1->Pt()*vp1->Pt();
      den=den+vp1->Pt();
    }
    return TMath::Sqrt(num)/den;
} 

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalQGTagging::GetJetpTD(AliEmcalJet *jet, Int_t jetContNb = 0) {
  //calc subtracted jet mass
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    return jet->GetSecondOrderSubtractedpTD();
  else
    return PTD(jet, jetContNb);
 
}

//_____________________________________________________________________________
Float_t AliAnalysisTaskEmcalQGTagging::Circularity(AliEmcalJet *jet, Int_t jetContNb = 0){

  AliJetContainer *jetCont = GetJetContainer(jetContNb);
  if (!jet->GetNumberOfTracks())
    return 0; 
  Double_t mxx    = 0.;
  Double_t myy    = 0.;
  Double_t mxy    = 0.;
  int  nc     = 0;
  Double_t sump2  = 0.;
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
Float_t AliAnalysisTaskEmcalQGTagging::GetJetCircularity(AliEmcalJet *jet, Int_t jetContNb =0 ) {
  //calc subtracted jet mass
 
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    return jet->GetSecondOrderSubtractedCircularity();
  else
    return Circularity(jet, jetContNb);
 
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalQGTagging::LeSub(AliEmcalJet *jet, Int_t jetContNb =0 ){

  AliJetContainer *jetCont = GetJetContainer(jetContNb);
  if (!jet->GetNumberOfTracks())
      return 0; 
    Double_t den=0.;
    Double_t num = 0.;
    AliVParticle *vp1 = 0x0;
    AliVParticle *vp2 = 0x0;
    std::vector<int> ordindex;
    ordindex=jet->SortConstituentsPt(jetCont->GetParticleContainer()->GetArray());
    
   vp1 = static_cast<AliVParticle*>(jet->TrackAt(ordindex[0], jetCont->GetParticleContainer()->GetArray()));  
   vp2 = static_cast<AliVParticle*>(jet->TrackAt(ordindex[1], jetCont->GetParticleContainer()->GetArray()));  
     
  num=vp1->Pt();
  den=vp2->Pt();
  
return num-den;
} 

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalQGTagging::GetJetLeSub(AliEmcalJet *jet, Int_t jetContNb =0) {
  //calc subtracted jet mass
 
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    return jet->GetSecondOrderSubtractedLeSub();
  else
    return LeSub(jet, jetContNb);
 
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalQGTagging::GetJetNumberOfConstituents(AliEmcalJet *jet,Int_t jetContNb=0) {
  //calc subtracted jet mass
  
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    return jet->GetSecondOrderSubtractedConstituent();
  else
    return jet->GetNumberOfTracks();
 
}
   

//______________________________________________________________________________
Float_t AliAnalysisTaskEmcalQGTagging::Sigma2(AliEmcalJet *jet, Int_t jetContNb=0){

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
       Double_t ppt=vp1->Pt();
       Double_t dphi = vp1->Phi()-jet->Phi();
       if(dphi<-1.*TMath::Pi()) dphi+=TMath::TwoPi();
       if(dphi>TMath::Pi()) dphi-=TMath::TwoPi();
       Double_t deta = vp1->Eta()-jet->Eta();
       mxx += ppt*ppt*deta*deta;
       myy += ppt*ppt*dphi*dphi;
       mxy -= ppt*ppt*deta*dphi;
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
Float_t AliAnalysisTaskEmcalQGTagging::GetSigma2(AliEmcalJet *jet, Int_t jetContNb=0) {
  //calc subtracted jet mass
 
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    return jet->GetSecondOrderSubtractedSigma2();
  else
    return Sigma2(jet, jetContNb);
 
}

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalQGTagging::SelectTrigger(Float_t minpT, Float_t maxpT){

  AliParticleContainer *partCont = GetParticleContainer(0);
  TClonesArray *tracksArray = partCont->GetArray();

  if(!partCont || !tracksArray) return 0; 
  AliAODTrack *track = 0x0;

  TList *trackList = new TList();
  Int_t triggers[100];
  for (Int_t iTrigger=0; iTrigger<100; iTrigger++) triggers[iTrigger] = 0; 
  Int_t iTT = 0; 
    
  for(Int_t iTrack=0; iTrack <= tracksArray->GetEntriesFast(); iTrack++){
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

  if (iTT == 0) return -9999;
  Int_t nbRn = 0, index = 0 ; 
  TRandom3* random = new TRandom3(0); 
  nbRn = random->Integer(iTT);

  index = triggers[nbRn];
  //Printf("iTT= %d, nbRn = %d, Index = %d",iTT, nbRn, index );
  return index; 
  
}

//__________________________________________________________________________________
Double_t AliAnalysisTaskEmcalQGTagging::RelativePhi(Double_t mphi,Double_t vphi){

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
Bool_t AliAnalysisTaskEmcalQGTagging::RetrieveEventObjects() {
  //
  // retrieve event objects
  //
  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  return kTRUE;
}

//_______________________________________________________________________
void AliAnalysisTaskEmcalQGTagging::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.

  // fTreeObservableTagging = dynamic_cast<TTree*>(GetOutputData(1));
  // if (!fTreeObservableTagging){
  //   Printf("ERROR: fTreeObservableTagging not available"); 
  //   return;
  // }

}

