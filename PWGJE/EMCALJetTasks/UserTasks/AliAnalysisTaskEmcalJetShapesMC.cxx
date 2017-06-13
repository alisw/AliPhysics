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
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
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
#include "AliEmcalPythiaInfo.h"
#include "TRandom3.h"
#include "AliEmcalJetFinder.h"
#include "AliAODEvent.h"
#include "AliAnalysisTaskEmcalJetShapesMC.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskEmcalJetShapesMC)

//________________________________________________________________________
AliAnalysisTaskEmcalJetShapesMC::AliAnalysisTaskEmcalJetShapesMC() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalJetShapesMC", kTRUE),
  fContainer(0),
  fMinFractionShared(0),
  fJetShapeType(kGenShapes),
  fJetShapeSub(kNoSub),
  fJetSelection(kInclusive),
  fPtThreshold(-9999.),
  fRMatching(0.2),
  fJetRadius(0.4),
  fSubjetRadius(0.2),
  fSelectedShapes(0),
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
  fTreeObservableTagging(0x0)

{
  for(Int_t i=0;i<33;i++){
    fShapesVar[i]=0;}
  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetShapesMC::AliAnalysisTaskEmcalJetShapesMC(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fContainer(0),
  fMinFractionShared(0),
  fJetShapeType(kGenShapes),
  fJetShapeSub(kNoSub),
  fJetSelection(kInclusive),
  fPtThreshold(-9999.),
  fRMatching(0.2),
  fSelectedShapes(0), 
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
  fTreeObservableTagging(0x0)
  
{
  // Standard constructor.
  
  
  for(Int_t i=0;i<33;i++){
    fShapesVar[i]=0;}
  
  SetMakeGeneralHistograms(kTRUE);

  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  

}

//________________________________________________________________________
AliAnalysisTaskEmcalJetShapesMC::~AliAnalysisTaskEmcalJetShapesMC()
{
  if(fTreeObservableTagging){
    delete fTreeObservableTagging;
    fTreeObservableTagging = 0;
  }

}

//________________________________________________________________________
 void AliAnalysisTaskEmcalJetShapesMC::UserCreateOutputObjects()
{
  // Create user output.
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(oldStatus);

  //fTreeObservableTagging = new TTree("fTreeJetShape", "fTreeJetShape");

  //TH1::AddDirectory(oldStatus);
  
  const char* nameoutput = GetOutputSlot(2)->GetContainer()->GetName();
  fTreeObservableTagging = new TTree(nameoutput, nameoutput);
  
  const Int_t nVar = 33;

  TString *fShapesVarNames = new TString [nVar];
  
  fShapesVarNames[0] = "partonCode";
  fShapesVarNames[1] = "ptJet"; 
  fShapesVarNames[2] = "ptDJet"; 
  fShapesVarNames[3] = "mJet";
  fShapesVarNames[4] = "nbOfConst";
  fShapesVarNames[5] = "angularity";
  fShapesVarNames[6] = "Nsubjet1kt";
  fShapesVarNames[7] = "Nsubjet2kt"; 
  fShapesVarNames[8] = "Nsubjet1Min"; 
  fShapesVarNames[9] = "Nsubjet2Min";
  fShapesVarNames[10] = "DeltaRkt";
  fShapesVarNames[11] = "DeltaRMin";
  fShapesVarNames[12] = "SDSymm";
   fShapesVarNames[13] = "SDDeltaR";
  fShapesVarNames[14] = "SDGroomedFrac"; 
  fShapesVarNames[15] = "SDGroomedN"; 
 fShapesVarNames[16] = "SDSymmkt";
   fShapesVarNames[17] = "SDDeltaRkt";
  fShapesVarNames[18] = "SDGroomedFrackt"; 
  fShapesVarNames[19] = "SDGroomedNkt";
   fShapesVarNames[20] = "SDSymmAkt";
   fShapesVarNames[21] = "SDDeltaRAkt";
  fShapesVarNames[22] = "SDGroomedFracAkt"; 
  fShapesVarNames[23] = "SDGroomedNAkt"; 
  fShapesVarNames[24] = "SDSymmBeta1";
   fShapesVarNames[25] = "SDDeltaRBeta1";
  fShapesVarNames[26] = "SDGroomedFracBeta1"; 
  fShapesVarNames[27] = "SDGroomedNBeta1"; 
   fShapesVarNames[28] = "SDSymmBeta2";
   fShapesVarNames[29] = "SDDeltaRBeta2";
  fShapesVarNames[30] = "SDGroomedFracBeta2"; 
  fShapesVarNames[31] = "SDGroomedNBeta2"; 
  fShapesVarNames[32] = "weightPythia"; 


  
   //fShapesVarNames[7] = "lesub";
  //fShapesVarNames[8] = "CoreFraction";
  //fShapesVarNames[9] = "Nsubjet1";
  //fShapesVarNames[10] = "Nsubjet2";
  //fShapesVarNames[11] = "DeltaR";
  //fShapesVarNames[12] = "OpenAngle";
  //fShapesVarNames[13] = "weightPythia";

  //fShapesVarNames[14] = "NT70";
  //fShapesVarNames[15] = "nConstNT70";
  //fShapesVarNames[16] = "NT80";
  //fShapesVarNames[17] = "nConstNT80";
  //fShapesVarNames[18] = "NT90";
  //fShapesVarNames[19] = "nConstNT90";
  //fShapesVarNames[20] = "NT95";
  //fShapesVarNames[21] = "nConstNT95";
  
  //fShapesVarNames[22] = "SubjetFraction";


   for(Int_t ivar=0; ivar < nVar; ivar++){
    cout<<"looping over variables"<<endl;
    fTreeObservableTagging->Branch(fShapesVarNames[ivar].Data(), &fShapesVar[ivar], Form("%s/F", fShapesVarNames[ivar].Data()));}

  

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
  
  //fOutput->Add(fTreeObservableTagging);
  
 
  PostData(1, fOutput); // Post data for ALL output slots > 0 here
  PostData(2, fTreeObservableTagging);
  
  delete [] fShapesVarNames;

}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetShapesMC::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetShapesMC::FillHistograms()
{
  // Fill histograms.
  //cout<<"IntoFillHistograms"<<endl;
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
      //Printf("jet1=%p", jet1);
      if (!jet1) continue;
      AliEmcalJet* jet2 = 0x0;
      AliEmcalJet* jet3 = 0x0;
      fPtJet->Fill(jet1->Pt());
      AliEmcalJet *jetUS = NULL;
      Int_t ifound=0, jfound=0;
      Int_t ilab=-1, jlab=-1;
      
      if(fSemigoodCorrect && (fJetSelection != kRecoil)){
        Double_t disthole=RelativePhi(jet1->Phi(),fHolePos);
        if(TMath::Abs(disthole)<fHoleWidth){
          continue;
        }
      }
      
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
      
      
      fShapesVar[0] = 0.;
      
      if (fJetShapeType == kGenShapes){
        const AliEmcalPythiaInfo *partonsInfo = 0x0;
        partonsInfo = GetPythiaInfo();
        //Printf("partonsInfo=%p",  partonsInfo);
        Double_t jp1=RelativePhi(jet1->Phi(),partonsInfo->GetPartonPhi6());
        Double_t detap1=(jet1->Eta())-(partonsInfo->GetPartonEta6());
        kWeight=partonsInfo->GetPythiaEventWeight();
        //Printf("kWeight=%f",  kWeight);
        fShapesVar[32] = kWeight;
        
        Float_t dRp1 = TMath::Sqrt(jp1 * jp1 + detap1 * detap1);
        fEtaJetCorr6->Fill(jet1->Eta(), partonsInfo->GetPartonEta6());
        fPhiJetCorr6->Fill(jet1->Phi(), partonsInfo->GetPartonPhi6());
        if(dRp1 < fRMatching) {
          fShapesVar[0] = partonsInfo->GetPartonFlag6();
          fPtJetCorr ->Fill(partonsInfo->GetPartonPt6(), jet1->Pt());
        }
        else {
          jp1=RelativePhi(jet1->Phi(),partonsInfo->GetPartonPhi7());
          detap1=(jet1->Eta())-(partonsInfo->GetPartonEta7());
          dRp1 = TMath::Sqrt(jp1 * jp1 + detap1 * detap1);
          fEtaJetCorr7->Fill(jet1->Eta(), partonsInfo->GetPartonEta7());
          fPhiJetCorr7->Fill(jet1->Phi(), partonsInfo->GetPartonPhi7());
          if(dRp1 < fRMatching) {
            fShapesVar[0] = partonsInfo->GetPartonFlag7();
            fPtJetCorr->Fill(partonsInfo->GetPartonPt7(), jet1->Pt());
          }
          else fShapesVar[0]=0;
        }
      }
      
      Double_t ptSubtracted = 0;
      ptSubtracted= jet1->Pt();
      //Printf("ptSubtracted=%f", ptSubtracted);
      
    
      if (ptSubtracted < fPtThreshold) continue;
    
      if ((fCentSelectOn == kFALSE) && (jet1->GetNumberOfTracks() < 2)) continue;
      
      //AliEmcalJetFinder *Reclusterer1; //Object containg Subjets from Subtracted Hybrid Jets
      //Reclusterer1 = Recluster(jet1, 0, fJetRadius, fSubjetRadius, 1, 0, "SubJetFinder_1");
      

  




      fShapesVar[1] = ptSubtracted;
      fShapesVar[2] = GetJetpTD(jet1,0);
      fShapesVar[3] = GetJetMass(jet1,0);
      fShapesVar[4] = 1.*GetJetNumberOfConstituents(jet1,0);
      fShapesVar[5] = GetJetAngularity(jet1,0);
      //nsub1 and nsub2 for kT
      fShapesVar[6] = fjNSubJettiness(jet1,0,1,0,1,0);
      fShapesVar[7] = fjNSubJettiness(jet1,0,2,0,1,0);
      //nsub1 and nsub2 for min_axis
      fShapesVar[8] =0;
	// 	fjNSubJettiness(jet1,0,1,10,1,0);
      fShapesVar[9] =0;
	// 	fjNSubJettiness(jet1,0,2,10,1,0);
      //deltaRkt
	fShapesVar[10] = 0;
      //fjNSubJettiness(jet1,0,2,0,1,2);
      //deltaRmin
      fShapesVar[11] = 0;
      //fjNSubJettiness(jet1,0,2,10,1,2);
    
      //SoftDropParameters for different reclustering strategies and beta values 
      SoftDrop(jet1,jetCont,0.1,0,0);
      SoftDrop(jet1,jetCont,0.1,0,1);
      SoftDrop(jet1,jetCont,0.1,0,2);
      SoftDrop(jet1,jetCont,0.1,1,0); 
      SoftDrop(jet1,jetCont,0.1,2,0); 
      
          
      // Float_t nTFractions[8]={0.,0.,0.,0.,0.,0.,0.,0.};
      //NTValues(jet1, 0, nTFractions);
      //shape 13 is pythia weight!
      //for (Int_t ishape=14; ishape<22; ishape++) fShapesVar[ishape] = nTFractions[ishape-14];
    
      //fShapesVar[22]= GetSubjetFraction(jet1,0,fJetRadius,Reclusterer1);
      
      fTreeObservableTagging->Fill();

    }

  }
  
  return kTRUE;
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetShapesMC::GetJetMass(AliEmcalJet *jet,Int_t jetContNb){
  //calc subtracted jet mass
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtracted();
      else return jet->GetShapeProperties()->GetSecondOrderSubtracted();
  else 
    return jet->M();
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetShapesMC::Angularity(AliEmcalJet *jet, Int_t jetContNb){

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
Float_t AliAnalysisTaskEmcalJetShapesMC::GetJetAngularity(AliEmcalJet *jet, Int_t jetContNb){

  if((fJetShapeSub==kDerivSub) && (jetContNb==0))
    if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedAngularity();
      else return jet->GetShapeProperties()->GetSecondOrderSubtractedAngularity();
  else
    return Angularity(jet, jetContNb);
 
}


//________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetShapesMC::PTD(AliEmcalJet *jet, Int_t jetContNb){

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
Float_t AliAnalysisTaskEmcalJetShapesMC::GetJetpTD(AliEmcalJet *jet, Int_t jetContNb){
  //calc subtracted jet mass
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedpTD();
      else return jet->GetShapeProperties()->GetSecondOrderSubtractedpTD();
  else
    return PTD(jet, jetContNb);
 
}

//_____________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetShapesMC::Circularity(AliEmcalJet *jet, Int_t jetContNb){

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
Float_t AliAnalysisTaskEmcalJetShapesMC::GetJetCircularity(AliEmcalJet *jet, Int_t jetContNb){
  //calc subtracted jet mass
 
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedCircularity();
      else return jet->GetShapeProperties()->GetSecondOrderSubtractedCircularity();
  else
    return Circularity(jet, jetContNb);
 
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetShapesMC::LeSub(AliEmcalJet *jet, Int_t jetContNb){

  AliJetContainer *jetCont = GetJetContainer(jetContNb);
  if (!jet->GetNumberOfTracks())
    return 0;
  Double_t den=0.;
  Double_t num = 0.;
  AliVParticle *vp1 = 0x0;
  AliVParticle *vp2 = 0x0;
  std::vector<int> ordindex;
  ordindex=jet->GetPtSortedTrackConstituentIndexes(jetCont->GetParticleContainer()->GetArray());
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
Float_t AliAnalysisTaskEmcalJetShapesMC::GetJetLeSub(AliEmcalJet *jet, Int_t jetContNb) {
  //calc subtracted jet mass
 
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedLeSub();
      else return jet->GetShapeProperties()->GetSecondOrderSubtractedLeSub();
  else
    return LeSub(jet, jetContNb);
 
}

//________________________________________________________________________
 Float_t AliAnalysisTaskEmcalJetShapesMC::GetJetNumberOfConstituents(AliEmcalJet *jet,Int_t jetContNb){
  //calc subtracted jet mass
  
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedConstituent();
      else return jet->GetShapeProperties()->GetSecondOrderSubtractedConstituent();
  else
    return jet->GetNumberOfTracks();
 
 }
   

//________________________________________________________________________
AliEmcalJetFinder *AliAnalysisTaskEmcalJetShapesMC::Recluster(AliEmcalJet *Jet, Int_t JetContNb, Double_t JetRadius, Double_t SubJetRadius, Double_t SubJetMinPt, Int_t Algorithm, const char* Name){
  
  AliJetContainer *JetCont = GetJetContainer(JetContNb);
  AliEmcalJetFinder *Reclusterer = new AliEmcalJetFinder(Name); //JetFinder Object for reclustered jets
  Reclusterer->SetRadius(SubJetRadius);
  Reclusterer->SetJetMinPt(SubJetMinPt);
  Reclusterer->SetJetAlgorithm(Algorithm); //0 for anti-kt     1 for kt
  Reclusterer->SetJetMaxEta(0.9-JetRadius);
  Reclusterer->SetRecombSheme(0);
  const AliVVertex *vert = InputEvent()->GetPrimaryVertex();

  //Double_t dVtx[3]={vert->GetX(),vert->GetY(),vert->GetZ()};
  Double_t dVtx[3]={0.,0.,0.};
  if(Reclusterer->AliEmcalJetFinder::Filter(Jet, JetCont, dVtx)){;}  //reclustering jet1 using the jetfinderobject Reclusterer
  return Reclusterer;
}




//----------------------------------------------------------------------
Double_t AliAnalysisTaskEmcalJetShapesMC::GetSubjetFraction(AliEmcalJet *Jet, Int_t JetContNb, Double_t JetRadius,  AliEmcalJetFinder *Reclusterer){
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
Float_t AliAnalysisTaskEmcalJetShapesMC::CoreFrac(AliEmcalJet *jet, Int_t jetContNb){
  
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
Float_t AliAnalysisTaskEmcalJetShapesMC::GetJetCoreFrac(AliEmcalJet *jet, Int_t jetContNb) {
  //calc subtracted jet mass
  
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedLeSub();
    else return jet->GetShapeProperties()->GetSecondOrderSubtractedLeSub();
    else
      return CoreFrac(jet, jetContNb);
  
}




//----------------------------------------------------------------------
Double_t AliAnalysisTaskEmcalJetShapesMC::SubJetOrdering(AliEmcalJet *Jet, AliEmcalJetFinder *Reclusterer, Int_t N, Int_t Type, Bool_t Index){
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





//______________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetShapesMC::Sigma2(AliEmcalJet *jet, Int_t jetContNb){

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
Float_t AliAnalysisTaskEmcalJetShapesMC::GetSigma2(AliEmcalJet *jet, Int_t jetContNb){
  //calc subtracted jet mass
 
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedSigma2();
    else return jet->GetShapeProperties()->GetSecondOrderSubtractedSigma2();
  else
    return Sigma2(jet, jetContNb);
 
}


//_________________________________________________________________________
void AliAnalysisTaskEmcalJetShapesMC::NTValues(AliEmcalJet *jet, Int_t jetContNb, Float_t* nTFractions){

  AliJetContainer *jetCont = GetJetContainer(jetContNb);
  if (!jet->GetNumberOfTracks())
    return;
  
  Double_t ptJet = jet->Pt();
  
  AliVParticle *vp1 = 0x0;
  std::vector<int> ordindex;
  ordindex=jet->GetPtSortedTrackConstituentIndexes(jetCont->GetParticleContainer()->GetArray());
  //Printf("Nbof const = %d", jet->GetNumberOfTracks());
  //Printf("ordindex[0] = %d, ordindex[1] = %d", ordindex[0], ordindex[1]);
  //if(ordindex.size()<2) return -1;
  
  for(Int_t iconst =0; iconst<jet->GetNumberOfTracks(); iconst++){
    
    vp1 = static_cast<AliVParticle*>(jet->TrackAt(ordindex[iconst], jetCont->GetParticleContainer()->GetArray()));
    if (!vp1){
      Printf("AliVParticle associated to Leading constituent not found");
      return;
    }

    if (nTFractions[0] <= 0.7*ptJet){
      nTFractions[0] += vp1->Pt();
      nTFractions[1] +=1;
    }
    
    if (nTFractions[2] <= 0.8*ptJet){
      nTFractions[2] += vp1->Pt();
      nTFractions[3] +=1;
    }
    
    if (nTFractions[4] <= 0.9*ptJet){
      nTFractions[4] += vp1->Pt();
      nTFractions[5] +=1;
    }
    
    if (nTFractions[6] <= 0.95*ptJet){
      nTFractions[6] += vp1->Pt();
      nTFractions[7] +=1;
    }
  }
}
//_________________________________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetShapesMC::fjNSubJettiness(AliEmcalJet *Jet, Int_t JetContNb, Int_t N, Int_t Algorithm, Double_t Beta, Int_t Option){
  
  //WARNING!!! Only works for parent jets that are clustered with Anti-Kt! To change go to AliEmcalJetFinder.cxx and look at the Nsubjettiness() function
  
  //Algorithm==0 -> kt_axes;
  // Algorithm==1 -> ca_axes;
  //Algorithm==2 -> antikt_0p2_axes;
  //Algorithm==3 -> wta_kt_axes;
  //Algorithm==4 -> wta_ca_axes;
  //Algorithm==5 -> onepass_kt_axes;
  //Algorithm==6 -> onepass_ca_axes;
  //Algorithm==7 -> onepass_antikt_0p2_axes;
  //Algorithm==8 -> onepass_wta_kt_axes;
  //Algorithm==9 -> onepass_wta_ca_axes;
  //Algorithm==10 -> min_axes;
  
  
  //Option==0 returns Nsubjettiness Value
  //Option==1 && N==2 returns opening angle between two subjet axes(Delta R?)
  //Option==2 && N==2 returns Delta R
  
  if (Jet->GetNumberOfTracks()>=N){
    AliJetContainer *JetCont = GetJetContainer(JetContNb);
      AliEmcalJetFinder *JetFinder=new AliEmcalJetFinder("Nsubjettiness");
      JetFinder->SetJetMaxEta(0.9-fJetRadius);
      JetFinder->SetRadius(fJetRadius); 
      JetFinder->SetJetAlgorithm(0); //0 for anti-kt     1 for kt  //this is for the JET!!!!!!!!!! Not the SubJets
      JetFinder->SetRecombSheme(0);
      JetFinder->SetJetMinPt(Jet->Pt());
    //Double_t dVtx[3]={vert->GetX(),vert->GetY(),vert->GetZ()};
    Double_t dVtx[3]={1,1,1};
    //Printf("JetFinder->Nsubjettiness =%f", JetFinder->Nsubjettiness(Jet,JetCont,dVtx,N,Algorithm,fSubjetRadius,Beta,Option));
    return JetFinder->Nsubjettiness(Jet,JetCont,dVtx,N,Algorithm,0.2,Beta,Option,0,0,0);
    
  }
  else return -2;
}



//________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetShapesMC::SelectTrigger(Float_t minpT, Float_t maxpT){

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

//__________________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetShapesMC::RelativePhi(Double_t mphi,Double_t vphi){

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
Bool_t AliAnalysisTaskEmcalJetShapesMC::RetrieveEventObjects() {
  //
  // retrieve event objects
  //
  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  return kTRUE;
}

//_______________________________________________________________________
void AliAnalysisTaskEmcalJetShapesMC::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.

  AliInfo("Terminate");
  AliAnalysisTaskSE::Terminate();
  
  fOutput = dynamic_cast<AliEmcalList*> (GetOutputData(1));
  if (!fOutput) {
    AliError("fOutput not available");
    return;
  }

  fTreeObservableTagging = dynamic_cast<TTree*>(GetOutputData(2));
  if (!fTreeObservableTagging){
    Printf("ERROR: fTreeObservableTagging not available");
    return;
  }

}

//_________________________________________________________________________
void AliAnalysisTaskEmcalJetShapesMC::SoftDrop(AliEmcalJet *fJet,AliJetContainer *fJetCont, double zcut, double beta, int ReclusterAlgo){
 
  std::vector<fastjet::PseudoJet>        fInputVectors;
  Double_t JetInvMass=0, PseudJetInvMass=0, TrackMom = 0, TrackEnergy = 0;
  
  AliParticleContainer *fTrackCont = fJetCont->GetParticleContainer();
 
  Double_t JetEta=fJet->Eta(),JetPhi=fJet->Phi();
  Double_t FJTrackEta[9999],FJTrackPhi[9999],FJTrackPt[9999],EmcalJetTrackEta[9999],EmcalJetTrackPhi[9999],EmcalJetTrackPt[9999];
  UShort_t FJNTracks=0,EmcalJetNTracks=0;
  
  if (fTrackCont) for (Int_t i=0; i<fJet->GetNumberOfTracks(); i++) {
      AliVParticle *fTrk = fJet->TrackAt(i, fTrackCont->GetArray());
      if (!fTrk) continue; 
      JetInvMass += fTrk->M();
      
      fastjet::PseudoJet PseudoTracks(fTrk->Px(), fTrk->Py(), fTrk->Pz(),fTrk->E());
      TrackMom += TMath::Sqrt(TMath::Power(fTrk->Px(),2)+TMath::Power(fTrk->Py(),2)+TMath::Power(fTrk->Pz(),2));
      TrackEnergy += fTrk->E();
      PseudoTracks.set_user_index(fJet->TrackAt(i)+100);
      PseudJetInvMass += PseudoTracks.m();
      fInputVectors.push_back(PseudoTracks);
      EmcalJetTrackEta[i]=fTrk->Eta();
      EmcalJetTrackPhi[i]=fTrk->Phi();
      EmcalJetTrackPt[i]=fTrk->Pt();
      EmcalJetNTracks++;
    }
  fastjet::JetDefinition                *fJetDef;         
  fastjet::ClusterSequence              *fClustSeqSA;
  


  fJetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, fJetRadius*2, static_cast<fastjet::RecombinationScheme>(0), fastjet::BestFJ30 ); 

  try {
    fClustSeqSA = new fastjet::ClusterSequence(fInputVectors, *fJetDef);
  } catch (fastjet::Error) {
    AliError(" [w] FJ Exception caught.");
    //return -1;
  }

  std::vector<fastjet::PseudoJet>       fOutputJets;
  fOutputJets.clear();
  fOutputJets=fClustSeqSA->inclusive_jets(0);
 

  std::vector<fastjet::PseudoJet> jet_constituents = fOutputJets[0].constituents();
   fastjet::contrib::SoftDrop softdrop(beta, zcut);
  //fastjet::contrib::SoftDrop softdrop_antikt(beta,zcut);
  softdrop.set_verbose_structure(kTRUE);
  //fastjet::JetDefinition jet_def_akt(fastjet::antikt_algorithm, 0.4);
  // fastjet::contrib::Recluster *antiKT_Recluster(jet_def_akt);
  fastjet::contrib::Recluster *recluster;
  if(ReclusterAlgo == 1) recluster = new fastjet::contrib::Recluster(fastjet::kt_algorithm,1,true);
  if(ReclusterAlgo == 2) recluster = new fastjet::contrib::Recluster(fastjet::antikt_algorithm,1,true);
  if(ReclusterAlgo == 0) recluster = new fastjet::contrib::Recluster(fastjet::cambridge_algorithm,1,true);  
  softdrop.set_reclustering(true,recluster);
  fastjet::PseudoJet finaljet = softdrop(fOutputJets[0]);
  // fastjet::PseudoJet finaljet_antikt = softdrop_antikt(fOutputJets[0]);
  //cout<< finaljet.structure_of<fastjet::contrib::SoftDrop>().symmetry()<<endl;
  //cout<< finaljet_antikt.structure_of<fastjet::contrib::SoftDrop>().symmetry()<<endl;


  //AliEmcalJet* jet = new AliEmcalJet(finaljet.perp(), finaljet.eta(), finaljet.phi(), finaljet.m());
  //std::vector<fastjet::PseudoJet> fSDTracks=finaljet.constituents();
  //Double_t FastjetTrackDelR,EmcalTrackDelR;
  //for(Int_t i=0;i<fJet->GetNumberOfConstituents();i++){
  //  if(i<=finaljet.constituents().size()){
  //    FastjetTrackDelR = TMath::Sqrt(TMath::Power(fSDTracks[i].eta()-JetEta,2)+TMath::Power(fSDTracks[i].phi()-JetPhi,2));
  //    FJTrackEta[i]=fSDTracks[i].eta();
  //    FJTrackPhi[i]=fSDTracks[i].phi();
  //    FJTrackPt[i]=fSDTracks[i].perp();
  //    FJNTracks++;
  //  }
  // AliVParticle *fTrk = fJet->TrackAt(i, fTrackCont->GetArray());
  // EmcalTrackDelR = TMath::Sqrt(TMath::Power(fTrk->Eta()-JetEta,2)+TMath::Power(fTrk->Phi()-JetPhi,2));       
  //}
  Int_t NDroppedTracks = fJet->GetNumberOfTracks()-finaljet.constituents().size();
  //Int_t nConstituents(fClustSeqSA->constituents(finaljet).size());
  //jet->SetNumberOfTracks(nConstituents);
  Double_t SymParam, Mu, DeltaR, GroomedPt;
  Int_t NGroomedBranches;
  SymParam=(finaljet.structure_of<fastjet::contrib::SoftDrop>().symmetry());
  Mu=(finaljet.structure_of<fastjet::contrib::SoftDrop>().mu());
  DeltaR=(finaljet.structure_of<fastjet::contrib::SoftDrop>().delta_R());
  NGroomedBranches=finaljet.structure_of<fastjet::contrib::SoftDrop>().dropped_count();
  GroomedPt=finaljet.perp();
  if(beta==0){
  if(ReclusterAlgo==0){
  fShapesVar[12]=SymParam;
  fShapesVar[13]=DeltaR;
  fShapesVar[14]=GroomedPt;
  fShapesVar[15]=NGroomedBranches;}
   if(ReclusterAlgo==1){
  fShapesVar[16]=SymParam;
  fShapesVar[17]=DeltaR;
  fShapesVar[18]=GroomedPt;
  fShapesVar[19]=NGroomedBranches;}

     if(ReclusterAlgo==2){
  fShapesVar[20]=SymParam;
  fShapesVar[21]=DeltaR;
  fShapesVar[22]=GroomedPt;
  fShapesVar[23]=NGroomedBranches;}}
  if(beta==1){
     fShapesVar[24]=SymParam;
  fShapesVar[25]=DeltaR;
  fShapesVar[26]=GroomedPt;
  fShapesVar[27]=NGroomedBranches;}

   if(beta==2){
  fShapesVar[28]=SymParam;
  fShapesVar[29]=DeltaR;
  fShapesVar[30]=GroomedPt;
  fShapesVar[31]=NGroomedBranches;}

  
  return;

  
}
