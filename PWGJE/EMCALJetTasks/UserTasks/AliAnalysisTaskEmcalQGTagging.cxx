//
// Jet QG tagging analysis task.
//
// Author: D. Caffarri, L. Cunqueiro 
//just testing
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
  fTreeObservableTagging(0)

{
   for(Int_t i=0;i<17;i++){
    fShapesVar[i]=0;}
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskEmcalQGTagging::AliAnalysisTaskEmcalQGTagging(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fContainer(0),
  fMinFractionShared(0),
  fJetShapeType(kData),
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
  fTreeObservableTagging(0)
  
{
  // Standard constructor.
  for(Int_t i=0;i<17;i++){
    fShapesVar[i]=0;}
  SetMakeGeneralHistograms(kTRUE);
  
 DefineOutput(1, TList::Class());
 DefineOutput(2, TTree::Class());
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
  
  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i=0; i<fOutput->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1*>(fOutput->At(i));
    if (h1){
      h1->Sumw2();
      continue;
    }
    THnSparse *hn = dynamic_cast<THnSparse*>(fOutput->At(i));
    if(hn)hn->Sumw2();
  }

 
  TH1::AddDirectory(oldStatus);
  const Int_t nVar = 19;
  const char* nameoutput = GetOutputSlot(2)->GetContainer()->GetName();
  fTreeObservableTagging = new TTree(nameoutput, nameoutput);
  
 
  TString *fShapesVarNames = new TString [nVar];

  fShapesVarNames[0] = "partonCode"; 
  fShapesVarNames[1] = "ptJet"; 
  fShapesVarNames[2] = "ptDJet"; 
  fShapesVarNames[3] = "mJet";
  // fShapesVarNames[4] = "nbOfConst";
  fShapesVarNames[4] = "angularity";
  fShapesVarNames[5] = "circularity";
  fShapesVarNames[6] = "lesub";
  fShapesVarNames[7] = "coronna";

  fShapesVarNames[8] = "ptJetMatch"; 
  fShapesVarNames[9] = "ptDJetMatch"; 
  fShapesVarNames[10] = "mJetMatch";
  // fShapesVarNames[12] = "nbOfConstMatch";
  fShapesVarNames[11] = "angularityMatch";
  fShapesVarNames[12] = "circularityMatch";
  fShapesVarNames[13] = "lesubMatch";
  fShapesVarNames[14] = "coronnaMatch";
  fShapesVarNames[15]="weightPythia";
  //fShapesVarNames[14]="ntrksEvt";
  fShapesVarNames[16]="rhoVal";
  fShapesVarNames[17]="rhoMassVal";
  fShapesVarNames[18]="ptUnsub";

   for(Int_t ivar=0; ivar < nVar; ivar++){
    cout<<"looping over variables"<<endl;
    fTreeObservableTagging->Branch(fShapesVarNames[ivar].Data(), &fShapesVar[ivar], Form("%s/F", fShapesVarNames[ivar].Data()));}
 
  PostData(1,fOutput);
  PostData(2,fTreeObservableTagging);

   delete [] fShapesVarNames;
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
  if (fCentSelectOn)
    if ((fCent>fCentMax) || (fCent<fCentMin)) return 0;
  
  AliAODTrack *triggerHadron = 0x0;
  
  if (fJetSelection == kRecoil) {
    //Printf("Recoil jets!!!, fminpTTrig = %f, fmaxpTTrig = %f", fminpTTrig, fmaxpTTrig);
    Int_t triggerHadronLabel = SelectTrigger(fminpTTrig, fmaxpTTrig);
     
    
    if (triggerHadronLabel==-99999) {
      //Printf ("Trigger Hadron not found, return");
      return 0;}

    AliTrackContainer *PartCont =NULL;
    AliParticleContainer *PartContMC=NULL;

    if (fJetShapeSub==kConstSub){
      if (fJetShapeType == AliAnalysisTaskEmcalQGTagging::kGenOnTheFly) PartContMC = GetParticleContainer(1);
      else PartCont = GetTrackContainer(1);
    }
    else{
      if (fJetShapeType == AliAnalysisTaskEmcalQGTagging::kGenOnTheFly) PartContMC = GetParticleContainer(0);
      else PartCont = GetTrackContainer(0);
    }
    TClonesArray *TrackArray = NULL;
    TClonesArray *TrackArrayMC = NULL;
    if (fJetShapeType == AliAnalysisTaskEmcalQGTagging::kGenOnTheFly) TrackArrayMC = PartContMC->GetArray();
    else TrackArray = PartCont->GetArray();    
    if (fJetShapeType == AliAnalysisTaskEmcalQGTagging::kGenOnTheFly) triggerHadron = static_cast<AliAODTrack*>(TrackArrayMC->At(triggerHadronLabel));
    else triggerHadron = static_cast<AliAODTrack*>(TrackArray->At(triggerHadronLabel));




    
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
  
  
  AliParticleContainer *partContAn = GetParticleContainer(0);
  TClonesArray *trackArrayAn = partContAn->GetArray();
  Int_t ntracksEvt = trackArrayAn->GetEntriesFast();
  
  Float_t rhoVal=0, rhoMassVal = 0.;
  if(jetCont) {

    jetCont->ResetCurrentID();
    if ((fJetShapeSub==kConstSub) || (fJetShapeSub==kDerivSub)){
      //rho                                                                                                   
      AliRhoParameter* rhoParam = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject("RhoSparseR020"));
      if (!rhoParam) {
	Printf("%s: Could not retrieve rho %s (some histograms will be filled with zero)!", GetName(), jetCont->GetRhoName().Data());
      } else rhoVal = rhoParam->GetVal();
      //rhom                                                                                                                                                          
      AliRhoParameter* rhomParam = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject("RhoMassSparseR020"));
      if (!rhomParam) {
	Printf("%s: Could not retrieve rho_m %s (some histograms will be filled with zero)!", GetName(), jetCont->GetRhoMassName().Data());	
      } else rhoMassVal = rhomParam->GetVal();
    }
    
    while((jet1 = jetCont->GetNextAcceptJet())) {
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
	  continue;}
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
      if(fJetShapeType == kDetEmbPartPythia){
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
        
        AliJetContainer *jetContPart=GetJetContainer(3);
        jet3=jet2->ClosestJet();
        
        if(!jet3){
          Printf("jet3 does not exist, returning");
          continue;
        }
        cout<<"jet 3 exists"<<jet3->Pt()<<endl;
        
        
        fh2ResponseUW->Fill(jet1->Pt(),jet2->Pt());
        
        Double_t fraction=0;
        if(!(fJetShapeSub==kConstSub))  fraction = jetCont->GetFractionSharedPt(jet1);
        if(fJetShapeSub==kConstSub) fraction = jetContUS->GetFractionSharedPt(jetUS);
        //if (fraction > 0.1) cout<<"***** hey a jet matched with fraction"<<fraction<<"  "<<jet1->Pt()<<" "<<jet2->Pt()<<" "<<fCent<<endl;
        
        if(fraction<fMinFractionShared) continue;
        //InputEvent()->Print();
        
      }
      
      

      if (fJetShapeType == kPythiaDef){
        
        AliJetContainer *jetContTrue = GetJetContainer(1);
        AliJetContainer *jetContUS = GetJetContainer(2);
        AliJetContainer *jetContPart = GetJetContainer(3);
        
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
        
          if (!jet2) {
            Printf("jet2 does not exist, returning");
            continue;
          }
          
          for(Int_t j=0; j<jetContPart->GetNJets(); j++) {
            
            jet3 = jetContPart->GetJet(j);
            if(!jet3) continue;
            if(jet3->GetLabel()==jet2->GetLabel()) {
              jfound++;
              if(jfound==1) jlab = j;
            }
          }
          if(jlab==-1) continue;
          jet3=jetContPart->GetJet(jlab);
          if(!jet3){
            Printf("jet3 does not exist, returning");
            continue;
          }
        }
        if(!(fJetShapeSub==kConstSub)) jet3 = jet1->ClosestJet();
        if (!jet3) {
          Printf("jet3 does not exist, returning");
          continue;
        }
        
      
        fh2ResponseUW->Fill(jet1->Pt(),jet3->Pt());
        
        
      }
      
      
      if (fJetShapeType == kGenOnTheFly){
        const AliEmcalPythiaInfo *partonsInfo = 0x0;
        partonsInfo = GetPythiaInfo();
        Double_t jp1=RelativePhi(jet1->Phi(),partonsInfo->GetPartonPhi6());
        Double_t detap1=(jet1->Eta())-(partonsInfo->GetPartonEta6());
        kWeight=partonsInfo->GetPythiaEventWeight();
        fh2ResponseW->Fill(jet1->Pt(),jet1->Pt(),kWeight);
        
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
      if (fJetShapeSub==kConstSub) ptSubtracted= jet1->Pt();
      
      else if (fJetShapeSub==kDerivSub)  {
	ptSubtracted=jet1->Pt()-GetRhoVal(0)*jet1->Area();
      }
      
      else if (fJetShapeSub==kNoSub) {
        if ((fJetShapeType==kData) || (fJetShapeType==kDetEmbPartPythia)) ptSubtracted=jet1->Pt()-GetRhoVal(0)*jet1->Area();
        else if ((fJetShapeType==kPythiaDef) || (fJetShapeType==kMCTrue) || (fJetShapeType==kGenOnTheFly)) ptSubtracted= jet1->Pt();
      }

      //Printf("ptSubtracted=%f,fPtThreshold =%f ", ptSubtracted, fPtThreshold);
      if (ptSubtracted < fPtThreshold) continue;
      
      if (fOneConstSelectOn == kTRUE) fNbOfConstvspT->Fill(GetJetNumberOfConstituents(jet1,0), ptSubtracted);
      if ((fCentSelectOn == kFALSE) && (jet1->GetNumberOfTracks() <= 1)) continue;

  
      fShapesVar[1] = ptSubtracted;
      fShapesVar[2] = GetJetpTD(jet1,0);
      fShapesVar[3] = GetJetMass(jet1,0);
      fShapesVar[4] = GetJetAngularity(jet1,0);
      fShapesVar[5] = GetJetCircularity(jet1,0);
      fShapesVar[6] = GetJetLeSub(jet1,0);
      fShapesVar[6] = GetJetCoronna(jet1,0);

      Float_t ptMatch=0., ptDMatch=0., massMatch=0., constMatch=0.,angulMatch=0.,circMatch=0., lesubMatch=0., sigma2Match=0., coronnaMatch=0;
      Int_t kMatched = 0;

       if (fJetShapeType==kPythiaDef) {
         kMatched =1;
         if(fJetShapeSub==kConstSub) kMatched = 3;
        
         ptMatch=jet3->Pt();
         ptDMatch=GetJetpTD(jet3, kMatched);
         massMatch=GetJetMass(jet3,kMatched);
         //constMatch=1.*GetJetNumberOfConstituents(jet2,kMatched);
         angulMatch=GetJetAngularity(jet3, kMatched);
        circMatch=GetJetCircularity(jet3, kMatched);
         lesubMatch=GetJetLeSub(jet3, kMatched);
	 coronnaMatch=GetJetCoronna(jet3,kMatched); 
         //sigma2Match = GetSigma2(jet2, kMatched);
       }
      
        if (fJetShapeType==kDetEmbPartPythia) {
        if(fJetShapeSub==kConstSub) kMatched = 3;
        if(fJetShapeSub==kDerivSub) kMatched = 2;
        ptMatch=jet3->Pt();
        ptDMatch=GetJetpTD(jet3, kMatched);
        massMatch=GetJetMass(jet3,kMatched);
        // constMatch=1.*GetJetNumberOfConstituents(jet3,kMatched);
        angulMatch=GetJetAngularity(jet3, kMatched);
        circMatch=GetJetCircularity(jet3, kMatched);
        lesubMatch=GetJetLeSub(jet3, kMatched);
        coronnaMatch = GetJetCoronna(jet3, kMatched);
        
      }


       
      if (fJetShapeType == kMCTrue || fJetShapeType == kData || fJetShapeType == kGenOnTheFly) {
        kMatched = 0;
        ptMatch=0.;
        ptDMatch=0.;
        massMatch=0.;
	//constMatch=0.;
        angulMatch=0.;
        circMatch=0.;
        lesubMatch=0.;
        coronnaMatch =0.;
        
      }
      
    

      fShapesVar[8] = ptMatch;
      fShapesVar[9] = ptDMatch;
      fShapesVar[10] = massMatch;
      fShapesVar[11] = angulMatch;
      fShapesVar[12] = circMatch;
      fShapesVar[13] = lesubMatch;
       fShapesVar[14] = coronnaMatch;
      fShapesVar[15] = kWeight;
      //fShapesVar[16] = ntracksEvt;
      fShapesVar[16] = rhoVal;
      fShapesVar[17] = rhoMassVal;
      fShapesVar[18] = jet1->Pt();


      fTreeObservableTagging->Fill();
      





    }
    
  }
  
  return kTRUE;
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalQGTagging::GetJetMass(AliEmcalJet *jet,Int_t jetContNb=0){
  //calc subtracted jet mass
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtracted();
      else return jet->GetShapeProperties()->GetSecondOrderSubtracted();
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
Float_t AliAnalysisTaskEmcalQGTagging::GetJetAngularity(AliEmcalJet *jet, Int_t jetContNb = 0){

  if((fJetShapeSub==kDerivSub) && (jetContNb==0))
    if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedAngularity();
      else return jet->GetShapeProperties()->GetSecondOrderSubtractedAngularity();
  else
    return Angularity(jet, jetContNb);
 
}

//____________________________________________________________________________

Float_t AliAnalysisTaskEmcalQGTagging::Coronna(AliEmcalJet *jet, Int_t jetContNb = 0){

  AliTrackContainer *PartCont = NULL;
  AliParticleContainer *PartContMC = NULL;


 if (fJetShapeSub==kConstSub){
    if (fJetShapeType == AliAnalysisTaskEmcalQGTagging::kGenOnTheFly) PartContMC = GetParticleContainer(1);
    else PartCont = GetTrackContainer(1);
  }
  else{
    if (fJetShapeType == AliAnalysisTaskEmcalQGTagging::kGenOnTheFly) PartContMC = GetParticleContainer(0);
    else PartCont = GetTrackContainer(0);
  }

  TClonesArray *TracksArray = NULL;
  TClonesArray *TracksArrayMC = NULL;

  if (fJetShapeType == AliAnalysisTaskEmcalQGTagging::kGenOnTheFly) TracksArrayMC = PartContMC->GetArray();
  else TracksArray = PartCont->GetArray();
 
  if (fJetShapeType == AliAnalysisTaskEmcalQGTagging::kGenOnTheFly){
    if(!PartContMC || !TracksArrayMC) return -2;
  }
  else {
    if(!PartCont || !TracksArray) return -2;
  }


  AliAODTrack *Track = 0x0;
  Float_t sumpt=0;
  Int_t NTracks=0;
  if (fJetShapeType == AliAnalysisTaskEmcalQGTagging::kGenOnTheFly) NTracks = TracksArrayMC->GetEntriesFast();
  else NTracks = TracksArray->GetEntriesFast();

  for(Int_t i=0; i < NTracks; i++){
    if (fJetShapeType == AliAnalysisTaskEmcalQGTagging::kGenOnTheFly){
      if((Track = static_cast<AliAODTrack*>(PartContMC->GetAcceptParticle(i)))){
	if (!Track) continue;
	if(TMath::Abs(Track->Eta())>0.9) continue;
	  Double_t dphi = RelativePhi(Track->Phi(),jet->Phi());
          Double_t dr2 = (Track->Eta()-jet->Eta())*(Track->Eta()-jet->Eta()) + dphi*dphi;
          Double_t dr = TMath::Sqrt(dr2);
	  if((dr>=0.8) && (dr<1)) sumpt=sumpt+Track->Pt();
      }
    }
    else{ 
      if((Track = static_cast<AliAODTrack*>(PartCont->GetAcceptTrack(i)))){
	if (!Track) continue;
	if(TMath::Abs(Track->Eta())>0.9) continue;
	if (Track->Pt()<0.15) continue;
	  Double_t dphi = RelativePhi(Track->Phi(),jet->Phi());
          Double_t dr2 = (Track->Eta()-jet->Eta())*(Track->Eta()-jet->Eta()) + dphi*dphi;
          Double_t dr = TMath::Sqrt(dr2);
	  if((dr>=0.8) && (dr<1)) sumpt=sumpt+Track->Pt();

      }
    } 
  }
 

  
  return sumpt; 
 


  
} 

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalQGTagging::GetJetCoronna(AliEmcalJet *jet, Int_t jetContNb = 0){

  if((fJetShapeSub==kDerivSub) && (jetContNb==0)) return -2;
  else
    return Coronna(jet, jetContNb);
 
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
Float_t AliAnalysisTaskEmcalQGTagging::GetJetpTD(AliEmcalJet *jet, Int_t jetContNb = 0){
  //calc subtracted jet mass
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedpTD();
      else return jet->GetShapeProperties()->GetSecondOrderSubtractedpTD();
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
Float_t AliAnalysisTaskEmcalQGTagging::GetJetCircularity(AliEmcalJet *jet, Int_t jetContNb =0 ){
  //calc subtracted jet mass
 
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedCircularity();
      else return jet->GetShapeProperties()->GetSecondOrderSubtractedCircularity();
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
Float_t AliAnalysisTaskEmcalQGTagging::GetJetLeSub(AliEmcalJet *jet, Int_t jetContNb =0) {
  //calc subtracted jet mass
 
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedLeSub();
      else return jet->GetShapeProperties()->GetSecondOrderSubtractedLeSub();
  else
    return LeSub(jet, jetContNb);
 
}

//________________________________________________________________________
 Float_t AliAnalysisTaskEmcalQGTagging::GetJetNumberOfConstituents(AliEmcalJet *jet,Int_t jetContNb=0){
  //calc subtracted jet mass
  
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedConstituent();
      else return jet->GetShapeProperties()->GetSecondOrderSubtractedConstituent();
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
Float_t AliAnalysisTaskEmcalQGTagging::GetSigma2(AliEmcalJet *jet, Int_t jetContNb=0){
  //calc subtracted jet mass
 
  if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
    if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedSigma2();
    else return jet->GetShapeProperties()->GetSecondOrderSubtractedSigma2();
  else
    return Sigma2(jet, jetContNb);
 
}

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalQGTagging::SelectTrigger(Float_t minpT, Float_t maxpT){

  AliTrackContainer *PartCont = NULL;
  AliParticleContainer *PartContMC = NULL;


 if (fJetShapeSub==kConstSub){
    if (fJetShapeType == AliAnalysisTaskEmcalQGTagging::kGenOnTheFly) PartContMC = GetParticleContainer(1);
    else PartCont = GetTrackContainer(1);
  }
  else{
    if (fJetShapeType == AliAnalysisTaskEmcalQGTagging::kGenOnTheFly) PartContMC = GetParticleContainer(0);
    else PartCont = GetTrackContainer(0);
  }

  TClonesArray *TracksArray = NULL;
  TClonesArray *TracksArrayMC = NULL;

  if (fJetShapeType == AliAnalysisTaskEmcalQGTagging::kGenOnTheFly) TracksArrayMC = PartContMC->GetArray();
  else TracksArray = PartCont->GetArray();
 
  if (fJetShapeType == AliAnalysisTaskEmcalQGTagging::kGenOnTheFly){
    if(!PartContMC || !TracksArrayMC) return -99999;
  }
  else {
    if(!PartCont || !TracksArray) return -99999;
  }


  AliAODTrack *Track = 0x0;

 
  
  TList *trackList = new TList();
  Int_t triggers[100];
  for (Int_t iTrigger=0; iTrigger<100; iTrigger++) triggers[iTrigger] = 0;
  Int_t iTT = 0;
  Int_t NTracks=0;
  if (fJetShapeType == AliAnalysisTaskEmcalQGTagging::kGenOnTheFly) NTracks = TracksArrayMC->GetEntriesFast();
  else NTracks = TracksArray->GetEntriesFast();

  for(Int_t i=0; i < NTracks; i++){
    if (fJetShapeType == AliAnalysisTaskEmcalQGTagging::kGenOnTheFly){
      if((Track = static_cast<AliAODTrack*>(PartContMC->GetAcceptParticle(i)))){
	if (!Track) continue;
	if(TMath::Abs(Track->Eta())>0.9) continue;
	if (Track->Pt()<0.15) continue;
	if ((Track->Pt() >= minpT) && (Track->Pt()< maxpT)) {
	  triggers[iTT] = i;
	  iTT++;
	}
      }
    }
    else{ 
      if((Track = static_cast<AliAODTrack*>(PartCont->GetAcceptTrack(i)))){
	if (!Track) continue;
	if(TMath::Abs(Track->Eta())>0.9) continue;
	if (Track->Pt()<0.15) continue;
	if ((Track->Pt() >= minpT) && (Track->Pt()< maxpT)) {
	  triggers[iTT] = i;
	  iTT++;
	}
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

