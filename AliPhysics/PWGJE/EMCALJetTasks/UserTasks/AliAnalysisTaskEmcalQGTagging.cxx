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
  fTrackCheckPlots(kFALSE),
  fCheckResolution(kFALSE),
  fSubjetCutoff(0.1),
  fMinPtConst(1),
  fHardCutoff(0),
  fDoTwoTrack(kFALSE),
  fDoAreaIterative(kTRUE),
  fPhiCutValue(0.02),
  fEtaCutValue(0.02),
  fMagFieldPolarity(1),
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
  fhTrackPhi(0x0),
  fHLundIterative(0x0),
  fHCheckResolutionSubjets(0x0),
  fNbOfConstvspT(0x0),
  fTreeObservableTagging(0)

{
   for(Int_t i=0;i<12;i++){
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
  fTrackCheckPlots(kFALSE),
  fCheckResolution(kFALSE),
  fSubjetCutoff(0.1),
  fMinPtConst(1),
  fHardCutoff(0),
  fDoTwoTrack(kFALSE),
  fDoAreaIterative(kTRUE),
  fPhiCutValue(0.02),
  fEtaCutValue(0.02),
  fMagFieldPolarity(1),
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
  fhTrackPhi(0x0),
  fHLundIterative(0x0),
  fHCheckResolutionSubjets(0x0),  
  fNbOfConstvspT(0x0),
  fTreeObservableTagging(0)
  
{
  // Standard constructor.
  for(Int_t i=0;i<12;i++){
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
  fhTrackPhi= new TH1F("fhTrackPhi", "fhTrackPhi", 100, 0, 2*TMath::Pi());
  fOutput->Add(fhTrackPhi);


  
   //log(1/theta),log(kt),jetpT,depth, algo, Eradiator// 
   const Int_t dimSpec   = 6;
   const Int_t nBinsSpec[6]     = {50,100,100,20,100,2};
   const Double_t lowBinSpec[6] = {0.,-10,0,0,0,0};
   const Double_t hiBinSpec[6]  = {5.,10.,100,20,100,2};
   fHLundIterative = new THnSparseF("fHLundIterative",
                   "LundIterativePlot [log(1/theta),log(z*theta),pTjet,algo]",
                   dimSpec,nBinsSpec,lowBinSpec,hiBinSpec);
  fOutput->Add(fHLundIterative);


  //// 
   const Int_t dimResol   = 5;
   const Int_t nBinsResol[5]     = {10,10,80,80,80};
   const Double_t lowBinResol[5] = {0,0,-1,-1,-1};
   const Double_t hiBinResol[5]  = {200,0.3,1,1,1};
   fHCheckResolutionSubjets = new THnSparseF("fHCheckResolutionSubjets",
                   "Mom.Resolution of Subjets vs opening angle",
					     dimResol,nBinsResol,lowBinResol,hiBinResol);
  fOutput->Add(fHCheckResolutionSubjets);


  
  
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
  const Int_t nVar = 12;
  const char* nameoutput = GetOutputSlot(2)->GetContainer()->GetName();
  fTreeObservableTagging = new TTree(nameoutput, nameoutput);
  
 
  TString *fShapesVarNames = new TString [nVar];

  fShapesVarNames[0] = "partonCode"; 
  fShapesVarNames[1] = "ptJet"; 
  fShapesVarNames[2] = "ktAv"; 
  fShapesVarNames[3] = "thetaAv";
  fShapesVarNames[4] = "zg";
  fShapesVarNames[5] = "deltaR";
  fShapesVarNames[6] = "ptJetMatch"; 
  fShapesVarNames[7] = "ktAvMatch"; 
  fShapesVarNames[8] = "thetaAvMatch";
   fShapesVarNames[9] = "zgMatch"; 
  fShapesVarNames[10] = "deltaRMatch";
  fShapesVarNames[11]="weightPythia";


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
    
    Int_t NTracks=0;
  if (fJetShapeType == AliAnalysisTaskEmcalQGTagging::kGenOnTheFly) NTracks = TrackArrayMC->GetEntriesFast();
  else NTracks = TrackArray->GetEntriesFast(); 




  if (fJetSelection == kRecoil) {
    //Printf("Recoil jets!!!, fminpTTrig = %f, fmaxpTTrig = %f", fminpTTrig, fmaxpTTrig);
    Int_t triggerHadronLabel = SelectTrigger(fminpTTrig, fmaxpTTrig);

   
    if (triggerHadronLabel==-99999) {
      //Printf ("Trigger Hadron not found, return");
      return 0;}

  if (fJetShapeType == AliAnalysisTaskEmcalQGTagging::kGenOnTheFly) triggerHadron = static_cast<AliAODTrack*>(TrackArrayMC->At(triggerHadronLabel));
    else triggerHadron = static_cast<AliAODTrack*>(TrackArray->At(triggerHadronLabel));
 
    /////////    
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
  

  if(fTrackCheckPlots){
      //here check tracks//
      AliAODTrack *Track = 0x0;
     for(Int_t i=0; i < NTracks; i++){
    if (fJetShapeType == AliAnalysisTaskEmcalQGTagging::kGenOnTheFly){
      if((Track = static_cast<AliAODTrack*>(PartContMC->GetAcceptParticle(i)))){
	if (!Track) continue;
	if(TMath::Abs(Track->Eta())>0.9) continue;
	if (Track->Pt()<0.15) continue;
	fhTrackPhi->Fill(Track->Phi());
      }
    }
    else{ 
      if((Track = static_cast<AliAODTrack*>(PartCont->GetAcceptTrack(i)))){
	if (!Track) continue;
	if(TMath::Abs(Track->Eta())>0.9) continue;
	if (Track->Pt()<0.15) continue;
	fhTrackPhi->Fill(Track->Phi());
      }
    } 
  }
  }
  
  
  // AliParticleContainer *partContAn = GetParticleContainer(0);
  //TClonesArray *trackArrayAn = partContAn->GetArray();
 
  
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
	// AliJetContainer *jetContTrue = GetJetContainer(1);
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
        
        //AliJetContainer *jetContPart=GetJetContainer(3);
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
        if(fCheckResolution) CheckSubjetResolution(jet1,jetCont,jet3,jetContTrue);
        
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
      RecursiveParents(jet1,jetCont);
     
      
      Float_t ptMatch=0.;
      Int_t kMatched = 0;
      Double_t ktAvMatch=0;;
      Double_t thetaAvMatch=0;
      Double_t zgMatch=0;
      Double_t deltaRMatch=0;
      Double_t aver1=0;
      Double_t aver2=0;
      Double_t aver3=0;
      Double_t aver4=0;
       if (fJetShapeType==kPythiaDef) {
         kMatched =1;
         if(fJetShapeSub==kConstSub) kMatched = 3;
        
         ptMatch=jet3->Pt();
	 RecursiveParentsMCAverage(jet3,kMatched, aver1, aver2,aver3,aver4);
	 ktAvMatch=aver1;
	 thetaAvMatch=aver2;
	 zgMatch=aver3;
	 deltaRMatch=aver4;
      
       }
      
        if (fJetShapeType==kDetEmbPartPythia) {
        if(fJetShapeSub==kConstSub) kMatched = 3;
        if(fJetShapeSub==kDerivSub) kMatched = 2;
        ptMatch=jet3->Pt();

        RecursiveParentsMCAverage(jet3,kMatched, aver1, aver2,aver3,aver4);
	 ktAvMatch=aver1;
	 thetaAvMatch=aver2;
	 zgMatch=aver3;
	 deltaRMatch=aver4;
        
      }


       
      if (fJetShapeType == kMCTrue || fJetShapeType == kData || fJetShapeType == kGenOnTheFly) {
        kMatched = 0;
        ptMatch=0.;
        ktAvMatch=0.;
        thetaAvMatch=0.;
	zgMatch=0;
	deltaRMatch=0;
	
        
      }
      
    

      fShapesVar[6] = ptMatch;
      fShapesVar[7] = ktAvMatch;
      fShapesVar[8] = thetaAvMatch;
      fShapesVar[9] = zgMatch;
      fShapesVar[10]=deltaRMatch;
      fShapesVar[11] = kWeight;
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

 
  
  TList trackList;
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
  TRandom3 random(0); 
  nbRn = random.Integer(iTT);
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


//_________________________________________________________________________
void AliAnalysisTaskEmcalQGTagging::RecursiveParents(AliEmcalJet *fJet,AliJetContainer *fJetCont){
 
  std::vector<fastjet::PseudoJet>  fInputVectors;
  fInputVectors.clear();
  fastjet::PseudoJet  PseudoTracks;
  
  AliParticleContainer *fTrackCont = fJetCont->GetParticleContainer();
  
    if (fTrackCont) for (Int_t i=0; i<fJet->GetNumberOfTracks(); i++) {
      AliVParticle *fTrk = fJet->TrackAt(i, fTrackCont->GetArray());
      if (!fTrk) continue;
      if(fDoTwoTrack==kTRUE && CheckClosePartner(i,fJet,fTrk,fTrackCont)) continue;
      PseudoTracks.reset(fTrk->Px(), fTrk->Py(), fTrk->Pz(),fTrk->E());
      PseudoTracks.set_user_index(fJet->TrackAt(i)+100);
      fInputVectors.push_back(PseudoTracks);
     
    }
    fastjet::JetAlgorithm jetalgo(fastjet::cambridge_algorithm);
     fastjet::GhostedAreaSpec ghost_spec(1, 1, 0.05);
   
  
  fastjet::JetDefinition fJetDef(jetalgo, 1., static_cast<fastjet::RecombinationScheme>(0), fastjet::BestFJ30 ); 
  fastjet::AreaDefinition fAreaDef(fastjet::passive_area,ghost_spec); 
  try {
    fastjet::ClusterSequenceArea fClustSeqSA(fInputVectors, fJetDef, fAreaDef);
    std::vector<fastjet::PseudoJet>   fOutputJets;
    fOutputJets.clear();
    fOutputJets=fClustSeqSA.inclusive_jets(0);
  
   fastjet::PseudoJet jj;
   fastjet::PseudoJet j1;
   fastjet::PseudoJet j2;
   jj=fOutputJets[0];
   double ktaverage=-20;
   double thetaverage=-20;
   double nall=0;
   double flagSubjet=0;
   double z=0.;
   double delta_R=-20;
   double area1=0;
   double area2=0;
    while(jj.has_parents(j1,j2) && z<fHardCutoff){
      nall=nall+1;
  
    flagSubjet=0;
      area1 = j1.area();
      area2 = j2.area();
      if(fJetShapeSub!=kNoSub) if(j1.perp() < j2.perp()) swap(j1,j2);
      
      if(fJetShapeSub==kNoSub && fDoAreaIterative== kTRUE){
	
	if((j1.perp()-area1*GetRhoVal(0)) < (j2.perp()-area2*GetRhoVal(0))) swap(j1,j2);
        area1 = j1.area();
        area2 = j2.area();
        if(j1.perp()-area1*GetRhoVal(0)<0 && j2.perp()-area2*GetRhoVal(0)<0) break;
       
      }

      
      if(j2.perp()-area2*GetRhoVal(0)>0){
       
    
     delta_R=j1.delta_R(j2);
   
     if(fJetShapeSub==kNoSub && fDoAreaIterative== kTRUE) {z = (j2.perp()-area2*GetRhoVal(0))/((j1.perp()-area1*GetRhoVal(0))+(j2.perp()-area2*GetRhoVal(0)));
      }
     if(fJetShapeSub!=kNoSub) z=j2.perp()/(j1.perp()+j2.perp());
  
    double y =log(1.0/delta_R);
    double lnpt_rel=log(j2.perp()*delta_R);
    double yh=j1.e()+j2.e();
    vector < fastjet::PseudoJet > constitj1 = sorted_by_pt(j1.constituents());
    if(constitj1[0].perp()>fMinPtConst) flagSubjet=1; 
    
      ktaverage=ktaverage+lnpt_rel;
      thetaverage=thetaverage+delta_R;
    Double_t LundEntries[6] = {y,lnpt_rel,fOutputJets[0].perp(),nall,yh,flagSubjet};  
    fHLundIterative->Fill(LundEntries);}
      
      jj=j1;} 

    cout<<"a z"<<z<<endl;
     fShapesVar[2]=ktaverage/nall;
     fShapesVar[3]=thetaverage/nall;
     fShapesVar[4]=z;
     fShapesVar[5]=delta_R;

  } catch (fastjet::Error) {
    AliError(" [w] FJ Exception caught.");
    //return -1;
  }




  return;

  
}
//_________________________________________________________________________
void AliAnalysisTaskEmcalQGTagging::RecursiveParentsMCAverage(AliEmcalJet *fJet,Int_t km, Double_t &average1, Double_t &average2, Double_t &average3, Double_t &average4){
  AliJetContainer *jetCont = GetJetContainer(km);
  std::vector<fastjet::PseudoJet>  fInputVectors;
  fInputVectors.clear();
  fastjet::PseudoJet  PseudoTracks;
  
  AliParticleContainer *fTrackCont = jetCont->GetParticleContainer();
  
    if (fTrackCont) for (Int_t i=0; i<fJet->GetNumberOfTracks(); i++) {
      AliVParticle *fTrk = fJet->TrackAt(i, fTrackCont->GetArray());
      if (!fTrk) continue;
     
      PseudoTracks.reset(fTrk->Px(), fTrk->Py(), fTrk->Pz(),fTrk->E());
      PseudoTracks.set_user_index(fJet->TrackAt(i)+100);
      fInputVectors.push_back(PseudoTracks);
     
    }
    fastjet::JetAlgorithm jetalgo(fastjet::cambridge_algorithm);

   
  
  fastjet::JetDefinition fJetDef(jetalgo, 1., static_cast<fastjet::RecombinationScheme>(0), fastjet::BestFJ30 ); 

  try {
    fastjet::ClusterSequence fClustSeqSA(fInputVectors, fJetDef);
    std::vector<fastjet::PseudoJet>   fOutputJets;
    fOutputJets.clear();
    fOutputJets=fClustSeqSA.inclusive_jets(0);
  
   fastjet::PseudoJet jj;
   fastjet::PseudoJet j1;
   fastjet::PseudoJet j2;
   jj=fOutputJets[0];
   double ktaverage=0;
   double thetaverage=0;
   double nall=0;
 
   double z=0;
   double delta_R=0;
    while(jj.has_parents(j1,j2) && z<fHardCutoff){
      nall=nall+1;
 

      if(j1.perp() < j2.perp()) swap(j1,j2);
    
   
    z=j2.perp()/(j1.perp()+j2.perp());
    delta_R=j1.delta_R(j2);
   
    double lnpt_rel=log(j2.perp()*delta_R);
      
      ktaverage=ktaverage+lnpt_rel;
      thetaverage=thetaverage+delta_R;
   
    jj=j1;} 
     
   
     average1=ktaverage/nall;
     average2=thetaverage/nall;
     average3=z;
     average4=delta_R;

  } catch (fastjet::Error) {
    AliError(" [w] FJ Exception caught.");
    //return -1;
  }




  return;

  
  }

//_________________________________________________________________________
void AliAnalysisTaskEmcalQGTagging::CheckSubjetResolution(AliEmcalJet *fJet,AliJetContainer *fJetCont,AliEmcalJet *fJetM,AliJetContainer *fJetContM){
 
  std::vector<fastjet::PseudoJet>  fInputVectors;
  fInputVectors.clear();
  fastjet::PseudoJet  PseudoTracks;

  std::vector<fastjet::PseudoJet>  fInputVectorsM;
  fInputVectorsM.clear();
  fastjet::PseudoJet  PseudoTracksM;
  
  AliParticleContainer *fTrackCont = fJetCont->GetParticleContainer();
  AliParticleContainer *fTrackContM = fJetContM->GetParticleContainer();
  
    if (fTrackCont) for (Int_t i=0; i<fJet->GetNumberOfTracks(); i++) {
      AliVParticle *fTrk = fJet->TrackAt(i, fTrackCont->GetArray());
      if (!fTrk) continue; 
      PseudoTracks.reset(fTrk->Px(), fTrk->Py(), fTrk->Pz(),fTrk->E());
      PseudoTracks.set_user_index(fJet->TrackAt(i)+100);
      fInputVectors.push_back(PseudoTracks);
     
    }
    fastjet::JetAlgorithm jetalgo(fastjet::cambridge_algorithm);
    fastjet::JetDefinition fJetDef(jetalgo, 1., static_cast<fastjet::RecombinationScheme>(0), fastjet::BestFJ30 ); 

   if (fTrackContM) for (Int_t i=0; i<fJetM->GetNumberOfTracks(); i++) {
      AliVParticle *fTrk = fJetM->TrackAt(i, fTrackContM->GetArray());
      if (!fTrk) continue; 
      PseudoTracksM.reset(fTrk->Px(), fTrk->Py(), fTrk->Pz(),fTrk->E());
      PseudoTracksM.set_user_index(fJetM->TrackAt(i)+100);
      fInputVectorsM.push_back(PseudoTracksM);
     
    }
    fastjet::JetAlgorithm jetalgoM(fastjet::cambridge_algorithm);
    fastjet::JetDefinition fJetDefM(jetalgoM, 1., static_cast<fastjet::RecombinationScheme>(0), fastjet::BestFJ30 ); 


    try {
    fastjet::ClusterSequence fClustSeqSA(fInputVectors, fJetDef);
    std::vector<fastjet::PseudoJet>   fOutputJets;
    fOutputJets.clear();
    fOutputJets=fClustSeqSA.inclusive_jets(0);

    fastjet::ClusterSequence fClustSeqSAM(fInputVectorsM, fJetDefM);
    std::vector<fastjet::PseudoJet>   fOutputJetsM;
    fOutputJetsM.clear();
    fOutputJetsM=fClustSeqSAM.inclusive_jets(0);

    
  
    fastjet::PseudoJet jj,jjM;
    fastjet::PseudoJet j1,j1M;
    fastjet::PseudoJet j2,j2M;
    jj=fOutputJets[0];
    jjM=fOutputJetsM[0];

   double z1=0;
   double z2=0;
   double zcut=0.1;
   while((jj.has_parents(j1,j2)) && (z1<zcut)){
    if(j1.perp() < j2.perp()) swap(j1,j2);
   
     z1=j2.perp()/(j1.perp()+j2.perp());
    jj=j1;} 
   if(z1<zcut) return;
    

     
   while((jjM.has_parents(j1M,j2M)) && (z2<zcut)){
    if(j1M.perp() < j2M.perp()) swap(j1M,j2M);
   
     z2=j2M.perp()/(j1M.perp()+j2M.perp());
    jjM=j1M;}
   if(z2<zcut) return;
        


   double delta_R1=j1.delta_R(j1M);
   double delta_R2=j2.delta_R(j2M);
   double delta_R=j1.delta_R(j2);
   double residz=(z1-z2)/z2;
   double resid1=(j1.perp()-j1M.perp())/j1M.perp(); 
   double resid2=(j2.perp()-j2M.perp())/j2M.perp(); 
    
   if((delta_R1<fSubjetCutoff) && (delta_R2<fSubjetCutoff)){
   Double_t ResolEntries[5] = {fOutputJets[0].perp(),delta_R,resid1,resid2,residz};  
   fHCheckResolutionSubjets->Fill(ResolEntries);}


   } catch (fastjet::Error) {
    AliError(" [w] FJ Exception caught.");
    //return -1;
  }




  return;





 

}


Bool_t AliAnalysisTaskEmcalQGTagging::CheckClosePartner(Int_t index, AliEmcalJet *fJet,AliVParticle *fTrk1, AliParticleContainer *fTrackCont){
      //check if tracks are close//
      for (Int_t i=0; i<fJet->GetNumberOfTracks(); i++) {
      AliVParticle *fTrk2 = fJet->TrackAt(i, fTrackCont->GetArray());
      if (!fTrk2) continue;
      if(i==index) continue;
      Double_t phi1 = fTrk1->Phi();
      Double_t phi2 = fTrk2->Phi();
      Double_t chg1 = fTrk1->Charge();
      Double_t chg2 = fTrk2->Charge();
      Double_t ptv1 = fTrk1->Pt();
      Double_t ptv2 = fTrk2->Pt();
      Double_t deta=fTrk2->Eta()-fTrk1->Eta();
      const Float_t kLimit = fPhiCutValue* 3;

      if (TMath::Abs(fTrk1->Eta()-fTrk2->Eta()) < fEtaCutValue*2.5*3){
      Float_t initdpsinner = (phi2 - TMath::ASin(0.075*chg2*fMagFieldPolarity*0.8/ptv2) -
      (phi1 - TMath::ASin(0.075*chg1*fMagFieldPolarity*0.8/ptv1)));
      
      Float_t initdpsouter = (phi2 - TMath::ASin(0.075*chg2*fMagFieldPolarity*2.5/ptv2) -
      (phi1 - TMath::ASin(0.075*chg1*fMagFieldPolarity*2.5/ptv1)));

      initdpsinner = TVector2::Phi_mpi_pi(initdpsinner);
      initdpsouter = TVector2::Phi_mpi_pi(initdpsouter);

      

     if (TMath::Abs(initdpsinner) < kLimit ||
      TMath::Abs(initdpsouter) < kLimit || initdpsinner * initdpsouter < 0 ) {
      Double_t mindps = 1e5;
    
   
   
        for (Double_t rad = 0.8; rad < 2.51; rad += 0.01) {
          Double_t dps = (phi2 - TMath::ASin(0.075*chg2*fMagFieldPolarity*rad/ptv2) - (phi1 - TMath::ASin(0.075*chg1*fMagFieldPolarity*rad/ptv1)));
          dps = TVector2::Phi_mpi_pi(dps);
          if (TMath::Abs(dps) < TMath::Abs(mindps))
            mindps = dps;
        }
	if(TMath::Abs(mindps)<fPhiCutValue && TMath::Abs(deta)<fEtaCutValue) return kTRUE;
     } }
	 return kFALSE;
      }}





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

