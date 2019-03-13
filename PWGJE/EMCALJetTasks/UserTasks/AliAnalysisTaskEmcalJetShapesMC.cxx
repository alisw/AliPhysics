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
#include "TF1.h"
#include "AliEmcalJetFinder.h"
#include "AliAODEvent.h"
#include "AliAnalysisTaskEmcalJetShapesMC.h"

#include <fastjet/config.h>
#if FASJET_VERSION_NUMBER >= 30302
#include <fastjet/tools/Recluster.hh>
#else 
#include <fastjet/contrib/Recluster.hh>
#endif

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
  fSwitchKtNSub(0),
  fSwitchMinNSub(0),
  fSwitchAktNSub(0),
  fSwitchSDKtNSub(0),
  fSwitchSDMinNSub(0),
  fAdditionalTracks(0),
  fHardCutoff(0),
  fOptionalPartonInfo(0),
  fminpTTrig(20.),
  fmaxpTTrig(50.),
  fangWindowRecoil(0.6),
  fSemigoodCorrect(0),
  fHolePos(0),
  fHoleWidth(0),
  fRandom(0),
  fqhat(1),
  fxlength(2),
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
  fHLundIterative(0x0),
  fHLundIterative_ktaxis(0x0),
  fHLundIterativeInject(0x0),
  fNbOfConstvspT(0x0),
  fTreeObservableTagging(0x0),
  fTf1SoftOmega(0x0),
  fTf1SoftKt(0x0),
  fTf1Omega(0x0),
  fTf1Kt(0x0),
  fScaleELoss(kFALSE),
  xfraction(1),
  fAddMedScat(kFALSE),
  fAddMedScatPtFrac(1),
  fAddMedScatN(100)


{
  for(Int_t i=0;i<13;i++){
    fShapesVar[i]=0;}
  SetMakeGeneralHistograms(kTRUE);
   DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
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
  fSwitchKtNSub(0),
  fSwitchMinNSub(0),
  fSwitchAktNSub(0),
  fSwitchSDKtNSub(0),
  fSwitchSDMinNSub(0),
  fAdditionalTracks(0),
  fHardCutoff(0),
  fOptionalPartonInfo(0),
  fminpTTrig(20.),
  fmaxpTTrig(50.),
  fangWindowRecoil(0.6),
  fSemigoodCorrect(0),
  fHolePos(0),
  fHoleWidth(0),
  fRandom(0),
  fqhat(1),
  fxlength(2),
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
  fHLundIterative(0x0),
  fHLundIterative_ktaxis(0x0),
  fHLundIterativeInject(0x0),
  fNbOfConstvspT(0x0),
  fTreeObservableTagging(0x0),
  fTf1SoftOmega(0x0),
  fTf1SoftKt(0x0),
  fTf1Omega(0x0),
  fTf1Kt(0x0),
  fScaleELoss(kFALSE),
  xfraction(1),
  fAddMedScat(kFALSE),
  fAddMedScatPtFrac(1),
  fAddMedScatN(100)
{
  // Standard constructor.


  for(Int_t i=0;i<13;i++){
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

   if(fRandom)      delete fRandom;
   if(fTf1SoftOmega)    delete fTf1SoftOmega;
   if(fTf1SoftKt)        delete fTf1SoftKt;
   if(fTf1Omega)    delete fTf1Omega;
   if(fTf1Kt)        delete fTf1Kt;

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

  const Int_t nVar = 13;

  TString *fShapesVarNames = new TString [nVar];

  fShapesVarNames[0] = "partonCode";
  fShapesVarNames[1] = "ptJet";
  fShapesVarNames[2] = "ptDJet";
  fShapesVarNames[3] = "mJet";
  fShapesVarNames[4] = "nbOfConst";
  fShapesVarNames[5] = "angularity";
  fShapesVarNames[6] = "nitersd";
  fShapesVarNames[7] = "niterall";
  fShapesVarNames[8] = "weightPythia";
  //fShapesVarNames[6] = "Nsubjet1kt";
  // fShapesVarNames[7] = "Nsubjet2kt";
  // fShapesVarNames[8] = "Nsubjet1Min";
  // fShapesVarNames[9] = "Nsubjet2Min";
  // fShapesVarNames[10] = "DeltaRkt";
  // fShapesVarNames[11] = "DeltaRMin";
  fShapesVarNames[9] = "SDSymm";
  fShapesVarNames[10] = "scaledptJet";
  fShapesVarNames[11] = "zg";
  fShapesVarNames[12] = "rg";
  // fShapesVarNames[15] = "SDGroomedN";
  // fShapesVarNames[16] = "SDMass";
  // fShapesVarNames[17] = "SDSymmkt";
  //  fShapesVarNames[18] = "SDDeltaRkt";
  // fShapesVarNames[19] = "SDGroomedFrackt";
  // fShapesVarNames[20] = "SDGroomedNkt";
  // fShapesVarNames[21] = "SDMasskt";
  //  fShapesVarNames[22] = "SDSymmAkt";
  //  fShapesVarNames[23] = "SDDeltaRAkt";
  // fShapesVarNames[24] = "SDGroomedFracAkt";
  // fShapesVarNames[25] = "SDGroomedNAkt";
  //  fShapesVarNames[26] = "SDMassAkt";
  // fShapesVarNames[27] = "SDSymmktForm";
  //  fShapesVarNames[28] = "SDDeltaRktForm";
  // fShapesVarNames[29] = "SDGroomedFracktForm";
  // fShapesVarNames[30] = "SDGroomedNktForm";
  // fShapesVarNames[31] = "SDMassktForm";
  // fShapesVarNames[32] = "SDSymmDemo";
  // fShapesVarNames[33] = "SDDeltaRDemo";
  // fShapesVarNames[34] = "SDGroomedFracDemo";
  // fShapesVarNames[35] = "SDGroomedNDemo";
  // fShapesVarNames[36] = "SDMassDemo";
  // fShapesVarNames[42] = "SDSymmForm";
  // fShapesVarNames[43] = "SDDeltaRForm";
  // fShapesVarNames[44] = "SDGroomedFracForm";
  // fShapesVarNames[45] = "SDGroomedNForm";
  // fShapesVarNames[46] = "SDMassForm";
  // fShapesVarNames[47] = "weightPythia";
  // fShapesVarNames[42] = "SDSymmNoCut";
  // fShapesVarNames[43] = "SDDeltaRNoCut";
  // fShapesVarNames[44] = "SDGroomedFracNoCut";
  // fShapesVarNames[45] = "SDGroomedNNoCut";
  // fShapesVarNames[46] = "SDMassNoCut";


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



    //log(1/theta),log(z*theta),jetpT,algo//
   const Int_t dimSpec   = 6;
   const Int_t nBinsSpec[6]     = {50,100,60,3,22,10};
   const Double_t lowBinSpec[6] = {0.9,-10,0.0,0.0,0.0,0.0};
   const Double_t hiBinSpec[6]  = {5.0,2.0,600.0,3.0,22.0,10.0};
   fHLundIterative = new THnSparseF("fHLundIterative",
                   "LundIterativePlot [log(1/theta),log(z*theta),pTjet,algo,partonFlavor,depth]",
                   dimSpec,nBinsSpec,lowBinSpec,hiBinSpec);
  fOutput->Add(fHLundIterative);
  const Int_t nBinsSpec_ktaxis[6]     = {50,100,60,3,22,10};
  const Double_t lowBinSpec_ktaxis[6] = {0.9,-4.0,0.0,0.0,0.0,0.0};
  const Double_t hiBinSpec_ktaxis[6]  = {5.0,8.0,600.0,3.0,22.0,10.0};
  fHLundIterative_ktaxis = new THnSparseF("fHLundIterative_ktaxis",
                  "LundIterativePlot [log(1/theta),log(z*theta),pTjet,algo,partonFlavor,depth]",
                  dimSpec,nBinsSpec_ktaxis,lowBinSpec_ktaxis,hiBinSpec_ktaxis);
  fOutput->Add(fHLundIterative_ktaxis);

  if(fAdditionalTracks>0){
  //log(1/theta),log(z*theta),jetpT of added tracks
   const Int_t dimSpecb   = 4;
   const Int_t nBinsSpecb[4]     = {50,50,20,20};
   const Double_t lowBinSpecb[4] = {0.0,-10,  0,0};
   const Double_t hiBinSpecb[4]  = {5.0,  0,200,200};
   fHLundIterativeInject = new THnSparseF("fHLundIterativeInject",
                   "LundIterativePlotInject [log(1/theta),log(z*theta),pTjet,algo]",
                   dimSpecb,nBinsSpecb,lowBinSpecb,hiBinSpecb);
   fOutput->Add(fHLundIterativeInject);}


  fNbOfConstvspT=new TH2F("fNbOfConstvspT", "fNbOfConstvspT", 100, 0, 100, 200, 0, 200);
  fOutput->Add(fNbOfConstvspT);

  //fOutput->Add(fTreeObservableTagging);

 fRandom = new TRandom3(0);
  PostData(1, fOutput); // Post data for ALL output slots > 0 here
  PostData(2, fTreeObservableTagging);

  delete [] fShapesVarNames;

   }

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetShapesMC::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().
  // if (gRandom) delete gRandom;
  //   gRandom = new TRandom3(0);
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

      fPtJet->Fill(jet1->Pt());




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

      if (fOptionalPartonInfo==1){
        const AliEmcalPythiaInfo *partonsInfo = 0x0;
        partonsInfo = GetPythiaInfo();
        //Printf("partonsInfo=%p",  partonsInfo);
        Double_t jp1=RelativePhi(jet1->Phi(),partonsInfo->GetPartonPhi6());
        Double_t detap1=(jet1->Eta())-(partonsInfo->GetPartonEta6());
        kWeight=partonsInfo->GetPythiaEventWeight();
        //Printf("kWeight=%f",  kWeight);
        fShapesVar[8] = kWeight;

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
      SoftDrop(jet1,jetCont,0.1,0,0);

      RecursiveParents(jet1,jetCont,0,fShapesVar[0]);
      RecursiveParents(jet1,jetCont,1,fShapesVar[0]);
      RecursiveParents(jet1,jetCont,2,fShapesVar[0]);


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


  //Double_t dVtx[3]={vert->GetX(),vert->GetY(),vert->GetZ()};
  Double_t dVtx[3]={0.,0.,0.};
  if(Reclusterer->AliEmcalJetFinder::Filter(Jet, JetCont, dVtx)){;}  //reclustering jet1 using the jetfinderobject Reclusterer
  return Reclusterer;
}




//----------------------------------------------------------------------
Double_t AliAnalysisTaskEmcalJetShapesMC::GetSubjetFraction(AliEmcalJet *Jet, Int_t JetContNb, Double_t JetRadius,  AliEmcalJetFinder *Reclusterer){
  AliJetContainer *JetCont = GetJetContainer(JetContNb);



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
Double_t AliAnalysisTaskEmcalJetShapesMC::FjNSubJettiness(AliEmcalJet *Jet, Int_t JetContNb, Int_t N, Int_t Algorithm, Double_t Beta, Int_t Option, Double_t Beta_SD, Double_t ZCut, Int_t SoftDropOn){

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
    return JetFinder->Nsubjettiness(Jet,JetCont,dVtx,N,Algorithm,0.2,Beta,Option,0,Beta_SD,ZCut,SoftDropOn);

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
void AliAnalysisTaskEmcalJetShapesMC::SoftDrop(AliEmcalJet *fJet,AliJetContainer *fJetCont, Double_t zcut, Double_t beta, Int_t ReclusterAlgo){

  std::vector<fastjet::PseudoJet>        fInputVectors;
  fInputVectors.clear();
  Double_t JetInvMass=0, PseudJetInvMass=0, TrackMom = 0, TrackEnergy = 0;
  fastjet::PseudoJet  PseudoTracks;
  fastjet::PseudoJet MyJet;
  fastjet::PseudoJet PseudoTracksCMS;
  AliParticleContainer *fTrackCont = fJetCont->GetParticleContainer();
  //cout<<"CALL TO SOFTDROP"<<endl;

  Double_t zeta=0;
  Double_t angle=0;
   if (fTrackCont) for (Int_t i=0; i<fJet->GetNumberOfTracks(); i++) {
      AliVParticle *fTrk = fJet->TrackAt(i, fTrackCont->GetArray());
      if (!fTrk) continue;
      JetInvMass += fTrk->M();

      PseudoTracks.reset(fTrk->Px(), fTrk->Py(), fTrk->Pz(),fTrk->E());
      TrackMom += TMath::Sqrt(TMath::Power(fTrk->Px(),2)+TMath::Power(fTrk->Py(),2)+TMath::Power(fTrk->Pz(),2));
      TrackEnergy += fTrk->E();
      PseudoTracks.set_user_index(fJet->TrackAt(i)+100);
      PseudJetInvMass += PseudoTracks.m();
      fInputVectors.push_back(PseudoTracks);

    }

   //fRandom->SetSeed(0);
  //here add N tracks with random phi and eta and theta according to bdmps distrib.

    MyJet.reset(fJet->Px(),fJet->Py(),fJet->Pz(),fJet->E());
    Double_t omegac=0.5*fqhat*fxlength*fxlength/0.2;
    Double_t thetac=TMath::Sqrt(12*0.2/(fqhat*TMath::Power(fxlength,3)));
    Double_t xQs=TMath::Sqrt(fqhat*fxlength);
    //cout<<"medium parameters "<<omegac<<" "<<thetac<<" "<<xQs<<endl;

   for(Int_t i=0;i<fAdditionalTracks;i++){

     Double_t ppx,ppy,ppz,kTscale,lim2o,lim1o;
     Double_t lim2=xQs;
    Double_t lim1=10000;

    //generation of kT according to 1/kT^4, with minimum QS=2 GeV and maximum ~sqrt(ptjet*T)
     fTf1Kt= new TF1("fTf1Kt","1/(x*x*x*x)",lim2,lim1);
     kTscale=fTf1Kt->GetRandom();
     //generation within the jet cone

     //generation of w according to 1/w, with minimum wc
     //omega needs to be larger than kT so to have well defined angles
     lim2o=kTscale;
     lim1o=kTscale/TMath::Sin(0.1);
     fTf1Omega= new TF1("fTf1Omega","1/x",lim2o,lim1o);
     Double_t omega=fTf1Omega->GetRandom();

     Double_t sinpptheta=kTscale/omega;
     Double_t pptheta=TMath::ASin(sinpptheta);
     //cout<<"angle_omega_kt"<<pptheta<<" "<<omega<<" "<<kTscale<<endl;
     if(pptheta>fJetRadius) continue;

     PseudoTracksCMS.reset(kTscale/TMath::Sqrt(2),kTscale/TMath::Sqrt(2),omega*TMath::Cos(pptheta),omega);
     //boost the particle in the rest frame of the jet to the lab frame
     fastjet::PseudoJet PseudoTracksLab=PseudoTracksCMS.boost(MyJet);
     PseudoTracksLab.set_user_index(i+fJet->GetNumberOfTracks()+100);
     fInputVectors.push_back(PseudoTracksLab);
     //in the frame of the jet
     zeta=omega/fJet->E();
     angle=pptheta;


   }









  fastjet::JetDefinition fJetDef(fastjet::antikt_algorithm, fJetRadius*2, static_cast<fastjet::RecombinationScheme>(0), fastjet::BestFJ30 );
  try {
    fastjet::ClusterSequence fClustSeqSA(fInputVectors, fJetDef);
    std::vector<fastjet::PseudoJet>       fOutputJets;
    fOutputJets.clear();
    fOutputJets=fClustSeqSA.inclusive_jets(0);

    //cout<<fOutputJets[0].perp()<<" "<<fJet->Pt()<<endl;

    fastjet::contrib::SoftDrop softdrop(beta, zcut);

    softdrop.set_verbose_structure(kTRUE);
    fastjet::JetAlgorithm jetalgo(fastjet::cambridge_algorithm);
    if(ReclusterAlgo==2) jetalgo=fastjet::antikt_algorithm;
  if(ReclusterAlgo==1) jetalgo=fastjet::kt_algorithm;
    if(ReclusterAlgo==0) jetalgo=fastjet::cambridge_algorithm;

#if FASTJET_VERSION_NUMBER >= 30302
    fastjet::Recluster recluster(jetalgo,1,fastjet::Recluster::keep_only_hardest);
#else
    fastjet::contrib::Recluster recluster(jetalgo,1,true);
#endif
    softdrop.set_reclustering(true, &recluster);
    fastjet::PseudoJet finaljet = softdrop(fOutputJets[0]);
    Int_t NDroppedTracks = fJet->GetNumberOfTracks()-finaljet.constituents().size();

    Double_t SymParam, Mu, DeltaR, GroomedPt,GroomedMass;
    Int_t NGroomedBranches;
    SymParam=(finaljet.structure_of<fastjet::contrib::SoftDrop>().symmetry());
    Mu=(finaljet.structure_of<fastjet::contrib::SoftDrop>().mu());
    DeltaR=(finaljet.structure_of<fastjet::contrib::SoftDrop>().delta_R());
    NGroomedBranches=finaljet.structure_of<fastjet::contrib::SoftDrop>().dropped_count();
    GroomedPt=finaljet.perp();
    GroomedMass=finaljet.m();


    fShapesVar[11]=SymParam;
    fShapesVar[12]=DeltaR;
  // fShapesVar[14]=zeta;
  // fShapesVar[15]=angle;
  // fShapesVar[16]=GroomedMass;}
  //  if(ReclusterAlgo==1){
  // fShapesVar[17]=SymParam;
  // fShapesVar[18]=DeltaR;
  // fShapesVar[19]=zeta;
  // fShapesVar[20]=angle;
  // fShapesVar[21]=GroomedMass; }

  //    if(ReclusterAlgo==2){
  // fShapesVar[22]=SymParam;
  // fShapesVar[23]=DeltaR;
  // fShapesVar[24]=zeta;
  // fShapesVar[25]=angle;
  // fShapesVar[26]=GroomedMass;
  //    }}
  // if(beta==1){
  //    fShapesVar[27]=SymParam;
  // fShapesVar[28]=DeltaR;
  // fShapesVar[29]=zeta;
  // fShapesVar[30]=angle;
  // fShapesVar[31]=GroomedMass;
  // }
  // //this one kills soft and large angle radiation
  // if((beta==1.5) && (zcut==0.5)){
  // fShapesVar[32]=SymParam;
  // fShapesVar[33]=DeltaR;
  // fShapesVar[34]=zeta;
  // fShapesVar[35]=angle;
  // fShapesVar[36]=GroomedMass; }
  //  //this option favour democratic branches at large kt
  //  if((beta==-1) && (zcut==0.005)){
  // fShapesVar[37]=SymParam;
  // fShapesVar[38]=DeltaR;
  // fShapesVar[39]=zeta;
  // fShapesVar[40]=angle;
  // fShapesVar[41]=GroomedMass; }

  // if((beta==-2) && (zcut==0.005)){
  // fShapesVar[42]=SymParam;
  // fShapesVar[43]=DeltaR;
  // fShapesVar[44]=zeta;
  // fShapesVar[45]=angle;
  // fShapesVar[46]=GroomedMass; }









  } catch (fastjet::Error) {
    AliError(" [w] FJ Exception caught.");
    //return -1;
  }
  if(fTf1Kt){ delete fTf1Kt;}
  if(fTf1Omega){ delete fTf1Omega;}
}


//_________________________________________________________________________
void AliAnalysisTaskEmcalJetShapesMC::RecursiveParents(AliEmcalJet *fJet,AliJetContainer *fJetCont, Int_t ReclusterAlgo, Float_t partonFlavor){

  std::vector<fastjet::PseudoJet>  fInputVectors;
  fInputVectors.clear();
  fastjet::PseudoJet  PseudoTracks;
  fastjet::PseudoJet  PseudoTracksLab;
  double xflagalgo=0;
  double lnpt_relinject=0;
  double yinject=0;
  int xflagAdded=0;
  double zinject,angleinject,pptheta,sinpptheta,omega,omega2,angle2;
  AliParticleContainer *fTrackCont = fJetCont->GetParticleContainer();
  Float_t pTscale=0., phiscale=0., thetascale=0., pXscale=0., pYscale=0., pZscale=0., pscale=0.;

  if (fTrackCont) for (Int_t i=0; i<fJet->GetNumberOfTracks(); i++) {
      AliVParticle *fTrk = fJet->TrackAt(i, fTrackCont->GetArray());
      if (!fTrk) continue;
      if (fScaleELoss){
	pTscale    = xfraction*sqrt(pow(fTrk->Px(),2)+pow(fTrk->Py(),2));
	phiscale   = fTrk->Phi();
	thetascale = 2.*TMath::ATan(TMath::Exp(-1.*(fTrk->Eta())));
	pXscale    = pTscale * TMath::Cos(phiscale);
	pYscale    = pTscale * TMath::Sin(phiscale);
	pZscale    = pTscale/TMath::Tan(thetascale);
	pscale     = TMath::Sqrt(pTscale*pTscale+pZscale*pZscale);
	PseudoTracks.reset(pXscale, pYscale, pZscale, pscale);
      }
      else PseudoTracks.reset(fTrk->Px(), fTrk->Py(), fTrk->Pz(),fTrk->E());
      PseudoTracks.set_user_index(fJet->TrackAt(i)+100);
      fInputVectors.push_back(PseudoTracks);

    }
  /*if(fAddMedScat){
    for(int i = 0; i < fAddMedScatN; i++){
      TRandom3 rand1(0),rand2(0),rand3(0); //set range +- jet R
      Double_t randN1 = 0.4*0.4*rand1.Rndm();
      Double_t randN2 = 2*TMath::Pi()*rand2.Rndm();
      Double_t phi_rand = (fJet->Phi())+TMath::Sqrt(randN1)*TMath::Sin(randN2);
      Double_t eta_rand = (fJet->Eta())+TMath::Sqrt(randN1)*TMath::Cos(randN2);
      Double_t fAddMedScatPt = (fAddMedScatPtFrac*fJet->Pt())/fAddMedScatN;
      PseudoTracks.reset(fAddMedScatPt*TMath::Cos(phi_rand),fAddMedScatPt*TMath::Sin(phi_rand),fAddMedScatPt/TMath::Tan(eta_rand),fAddMedScatPt);
      PseudoTracks.set_user_index(i+fJet->GetNumberOfTracks()+100);
      fInputVectors.push_back(PseudoTracks);
    }
    }*/
  if(fAddMedScat){
    for(int i = 0; i < fAddMedScatN; i++){
      Double_t ppx,ppy,ppz,SoftkTscale,lim2o,lim1o;
      Double_t lim1=0.1;
      Double_t lim2=0.5;
      fTf1SoftKt= new TF1("fTf1SoftKt","1/(x)",lim1,lim2);
      SoftkTscale=fTf1SoftKt->GetRandom();

      lim2o=SoftkTscale;
      lim1o=SoftkTscale/TMath::Sin(0.1);
      fTf1SoftOmega= new TF1("fTf1SoftOmega","1/x",lim2o,lim1o);
      omega=fTf1SoftOmega->GetRandom();
      sinpptheta=SoftkTscale/omega;
      pptheta=TMath::ASin(sinpptheta);
      if(pptheta>fJetRadius) continue;

      TLorentzVector pTrackCMS(SoftkTscale/TMath::Sqrt(2),SoftkTscale/TMath::Sqrt(2),omega*TMath::Cos(pptheta),omega);
      TVector3 MyJet(fJet->Px(),fJet->Py(),fJet->Pz());
      TVector3 direction = MyJet.Unit();
      //rotate the track to the jet frame
      pTrackCMS.RotateUz(direction);

      //add the rotated track to the jet
      PseudoTracksLab.reset(pTrackCMS.Px(),pTrackCMS.Py(),pTrackCMS.Pz(),pTrackCMS.E());

      PseudoTracksLab.set_user_index(i+fJet->GetNumberOfTracks()+100);

      omega2=PseudoTracksLab.perp();
      angle2=pTrackCMS.Angle(MyJet);

      fInputVectors.push_back(PseudoTracksLab);
    }
  }


    //add tracks to the jet prior to the reclusterer in case of iterative mapping of splittings

    Double_t omegac=0.5*fqhat*fxlength*fxlength/0.2;
    Double_t thetac=TMath::Sqrt(12*0.2/(fqhat*TMath::Power(fxlength,3)));
    Double_t xQs=TMath::Sqrt(fqhat*fxlength);


   for(Int_t i=0;i<fAdditionalTracks;i++){

    Double_t ppx,ppy,ppz,kTscale,lim2o,lim1o;
    Double_t lim2=xQs;
    Double_t lim1=10000;

    //generation of kT according to 1/kT^4, with minimum QS=2 GeV and maximum ~sqrt(ptjet*T)
     fTf1Kt= new TF1("fTf1Kt","1/(x*x*x*x)",lim2,lim1);
     kTscale=fTf1Kt->GetRandom();
     //generation within the jet cone

     //generation of w according to 1/w, with minimum wc
     //omega needs to be larger than kT so to have well defined angles
     lim2o=kTscale;
     lim1o=kTscale/TMath::Sin(0.1);
     fTf1Omega= new TF1("fTf1Omega","1/x",lim2o,lim1o);
     omega=fTf1Omega->GetRandom();
     sinpptheta=kTscale/omega;
     pptheta=TMath::ASin(sinpptheta);
     if(pptheta>fJetRadius) continue;

     //Lorentz vector in the frame where the jet moves along z axis
     TLorentzVector pTrackCMS(kTscale/TMath::Sqrt(2),kTscale/TMath::Sqrt(2),omega*TMath::Cos(pptheta),omega);
     TVector3 MyJet(fJet->Px(),fJet->Py(),fJet->Pz());
     TVector3 direction = MyJet.Unit();
     //rotate the track to the jet frame
     pTrackCMS.RotateUz(direction);

     //add the rotated track to the jet
     PseudoTracksLab.reset(pTrackCMS.Px(),pTrackCMS.Py(),pTrackCMS.Pz(),pTrackCMS.E());

     PseudoTracksLab.set_user_index(i+fJet->GetNumberOfTracks()+100);

     omega2=PseudoTracksLab.perp();
     angle2=pTrackCMS.Angle(MyJet);

     fInputVectors.push_back(PseudoTracksLab);
     xflagAdded=1;
   }




    fastjet::JetAlgorithm jetalgo(fastjet::antikt_algorithm);
    if(ReclusterAlgo==0){ xflagalgo=0.5;
      jetalgo=fastjet::kt_algorithm ;}

      if(ReclusterAlgo==1){ xflagalgo=1.5;
	jetalgo=fastjet::cambridge_algorithm;}
	if(ReclusterAlgo==2){ xflagalgo=2.5;
	  jetalgo=fastjet::antikt_algorithm;}

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
   double ndepth=0;
   double nsd=0;
   double nall=0;
   fShapesVar[10]=fOutputJets[0].perp();
    while(jj.has_parents(j1,j2)){

    if(j1.perp() < j2.perp()) swap(j1,j2);
    double delta_R=j1.delta_R(j2);
    double z=j2.perp()/(j1.perp()+j2.perp());
    double y =log(1.0/delta_R);
    double lnpt_rel=log(z*delta_R);
    double lnkt=log(jj.perp()*z*delta_R);

    nall=nall+1;
    if(z>fHardCutoff){
    ndepth=ndepth+1;
    nsd=nsd+1;
    if(nsd==1 && ReclusterAlgo==1) fShapesVar[9]=z;
    Double_t LundEntries[6] = {y,lnpt_rel,fOutputJets[0].perp(),xflagalgo,partonFlavor,ndepth};
    Double_t LundEntries_kt[6] = {y,lnkt,fOutputJets[0].perp(),xflagalgo,partonFlavor,ndepth};
    fHLundIterative->Fill(LundEntries);
    fHLundIterative_ktaxis->Fill(LundEntries_kt);
    }
    jj=j1;}

    if(ReclusterAlgo==1){
    fShapesVar[6]=nsd;
    fShapesVar[7]=nall;}

    if(fAdditionalTracks>0 && xflagAdded>0){
     zinject=omega2/fOutputJets[0].perp();
     angleinject=angle2;

     yinject =log(1.0/angleinject);
     lnpt_relinject=log(zinject*angleinject);
     Double_t LundEntriesInject[4] = {yinject,lnpt_relinject,fOutputJets[0].perp(),fJet->Pt()};
     fHLundIterativeInject->Fill(LundEntriesInject);}


  } catch (fastjet::Error) {
    AliError(" [w] FJ Exception caught.");
    //return -1;
  }


  if(fTf1SoftKt){ delete fTf1SoftKt;}
  if(fTf1SoftOmega){ delete fTf1SoftOmega;}
  if(fTf1Kt){ delete fTf1Kt;}
  if(fTf1Omega){ delete fTf1Omega;}


  return;


}


AliAnalysisTaskEmcalJetShapesMC* AliAnalysisTaskEmcalJetShapesMC::AddTaskJetShapesMC(const char * njetsBase,
                                                    const Double_t jetradius,
                                                    const Double_t subjetradius,
                                                    const char *ntracksPartLevel,
                                                    const char *type,
                                                    const char *CentEst,
                                                    Int_t       pSel,
                                                    TString     trigClass,
                                                    TString     kEmcalTriggers,
                                                    TString     tag,
                                                    const char *rhoName,
                                                    AliAnalysisTaskEmcalJetShapesMC::JetShapeType jetShapeType,
                                                    AliAnalysisTaskEmcalJetShapesMC::JetShapeSub jetShapeSub,
                                                    AliAnalysisTaskEmcalJetShapesMC::JetSelectionType jetSelection,
                                                    Float_t minpTHTrigger,  Float_t maxpTHTrigger,
                                                    AliAnalysisTaskEmcalJetShapesMC::DerivSubtrOrder derivSubtrOrder) {



  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    {
      ::Error("AliAnalysisTaskEmcalJetShapesMC","No analysis manager found.");
      return 0;
    }
  Bool_t ismc=kFALSE;
  ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
    {
      ::Error("AddTaskJetShapesMC", "This task requires an input event handler");
      return NULL;
    }

  TString wagonName1 = Form("JetShapesMC_%s_Histos%s%s",njetsBase,trigClass.Data(),tag.Data());
  TString wagonName2 = Form("JetShapesMC_%s_Tree%s%s",njetsBase,trigClass.Data(),tag.Data());

  //Configure jet tagger task
  AliAnalysisTaskEmcalJetShapesMC *task = new AliAnalysisTaskEmcalJetShapesMC(wagonName1.Data());

  //task->SetNCentBins(4);
  task->SetJetShapeType(jetShapeType);
  task->SetJetShapeSub(jetShapeSub);
  task->SetJetSelection(jetSelection);
  task->SetDerivativeSubtractionOrder(derivSubtrOrder);
  task->SetJetRadius(jetradius);
  task->SetSubjetRadius(subjetradius);

  if (jetSelection == AliAnalysisTaskEmcalJetShapesMC::kRecoil) task->SetPtTriggerSelections(minpTHTrigger, maxpTHTrigger);

  TString thename(njetsBase);

  AliParticleContainer *trackContPartLevel = task->AddMCParticleContainer(ntracksPartLevel);

  AliJetContainer *jetContBase=0x0;
  TString strType(type);

  jetContBase = task->AddJetContainer(njetsBase,strType,jetradius);
  if(jetContBase) {
    jetContBase->SetRhoName(rhoName);
    jetContBase->ConnectParticleContainer(trackContPartLevel);
    //jetContBase->ConnectClusterContainer(clusterCont);
    jetContBase->SetPercAreaCut(0.6);
  }

  task->SetCaloTriggerPatchInfoName(kEmcalTriggers.Data());
  task->SetCentralityEstimator(CentEst);
  task->SelectCollisionCandidates(pSel);
  task->SetUseAliAnaUtils(kFALSE);

  mgr->AddTask(task);

  //Connnect input
  mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer() );

  //Connect output
  TString contName1(wagonName1);
  TString contName2(wagonName2);

  if (jetShapeType == AliAnalysisTaskEmcalJetShapesMC::kGenShapes) {
    contName1 += "_GenShapes";
    contName2 += "_GenShapes";
  }

  switch (jetShapeSub) {

    case AliAnalysisTaskEmcalJetShapesMC::kNoSub:
      contName1 += "_NoSub";
      contName2 += "_NoSub";
      break;

    case AliAnalysisTaskEmcalJetShapesMC::kConstSub:
      contName1 += "_ConstSub";
      contName2 += "_ConstSub";
      break;

    case AliAnalysisTaskEmcalJetShapesMC::kDerivSub:
      contName1 += "_DerivSub";
      contName2 += "_DerivSub";
      break;


  }

  switch (jetSelection) {
    case AliAnalysisTaskEmcalJetShapesMC::kInclusive:
      contName1 += "_Incl";
      contName2 += "_Incl";
      break;


    case AliAnalysisTaskEmcalJetShapesMC::kRecoil:
      TString recoilTriggerString = Form("_Recoil_%.0f_%0.f", minpTHTrigger, maxpTHTrigger);
      contName1 += recoilTriggerString;
      contName2 += recoilTriggerString;

      break;

  }


  TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());


  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName1.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(contName2.Data(),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());

  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);

  return task;

}
