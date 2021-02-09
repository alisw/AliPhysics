//
// Subjet fraction analysis task.
//
//Nima Zardoshti

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TCanvas.h>
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
//#include "AliPythiaInfo.h"
#include "TRandom3.h"
#include "AliPicoTrack.h"
#include "AliEmcalJetFinder.h"
#include "AliAODEvent.h"
#include "AliAnalysisTaskSubJetFraction.h"
#include "AliEmcalPythiaInfo.h"

#include "FJ_includes.h"

//Globals






using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskSubJetFraction)

//________________________________________________________________________
AliAnalysisTaskSubJetFraction::AliAnalysisTaskSubJetFraction() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskSubJetFraction", kTRUE),
  fContainer(0),
  fMinFractionShared(0),
  fJetShapeType(kData),
  fJetShapeSub(kNoSub),
  fJetSelection(kInclusive),
  fTreeSize(nVar),
  fPtThreshold(-9999.),
  fRMatching(0.2),
  fPtMinTriggerHadron(20.),
  fPtMaxTriggerHadron(50.),
  fRecoilAngularWindow(0.6),
  fSemigoodCorrect(0),
  fHolePos(0),
  fHoleWidth(0),
  fRandom(0x0),
  fCentSelectOn(kTRUE),
  fCentMin(0),
  fCentMax(10),
  fSubJetAlgorithm(0),
  fSubJetRadius(0.1),
  fSubJetMinPt(1),
  fJetRadius(0.4),
  fRMatched(0.2),
  fSharedFractionPtMin(0.5),
  fDerivSubtrOrder(0),
  fFullTree(kFALSE),
  fBeta_SD(0),
  fZCut(0.1),
  fReclusteringAlgorithm(0),
  fSoftDropOn(0),
  fMLOn(0),
  fNsubMeasure(kFALSE),
  fRandomisationEqualPt(kFALSE),
  fhPtTriggerHadron(0x0),
  fhJetPt(0x0),
  fhJetPt_1(0x0),
  fhJetPt_2(0x0),
  fhJetPhi(0x0),
  fhJetPhi_1(0x0),
  fhJetPhi_2(0x0),
  fhJetEta(0x0),
  fhJetEta_1(0x0),
  fhJetEta_2(0x0),
  fhJetMass(0x0),
  fhJetMass_1(0x0),
  fhJetMass_2(0x0),
  fhJetRadius(0x0),
  fhJetRadius_1(0x0),
  fhJetRadius_2(0x0),
  fhJetCounter(0x0),
  fhJetCounter_1(0x0),
  fhJetCounter_2(0x0),
  fhNumberOfJetTracks(0x0),
  fhNumberOfJetTracks_1(0x0),
  fhNumberOfJetTracks_2(0x0),
  fhSubJetPt(0x0),
  fhSubJetPt_1(0x0),
  fhSubJetPt_2(0x0),
  fhSubJetMass(0x0),
  fhSubJetMass_1(0x0),
  fhSubJetMass_2(0x0),
  fhSubJetRadius(0x0),
  fhSubJetRadius_1(0x0),
  fhSubJetRadius_2(0x0),
  fhSubJetCounter(0x0),
  fhSubJetCounter_1(0x0),
  fhSubJetCounter_2(0x0),
  fhNumberOfSubJetTracks(0x0),
  fhNumberOfSubJetTracks_1(0x0),
  fhNumberOfSubJetTracks_2(0x0),
  fh2PtTriggerHadronJet(0x0),
  fhPhiTriggerHadronJet(0x0),
  fhPhiTriggerHadronEventPlane(0x0),
  fhPhiTriggerHadronEventPlaneTPC(0x0),
  fhTrackPhi(0x0),
  fhTrackPhi_Cut(0x0),
  fh2PtRatio(0x0),
  fhEventCounter(0x0),
  fhEventCounter_1(0x0),
  fhEventCounter_2(0x0),
  fhSubJettiness1CheckRatio_FJ_AKT(0x0),
  fhSubJettiness1CheckRatio_FJ_KT(0x0),
  fhSubJettiness1CheckRatio_FJ_CA(0x0),
  fhSubJettiness1CheckRatio_FJ_WTA_KT(0x0),
  fhSubJettiness1CheckRatio_FJ_WTA_CA(0x0),
  fhSubJettiness1CheckRatio_FJ_OP_AKT(0x0),
  fhSubJettiness1CheckRatio_FJ_OP_KT(0x0),
  fhSubJettiness1CheckRatio_FJ_OP_CA(0x0),
  fhSubJettiness1CheckRatio_FJ_OP_WTA_KT(0x0),
  fhSubJettiness1CheckRatio_FJ_OP_WTA_CA(0x0),
  fhSubJettiness1CheckRatio_FJ_MIN(0x0),
  fhSubJettiness2CheckRatio_FJ_AKT(0x0),
  fhSubJettiness2CheckRatio_FJ_KT(0x0),
  fhSubJettiness2CheckRatio_FJ_CA(0x0),
  fhSubJettiness2CheckRatio_FJ_WTA_KT(0x0),
  fhSubJettiness2CheckRatio_FJ_WTA_CA(0x0),
  fhSubJettiness2CheckRatio_FJ_OP_AKT(0x0),
  fhSubJettiness2CheckRatio_FJ_OP_KT(0x0),
  fhSubJettiness2CheckRatio_FJ_OP_CA(0x0),
  fhSubJettiness2CheckRatio_FJ_OP_WTA_KT(0x0),
  fhSubJettiness2CheckRatio_FJ_OP_WTA_CA(0x0),
  fhSubJettiness2CheckRatio_FJ_MIN(0x0),
  fhSubJettiness2to1CheckRatio_FJ_AKT(0x0),
  fhSubJettiness2to1CheckRatio_FJ_KT(0x0),
  fhSubJettiness2to1CheckRatio_FJ_CA(0x0),
  fhSubJettiness2to1CheckRatio_FJ_WTA_KT(0x0),
  fhSubJettiness2to1CheckRatio_FJ_WTA_CA(0x0),
  fhSubJettiness2to1CheckRatio_FJ_OP_AKT(0x0),
  fhSubJettiness2to1CheckRatio_FJ_OP_KT(0x0),
  fhSubJettiness2to1CheckRatio_FJ_OP_CA(0x0),
  fhSubJettiness2to1CheckRatio_FJ_OP_WTA_KT(0x0),
  fhSubJettiness2to1CheckRatio_FJ_OP_WTA_CA(0x0),
  fhSubJettiness2to1CheckRatio_FJ_MIN(0x0),
  fhSubJettiness1_FJ_AKT(0x0),
  fhSubJettiness1_FJ_KT(0x0),
  fhSubJettiness1_FJ_CA(0x0),
  fhSubJettiness1_FJ_WTA_KT(0x0),
  fhSubJettiness1_FJ_WTA_CA(0x0),
  fhSubJettiness1_FJ_OP_AKT(0x0),
  fhSubJettiness1_FJ_OP_KT(0x0),
  fhSubJettiness1_FJ_OP_CA(0x0),
  fhSubJettiness1_FJ_OP_WTA_KT(0x0),
  fhSubJettiness1_FJ_OP_WTA_CA(0x0),
  fhSubJettiness1_FJ_MIN(0x0),
  fhSubJettiness2_FJ_AKT(0x0),
  fhSubJettiness2_FJ_KT(0x0),
  fhSubJettiness2_FJ_CA(0x0),
  fhSubJettiness2_FJ_WTA_KT(0x0),
  fhSubJettiness2_FJ_WTA_CA(0x0),
  fhSubJettiness2_FJ_OP_AKT(0x0),
  fhSubJettiness2_FJ_OP_KT(0x0),
  fhSubJettiness2_FJ_OP_CA(0x0),
  fhSubJettiness2_FJ_OP_WTA_KT(0x0),
  fhSubJettiness2_FJ_OP_WTA_CA(0x0),
  fhSubJettiness2_FJ_MIN(0x0),
  fhSubJettiness2to1_FJ_AKT(0x0),
  fhSubJettiness2to1_FJ_KT(0x0),
  fhSubJettiness2to1_FJ_CA(0x0),
  fhSubJettiness2to1_FJ_WTA_KT(0x0),
  fhSubJettiness2to1_FJ_WTA_CA(0x0),
  fhSubJettiness2to1_FJ_OP_AKT(0x0),
  fhSubJettiness2to1_FJ_OP_KT(0x0),
  fhSubJettiness2to1_FJ_OP_CA(0x0),
  fhSubJettiness2to1_FJ_OP_WTA_KT(0x0),
  fhSubJettiness2to1_FJ_OP_WTA_CA(0x0),
  fhSubJettiness2to1_FJ_MIN(0x0),
  fShapesVar_Tracks_Rec(0),
  fShapesVar_Tracks_Truth(0),
  fTreeResponseMatrixAxis(0),
  fTreeTracks(0)

{
  for(Int_t i=0;i<nVar;i++){
    fShapesVar[i]=0;
  }
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskSubJetFraction::AliAnalysisTaskSubJetFraction(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fContainer(0),
  fMinFractionShared(0),
  fJetShapeType(kData),
  fJetShapeSub(kNoSub),
  fJetSelection(kInclusive),
  fTreeSize(nVar),
  fPtThreshold(-9999.),
  fRMatching(0.2),
  fPtMinTriggerHadron(20.),
  fPtMaxTriggerHadron(50.),
  fRecoilAngularWindow(0.6),
  fSemigoodCorrect(0),
  fHolePos(0),
  fHoleWidth(0),
  fRandom(0x0),
  fCentSelectOn(kTRUE),
  fCentMin(0),
  fCentMax(10),
  fSubJetAlgorithm(0),
  fSubJetRadius(0.1),
  fSubJetMinPt(1),
  fJetRadius(0.4),
  fRMatched(0.2),
  fSharedFractionPtMin(0.5),
  fDerivSubtrOrder(0),
  fFullTree(kFALSE),
  fBeta_SD(0),
  fZCut(0.1),
  fReclusteringAlgorithm(0),
  fSoftDropOn(0),
  fMLOn(0),
  fNsubMeasure(kFALSE),
  fRandomisationEqualPt(kFALSE),
  fhPtTriggerHadron(0x0),
  fhJetPt(0x0),
  fhJetPt_1(0x0),
  fhJetPt_2(0x0),
  fhJetPhi(0x0),
  fhJetPhi_1(0x0),
  fhJetPhi_2(0x0),
  fhJetEta(0x0),
  fhJetEta_1(0x0),
  fhJetEta_2(0x0),
  fhJetMass(0x0),
  fhJetMass_1(0x0),
  fhJetMass_2(0x0),
  fhJetRadius(0x0),
  fhJetRadius_1(0x0),
  fhJetRadius_2(0x0),
  fhJetCounter(0x0),
  fhJetCounter_1(0x0),
  fhJetCounter_2(0x0),
  fhNumberOfJetTracks(0x0),
  fhNumberOfJetTracks_1(0x0),
  fhNumberOfJetTracks_2(0x0),
  fhSubJetPt(0x0),
  fhSubJetPt_1(0x0),
  fhSubJetPt_2(0x0),
  fhSubJetMass(0x0),
  fhSubJetMass_1(0x0),
  fhSubJetMass_2(0x0),
  fhSubJetRadius(0x0),
  fhSubJetRadius_1(0x0),
  fhSubJetRadius_2(0x0),
  fhSubJetCounter(0x0),
  fhSubJetCounter_1(0x0),
  fhSubJetCounter_2(0x0),
  fhNumberOfSubJetTracks(0x0),
  fhNumberOfSubJetTracks_1(0x0),
  fhNumberOfSubJetTracks_2(0x0),
  fh2PtTriggerHadronJet(0x0),
  fhPhiTriggerHadronJet(0x0),
  fhPhiTriggerHadronEventPlane(0x0),
  fhPhiTriggerHadronEventPlaneTPC(0x0),
  fhTrackPhi(0x0),
  fhTrackPhi_Cut(0x0),
  fh2PtRatio(0x0),
  fhEventCounter(0x0),
  fhEventCounter_1(0x0),
  fhEventCounter_2(0x0),
  fhSubJettiness1CheckRatio_FJ_AKT(0x0),
  fhSubJettiness1CheckRatio_FJ_KT(0x0),
  fhSubJettiness1CheckRatio_FJ_CA(0x0),
  fhSubJettiness1CheckRatio_FJ_WTA_KT(0x0),
  fhSubJettiness1CheckRatio_FJ_WTA_CA(0x0),
  fhSubJettiness1CheckRatio_FJ_OP_AKT(0x0),
  fhSubJettiness1CheckRatio_FJ_OP_KT(0x0),
  fhSubJettiness1CheckRatio_FJ_OP_CA(0x0),
  fhSubJettiness1CheckRatio_FJ_OP_WTA_KT(0x0),
  fhSubJettiness1CheckRatio_FJ_OP_WTA_CA(0x0),
  fhSubJettiness1CheckRatio_FJ_MIN(0x0),
  fhSubJettiness2CheckRatio_FJ_AKT(0x0),
  fhSubJettiness2CheckRatio_FJ_KT(0x0),
  fhSubJettiness2CheckRatio_FJ_CA(0x0),
  fhSubJettiness2CheckRatio_FJ_WTA_KT(0x0),
  fhSubJettiness2CheckRatio_FJ_WTA_CA(0x0),
  fhSubJettiness2CheckRatio_FJ_OP_AKT(0x0),
  fhSubJettiness2CheckRatio_FJ_OP_KT(0x0),
  fhSubJettiness2CheckRatio_FJ_OP_CA(0x0),
  fhSubJettiness2CheckRatio_FJ_OP_WTA_KT(0x0),
  fhSubJettiness2CheckRatio_FJ_OP_WTA_CA(0x0),
  fhSubJettiness2CheckRatio_FJ_MIN(0x0),
  fhSubJettiness2to1CheckRatio_FJ_AKT(0x0),
  fhSubJettiness2to1CheckRatio_FJ_KT(0x0),
  fhSubJettiness2to1CheckRatio_FJ_CA(0x0),
  fhSubJettiness2to1CheckRatio_FJ_WTA_KT(0x0),
  fhSubJettiness2to1CheckRatio_FJ_WTA_CA(0x0),
  fhSubJettiness2to1CheckRatio_FJ_OP_AKT(0x0),
  fhSubJettiness2to1CheckRatio_FJ_OP_KT(0x0),
  fhSubJettiness2to1CheckRatio_FJ_OP_CA(0x0),
  fhSubJettiness2to1CheckRatio_FJ_OP_WTA_KT(0x0),
  fhSubJettiness2to1CheckRatio_FJ_OP_WTA_CA(0x0),
  fhSubJettiness2to1CheckRatio_FJ_MIN(0x0),
  fhSubJettiness1_FJ_AKT(0x0),
  fhSubJettiness1_FJ_KT(0x0),
  fhSubJettiness1_FJ_CA(0x0),
  fhSubJettiness1_FJ_WTA_KT(0x0),
  fhSubJettiness1_FJ_WTA_CA(0x0),
  fhSubJettiness1_FJ_OP_AKT(0x0),
  fhSubJettiness1_FJ_OP_KT(0x0),
  fhSubJettiness1_FJ_OP_CA(0x0),
  fhSubJettiness1_FJ_OP_WTA_KT(0x0),
  fhSubJettiness1_FJ_OP_WTA_CA(0x0),
  fhSubJettiness1_FJ_MIN(0x0),
  fhSubJettiness2_FJ_AKT(0x0),
  fhSubJettiness2_FJ_KT(0x0),
  fhSubJettiness2_FJ_CA(0x0),
  fhSubJettiness2_FJ_WTA_KT(0x0),
  fhSubJettiness2_FJ_WTA_CA(0x0),
  fhSubJettiness2_FJ_OP_AKT(0x0),
  fhSubJettiness2_FJ_OP_KT(0x0),
  fhSubJettiness2_FJ_OP_CA(0x0),
  fhSubJettiness2_FJ_OP_WTA_KT(0x0),
  fhSubJettiness2_FJ_OP_WTA_CA(0x0),
  fhSubJettiness2_FJ_MIN(0x0),
  fhSubJettiness2to1_FJ_AKT(0x0),
  fhSubJettiness2to1_FJ_KT(0x0),
  fhSubJettiness2to1_FJ_CA(0x0),
  fhSubJettiness2to1_FJ_WTA_KT(0x0),
  fhSubJettiness2to1_FJ_WTA_CA(0x0),
  fhSubJettiness2to1_FJ_OP_AKT(0x0),
  fhSubJettiness2to1_FJ_OP_KT(0x0),
  fhSubJettiness2to1_FJ_OP_CA(0x0),
  fhSubJettiness2to1_FJ_OP_WTA_KT(0x0),
  fhSubJettiness2to1_FJ_OP_WTA_CA(0x0),
  fhSubJettiness2to1_FJ_MIN(0x0),
  fShapesVar_Tracks_Rec(0),
  fShapesVar_Tracks_Truth(0),
  fTreeResponseMatrixAxis(0),
  fTreeTracks(0)
{
  // Standard constructor.
  for(Int_t i=0;i<nVar;i++){
    fShapesVar[i]=0;
  }
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskSubJetFraction::~AliAnalysisTaskSubJetFraction()
{
  // Destructor.
}

//________________________________________________________________________
 void AliAnalysisTaskSubJetFraction::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  TH1::AddDirectory(oldStatus);
  //create a tree used for the MC data and making a 4D response matrix


  

  const char* nameoutput_1 = GetOutputSlot(3)->GetContainer()->GetName();
  fTreeTracks = new TTree(nameoutput_1, nameoutput_1);
  TString *fShapesVarNames_Tracks=new TString[2];
  fShapesVarNames_Tracks[0] = "Jet_Constituents_Det";
  fShapesVarNames_Tracks[1] = "Jet_Constituents_Truth";
  fTreeTracks->Branch(fShapesVarNames_Tracks[0].Data(), &fShapesVar_Tracks_Rec, 0,1);
  fTreeTracks->Branch(fShapesVarNames_Tracks[1].Data(), &fShapesVar_Tracks_Truth, 0,1);
  //fTreeTracks->Branch(fShapesVarNames_Tracks[0].Data(), &fShapesVar_Tracks[0], "dfdf");
  //fTreeTracks->Branch(fShapesVarNames_Tracks[1].Data(), &fShapesVar_Tracks[1], "fdss");
  

  
  const char* nameoutput = GetOutputSlot(2)->GetContainer()->GetName();
  fTreeResponseMatrixAxis = new TTree(nameoutput, nameoutput);
  if (fFullTree){
    TString *fShapesVarNames = new TString [nVar];
    if (fMLOn==0){
      fShapesVarNames[0] = "Pt";
      fShapesVarNames[1] = "Pt_Truth";
      fShapesVarNames[2] = "Tau1";
      fShapesVarNames[3] = "Tau1_Truth";
      fShapesVarNames[4] = "Tau2";
      fShapesVarNames[5] = "Tau2_Truth";
      fShapesVarNames[6] = "Tau3";
      fShapesVarNames[7] = "Tau3_Truth";
      fShapesVarNames[8] = "OpeningAngle";
      fShapesVarNames[9] = "OpeningAngle_Truth";
      fShapesVarNames[10] = "JetMultiplicity";
      fShapesVarNames[11] = "JetMultiplicity_Truth";
      fShapesVarNames[12] = "OpeningAngleSD";
      fShapesVarNames[13] = "OpeningAngleSD_Truth";
      fShapesVarNames[14] = "Zg";
      fShapesVarNames[15] = "Zg_Truth";
      fShapesVarNames[16] = "LeadingTrackPt";
      fShapesVarNames[17] = "LeadingTrackPt_Truth";
      fShapesVarNames[18] = "EventPlaneTriggerHadron";
      fShapesVarNames[19] = "EventPlaneTriggerHadron_Truth";
      fShapesVarNames[20] = "EventPlaneTPCTriggerHadron";
      fShapesVarNames[21] = "EventPlaneTPCTriggerHadron_Truth";
      fShapesVarNames[22] = "DeltaR";
      fShapesVarNames[23] = "DeltaR_Truth";
      fShapesVarNames[24] = "Frac1";
      fShapesVarNames[25] = "Frac1_Truth";
      fShapesVarNames[26] = "Frac2";
      fShapesVarNames[27] = "Frac2_Truth";
    }
    if (fMLOn==1){
      fShapesVarNames[0] = "Pt_Rec";
      fShapesVarNames[1] = "Pt_Truth";
      fShapesVarNames[2] = "Eta_Rec";
      fShapesVarNames[3] = "Eta_Truth";
      fShapesVarNames[4] = "Phi_Rec";
      fShapesVarNames[5] = "Phi_Truth";
      fShapesVarNames[6] = "Mass_Rec";
      fShapesVarNames[7] = "Mass_Truth";
      fShapesVarNames[8] = "JetMultiplicity_Rec";
      fShapesVarNames[9] = "JetMultiplicity_Truth";
      fShapesVarNames[10] = "Parton_1_Flag";
      fShapesVarNames[11] = "Parton_2_Flag";
      fShapesVarNames[12] = "Parton_1_Eta";
      fShapesVarNames[13] = "Parton_2_Eta";
      fShapesVarNames[14] = "Parton_1_Phi";
      fShapesVarNames[15] = "Parton_2_Phi";
      fShapesVarNames[16] = "Angularity";
      fShapesVarNames[17] = "Angularity_Truth";
      fShapesVarNames[18] = "PTD";
      fShapesVarNames[19] = "PTD_Truth";
      fShapesVarNames[20] = "Blank_1";
      fShapesVarNames[21] = "Blank_2";
      fShapesVarNames[22] = "Blank_3";
      fShapesVarNames[23] = "Blank_4";
      fShapesVarNames[24] = "Blank_5";
      fShapesVarNames[25] = "Blank_6";
      fShapesVarNames[26] = "Blank_7";
      fShapesVarNames[27] = "Blank_8";
    }
    for(Int_t ivar=0; ivar < nVar; ivar++){
      cout<<"looping over variables"<<endl;
      fTreeResponseMatrixAxis->Branch(fShapesVarNames[ivar].Data(), &fShapesVar[ivar], Form("%s/D", fShapesVarNames[ivar].Data()));
    }
  }

  Double_t Eta_Up=1.00;
  Double_t Eta_Low=-1.00;
  Int_t Eta_Bins=100;
  Double_t Phi_Up=2*(TMath::Pi());
  Double_t Phi_Low=(-1*(TMath::Pi()));
  Int_t Phi_Bins=100;

  
  if (!fFullTree){
    const Int_t nVarMin = 22; 
    TString *fShapesVarNames = new TString [nVarMin];
    if (fMLOn==0){
      fShapesVarNames[0] = "Pt";
      fShapesVarNames[1] = "Pt_Truth";
      if(fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) fShapesVarNames[1] = "Tau1_Randomisation";
      fShapesVarNames[2] = "Tau1";
      fShapesVarNames[3] = "Tau1_Truth";
      if(fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) fShapesVarNames[3] = "Tau1_ExtraProng_01_10";
      fShapesVarNames[4] = "Tau2";
      fShapesVarNames[5] = "Tau2_Truth";
      if(fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) fShapesVarNames[5] = "Tau1_ExtraProng_01_15";
      fShapesVarNames[6] = "SubJet1LeadingTrackPt";
      if(fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) fShapesVarNames[6] = "Tau1_ExtraProng_01_20";
      fShapesVarNames[7] = "SubJet1LeadingTrackPt_Truth";
      if(fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) fShapesVarNames[7] = "Tau1_ExtraProng_01_25";
      fShapesVarNames[8] = "OpeningAngle";
      fShapesVarNames[9] = "OpeningAngle_Truth";
      //if(fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) fShapesVarNames[9] = "Tau1_ExtraProng_01_30";
      if(fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) fShapesVarNames[9] = "Tau1_kTTracks_1_2_1";
      fShapesVarNames[10] = "SubJet2LeadingTrackPt";
      if(fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) fShapesVarNames[10] = "Tau1_ExtraProng_03_10";
      fShapesVarNames[11] = "SubJet2LeadingTrackPt_Truth";
      if(fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) fShapesVarNames[11] = "Tau1_ExtraProng_03_15";
      fShapesVarNames[12] = "OpeningAngleSD";
      if(fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) fShapesVarNames[12] = "Tau1_ExtraProng_03_20";
      fShapesVarNames[13] = "OpeningAngleSD_Truth";
      if(fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) fShapesVarNames[13] = "Tau2_Randomisation";
      fShapesVarNames[14] = "SubJet1Pt";
      if(fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) fShapesVarNames[14] = "Tau2_ExtraProng_01_10";
      fShapesVarNames[15] = "SubJet1Pt_Truth";
      if(fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) fShapesVarNames[15] = "Tau2_ExtraProng_01_15";
      fShapesVarNames[16] = "LeadingTrackPt";
      if(fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) fShapesVarNames[16] = "Tau2_ExtraProng_01_20";
      fShapesVarNames[17] = "LeadingTrackPt_Truth";
      if(fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) fShapesVarNames[17] = "Tau2_ExtraProng_01_25";
      fShapesVarNames[18] = "EventPlaneTriggerHadron";
      // if(fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) fShapesVarNames[18] = "Tau2_ExtraProng_01_30";
      if(fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) fShapesVarNames[18] = "Tau2_kTTracks_1_2_1";
      fShapesVarNames[19] = "EventPlaneTriggerHadron_Truth";
      if(fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) fShapesVarNames[19] = "Tau2_ExtraProng_03_10";
      fShapesVarNames[20] = "SubJet2Pt";
      if(fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) fShapesVarNames[20] = "Tau2_ExtraProng_03_15";
      fShapesVarNames[21] = "SubJet2Pt_Truth";
      if(fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) fShapesVarNames[21] = "Tau2_ExtraProng_03_20";
    }
    if (fMLOn==1){
      fShapesVarNames[0] = "Pt_Rec";
      fShapesVarNames[1] = "Pt_Truth";
      fShapesVarNames[2] = "Eta_Rec";
      fShapesVarNames[3] = "Eta_Truth";
      fShapesVarNames[4] = "Phi_Rec";
      fShapesVarNames[5] = "Phi_Truth";
      fShapesVarNames[6] = "Mass_Rec";
      fShapesVarNames[7] = "Mass_Truth";
      fShapesVarNames[8] = "JetMultiplicity_Rec";
      fShapesVarNames[9] = "JetMultiplicity_Truth";
      fShapesVarNames[10] = "Parton_1_Flag";
      fShapesVarNames[11] = "Parton_2_Flag";
      fShapesVarNames[12] = "Parton_1_Eta";
      fShapesVarNames[13] = "Parton_2_Eta";
      fShapesVarNames[14] = "Parton_1_Phi";
      fShapesVarNames[15] = "Parton_2_Phi";
      fShapesVarNames[16] = "Angularity";
      fShapesVarNames[17] = "Angularity_Truth";
      fShapesVarNames[18] = "PTD";
      fShapesVarNames[19] = "PTD_Truth";
      fShapesVarNames[20] = "Blank_1";
      fShapesVarNames[21] = "Blank_2";
    }
    
    for(Int_t ivar=0; ivar < nVarMin; ivar++){
      cout<<"looping over variables"<<endl;
      fTreeResponseMatrixAxis->Branch(fShapesVarNames[ivar].Data(), &fShapesVar[ivar], Form("%s/D", fShapesVarNames[ivar].Data()));
    }
  }
  
  if (fJetSelection == AliAnalysisTaskSubJetFraction::kRecoil){
    fhPtTriggerHadron= new TH1F("fhPtTriggerHadron", "fhPtTriggerHadron",1500,-0.5,149.5);  
    fOutput->Add(fhPtTriggerHadron);
    fh2PtTriggerHadronJet= new TH2F("fh2PtTriggerHadronJet", "fh2PtTriggerHadronJet",1500,-0.5,149.5,1500,-0.5,149.5);  
    fOutput->Add(fh2PtTriggerHadronJet);
    fhPhiTriggerHadronJet= new TH1F("fhPhiTriggerHadronJet", "fhPhiTriggerHadronJet",360 , -1.5*(TMath::Pi()), 1.5*(TMath::Pi()));  
    fOutput->Add(fhPhiTriggerHadronJet);
    fhPhiTriggerHadronEventPlane= new TH1F("fhPhiTriggerHadronEventPlane", "fhPhiTriggerHadronEventPlane",360 , -1.5*(TMath::Pi()), 1.5*(TMath::Pi()));  
    fOutput->Add(fhPhiTriggerHadronEventPlane);
    fhPhiTriggerHadronEventPlaneTPC= new TH1F("fhPhiTriggerHadronEventPlaneTPC", "fhPhiTriggerHadronEventPlaneTPC",360 , -1.5*(TMath::Pi()), 1.5*(TMath::Pi()));  
    fOutput->Add(fhPhiTriggerHadronEventPlaneTPC);
  }
  if (fJetShapeType==AliAnalysisTaskSubJetFraction::kData || fJetShapeType==AliAnalysisTaskSubJetFraction::kSim || fJetShapeType==AliAnalysisTaskSubJetFraction::kGenOnTheFly){
    
    //  fhJetPt= new TH1F("fhJetPt", "Jet Pt", (XBinsJetPtSize)-1, XBinsJetPt);
    fhJetPt= new TH1F("fhJetPt", "Jet Pt",1500,-0.5,149.5 );   
    fOutput->Add(fhJetPt);
    //fhJetPhi= new TH1F("fhJetPhi", "Jet Phi", Phi_Bins, Phi_Low, Phi_Up);
    // fhJetPhi= new TH1F("fhJetPhi", "Jet Phi",360 , -1.5*(TMath::Pi()), 1.5*(TMath::Pi()));
    fhJetPhi= new TH1F("fhJetPhi", "Jet Phi",780 , -7, 7);
    fOutput->Add(fhJetPhi);
    fhJetEta= new TH1F("fhJetEta", "Jet Eta", Eta_Bins, Eta_Low, Eta_Up);
    fOutput->Add(fhJetEta);
    fhJetMass= new TH1F("fhJetMass", "Jet Mass", 4000,-0.5, 39.5);
    fOutput->Add(fhJetMass);
    fhJetRadius= new TH1F("fhJetRadius", "Jet Radius", 100, -0.05,0.995);
    fOutput->Add(fhJetRadius);
    fhNumberOfJetTracks= new TH1F("fhNumberOfJetTracks", "Number of Tracks within a Jet", 300, -0.5,299.5);
    fOutput->Add(fhNumberOfJetTracks);
    fhSubJetRadius= new TH1F("fhSubJetRadius", "SubJet Radius", 100, -0.05,0.995);
    fOutput->Add(fhSubJetRadius);
    fhSubJetPt= new TH1F("fhSubJetPt", "SubJet Pt", 1500, -0.5,149.5);
    fOutput->Add(fhSubJetPt);
    fhSubJetMass= new TH1F("fhSubJetMass", "Sub Jet Mass", 4000,-0.5, 39.5);
    fOutput->Add(fhSubJetMass);
    fhNumberOfSubJetTracks= new TH1F("fhNumberOfSubJetTracks", "Number of Tracks within a Sub Jet", 300, -0.5,299.5);
    fOutput->Add(fhNumberOfSubJetTracks);
    fhJetCounter= new TH1F("fhJetCounter", "Jet Counter", 150, -0.5, 149.5);
    fOutput->Add(fhJetCounter);
    fhSubJetCounter = new TH1F("fhSubJetCounter", "SubJet Counter",50, -0.5,49.5);
    fOutput->Add(fhSubJetCounter);
    fhEventCounter= new TH1F("fhEventCounter", "Event Counter", 15,0.5,15.5);
    fOutput->Add(fhEventCounter);
    fhTrackPhi= new TH1F("fhTrackPhi", "fhTrackPhi",780 , -7, 7);   
    fOutput->Add(fhTrackPhi);
    fhTrackPhi_Cut= new TH1F("fhTrackPhi_Cut", "fhTrackPhi_Cut",780 , -7, 7);   
    fOutput->Add(fhTrackPhi_Cut);
  }
  if(fJetShapeType==AliAnalysisTaskSubJetFraction::kSim || fJetShapeType==AliAnalysisTaskSubJetFraction::kGenOnTheFly){
    fhSubJettiness1_FJ_KT= new TH1D("fhSubJettiness1_FJ_KT","fhSubJettiness1_FJ_KT",400,-2,2);
    fOutput->Add(fhSubJettiness1_FJ_KT);
    fhSubJettiness1_FJ_MIN= new TH1D("fhSubJettiness1_FJ_MIN","fhSubJettiness1_FJ_MIN",400,-2,2);
    fOutput->Add(fhSubJettiness1_FJ_MIN);
    fhSubJettiness2_FJ_KT= new TH1D("fhSubJettiness2_FJ_KT","fhSubJettiness2_FJ_KT",400,-2,2);
    fOutput->Add(fhSubJettiness2_FJ_KT);
    fhSubJettiness2_FJ_MIN= new TH1D("fhSubJettiness2_FJ_MIN","fhSubJettiness2_FJ_MIN",400,-2,2);
    fOutput->Add(fhSubJettiness2_FJ_MIN);
    fhSubJettiness1CheckRatio_FJ_AKT = new TH2D("fhSubJettiness1CheckRatio_FJ_AKT","fhSubJettiness1CheckRatio_FJ_AKT",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness1CheckRatio_FJ_AKT);
    fhSubJettiness1CheckRatio_FJ_KT= new TH2D("fhSubJettiness1CheckRatio_FJ_KT","fhSubJettiness1CheckRatio_FJ_KT",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness1CheckRatio_FJ_KT);
    fhSubJettiness1CheckRatio_FJ_CA= new TH2D("fhSubJettiness1CheckRatio_FJ_CA","fhSubJettiness1CheckRatio_FJ_CA",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness1CheckRatio_FJ_CA);
    fhSubJettiness1CheckRatio_FJ_WTA_KT= new TH2D("fhSubJettiness1CheckRatio_FJ_WTA_KT","fhSubJettiness1CheckRatio_FJ_WTA_KT",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness1CheckRatio_FJ_WTA_KT);
    fhSubJettiness1CheckRatio_FJ_WTA_CA= new TH2D("fhSubJettiness1CheckRatio_FJ_WTA_CA","fhSubJettiness1CheckRatio_FJ_WTA_CA",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness1CheckRatio_FJ_WTA_CA);
    fhSubJettiness1CheckRatio_FJ_OP_AKT= new TH2D("fhSubJettiness1CheckRatio_FJ_OP_AKT","fhSubJettiness1CheckRatio_FJ_OP_AKT",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness1CheckRatio_FJ_OP_AKT);
    fhSubJettiness1CheckRatio_FJ_OP_KT= new TH2D("fhSubJettiness1CheckRatio_FJ_OP_KT","fhSubJettiness1CheckRatio_FJ_OP_KT",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness1CheckRatio_FJ_OP_KT);
    fhSubJettiness1CheckRatio_FJ_OP_CA= new TH2D("fhSubJettiness1CheckRatio_FJ_OP_CA","fhSubJettiness1CheckRatio_FJ_OP_CA",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness1CheckRatio_FJ_OP_CA);
    fhSubJettiness1CheckRatio_FJ_OP_WTA_KT= new TH2D("fhSubJettiness1CheckRatio_FJ_OP_WTA_KT","fhSubJettiness1CheckRatio_FJ_OP_WTA_KT",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness1CheckRatio_FJ_OP_WTA_KT);
    fhSubJettiness1CheckRatio_FJ_OP_WTA_CA= new TH2D("fhSubJettiness1CheckRatio_FJ_OP_WTA_CA","fhSubJettiness1CheckRatio_FJ_OP_WTA_CA",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness1CheckRatio_FJ_OP_WTA_CA);
    fhSubJettiness1CheckRatio_FJ_MIN= new TH2D("fhSubJettiness1CheckRatio_FJ_MIN","fhSubJettiness1CheckRatio_FJ_MIN",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness1CheckRatio_FJ_MIN);
    fhSubJettiness2CheckRatio_FJ_AKT = new TH2D("fhSubJettiness2CheckRatio_FJ_AKT","fhSubJettiness2CheckRatio_FJ_AKT",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness2CheckRatio_FJ_AKT);
    fhSubJettiness2CheckRatio_FJ_KT= new TH2D("fhSubJettiness2CheckRatio_FJ_KT","fhSubJettiness2CheckRatio_FJ_KT",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness2CheckRatio_FJ_KT);
    fhSubJettiness2CheckRatio_FJ_CA= new TH2D("fhSubJettiness2CheckRatio_FJ_CA","fhSubJettiness2CheckRatio_FJ_CA",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness2CheckRatio_FJ_CA);
    fhSubJettiness2CheckRatio_FJ_WTA_KT= new TH2D("fhSubJettiness2CheckRatio_FJ_WTA_KT","fhSubJettiness2CheckRatio_FJ_WTA_KT",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness2CheckRatio_FJ_WTA_KT);
    fhSubJettiness2CheckRatio_FJ_WTA_CA= new TH2D("fhSubJettiness2CheckRatio_FJ_WTA_CA","fhSubJettiness2CheckRatio_FJ_WTA_CA",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness2CheckRatio_FJ_WTA_CA);
    fhSubJettiness2CheckRatio_FJ_OP_AKT= new TH2D("fhSubJettiness2CheckRatio_FJ_OP_AKT","fhSubJettiness2CheckRatio_FJ_OP_AKT",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness2CheckRatio_FJ_OP_AKT);
    fhSubJettiness2CheckRatio_FJ_OP_KT= new TH2D("fhSubJettiness2CheckRatio_FJ_OP_KT","fhSubJettiness2CheckRatio_FJ_OP_KT",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness2CheckRatio_FJ_OP_KT);
    fhSubJettiness2CheckRatio_FJ_OP_CA= new TH2D("fhSubJettiness2CheckRatio_FJ_OP_CA","fhSubJettiness2CheckRatio_FJ_OP_CA",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness2CheckRatio_FJ_OP_CA);
    fhSubJettiness2CheckRatio_FJ_OP_WTA_KT= new TH2D("fhSubJettiness2CheckRatio_FJ_OP_WTA_KT","fhSubJettiness2CheckRatio_FJ_OP_WTA_KT",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness2CheckRatio_FJ_OP_WTA_KT);
    fhSubJettiness2CheckRatio_FJ_OP_WTA_CA= new TH2D("fhSubJettiness2CheckRatio_FJ_OP_WTA_CA","fhSubJettiness2CheckRatio_FJ_OP_WTA_CA",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness2CheckRatio_FJ_OP_WTA_CA);
    fhSubJettiness2CheckRatio_FJ_MIN= new TH2D("fhSubJettiness2CheckRatio_FJ_MIN","fhSubJettiness2CheckRatio_FJ_MIN",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness2CheckRatio_FJ_MIN);
    fhSubJettiness2to1CheckRatio_FJ_AKT = new TH2D("fhSubJettiness2to1CheckRatio_FJ_AKT","fhSubJettiness2to1CheckRatio_FJ_AKT",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness2to1CheckRatio_FJ_AKT);
    fhSubJettiness2to1CheckRatio_FJ_KT= new TH2D("fhSubJettiness2to1CheckRatio_FJ_KT","fhSubJettiness2to1CheckRatio_FJ_KT",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness2to1CheckRatio_FJ_KT);
    fhSubJettiness2to1CheckRatio_FJ_CA= new TH2D("fhSubJettiness2to1CheckRatio_FJ_CA","fhSubJettiness2to1CheckRatio_FJ_CA",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness2to1CheckRatio_FJ_CA);
    fhSubJettiness2to1CheckRatio_FJ_WTA_KT= new TH2D("fhSubJettiness2to1CheckRatio_FJ_WTA_KT","fhSubJettiness2to1CheckRatio_FJ_WTA_KT",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness2to1CheckRatio_FJ_WTA_KT);
    fhSubJettiness2to1CheckRatio_FJ_WTA_CA= new TH2D("fhSubJettiness2to1CheckRatio_FJ_WTA_CA","fhSubJettiness2to1CheckRatio_FJ_WTA_CA",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness2to1CheckRatio_FJ_WTA_CA);
    fhSubJettiness2to1CheckRatio_FJ_OP_AKT= new TH2D("fhSubJettiness2to1CheckRatio_FJ_OP_AKT","fhSubJettiness2to1CheckRatio_FJ_OP_AKT",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness2to1CheckRatio_FJ_OP_AKT);
    fhSubJettiness2to1CheckRatio_FJ_OP_KT= new TH2D("fhSubJettiness2to1CheckRatio_FJ_OP_KT","fhSubJettiness2to1CheckRatio_FJ_OP_KT",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness2to1CheckRatio_FJ_OP_KT);
    fhSubJettiness2to1CheckRatio_FJ_OP_CA= new TH2D("fhSubJettiness2to1CheckRatio_FJ_OP_CA","fhSubJettiness2to1CheckRatio_FJ_OP_CA",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness2to1CheckRatio_FJ_OP_CA);
    fhSubJettiness2to1CheckRatio_FJ_OP_WTA_KT= new TH2D("fhSubJettiness2to1CheckRatio_FJ_OP_WTA_KT","fhSubJettiness2to1CheckRatio_FJ_OP_WTA_KT",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness2to1CheckRatio_FJ_OP_WTA_KT);
    fhSubJettiness2to1CheckRatio_FJ_OP_WTA_CA= new TH2D("fhSubJettiness2to1CheckRatio_FJ_OP_WTA_CA","fhSubJettiness2to1CheckRatio_FJ_OP_WTA_CA",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness2to1CheckRatio_FJ_OP_WTA_CA);
    fhSubJettiness2to1CheckRatio_FJ_MIN= new TH2D("fhSubJettiness2to1CheckRatio_FJ_MIN","fhSubJettiness2to1CheckRatio_FJ_MIN",400,-2,2,300,-1,2);
    fOutput->Add(fhSubJettiness2to1CheckRatio_FJ_MIN);
    fhSubJettiness1_FJ_AKT= new TH1D("fhSubJettiness1_FJ_AKT","fhSubJettiness1_FJ_AKT",400,-2,2);
    fOutput->Add(fhSubJettiness1_FJ_AKT);
    fhSubJettiness1_FJ_CA= new TH1D("fhSubJettiness1_FJ_CA","fhSubJettiness1_FJ_CA",400,-2,2);
    fOutput->Add(fhSubJettiness1_FJ_CA);
    fhSubJettiness1_FJ_WTA_KT= new TH1D("fhSubJettiness1_FJ_WTA_KT","fhSubJettiness1_FJ_WTA_KT",400,-2,2);
    fOutput->Add(fhSubJettiness1_FJ_WTA_KT);
    fhSubJettiness1_FJ_WTA_CA= new TH1D("fhSubJettiness1_FJ_WTA_CA","fhSubJettiness1_FJ_WTA_CA",400,-2,2);
    fOutput->Add(fhSubJettiness1_FJ_WTA_CA);
    fhSubJettiness1_FJ_OP_AKT= new TH1D("fhSubJettiness1_FJ_OP_AKT","fhSubJettiness1_FJ_OP_AKT",400,-2,2);
    fOutput->Add(fhSubJettiness1_FJ_OP_AKT);
    fhSubJettiness1_FJ_OP_KT= new TH1D("fhSubJettiness1_FJ_OP_KT","fhSubJettiness1_FJ_OP_KT",400,-2,2);
    fOutput->Add(fhSubJettiness1_FJ_OP_KT);
    fhSubJettiness1_FJ_OP_CA= new TH1D("fhSubJettiness1_FJ_OP_CA","fhSubJettiness1_FJ_OP_CA",400,-2,2);
    fOutput->Add(fhSubJettiness1_FJ_OP_CA);
    fhSubJettiness1_FJ_OP_WTA_KT= new TH1D("fhSubJettiness1_FJ_OP_WTA_KT","fhSubJettiness1_FJ_OP_WTA_KT",400,-2,2);
    fOutput->Add(fhSubJettiness1_FJ_OP_WTA_KT);
    fhSubJettiness1_FJ_OP_WTA_CA= new TH1D("fhSubJettiness1_FJ_OP_WTA_CA","fhSubJettiness1_FJ_OP_WTA_CA",400,-2,2);
    fOutput->Add(fhSubJettiness1_FJ_OP_WTA_CA);
    fhSubJettiness2_FJ_AKT= new TH1D("fhSubJettiness2_FJ_AKT","fhSubJettiness2_FJ_AKT",400,-2,2);
    fOutput->Add(fhSubJettiness2_FJ_AKT);
    fhSubJettiness2_FJ_CA= new TH1D("fhSubJettiness2_FJ_CA","fhSubJettiness2_FJ_CA",400,-2,2);
    fOutput->Add(fhSubJettiness2_FJ_CA);
    fhSubJettiness2_FJ_WTA_KT= new TH1D("fhSubJettiness2_FJ_WTA_KT","fhSubJettiness2_FJ_WTA_KT",400,-2,2);
    fOutput->Add(fhSubJettiness2_FJ_WTA_KT);
    fhSubJettiness2_FJ_WTA_CA= new TH1D("fhSubJettiness2_FJ_WTA_CA","fhSubJettiness2_FJ_WTA_CA",400,-2,2);
    fOutput->Add(fhSubJettiness2_FJ_WTA_CA);
    fhSubJettiness2_FJ_OP_AKT= new TH1D("fhSubJettiness2_FJ_OP_AKT","fhSubJettiness2_FJ_OP_AKT",400,-2,2);
    fOutput->Add(fhSubJettiness2_FJ_OP_AKT);
    fhSubJettiness2_FJ_OP_KT= new TH1D("fhSubJettiness2_FJ_OP_KT","fhSubJettiness2_FJ_OP_KT",400,-2,2);
    fOutput->Add(fhSubJettiness2_FJ_OP_KT);
    fhSubJettiness2_FJ_OP_CA= new TH1D("fhSubJettiness2_FJ_OP_CA","fhSubJettiness2_FJ_OP_CA",400,-2,2);
    fOutput->Add(fhSubJettiness2_FJ_OP_CA);
    fhSubJettiness2_FJ_OP_WTA_KT= new TH1D("fhSubJettiness2_FJ_OP_WTA_KT","fhSubJettiness2_FJ_OP_WTA_KT",400,-2,2);
    fOutput->Add(fhSubJettiness2_FJ_OP_WTA_KT);
    fhSubJettiness2_FJ_OP_WTA_CA= new TH1D("fhSubJettiness2_FJ_OP_WTA_CA","fhSubJettiness2_FJ_OP_WTA_CA",400,-2,2);
    fOutput->Add(fhSubJettiness2_FJ_OP_WTA_CA);
    fhSubJettiness2to1_FJ_KT= new TH1D("fhSubJettiness2to1_FJ_KT","fhSubJettiness2to1_FJ_KT",400,-2,2);
    fOutput->Add(fhSubJettiness2to1_FJ_KT);
    fhSubJettiness2to1_FJ_AKT= new TH1D("fhSubJettiness2to1_FJ_AKT","fhSubJettiness2to1_FJ_AKT",400,-2,2);
    fOutput->Add(fhSubJettiness2to1_FJ_AKT);
    fhSubJettiness2to1_FJ_CA= new TH1D("fhSubJettiness2to1_FJ_CA","fhSubJettiness2to1_FJ_CA",400,-2,2);
    fOutput->Add(fhSubJettiness2to1_FJ_CA);
    fhSubJettiness2to1_FJ_WTA_KT= new TH1D("fhSubJettiness2to1_FJ_WTA_KT","fhSubJettiness2to1_FJ_WTA_KT",400,-2,2);
    fOutput->Add(fhSubJettiness2to1_FJ_WTA_KT);
    fhSubJettiness2to1_FJ_WTA_CA= new TH1D("fhSubJettiness2to1_FJ_WTA_CA","fhSubJettiness2to1_FJ_WTA_CA",400,-2,2);
    fOutput->Add(fhSubJettiness2to1_FJ_WTA_CA);
    fhSubJettiness2to1_FJ_OP_AKT= new TH1D("fhSubJettiness2to1_FJ_OP_AKT","fhSubJettiness2to1_FJ_OP_AKT",400,-2,2);
    fOutput->Add(fhSubJettiness2to1_FJ_OP_AKT);
    fhSubJettiness2to1_FJ_OP_KT= new TH1D("fhSubJettiness2to1_FJ_OP_KT","fhSubJettiness2to1_FJ_OP_KT",400,-2,2);
    fOutput->Add(fhSubJettiness2to1_FJ_OP_KT);
    fhSubJettiness2to1_FJ_OP_CA= new TH1D("fhSubJettiness2to1_FJ_OP_CA","fhSubJettiness2to1_FJ_OP_CA",400,-2,2);
    fOutput->Add(fhSubJettiness2to1_FJ_OP_CA);
    fhSubJettiness2to1_FJ_OP_WTA_KT= new TH1D("fhSubJettiness2to1_FJ_OP_WTA_KT","fhSubJettiness2to1_FJ_OP_WTA_KT",400,-2,2);
    fOutput->Add(fhSubJettiness2to1_FJ_OP_WTA_KT);
    fhSubJettiness2to1_FJ_OP_WTA_CA= new TH1D("fhSubJettiness2to1_FJ_OP_WTA_CA","fhSubJettiness2to1_FJ_OP_WTA_CA",400,-2,2);
    fOutput->Add(fhSubJettiness2to1_FJ_OP_WTA_CA);
    fhSubJettiness2to1_FJ_MIN= new TH1D("fhSubJettiness2to1_FJ_MIN","fhSubJettiness2to1_FJ_MIN",400,-2,2);
    fOutput->Add(fhSubJettiness2to1_FJ_MIN);
  }
  if(fJetShapeType==AliAnalysisTaskSubJetFraction::kTrueDet){
    fhJetPt_1= new TH1F("fhJetPt_1", "Jet Pt Detector Level",1500,-0.5,149.5 );
    fOutput->Add(fhJetPt_1);
    fhJetPt_2= new TH1F("fhJetPt_2", "Jet Pt Particle Level",1500,-0.5,149.5 );
    fOutput->Add(fhJetPt_2);
    fhJetPhi_1= new TH1F("fhJetPhi_1", "Jet Phi Detector Level",360 , -1.5*(TMath::Pi()), 1.5*(TMath::Pi()));
    fOutput->Add(fhJetPhi_1);
    fhJetPhi_2= new TH1F("fhJetPhi_2", "Jet Phi Particle Level",360 , -1.5*(TMath::Pi()), 1.5*(TMath::Pi()));
    fOutput->Add(fhJetPhi_2);
    fhJetEta_1= new TH1F("fhJetEta_1", "Jet Eta Detector Level", Eta_Bins, Eta_Low, Eta_Up);
    fOutput->Add(fhJetEta_1);
    fhJetEta_2= new TH1F("fhJetEta_2", "Jet Eta Particle Level", Eta_Bins, Eta_Low, Eta_Up);
    fOutput->Add(fhJetEta_2);
    fhJetMass_1= new TH1F("fhJetMass_1", "Jet Mass Detector Level", 4000,-0.5, 39.5);
    fOutput->Add(fhJetMass_1);
    fhJetMass_2= new TH1F("fhJetMass_2", "Jet Mass Particle Level", 4000,-0.5, 39.5);
    fOutput->Add(fhJetMass_2);
    fhJetRadius_1= new TH1F("fhJetRadius_1", "Jet Radius Detector Level", 100, -0.05,0.995);
    fOutput->Add(fhJetRadius_1);
    fhJetRadius_2= new TH1F("fhJetRadius_2", "Jet Radius Particle Level", 100, -0.05,0.995);
    fOutput->Add(fhJetRadius_2);
    fhNumberOfJetTracks_1= new TH1F("fhNumberOfJetTracks_1", "Number of Tracks within a Jet Detector Level", 300, -0.5,299.5);
    fOutput->Add(fhNumberOfJetTracks_1);
    fhNumberOfJetTracks_2= new TH1F("fhNumberOfJetTracks_2", "Number of Tracks within a Jet Particle Level", 300, -0.5,299.5);
    fOutput->Add(fhNumberOfJetTracks_2);
    fhSubJetRadius_1= new TH1F("fhSubJetRadius_1", "SubJet Radius Detector Level", 100, -0.05,0.995);
    fOutput->Add(fhSubJetRadius_1);
    fhSubJetRadius_2= new TH1F("fhSubJetRadius_2", "SubJet Radius Particle Level", 100, -0.05,0.995);
    fOutput->Add(fhSubJetRadius_2);
    fhSubJetPt_1= new TH1F("fhSubJetPt_1", "SubJet Pt Detector Level", 1500, -0.5,149.5);
    fOutput->Add(fhSubJetPt_1);
    fhSubJetPt_2= new TH1F("fhSubJetPt_2", "SubJet Pt Particle Level", 1500, -0.5,149.5);
    fOutput->Add(fhSubJetPt_2);
    fhSubJetMass_1= new TH1F("fhSubJetMass_1", "Sub Jet Mass Detector Level", 4000,-0.5, 39.5);
    fOutput->Add(fhSubJetMass_1);
    fhSubJetMass_2= new TH1F("fhSubJetMass_2", "Sub Jet Mass Particle Level", 4000,-0.5, 39.5);
    fOutput->Add(fhSubJetMass_2);
    fhNumberOfSubJetTracks_1= new TH1F("fhNumberOfSubJetTracks_1", "Number of Tracks within a Sub Jet Detector Level", 300, -0.5,299.5);
    fOutput->Add(fhNumberOfSubJetTracks_1);
    fhNumberOfSubJetTracks_2= new TH1F("fhNumberOfSubJetTracks_2", "Number of Tracks within a Sub Jet Particle Level", 300, -0.5,299.5);
    fOutput->Add(fhNumberOfSubJetTracks_2);
    fhJetCounter_1= new TH1F("fhJetCounter_1", "Jet Counter Detector Level", 150, -0.5, 149.5);
    fOutput->Add(fhJetCounter_1);
    fhJetCounter_2= new TH1F("fhJetCounter_2", "Jet Counter Particle Level", 150, -0.5, 149.5);
    fOutput->Add(fhJetCounter_2);
    fhSubJetCounter_1 = new TH1F("fhSubJetCounter_1", "SubJet Counter Detector Level",50, -0.5,49.5);
    fOutput->Add(fhSubJetCounter_1);
    fhSubJetCounter_2 = new TH1F("fhSubJetCounter_2", "SubJet Counter Particle Level",50, -0.5,49.5);
    fOutput->Add(fhSubJetCounter_2);
    fh2PtRatio= new TH2F("fhPtRatio", "MC pt for detector level vs particle level jets",1500,-0.5,149.5,1500,-0.5,149.5);
    fOutput->Add(fh2PtRatio);
    fhEventCounter_1= new TH1F("fhEventCounter_1", "Event Counter Detector Level", 15,0.5,15.5);
    fOutput->Add(fhEventCounter_1);
    fhEventCounter_2= new TH1F("fhEventCounter_2", "Event Counter Particle Level", 15,0.5,15.5);
    fOutput->Add(fhEventCounter_2);
    fhTrackPhi= new TH1F("fhTrackPhi", "fhTrackPhi",780 , -7, 7);   
    fOutput->Add(fhTrackPhi);
    fhTrackPhi_Cut= new TH1F("fhTrackPhi_Cut", "fhTrackPhi_Cut",780 , -7, 7);   
    fOutput->Add(fhTrackPhi_Cut);
  }
  if (fJetShapeType == AliAnalysisTaskSubJetFraction::kDetEmbPart){
    fhEventCounter= new TH1F("fhEventCounter", "Event Counter", 15,0.5,15.5);
    fOutput->Add(fhEventCounter);
  }
  
  PostData(1,fOutput);
  PostData(2,fTreeResponseMatrixAxis);
  PostData(3,fTreeTracks);
  // delete [] fShapesVarNames;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSubJetFraction::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().


  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSubJetFraction::FillHistograms()
{

  //fhEventCounter Key:
  // 1: Number of events with a Jet Container 
  // 2: Number of Jets without a Jet Container
  // 3:
  // 4: Number of Jets found in all events
  // 5: Number of Jets that were reclustered in all events
  // 6: Number of SubJets found in all events 
  // 7: Number of events were primary vertext was found
  // 8: Number of Jets with more than one SubJet
  // 9:Number of Jets with more than two SubJets
  // 12:Number of SubJetinessEvents in kData


  //  const AliVVertex *vert = InputEvent()->GetPrimaryVertex();
  //  Double_t dVtx[3]={vert->GetX(),vert->GetY(),vert->GetZ()};
  // if(vert) fhEventCounter->Fill(7);

  //cout << ((AliVAODHeader*)(InputEvent()->GetHeader()))->GetEventplane()<< "    "<<fEPV0<<endl;
  // cout << InputEvent()->GetEventplane()->GetEventplane("Q")<< "    "<<fEPV0<<endl;   
  if (fCentSelectOn){
    if ((fCent>fCentMax) || (fCent<fCentMin)) return 0;
  }

  ////Filling Track Phi
  AliTrackContainer *PartCont_Particles = NULL;
  AliParticleContainer *PartContMC_Particles = NULL;
  
  if (fJetShapeSub==kConstSub){
    if (fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) PartContMC_Particles = GetParticleContainer(1);
    else PartCont_Particles = GetTrackContainer(1);
  }
  else{
    if (fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) PartContMC_Particles = GetParticleContainer(0);
    else PartCont_Particles = GetTrackContainer(0);
  }
  
  TClonesArray *TracksArray_Particles = NULL;
  TClonesArray *TracksArrayMC_Particles = NULL;
  if (fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly && PartContMC_Particles) TracksArrayMC_Particles = PartContMC_Particles->GetArray();
  else if(PartCont_Particles) TracksArray_Particles = PartCont_Particles->GetArray();

  AliAODTrack *Track_Particles = 0x0;
  Int_t NTracks_Particles=0;
  if (fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly && TracksArrayMC_Particles) NTracks_Particles = TracksArrayMC_Particles->GetEntriesFast();
  else if(TracksArray_Particles) NTracks_Particles = TracksArray_Particles->GetEntriesFast();

  if (fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly){
    if (PartContMC_Particles && TracksArrayMC_Particles){
      for(Int_t i=0; i < NTracks_Particles; i++){
	if((Track_Particles = static_cast<AliAODTrack*>(PartContMC_Particles->GetAcceptParticle(i)))){
	  if (!Track_Particles) continue;
	  if(TMath::Abs(Track_Particles->Eta())>0.9) continue;
	  if (Track_Particles->Pt()<0.15) continue;
	  fhTrackPhi->Fill(Track_Particles->Phi());
	  if (Track_Particles->Pt()>=4.0) fhTrackPhi_Cut->Fill(Track_Particles->Phi());
	}
      }
    }
  }
  else{
    if (PartCont_Particles && TracksArray_Particles){
      for(Int_t i=0; i < NTracks_Particles; i++){
	if((Track_Particles = static_cast<AliAODTrack*>(PartCont_Particles->GetAcceptTrack(i)))){
	  if (!Track_Particles) continue;
	  if(TMath::Abs(Track_Particles->Eta())>0.9) continue;
	  if (Track_Particles->Pt()<0.15) continue;
	  fhTrackPhi->Fill(Track_Particles->Phi());
	  if (Track_Particles->Pt()>=4.0) fhTrackPhi_Cut->Fill(Track_Particles->Phi());
	}
      }
    }
  }
  ////

  AliAODTrack *TriggerHadron = 0x0;
  if (fJetSelection == kRecoil) {
    //you have to set a flag and the limits of the pT interval for your trigger
    Int_t TriggerHadronLabel = SelectTriggerHadron(fPtMinTriggerHadron, fPtMaxTriggerHadron);    
    if (TriggerHadronLabel==-99999) return 0;  //Trigger Hadron Not Found
    AliTrackContainer *PartCont =NULL;
    AliParticleContainer *PartContMC=NULL;
    if (fJetShapeSub==kConstSub){
      if (fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) PartContMC = GetParticleContainer(1);
      else PartCont = GetTrackContainer(1);
    }
    else{
      if (fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) PartContMC = GetParticleContainer(0);
      else PartCont = GetTrackContainer(0);
    }
    TClonesArray *TrackArray = NULL;
    TClonesArray *TrackArrayMC = NULL;
    if (fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) TrackArrayMC = PartContMC->GetArray();
    else TrackArray = PartCont->GetArray();    
    if (fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) TriggerHadron = static_cast<AliAODTrack*>(TrackArrayMC->At(TriggerHadronLabel));
    else TriggerHadron = static_cast<AliAODTrack*>(TrackArray->At(TriggerHadronLabel));
    if (!TriggerHadron) return 0;//No trigger hadron with label found   
    if(fSemigoodCorrect){
      Double_t HoleDistance=RelativePhi(TriggerHadron->Phi(),fHolePos);
      if(TMath::Abs(HoleDistance)+fHoleWidth+fJetRadius>TMath::Pi()-fRecoilAngularWindow) return 0;
    }
    fhPtTriggerHadron->Fill(TriggerHadron->Pt()); //Needed for per trigger Normalisation
    if (fJetShapeType != AliAnalysisTaskSubJetFraction::kGenOnTheFly) fhPhiTriggerHadronEventPlane->Fill(TMath::Abs(RelativePhiEventPlane(fEPV0,TriggerHadron->Phi()))); //fEPV0 is the event plane from AliAnalysisTaskEmcal
    if (fJetShapeType != AliAnalysisTaskSubJetFraction::kGenOnTheFly) fhPhiTriggerHadronEventPlaneTPC->Fill(TMath::Abs(RelativePhiEventPlane(((AliVAODHeader*)(InputEvent()->GetHeader()))->GetEventplane(),TriggerHadron->Phi()))); //TPC event plane 
  }

  
  ////////////////////////////////////Embedding////////////////////////////////////////
  if (fJetShapeType == AliAnalysisTaskSubJetFraction::kDetEmbPart){
    AliEmcalJet *Jet1 = NULL; //Subtracted hybrid Jet  
    AliEmcalJet *Jet2 = NULL; //Unsubtracted Hybrid Jet                                                                                                                     
    AliEmcalJet *Jet3 = NULL; //Detector Level Pythia Jet
    AliEmcalJet *Jet4 = NULL; //Particle Level Pyhtia Jet
    AliJetContainer *JetCont1= GetJetContainer(0); //Jet Container for Subtracted Hybrid Jets 
    AliJetContainer *JetCont2= GetJetContainer(1); //Jet Container for Unsubtracted Hybrid Jets                                                                                  
    AliJetContainer *JetCont3= GetJetContainer(2); //Jet Container for Detector Level Pyhtia Jets 
    AliJetContainer *JetCont4= GetJetContainer(3); //Jet Container for Particle Level Pythia Jets
    AliEmcalJetFinder *Reclusterer1; //Object containg Subjets from Subtracted Hybrid Jets     
    AliEmcalJetFinder *Reclusterer4; //Object containg Subjets from Particle Level Pythia Jets  
    Bool_t JetsMatched=kFALSE;
    Bool_t EventCounter=kFALSE;
    Int_t JetNumber=-1;
    Double_t JetPtThreshold=-2;
    const AliEmcalPythiaInfo *Parton_Info = 0x0;
    Parton_Info=GetPythiaInfo();
    fhEventCounter->Fill(1);
    if(JetCont1) {
      fhEventCounter->Fill(2);
      JetCont1->ResetCurrentID();
      while((Jet1=JetCont1->GetNextAcceptJet())) {
	if (fJetShapeSub==kConstSub) JetPtThreshold=Jet1->Pt();
	else JetPtThreshold=Jet1->Pt()-(GetRhoVal(0)*Jet1->Area());
        if ( (!Jet1) || (JetPtThreshold<fPtThreshold)) continue;
        else {
	  /*  if(fSemigoodCorrect){
	    Double_t HoleDistance=RelativePhi(Jet1->Phi(),fHolePos);
	    if(TMath::Abs(HoleDistance)<fHoleWidth) continue;
	    }*/
	  Float_t RecoilDeltaPhi = 0.;
	  if (fJetSelection == kRecoil){
	    RecoilDeltaPhi = RelativePhi(TriggerHadron->Phi(), Jet1->Phi());
	    if (TMath::Abs(RecoilDeltaPhi) < (TMath::Pi() - fRecoilAngularWindow)) continue;  //accept the jet only if it overlaps with the recoil phi area of the trigger
	    fh2PtTriggerHadronJet->Fill(TriggerHadron->Pt(), Jet1->Pt());
	    fhPhiTriggerHadronJet->Fill(RelativePhi(TriggerHadron->Phi(), Jet1->Phi()));
	  }
   	  if (!EventCounter){
	    fhEventCounter->Fill(3);
	    EventCounter=kTRUE;
	  }
	  if(fJetShapeSub==kConstSub){
	    JetNumber=-1;
	    for(Int_t i = 0; i<JetCont2->GetNJets(); i++) {
	      Jet2 = JetCont2->GetJet(i);
	      if(Jet2->GetLabel()==Jet1->GetLabel()) {
		JetNumber=i;
	      }
	    }
	    if(JetNumber==-1) continue;
	    Jet2=JetCont2->GetJet(JetNumber);
	    if (JetCont2->AliJetContainer::GetFractionSharedPt(Jet2)<fSharedFractionPtMin) continue;
	    Jet3=Jet2->ClosestJet();
	  }
	  if(!(fJetShapeSub==kConstSub)){
	    if (JetCont1->AliJetContainer::GetFractionSharedPt(Jet1)<fSharedFractionPtMin) continue;
	    Jet3 = Jet1->ClosestJet();   //Note for NoSub and Deriv Sub cases you must fill both the Unsubtracted and Subtracted Hybrid jet containers with the same jet branch
	  }
	  if (!Jet3) continue;
	  Jet4=Jet3->ClosestJet();
	  if(!Jet4) continue;
	  JetsMatched=kTRUE;

	  std::vector <Double_t> Jet_Constituents_Emb;
	  std::vector <Double_t> Jet_Constituents_Truth;
	  AliVParticle *JetConstituent_Emb=NULL;
	  AliVParticle *JetConstituent_Truth=NULL;
	  if (fMLOn==1){
	    for (Int_t Tracks_Emb=0; Tracks_Emb<Jet1->GetNumberOfTracks(); Tracks_Emb++){
	      JetConstituent_Emb = static_cast<AliVParticle*>(Jet1->TrackAt(Tracks_Emb, JetCont1->GetParticleContainer()->GetArray()));
	      if(!JetConstituent_Emb) continue;
	      Jet_Constituents_Emb.push_back(JetConstituent_Emb->Px());
	      Jet_Constituents_Emb.push_back(JetConstituent_Emb->Py());
	      Jet_Constituents_Emb.push_back(JetConstituent_Emb->Pz());
	      Jet_Constituents_Emb.push_back(JetConstituent_Emb->E());
	    }
	    for (Int_t Tracks_Truth=0; Tracks_Truth<Jet4->GetNumberOfTracks(); Tracks_Truth++){
	      JetConstituent_Truth = static_cast<AliVParticle*>(Jet4->TrackAt(Tracks_Truth, JetCont4->GetParticleContainer()->GetArray()));
	      if(!JetConstituent_Truth) continue;
	      Jet_Constituents_Truth.push_back(JetConstituent_Truth->Px());
	      Jet_Constituents_Truth.push_back(JetConstituent_Truth->Py());
	      Jet_Constituents_Truth.push_back(JetConstituent_Truth->Pz());
	      Jet_Constituents_Truth.push_back(JetConstituent_Truth->E());
	    }
	  }

	  if (fMLOn==1) fShapesVar_Tracks_Rec.push_back(Jet_Constituents_Emb);
	  if (fMLOn==1)fShapesVar_Tracks_Truth.push_back(Jet_Constituents_Truth);
	  
	  ///////////So at the moment tree is only filled when we have matched....is this ok?
	  if (fJetShapeSub==kConstSub) fShapesVar[0]=Jet1->Pt();
	  else fShapesVar[0]=Jet1->Pt()-(GetRhoVal(0)*Jet1->Area());
	  if (fMLOn==0) fShapesVar[2]=FjNSubJettiness(Jet1,0,1,0,1,0);
	  else fShapesVar[2]=Jet1->Eta();
	  if (fMLOn==0) fShapesVar[4]=FjNSubJettiness(Jet1,0,2,0,1,0);
	  else fShapesVar[4]=Jet1->Phi();
	  if (fMLOn==0) fShapesVar[6]=FjNSubJettiness(Jet1,0,2,0,1,8);
	  else fShapesVar[6]=Jet1->M();
	  if (fMLOn==0) fShapesVar[8]=FjNSubJettiness(Jet1,0,2,0,1,1);
	  else fShapesVar[8]=Jet1->GetNumberOfTracks();
	  if (fMLOn==0) fShapesVar[10]=FjNSubJettiness(Jet1,0,2,0,1,9);
	  else fShapesVar[10]=Parton_Info->GetPartonFlag6();
	  if (fMLOn==0) fShapesVar[12]=FjNSubJettiness(Jet1,0,2,0,1,3,fBeta_SD,fZCut);
	  else fShapesVar[12]=Parton_Info->GetPartonEta6();
	  if (fMLOn==0) fShapesVar[14]=FjNSubJettiness(Jet1,0,2,0,1,5,fBeta_SD,fZCut);
	  else fShapesVar[14]=Parton_Info->GetPartonPhi6();
	  if (fMLOn==0) fShapesVar[16]=Jet1->GetLeadingTrack(JetCont1->GetParticleContainer()->GetArray())->Pt();
	  else fShapesVar[16]=Angularity(Jet1,0);
	  if (fMLOn==0) fShapesVar[18]=RelativePhiEventPlane(fEPV0,Jet1->Phi());
	  else fShapesVar[18]=PTD(Jet1,0);
	  //fShapesVar[20]=RelativePhiEventPlane(((AliVAODHeader*)(InputEvent()->GetHeader()))->GetEventplane(),Jet1->Phi());
	  if (fMLOn==0) fShapesVar[20]=FjNSubJettiness(Jet1,0,2,0,1,6,fBeta_SD,fZCut);
	  else fShapesVar[20]=0.0;
	  if (fFullTree){
	    if (fMLOn==0) fShapesVar[22]=FjNSubJettiness(Jet1,0,2,0,1,2);
	    else fShapesVar[22]=0.0;
	    Reclusterer1 = Recluster(Jet1, 0, fSubJetRadius, fSubJetMinPt, fSubJetAlgorithm, "SubJetFinder_1");
	    if (fMLOn==0) fShapesVar[24]=SubJetFraction(Jet1, Reclusterer1, 1, 0, kTRUE, kFALSE);
	    else fShapesVar[24]=0.0;
	    if (fMLOn==0) fShapesVar[26]=SubJetFraction(Jet1, Reclusterer1, 2, 0, kTRUE, kFALSE);
	    else fShapesVar[26]=0.0;
	  }
	  if (JetsMatched){ //even needed? Not now but might be if you want to fill trees when jets aren't matched too
	    fShapesVar[1]=Jet4->Pt();
	    if (fMLOn==0) fShapesVar[3]=FjNSubJettiness(Jet4,3,1,0,1,0);
	    else fShapesVar[3]=Jet4->Eta();
	    if (fMLOn==0) fShapesVar[5]=FjNSubJettiness(Jet4,3,2,0,1,0);
	    else fShapesVar[5]=Jet4->Phi();
	    if (fMLOn==0) fShapesVar[7]=FjNSubJettiness(Jet4,3,2,0,1,8);
	    else fShapesVar[7]=Jet4->M();
	    if (fMLOn==0) fShapesVar[9]=FjNSubJettiness(Jet4,3,2,0,1,1);
	    else fShapesVar[9]=Jet4->GetNumberOfTracks();
	    if (fMLOn==0) fShapesVar[11]=FjNSubJettiness(Jet4,3,2,0,1,9);
	    else fShapesVar[11]=Parton_Info->GetPartonFlag7();
	    if (fMLOn==0) fShapesVar[13]=FjNSubJettiness(Jet4,3,2,0,1,3,fBeta_SD,fZCut);
	    else fShapesVar[13]=Parton_Info->GetPartonEta7();
	    if (fMLOn==0) fShapesVar[15]=FjNSubJettiness(Jet4,3,2,0,1,5,fBeta_SD,fZCut);
	    else fShapesVar[15]=Parton_Info->GetPartonPhi7();
	    if (fMLOn==0) fShapesVar[17]=Jet4->GetLeadingTrack(JetCont4->GetParticleContainer()->GetArray())->Pt();
	    else fShapesVar[17]=Angularity(Jet4,3);
	    if (fMLOn==0) fShapesVar[19]=RelativePhiEventPlane(fEPV0,Jet4->Phi());
	    else fShapesVar[19]=PTD(Jet4,3);
	    //fShapesVar[21]=RelativePhiEventPlane(((AliVAODHeader*)(InputEvent()->GetHeader()))->GetEventplane(),Jet4->Phi());
	    if (fMLOn==0) fShapesVar[21]=FjNSubJettiness(Jet4,3,2,0,1,6,fBeta_SD,fZCut);
	    else fShapesVar[21]=0.0;
	    if (fFullTree){
	      if (fMLOn==0) fShapesVar[23]=FjNSubJettiness(Jet4,3,2,0,1,2);
	      else fShapesVar[23]=0.0;
	      Reclusterer4=Recluster(Jet4, 3, fSubJetRadius, 0, fSubJetAlgorithm, "SubJetFinder_4");
	      if (fMLOn==0) fShapesVar[25]=SubJetFraction(Jet4, Reclusterer4, 1, 0, kTRUE, kFALSE);
	      else fShapesVar[25]=0.0;
	      if (fMLOn==0) fShapesVar[27]=SubJetFraction(Jet4, Reclusterer4, 2, 0, kTRUE, kFALSE);
	      else fShapesVar[27]=0.0;
	    } 
	  }
	  else{
	    fShapesVar[1]=-2;
	    fShapesVar[3]=-2;
	    fShapesVar[5]=-2;
	    fShapesVar[7]=-2;
	    fShapesVar[9]=-2;
	    fShapesVar[11]=-2;
	    fShapesVar[13]=-2;
	    fShapesVar[15]=-2;
	    fShapesVar[17]=-2;
	    fShapesVar[19]=-2;
	    fShapesVar[21]=-2;
	    if (fFullTree){
	      fShapesVar[23]=-2;
	      fShapesVar[25]=-2;
	      fShapesVar[27]=-2;
	    }
	  }
	  fTreeResponseMatrixAxis->Fill();
	  if (fMLOn==1) fTreeTracks->Fill();
	  JetsMatched=kFALSE;
	}
      }
    }
  }
  
////////////////////////////////////MC Det vs Part Level////////////////////////////////////////
  if (fJetShapeType == AliAnalysisTaskSubJetFraction::kTrueDet){	        
    AliEmcalJet *Jet1 = NULL; //Detector Level Jet                                                                                                                  
    AliEmcalJet *Jet2 = NULL; //Particle Level Jet  
    AliJetContainer *JetCont1= GetJetContainer(0); //Jet Container for Detector Level Pythia
    AliJetContainer *JetCont2= GetJetContainer(1); //Jet Container for Particle Level Pythia
    AliEmcalJetFinder *Reclusterer1; //Object containg Subjets Detector Level 
    AliEmcalJetFinder *Reclusterer2; //Object containg Subjets Particle Level 
    Int_t JetCounter1=0; //Counts number of jets in event
    Int_t JetCounter2=0; //Counts number of jets in event                                                                                                      
    Double_t JetPhi1=0;
    Double_t JetPhi2=0;
    Bool_t JetsMatched=kFALSE;
    Double_t Pythia_Event_Weight=1;
    Bool_t EventCounter=kFALSE;
    const AliEmcalPythiaInfo *Parton_Info = 0x0;
    Parton_Info=GetPythiaInfo();
    fhEventCounter_1->Fill(1);
    if(JetCont1) {
      fhEventCounter_1->Fill(2); //Number of events with a jet container                                                                                               
      JetCont1->ResetCurrentID();
      while((Jet1=JetCont1->GetNextAcceptJet())) {
    std::vector <Double_t> Jet_Constituents_Det;
    std::vector <Double_t> Jet_Constituents_Truth;
	if( (!Jet1) || ((Jet1->Pt())<fPtThreshold)) {
	  // fhEventCounter_1->Fill(3); //events where the jet had a problem                                                                                   
	  continue;
	}
	else {
	  /* if(fSemigoodCorrect){
	    Double_t HoleDistance=RelativePhi(Jet1->Phi(),fHolePos);
	    if(TMath::Abs(HoleDistance)<fHoleWidth) continue;
	    }*/
	  Float_t RecoilDeltaPhi = 0.;
	  if (fJetSelection == kRecoil){
	    RecoilDeltaPhi = RelativePhi(TriggerHadron->Phi(), Jet1->Phi());
	    if (TMath::Abs(RecoilDeltaPhi) < (TMath::Pi() - fRecoilAngularWindow)) continue;  //accept the jet only if it overlaps with the recoil phi area of the trigger
	    fh2PtTriggerHadronJet->Fill(TriggerHadron->Pt(), Jet1->Pt());
	    fhPhiTriggerHadronJet->Fill(RelativePhi(TriggerHadron->Phi(), Jet1->Phi()));
	  }
	  if (!EventCounter){
	    fhEventCounter_1->Fill(3);
	    EventCounter=kTRUE;
	  }
	  if((Jet1->GetNumberOfTracks())==0){
	    fhEventCounter_1->Fill(10); //zero track jets                                                                                                         
	  }
	  if((Jet1->GetNumberOfTracks())==1){
	    fhEventCounter_1->Fill(11); //one track jets                                                                                                              
	  }
	  fhEventCounter_1->Fill(4); //Number of Jets found in all events                                                                                         
	  JetCounter1++;
	  fhJetPt_1->Fill(Jet1->Pt());
	  JetPhi1=Jet1->Phi();
	  if(JetPhi1 < -1*TMath::Pi()) JetPhi1 += (2*TMath::Pi());
	  else if (JetPhi1 > TMath::Pi()) JetPhi1 -= (2*TMath::Pi());
	  fhJetPhi_1->Fill(JetPhi1);
	  fhJetEta_1->Fill(Jet1->Eta());
	  fhJetMass_1->Fill(Jet1->M());
	  fhJetRadius_1->Fill(TMath::Sqrt((Jet1->Area()/TMath::Pi()))); //Radius of Jets per event                                                            
	  fhNumberOfJetTracks_1->Fill(Jet1->GetNumberOfTracks());
	  if((Jet2 = Jet1->ClosestJet())){
	    JetsMatched=kTRUE;
	    if((Jet2->GetNumberOfTracks())==0){
	      fhEventCounter_2->Fill(10); //zero track jets                                                                                              
	    }
	    if((Jet2->GetNumberOfTracks())==1){
	      fhEventCounter_2->Fill(11); //one track jets                                                                                              
	    }
	    fhEventCounter_2->Fill(4); //Number of Jets found in all events                                                                                                  
	    JetCounter2++;
	    fhJetPt_2->Fill(Jet2->Pt());
	    JetPhi2=Jet2->Phi();
	    if(JetPhi2 < -1*TMath::Pi()) JetPhi2 += (2*TMath::Pi());
	    else if (JetPhi2 > TMath::Pi()) JetPhi2 -= (2*TMath::Pi());
	    fhJetPhi_2->Fill(JetPhi2);
	    fhJetEta_2->Fill(Jet2->Eta());
	    fhJetMass_2->Fill(Jet2->M());
	    fhJetRadius_2->Fill(TMath::Sqrt((Jet2->Area()/TMath::Pi()))); //Radius of Jets per event                                                                         
	    fhNumberOfJetTracks_2->Fill(Jet2->GetNumberOfTracks());
	    fh2PtRatio->Fill(Jet1->Pt(),Jet2->Pt(),Pythia_Event_Weight);

	    
	    AliVParticle *JetConstituent_Det=NULL;
	    AliVParticle *JetConstituent_Truth=NULL;
	    if (fMLOn==1){
	      for (Int_t Tracks_Det=0; Tracks_Det<Jet1->GetNumberOfTracks(); Tracks_Det++){
		JetConstituent_Det = static_cast<AliVParticle*>(Jet1->TrackAt(Tracks_Det, JetCont1->GetParticleContainer()->GetArray()));
		if(!JetConstituent_Det) continue;
		Jet_Constituents_Det.push_back(JetConstituent_Det->Px());
		Jet_Constituents_Det.push_back(JetConstituent_Det->Py());
		Jet_Constituents_Det.push_back(JetConstituent_Det->Pz());
		Jet_Constituents_Det.push_back(JetConstituent_Det->E());
	      }
	      for (Int_t Tracks_Truth=0; Tracks_Truth<Jet2->GetNumberOfTracks(); Tracks_Truth++){
		JetConstituent_Truth = static_cast<AliVParticle*>(Jet2->TrackAt(Tracks_Truth, JetCont2->GetParticleContainer()->GetArray()));
		if(!JetConstituent_Truth) continue;
		Jet_Constituents_Truth.push_back(JetConstituent_Truth->Px());
		Jet_Constituents_Truth.push_back(JetConstituent_Truth->Py());
		Jet_Constituents_Truth.push_back(JetConstituent_Truth->Pz());
		Jet_Constituents_Truth.push_back(JetConstituent_Truth->E());
	      }
	    }
	    
	  }
          else {
            fhEventCounter_2->Fill(1);
	    //continue;
          }
	  /////////////////How do you do this???? Do you only fill them for when both sets work or only for when one works??????///////////////////////
	  
	  if (fMLOn==1) fShapesVar_Tracks_Rec.push_back(Jet_Constituents_Det);
	  if (fMLOn==1)fShapesVar_Tracks_Truth.push_back(Jet_Constituents_Truth);
	  
	  fShapesVar[0]=Jet1->Pt();
	  if (fMLOn==0) fShapesVar[2]=FjNSubJettiness(Jet1,0,1,0,1,0);
	  else fShapesVar[2]=Jet1->Eta();
	  if (fMLOn==0) fShapesVar[4]=FjNSubJettiness(Jet1,0,2,0,1,0);
	  else fShapesVar[4]=Jet1->Phi();
	  if (fMLOn==0) fShapesVar[6]=FjNSubJettiness(Jet1,0,2,0,1,8);
	  else fShapesVar[6]=Jet1->M();
	  if (fMLOn==0) fShapesVar[8]=FjNSubJettiness(Jet1,0,2,0,1,1);
	  else fShapesVar[8]=Jet1->GetNumberOfTracks();
	  if (fMLOn==0)  fShapesVar[10]=FjNSubJettiness(Jet1,0,2,0,1,9);
	  else fShapesVar[10]=Parton_Info->GetPartonFlag6();
	  if (fMLOn==0) fShapesVar[12]=FjNSubJettiness(Jet1,0,2,0,1,3,fBeta_SD,fZCut);
	  else fShapesVar[12]=Parton_Info->GetPartonEta6();
	  if (fMLOn==0) fShapesVar[14]=FjNSubJettiness(Jet1,0,2,0,1,5,fBeta_SD,fZCut);
	  else fShapesVar[14]=Parton_Info->GetPartonPhi6();
	  if (fMLOn==0) fShapesVar[16]=Jet1->GetLeadingTrack(JetCont1->GetParticleContainer()->GetArray())->Pt();
	  else fShapesVar[16]=Angularity(Jet1,0);
	  if (fMLOn==0) fShapesVar[18]=-2; //event plane calculation only needed for PbPb recoils
	  else fShapesVar[18]=PTD(Jet1,0);
	  if (fMLOn==0) fShapesVar[20]=FjNSubJettiness(Jet1,0,2,0,1,6,fBeta_SD,fZCut);
	  else fShapesVar[20]=0.0;
	  Reclusterer1 = Recluster(Jet1, 0, fSubJetRadius, fSubJetMinPt, fSubJetAlgorithm, "SubJetFinder_1");
	  if (fFullTree){
	    if (fMLOn==0)  fShapesVar[22]=FjNSubJettiness(Jet1,0,2,0,1,2);
	    else fShapesVar[22]=0.0;
	    if (fMLOn==0) fShapesVar[24]=SubJetFraction(Jet1, Reclusterer1, 1, 0, kTRUE, kFALSE);
	    else fShapesVar[24]=0.0;
	    if (fMLOn==0) fShapesVar[26]=SubJetFraction(Jet1, Reclusterer1, 2, 0, kTRUE, kFALSE);
	    else fShapesVar[26]=0.0;
	  }
	  if (JetsMatched){ //even needed? Not now but might be if you want to fill trees when jets aren't matched too
	    fShapesVar[1]=Jet2->Pt();
	    if (fMLOn==0) fShapesVar[3]=FjNSubJettiness(Jet2,1,1,0,1,0);
	    else fShapesVar[3]=Jet2->Eta();
	    if (fMLOn==0) fShapesVar[5]=FjNSubJettiness(Jet2,1,2,0,1,0);
	    else fShapesVar[5]=Jet2->Phi();
	    if (fMLOn==0) fShapesVar[7]=FjNSubJettiness(Jet2,1,2,0,1,8);
	    else fShapesVar[7]=Jet2->M();
	    if (fMLOn==0) fShapesVar[9]=FjNSubJettiness(Jet2,1,2,0,1,1);
	    else fShapesVar[9]=Jet2->GetNumberOfTracks();
	    if (fMLOn==0) fShapesVar[11]=FjNSubJettiness(Jet2,1,2,0,1,9);
	    else fShapesVar[11]=Parton_Info->GetPartonFlag7();
	    if (fMLOn==0) fShapesVar[13]=FjNSubJettiness(Jet2,1,2,0,1,3,fBeta_SD,fZCut);
	    else fShapesVar[13]=Parton_Info->GetPartonEta7();
	    if (fMLOn==0) fShapesVar[15]=FjNSubJettiness(Jet2,1,2,0,1,5,fBeta_SD,fZCut);
	    else fShapesVar[15]=Parton_Info->GetPartonPhi7();
	    if (fMLOn==0) fShapesVar[17]=Jet2->GetLeadingTrack(JetCont2->GetParticleContainer()->GetArray())->Pt();
	    else fShapesVar[17]=Angularity(Jet2,1);
	    if (fMLOn==0) fShapesVar[19]=-2;
	    else fShapesVar[19]=PTD(Jet2,1);
	    if (fMLOn==0) fShapesVar[21]=FjNSubJettiness(Jet2,1,2,0,1,6,fBeta_SD,fZCut);
	    else fShapesVar[21]=0.0;
	    Reclusterer2 = Recluster(Jet2, 1, fSubJetRadius, 0, fSubJetAlgorithm, "SubJetFinder_2");
	    if (fFullTree){
	      if (fMLOn==0) fShapesVar[23]=FjNSubJettiness(Jet2,1,2,0,1,2);
	      else fShapesVar[23]=0.0;
	      if (fMLOn==0) fShapesVar[25]=SubJetFraction(Jet2, Reclusterer2, 1, 0, kTRUE, kFALSE);
	      else fShapesVar[25]=0.0;
	      if (fMLOn==0) fShapesVar[27]=SubJetFraction(Jet2, Reclusterer2, 2, 0, kTRUE, kFALSE);
	      else fShapesVar[27]=0.0;
	    } 
	  }
	  else{
	    fShapesVar[1]=-2;
	    fShapesVar[3]=-2;
	    fShapesVar[5]=-2;
	    fShapesVar[7]=-2;
	    fShapesVar[9]=-2;
	    fShapesVar[11]=-2;
	    fShapesVar[13]=-2;
	    fShapesVar[15]=-2;
	    fShapesVar[17]=-2;
	    fShapesVar[19]=-2;
	    fShapesVar[21]=-2;
	    if (fFullTree){
	      fShapesVar[23]=-2;
	      fShapesVar[25]=-2;
	      fShapesVar[27]=-2;
	    }
	  }
	  fTreeResponseMatrixAxis->Fill();
	  if (fMLOn==1)fTreeTracks->Fill();
	  
	  fhSubJetCounter_1->Fill(Reclusterer1->GetNumberOfJets());
	  for (Int_t i= 0; i<Reclusterer1->GetNumberOfJets(); i++){
	    fhEventCounter_1->Fill(6); //Number of overall subjets in all events                                                                                           
	    fhSubJetPt_1->Fill(Reclusterer1->GetJet(i)->Pt());
	    fhSubJetMass_1->Fill(Reclusterer1->GetJet(i)->M());
	    fhNumberOfSubJetTracks_1->Fill(Reclusterer1->GetJet(i)->GetNumberOfTracks());
	    fhSubJetRadius_1->Fill(TMath::Sqrt((Reclusterer1->GetJet(i)->Area()/TMath::Pi()))); //Radius of SubJets per event                                                 
	  }
	  if(JetsMatched){
	    fhSubJetCounter_2->Fill(Reclusterer2->GetNumberOfJets());
	    for (Int_t i= 0; i<Reclusterer2->GetNumberOfJets(); i++){
	      fhEventCounter_2->Fill(6); //Number of overall subjets in all events                                                                                          
	      fhSubJetPt_2->Fill(Reclusterer2->GetJet(i)->Pt());
	      fhSubJetMass_2->Fill(Reclusterer2->GetJet(i)->M());
	      fhNumberOfSubJetTracks_2->Fill(Reclusterer2->GetJet(i)->GetNumberOfTracks());
	      fhSubJetRadius_2->Fill(TMath::Sqrt((Reclusterer2->GetJet(i)->Area()/TMath::Pi()))); //Radius of SubJets per event                                            
	    }
	  }
	  JetsMatched=kFALSE;
	}
      }
      fhJetCounter_1->Fill(JetCounter1); //Number of Jets in Each Event Particle Level                                                                 
      fhJetCounter_2->Fill(JetCounter2); //Number of Jets in Each Event Detector Level        
    }
    //else {fhEventCounter_1->Fill(3);} //Events with no jet container 
  }


  
  if (fJetShapeType == AliAnalysisTaskSubJetFraction::kData){
    AliEmcalJet *Jet1 = NULL; //Original Jet in the event                                                                                                                
    AliJetContainer *JetCont= GetJetContainer(0); //Jet Container for event 
    Int_t JetCounter=0; //Counts number of jets in event  
    Double_t JetPhi=0;
    Bool_t EventCounter=kFALSE;
    Double_t JetPt_ForThreshold=0;
    fhEventCounter->Fill(1);
    if(JetCont) {
      fhEventCounter->Fill(2); //Number of events with a jet container
      JetCont->ResetCurrentID();
      while((Jet1=JetCont->GetNextAcceptJet())) {
	if(!Jet1) continue;
	if (fJetShapeSub==kNoSub || fJetShapeSub==kDerivSub) JetPt_ForThreshold = Jet1->Pt()-(GetRhoVal(0)*Jet1->Area());
	else JetPt_ForThreshold = Jet1->Pt();
	if(JetPt_ForThreshold<fPtThreshold) {
	  //fhEventCounter->Fill(3); //events where the jet had a problem
	  continue;
	}
	else {
	  /* if(fSemigoodCorrect){
	    Double_t HoleDistance=RelativePhi(Jet1->Phi(),fHolePos);
	    if(TMath::Abs(HoleDistance)<fHoleWidth) continue;
	  }*/
	  Float_t RecoilDeltaPhi = 0.;
	  if (fJetSelection == kRecoil){
	    RecoilDeltaPhi = RelativePhi(TriggerHadron->Phi(), Jet1->Phi());
	    if (TMath::Abs(RecoilDeltaPhi) < (TMath::Pi() - fRecoilAngularWindow)) continue;  //accept the jet only if it overlaps with the recoil phi area of the trigger
	    fh2PtTriggerHadronJet->Fill(TriggerHadron->Pt(), Jet1->Pt());
	    fhPhiTriggerHadronJet->Fill(RelativePhi(TriggerHadron->Phi(), Jet1->Phi()));
	  }
          if (!EventCounter){
            fhEventCounter->Fill(3);
            EventCounter=kTRUE;
          }
	  if((Jet1->GetNumberOfTracks())==0){
	    fhEventCounter->Fill(10); //zero track jets
	  }
	  if((Jet1->GetNumberOfTracks())==1){
	    fhEventCounter->Fill(11); //one track jets
	  }
          fhEventCounter->Fill(4); //Number of Jets found in all events  
	  JetCounter++;
	  fhJetPt->Fill(Jet1->Pt());    
	  JetPhi=Jet1->Phi();
	  //if(JetPhi < -1*TMath::Pi()) JetPhi += (2*TMath::Pi());
	  //else if (JetPhi > TMath::Pi()) JetPhi -= (2*TMath::Pi());
	  fhJetPhi->Fill(JetPhi);
	  fhJetEta->Fill(Jet1->Eta());
	  fhJetMass->Fill(Jet1->M());
	  fhJetRadius->Fill(TMath::Sqrt((Jet1->Area()/TMath::Pi()))); //Radius of Jets per event
          fhNumberOfJetTracks->Fill(Jet1->GetNumberOfTracks());
	  std::vector <Double_t> Jet_Constituents_Data;
	  AliVParticle *JetConstituent_Data=NULL;
	  AliVParticle *JetConstituent_Truth=NULL;
	  if (fMLOn==1){
	    for (Int_t Tracks_Data=0; Tracks_Data<Jet1->GetNumberOfTracks(); Tracks_Data++){
	      JetConstituent_Data = static_cast<AliVParticle*>(Jet1->TrackAt(Tracks_Data, JetCont->GetParticleContainer()->GetArray()));
	      if(!JetConstituent_Data) continue;
	      Jet_Constituents_Data.push_back(JetConstituent_Data->Px());
	      Jet_Constituents_Data.push_back(JetConstituent_Data->Py());
	      Jet_Constituents_Data.push_back(JetConstituent_Data->Pz());
	      Jet_Constituents_Data.push_back(JetConstituent_Data->E());
	    }
	    fShapesVar_Tracks_Rec.push_back(Jet_Constituents_Data);
	    fShapesVar_Tracks_Truth.push_back(Jet_Constituents_Data);
	  }
	  if(fJetShapeSub==kNoSub || fJetShapeSub==kDerivSub) fShapesVar[0]= Jet1->Pt()-(GetRhoVal(0)*Jet1->Area());
	  else fShapesVar[0]=Jet1->Pt();
	  if (fMLOn==0) fShapesVar[2]=FjNSubJettiness(Jet1,0,1,0,1,0);
	  else fShapesVar[2]=Jet1->Eta();
	  if (fMLOn==0) fShapesVar[4]=FjNSubJettiness(Jet1,0,2,0,1,0);
	  else fShapesVar[4]=Jet1->Phi();
	  if (fMLOn==0) fShapesVar[6]=FjNSubJettiness(Jet1,0,2,0,1,8);
	  else fShapesVar[6]=Jet1->M();
	  if (fMLOn==0) fShapesVar[8]=FjNSubJettiness(Jet1,0,2,0,1,1);
	  else fShapesVar[8]=Jet1->GetNumberOfTracks();
	  if (fMLOn==0) fShapesVar[10]=FjNSubJettiness(Jet1,0,2,0,1,9);
	  else fShapesVar[10]=0.0;
	  if (fMLOn==0) fShapesVar[12]=FjNSubJettiness(Jet1,0,2,0,1,3,fBeta_SD,fZCut);
	  else fShapesVar[12]=0.0;
	  if (fMLOn==0) fShapesVar[14]=FjNSubJettiness(Jet1,0,2,0,1,5,fBeta_SD,fZCut);
	  else fShapesVar[14]=0.0;
	  if (fMLOn==0) fShapesVar[16]=Jet1->GetLeadingTrack(JetCont->GetParticleContainer()->GetArray())->Pt();
	  else fShapesVar[16]=Angularity(Jet1,0);
	  if (fMLOn==0) fShapesVar[18]=-2; //event plane calculation not needed for data
	  else fShapesVar[18]=PTD(Jet1,0);
	  if (fMLOn==0) fShapesVar[20]=FjNSubJettiness(Jet1,0,2,0,1,6,fBeta_SD,fZCut);
	  else fShapesVar[20]=0.0;
	  AliEmcalJetFinder *Reclusterer1 = Recluster(Jet1, 0, fSubJetRadius, fSubJetMinPt, fSubJetAlgorithm, "SubJetFinder");
	  if (fFullTree){
	    if (fMLOn==0) fShapesVar[22]=FjNSubJettiness(Jet1,0,2,0,1,2);
	    else fShapesVar[22]=0.0;
	    if (fMLOn==0) fShapesVar[24]=SubJetFraction(Jet1, Reclusterer1, 1, 0, kTRUE, kFALSE);
	    else fShapesVar[24]=0.0;
	    if (fMLOn==0) fShapesVar[26]=SubJetFraction(Jet1, Reclusterer1, 2, 0, kTRUE, kFALSE);
	    else fShapesVar[26]=0.0;
	  }
	  fShapesVar[1]=-2;
	  fShapesVar[3]=-2;
	  fShapesVar[5]=-2;
	  fShapesVar[7]=-2;
	  fShapesVar[9]=-2;
	  fShapesVar[11]=-2;
	  fShapesVar[13]=-2;
	  fShapesVar[15]=-2;
	  fShapesVar[17]=-2;
	  fShapesVar[19]=-2;
	  fShapesVar[21]=-2;
	  if (fFullTree){
	    fShapesVar[23]=-2;
	    fShapesVar[25]=-2;
	    fShapesVar[27]=-2;
	  }
	  fTreeResponseMatrixAxis->Fill();
	  if (fMLOn==1) fTreeTracks->Fill();
	  fhSubJetCounter->Fill(Reclusterer1->GetNumberOfJets());
	  for (Int_t i= 0; i<Reclusterer1->GetNumberOfJets(); i++){
	    fhEventCounter->Fill(6); //Number of overall subjets in all events
	    fhSubJetPt->Fill(Reclusterer1->GetJet(i)->Pt());
	    fhSubJetMass->Fill(Reclusterer1->GetJet(i)->M());
	    fhNumberOfSubJetTracks->Fill(Reclusterer1->GetJet(i)->GetNumberOfTracks());
	    fhSubJetRadius->Fill(TMath::Sqrt((Reclusterer1->GetJet(i)->Area()/TMath::Pi()))); //Radius of SubJets per event
	  }
	}
      }
      fhJetCounter->Fill(JetCounter); //Number of Jets in Each Event
    }
    //else {fhEventCounter->Fill(2);} //Events with no jet container
  }

  
  if (fJetShapeType==AliAnalysisTaskSubJetFraction::kSim || fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly){
    
    AliEmcalJet *Jet1 = NULL; //Original Jet in the event                                                                                                                
    AliJetContainer *JetCont= GetJetContainer(0); //Jet Container for event 
    Int_t JetCounter=0; //Counts number of jets in event  
    Double_t JetPhi=0;
    Bool_t EventCounter=kFALSE;
    Double_t JetPt_ForThreshold=0;
    const AliEmcalPythiaInfo *Parton_Info = 0x0;
    Parton_Info=GetPythiaInfo();
    fhEventCounter->Fill(1);
    if(JetCont) {
      fhEventCounter->Fill(2); //Number of events with a jet container
      JetCont->ResetCurrentID();
      while((Jet1=JetCont->GetNextAcceptJet())) {
	//	if (((Jet1=JetCont->GetLeadingJet("rho")) && (fPtThreshold<=Jet1->Pt()-(GetRhoVal(0)*Jet1->Area())))){ 
	if(!Jet1) continue;
	if(fJetShapeSub==kNoSub || fJetShapeSub==kDerivSub) JetPt_ForThreshold = Jet1->Pt()-(GetRhoVal(0)*Jet1->Area());
	else JetPt_ForThreshold = Jet1->Pt(); 
	if(JetPt_ForThreshold<fPtThreshold) {
	  //fhEventCounter->Fill(3); //events where the jet had a problem
	  continue;
	}
	else {
	  /* if(fSemigoodCorrect){
	    Double_t HoleDistance=RelativePhi(Jet1->Phi(),fHolePos);
	    if(TMath::Abs(HoleDistance)<fHoleWidth) continue;
	  }*/
	  Float_t RecoilDeltaPhi = 0.;
	  if (fJetSelection == kRecoil){
	    RecoilDeltaPhi = RelativePhi(TriggerHadron->Phi(), Jet1->Phi());
	    if (TMath::Abs(RecoilDeltaPhi) < (TMath::Pi() - fRecoilAngularWindow)) continue;  //accept the jet only if it overlaps with the recoil phi area of the trigger
	    fh2PtTriggerHadronJet->Fill(TriggerHadron->Pt(), Jet1->Pt());
	    fhPhiTriggerHadronJet->Fill(RelativePhi(TriggerHadron->Phi(), Jet1->Phi()));
	  }
          if (!EventCounter){
            fhEventCounter->Fill(3);
            EventCounter=kTRUE;
          }
	  if((Jet1->GetNumberOfTracks())==0){
	    fhEventCounter->Fill(10); //zero track jets
	  }
	  if((Jet1->GetNumberOfTracks())==1){
	    fhEventCounter->Fill(11); //one track jets
	  }
          fhEventCounter->Fill(4); //Number of Jets found in all events  
	  JetCounter++;
	  fhJetPt->Fill(Jet1->Pt());    
	  JetPhi=Jet1->Phi();
	  if(JetPhi < -1*TMath::Pi()) JetPhi += (2*TMath::Pi());
	  else if (JetPhi > TMath::Pi()) JetPhi -= (2*TMath::Pi());
	  fhJetPhi->Fill(JetPhi);
	  fhJetEta->Fill(Jet1->Eta());
	  fhJetMass->Fill(Jet1->M());
	  fhJetRadius->Fill(TMath::Sqrt((Jet1->Area()/TMath::Pi()))); //Radius of Jets per event
          fhNumberOfJetTracks->Fill(Jet1->GetNumberOfTracks());
	  std::pair<fastjet::PseudoJet,fastjet::ClusterSequence *> Randomised_Jet=ModifyJet(Jet1,0,"Randomise");
	  std::pair<fastjet::PseudoJet,fastjet::ClusterSequence *> ExtraProng_Jet_02_15=ModifyJet(Jet1,0,"AddExtraProng_02_15");
	  std::pair<fastjet::PseudoJet,fastjet::ClusterSequence *> ExtraProng_Jet_02_30=ModifyJet(Jet1,0,"AddExtraProng_02_30");
	  std::pair<fastjet::PseudoJet,fastjet::ClusterSequence *> ExtraProng_Jet_02_45=ModifyJet(Jet1,0,"AddExtraProng_02_45");
	  std::pair<fastjet::PseudoJet,fastjet::ClusterSequence *> ExtraProng_Jet_02_60=ModifyJet(Jet1,0,"AddExtraProng_02_60");
	  std::pair<fastjet::PseudoJet,fastjet::ClusterSequence *> ExtraProng_Jet_01_10=ModifyJet(Jet1,0,"AddExtraProng_01_10");
	  std::pair<fastjet::PseudoJet,fastjet::ClusterSequence *> ExtraProng_Jet_01_15=ModifyJet(Jet1,0,"AddExtraProng_01_15");
	  std::pair<fastjet::PseudoJet,fastjet::ClusterSequence *> ExtraProng_Jet_01_20=ModifyJet(Jet1,0,"AddExtraProng_01_20");
	  std::pair<fastjet::PseudoJet,fastjet::ClusterSequence *> ExtraProng_Jet_01_25=ModifyJet(Jet1,0,"AddExtraProng_01_25");
	  std::pair<fastjet::PseudoJet,fastjet::ClusterSequence *> ExtraProng_Jet_01_30=ModifyJet(Jet1,0,"AddExtraProng_01_30");
	  std::pair<fastjet::PseudoJet,fastjet::ClusterSequence *> ExtraProng_Jet_01_45=ModifyJet(Jet1,0,"AddExtraProng_01_45");
	  std::pair<fastjet::PseudoJet,fastjet::ClusterSequence *> ExtraProng_Jet_03_10=ModifyJet(Jet1,0,"AddExtraProng_03_10");
	  std::pair<fastjet::PseudoJet,fastjet::ClusterSequence *> ExtraProng_Jet_03_15=ModifyJet(Jet1,0,"AddExtraProng_03_15");
	  std::pair<fastjet::PseudoJet,fastjet::ClusterSequence *> ExtraProng_Jet_03_20=ModifyJet(Jet1,0,"AddExtraProng_03_20");
	  std::pair<fastjet::PseudoJet,fastjet::ClusterSequence *> ExtraProng_Jet_03_30=ModifyJet(Jet1,0,"AddExtraProng_03_30");
	  std::pair<fastjet::PseudoJet,fastjet::ClusterSequence *> ExtraProng_Jet_03_45=ModifyJet(Jet1,0,"AddExtraProng_03_45");

	  
	  std::pair<fastjet::PseudoJet,fastjet::ClusterSequence *> kTTrack_1_2_1=ModifyJet(Jet1,0,"AddkTTrack_1_2_1");

	  std::vector <Double_t> Jet_Constituents_Data;
	  AliVParticle *JetConstituent_Data=NULL;
	  AliVParticle *JetConstituent_Truth=NULL;
	  if (fMLOn==1){
	    for (Int_t Tracks_Data=0; Tracks_Data<Jet1->GetNumberOfTracks(); Tracks_Data++){
	      JetConstituent_Data = static_cast<AliVParticle*>(Jet1->TrackAt(Tracks_Data, JetCont->GetParticleContainer()->GetArray()));
	      if(!JetConstituent_Data) continue;
	      Jet_Constituents_Data.push_back(JetConstituent_Data->Px());
	      Jet_Constituents_Data.push_back(JetConstituent_Data->Py());
	      Jet_Constituents_Data.push_back(JetConstituent_Data->Pz());
	      Jet_Constituents_Data.push_back(JetConstituent_Data->E());
	    }
	    fShapesVar_Tracks_Rec.push_back(Jet_Constituents_Data);
	    fShapesVar_Tracks_Truth.push_back(Jet_Constituents_Data);
	  }

	  
	  if(fJetShapeSub==kNoSub || fJetShapeSub==kDerivSub) fShapesVar[0]= Jet1->Pt()-(GetRhoVal(0)*Jet1->Area());
	  else fShapesVar[0]=Jet1->Pt(); 
	  if (fMLOn==0) fShapesVar[2]=FjNSubJettiness(Jet1,0,1,0,1,0);
	  else fShapesVar[2]=Jet1->Eta();
	  if (fMLOn==0) fShapesVar[4]=FjNSubJettiness(Jet1,0,2,0,1,0);
	  else fShapesVar[4]=Jet1->Phi();
	  if (fMLOn==0) fShapesVar[6]=FjNSubJettinessFastJet(ExtraProng_Jet_01_20,0,1,0,1,0);
	  else fShapesVar[6]=Jet1->M();
	  if (fMLOn==0) fShapesVar[8]=FjNSubJettiness(Jet1,0,2,0,1,1);
	  else fShapesVar[8]=Jet1->GetNumberOfTracks();
	  if (fMLOn==0) fShapesVar[10]=FjNSubJettinessFastJet(ExtraProng_Jet_03_10,0,1,0,1,0);
	  else fShapesVar[10]=Parton_Info->GetPartonFlag6();
	  if (fMLOn==0) fShapesVar[12]=FjNSubJettinessFastJet(ExtraProng_Jet_03_20,0,1,0,1,0);
	  else fShapesVar[12]=Parton_Info->GetPartonEta6();
	  if (fMLOn==0) fShapesVar[14]=FjNSubJettinessFastJet(ExtraProng_Jet_01_10,0,2,0,1,0);
	  else fShapesVar[14]=Parton_Info->GetPartonPhi6();
	  if (fMLOn==0) fShapesVar[16]=FjNSubJettinessFastJet(ExtraProng_Jet_01_20,0,2,0,1,0);
	  else fShapesVar[16]=Angularity(Jet1,0);
	  // fShapesVar[18]=FjNSubJettinessFastJet(ExtraProng_Jet_01_30,0,2,0,1,0);
	  if (fMLOn==0) fShapesVar[18]=FjNSubJettinessFastJet(kTTrack_1_2_1,0,2,0,1,0);
	  else fShapesVar[18]=PTD(Jet1,0);
	  if (fMLOn==0) fShapesVar[20]=FjNSubJettinessFastJet(ExtraProng_Jet_03_15,0,2,0,1,0);
	  else fShapesVar[20]=0.0;
	  AliEmcalJetFinder *Reclusterer1 = Recluster(Jet1, 0, fSubJetRadius, fSubJetMinPt, fSubJetAlgorithm, "SubJetFinder");
	  if (fFullTree){
	    if (fMLOn==0) fShapesVar[22]=FjNSubJettiness(Jet1,0,2,0,1,2);
	    else fShapesVar[22]=0.0;
	    if (fMLOn==0) fShapesVar[24]=SubJetFraction(Jet1, Reclusterer1, 1, 0, kTRUE, kFALSE);
	    else fShapesVar[24]=0.0;
	    if (fMLOn==0) fShapesVar[26]=SubJetFraction(Jet1, Reclusterer1, 2, 0, kTRUE, kFALSE);
	    else fShapesVar[26]=0.0;
	  }
	  if (fMLOn==0) fShapesVar[1]=FjNSubJettinessFastJet(Randomised_Jet,0,1,0,1,0);
	  else fShapesVar[1]=0.0;
	  // fShapesVar[3]=FjNSubJettiness(std::unique_ptr<AliEmcalJet>(ModifyJet(Jet1,0,"Randomise")).get(),0,1,0,1,0); 
	  //fShapesVar[5]=FjNSubJettiness(std::unique_ptr<AliEmcalJet>(ModifyJet(Jet1,0,"AddExtraProng_02_30")).get(),0,1,0,1,0);
	  if (fMLOn==0) fShapesVar[3]=FjNSubJettinessFastJet(ExtraProng_Jet_01_10,0,1,0,1,0);
	  else fShapesVar[3]=0.0;
	  if (fMLOn==0) fShapesVar[5]=FjNSubJettinessFastJet(ExtraProng_Jet_01_15,0,1,0,1,0);
	  else fShapesVar[5]=0.0;
	  if (fMLOn==0) fShapesVar[7]=FjNSubJettinessFastJet(ExtraProng_Jet_01_25,0,1,0,1,0);
	  else fShapesVar[7]=0.0;
	  //fShapesVar[9]=FjNSubJettinessFastJet(ExtraProng_Jet_01_30,0,1,0,1,0);
	  if (fMLOn==0) fShapesVar[9]=FjNSubJettinessFastJet(kTTrack_1_2_1,0,1,0,1,0);
	  else fShapesVar[9]=0.0;
	  if (fMLOn==0) fShapesVar[11]=FjNSubJettinessFastJet(ExtraProng_Jet_03_15,0,1,0,1,0);
	  else fShapesVar[11]=Parton_Info->GetPartonFlag7();
	  if (fMLOn==0) fShapesVar[13]=FjNSubJettinessFastJet(Randomised_Jet,0,2,0,1,0);
	  else fShapesVar[13]=Parton_Info->GetPartonEta7();
	  if (fMLOn==0) fShapesVar[15]=FjNSubJettinessFastJet(ExtraProng_Jet_01_15,0,2,0,1,0);
	  else fShapesVar[15]=Parton_Info->GetPartonPhi7();
	  if (fMLOn==0) fShapesVar[17]=FjNSubJettinessFastJet(ExtraProng_Jet_01_25,0,2,0,1,0);
	  else fShapesVar[17]=0.0;
	  if (fMLOn==0) fShapesVar[19]=FjNSubJettinessFastJet(ExtraProng_Jet_03_10,0,2,0,1,0);
	  else fShapesVar[19]=0.0;
	  if (fMLOn==0) fShapesVar[21]=FjNSubJettinessFastJet(ExtraProng_Jet_03_20,0,2,0,1,0);
	  else fShapesVar[21]=0.0;
	  if (fFullTree){
	    if (fMLOn==0) fShapesVar[23]=-2;
	    else fShapesVar[23]=0.0;
	    if (fMLOn==0) fShapesVar[25]=-2;
	    else fShapesVar[25]=0.0;
	    if (fMLOn==0) fShapesVar[27]=-2;
	    else fShapesVar[27]=0.0;
	  }
	  fTreeResponseMatrixAxis->Fill();
	  if (fMLOn==0) fTreeTracks->Fill();
	  
	  fhSubJetCounter->Fill(Reclusterer1->GetNumberOfJets());
	  for (Int_t i= 0; i<Reclusterer1->GetNumberOfJets(); i++){
	    fhEventCounter->Fill(6); //Number of overall subjets in all events 
	    fhSubJetPt->Fill(Reclusterer1->GetJet(i)->Pt());
	    fhSubJetMass->Fill(Reclusterer1->GetJet(i)->M());
	    fhNumberOfSubJetTracks->Fill(Reclusterer1->GetJet(i)->GetNumberOfTracks());                                                                            
	    fhSubJetRadius->Fill(TMath::Sqrt((Reclusterer1->GetJet(i)->Area()/TMath::Pi()))); //Radius of SubJets per event   
	  }
	  delete Randomised_Jet.second;
	  delete ExtraProng_Jet_02_15.second;
	  delete ExtraProng_Jet_02_30.second;
	  delete ExtraProng_Jet_02_45.second;
	  delete ExtraProng_Jet_02_60.second;
	  delete ExtraProng_Jet_01_30.second;
	  delete ExtraProng_Jet_01_45.second;
	  delete ExtraProng_Jet_03_30.second;
	  delete ExtraProng_Jet_03_45.second;
	}    
      }
      fhJetCounter->Fill(JetCounter); //Number of Jets in Each Event
    }
    //else {fhEventCounter->Fill(2);} //Events with no jet container
  }
  return kTRUE;
}
//________________________________________________________________________
Double_t AliAnalysisTaskSubJetFraction::RelativePhiEventPlane(Double_t EventPlane, Double_t Phi){

  if(Phi < -1*TMath::Pi()) Phi += (2*TMath::Pi());
  else if (Phi > TMath::Pi()) Phi -= (2*TMath::Pi());
  Double_t DeltaPhi=Phi-EventPlane;
  if(DeltaPhi < -1*TMath::Pi()) DeltaPhi += (2*TMath::Pi());
  else if (DeltaPhi > TMath::Pi()) DeltaPhi -= (2*TMath::Pi());
  return DeltaPhi;
}
//________________________________________________________________________
Double_t AliAnalysisTaskSubJetFraction::RelativePhi(Double_t Phi1, Double_t Phi2){

  if(Phi1 < -1*TMath::Pi()) Phi1 += (2*TMath::Pi()); // Turns the range of 0to2Pi into -PitoPi ???????????                                                             
  else if (Phi1 > TMath::Pi()) Phi1 -= (2*TMath::Pi());
  if(Phi2 < -1*TMath::Pi()) Phi2 += (2*TMath::Pi());
  else if (Phi2 > TMath::Pi()) Phi2 -= (2*TMath::Pi());
  Double_t DeltaPhi=Phi2-Phi1;
  if(DeltaPhi < -1*TMath::Pi()) DeltaPhi += (2*TMath::Pi());
  else if (DeltaPhi > TMath::Pi()) DeltaPhi -= (2*TMath::Pi());
  return DeltaPhi;
}


//--------------------------------------------------------------------------
Int_t AliAnalysisTaskSubJetFraction::SelectTriggerHadron(Float_t PtMin, Float_t PtMax){

  AliTrackContainer *PartCont = NULL;
  AliParticleContainer *PartContMC = NULL;


  if (fJetShapeSub==kConstSub){
    if (fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) PartContMC = GetParticleContainer(1);
    else PartCont = GetTrackContainer(1);
  }
  else{
    if (fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) PartContMC = GetParticleContainer(0);
    else PartCont = GetTrackContainer(0);
  }
  
  TClonesArray *TracksArray = NULL;
  TClonesArray *TracksArrayMC = NULL;
  if (fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) TracksArrayMC = PartContMC->GetArray();
  else TracksArray = PartCont->GetArray();
 
  if (fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly){
    if(!PartContMC || !TracksArrayMC) return -99999;
  }
  else {
    if(!PartCont || !TracksArray) return -99999;
  }
    
  AliAODTrack *Track = 0x0;
  Int_t Trigger_Index[100];
  for (Int_t i=0; i<100; i++) Trigger_Index[i] = 0;
  Int_t Trigger_Counter = 0;
  Int_t NTracks=0;
  if (fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly) NTracks = TracksArrayMC->GetEntriesFast();
  else NTracks = TracksArray->GetEntriesFast();
  for(Int_t i=0; i < NTracks; i++){
    if (fJetShapeType == AliAnalysisTaskSubJetFraction::kGenOnTheFly){
      if((Track = static_cast<AliAODTrack*>(PartContMC->GetAcceptParticle(i)))){
	if (!Track) continue;
	if(TMath::Abs(Track->Eta())>0.9) continue;
	if (Track->Pt()<0.15) continue;
	if ((Track->Pt() >= PtMin) && (Track->Pt()< PtMax)) {
	  Trigger_Index[Trigger_Counter] = i;
	  Trigger_Counter++;
	}
      }
    }
    else{ 
      if((Track = static_cast<AliAODTrack*>(PartCont->GetAcceptTrack(i)))){
	if (!Track) continue;
	if(TMath::Abs(Track->Eta())>0.9) continue;
	if (Track->Pt()<0.15) continue;
	if ((Track->Pt() >= PtMin) && (Track->Pt()< PtMax)) {
	  Trigger_Index[Trigger_Counter] = i;
	  Trigger_Counter++;
	}
      }
    } 
  }
  if (Trigger_Counter == 0) return -99999;
  Int_t RandomNumber = 0, Index = 0 ; 
  TRandom3* Random = new TRandom3(0); 
  RandomNumber = Random->Integer(Trigger_Counter);
  Index = Trigger_Index[RandomNumber];
  return Index; 
}


//--------------------------------------------------------------------------
Double_t AliAnalysisTaskSubJetFraction::Angularity(AliEmcalJet *Jet, Int_t JetContNb){
  
  AliJetContainer *JetCont = GetJetContainer(JetContNb);
  Double_t Angularity_Numerator=0;  //Reset these values                                                                            
  Double_t Angularity_Denominator=0;
  AliVParticle *Particle=0x0;
  Double_t DeltaPhi=0;


  for (Int_t i=0; i< Jet->GetNumberOfTracks(); i++){  //loops through all tracks (particles in the jet                                                                              
    Particle = static_cast<AliVParticle*>(Jet->TrackAt(i, JetCont->GetParticleContainer()->GetArray()));                                                                 
    if(!Particle) continue;
    DeltaPhi=RelativePhi(Jet->Phi(),Particle->Phi());
    Angularity_Numerator=Angularity_Numerator+(Particle->Pt()*TMath::Sqrt(((Particle->Eta()-Jet->Eta())*(Particle->Eta()-Jet->Eta()))+(DeltaPhi*DeltaPhi)));
    Angularity_Denominator= Angularity_Denominator+Particle->Pt();
  }
  if(Angularity_Denominator!=0) return Angularity_Numerator/Angularity_Denominator;
  else return -1;

}



//--------------------------------------------------------------------------                                                                                                         
Double_t AliAnalysisTaskSubJetFraction::PTD(AliEmcalJet *Jet, Int_t JetContNb){

  AliJetContainer *JetCont = GetJetContainer(JetContNb);
  Double_t PTD_Numerator=0;  //Reset these values                                                                                                                            
  Double_t PTD_Denominator=0;
  AliVParticle *Particle=0x0;
  Double_t DeltaPhi=0;
  for (Int_t i=0; i< Jet->GetNumberOfTracks(); i++){  //loops through all tracks (particles in the jet                                                                            
    Particle = static_cast<AliVParticle*>(Jet->TrackAt(i, JetCont->GetParticleContainer()->GetArray()));
    if(!Particle) continue;
    PTD_Numerator=PTD_Numerator+(Particle->Pt()*Particle->Pt());
    PTD_Denominator=PTD_Denominator+Particle->Pt();
  }
  if(PTD_Denominator!=0) return TMath::Sqrt(PTD_Numerator)/PTD_Denominator;
  else return -1;

}


//----------------------------------------------------------------------
///////////////returns jet finder object containg subjets
AliEmcalJetFinder *AliAnalysisTaskSubJetFraction::Recluster(AliEmcalJet *Jet, Int_t JetContNb, Double_t SubJetRadius, Double_t SubJetMinPt, Int_t Algorithm, const char* Name){

  AliJetContainer *JetCont = GetJetContainer(JetContNb);
  AliEmcalJetFinder *Reclusterer = new AliEmcalJetFinder(Name); //JetFinder Object for reclustered jets                                                                 
  Reclusterer->SetRadius(SubJetRadius); 
  Reclusterer->SetJetMinPt(SubJetMinPt);
  Reclusterer->SetJetAlgorithm(Algorithm); //0 for anti-kt     1 for kt
  Reclusterer->SetJetMaxEta(0.9);
  Reclusterer->SetRecombSheme(0);
  if(fJetShapeType != AliAnalysisTaskSubJetFraction::kGenOnTheFly){
    const AliVVertex *vert = InputEvent()->GetPrimaryVertex();
    Double_t dVtx[3]={vert->GetX(),vert->GetY(),vert->GetZ()};
    if(Reclusterer->AliEmcalJetFinder::Filter(Jet, JetCont, dVtx)){;}  //reclustering jet1 using the jetfinderobject Reclusterer
  }
  else{
    Double_t dVtx[3]={1,1,1};
    if(Reclusterer->AliEmcalJetFinder::Filter(Jet, JetCont, dVtx)){;}  //reclustering jet1 using the jetfinderobject Reclusterer
  }
  return Reclusterer;
  
}



//----------------------------------------------------------------------
Double_t AliAnalysisTaskSubJetFraction::SubJetOrdering(AliEmcalJet *Jet, AliEmcalJetFinder *Reclusterer, Int_t N, Int_t Type, Bool_t Index){
  AliEmcalJet *SubJet=NULL;
  Double_t SortingVariable;
  Int_t ArraySize =N+1;
  TArrayD JetSorter(ArraySize);
  TArrayD JetIndexSorter(ArraySize);
  for (Int_t i=0; i<ArraySize; i++){
    JetSorter[i]=0;
  }
  for (Int_t i=0; i<ArraySize; i++){
    JetIndexSorter[i]=0;
  }
  if(Reclusterer->GetNumberOfJets()<N) return -999;
  for (Int_t i=0; i<Reclusterer->GetNumberOfJets(); i++){
    SubJet=Reclusterer->GetJet(i);
    if (Type==0) SortingVariable=SubJet->Pt();
    else if (Type==1) SortingVariable=SubJet->E();
    else if (Type==2) SortingVariable=SubJet->M();
    for (Int_t j=0; j<N; j++){
      if (SortingVariable>JetSorter[j]){
	for (Int_t k=N-1; k>=j; k--){
	  JetSorter[k+1]=JetSorter[k];
	  JetIndexSorter[k+1]=JetIndexSorter[k];
	}
	JetSorter[j]=SortingVariable;
	JetIndexSorter[j]=i;
	break;
      }
    }
  }
  if (!Index) return JetSorter[N-1];
  else return JetIndexSorter[N-1];
}



//returns -1 if the Nth hardest jet is requested where N>number of available jets
//type:  0=Pt  1=E  2=M 
//Index TRUE=returns index   FALSE=returns value of quantatiy in question







//----------------------------------------------------------------------                                                                                                        
Double_t AliAnalysisTaskSubJetFraction::SubJetFraction(AliEmcalJet *Jet, AliEmcalJetFinder *Reclusterer, Int_t N, Int_t Type, Bool_t Add, Bool_t Loss){
  AliEmcalJet *SubJet=NULL;
  Double_t Observable=0;
  Double_t Fraction_Numerator=0;
  Bool_t Error=kFALSE;
  if (!Jet->GetNumberOfTracks())    return -2;
  if (Add){
    for (Int_t i=1; i<=N; i++){
      Observable=SubJetOrdering(Jet,Reclusterer,i,Type,kFALSE);
      if(Observable==-999){
	Error = kTRUE;
	return -2;
      }
      Fraction_Numerator=Fraction_Numerator+Observable;
	}
  }
  else {
    Fraction_Numerator=SubJetOrdering(Jet,Reclusterer,N,Type,kFALSE);
    if (Fraction_Numerator==-999) return -2;
  }
  if (Type==0){
    if(Loss) return (Jet->Pt()-Fraction_Numerator)/Jet->Pt();
    else return Fraction_Numerator/Jet->Pt();
  }
  else if (Type==1){
    if(Loss) return (Jet->E()-Fraction_Numerator)/Jet->E();
    else return Fraction_Numerator/Jet->E();
  }
  else { //change to else if if you want to add more later 
    if(Loss) return (Jet->M()-Fraction_Numerator)/Jet->M();
    else return Fraction_Numerator/Jet->M();
  }
}
//N number of hardest subjets involved
//Type 0=Pt 1=E 2=M
// Add TRUE: Add 1 to Nth hardest subjets together/total jet   False: Nth hardest Jet/toal jet
//Loss TRUE: Jet-Subjet(s)/Jet FALSE: Subjet(s)/Jet 


//----------------------------------------------------------------------
Double_t AliAnalysisTaskSubJetFraction::NSubJettiness(AliEmcalJet *Jet, Int_t JetContNb, Double_t JetRadius,  AliEmcalJetFinder *Reclusterer, Int_t N, Int_t A, Int_t B){
  AliJetContainer *JetCont = GetJetContainer(JetContNb);
  AliEmcalJet *SubJet=NULL;
  Double_t DeltaR1=0;
  Double_t DeltaR2=0;
  AliVParticle *JetParticle=0x0;
  Double_t SubJetiness_Numerator = 0;
  Double_t SubJetiness_Denominator = 0;
  Double_t Index=-2;
  Bool_t Error=kFALSE;
  //  JetRadius=TMath::Sqrt((Jet->Area()/TMath::Pi())); //comment out later
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
  if (SubJetiness_Denominator!=0 && !Error){
    // if (SubJetiness_Numerator/SubJetiness_Denominator>1.0) cout << JetRadius<<endl;
    return SubJetiness_Numerator/SubJetiness_Denominator;                                                                                  }
  else return -2;
}


//______________________________________________________________________________________
Double_t AliAnalysisTaskSubJetFraction::FjNSubJettiness(AliEmcalJet *Jet, Int_t JetContNb,Int_t N, Int_t Algorithm, Double_t Beta, Int_t Option, Double_t Beta_SD, Double_t ZCut){

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


  //fSoftDropOn==1    All the below are done on the jet after it has undergone softdrop with CA. Option 3 and 4 are done whether this is turned on or not.

  
  //Option==0 returns Nsubjettiness Value
  //Option==1 && N==2 returns opening angle between two subjet axes(Delta R?)
  //Option==2 && N==2 returns Delta R
  //Option==3 returns first splitting distance for soft dropped jet (CA)
  //Option==4 returns Symmetry measure (Zg) for soft dropped jet  (CA)
  //Option==5 returns Pt of Subjet1
  //Option==6 returns Pt of Subjet2
  //Options==7 trutns deltaR of subjets...Is this different to before??
  //Option==8 Subjet1 Leading track Pt
  //Option==9 Subjet1 Leading track Pt

  
  Algorithm=fReclusteringAlgorithm;
  if (Jet->GetNumberOfTracks()>=N){
    if((fJetShapeSub==kDerivSub) && (JetContNb==0) && (N==1) && (Algorithm==0) && (Beta==1.0) && (Option==0)){
      if (fDerivSubtrOrder == kFirstOrder) return Jet->GetShapeProperties()->GetFirstOrderSubtracted1subjettiness_kt();
      else return Jet->GetShapeProperties()->GetSecondOrderSubtracted1subjettiness_kt();
    }
    else if((fJetShapeSub==kDerivSub) && (JetContNb==0) && (N==2) && (Algorithm==0) && (Beta==1.0) && (Option==0)){
      if (fDerivSubtrOrder == kFirstOrder) return Jet->GetShapeProperties()->GetFirstOrderSubtracted2subjettiness_kt();
      else return Jet->GetShapeProperties()->GetSecondOrderSubtracted2subjettiness_kt();
    }
    else if((fJetShapeSub==kDerivSub) && (JetContNb==0) && (N==3) && (Algorithm==0) && (Beta==1.0) && (Option==0)){
      if (fDerivSubtrOrder == kFirstOrder) return Jet->GetShapeProperties()->GetFirstOrderSubtracted3subjettiness_kt();
      else return Jet->GetShapeProperties()->GetSecondOrderSubtracted3subjettiness_kt();
    }
    else if((fJetShapeSub==kDerivSub) && (JetContNb==0) && (N==2) && (Algorithm==0) && (Beta==1.0) && (Option==1)){
      if (fDerivSubtrOrder == kFirstOrder) return Jet->GetShapeProperties()->GetFirstOrderSubtractedOpeningAngle_kt();
      else return Jet->GetShapeProperties()->GetSecondOrderSubtractedOpeningAngle_kt();
    }
    else if((fJetShapeSub==kDerivSub) && (JetContNb==0) && (N==1) && (Algorithm==1) && (Beta==1.0) && (Option==0)){
      if (fDerivSubtrOrder == kFirstOrder) return Jet->GetShapeProperties()->GetFirstOrderSubtracted1subjettiness_ca();
      else return Jet->GetShapeProperties()->GetSecondOrderSubtracted1subjettiness_ca();
    }
    else if((fJetShapeSub==kDerivSub) && (JetContNb==0) && (N==2) && (Algorithm==1) && (Beta==1.0) && (Option==0)){
      if (fDerivSubtrOrder == kFirstOrder) return Jet->GetShapeProperties()->GetFirstOrderSubtracted2subjettiness_ca();
      else return Jet->GetShapeProperties()->GetSecondOrderSubtracted2subjettiness_ca();
    }
    else if((fJetShapeSub==kDerivSub) && (JetContNb==0) && (N==2) && (Algorithm==1) && (Beta==1.0) && (Option==1)){
      if (fDerivSubtrOrder == kFirstOrder) return Jet->GetShapeProperties()->GetFirstOrderSubtractedOpeningAngle_ca();
      else return Jet->GetShapeProperties()->GetSecondOrderSubtractedOpeningAngle_ca();
    }   
    else if((fJetShapeSub==kDerivSub) && (JetContNb==0) && (N==1) && (Algorithm==2) && (Beta==1.0) && (Option==0)){
      if (fDerivSubtrOrder == kFirstOrder) return Jet->GetShapeProperties()->GetFirstOrderSubtracted1subjettiness_akt02();
      else return Jet->GetShapeProperties()->GetSecondOrderSubtracted1subjettiness_akt02();
    }
    else if((fJetShapeSub==kDerivSub) && (JetContNb==0) && (N==2) && (Algorithm==2) && (Beta==1.0) && (Option==0)){
      if (fDerivSubtrOrder == kFirstOrder) return Jet->GetShapeProperties()->GetFirstOrderSubtracted2subjettiness_akt02();
      else return Jet->GetShapeProperties()->GetSecondOrderSubtracted2subjettiness_akt02();
    }
    else if((fJetShapeSub==kDerivSub) && (JetContNb==0) && (N==2) && (Algorithm==2) && (Beta==1.0) && (Option==1)){
      if (fDerivSubtrOrder == kFirstOrder) return Jet->GetShapeProperties()->GetFirstOrderSubtractedOpeningAngle_akt02();
      else return Jet->GetShapeProperties()->GetSecondOrderSubtractedOpeningAngle_akt02();
    }
    else if((fJetShapeSub==kDerivSub) && (JetContNb==0) && (N==1) && (Algorithm==6) && (Beta==1.0) && (Option==0)){
      if (fDerivSubtrOrder == kFirstOrder) return Jet->GetShapeProperties()->GetFirstOrderSubtracted1subjettiness_onepassca();
      else return Jet->GetShapeProperties()->GetSecondOrderSubtracted1subjettiness_onepassca();
    }
    else if((fJetShapeSub==kDerivSub) && (JetContNb==0) && (N==2) && (Algorithm==6) && (Beta==1.0) && (Option==0)){
      if (fDerivSubtrOrder == kFirstOrder) return Jet->GetShapeProperties()->GetFirstOrderSubtracted2subjettiness_onepassca();
      else return Jet->GetShapeProperties()->GetSecondOrderSubtracted2subjettiness_onepassca();
    }
    else if((fJetShapeSub==kDerivSub) && (JetContNb==0) && (N==2) && (Algorithm==6) && (Beta==1.0) && (Option==1)){
      if (fDerivSubtrOrder == kFirstOrder) return Jet->GetShapeProperties()->GetFirstOrderSubtractedOpeningAngle_onepassca();
      else return Jet->GetShapeProperties()->GetSecondOrderSubtractedOpeningAngle_onepassca();
    } 
    else{
      //Algorithm=fReclusteringAlgorithm;
      AliJetContainer *JetCont = GetJetContainer(JetContNb);
      AliEmcalJetFinder JetFinder("Nsubjettiness");
      JetFinder.SetJetMaxEta(0.9-fJetRadius);
      JetFinder.SetRadius(fJetRadius); 
      JetFinder.SetJetAlgorithm(0); //0 for anti-kt     1 for kt  //this is for the JET!!!!!!!!!! Not the SubJets
      JetFinder.SetRecombSheme(0);
      JetFinder.SetJetMinPt(Jet->Pt());
      if(fJetShapeType != AliAnalysisTaskSubJetFraction::kGenOnTheFly){
	const AliVVertex *vert = InputEvent()->GetPrimaryVertex();
	Double_t dVtx[3]={vert->GetX(),vert->GetY(),vert->GetZ()};
	return JetFinder.Nsubjettiness(Jet,JetCont,dVtx,N,Algorithm,fSubJetRadius,Beta,Option,0,Beta_SD,ZCut,fSoftDropOn);
      }
      else{
	Double_t dVtx[3]={1,1,1};
	if (!fNsubMeasure){
	  return JetFinder.Nsubjettiness(Jet,JetCont,dVtx,N,Algorithm,fSubJetRadius,Beta,Option,0,Beta_SD,ZCut,fSoftDropOn);
	}
	else{
	  if (Option==3) return JetFinder.Nsubjettiness(Jet,JetCont,dVtx,N,Algorithm,fSubJetRadius,Beta,Option,0,Beta_SD,ZCut,fSoftDropOn);
	  else return JetFinder.Nsubjettiness(Jet,JetCont,dVtx,N,Algorithm,fSubJetRadius,Beta,Option,1);
	}
      }
    }
  }
  else return -2;
}


//________________________________________________________________________
Bool_t AliAnalysisTaskSubJetFraction::RetrieveEventObjects() {
  //
  // retrieve event objects
  //
  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  return kTRUE;
}


//_______________________________________________________________________
void AliAnalysisTaskSubJetFraction::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
  if (fJetShapeType==AliAnalysisTaskSubJetFraction::kData){

    /*
    for (int i=1; i<=(XBinsJetPtSize)-1; i++){ //Rescales the fhJetPt Histograms according to the with of the variable bins and number of events
      fhJetPt->SetBinContent(i,((1.0/(XBinsJetPt[i]-XBinsJetPt[i-1]))*((fhJetPt->GetBinContent(i))*(1.0/(fhEventCounter->GetBinContent(1))))));      
      fhSubJetPt->SetBinContent(i,((1.0/(XBinsJetPt[i]-XBinsJetPt[i-1]))*((fhSubJetPt->GetBinContent(i))*(1.0/(fhEventCounter->GetBinContent(5))))));
      
      //fhJetPt->SetBinContent(i,((1.0/(fhPt->GetBinWidth(i)))*((fhJetPt->GetBinContent(i))*(fhEventCounter->GetBinContent(1)))));   
      // fhJetPt->SetBinContent(i,((1.0/((fhPt->GetLowEdge(i+1))-(fhJetPt->GetLowEdge(i))))*((fhPt->GetBinContent(i))*(fhEventCounter->GetBinContent(1)))));
    }
    for (int i=1; i<=(XBinsJetMassSize)-1; i++){ //Rescales the fhPt Histograms according to the with of the variable bins and number of events                     
      fhJetMass->SetBinContent(i,((1.0/(XBinsJetMass[i]-XBinsJetMass[i-1]))*((fhJetMass->GetBinContent(i))*(1.0/(fhEventCounter->GetBinContent(1))))));                                                                            
    }
    
    fhJetPhi->Scale(Phi_Bins/((Phi_Up-Phi_Low)*((fhEventCounter->GetBinContent(1)))));
    fhJetEta->Scale(Eta_Bins/((Eta_Up-Eta_Low)*((fhEventCounter->GetBinContent(1)))));
    fhJetRadius->Scale(100/(fhEventCounter->GetBinContent(4)));  //should this and JetAngularity be divided by Bin 1 or 4???? 
    fhNumberOfJetTracks->Scale(1.0/(fhEventCounter->GetBinContent(4)));
    fhJetCounter->Scale(1.0/(fhEventCounter->GetBinContent(1)));  //is the first bin the correct one to look at?
   
 fhSubJetCounter->Scale(1.0/(fhEventCounter->GetBinContent(5)));
    */
    /*  
  fhJetPt->Scale(1.0/(fhEventCounter->GetBinContent(1)));                                      
  fhSubJetPt->Scale(1.0/(fhEventCounter->GetBinContent(5)));                                
  fhJetMass->Scale(1.0/(fhEventCounter->GetBinContent(1)));                                                                                                                     
  fhJetPhi->Scale(1.0/(fhEventCounter->GetBinContent(1)));                                                                                               
  fhJetEta->Scale(1.0/(fhEventCounter->GetBinContent(1)));                                                                                               
  fhJetRadius->Scale(1.0/(fhEventCounter->GetBinContent(4)));  //should this and JetAngularity be divided by Bin 1 or 4????                                                                                                                                                                                                                                           
  fhNumberOfJetTracks->Scale(1.0/(fhEventCounter->GetBinContent(4)));                                                                                                                                                                                                                           
  fhJetCounter->Scale(1.0/(fhEventCounter->GetBinContent(1)));  //is the first bin the correct one to look at?                                                                     
  fhSubJetCounter->Scale(1.0/(fhEventCounter->GetBinContent(5)));            
    */
  }
  
}

///////////Function to modify jets based on Toy Models/////////
std::pair<fastjet::PseudoJet,fastjet::ClusterSequence *> AliAnalysisTaskSubJetFraction::ModifyJet(AliEmcalJet* Jet, Int_t JetContNb, TString Modification){
  std::pair<fastjet::PseudoJet,fastjet::ClusterSequence *> Jet_ClusterSequence;
  std::vector<fastjet::PseudoJet> fInputVectors;
  AliJetContainer *fJetCont = GetJetContainer(JetContNb);
  AliParticleContainer *fTrackCont = fJetCont->GetParticleContainer();
  std::vector<fastjet::PseudoJet> Modified_Jet;
  fastjet::ClusterSequence *fClustSeqSA;
  Int_t NJetTracksTest = Jet->GetNumberOfTracks();
  if (fTrackCont) for (Int_t i=0; i<Jet->GetNumberOfTracks(); i++) {
      AliVParticle *fTrk = Jet->TrackAt(i, fTrackCont->GetArray());
      if (!fTrk) continue;    
      fastjet::PseudoJet PseudoTracks(fTrk->Px(), fTrk->Py(), fTrk->Pz(),fTrk->E());
      PseudoTracks.set_user_index(Jet->TrackAt(i)+100); //100 is very arbitary....why is it here anyway?
      fInputVectors.push_back(PseudoTracks);
    }
  if (Modification=="Randomise") fInputVectors=RandomiseTracks(Jet,fInputVectors);
  if (Modification=="AddExtraProng_02_15") fInputVectors=AddExtraProng(fInputVectors,0.2,0.15);
  if (Modification=="AddExtraProng_02_30") fInputVectors=AddExtraProng(fInputVectors,0.2,0.30);
  if (Modification=="AddExtraProng_02_45") fInputVectors=AddExtraProng(fInputVectors,0.2,0.45);
  if (Modification=="AddExtraProng_02_60") fInputVectors=AddExtraProng(fInputVectors,0.2,0.60);
  if (Modification=="AddExtraProng_01_10") fInputVectors=AddExtraProng(fInputVectors,0.1,0.10);
  if (Modification=="AddExtraProng_01_15") fInputVectors=AddExtraProng(fInputVectors,0.1,0.15);
  if (Modification=="AddExtraProng_01_20") fInputVectors=AddExtraProng(fInputVectors,0.1,0.20);
  if (Modification=="AddExtraProng_01_25") fInputVectors=AddExtraProng(fInputVectors,0.1,0.25); 
  if (Modification=="AddExtraProng_01_30") fInputVectors=AddExtraProng(fInputVectors,0.1,0.30);
  if (Modification=="AddExtraProng_01_45") fInputVectors=AddExtraProng(fInputVectors,0.1,0.45);
  if (Modification=="AddExtraProng_03_10") fInputVectors=AddExtraProng(fInputVectors,0.3,0.10);
  if (Modification=="AddExtraProng_03_15") fInputVectors=AddExtraProng(fInputVectors,0.3,0.15);
  if (Modification=="AddExtraProng_03_20") fInputVectors=AddExtraProng(fInputVectors,0.3,0.20);
  if (Modification=="AddExtraProng_03_30") fInputVectors=AddExtraProng(fInputVectors,0.3,0.30);
  if (Modification=="AddExtraProng_03_45") fInputVectors=AddExtraProng(fInputVectors,0.3,0.45);
  if (Modification=="AddkTTrack_1_2_1") fInputVectors=AddkTTracks(Jet,fInputVectors,1.0,2.0,1);
  fJetRadius=fJetRadius*2.0;
  try {
    fastjet::JetDefinition fJetDef(fastjet::antikt_algorithm, fJetRadius, static_cast<fastjet::RecombinationScheme>(0), fastjet::BestFJ30 );         
    fClustSeqSA=new fastjet::ClusterSequence(fInputVectors, fJetDef);
    Modified_Jet=fClustSeqSA->inclusive_jets(0);
    Jet_ClusterSequence.second=fClustSeqSA; 
  } catch (fastjet::Error) {
    AliError(" [w] FJ Exception caught.");
    //return -1;
  }
  fJetRadius=fJetRadius/2.0;
  Jet_ClusterSequence.first=Modified_Jet[0];
  return Jet_ClusterSequence;
}


std::vector<fastjet::PseudoJet> AliAnalysisTaskSubJetFraction::RandomiseTracks(AliEmcalJet *Jet,std::vector<fastjet::PseudoJet> fInputVectors){
  TRandom3 fRandom1;
  fRandom1.SetSeed(0);
  Double_t Random_Phi=0.0;
  // Double_t Random_Phi_Change=0.0;
  Double_t Random_Eta=0.0;
  // Double_t Random_Eta_Change=0.0;
  Double_t Track_Pt=0.0;
  std::vector<fastjet::PseudoJet> Random_Track_Vector;
  Random_Track_Vector.clear();
  for (Int_t i=0; i< fInputVectors.size(); i++){
    do{
    Random_Phi=fRandom1.Uniform(Jet->Phi()-fJetRadius,Jet->Phi()+fJetRadius);
    Random_Eta=fRandom1.Uniform(Jet->Eta()-fJetRadius,Jet->Eta()+fJetRadius);
    }while(TMath::Sqrt(TMath::Power(RelativePhi(Jet->Phi(),Random_Phi),2)+TMath::Power(Jet->Eta()-Random_Eta,2))>fJetRadius);
    /*    Random_Phi_Change=fRandom1.Uniform(-1*fJetRadius,fJetRadius);
    Random_Phi=Jet->Phi()+Random_Phi_Change;
    Random_Eta_Change=fRandom1.Uniform(-1*(TMath::Sqrt((fJetRadius*fJetRadius)-(Random_Phi_Change*Random_Phi_Change))),TMath::Sqrt((fJetRadius*fJetRadius)-(Random_Phi_Change*Random_Phi_Change)));
    Random_Eta=Jet->Eta()+Random_Eta_Change;
    */
    if (fRandomisationEqualPt) Track_Pt=Jet->Pt()/fInputVectors.size();
    else Track_Pt=fInputVectors[i].perp();
    fastjet::PseudoJet Random_Track(Track_Pt*TMath::Cos(Random_Phi),Track_Pt*TMath::Sin(Random_Phi),Track_Pt*TMath::SinH(Random_Eta),TMath::Sqrt((Track_Pt*Track_Pt)+(Track_Pt*TMath::SinH(Random_Eta)*Track_Pt*TMath::SinH(Random_Eta))));
    Random_Track.set_user_index(i);
    Random_Track_Vector.push_back(Random_Track);
  }
  fInputVectors.clear();
  return Random_Track_Vector;
}

std::vector<fastjet::PseudoJet> AliAnalysisTaskSubJetFraction::AddExtraProng(std::vector<fastjet::PseudoJet> fInputVectors, Double_t Distance, Double_t PtFrac){
  TRandom3 fRandom1;
  fRandom1.SetSeed(0);
  std::vector<fastjet::PseudoJet> Jet;
  Jet.clear();
  try {
    fastjet::JetDefinition fJetDef(fastjet::antikt_algorithm, fJetRadius*2, static_cast<fastjet::RecombinationScheme>(0), fastjet::BestFJ30 ); 
    fastjet::ClusterSequence fClustSeqSA_Prong(fInputVectors, fJetDef);
    Jet= fClustSeqSA_Prong.inclusive_jets(0);
  } catch (fastjet::Error) {
    AliError(" [w] FJ Exception caught.");
  }
  Double_t Extra_Track_Phi_Change=fRandom1.Uniform(-1*Distance,Distance);
  Double_t Extra_Track_Phi=Jet[0].phi()+Extra_Track_Phi_Change;
  Double_t Extra_Track_Eta_Change=fRandom1.Uniform(-1*(TMath::Sqrt((Distance*Distance)-(Extra_Track_Phi_Change*Extra_Track_Phi_Change))),TMath::Sqrt((Distance*Distance)-(Extra_Track_Phi_Change*Extra_Track_Phi_Change)));
  Double_t Extra_Track_Eta=Jet[0].pseudorapidity()+Extra_Track_Eta_Change;
  Double_t Extra_Track_Pt=Jet[0].perp()*PtFrac;
  fastjet::PseudoJet ExtraTrack(Extra_Track_Pt*TMath::Cos(Extra_Track_Phi),Extra_Track_Pt*TMath::Sin(Extra_Track_Phi),Extra_Track_Pt*TMath::SinH(Extra_Track_Eta),TMath::Sqrt((Extra_Track_Pt*Extra_Track_Pt)+(Extra_Track_Pt*TMath::SinH(Extra_Track_Eta)*Extra_Track_Pt*TMath::SinH(Extra_Track_Eta))));
  ExtraTrack.set_user_index(150);
  fInputVectors.push_back(ExtraTrack); 
  return fInputVectors;
}

std::vector<fastjet::PseudoJet> AliAnalysisTaskSubJetFraction::AddkTTracks(AliEmcalJet *Jet,std::vector<fastjet::PseudoJet> fInputVectors, Double_t QHat,Double_t XLength, Int_t NAdditionalTracks)
{
//here add N tracks with random phi and eta and theta according to bdmps distrib.
  fastjet::PseudoJet Lab_Jet;
  Lab_Jet.reset(Jet->Px(),Jet->Py(),Jet->Pz(),Jet->E());
  // Double_t Omega_C=0.5*QHat*XLength*XLength/0.2;
  // Double_t Theta_C=TMath::Sqrt(12*0.2/(QHat*TMath::Power(XLength,3)));
  Double_t xQs=TMath::Sqrt(QHat*XLength); //QHat=1,XLength=2				
    //cout<<"medium parameters "<<Omega_C<<" "<<Theta_C<<" "<<xQs<<endl;

   for(Int_t i=0;i<NAdditionalTracks;i++){

     Double_t kT_Scale,Limit_Min,Limit_Max;
    
    //generation of kT according to 1/kT^4, with minimum QS^2=2 GeV and maximum ~sqrt(ptjet*T)	
     TF1 *kT_Func= new TF1("kT_Func","1/(x*x*x*x)",xQs*xQs,10000);
     kT_Scale=kT_Func->GetRandom();
     //generation within the jet cone
    
     //generation of w according to 1/w, with minimum wc
     //omega needs to be larger than kT so to have well defined angles
     Limit_Min=kT_Scale;
     Limit_Max=kT_Scale/TMath::Sin(0.01);
     TF1 *Omega_Func= new TF1("Omega_Func","1/x",Limit_Min,Limit_Max);
     Double_t Omega=Omega_Func->GetRandom();
     
     Double_t Theta=TMath::ASin(kT_Scale/Omega);
     //cout<<"angle_omega_kt"<<Theta<<" "<<Omega<<" "<<kT_Scale<<endl;
     if(Theta>fJetRadius) continue;

     fastjet::PseudoJet ExtraTrack(kT_Scale/TMath::Sqrt(2),kT_Scale/TMath::Sqrt(2),Omega*TMath::Cos(Theta),Omega);
     //boost the particle in the rest frame of the jet to the lab frame
     fastjet::PseudoJet ExtraTrack_LabFrame=ExtraTrack.boost(Lab_Jet);
     ExtraTrack_LabFrame.set_user_index(i+Jet->GetNumberOfTracks()+100);											 
     fInputVectors.push_back(ExtraTrack_LabFrame);
     delete kT_Func;
     delete Omega_Func;
   }
     return fInputVectors;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________________
Double_t AliAnalysisTaskSubJetFraction::FjNSubJettinessFastJet(std::pair<fastjet::PseudoJet,fastjet::ClusterSequence *> Jet_ClusterSequence, Int_t JetContNb,Int_t N, Int_t Algorithm, Double_t Beta, Int_t Option, Double_t Beta_SD, Double_t ZCut){

  //Only Works for kGenOnTheFly
  
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


  //fSoftDropOn==1    All the below are done on the jet after it has undergone softdrop with CA. Option 3 and 4 are done whether this is turned on or not.

  
  //Option==0 returns Nsubjettiness Value
  //Option==1 && N==2 returns opening angle between two subjet axes(Delta R?)
  //Option==2 && N==2 returns Delta R
  //Option==3 returns first splitting distance for soft dropped jet (CA)
  //Option==4 returns Symmetry measure (Zg) for soft dropped jet  (CA)
  //Option==5 returns Pt of Subjet1
  //Option==6 returns Pt of Subjet2
  //Options==7 trutns deltaR of subjets...Is this different to before??
  //Option==8 Subjet1 Leading track Pt
  //Option==9 Subjet1 Leading track Pt

  fastjet::PseudoJet Jet=Jet_ClusterSequence.first;
  AliFJWrapper FJWrapper("FJWrapper", "FJWrapper");
  FJWrapper.SetR(fJetRadius);
  FJWrapper.SetMinJetPt(Jet.perp());
  FJWrapper.SetAlgorithm(fastjet::antikt_algorithm);  //this is for the jet clustering not the subjet reclustering. 
  FJWrapper.SetRecombScheme(static_cast<fastjet::RecombinationScheme>(0));
  if (Jet.constituents().size()>=N){
    Algorithm=fReclusteringAlgorithm;
    Double_t dVtx[3]={1,1,1};
    return FJWrapper.NSubjettinessDerivativeSub(N,Algorithm,fSubJetRadius,Beta,fJetRadius,Jet,Option,0,Beta_SD,ZCut,fSoftDropOn); //change this
  }
  else return -2;
}
