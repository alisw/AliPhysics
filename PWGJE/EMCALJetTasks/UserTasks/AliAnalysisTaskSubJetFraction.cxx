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

//Globals


Double_t XBinsPt[66]={0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.9,0.95,1.00,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00,2.20,2.40,2.60,2.80,3.00,3.20,3.40,3.60,3.80,4.00,4.50,5.00,5.50,6.00,6.50,7.00,8.00,9.00,10.00,11.00,12.00,13.00,14.00,15.00,16.00,18.00,20.00,22.00,24.00,26.00,28.00,30.00,32.00,34.00,36.00,40.00,45.00,50.00}; //Lower edges of new bins for rebbing purposes for track pt                                                                                                           
Double_t XBinsPtSize=66;
Double_t XBinsJetPt[35]={0,0.50,1.00,2.00,3.00,4.00,5.00,6.00,7.00,8.00,9.00,10.00,12.00,14.00,16.00,18.00,20.00,25.00,30.00,35.00,40.00,45.00,50.00,60.00,70.00,80.00,90.00,100.00,120.00,140.00,160.00,180.00,200.00,250.00,300.00}; //for jet pt
Int_t XBinsJetPtSize=35;
Double_t XBinsJetMass[28]={0,0.50,1.00,2.00,3.00,4.00,5.00,6.00,7.00,8.00,9.00,10.00,12.00,14.00,16.00,18.00,20.00,25.00,30.00,35.00,40.00,45.00,50.00,60.00,70.00,80.00,90.00,100.00}; //for jet mass
Int_t XBinsJetMassSize=28;
Double_t Eta_Up=1.00;
Double_t Eta_Low=-1.00;
Int_t Eta_Bins=100;
Double_t Phi_Up=2*(TMath::Pi());
Double_t Phi_Low=(-1*(TMath::Pi()));
Int_t Phi_Bins=100;





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
  fSubJetAlgorithm(0),
  fSubJetRadius(0.1),
  fSubJetMinPt(1),
  fJetRadius(0.4),
  fRMatched(0.2),
  fSharedFractionPtMin(0.5),
  fDerivSubtrOrder(0),
  fFullTree(kFALSE),
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
  fhJetAngularity(0x0),
  fhJetAngularityJetPt(0x0),
  fhJetPTD(0x0),
  fhJetPTDJetPt(0x0),
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
  fhSubJetPtFrac(0x0),
  fhSubJetPtFrac2(0x0),
  fhSubJetPtLoss(0x0),
  fhSubJetPtLoss2(0x0),
  fhSubJetEnergyFrac(0x0),
  fhSubJetEnergyFrac2(0x0),
  fhSubJetEnergyLoss(0x0),
  fhSubJetEnergyLoss2(0x0),
  fh2PtRatio(0x0),
  fhSubJetiness1(0x0),
  fhSubJetiness1JetPt(0x0),
  fhSubJetiness2(0x0),
  fhSubJetiness2JetPt(0x0),
  fh2to1SubJetinessRatio(0x0),
  fh2to1SubJetinessRatioJetPt(0x0),
  fhEventCounter(0x0),
  fhEventCounter_1(0x0),
  fhEventCounter_2(0x0),
  fhJetPtJetEta(0x0),
  fhSubJetiness2Distance(0x0),
  fh2JetTracksEtaPhiPt(0x0),
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
  fTreeResponseMatrixAxis(0)

{
  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskSubJetFraction::AliAnalysisTaskSubJetFraction(const char *name) : 
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
  fSubJetAlgorithm(0),
  fSubJetRadius(0.1),
  fSubJetMinPt(1),
  fJetRadius(0.4),
  fRMatched(0.2),
  fSharedFractionPtMin(0.5),
  fDerivSubtrOrder(0),
  fFullTree(kFALSE),
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
  fhJetAngularity(0x0),
  fhJetAngularityJetPt(0x0),
  fhJetPTD(0x0),
  fhJetPTDJetPt(0x0),
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
  fhSubJetPtFrac(0x0),
  fhSubJetPtFrac2(0x0),
  fhSubJetPtLoss(0x0),
  fhSubJetPtLoss2(0x0),
  fhSubJetEnergyFrac(0x0),
  fhSubJetEnergyFrac2(0x0),
  fhSubJetEnergyLoss(0x0),
  fhSubJetEnergyLoss2(0x0),
  fh2PtRatio(0x0),
  fhSubJetiness1(0x0),
  fhSubJetiness1JetPt(0x0),
  fhSubJetiness2(0x0),
  fhSubJetiness2JetPt(0x0),
  fh2to1SubJetinessRatio(0x0),
  fh2to1SubJetinessRatioJetPt(0x0),
  fhEventCounter(0x0),
  fhEventCounter_1(0x0),
  fhEventCounter_2(0x0),
  fhJetPtJetEta(0x0),
  fhSubJetiness2Distance(0x0),
  fh2JetTracksEtaPhiPt(0x0),
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
  fTreeResponseMatrixAxis(0)
  
{
  // Standard constructor.
  
  SetMakeGeneralHistograms(kTRUE);

  DefineOutput(1, TTree::Class());

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

  //create a tree used for the MC data and making a 4D response matrix
  fTreeResponseMatrixAxis = new TTree("fTreeJetShape", "fTreeJetShape");
  if (fFullTree){
    const Int_t nVar = 18;
    fShapesVar = new Double_t [nVar]; //shapes used for tagging   
    TString *fShapesVarNames = new TString [nVar];
  
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
    fShapesVarNames[12] = "DeltaR";
    fShapesVarNames[13] = "DeltaR_Truth";
    fShapesVarNames[14] = "Frac1";
    fShapesVarNames[15] = "Frac1_Truth";
    fShapesVarNames[16] = "Frac2";
    fShapesVarNames[17] = "Frac2_Truth";
    for(Int_t ivar=0; ivar < nVar; ivar++){
      cout<<"looping over variables"<<endl;
      fTreeResponseMatrixAxis->Branch(fShapesVarNames[ivar].Data(), &fShapesVar[ivar], Form("%s/D", fShapesVarNames[ivar].Data()));
    }
  }
  
    if (!fFullTree){
    const Int_t nVar = 12;
    fShapesVar = new Double_t [nVar]; //shapes used for tagging   
    TString *fShapesVarNames = new TString [nVar];
  
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
    for(Int_t ivar=0; ivar < nVar; ivar++){
      cout<<"looping over variables"<<endl;
      fTreeResponseMatrixAxis->Branch(fShapesVarNames[ivar].Data(), &fShapesVar[ivar], Form("%s/D", fShapesVarNames[ivar].Data()));
    }
  }
  
  
  if (fJetShapeType==AliAnalysisTaskSubJetFraction::kData || fJetShapeType==AliAnalysisTaskSubJetFraction::kSim || fJetShapeType==AliAnalysisTaskSubJetFraction::kGenOnTheFly){
    
    //  fhJetPt= new TH1F("fhJetPt", "Jet Pt", (XBinsJetPtSize)-1, XBinsJetPt);
    fhJetPt= new TH1F("fhJetPt", "Jet Pt",1500,-0.5,149.5 );   
    fOutput->Add(fhJetPt);
    //fhJetPhi= new TH1F("fhJetPhi", "Jet Phi", Phi_Bins, Phi_Low, Phi_Up);
    fhJetPhi= new TH1F("fhJetPhi", "Jet Phi",360 , -1.5*(TMath::Pi()), 1.5*(TMath::Pi()));
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
  }
  if(fJetShapeType==AliAnalysisTaskSubJetFraction::kSim || fJetShapeType==AliAnalysisTaskSubJetFraction::kGenOnTheFly){
    fhJetAngularity= new TH1F("fhJetAngularity", "Jet Angularity", 400, -2,2);
    fOutput->Add(fhJetAngularity);
    fhJetAngularityJetPt= new TH2F("fhJetAngularityJetPt", "Jet Angularity vs Jet Pt", 1500, -0.5, 149.50, 400, -2,2);
    fOutput->Add(fhJetAngularityJetPt);
    fhJetPTD= new TH1F("fhJetPTD", "Jet PTD", 400, -2,2);
    fOutput->Add(fhJetPTD);
    fhJetPTDJetPt= new TH2F("fhJetPTDJetPt", "Jet PTD vs Jet Pt", 1500, -0.5, 149.50, 400, -2,2);
    fOutput->Add(fhJetPTDJetPt);
    fhSubJetPtFrac= new TH1F("fhSubJetPtFrac", "Pt Fraction of Highest Pt Subjet compared to original Jet",101, -0.05,1.05);
    fOutput->Add(fhSubJetPtFrac);
    fhSubJetPtFrac2= new TH1F("fhSubJetPtFrac2", "Pt Fraction of Two Highest Pt Subjets compared to original Jet",101, -0.05,1.05);
    fOutput->Add(fhSubJetPtFrac2);
    fhSubJetPtLoss= new TH1F("fhSubJetPtLoss", "Pt Difference of Highest Pt Subjet compared to original Jet",101, -0.05,1.05);
    fOutput->Add(fhSubJetPtLoss);
    fhSubJetPtLoss2= new TH1F("fhSubJetPtLoss2", "Pt Difference of Two Highest Pt Subjets compared to original Jet",101, -0.05,1.05);
    fOutput->Add(fhSubJetPtLoss2);
    fhSubJetEnergyFrac= new TH1F("fhSubJetEnergyFrac", "Energy Fraction of Most Energetic Subjet compared to original Jet",101, -0.05,1.05);
    fOutput->Add(fhSubJetEnergyFrac);
    fhSubJetEnergyFrac2= new TH1F("fhSubJetEnergyFrac2", "Energy Fraction of Two Most Energetic Subjets compared to original Jet",101, -0.05,1.05);
    fOutput->Add(fhSubJetEnergyFrac2);
    fhSubJetEnergyLoss= new TH1F("fhSubJetEnergyLoss", "Energy Difference of Most Energetic Subjet compared to original Jet",101, -0.05,1.05);
    fOutput->Add(fhSubJetEnergyLoss);
    fhSubJetEnergyLoss2= new TH1F("fhSubJetEnergyLoss2", "Pt Difference of Two Most Energetic Subjets compared to original Jet",101, -0.05,1.05);
    fOutput->Add(fhSubJetEnergyLoss2);
    fhSubJetiness1= new TH1F("fhSubjetiness1", "Tau 1 value",101, -0.05,1.05);
    fOutput->Add(fhSubJetiness1);
    fhSubJetiness1JetPt= new TH2F("fhSubjetiness1JetPt", "Tau 1 value vs Jet Pt",1500, -0.5, 149.5, 101, -0.05,1.05);
    fOutput->Add(fhSubJetiness1JetPt);
    fhSubJetiness2= new TH1F("fhSubjetiness2", "Tau 2 value",101, -0.05,1.05);
    fOutput->Add(fhSubJetiness2);
    fhSubJetiness2JetPt= new TH2F("fhSubjetiness2JetPt", "Tau 2 value vs Jet Pt",1500, -0.5, 149.5, 101, -0.05,1.05);
    fOutput->Add(fhSubJetiness2JetPt);
    fh2to1SubJetinessRatio= new TH1F("fh2to1SubJetinessRatio", "Ratio of #tau 1 to #tau 2",200, -0.5,1.5);
    fOutput->Add(fh2to1SubJetinessRatio);
    fh2to1SubJetinessRatioJetPt= new TH2F("fh2to1SubJetinessRatioJetPt", "Ratio of #tau 1 to #tau 2 vs Jet Pt", 1500, -0.5, 149.5, 200, -0.5,1.5);
    fOutput->Add(fh2to1SubJetinessRatioJetPt);
    fhJetPtJetEta = new TH2F("fhJetPtJetEta", "Jet Pt vs Jet Eta", 1500,-0.5,150,Eta_Bins, Eta_Low, Eta_Up);
    fOutput->Add(fhJetPtJetEta);
    fhSubJetiness2Distance = new TH2F("fhSubJetiness2Distance", "#tau 2 as a function of distance between subject axis",100,0,10,101,-0.05,1.05);
    fOutput->Add(fhSubJetiness2Distance);
    fhSubJettiness1_FJ_KT= new TH1D("fhSubJettiness1_FJ_KT","fhSubJettiness1_FJ_KT",400,-2,2);
    fOutput->Add(fhSubJettiness1_FJ_KT);
    fhSubJettiness1_FJ_MIN= new TH1D("fhSubJettiness1_FJ_MIN","fhSubJettiness1_FJ_MIN",400,-2,2);
    fOutput->Add(fhSubJettiness1_FJ_MIN);
    fhSubJettiness2_FJ_KT= new TH1D("fhSubJettiness2_FJ_KT","fhSubJettiness2_FJ_KT",400,-2,2);
    fOutput->Add(fhSubJettiness2_FJ_KT);
    fhSubJettiness2_FJ_MIN= new TH1D("fhSubJettiness2_FJ_MIN","fhSubJettiness2_FJ_MIN",400,-2,2);
    fOutput->Add(fhSubJettiness2_FJ_MIN);
    fh2JetTracksEtaPhiPt = new TH2D("fh3JetTracksEtaPhiPt","fh3JetTracksEtaPgiPt",Eta_Bins, Eta_Low, Eta_Up,360 , -1.5*(TMath::Pi()), 1.5*(TMath::Pi()));
    fOutput->Add(fh2JetTracksEtaPhiPt);
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
  }
  if (fJetShapeType == AliAnalysisTaskSubJetFraction::kDetEmbPart){
    fhEventCounter= new TH1F("fhEventCounter", "Event Counter", 15,0.5,15.5);
    fOutput->Add(fhEventCounter);
  }
  fOutput->Add(fTreeResponseMatrixAxis);
  TH1::AddDirectory(oldStatus);
  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
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
  if (fCentSelectOn){
    if ((fCent>fCentMax) || (fCent<fCentMin)) return 0;
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
    fhEventCounter->Fill(1);
    if(JetCont1) {
      fhEventCounter->Fill(2);
      JetCont1->ResetCurrentID();
      while((Jet1=JetCont1->GetNextAcceptJet())) {
	if (fJetShapeSub==kConstSub) JetPtThreshold=Jet1->Pt();
	else JetPtThreshold=Jet1->Pt()-(GetRhoVal(0)*Jet1->Area());
        if ( (!Jet1) || (JetPtThreshold<fPtThreshold)) continue;
        else {
	  if(fSemigoodCorrect){
	    Double_t disthole=RelativePhi(Jet1->Phi(),fHolePos);
	    if(TMath::Abs(disthole)<fHoleWidth){
	      continue;}
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
	    if (!(JetCont1->AliJetContainer::GetFractionSharedPt(Jet1)<fSharedFractionPtMin)) continue;
	    Jet3 = Jet1->ClosestJet();
	  }
	  if (!Jet3) continue;
	  Jet4=Jet3->ClosestJet();
	  if(!Jet4) continue;
	  JetsMatched=kTRUE;  
	  ///////////So at the moment tree is only filled when we have matched....is this ok?
	  if (fJetShapeSub==kConstSub) fShapesVar[0]=Jet1->Pt();
	  else fShapesVar[0]=Jet1->Pt()-(GetRhoVal(0)*Jet1->Area());
	  fShapesVar[2]=fjNSubJettiness(Jet1,0,1,0,1,0);
	  fShapesVar[4]=fjNSubJettiness(Jet1,0,2,0,1,0);
	  fShapesVar[6]=fjNSubJettiness(Jet1,0,3,0,1,0);
	  fShapesVar[8]=fjNSubJettiness(Jet1,0,2,0,1,1);
	  fShapesVar[10]=Jet1->GetNumberOfTracks();
	  if (fFullTree){
	    fShapesVar[12]=fjNSubJettiness(Jet1,0,2,0,1,2);
	    Reclusterer1 = Recluster(Jet1, 0, fSubJetRadius, fSubJetMinPt, fSubJetAlgorithm, "SubJetFinder_1");
	    fShapesVar[14]=SubJetFraction(Jet1, Reclusterer1, 1, 0, kTRUE, kFALSE);
	    fShapesVar[16]=SubJetFraction(Jet1, Reclusterer1, 2, 0, kTRUE, kFALSE);
	  }
	  if (JetsMatched){ //even needed? Not now but might be if you want to fill trees when jets aren't matched too
	    fShapesVar[1]=Jet4->Pt();
	    fShapesVar[3]=fjNSubJettiness(Jet4,3,1,0,1,0);
	    fShapesVar[5]=fjNSubJettiness(Jet4,3,2,0,1,0);
	    fShapesVar[7]=fjNSubJettiness(Jet4,3,3,0,1,0);
	    fShapesVar[9]=fjNSubJettiness(Jet4,3,2,0,1,1);
	    fShapesVar[11]=Jet4->GetNumberOfTracks();
	    if (fFullTree){
	      fShapesVar[13]=fjNSubJettiness(Jet4,3,2,0,1,2);
	      Reclusterer4=Recluster(Jet4, 3, fSubJetRadius, 0, fSubJetAlgorithm, "SubJetFinder_4");
	      fShapesVar[15]=SubJetFraction(Jet4, Reclusterer4, 1, 0, kTRUE, kFALSE);
	      fShapesVar[17]=SubJetFraction(Jet4, Reclusterer4, 2, 0, kTRUE, kFALSE);
	    } 
	  }
	  else{
	    fShapesVar[1]=-2;
	    fShapesVar[3]=-2;
	    fShapesVar[5]=-2;
	    fShapesVar[7]=-2;
	    fShapesVar[9]=-2;
	    fShapesVar[11]=-2;
	    if (fFullTree){
	      fShapesVar[13]=-2;
	      fShapesVar[15]=-2;
	      fShapesVar[17]=-2;
	    }
	  }
	  fTreeResponseMatrixAxis->Fill();
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
    fhEventCounter_1->Fill(1);
    if(JetCont1) {
      fhEventCounter_1->Fill(2); //Number of events with a jet container                                                                                               
      JetCont1->ResetCurrentID();
      while((Jet1=JetCont1->GetNextAcceptJet())) {
	if( (!Jet1) || ((Jet1->Pt())<fPtThreshold)) {
	  // fhEventCounter_1->Fill(3); //events where the jet had a problem                                                                                   
	  continue;
	}
	else {
	  if(fSemigoodCorrect){
	    Double_t disthole=RelativePhi(Jet1->Phi(),fHolePos);
	    if(TMath::Abs(disthole)<fHoleWidth){
	      continue;}
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
	  }
          else {
            fhEventCounter_2->Fill(1);
	    //continue;
          }
	  /////////////////How do you do this???? Do you only fill them for when both sets work or only for when one works??????///////////////////////
	  

	  fShapesVar[0]=Jet1->Pt();
	  fShapesVar[2]=fjNSubJettiness(Jet1,0,1,0,1,0);
	  fShapesVar[4]=fjNSubJettiness(Jet1,0,2,0,1,0);
	  fShapesVar[6]=fjNSubJettiness(Jet1,0,3,0,1,0);
	  fShapesVar[8]=fjNSubJettiness(Jet1,0,2,0,1,1);
	  fShapesVar[10]=Jet1->GetNumberOfTracks();
	  Reclusterer1 = Recluster(Jet1, 0, fSubJetRadius, fSubJetMinPt, fSubJetAlgorithm, "SubJetFinder_1");
	  if (fFullTree){
	    fShapesVar[12]=fjNSubJettiness(Jet1,0,2,0,1,2);
	    fShapesVar[14]=SubJetFraction(Jet1, Reclusterer1, 1, 0, kTRUE, kFALSE);
	    fShapesVar[16]=SubJetFraction(Jet1, Reclusterer1, 2, 0, kTRUE, kFALSE);
	  }
	  if (JetsMatched){ //even needed? Not now but might be if you want to fill trees when jets aren't matched too
	    fShapesVar[1]=Jet2->Pt();
	    fShapesVar[3]=fjNSubJettiness(Jet2,1,1,0,1,0);
	    fShapesVar[5]=fjNSubJettiness(Jet2,1,2,0,1,0);
	    fShapesVar[7]=fjNSubJettiness(Jet2,1,3,0,1,0);
	    fShapesVar[9]=fjNSubJettiness(Jet2,1,2,0,1,1);
	    fShapesVar[11]=Jet2->GetNumberOfTracks();
	    Reclusterer2 = Recluster(Jet2, 1, fSubJetRadius, 0, fSubJetAlgorithm, "SubJetFinder_2");
	    if (fFullTree){
	      fShapesVar[13]=fjNSubJettiness(Jet2,1,2,0,1,2);
	      fShapesVar[15]=SubJetFraction(Jet2, Reclusterer2, 1, 0, kTRUE, kFALSE);
	      fShapesVar[17]=SubJetFraction(Jet2, Reclusterer2, 2, 0, kTRUE, kFALSE);
	    } 
	  }
	  else{
	    fShapesVar[1]=-2;
	    fShapesVar[3]=-2;
	    fShapesVar[5]=-2;
	    fShapesVar[7]=-2;
	    fShapesVar[9]=-2;
	    fShapesVar[11]=-2;
	    if (fFullTree){
	      fShapesVar[13]=-2;
	      fShapesVar[15]=-2;
	      fShapesVar[17]=-2;
	    }
	  }
	  fTreeResponseMatrixAxis->Fill();

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
	  if(fSemigoodCorrect){
	    Double_t disthole=RelativePhi(Jet1->Phi(),fHolePos);
	    if(TMath::Abs(disthole)<fHoleWidth){
	      continue;}
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
	  if(fJetShapeSub==kNoSub || fJetShapeSub==kDerivSub) fShapesVar[0]= Jet1->Pt()-(GetRhoVal(0)*Jet1->Area());
	  else fShapesVar[0]=Jet1->Pt();
	  fShapesVar[2]=fjNSubJettiness(Jet1,0,1,0,1,0);
	  fShapesVar[4]=fjNSubJettiness(Jet1,0,2,0,1,0);
	  fShapesVar[6]=fjNSubJettiness(Jet1,0,3,0,1,0);
	  fShapesVar[8]=fjNSubJettiness(Jet1,0,2,0,1,1);
	  fShapesVar[10]=Jet1->GetNumberOfTracks();
	  AliEmcalJetFinder *Reclusterer1 = Recluster(Jet1, 0, fSubJetRadius, fSubJetMinPt, fSubJetAlgorithm, "SubJetFinder");
	  if (fFullTree){
	    fShapesVar[12]=fjNSubJettiness(Jet1,0,2,0,1,2);
	    fShapesVar[14]=SubJetFraction(Jet1, Reclusterer1, 1, 0, kTRUE, kFALSE);
	    fShapesVar[16]=SubJetFraction(Jet1, Reclusterer1, 2, 0, kTRUE, kFALSE);
	  }
	  fShapesVar[1]=-2;
	  fShapesVar[3]=-2;
	  fShapesVar[5]=-2;
	  fShapesVar[7]=-2;
	  fShapesVar[9]=-2;
	  fShapesVar[11]=-2;
	  if (fFullTree){
	    fShapesVar[13]=-2;
	    fShapesVar[15]=-2;
	    fShapesVar[17]=-2;
	  }
	  fTreeResponseMatrixAxis->Fill();
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
	  if(fSemigoodCorrect){
	    Double_t disthole=RelativePhi(Jet1->Phi(),fHolePos);
	    if(TMath::Abs(disthole)<fHoleWidth){
	      continue;}
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

	  Int_t FJ_Algorithm=9; //KT=0, CA=1, AKT_R02=2, WTA_KT=3, WTA_CA=4, OP_KT=5, OP_CA=6, OP_AKT_R02=7, OP_WTA_KT=8, OP_WTA_CA=9, MIN=10
	  Double_t FJ_Beta=1.0;
	  fhSubJettiness1CheckRatio_FJ_KT->Fill(fjNSubJettiness(Jet1,0,1,0,FJ_Beta,0)-fjNSubJettiness(Jet1,0,1,0,FJ_Beta,0),fjNSubJettiness(Jet1,0,1,0,FJ_Beta,0));
	  fhSubJettiness1CheckRatio_FJ_CA->Fill(fjNSubJettiness(Jet1,0,1,0,FJ_Beta,0)-fjNSubJettiness(Jet1,0,1,1,FJ_Beta,0),fjNSubJettiness(Jet1,0,1,0,FJ_Beta,0));
	  fhSubJettiness1CheckRatio_FJ_AKT->Fill(fjNSubJettiness(Jet1,0,1,0,FJ_Beta,0)-fjNSubJettiness(Jet1,0,1,2,FJ_Beta,0),fjNSubJettiness(Jet1,0,1,0,FJ_Beta,0));
	  fhSubJettiness1CheckRatio_FJ_WTA_KT->Fill(fjNSubJettiness(Jet1,0,1,0,FJ_Beta,0)-fjNSubJettiness(Jet1,0,1,3,FJ_Beta,0),fjNSubJettiness(Jet1,0,1,0,FJ_Beta,0));
	  fhSubJettiness1CheckRatio_FJ_WTA_CA->Fill(fjNSubJettiness(Jet1,0,1,0,FJ_Beta,0)-fjNSubJettiness(Jet1,0,1,4,FJ_Beta,0),fjNSubJettiness(Jet1,0,1,0,FJ_Beta,0));
	  fhSubJettiness1CheckRatio_FJ_OP_KT->Fill(fjNSubJettiness(Jet1,0,1,0,FJ_Beta,0)-fjNSubJettiness(Jet1,0,1,5,FJ_Beta,0),fjNSubJettiness(Jet1,0,1,0,FJ_Beta,0));
	  fhSubJettiness1CheckRatio_FJ_OP_CA->Fill(fjNSubJettiness(Jet1,0,1,0,FJ_Beta,0)-fjNSubJettiness(Jet1,0,1,6,FJ_Beta,0),fjNSubJettiness(Jet1,0,1,0,FJ_Beta,0));
	  fhSubJettiness1CheckRatio_FJ_OP_AKT->Fill(fjNSubJettiness(Jet1,0,1,0,FJ_Beta,0)-fjNSubJettiness(Jet1,0,1,7,FJ_Beta,0),fjNSubJettiness(Jet1,0,1,0,FJ_Beta,0));
	  fhSubJettiness1CheckRatio_FJ_OP_WTA_KT->Fill(fjNSubJettiness(Jet1,0,1,0,FJ_Beta,0)-fjNSubJettiness(Jet1,0,1,8,FJ_Beta,0),fjNSubJettiness(Jet1,0,1,0,FJ_Beta,0));
	  fhSubJettiness1CheckRatio_FJ_OP_WTA_CA->Fill(fjNSubJettiness(Jet1,0,1,0,FJ_Beta,0)-fjNSubJettiness(Jet1,0,1,9,FJ_Beta,0),fjNSubJettiness(Jet1,0,1,0,FJ_Beta,0));
	  fhSubJettiness1CheckRatio_FJ_MIN->Fill(fjNSubJettiness(Jet1,0,1,0,FJ_Beta,0)-fjNSubJettiness(Jet1,0,1,10,FJ_Beta,0),fjNSubJettiness(Jet1,0,1,0,FJ_Beta,0));
	  fhSubJettiness2CheckRatio_FJ_KT->Fill(fjNSubJettiness(Jet1,0,2,0,FJ_Beta,0)-fjNSubJettiness(Jet1,0,2,0,FJ_Beta,0),fjNSubJettiness(Jet1,0,2,0,FJ_Beta,0));
	  fhSubJettiness2CheckRatio_FJ_CA->Fill(fjNSubJettiness(Jet1,0,2,0,FJ_Beta,0)-fjNSubJettiness(Jet1,0,2,1,FJ_Beta,0),fjNSubJettiness(Jet1,0,2,0,FJ_Beta,0));
	  fhSubJettiness2CheckRatio_FJ_AKT->Fill(fjNSubJettiness(Jet1,0,2,0,FJ_Beta,0)-fjNSubJettiness(Jet1,0,2,2,FJ_Beta,0),fjNSubJettiness(Jet1,0,2,0,FJ_Beta,0));
	  fhSubJettiness2CheckRatio_FJ_WTA_KT->Fill(fjNSubJettiness(Jet1,0,2,0,FJ_Beta,0)-fjNSubJettiness(Jet1,0,2,3,FJ_Beta,0),fjNSubJettiness(Jet1,0,2,0,FJ_Beta,0));
	  fhSubJettiness2CheckRatio_FJ_WTA_CA->Fill(fjNSubJettiness(Jet1,0,2,0,FJ_Beta,0)-fjNSubJettiness(Jet1,0,2,4,FJ_Beta,0),fjNSubJettiness(Jet1,0,2,0,FJ_Beta,0));
	  fhSubJettiness2CheckRatio_FJ_OP_KT->Fill(fjNSubJettiness(Jet1,0,2,0,FJ_Beta,0)-fjNSubJettiness(Jet1,0,2,5,FJ_Beta,0),fjNSubJettiness(Jet1,0,2,0,FJ_Beta,0));
	  fhSubJettiness2CheckRatio_FJ_OP_CA->Fill(fjNSubJettiness(Jet1,0,2,0,FJ_Beta,0)-fjNSubJettiness(Jet1,0,2,6,FJ_Beta,0),fjNSubJettiness(Jet1,0,2,0,FJ_Beta,0));
	  fhSubJettiness2CheckRatio_FJ_OP_AKT->Fill(fjNSubJettiness(Jet1,0,2,0,FJ_Beta,0)-fjNSubJettiness(Jet1,0,2,7,FJ_Beta,0),fjNSubJettiness(Jet1,0,2,0,FJ_Beta,0));
	  fhSubJettiness2CheckRatio_FJ_OP_WTA_KT->Fill(fjNSubJettiness(Jet1,0,2,0,FJ_Beta,0)-fjNSubJettiness(Jet1,0,2,8,FJ_Beta,0),fjNSubJettiness(Jet1,0,2,0,FJ_Beta,0));
	  fhSubJettiness2CheckRatio_FJ_OP_WTA_CA->Fill(fjNSubJettiness(Jet1,0,2,0,FJ_Beta,0)-fjNSubJettiness(Jet1,0,2,9,FJ_Beta,0),fjNSubJettiness(Jet1,0,2,0,FJ_Beta,0));
	  fhSubJettiness2CheckRatio_FJ_MIN->Fill(fjNSubJettiness(Jet1,0,2,0,FJ_Beta,0)-fjNSubJettiness(Jet1,0,2,10,FJ_Beta,0),fjNSubJettiness(Jet1,0,2,0,FJ_Beta,0));
 
	  fhSubJettiness1_FJ_KT->Fill(fjNSubJettiness(Jet1,0,1,0,FJ_Beta,0));
	  fhSubJettiness1_FJ_CA->Fill(fjNSubJettiness(Jet1,0,1,1,FJ_Beta,0));
	  fhSubJettiness1_FJ_AKT->Fill(fjNSubJettiness(Jet1,0,1,2,FJ_Beta,0));
	  fhSubJettiness1_FJ_WTA_KT->Fill(fjNSubJettiness(Jet1,0,1,3,FJ_Beta,0));
	  fhSubJettiness1_FJ_WTA_CA->Fill(fjNSubJettiness(Jet1,0,1,4,FJ_Beta,0));
	  fhSubJettiness1_FJ_OP_KT->Fill(fjNSubJettiness(Jet1,0,1,5,FJ_Beta,0));
	  fhSubJettiness1_FJ_OP_CA->Fill(fjNSubJettiness(Jet1,0,1,6,FJ_Beta,0));
	  fhSubJettiness1_FJ_OP_AKT->Fill(fjNSubJettiness(Jet1,0,1,7,FJ_Beta,0));
	  fhSubJettiness1_FJ_OP_WTA_KT->Fill(fjNSubJettiness(Jet1,0,1,8,FJ_Beta,0));
	  fhSubJettiness1_FJ_OP_WTA_CA->Fill(fjNSubJettiness(Jet1,0,1,9,FJ_Beta,0));
	  fhSubJettiness1_FJ_MIN->Fill(fjNSubJettiness(Jet1,0,1,10,FJ_Beta,0));
	  fhSubJettiness2_FJ_KT->Fill(fjNSubJettiness(Jet1,0,2,0,FJ_Beta,0));
	  fhSubJettiness2_FJ_CA->Fill(fjNSubJettiness(Jet1,0,2,1,FJ_Beta,0));
	  fhSubJettiness2_FJ_AKT->Fill(fjNSubJettiness(Jet1,0,2,2,FJ_Beta,0)); //because in this case we aren't garantueed to get 2 subjets
	  fhSubJettiness2_FJ_WTA_KT->Fill(fjNSubJettiness(Jet1,0,2,3,FJ_Beta,0));
	  fhSubJettiness2_FJ_WTA_CA->Fill(fjNSubJettiness(Jet1,0,2,4,FJ_Beta,0));
	  fhSubJettiness2_FJ_OP_KT->Fill(fjNSubJettiness(Jet1,0,2,5,FJ_Beta,0));
	  fhSubJettiness2_FJ_OP_CA->Fill(fjNSubJettiness(Jet1,0,2,6,FJ_Beta,0));
	  fhSubJettiness2_FJ_OP_AKT->Fill(fjNSubJettiness(Jet1,0,2,7,FJ_Beta,0));
	  fhSubJettiness2_FJ_OP_WTA_KT->Fill(fjNSubJettiness(Jet1,0,2,8,FJ_Beta,0));
	  fhSubJettiness2_FJ_OP_WTA_CA->Fill(fjNSubJettiness(Jet1,0,2,9,FJ_Beta,0));
	  fhSubJettiness2_FJ_MIN->Fill(fjNSubJettiness(Jet1,0,2,10,FJ_Beta,0));

	  if(fJetShapeSub==kNoSub || fJetShapeSub==kDerivSub) fShapesVar[0]= Jet1->Pt()-(GetRhoVal(0)*Jet1->Area());
	  else fShapesVar[0]=Jet1->Pt(); 
	  fShapesVar[2]=fjNSubJettiness(Jet1,0,1,0,1,0);
	  fShapesVar[4]=fjNSubJettiness(Jet1,0,2,0,1,0);
	  fShapesVar[6]=fjNSubJettiness(Jet1,0,3,0,1,0);
	  fShapesVar[8]=fjNSubJettiness(Jet1,0,2,0,1,1);
	  fShapesVar[10]=Jet1->GetNumberOfTracks();
	  AliEmcalJetFinder *Reclusterer1 = Recluster(Jet1, 0, fSubJetRadius, fSubJetMinPt, fSubJetAlgorithm, "SubJetFinder");
	  if (fFullTree){
	    fShapesVar[12]=fjNSubJettiness(Jet1,0,2,0,1,2);
	    fShapesVar[14]=SubJetFraction(Jet1, Reclusterer1, 1, 0, kTRUE, kFALSE);
	    fShapesVar[16]=SubJetFraction(Jet1, Reclusterer1, 2, 0, kTRUE, kFALSE);
	  }
	  fShapesVar[1]=-2;
	  fShapesVar[3]=-2;
	  fShapesVar[5]=-2;
	  fShapesVar[7]=-2;
	  fShapesVar[9]=-2;
	  fShapesVar[11]=-2;
	  if (fFullTree){
	    fShapesVar[13]=-2;
	    fShapesVar[15]=-2;
	    fShapesVar[17]=-2;
	  }
	  fTreeResponseMatrixAxis->Fill();
	  
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
  return kTRUE;
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
  const AliVVertex *vert = InputEvent()->GetPrimaryVertex();
  Double_t dVtx[3]={vert->GetX(),vert->GetY(),vert->GetZ()};
  if(Reclusterer->AliEmcalJetFinder::Filter(Jet, JetCont, dVtx)){;}  //reclustering jet1 using the jetfinderobject Reclusterer                            
  return Reclusterer;
}









//----------------------------------------------------------------------
Double_t AliAnalysisTaskSubJetFraction::SubJetOrdering(AliEmcalJet *Jet, AliEmcalJetFinder *Reclusterer, Int_t N, Int_t Type, Bool_t Index){
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
Double_t AliAnalysisTaskSubJetFraction::fjNSubJettiness(AliEmcalJet *Jet, Int_t JetContNb,Int_t N, Int_t Algorithm, Double_t Beta, Int_t Option){

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
    else{
      AliJetContainer *JetCont = GetJetContainer(JetContNb);
      AliEmcalJetFinder *JetFinder=new AliEmcalJetFinder("Nsubjettiness");
      JetFinder->SetJetMaxEta(0.9-fJetRadius);
      JetFinder->SetRadius(fJetRadius); 
      JetFinder->SetJetAlgorithm(0); //0 for anti-kt     1 for kt  //this is for the JET!!!!!!!!!! Not the SubJets
      JetFinder->SetRecombSheme(0);
      JetFinder->SetJetMinPt(Jet->Pt());
      const AliVVertex *vert = InputEvent()->GetPrimaryVertex();
      Double_t dVtx[3]={vert->GetX(),vert->GetY(),vert->GetZ()};
      return JetFinder->Nsubjettiness(Jet,JetCont,dVtx,N,Algorithm,fSubJetRadius,Beta,Option);
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
    fhJetAngularity->Scale(100/(fhEventCounter->GetBinContent(4)));
    fhJetPTD->Scale(100/(fhEventCounter->GetBinContent(4)));
    fhNumberOfJetTracks->Scale(1.0/(fhEventCounter->GetBinContent(4)));
    fhSubJetPtFrac->Scale(100/(fhEventCounter->GetBinContent(8)));
    fhSubJetPtFrac2->Scale(100/(fhEventCounter->GetBinContent(9)));
    fhSubJetEnergyFrac->Scale(100/(fhEventCounter->GetBinContent(8)));
    fhSubJetEnergyFrac2->Scale(100/(fhEventCounter->GetBinContent(9)));
    fhJetCounter->Scale(1.0/(fhEventCounter->GetBinContent(1)));  //is the first bin the correct one to look at?
   

 fhSubJetCounter->Scale(1.0/(fhEventCounter->GetBinContent(5)));
    fhSubJetiness1->Scale(100/(fhEventCounter->GetBinContent(12)));
    fhSubJetiness2->Scale(100/(fhEventCounter->GetBinContent(12)));
    fh2to1SubJetinessRatio->Scale(100/(fhEventCounter->GetBinContent(12)));

    */
    /*  

  fhJetPt->Scale(1.0/(fhEventCounter->GetBinContent(1)));                                      
  fhSubJetPt->Scale(1.0/(fhEventCounter->GetBinContent(5)));                                
  fhJetMass->Scale(1.0/(fhEventCounter->GetBinContent(1)));                                                                                                                     
  fhJetPhi->Scale(1.0/(fhEventCounter->GetBinContent(1)));                                                                                               
  fhJetEta->Scale(1.0/(fhEventCounter->GetBinContent(1)));                                                                                               
  fhJetRadius->Scale(1.0/(fhEventCounter->GetBinContent(4)));  //should this and JetAngularity be divided by Bin 1 or 4????                                                        
  fhJetAngularity->Scale(1.0/(fhEventCounter->GetBinContent(4)));                                                                                                                  
  fhJetPTD->Scale(1.0/(fhEventCounter->GetBinContent(4)));                                                                                                                         
  fhNumberOfJetTracks->Scale(1.0/(fhEventCounter->GetBinContent(4)));                                                                                                              
  fhSubJetPtFrac->Scale(1.0/(fhEventCounter->GetBinContent(8)));                                                                                                                   
  fhSubJetPtFrac2->Scale(1.0/(fhEventCounter->GetBinContent(9)));                                                                                                                  
  fhSubJetEnergyFrac->Scale(1.0/(fhEventCounter->GetBinContent(8)));                                                                                                               
  fhSubJetEnergyFrac2->Scale(1.0/(fhEventCounter->GetBinContent(9)));                                                                                                              
  fhJetCounter->Scale(1.0/(fhEventCounter->GetBinContent(1)));  //is the first bin the correct one to look at?                                                                     
  fhSubJetCounter->Scale(1.0/(fhEventCounter->GetBinContent(5)));                                                                                                            
  fhSubJetiness1->Scale(1.0/(fhEventCounter->GetBinContent(12)));                                                                                                                  
  fhSubJetiness2->Scale(1.0/(fhEventCounter->GetBinContent(12)));                                                                                                                  
  fh2to1SubJetinessRatio->Scale(1.0/(fhEventCounter->GetBinContent(12)));            
    */
  }
  
}
