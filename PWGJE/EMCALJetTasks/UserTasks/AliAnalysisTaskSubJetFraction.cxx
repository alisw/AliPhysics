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
#include "AliPythiaInfo.h"
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
  fhJetPt(0x0),
  fhJetPhi(0x0),
  fhJetEta(0x0),
  fhJetMass(0x0),
  fhJetRadius(0x0),
  fhJetAngularity(0x0),
  fhJetPTD(0x0),
  fhJetCounter(0x0),
  fhNumberOfJetTracks(0x0),
  fhSubJetPt(0x0),
  fhSubJetMass(0x0),
  fhSubJetRadius(0x0),
  fhSubJetCounter(0x0),
  fhNumberOfSubJetTracks(0x0),
  fhSubJetPtFrac(0x0),
  fhSubJetPtFrac2(0x0),
  fhSubJetPtLoss(0x0),
  fhSubJetPtLoss2(0x0),
  fhSubJetEnergyFrac(0x0),
  fhSubJetEnergyFrac2(0x0),
  fhSubJetEnergyLoss(0x0),
  fhSubJetEnergyLoss2(0x0),
  fhPtRatio(0x0),
  fhParticleSubJetPtFrac(0x0),
  fhDetectorSubJetPtFrac(0x0),
  fhParticleSubJetPtFrac2(0x0),
  fhDetectorSubJetPtFrac2(0x0),
  fhSubJetPtFracRatio(0x0),
  fhSubJetPtFrac2Ratio(0x0),
  fhParticleSubJetiness1(0x0),
  fhDetectorSubJetiness1(0x0),
  fhSubJetiness1Ratio(0x0),
  fhSubJetiness1(0x0),
  fhParticleSubJetiness2(0x0),
  fhDetectorSubJetiness2(0x0),
  fhSubJetiness2Ratio(0x0),
  fhSubJetiness2(0x0),
  fh2to1SubJetinessRatio(0x0),
  fhEventCounter(0x0),
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
  fhJetPt(0x0),
  fhJetPhi(0x0),
  fhJetEta(0x0),
  fhJetMass(0x0),
  fhJetRadius(0x0),
  fhJetAngularity(0x0),
  fhJetPTD(0x0),
  fhJetCounter(0x0),
  fhNumberOfJetTracks(0x0),
  fhSubJetPt(0x0),
  fhSubJetMass(0x0),
  fhSubJetRadius(0x0),
  fhSubJetCounter(0x0),
  fhNumberOfSubJetTracks(0x0),
  fhSubJetPtFrac(0x0),
  fhSubJetPtFrac2(0x0),
  fhSubJetPtLoss(0x0),
  fhSubJetPtLoss2(0x0),
  fhSubJetEnergyFrac(0x0),
  fhSubJetEnergyFrac2(0x0),
  fhSubJetEnergyLoss(0x0),
  fhSubJetEnergyLoss2(0x0),
  fhPtRatio(0x0),
  fhParticleSubJetPtFrac(0x0),
  fhDetectorSubJetPtFrac(0x0),
  fhParticleSubJetPtFrac2(0x0),
  fhDetectorSubJetPtFrac2(0x0),
  fhSubJetPtFracRatio(0x0),
  fhSubJetPtFrac2Ratio(0x0),
  fhParticleSubJetiness1(0x0),
  fhDetectorSubJetiness1(0x0),
  fhSubJetiness1Ratio(0x0),
  fhSubJetiness1(0x0),
  fhParticleSubJetiness2(0x0),
  fhDetectorSubJetiness2(0x0),
  fhSubJetiness2Ratio(0x0),
  fhSubJetiness2(0x0),
  fh2to1SubJetinessRatio(0x0),
  fhEventCounter(0x0),
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
  const Int_t nVar = 10;
  fShapesVar = new Double_t [nVar]; //shapes used for tagging   
  fTreeResponseMatrixAxis = new TTree("fTreeJetShape", "fTreeJetShape");
  TString *fShapesVarNames = new TString [nVar];
  
  fShapesVarNames[0] = "Pt_Particle_Level";
  fShapesVarNames[1] = "Frac_Particle_Level";
  fShapesVarNames[2] = "Frac2_Particle_Level";
  fShapesVarNames[3] = "Pt_Detector_Level";
  fShapesVarNames[4] = "Frac_Detector_Level";
  fShapesVarNames[5] = "Frac2_Detector_Level";
  fShapesVarNames[6] = "1Subjetiness_Particle_Level";
  fShapesVarNames[7] = "1Subjetiness_Detector_Level";
  fShapesVarNames[8] = "2Subjetiness_Particle_Level";
  fShapesVarNames[9] = "2Subjetiness_Detector_Level";
  for(Int_t ivar=0; ivar < nVar; ivar++){
    cout<<"looping over variables"<<endl;
    fTreeResponseMatrixAxis->Branch(fShapesVarNames[ivar].Data(), &fShapesVar[ivar], Form("%s/D", fShapesVarNames[ivar].Data()));
  }
  
  
  
  if (fJetShapeType==AliAnalysisTaskSubJetFraction::kData){
    
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
    fhJetAngularity= new TH1F("fhJetAngularity", "Jet Angularity", 400, -2,2);
    fOutput->Add(fhJetAngularity);
    fhJetPTD= new TH1F("fhJetPTD", "Jet PTD", 400, -2,2);
    fOutput->Add(fhJetPTD);
    fhNumberOfJetTracks= new TH1F("fhNumberOfJetTracks", "Number of Tracks within a Jet", 30, -0.5,29.5);
    fOutput->Add(fhNumberOfJetTracks);
    fhSubJetRadius= new TH1F("fhSubJetRadius", "SubJet Radius", 100, -0.05,0.995);
    fOutput->Add(fhSubJetRadius);
    fhSubJetPt= new TH1F("fhSubJetPt", "SubJet Pt", 1500, -0.5,149.5);
    fOutput->Add(fhSubJetPt);
    fhSubJetMass= new TH1F("fhSubJetMass", "Sub Jet Mass", 4000,-0.5, 39.5);
    fOutput->Add(fhSubJetMass);
    fhNumberOfSubJetTracks= new TH1F("fhNumberOfSubJetTracks", "Number of Tracks within a Sub Jet", 30, -0.5,29.5);
    fOutput->Add(fhNumberOfSubJetTracks);
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
    fhSubJetiness2= new TH1F("fhSubjetiness2", "Tau 2 value",101, -0.05,1.05);
    fOutput->Add(fhSubJetiness2);
    fhJetCounter= new TH1F("fhJetCounter", "Jet Counter", 100, -0.5, 99.5);
    fOutput->Add(fhJetCounter);
    fhSubJetCounter = new TH1F("fhSubJetCounter", "SubJet Counter",50, -0.5,49.5);
    fOutput->Add(fhSubJetCounter);
    fh2to1SubJetinessRatio= new TH1F("fh2to1SubJetinessRatio", "Ratio of #tau 1 to #tau 2",200, -0.5,1.5);
    fOutput->Add(fh2to1SubJetinessRatio);
  }
  if(fJetShapeType==AliAnalysisTaskSubJetFraction::kTrueDet){
    fhPtRatio= new TH1F("fhPtRatio", "MC pt ratio for particle and detector level jets",100, -0.5,9.5);
    fOutput->Add(fhPtRatio);
    fhParticleSubJetPtFrac= new TH1F("fhParticleSubJetPtFrac", "Pt Fraction of Highest Pt Subjet compared to original Jet for MC Particle Level data",101, -0.05,1.05);
    fOutput->Add(fhParticleSubJetPtFrac);
    fhParticleSubJetPtFrac2= new TH1F("fhParticleSubJetPtFrac2", "Pt Fraction of Two Highest Pt Subjets compared to original Jet for MC Particle Level data",101, -0.05,1.05);
    fOutput->Add(fhParticleSubJetPtFrac2);
    fhDetectorSubJetPtFrac= new TH1F("fhDetectorSubJetPtFrac", "Pt Fraction of Highest Pt Subjet compared to original Jet for MC Detector Level data",101, -0.05,1.05);
    fOutput->Add(fhDetectorSubJetPtFrac);
    fhDetectorSubJetPtFrac2= new TH1F("fhDetectorSubJetPtFrac2", "Pt Fraction of Two Highest Pt Subjets compared to original Jet for MC Detector Level data",101, -0.05,1.05);
    fOutput->Add(fhDetectorSubJetPtFrac2);
    fhSubJetPtFracRatio= new TH1F("fhSubJetPtFracRatio", "Ratio of Pt Fraction of Highest Pt Subjet compared to original Jet for MC particle and Detector Level data", 1010,-0.05,10.5); 
    fOutput->Add(fhSubJetPtFracRatio);
    fhSubJetPtFrac2Ratio= new TH1F("fhSubJetPtFrac2Ratio", "Ratio of Pt Fraction of Two Highest Pt Subjets compared to original Jet for MC particle and Detector Level data", 1010,-0.05,10.5);
    fOutput->Add(fhSubJetPtFrac2Ratio);
    fhParticleSubJetiness1= new TH1F("fhParticleSubJetiness1", "1Tau  for MC Particle Level data",101, -0.05,1.05);
    fOutput->Add(fhParticleSubJetiness1);
    fhDetectorSubJetiness1= new TH1F("fhDetectorSubJetiness1", "1Tau for MC Detector Level data",101, -0.05,1.05);
    fOutput->Add(fhDetectorSubJetiness1);
    fhSubJetiness1Ratio= new TH1F("fhSubJetiness1Ratio", "1Tau Ratio for MC Particle Level and Detector Level data",2010, -10.05,10.05);
    fOutput->Add(fhSubJetiness1Ratio);
    fhParticleSubJetiness2= new TH1F("fhParticleSubJetiness2", "2Tau  for MC Particle Level data",101, -0.05,1.05);
    fOutput->Add(fhParticleSubJetiness2);
    fhDetectorSubJetiness2= new TH1F("fhDetectorSubJetiness2", "2Tau for MC Detector Level data",101, -0.05,1.05);
    fOutput->Add(fhDetectorSubJetiness2);
    fhSubJetiness2Ratio= new TH1F("fhSubJetiness2Ratio", "2Tau Ratio for MC Particle Level and Detector Level data",2010, -10.05,10.05);
    fOutput->Add(fhSubJetiness2Ratio);
  }
  fOutput->Add(fTreeResponseMatrixAxis);
  fhEventCounter= new TH1F("fhEventCounter", "Event Counter", 15,0.5,15.5);
  fOutput->Add(fhEventCounter);
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


  //Jet1 -> KData Jet
  //Jet2 -> KData SubJet
  //Jet3 -> KTrueDet MC Particle level Jet
  //Jet4 -> KTrueDet MC Detector level Jet
  //Jet5 -> KTrueDet MC Particle level SubJet
  //Jet6 -> KTrueDet MC Detector level SubJet
  const AliVVertex *vert = InputEvent()->GetPrimaryVertex();
  Double_t dVtx[3]={vert->GetX(),vert->GetY(),vert->GetZ()};
  if(vert) fhEventCounter->Fill(7);


  if (fJetShapeType == AliAnalysisTaskSubJetFraction::kTrueDet){
    AliJetContainer *JetContParticle = GetJetContainer(0);
    AliJetContainer *JetContDetector = GetJetContainer(1);
    AliEmcalJet *Jet3=NULL; //MC particle level jet
    AliEmcalJet *Jet4=NULL; //MC detector level jet
    AliEmcalJet *Jet5=NULL; //MC particle level subjet                                                                                                                             
    AliEmcalJet *Jet6=NULL; //MC detector level subjet  
    AliEmcalJetFinder *ReclustererParticle = new AliEmcalJetFinder("SubJetFinder2"); //JetFinder Object for reclustered jets from MC particle Level data                   
    ReclustererParticle->SetRadius(fSubJetRadius);
    ReclustererParticle->SetJetMinPt(fSubJetMinPt);
    ReclustererParticle->SetJetAlgorithm(fSubJetAlgorithm); //0 for anti-kt     1 for kt  
    AliEmcalJetFinder *ReclustererDetector = new AliEmcalJetFinder("SubJetFinder3"); //JetFinder Object for reclustered jets from MC detector Level data                        
    ReclustererDetector->SetRadius(fSubJetRadius);
    ReclustererDetector->SetJetMinPt(fSubJetMinPt);
    ReclustererDetector->SetJetAlgorithm(fSubJetAlgorithm); //0 for anti-kt     1 for kt           
    AliVParticle *ParticleJetParticle = 0x0; //Individual constituent tracks (particles) in a jet   
    AliVParticle *DetectorJetParticle = 0x0; //Individual constituent tracks (particles) in a jet   
    Int_t SubJetCounterParticle=0;
    Int_t SubJetCounterDetector=0;
    Int_t ParticleReclusterOk=0; //makes sure clustering has happened before commiting values to the response matrix
    Int_t DetectorReclusterOk=0;
    Double_t HighestParticleSubJetPt=-1;
    Double_t NextHighestParticleSubJetPt=-1;
    Double_t HighestDetectorSubJetPt=-1;
    Double_t NextHighestDetectorSubJetPt=-1;
    Double_t ParticleSubJetPtFrac=0;
    Double_t DetectorSubJetPtFrac=0;
    Double_t ParticleSubJetPtFrac2=0;
    Double_t DetectorSubJetPtFrac2=0;

    Int_t HardestParticleSubJetIndex=-1;
    Int_t NextHardestParticleSubJetIndex=-1;
    Double_t HardestParticleSubJetEta=-10;
    Double_t NextHardestParticleSubJetEta=-10;
    Double_t HardestParticleSubJetPhi=-10;
    Double_t NextHardestParticleSubJetPhi=-10;
    Double_t SubJetiness1_Numerator_Particle=0;
    Double_t SubJetiness1_Denominator_Particle=0;
    Double_t SubJetiness2_Numerator_Particle=0;
    Double_t SubJetiness2_Denominator_Particle=0;
    Double_t DeltaR_Particle=0;
    Double_t DeltaR1_Particle=0;
    Double_t DeltaR2_Particle=0;
    Double_t TempSubJetCounterParticle=-1;

    Int_t HardestDetectorSubJetIndex=-1;
    Int_t NextHardestDetectorSubJetIndex=-1;
    Double_t HardestDetectorSubJetEta=-10;
    Double_t NextHardestDetectorSubJetEta=-10;
    Double_t HardestDetectorSubJetPhi=-10;
    Double_t NextHardestDetectorSubJetPhi=-10;
    Double_t SubJetiness1_Numerator_Detector=0;
    Double_t SubJetiness1_Denominator_Detector=0;
    Double_t SubJetiness2_Numerator_Detector=0;
    Double_t SubJetiness2_Denominator_Detector=0;
    Double_t DeltaR_Detector=0;
    Double_t DeltaR1_Detector=0;
    Double_t DeltaR2_Detector=0;
    Double_t TempSubJetCounterDetector=-1;
    if(JetContParticle){
      JetContParticle->ResetCurrentID();
      //      while((Jet3=JetContParticle->GetNextAcceptJet()) && ((Jet3->GetNumberOfTracks())>1) && (Jet3->Pt()>=10)) {      //exludes one or zero(?) track jets           
      while((Jet3=JetContParticle->GetNextAcceptJet()) &&  (Jet3->Pt()>=10)) {      //exludes one or zero(?) track jets
	if(!Jet3) {
	  continue;
	}
	else {
	  if((Jet4 = (Jet3->ClosestJet()))){
	    if ((Jet4->Pt())!=0) fhPtRatio->Fill(Jet3->Pt()/Jet4->Pt());
	    HighestParticleSubJetPt=-1;
	    HighestDetectorSubJetPt=-1;
	    TempSubJetCounterParticle=0;
	    ParticleReclusterOk=0;
	    if(ReclustererParticle->AliEmcalJetFinder::Filter(Jet3, JetContParticle, dVtx)){  //reclustering jet1 using the jetfinderobject Reclusterer                
	      ParticleReclusterOk=1;
	      SubJetCounterParticle=ReclustererParticle->GetNumberOfJets(); // Number of reclustered SubJets in each original jet                                            
	      for (Int_t i=0; i<SubJetCounterParticle; i++){ //Loops through each SubJet in a reclustered jet                                                              
		//	if((Jet5=ReclustererParticle->GetJet(i)) && Jet5->GetNumberOfTracks()>1){ //jet2 is now set to the Subjet                                              
		if((Jet5=ReclustererParticle->GetJet(i))){  
		  if((Jet5->Pt())>HighestParticleSubJetPt){ //This finds the highest pt subjet in each jet
		    NextHighestParticleSubJetPt=HighestParticleSubJetPt;                                                                       
		    HighestParticleSubJetPt=Jet5->Pt();
		    NextHardestParticleSubJetIndex=HardestParticleSubJetIndex;
		    HardestParticleSubJetIndex=i; //NSUBJETINESS              
		    NextHardestParticleSubJetEta=HardestParticleSubJetEta;                                                                                                      
		    HardestParticleSubJetEta=Jet5->Eta(); //NSUBJETINESS  
		    NextHardestParticleSubJetPhi=HardestParticleSubJetPhi;                                                                                                          
		    HardestParticleSubJetPhi=Jet5->Phi(); //NSUBJETINESS 
		  }
		  else if((Jet5->Pt())>NextHighestParticleSubJetPt){ //This finds the 2nd highest pt subjet in each jet                                                       
		    NextHighestParticleSubJetPt=Jet5->Pt();
		    NextHardestParticleSubJetIndex=i; //NSUBJETINESS                                                                                       
		    NextHardestParticleSubJetEta=Jet5->Eta(); //NSUBJETINESS                                                                                                       
		    NextHardestParticleSubJetPhi=Jet5->Phi(); //NSUBJETINESS                                                                                                   
		  }
		}
		else{TempSubJetCounterParticle++;}
	      }
	      SubJetCounterParticle=SubJetCounterParticle-TempSubJetCounterParticle;
	      if (SubJetCounterParticle>=1){
		ParticleSubJetPtFrac=HighestParticleSubJetPt/(Jet3->Pt());
		fhParticleSubJetPtFrac->Fill(ParticleSubJetPtFrac); //Pt fraction of highest Pt subjet compared to original jet                            
	      }
	      if (SubJetCounterParticle>=2){
		ParticleSubJetPtFrac2=(HighestParticleSubJetPt+NextHighestParticleSubJetPt)/(Jet3->Pt());
		fhParticleSubJetPtFrac2->Fill(ParticleSubJetPtFrac2); //Pt fraction of two highest Pt subjets compared to original jet
	      }
	      if(SubJetCounterParticle>=2){
		SubJetiness1_Numerator_Particle=0;
		SubJetiness1_Denominator_Particle=0;	
		SubJetiness2_Numerator_Particle=0;
		SubJetiness2_Denominator_Particle=0;
		DeltaR_Particle=0;
		DeltaR1_Particle=0;
		DeltaR2_Particle=0;
		for (Int_t i=0; i< (Jet3->GetNumberOfTracks()); i++){  //loops through all tracks (particles in the jet                                                          
		  ParticleJetParticle = static_cast<AliVParticle*>(Jet3->TrackAt(i, JetContParticle->GetParticleContainer()->GetArray()));
		  DeltaR1_Particle=TMath::Sqrt((((ParticleJetParticle->Eta())-HardestParticleSubJetEta)*((ParticleJetParticle->Eta())- HardestParticleSubJetEta))+((RelativePhi(HardestParticleSubJetPhi,ParticleJetParticle->Phi()))*(RelativePhi(HardestParticleSubJetPhi,ParticleJetParticle->Phi()))));
		  DeltaR2_Particle=TMath::Sqrt((((ParticleJetParticle->Eta())-NextHardestParticleSubJetEta)*((ParticleJetParticle->Eta())- NextHardestParticleSubJetEta))+((RelativePhi(NextHardestParticleSubJetPhi,ParticleJetParticle->Phi()))*(RelativePhi(NextHardestParticleSubJetPhi,ParticleJetParticle->Phi()))));
		  
		  if(DeltaR1_Particle<=DeltaR2_Particle) {DeltaR_Particle=DeltaR1_Particle;}
		  else {DeltaR_Particle=DeltaR2_Particle;}
                  SubJetiness1_Numerator_Particle=SubJetiness1_Numerator_Particle+(ParticleJetParticle->Pt()*DeltaR1_Particle);
                  SubJetiness1_Denominator_Particle=SubJetiness1_Denominator_Particle+(ParticleJetParticle->Pt()*(fJetRadius));
		  SubJetiness2_Numerator_Particle=SubJetiness2_Numerator_Particle+(ParticleJetParticle->Pt()*DeltaR_Particle);
		  SubJetiness2_Denominator_Particle=SubJetiness2_Denominator_Particle+(ParticleJetParticle->Pt()*(fJetRadius));
		}
	      }
	    }
	    TempSubJetCounterDetector=0;
	    DetectorReclusterOk=0;
	    if(ReclustererDetector->AliEmcalJetFinder::Filter(Jet4, JetContDetector, dVtx)){  //reclustering jet1 using the jetfinderobject Reclusterer
	      DetectorReclusterOk=1;
	      SubJetCounterDetector=ReclustererDetector->GetNumberOfJets(); // Number of reclustered SubJets in each original jet                                       
	      for (Int_t i=0; i<SubJetCounterDetector; i++){ //Loops through each SubJet in a reclustered jet                                                           
		//if((Jet6=ReclustererDetector->GetJet(i)) && Jet6->GetNumberOfTracks()>1){; //jet2 is now set to the Subjet
		if((Jet6=ReclustererDetector->GetJet(i))){
		  if((Jet6->Pt())>HighestDetectorSubJetPt){ //This finds the highest pt subjet in each jet                                                                     
		    NextHighestDetectorSubJetPt=HighestDetectorSubJetPt;
		    HighestDetectorSubJetPt=Jet6->Pt();
		    NextHardestDetectorSubJetIndex=HardestDetectorSubJetIndex;
		    HardestDetectorSubJetIndex=i; //NSUBJETINESS              
		    NextHardestDetectorSubJetEta=HardestDetectorSubJetEta;                                                                                                
		    HardestDetectorSubJetEta=Jet6->Eta(); //NSUBJETINESS  
		    NextHardestDetectorSubJetPhi=HardestDetectorSubJetPhi;
		    HardestDetectorSubJetPhi=Jet6->Phi(); //NSUBJETINESS            
		  }
		  else if((Jet6->Pt())>NextHighestDetectorSubJetPt){ //This finds the 2nd highest pt subjet in each jet                                                    
		    NextHighestDetectorSubJetPt=Jet6->Pt();
		    NextHardestDetectorSubJetIndex=i; //NSUBJETINESS                                                                                                 
		    NextHardestDetectorSubJetEta=Jet6->Eta(); //NSUBJETINESS                                                                                   
		    NextHardestDetectorSubJetPhi=Jet6->Phi(); //NSUBJETINESS                                                                                                    
		  }
		}
		else{TempSubJetCounterDetector++;}
	      }
	      SubJetCounterDetector=SubJetCounterDetector-TempSubJetCounterDetector;
	      if (SubJetCounterDetector>=1){
		DetectorSubJetPtFrac=HighestDetectorSubJetPt/(Jet4->Pt());
		fhDetectorSubJetPtFrac->Fill(DetectorSubJetPtFrac); //Pt fraction of highest Pt subjet compared to original jet                                              
	      }
	      if (SubJetCounterDetector>=2){
                DetectorSubJetPtFrac2=(HighestDetectorSubJetPt+NextHighestDetectorSubJetPt)/(Jet4->Pt());
                fhDetectorSubJetPtFrac2->Fill(DetectorSubJetPtFrac2); //Pt fraction of two highest Pt subjets compared to original jet
	      }
	      if(SubJetCounterDetector>=2){
                SubJetiness1_Numerator_Detector=0;
                SubJetiness1_Denominator_Detector=0;	
		SubJetiness2_Numerator_Detector=0;
                SubJetiness2_Denominator_Detector=0;
                DeltaR_Detector=0;
                DeltaR1_Detector=0;
                DeltaR2_Detector=0;
                for (Int_t i=0; i< (Jet4->GetNumberOfTracks()); i++){  //loops through all tracks (particles in the jet                                                          
                  DetectorJetParticle = static_cast<AliVParticle*>(Jet4->TrackAt(i, JetContDetector->GetParticleContainer()->GetArray()));
                  DeltaR1_Detector=TMath::Sqrt((((DetectorJetParticle->Eta())-HardestDetectorSubJetEta)*((DetectorJetParticle->Eta())- HardestDetectorSubJetEta))+((RelativePhi(HardestDetectorSubJetPhi,DetectorJetParticle->Phi()))*(RelativePhi(HardestDetectorSubJetPhi,DetectorJetParticle->Phi()))));
		  
                  DeltaR2_Detector=TMath::Sqrt((((DetectorJetParticle->Eta())-NextHardestDetectorSubJetEta)*((DetectorJetParticle->Eta())- NextHardestDetectorSubJetEta))+((RelativePhi(NextHardestDetectorSubJetPhi,DetectorJetParticle->Phi()))*(RelativePhi(NextHardestDetectorSubJetPhi,DetectorJetParticle->Phi()))));
                  if(DeltaR1_Detector<=DeltaR2_Detector) {DeltaR_Detector=DeltaR1_Detector;}
                  else {DeltaR_Detector=DeltaR2_Detector;}
                  SubJetiness1_Numerator_Detector=SubJetiness1_Numerator_Detector+(DetectorJetParticle->Pt()*DeltaR1_Detector);
                  SubJetiness1_Denominator_Detector=SubJetiness1_Denominator_Detector+(DetectorJetParticle->Pt()*(fJetRadius));
                  SubJetiness2_Numerator_Detector=SubJetiness2_Numerator_Detector+(DetectorJetParticle->Pt()*DeltaR_Detector);
                  SubJetiness2_Denominator_Detector=SubJetiness2_Denominator_Detector+(DetectorJetParticle->Pt()*(fJetRadius));
		}
	      }
	    }
	    if((ParticleReclusterOk==1) && (DetectorReclusterOk==1)){
	      if((SubJetCounterParticle>=1) && (SubJetCounterDetector>=1)){
	      if (DetectorSubJetPtFrac>0){
		fhSubJetPtFracRatio->Fill(ParticleSubJetPtFrac/DetectorSubJetPtFrac);
	      }
	      fShapesVar[0] = (Jet3->Pt());
	      fShapesVar[1] = ParticleSubJetPtFrac;
	      fShapesVar[3] = (Jet4->Pt());
	      fShapesVar[4] = DetectorSubJetPtFrac;
	      }
	      if((SubJetCounterParticle>=2) && (SubJetCounterDetector>=2)){
		fShapesVar[2]=ParticleSubJetPtFrac2;
		fShapesVar[5]=DetectorSubJetPtFrac2;
		fhSubJetPtFrac2Ratio->Fill(ParticleSubJetPtFrac2/DetectorSubJetPtFrac2);
	      }
	      else{
                fShapesVar[2]=0;
                fShapesVar[5]=0;
	      }
	      if(SubJetCounterParticle>=2 && SubJetCounterDetector>=2 && SubJetiness1_Denominator_Particle!=0 && SubJetiness1_Denominator_Detector!=0 && SubJetiness2_Denominator_Particle!=0 && SubJetiness2_Denominator_Detector!=0){
		fhParticleSubJetiness1->Fill(SubJetiness1_Numerator_Particle/SubJetiness1_Denominator_Particle);
		fhDetectorSubJetiness1->Fill(SubJetiness1_Numerator_Detector/SubJetiness1_Denominator_Detector);
		fShapesVar[6]=SubJetiness1_Numerator_Particle/SubJetiness1_Denominator_Particle;
                  fShapesVar[7]=SubJetiness1_Numerator_Detector/SubJetiness1_Denominator_Detector;
                  fhSubJetiness1Ratio->Fill((SubJetiness1_Numerator_Particle/SubJetiness1_Denominator_Particle)/(SubJetiness1_Numerator_Detector/SubJetiness1_Denominator_Detector));		 
		  fhParticleSubJetiness2->Fill(SubJetiness2_Numerator_Particle/SubJetiness2_Denominator_Particle);
                  fhDetectorSubJetiness2->Fill(SubJetiness2_Numerator_Detector/SubJetiness2_Denominator_Detector);	      
		  fShapesVar[8]=SubJetiness2_Numerator_Particle/SubJetiness2_Denominator_Particle;
		  fShapesVar[9]=SubJetiness2_Numerator_Detector/SubJetiness2_Denominator_Detector;
		  fhSubJetiness2Ratio->Fill((SubJetiness2_Numerator_Particle/SubJetiness2_Denominator_Particle)/(SubJetiness2_Numerator_Detector/SubJetiness2_Denominator_Detector));
	      }
 
	      
	      else{
		fShapesVar[6]=0;
		fShapesVar[7]=0;
		fhParticleSubJetiness2->Fill(0);
		fShapesVar[8]=0;
		fShapesVar[9]=0;
	      }
	      fTreeResponseMatrixAxis->Fill();
	      
	      ParticleReclusterOk=0;
	      DetectorReclusterOk=0;  
	    }
	  }
	}
      }
    }
  }
  
  
  





  if (fJetShapeType == AliAnalysisTaskSubJetFraction::kData){

    AliEmcalJet *Jet1 = NULL; //Original Jet in the event                                                                                                                           
    AliEmcalJet *Jet2 = NULL; //Reclustered SubJet                                                                                                                               
    AliEmcalJetFinder *Reclusterer = new AliEmcalJetFinder("SubjetFinder3"); //JetFinder Object for reclustered jets                                 
    Reclusterer->SetRadius(fSubJetRadius);
    Reclusterer->SetJetMinPt(fSubJetMinPt);
    Reclusterer->SetJetAlgorithm(fSubJetAlgorithm); //0 for anti-kt     1 for kt                                                                                     
    Double_t Jet1_Pt_Subtracted = 0;
    AliJetContainer *JetCont= GetJetContainer(0); //Jet Container for event                                                                                                       
    Double_t JetPhi=0;
    Double_t JetEta=0;
    Int_t JetCounter=0; //Counts number of jets in event                                                                                                                           
    Int_t SubJetCounter=0; //Counts number of reclustered jets in a jet                                                                                                            
    Double_t HighestSubJetPt;
    Double_t NextHighestSubJetPt;
    Double_t HighestSubJetEnergy;
    Double_t NextHighestSubJetEnergy;
    Int_t NumberOfJetTracks;
    AliVParticle *JetParticle = 0x0; //Individual constituent tracks (particles) in a jet                                                                                           
    Double_t DeltaPhi=0;
    Double_t JetParticlePhi=0;
    Double_t JetParticlePt=0;
    Double_t Angularity_Numerator=0;
    Double_t Angularity_Denominator=0;
    Double_t PTD_Numerator=0;
    Double_t PTD_Denominator=0;

    Int_t HardestSubJetIndex=-1;
    Int_t NextHardestSubJetIndex=-1;
    Double_t HardestSubJetEta=-10;
    Double_t NextHardestSubJetEta=-10;
    Double_t HardestSubJetPhi=-10;
    Double_t NextHardestSubJetPhi=-10;
    Double_t SubJetiness1_Numerator=0;
    Double_t SubJetiness1_Denominator=0;
    Double_t SubJetiness2_Numerator=0;
    Double_t SubJetiness2_Denominator=0;
    Double_t DeltaR=0;
    Double_t DeltaR1=0;
    Double_t DeltaR2=0;
    Double_t TempSubJetCounter=-1;




    if(JetCont) {
      //    cout << "Jet Container Size" << JetCont->GetNEntries()<<endl;
      fhEventCounter->Fill(1); //Number of events with a jet container
      JetCont->ResetCurrentID();
      // while((Jet1=JetCont->GetNextAcceptJet()) && ((Jet1->GetNumberOfTracks())>1) && (Jet1->Pt()>=10)) {      //exludes one or zero(?) track jets
      while((Jet1=JetCont->GetNextAcceptJet()) && (Jet1->Pt()>=10)) {
	if(!Jet1) {
	  fhEventCounter->Fill(3);
	  continue;
	}
	else {
	  if((Jet1->GetNumberOfTracks())==0){
	    fhEventCounter->Fill(10);
	  }
	  if((Jet1->GetNumberOfTracks())==1){
	    fhEventCounter->Fill(11);
	  }
	  JetCounter++;
	  fhJetPt->Fill(Jet1->Pt());
	  JetPhi=Jet1->Phi();
	  if(JetPhi < -1*TMath::Pi()) JetPhi += (2*TMath::Pi());
	  else if (JetPhi > TMath::Pi()) JetPhi -= (2*TMath::Pi());
	  fhJetPhi->Fill(JetPhi);
	  JetEta=Jet1->Eta();
	  fhJetEta->Fill(JetEta);
	  fhJetMass->Fill(Jet1->M());
	  fhEventCounter->Fill(4); //Number of Jets found in all events
	  HighestSubJetPt=0.0;
	  NextHighestSubJetPt=0.0;
	  HighestSubJetEnergy=0.0;
	  NextHighestSubJetEnergy=0.0;
	  fhJetRadius->Fill(TMath::Sqrt((Jet1->Area()/TMath::Pi()))); //Radius of Jets per event
	  
	  ///////////////////////////////////////////////Jet Angularity and PTD/////////////////////////////////////////////////////////////////////////////////////
	  JetPhi=Jet1->Phi();//Done again becasue previously we changed JetPhi to put it between certain limits
	  NumberOfJetTracks = Jet1->GetNumberOfTracks();
	  fhNumberOfJetTracks->Fill(NumberOfJetTracks);
	  Angularity_Numerator=0;  //Reset these values
	  Angularity_Denominator=0;
	  PTD_Numerator=0;
	  PTD_Denominator=0;
	  for (Int_t i=0; i< NumberOfJetTracks; i++){  //loops through all tracks (particles in the jet
	    JetParticle = static_cast<AliVParticle*>(Jet1->TrackAt(i, JetCont->GetParticleContainer()->GetArray()));  //TrackAt(i,TClonesArray) returns the ith particle in the TClones Array. This sets JetParticle to each particle in the Jet
	    if(!JetParticle) continue;
	    JetParticlePt=JetParticle->Pt();
	    JetParticlePhi=JetParticle->Phi();
	    /*
	    if(JetPhi < -1*TMath::Pi()) JetPhi += (2*TMath::Pi()); // Turns the range of 0to2Pi into -PitoPi ???????????
	    else if (JetPhi > TMath::Pi()) JetPhi -= (2*TMath::Pi());
	    if(JetParticlePhi < -1*TMath::Pi()) JetParticlePhi += (2*TMath::Pi()); 
	    else if (JetParticlePhi > TMath::Pi()) JetParticlePhi -= (2*TMath::Pi());
	    DeltaPhi=JetParticlePhi-JetPhi;
	    if(DeltaPhi < -1*TMath::Pi()) DeltaPhi += (2*TMath::Pi());
	    else if (DeltaPhi > TMath::Pi()) DeltaPhi -= (2*TMath::Pi());
	    */
	    DeltaPhi=RelativePhi(JetPhi,JetParticlePhi);
	    Angularity_Numerator=Angularity_Numerator+(JetParticlePt*TMath::Sqrt(((JetParticle->Eta()-JetEta)*(JetParticle->Eta()-JetEta))+(DeltaPhi*DeltaPhi)));
	    Angularity_Denominator= Angularity_Denominator+JetParticlePt;
	    PTD_Numerator=PTD_Numerator+(JetParticlePt*JetParticlePt);
	    PTD_Denominator=PTD_Denominator+JetParticlePt;
	    
	  }
	  if(Angularity_Denominator!=0)	fhJetAngularity->Fill(Angularity_Numerator/Angularity_Denominator);
	  if(PTD_Denominator!=0) fhJetPTD->Fill((TMath::Sqrt(PTD_Numerator))/PTD_Denominator);
	  //////////////////////////////////////////////////Jet Reclustering and NSubJetiness//////////////////////////////////////////////////////////////////////////////
	  HardestSubJetIndex=-1;
	  NextHardestSubJetIndex=-1; 
	  HardestSubJetEta=-10;    
	  NextHardestSubJetEta=-10;    
          HardestSubJetPhi=-10;                                                                                                                       
	  NextHardestSubJetPhi=-10;
          SubJetiness1_Numerator=0;
          SubJetiness1_Denominator=0;
	  SubJetiness2_Numerator=0;
	  SubJetiness2_Denominator=0;
	  DeltaR=0;
	  DeltaR1=0;
	  DeltaR2=0;
	  TempSubJetCounter=0;
 
	  if(Reclusterer->AliEmcalJetFinder::Filter(Jet1, JetCont, dVtx)){  //reclustering jet1 using the jetfinderobject Reclusterer
	    fhEventCounter->Fill(5); //Number of times jets were reclustered
	    SubJetCounter=Reclusterer->GetNumberOfJets(); // Number of reclustered SubJets in each original jet
	    for (Int_t i=0; i<SubJetCounter; i++){ //Loops through each SubJet in a reclustered jet
	      //if((Jet2=Reclusterer->GetJet(i)) && ((Jet2->GetNumberOfTracks())>1)){ //jet2 is now set to the Subjet
	      if((Jet2=Reclusterer->GetJet(i))){
		fhEventCounter->Fill(6); //Number of overall subjets in all events 		
		fhSubJetPt->Fill(Jet2->Pt()); 
		fhSubJetMass->Fill(Jet2->M());
		fhNumberOfSubJetTracks->Fill(Jet2->GetNumberOfTracks());
	      // cout << Jet2->Area()<<endl;
		fhSubJetRadius->Fill(TMath::Sqrt((Jet2->Area()/TMath::Pi()))); //Radius of SubJets per event  
		if((Jet2->Pt())>HighestSubJetPt){ //This finds the highest pt subjet in each jet
		  NextHighestSubJetPt=HighestSubJetPt;	      
		  HighestSubJetPt=Jet2->Pt();
		  NextHardestSubJetIndex=HardestSubJetIndex;
		  HardestSubJetIndex=i; //NSUBJETINESS
		  NextHardestSubJetEta=HardestSubJetEta;
		  HardestSubJetEta=Jet2->Eta(); //NSUBJETINESS
		  NextHardestSubJetPhi=HardestSubJetPhi;
		  HardestSubJetPhi=Jet2->Phi(); //NSUBJETINESS
		}
		else if((Jet2->Pt())>NextHighestSubJetPt){ //This finds the 2nd highest pt subjet in each jet                                                                    
		  NextHighestSubJetPt=Jet2->Pt();
		  NextHardestSubJetIndex=i; //NSUBJETINESS
		  NextHardestSubJetEta=Jet2->Eta(); //NSUBJETINESS                                                                                                      
		  NextHardestSubJetPhi=Jet2->Phi(); //NSUBJETINESS 
		}
		if((Jet2->E())>HighestSubJetEnergy){ //This finds the highest Energy subjet in each jet                                                                         
		  NextHighestSubJetEnergy=HighestSubJetEnergy;
		  HighestSubJetEnergy=Jet2->E();
		}
		else if((Jet2->E())>NextHighestSubJetEnergy){ //This finds the 2nd highest Energy subjet in each jet                                                               
		  NextHighestSubJetEnergy=Jet2->E();
		}
	      }
	      else {TempSubJetCounter++;}
	    }
	    SubJetCounter=SubJetCounter-TempSubJetCounter;
	    fhSubJetCounter->Fill(SubJetCounter); //Counts Number of SubJets in each jet     
	    
	    if (SubJetCounter>=1){
	      fhSubJetPtFrac->Fill((HighestSubJetPt)/(Jet1->Pt())); //Pt fraction of highest Pt subjet compared to original jet
	      fhSubJetPtLoss->Fill(((Jet1->Pt())-HighestSubJetPt)/(Jet1->Pt()));  //Pt difference of jet and its highest Pt subjet
	      fhSubJetEnergyFrac->Fill((HighestSubJetEnergy)/(Jet1->E())); //Energy fraction of most energetic subjet compared to original jet       
	      fhSubJetEnergyLoss->Fill(((Jet1->E())-HighestSubJetEnergy)/(Jet1->E()));	//Energy differance of jet and its most energetic subjet  
	      fhEventCounter->Fill(8);
	    }
	    if (SubJetCounter>=2){
	      fhSubJetPtFrac2->Fill((HighestSubJetPt+NextHighestSubJetPt)/(Jet1->Pt())); //Pt fraction of two highest Pt subjets compared to original jet  
	      fhSubJetPtLoss2->Fill(((Jet1->Pt())-(HighestSubJetPt+NextHighestSubJetPt))/(Jet1->Pt())); //Pt difference of jet and its two highest Pt subjets
	      fhSubJetEnergyFrac2->Fill((HighestSubJetEnergy+NextHighestSubJetEnergy)/(Jet1->E())); //Energy fraction of the two most energetic subjets compared to original jet 
	      fhSubJetEnergyLoss2->Fill(((Jet1->E())-(HighestSubJetEnergy+NextHighestSubJetEnergy))/(Jet1->Pt())); //Energy difference of jet and its two most energetic subjets  
	      fhEventCounter->Fill(9);
	    }
	    if (SubJetCounter>=2){
              SubJetiness1_Numerator=0;
              SubJetiness1_Denominator=0;	      
	      SubJetiness2_Numerator=0;
	      SubJetiness2_Denominator=0;
	      DeltaR=0;
	      DeltaR1=0;
	      DeltaR2=0;
	      for (Int_t i=0; i< NumberOfJetTracks; i++){  //loops through all tracks (particles in the jet                                                               
		DeltaR=0;
		DeltaR1=0;
		DeltaR2=0;
		JetParticle = static_cast<AliVParticle*>(Jet1->TrackAt(i, JetCont->GetParticleContainer()->GetArray()));
		DeltaR1=TMath::Sqrt((((JetParticle->Eta())-HardestSubJetEta)*((JetParticle->Eta())- HardestSubJetEta))+((RelativePhi(HardestSubJetPhi,JetParticle->Phi()))*(RelativePhi(HardestSubJetPhi,JetParticle->Phi()))));
		DeltaR2=TMath::Sqrt((((JetParticle->Eta())-NextHardestSubJetEta)*((JetParticle->Eta())- NextHardestSubJetEta))+((RelativePhi(NextHardestSubJetPhi,JetParticle->Phi()))*(RelativePhi(NextHardestSubJetPhi,JetParticle->Phi()))));
		if(DeltaR1<=DeltaR2) {DeltaR=DeltaR1;}
		else {DeltaR=DeltaR2;}
		SubJetiness1_Numerator=SubJetiness1_Numerator+(JetParticle->Pt()*DeltaR1);
                SubJetiness1_Denominator=SubJetiness1_Denominator+(JetParticle->Pt()*fJetRadius);
		SubJetiness2_Numerator=SubJetiness2_Numerator+(JetParticle->Pt()*DeltaR);
		SubJetiness2_Denominator=SubJetiness2_Denominator+(JetParticle->Pt()*fJetRadius);
	      }
	    }
            if(SubJetCounter>=1){
              fShapesVar[0]=Jet1->Pt();
              fShapesVar[1]=(HighestSubJetPt)/(Jet1->Pt());
	    }
	    else{             
	      fShapesVar[0]=0;
              fShapesVar[1]=0;
	    }
	    if(SubJetCounter>=2){
	      fhEventCounter->Fill(12);
	      fShapesVar[2]=(HighestSubJetPt+NextHighestSubJetPt)/(Jet1->Pt());
	    }
	    else{fShapesVar[2]=0;}
	    if (SubJetCounter>=2 && SubJetiness1_Denominator!=0 && SubJetiness2_Denominator!=0){
	      fhSubJetiness1->Fill(SubJetiness1_Numerator/SubJetiness1_Denominator);
	      fShapesVar[6]=SubJetiness1_Numerator/SubJetiness1_Denominator;
	      fhSubJetiness2->Fill(SubJetiness2_Numerator/SubJetiness2_Denominator);
	      fShapesVar[8]=SubJetiness2_Numerator/SubJetiness2_Denominator;
		  fh2to1SubJetinessRatio->Fill((SubJetiness2_Numerator/SubJetiness2_Denominator)/(SubJetiness1_Numerator/SubJetiness1_Denominator));
	    }
	    else{
	      fShapesVar[6]=0;
	      fShapesVar[8]=0;
	    }
	    fShapesVar[3]=0;
	    fShapesVar[4]=0;
	    fShapesVar[5]=0;
	    fShapesVar[7]=0;
	    fShapesVar[9]=0;
	    fTreeResponseMatrixAxis->Fill();   
	  }
	}
      }
      fhJetCounter->Fill(JetCounter); //Number of Jets in Each Event
      
    }
    else {fhEventCounter->Fill(2);} //Events with no jet container
  }
  // if (JetCounter>40) cout << "Too Big!" << JetCounter <<endl;
  
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
