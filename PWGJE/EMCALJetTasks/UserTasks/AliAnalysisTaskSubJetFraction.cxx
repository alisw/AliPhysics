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
Double_t XBinsJetPt[28]={0,0.50,1.00,2.00,3.00,4.00,5.00,6.00,7.00,8.00,9.00,10.00,12.00,14.00,16.00,18.00,20.00,25.00,30.00,35.00,40.00,45.00,50.00,60.00,70.00,80.00,90.00,100.00}; //for jet pt
Int_t XBinsJetPtSize=28;
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
  fhJetPt(0x0),
  fhJetPhi(0x0),
  fhJetEta(0x0),
  fhJetMass(0x0),
  fhJetRadius(0x0),
  fhJetAngularity(0x0),
  fhJetPTD(0x0),
  fhNumberOfJetTracks(0x0),
  fhSubJetPt(0x0),
  fhSubJetRadius(0x0),
  fhSubJetPtFrac(0x0),
  fhSubJetPtFrac2(0x0),
  fhSubJetPtLoss(0x0),
  fhSubJetPtLoss2(0x0),
  fhSubJetEnergyFrac(0x0),
  fhSubJetEnergyFrac2(0x0),
  fhSubJetEnergyLoss(0x0),
  fhSubJetEnergyLoss2(0x0),
  fhJetiness(0x0),
  fhEventCounter(0x0),
  fhJetCounter(0x0),
  fhSubJetCounter(0x0),
  fhPtRatio(0x0),
  fhParticleSubJetPtFrac(0x0),
  fhDetectorSubJetPtFrac(0x0),
  fhSubJetPtFracRatio(0x0),
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
  fhJetPt(0x0),
  fhJetPhi(0x0),
  fhJetEta(0x0),
  fhJetMass(0x0),
  fhJetRadius(0x0),
  fhJetAngularity(0x0),
  fhJetPTD(0x0),
  fhNumberOfJetTracks(0x0),
  fhSubJetPt(0x0),
  fhSubJetRadius(0x0),
  fhSubJetPtFrac(0x0),
  fhSubJetPtFrac2(0x0),
  fhSubJetPtLoss(0x0),
  fhSubJetPtLoss2(0x0),
  fhSubJetEnergyFrac(0x0),
  fhSubJetEnergyFrac2(0x0),
  fhSubJetEnergyLoss(0x0),
  fhSubJetEnergyLoss2(0x0),
  fhJetiness(0x0),
  fhEventCounter(0x0),
  fhJetCounter(0x0),
  fhSubJetCounter(0x0),
  fhPtRatio(0x0),
  fhParticleSubJetPtFrac(0x0),
  fhDetectorSubJetPtFrac(0x0),
  fhSubJetPtFracRatio(0x0),
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

  const Int_t nVar = 4;
  fShapesVar = new Double_t [nVar]; //shapes used for tagging   
  if(fJetShapeType==AliAnalysisTaskSubJetFraction::kTrueDet){
  //create a tree used for the MC data and making a 4D response matrix
  
    fTreeResponseMatrixAxis = new TTree("fTreeJetShape", "fTreeJetShape");
    TString *fShapesVarNames = new TString [nVar];
    
    fShapesVarNames[0] = "Pt_Particle_Level";
    fShapesVarNames[1] = "Frac_Particle_Level";
    fShapesVarNames[2] = "Pt_Detector_Level";
    fShapesVarNames[3] = "Frac_Detector_Level";
    
    for(Int_t ivar=0; ivar < nVar; ivar++){
      cout<<"looping over variables"<<endl;
      fTreeResponseMatrixAxis->Branch(fShapesVarNames[ivar].Data(), &fShapesVar[ivar], Form("%s/F", fShapesVarNames[ivar].Data()));
    }
    
  }
  

  if (fJetShapeType==AliAnalysisTaskSubJetFraction::kData){
    
    fhJetPt= new TH1F("fhJetPt", "Jet Pt", (XBinsJetPtSize)-1, XBinsJetPt);
    fOutput->Add(fhJetPt);
    fhJetPhi= new TH1F("fhJetPhi", "Jet Phi", Phi_Bins, Phi_Low, Phi_Up);
    fOutput->Add(fhJetPhi);
    fhJetEta= new TH1F("fhJetEta", "Jet Eta", Eta_Bins, Eta_Low, Eta_Up);
    fOutput->Add(fhJetEta);
    fhJetMass= new TH1F("fhJetMass", "Jet Mass", (XBinsJetMassSize)-1, XBinsJetMass);
    fOutput->Add(fhJetMass);
    fhJetRadius= new TH1F("fhJetRadius", "Jet Radius", 100, -0.05,0.995);
    fOutput->Add(fhJetRadius);
    fhJetAngularity= new TH1F("fhJetAngularity", "Jet Angularity", 2000, -10,10);
    fOutput->Add(fhJetAngularity);
    fhJetPTD= new TH1F("fhJetPTD", "Jet PTD", 2000, -10,10);
    fOutput->Add(fhJetPTD);
    fhNumberOfJetTracks= new TH1F("fhNumberOfJetTracks", "Number of Tracks within a Jet", 30, -0.5,29.5);
    fOutput->Add(fhNumberOfJetTracks);
    fhSubJetRadius= new TH1F("fhSubJetRadius", "SubJet Radius", 100, -0.05,0.995);
    fOutput->Add(fhSubJetRadius);
    fhSubJetPt= new TH1F("fhSubJetPt", "SubJet Pt", (XBinsJetPtSize)-1, XBinsJetPt);
    fOutput->Add(fhSubJetPt);
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
    fhJetiness= new TH1F("fhJetiness", "Jettiness",10100, -0.05,500.05);
    fOutput->Add(fhJetiness);
    fhJetCounter= new TH1F("fhJetCounter", "Jet Counter", 100, -0.5, 99.5);
    fOutput->Add(fhJetCounter);
    fhSubJetCounter = new TH1F("fhSubJetCounter", "SubJet Counter",50, -0.5,49.5);
    fOutput->Add(fhSubJetCounter);
  }
  if(fJetShapeType==AliAnalysisTaskSubJetFraction::kTrueDet){
    fhPtRatio= new TH1F("fhPtRatio", "MC pt ratio for particle and detector level jets",100, -0.5,9.5);
    fOutput->Add(fhPtRatio);
    fhParticleSubJetPtFrac= new TH1F("fhParticleSubJetPtFrac", "Pt Fraction of Highest Pt Subjet compared to original Jet for MC Particle Level data",101, -0.05,1.05);
    fOutput->Add(fhParticleSubJetPtFrac);
    fhDetectorSubJetPtFrac= new TH1F("fhDetectorSubJetPtFrac", "Pt Fraction of Highest Pt Subjet compared to original Jet for MC Detector Level data",101, -0.05,1.05);
    fOutput->Add(fhDetectorSubJetPtFrac);
    fhSubJetPtFracRatio= new TH1F("fhSubJetPtFracRatio", "Ratio of Pt Fraction of Highest Pt Subjet compared to original Jet for MC particle and Detector Level data", 1010,-0.05,10.5); 
    fOutput->Add(fhSubJetPtFracRatio);
    fOutput->Add(fTreeResponseMatrixAxis);
  }
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
    AliEmcalJetFinder *ReclustererParticle =NULL; //JetFinder Object for reclustered jets from MC particle Level data                                                          
    ReclustererParticle->SetRadius(0.1);
    ReclustererParticle->SetJetMinPt(0.5);
    ReclustererParticle->SetJetAlgorithm(0); //0 for anti-kt     1 for kt  
    AliEmcalJetFinder *ReclustererDetector =NULL; //JetFinder Object for reclustered jets from MC detector Level data                                                         
    ReclustererDetector->SetRadius(0.1);
    ReclustererDetector->SetJetMinPt(0.5);
    ReclustererDetector->SetJetAlgorithm(0); //0 for anti-kt     1 for kt           
    Int_t SubJetCounterParticle=0;
    Int_t SubJetCounterDetector=0;
    Int_t ParticleReclusterOk=0; //makes sure clustering has happened before commiting values to the response matrix
    Int_t DetectorReclusterOk=0;
    Double_t HighestParticleSubJetPt=-1;
    Double_t HighestDetectorSubJetPt=-1;
    Double_t ParticleSubJetPtFrac=0;
    Double_t DetectorSubJetPtFrac=0;
    if(JetContParticle){
      JetContParticle->ResetCurrentID();
      while((Jet3=JetContParticle->GetNextAcceptJet()) && ((Jet3->GetNumberOfTracks())>1) && (Jet3->Pt()>=10)) {      //exludes one or zero(?) track jets                           
	if(!Jet3) {
	  continue;
	}
	else {
	  if((Jet4 = (Jet3->ClosestJet()))){
	    fhPtRatio->Fill(Jet3->Pt()/Jet4->Pt());
	    HighestParticleSubJetPt=-1;
	    HighestDetectorSubJetPt=-1;
	    if(ReclustererParticle->AliEmcalJetFinder::Filter(Jet3, JetContParticle, dVtx)){  //reclustering jet1 using the jetfinderobject Reclusterer                
	      ParticleReclusterOk=1;
	      SubJetCounterParticle=ReclustererParticle->GetNumberOfJets(); // Number of reclustered SubJets in each original jet                                            
	      for (Int_t i=0; i<SubJetCounterParticle; i++){ //Loops through each SubJet in a reclustered jet                                                              
		Jet5=ReclustererParticle->GetJet(i); //jet2 is now set to the Subjet                                                                                           
		if((Jet5->Pt())>HighestParticleSubJetPt){ //This finds the highest pt subjet in each jet                                                                       
		  HighestParticleSubJetPt=Jet5->Pt();
		}
	      }
	      if (SubJetCounterParticle>=1){
	      ParticleSubJetPtFrac=HighestParticleSubJetPt/(Jet3->Pt());
	      fhParticleSubJetPtFrac->Fill(ParticleSubJetPtFrac); //Pt fraction of highest Pt subjet compared to original jet                            
	      }
	    }
	    
	    if(ReclustererDetector->AliEmcalJetFinder::Filter(Jet4, JetContDetector, dVtx)){  //reclustering jet1 using the jetfinderobject Reclusterer
	      DetectorReclusterOk=1;
	      SubJetCounterDetector=ReclustererDetector->GetNumberOfJets(); // Number of reclustered SubJets in each original jet                                       
	      for (Int_t i=0; i<SubJetCounterDetector; i++){ //Loops through each SubJet in a reclustered jet                                                            
		Jet6=ReclustererDetector->GetJet(i); //jet2 is now set to the Subjet                                                                                             
		if((Jet6->Pt())>HighestDetectorSubJetPt){ //This finds the highest pt subjet in each jet                                                                     
		  HighestDetectorSubJetPt=Jet6->Pt();
		}
	      }
	      if (SubJetCounterDetector>=1){
		DetectorSubJetPtFrac=HighestDetectorSubJetPt/(Jet4->Pt());
		fhDetectorSubJetPtFrac->Fill(DetectorSubJetPtFrac); //Pt fraction of highest Pt subjet compared to original jet                                              
	      }
	    }
	    if((ParticleReclusterOk==1) && (DetectorReclusterOk==1)){
	      if (DetectorSubJetPtFrac>0) fhSubJetPtFracRatio->Fill(ParticleSubJetPtFrac/DetectorSubJetPtFrac);
	      fShapesVar[0] = (Jet3->Pt());
	      fShapesVar[1] = ParticleSubJetPtFrac;
	      fShapesVar[2] = (Jet4->Pt());
	      fShapesVar[3] = DetectorSubJetPtFrac;
	      fTreeResponseMatrixAxis->Fill();
	    }
	    ParticleReclusterOk=0;
	    DetectorReclusterOk=0;
	    
	  }
	}
      }
    }
    
  }
  






  if (fJetShapeType == AliAnalysisTaskSubJetFraction::kData){

    AliEmcalJet *Jet1 = NULL; //Original Jet in the event                                                                                                                           
    AliEmcalJet *Jet2 = NULL; //Reclustered SubJet                                                                                                                               
    AliEmcalJetFinder *Reclusterer =NULL; //JetFinder Object for reclustered jets                                                                                                  
    Reclusterer->SetRadius(0.1);
    Reclusterer->SetJetMinPt(0.5);
    Reclusterer->SetJetAlgorithm(0); //0 for anti-kt     1 for kt                                                                                                                  
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


    if(JetCont) {
      //    cout << "Jet Container Size" << JetCont->GetNEntries()<<endl;
      fhEventCounter->Fill(1); //Number of events with a jet container
      JetCont->ResetCurrentID();
      while((Jet1=JetCont->GetNextAcceptJet()) && ((Jet1->GetNumberOfTracks())>1) && (Jet1->Pt()>=10)) {      //exludes one or zero(?) track jets
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
	  if (JetPhi>((3*TMath::Pi())/2.0)) JetPhi=JetPhi-(2*TMath::Pi());
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
	    if(JetPhi < -1*TMath::Pi()) JetPhi += (2*TMath::Pi()); // Turns the range of 0to2Pi into -PitoPi ???????????
	    else if (JetPhi > TMath::Pi()) JetPhi -= (2*TMath::Pi());
	    if(JetParticlePhi < -1*TMath::Pi()) JetParticlePhi += (2*TMath::Pi()); 
	    else if (JetParticlePhi > TMath::Pi()) JetParticlePhi -= (2*TMath::Pi());
	    DeltaPhi=JetParticlePhi-JetPhi;
	    if(DeltaPhi < -1*TMath::Pi()) DeltaPhi += (2*TMath::Pi());
	    else if (DeltaPhi > TMath::Pi()) DeltaPhi -= (2*TMath::Pi());
	    Angularity_Numerator=Angularity_Numerator+(JetParticlePt*TMath::Sqrt(((JetParticle->Eta()-JetEta)*(JetParticle->Eta()-JetEta))+(DeltaPhi*DeltaPhi)));
	    Angularity_Denominator= Angularity_Denominator+JetParticlePt;
	    PTD_Numerator=PTD_Numerator+(JetParticlePt*JetParticlePt);
	    PTD_Denominator=PTD_Denominator+JetParticlePt;
	    
	  }
	  if(Angularity_Denominator!=0)	fhJetAngularity->Fill(Angularity_Numerator/Angularity_Denominator);
	  if(PTD_Denominator!=0) fhJetPTD->Fill((TMath::Sqrt(PTD_Numerator))/PTD_Denominator);
	  //////////////////////////////////////////////////Jet Reclustering//////////////////////////////////////////////////////////////////////////////
	  if(Reclusterer->AliEmcalJetFinder::Filter(Jet1, JetCont, dVtx)){  //reclustering jet1 using the jetfinderobject Reclusterer
	    fhEventCounter->Fill(5); //Number of times jets were reclustered
	    SubJetCounter=Reclusterer->GetNumberOfJets(); // Number of reclustered SubJets in each original jet
	    fhSubJetCounter->Fill(SubJetCounter); //Counts Number of SubJets in each jet
	    for (Int_t i=0; i<SubJetCounter; i++){ //Loops through each SubJet in a reclustered jet
	      fhEventCounter->Fill(6); //Number of overall subjects in all events
	      Jet2=Reclusterer->GetJet(i); //jet2 is now set to the Subjet
	      fhSubJetPt->Fill(Jet2->Pt()); 
	      // cout << Jet2->Area()<<endl;
	      fhSubJetRadius->Fill(TMath::Sqrt((Jet2->Area()/TMath::Pi()))); //Radius of SubJets per event  
	      if((Jet2->Pt())>HighestSubJetPt){ //This finds the highest pt subjet in each jet
		NextHighestSubJetPt=HighestSubJetPt;	      
		HighestSubJetPt=Jet2->Pt();
	      }
	      else if((Jet2->Pt())>NextHighestSubJetPt){ //This finds the 2nd highest pt subjet in each jet                                                                    
		NextHighestSubJetPt=Jet2->Pt();
	      }
	      if((Jet2->E())>HighestSubJetEnergy){ //This finds the highest Energy subjet in each jet                                                                         
		NextHighestSubJetEnergy=HighestSubJetEnergy;
		HighestSubJetEnergy=Jet2->E();
	      }
	      else if((Jet2->E())>NextHighestSubJetEnergy){ //This finds the 2nd highest Energy subjet in each jet                                                               
		NextHighestSubJetEnergy=Jet2->E();
	      }     
	    }
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
	  /*
	                                               ////////////////Jetiness///////////
	  Double_t DeltaR_temp=0;
	  Double_t DeltaR=100;
	  Double_t Jetiness_Numerator=0;
	  Double_t Jetiness_Denominator=0;
	  for (Int_t i=0; i< NumberOfJetTracks; i++){  //loops through all tracks (particles in the jet                                                                          
	    JetParticle = static_cast<AliVParticle*>(Jet1->TrackAt(i, JetCont->GetParticleContainer()->GetArray()));
	    for (Int_t i=0; i<SubJetCounter; i++){ //Loops through each SubJet in a reclustered jet                                             
	      Jet2=Reclusterer->GetJet(i); //jet2 is now set to the Subjet     
	      DeltaR_temp=TMath::Sqrt(((JetParticle->Eta()-(Jet2->Eta()))*(JetParticle->Eta()-(Jet2->Eta())))+((JetParticle->Phi()-Jet2->Phi())*(JetParticle->Phi()-Jet2->Phi())));
	      if (DeltaR_temp<DeltaR){
		DeltaR=DeltaR_temp;
	      }
	    }
	    Jetiness_Numerator=Jetiness_Numerator+(JetParticle->Pt()*DeltaR);
	    Jetiness_Denominator=Jetiness_Denominator+(JetParticle->Pt()*TMath::Sqrt((Jet1->Area()/TMath::Pi())));				 
	    DeltaR_temp=0;
	    DeltaR=100;	  
	  }
	  // cout << Jetiness_Numerator/Jetiness_Denominator<<endl;
	  fhJetiness->Fill(Jetiness_Numerator/Jetiness_Denominator);
	 */ 
	  }  
	}
      }
      fhJetCounter->Fill(JetCounter); //Number of Jets in Each Event
      
    }
    else fhEventCounter->Fill(2); //Events with no jet container
  } 
  // if (JetCounter>40) cout << "Too Big!" << JetCounter <<endl;
  
  return kTRUE;
  
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
  }

}

