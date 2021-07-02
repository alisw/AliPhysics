//
// Basic analysis task.
//
// Basic analysis task template for analysis jets storing information in both tree
// branches and histograms


#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TTree.h>
#include <TList.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include <TProfile.h>
#include <TChain.h>

// aliroot Headers
#include "AliAODEvent.h"
#include "AliParticleContainer.h"
#include "AliJetContainer.h"
#include "AliMCEvent.h"
#include "AliEmcalJet.h"

//My Header
#include "AliAnalysisTaskJetChargeFlavourTemplates.h"

//Globals
using std::cout;
using std::endl;



ClassImp(AliAnalysisTaskJetChargeFlavourTemplates)

//________________________________________________________________________
AliAnalysisTaskJetChargeFlavourTemplates::AliAnalysisTaskJetChargeFlavourTemplates() :
  AliAnalysisTaskEmcalJet("AliAnalysisTaskJetChargeFlavourTemplates", kTRUE),
  fContainer(0),
  pChain(0),
  fPtThreshold(-9999.),
  fCentSelectOn(kTRUE),
  fCentMin(0),
  fCentMax(10),
  fJetRadius(0),
  JetChargeK(0.5),
  MotherFraction(0.8),
  JetMidPt(40),
  JetHighPt(80),

  fhJetPt(0x0),
  fhJetPhi(0x0),
  fhJetEta(0x0),

  JC(0x0),

  JCUp(0x0),
  JCDown(0x0),
  JCGluon(0x0),
  JCOther(0x0),
  JCUnmatched(0x0),

  JCLow(0x0),

  JCUpLow(0x0),
  JCDownLow(0x0),
  JCGluonLow(0x0),
  JCOtherLow(0x0),
  JCUnmatchedLow(0x0),


  JCMid(0x0),

  JCUpMid(0x0),
  JCDownMid(0x0),
  JCGluonMid(0x0),
  JCOtherMid(0x0),
  JCUnmatchedMid(0x0),

  JCHigh(0x0),

  JCUpHigh(0x0),
  JCDownHigh(0x0),
  JCGluonHigh(0x0),
  JCOtherHigh(0x0),
  JCUnmatchedHigh(0x0),

  fhParticleJetPt(0x0),
  fhParticleJetPhi(0x0),
  fhParticleJetEta(0x0),

  JCParticle(0x0),

  JCParticleUp(0x0),
  JCParticleDown(0x0),
  JCParticleGluon(0x0),
  JCParticleOther(0x0),
  JCParticleUnmatched(0x0),

  JCParticleLow(0x0),

  JCParticleUpLow(0x0),
  JCParticleDownLow(0x0),
  JCParticleGluonLow(0x0),
  JCParticleOtherLow(0x0),
  JCParticleUnmatchedLow(0x0),


  JCParticleMid(0x0),

  JCParticleUpMid(0x0),
  JCParticleDownMid(0x0),
  JCParticleGluonMid(0x0),
  JCParticleOtherMid(0x0),
  JCParticleUnmatchedMid(0x0),

  JCParticleHigh(0x0),

  JCParticleUpHigh(0x0),
  JCParticleDownHigh(0x0),
  JCParticleGluonHigh(0x0),
  JCParticleOtherHigh(0x0),
  JCParticleUnmatchedHigh(0x0),

  PtComparison(0x0),
  JCComparison(0x0),
  JCComparisonUp(0x0),
  JCComparisonDown(0x0),
  JCComparisonGluon(0x0),
  JCComparisonOther(0x0),
  JCComparisonUnmatched(0x0),

  Pt2DCompare(0x0),
  JC2DCompare(0x0),

  PtComparisonVsJCDiff(0x0),
  JCComparisonVsPtDiff(0x0),
  JCComparisonVsPtDiffUp(0x0),
  JCComparisonVsPtDiffDown(0x0),
  JCComparisonVsPtDiffGluon(0x0),
  JCComparisonVsPtDiffOther(0x0),
  JCComparisonVsPtDiffUnmatched(0x0),


  //ParticlePtAndJC(0x0),
  //ParticlePtAndJCUp(0x0),
  //ParticlePtAndJCDown(0x0),
  //ParticlePtAndJCGluon(0x0),
  //ParticlePtAndJCOther(0x0),
  //ParticlePtAndJCUnmatched(0x0),

  fTreeJets(0)
{

}

//________________________________________________________________________
AliAnalysisTaskJetChargeFlavourTemplates::AliAnalysisTaskJetChargeFlavourTemplates(const char *name) :
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fContainer(0),
  pChain(0),
  fPtThreshold(-9999.),
  fCentSelectOn(kTRUE),
  fCentMin(0),
  fCentMax(10),
  fJetRadius(0),
  JetChargeK(0.5),
  MotherFraction(0.8),
  JetMidPt(40),
  JetHighPt(80),


  fhJetPt(0x0),
  fhJetPhi(0x0),
  fhJetEta(0x0),

  JC(0x0),

  JCUp(0x0),
  JCDown(0x0),
  JCGluon(0x0),
  JCOther(0x0),
  JCUnmatched(0x0),

  JCLow(0x0),

  JCUpLow(0x0),
  JCDownLow(0x0),
  JCGluonLow(0x0),
  JCOtherLow(0x0),
  JCUnmatchedLow(0x0),


  JCMid(0x0),

  JCUpMid(0x0),
  JCDownMid(0x0),
  JCGluonMid(0x0),
  JCOtherMid(0x0),
  JCUnmatchedMid(0x0),

  JCHigh(0x0),

  JCUpHigh(0x0),
  JCDownHigh(0x0),
  JCGluonHigh(0x0),
  JCOtherHigh(0x0),
  JCUnmatchedHigh(0x0),

  fhParticleJetPt(0x0),
  fhParticleJetPhi(0x0),
  fhParticleJetEta(0x0),

  JCParticle(0x0),

  JCParticleUp(0x0),
  JCParticleDown(0x0),
  JCParticleGluon(0x0),
  JCParticleOther(0x0),
  JCParticleUnmatched(0x0),

  JCParticleLow(0x0),

  JCParticleUpLow(0x0),
  JCParticleDownLow(0x0),
  JCParticleGluonLow(0x0),
  JCParticleOtherLow(0x0),
  JCParticleUnmatchedLow(0x0),


  JCParticleMid(0x0),

  JCParticleUpMid(0x0),
  JCParticleDownMid(0x0),
  JCParticleGluonMid(0x0),
  JCParticleOtherMid(0x0),
  JCParticleUnmatchedMid(0x0),

  JCParticleHigh(0x0),

  JCParticleUpHigh(0x0),
  JCParticleDownHigh(0x0),
  JCParticleGluonHigh(0x0),
  JCParticleOtherHigh(0x0),
  JCParticleUnmatchedHigh(0x0),

  PtComparison(0x0),
  JCComparison(0x0),
  JCComparisonUp(0x0),
  JCComparisonDown(0x0),
  JCComparisonGluon(0x0),
  JCComparisonOther(0x0),
  JCComparisonUnmatched(0x0),

  Pt2DCompare(0x0),
  JC2DCompare(0x0),

  PtComparisonVsJCDiff(0x0),
  PtComparisonVsJCDiffUp(0x0),
  PtComparisonVsJCDiffDown(0x0),
  PtComparisonVsJCDiffGluon(0x0),
  PtComparisonVsJCDiffOther(0x0),
  PtComparisonVsJCDiffUnmatched(0x0),

  JCComparisonVsPtDiff(0x0),
  JCComparisonVsPtDiffUp(0x0),
  JCComparisonVsPtDiffDown(0x0),
  JCComparisonVsPtDiffGluon(0x0),
  JCComparisonVsPtDiffOther(0x0),
  JCComparisonVsPtDiffUnmatched(0x0),


  //ParticlePtAndJC(0x0),
  //ParticlePtAndJCUp(0x0),
  //ParticlePtAndJCDown(0x0),
  //ParticlePtAndJCGluon(0x0),
  //ParticlePtAndJCOther(0x0),
  //ParticlePtAndJCUnmatched(0x0),


  fTreeJets(0)
{
  // Standard constructor.
  for(Int_t i=0;i<nBranchesJetChargeFlavourTemplates;i++){
    fTreeBranch[i]=0;
  }
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskJetChargeFlavourTemplates::~AliAnalysisTaskJetChargeFlavourTemplates()
{
  // Destructor.
}

//________________________________________________________________________
 void AliAnalysisTaskJetChargeFlavourTemplates::UserCreateOutputObjects()
{
  // Echo jet radius
  //Info("TaskJets","Using jet radius R=%f",fJetRadius);

  // Create user output.
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  TH1::AddDirectory(oldStatus);
  //create output TTree
  const char* nameoutput = GetOutputSlot(2)->GetContainer()->GetName();
  fTreeJets = new TTree(nameoutput, nameoutput);
  // Names for the branches
  TString *fTreeBranchName = new TString [nBranchesJetChargeFlavourTemplates];

  // Name the branches of your TTree here
  fTreeBranchName[0]  = "Pt";
  fTreeBranchName[1]  = "Phi";
  fTreeBranchName[2]  = "Eta";

  fTreeBranchName[3]  = "JetCharge";

  fTreeBranchName[4] = "JCUp";
  fTreeBranchName[5] = "JCDown";
  fTreeBranchName[6] = "JCGluon";
  fTreeBranchName[7] = "JCOther";
  fTreeBranchName[8] = "JCUnmatched";

  fTreeBranchName[9]  = "ParticlePt";
  fTreeBranchName[10]  = "ParticlePhi";
  fTreeBranchName[11]  = "ParticleEta";

  fTreeBranchName[12]  = "ParticleJetCharge";

  fTreeBranchName[13] = "ParticleJCUp";
  fTreeBranchName[14] = "ParticleJCDown";
  fTreeBranchName[15] = "ParticleJCGluon";
  fTreeBranchName[16] = "ParticleJCOther";
  fTreeBranchName[17] = "ParticleJCUnmatched";

  fTreeBranchName[18] = "PtComparison";
  fTreeBranchName[19] = "JCComparison";
  fTreeBranchName[20] = "JCComparisonUp";
  fTreeBranchName[21] = "JCComparisonDown";
  fTreeBranchName[22] = "JCComparisonGluon";
  fTreeBranchName[23] = "JCComparisonOther";
  fTreeBranchName[24] = "JCComparisonUnmatched";




  // Associate the branches
  for(Int_t iBranch=0; iBranch < nBranchesJetChargeFlavourTemplates; iBranch++){
    cout<<"looping over variables"<<endl;
    fTreeJets->Branch(fTreeBranchName[iBranch].Data(), &fTreeBranch[iBranch], Form("%s/D", fTreeBranchName[iBranch].Data()));
  }

  // Define histograms
  fhJetPt= new TH1F("fhJetPt", "Jet Pt",1500,-0.5,149.5 );
  fOutput->Add(fhJetPt);
  fhJetPhi= new TH1F("fhJetPhi", "Jet Phi",360 , -1.5*(TMath::Pi()), 1.5*(TMath::Pi()));
  fOutput->Add(fhJetPhi);
  fhJetEta= new TH1F("fhJetEta", "Jet Eta",100,-2,2);
  fOutput->Add(fhJetEta);

  /*
  fhEventCounter= new TH1F("fhEventCounter", "Event Counter",10,10,10);
  fOutput->Add(fhEventCounter);

  fhRunNumberCounter= new TH1F("fhRunNumberCounter", "Event Counter",10,10,10);
  fOutput->Add(fhRunNumberCounter);
  */

  JC= new TH1F("JC", "Jet Charge", 25, -3, 3);
  fOutput->Add(JC);

  JCUp= new TH1F("JCUp", "Jet Charge Up", 25, -3, 3);
  fOutput->Add(JCUp);
  JCDown= new TH1F("JCDown", "Jet Charge Down", 25, -3, 3);
  fOutput->Add(JCDown);
  JCGluon= new TH1F("JCGluon", "Jet Charge Gluon", 25, -3, 3);
  fOutput->Add(JCGluon);
  JCOther= new TH1F("JCOther", "Jet Charge Other", 25, -3, 3);
  fOutput->Add(JCOther);
  JCUnmatched= new TH1F("JCUnmatched", "Jet Charge Unmatched", 25, -3, 3);
  fOutput->Add(JCUnmatched);




  JCLow= new TH1F("JCLow", "Jet Charge Low Pt ", 25, -3, 3);
  fOutput->Add(JCLow);

  JCUpLow= new TH1F("JCUpLow", "Jet Charge Up Low Pt ", 25, -3, 3);
  fOutput->Add(JCUpLow);
  JCDownLow= new TH1F("JCDownLow", "Jet Charge Down Low Pt", 25, -3, 3);
  fOutput->Add(JCDownLow);
  JCGluonLow= new TH1F("JCGluonLow", "Jet Charge Gluon Low Pt", 25, -3, 3);
  fOutput->Add(JCGluonLow);
  JCOtherLow= new TH1F("JCOtherLow", "Jet Charge Other Low Pt", 25, -3, 3);
  fOutput->Add(JCOtherLow);
  JCUnmatchedLow= new TH1F("JCUnmatchedLow", "Jet Charge Unmatched Low Pt", 25, -3, 3);
  fOutput->Add(JCUnmatchedLow);



  JCMid= new TH1F("JCMid", "Jet Charge Mid Pt ", 25, -3, 3);
  fOutput->Add(JCMid);

  JCUpMid= new TH1F("JCUpMid", "Jet Charge Up Mid Pt ", 25, -3, 3);
  fOutput->Add(JCUpMid);
  JCDownMid= new TH1F("JCDownMid", "Jet Charge Down Mid Pt", 25, -3, 3);
  fOutput->Add(JCDownMid);
  JCGluonMid= new TH1F("JCGluonMid", "Jet Charge Gluon Mid Pt", 25, -3, 3);
  fOutput->Add(JCGluonMid);
  JCOtherMid= new TH1F("JCOtherMid", "Jet Charge Other Mid Pt", 25, -3, 3);
  fOutput->Add(JCOtherMid);
  JCUnmatchedMid= new TH1F("JCUnmatchedMid", "Jet Charge Unmatched Mid Pt", 25, -3, 3);
  fOutput->Add(JCUnmatchedMid);

  JCHigh= new TH1F("JCHigh", "Jet Charge High Pt ", 25, -3, 3);
  fOutput->Add(JCHigh);

  JCUpHigh= new TH1F("JCUpHigh", "Jet Charge Up High Pt ", 25, -3, 3);
  fOutput->Add(JCUpHigh);
  JCDownHigh= new TH1F("JCDownHigh", "Jet Charge Down High Pt", 25, -3, 3);
  fOutput->Add(JCDownHigh);
  JCGluonHigh= new TH1F("JCGluonHigh", "Jet Charge Gluon High Pt", 25, -3, 3);
  fOutput->Add(JCGluonHigh);
  JCOtherHigh= new TH1F("JCOtherHigh", "Jet Charge Other High Pt", 25, -3, 3);
  fOutput->Add(JCOtherHigh);
  JCUnmatchedHigh= new TH1F("JCUnmatchedHigh", "Jet Charge Unmatched High Pt", 25, -3, 3);
  fOutput->Add(JCUnmatchedHigh);

  fhParticleJetPt= new TH1F("fhParticleJetPt", "Jet Pt",1500,-0.5,149.5 );
  fOutput->Add(fhParticleJetPt);
  fhParticleJetPhi= new TH1F("fhParticleJetPhi", "Jet Phi",360 , -1.5*(TMath::Pi()), 1.5*(TMath::Pi()));
  fOutput->Add(fhParticleJetPhi);
  fhParticleJetEta= new TH1F("fhParticleJetEta", "Jet Eta",100,-2,2);
  fOutput->Add(fhParticleJetEta);


// Add Particle Jet Charge

  JCParticle= new TH1F("JCParticle", "Jet Charge", 25, -3, 3);
  fOutput->Add(JCParticle);

  JCParticleUp= new TH1F("JCParticleUp", "Jet Charge Up", 25, -3, 3);
  fOutput->Add(JCParticleUp);
  JCParticleDown= new TH1F("JCParticleDown", "Jet Charge Down", 25, -3, 3);
  fOutput->Add(JCParticleDown);
  JCParticleGluon= new TH1F("JCParticleGluon", "Jet Charge Gluon", 25, -3, 3);
  fOutput->Add(JCParticleGluon);
  JCParticleOther= new TH1F("JCParticleOther", "Jet Charge Other", 25, -3, 3);
  fOutput->Add(JCParticleOther);
  JCParticleUnmatched= new TH1F("JCParticleUnmatched", "Jet Charge Unmatched", 25, -3, 3);
  fOutput->Add(JCParticleUnmatched);



  JCParticleLow= new TH1F("JCParticleLow", "Jet Charge Low Pt ", 25, -3, 3);
  fOutput->Add(JCParticleLow);

  JCParticleUpLow= new TH1F("JCParticleUpLow", "Jet Charge Up Low Pt ", 25, -3, 3);
  fOutput->Add(JCParticleUpLow);
  JCParticleDownLow= new TH1F("JCParticleDownLow", "Jet Charge Down Low Pt", 25, -3, 3);
  fOutput->Add(JCParticleDownLow);
  JCParticleGluonLow= new TH1F("JCParticleGluonLow", "Jet Charge Gluon Low Pt", 25, -3, 3);
  fOutput->Add(JCParticleGluonLow);
  JCParticleOtherLow= new TH1F("JCParticleOtherLow", "Jet Charge Other Low Pt", 25, -3, 3);
  fOutput->Add(JCParticleOtherLow);
  JCParticleUnmatchedLow= new TH1F("JCParticleUnmatchedLow", "Jet Charge Unmatched Low Pt", 25, -3, 3);
  fOutput->Add(JCParticleUnmatchedLow);



  JCParticleMid= new TH1F("JCParticleMid", "Jet Charge Mid Pt ", 25, -3, 3);
  fOutput->Add(JCParticleMid);

  JCParticleUpMid= new TH1F("JCParticleUpMid", "Jet Charge Up Mid Pt ", 25, -3, 3);
  fOutput->Add(JCParticleUpMid);
  JCParticleDownMid= new TH1F("JCParticleDownMid", "Jet Charge Down Mid Pt", 25, -3, 3);
  fOutput->Add(JCParticleDownMid);
  JCParticleGluonMid= new TH1F("JCParticleGluonMid", "Jet Charge Gluon Mid Pt", 25, -3, 3);
  fOutput->Add(JCParticleGluonMid);
  JCParticleOtherMid= new TH1F("JCParticleOtherMid", "Jet Charge Other Mid Pt", 25, -3, 3);
  fOutput->Add(JCParticleOtherMid);
  JCParticleUnmatchedMid= new TH1F("JCParticleUnmatchedMid", "Jet Charge Unmatched Mid Pt", 25, -3, 3);
  fOutput->Add(JCParticleUnmatchedMid);

  JCParticleHigh= new TH1F("JCParticleHigh", "Jet Charge High Pt ", 25, -3, 3);
  fOutput->Add(JCParticleHigh);

  JCParticleUpHigh= new TH1F("JCParticleUpHigh", "Jet Charge Up High Pt ", 25, -3, 3);
  fOutput->Add(JCParticleUpHigh);
  JCParticleDownHigh= new TH1F("JCParticleDownHigh", "Jet Charge Down High Pt", 25, -3, 3);
  fOutput->Add(JCParticleDownHigh);
  JCParticleGluonHigh= new TH1F("JCParticleGluonHigh", "Jet Charge Gluon High Pt", 25, -3, 3);
  fOutput->Add(JCParticleGluonHigh);
  JCParticleOtherHigh= new TH1F("JCParticleOtherHigh", "Jet Charge Other High Pt", 25, -3, 3);
  fOutput->Add(JCParticleOtherHigh);
  JCParticleUnmatchedHigh= new TH1F("JCParticleUnmatchedHigh", "Jet Charge Unmatched High Pt", 25, -3, 3);
  fOutput->Add(JCParticleUnmatchedHigh);

  //Add Comparison Histograms

  PtComparison= new TH1F("PtComparison", "Pt Diffrence",100,-2,2 );
  fOutput->Add(PtComparison);
  JCComparison = new TH1F("JCComparison", "Jet Charge Comparison", 100, -3, 3);
  fOutput->Add(JCComparison);

  JCComparisonUp = new TH1F("JCComparisonUp", "Jet Charge Comparison Up", 100, -3, 3);
  fOutput->Add(JCComparisonUp);
  JCComparisonDown = new TH1F("JCComparisonDown", "Jet Charge Comparison Down", 100, -3, 3);
  fOutput->Add(JCComparisonDown);
  JCComparisonGluon = new TH1F("JCComparisonGluon", "Jet Charge Comparison Gluon", 100, -3, 3);
  fOutput->Add(JCComparisonGluon);
  JCComparisonOther = new TH1F("JCComparisonOther", "Jet Charge Comparison Other", 100, -3, 3);
  fOutput->Add(JCComparisonOther);
  JCComparisonUnmatched = new TH1F("JCComparisonUnmatched", "Jet Charge Comparison Unmatched", 100, -3, 3);
  fOutput->Add(JCComparisonUnmatched);

  Pt2DCompare = new TH2F("Pt2DCompare", "Pt 2D Particle vs Detector",1500,-0.5,149.5 ,1500,-0.5,149.5);
  fOutput->Add(Pt2DCompare);
  JC2DCompare = new TH2F("JC2DCompare", "JC 2D Particle vs Detector",150, -4, 4, 50, -4, 4);
  fOutput->Add(JC2DCompare);


  PtComparisonVsJCDiff = new TH2F("PtComparisonVsJCDiff", "Pt Det Vs Part Compared to JC",50,-0.5,149.5,50,-0.5,149.5);
  fOutput->Add(PtComparisonVsJCDiff);

  PtComparisonVsJCDiffUp = new TH2F("PtComparisonVsJCDiffUp", "Pt Det Vs Part Compared to JC Up",50,-0.5,149.5,50,-0.5,149.5);
  fOutput->Add(PtComparisonVsJCDiffUp);
  PtComparisonVsJCDiffDown = new TH2F("PtComparisonVsJCDiffDown", "Pt Det Vs Part Compared to JC Down",50,-0.5,149.5,50,-0.5,149.5);
  fOutput->Add(PtComparisonVsJCDiffDown);
  PtComparisonVsJCDiffGluon = new TH2F("PtComparisonVsJCDiffGluon", "Pt Det Vs Part Compared to JC Gluon",50,-0.5,149.5,50,-0.5,149.5);
  fOutput->Add(PtComparisonVsJCDiffGluon);
  PtComparisonVsJCDiffOther = new TH2F("PtComparisonVsJCDiffOther", "Pt Det Vs Part Compared to JC Other",50,-0.5,149.5,50,-0.5,149.5);
  fOutput->Add(PtComparisonVsJCDiffOther);
  PtComparisonVsJCDiffUnmatched = new TH2F("PtComparisonVsJCDiffUnmatched", "Pt Det Vs Part Compared to JC Unmatched",50,-0.5,149.5,50,-0.5,149.5);
  fOutput->Add(PtComparisonVsJCDiffUnmatched);

  JCComparisonVsPtDiff = new TH2F("JCComparisonVsPtDiff", "JC Det vs Part Compared to Pt", 50, -4, 4, 50, -4, 4);
  fOutput->Add(JCComparisonVsPtDiff);

  JCComparisonVsPtDiffUp = new TH2F("JCComparisonVsPtDiffUp", "JC Det vs Part Compared to Pt Up", 50, -4, 4, 50, -4, 4);
  fOutput->Add(JCComparisonVsPtDiffUp);
  JCComparisonVsPtDiffDown = new TH2F("JCComparisonVsPtDiffDown", "JC Det vs Part Compared to Pt Down", 50, -4, 4, 50, -4, 4);
  fOutput->Add(JCComparisonVsPtDiffDown);
  JCComparisonVsPtDiffGluon = new TH2F("JCComparisonVsPtDiffGluon", "JC Det vs Part Compared to Pt Gluon", 50, -4, 4, 50, -4, 4);
  fOutput->Add(JCComparisonVsPtDiffGluon);
  JCComparisonVsPtDiffOther = new TH2F("JCComparisonVsPtDiffOther", "JC Det vs Part Compared to Pt Other", 50, -4, 4, 50, -4, 4);
  fOutput->Add(JCComparisonVsPtDiffOther);
  JCComparisonVsPtDiffUnmatched = new TH2F("JCComparisonVsPtDiffUnmatched", "JC Det vs Part Compared to Pt Unmatched", 30, -4, 4, 30, -4, 4);
  fOutput->Add(JCComparisonVsPtDiffUnmatched);

/*
  ParticlePtAndJC = new TH3F("ParticlePtAndJC", "Pt Det vs Pt Part Vs JC Part Plot to JC Det",1500,-0.5,149.5,1500,-0.5,149.5, 30, -4, 4);
  fOutput->Add(ParticlePtAndJC);

  ParticlePtAndJCUp = new TH3F("ParticlePtAndJCUp", "Pt Det vs Pt Part Vs JC Part Plot to JC Det - Up",1500,-0.5,149.5,1500,-0.5,149.5, 30, -4, 4);
  fOutput->Add(ParticlePtAndJCUp);
  ParticlePtAndJCDown = new TH3F("ParticlePtAndJCDown", "Pt Det vs Pt Part Vs JC Part Plot to JC Det - Down",1500,-0.5,149.5,1500,-0.5,149.5, 30, -4, 4);
  fOutput->Add(ParticlePtAndJCDown);
  ParticlePtAndJCGluon = new TH3F("ParticlePtAndJCGluon", "Pt Det vs Pt Part Vs JC Part Plot to JC Det - Gluon",1500,-0.5,149.5,1500,-0.5,149.5, 30, -4, 4);
  fOutput->Add(ParticlePtAndJCGluon);
  ParticlePtAndJCOther = new TH3F("ParticlePtAndJCOther", "Pt Det vs Pt Part Vs JC Part Plot to JC Det - Other",1500,-0.5,149.5,1500,-0.5,149.5, 30, -4, 4);
  fOutput->Add(ParticlePtAndJCOther);
  ParticlePtAndJCUnmatched = new TH3F("ParticlePtAndJCUnmatched", "Pt Det vs Pt Part Vs JC Part Plot to JC Det - Unmatched",1500,-0.5,149.5,1500,-0.5,149.5, 30, -4, 4);
  fOutput->Add(ParticlePtAndJCUnmatched);
*/





  // Make sure that the outputs get written out
  PostData(1,fOutput);
  PostData(2,fTreeJets);
  // delete [] fShapesVarNames;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetChargeFlavourTemplates::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().




  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetChargeFlavourTemplates::FillHistograms()
{
  // Centrality selection, if enabled
  if (fCentSelectOn){
    if ((fCent>fCentMax) || (fCent<fCentMin))
      return 0;
  }
  // Initialise jet pointer
  //cout << "Running Fill Histograms" << endl;
  AliEmcalJet *Jet1 = NULL; //Original Jet in the event                                                                                                         // Get jet container (0 = ?)
  AliJetContainer *JetCont= GetJetContainer(0); //Jet Container for event
  AliEmcalJet *TruthJet = NULL; //Original Jet in the event                                                                                                         // Get jet container (0 = ?)
  AliJetContainer *JetGen= GetJetContainer(1); //Jet Container for event


  AliParticleContainer *MCParticleContainer = JetGen->GetParticleContainer();
  Int_t nAcceptedJets = JetCont->GetNAcceptedJets();
  //TClonesArray *trackArr = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("HybridTracks"));
  Double_t JetPhi=0;
  Double_t JetParticlePhi=0;
  Double_t JetPt_ForThreshold=0;

  if(JetCont) {
    // Technical detail; fix possibly corrupted jet container ID
    JetCont->ResetCurrentID();
    // Jet is acceptable?
    while((Jet1=JetCont->GetNextAcceptJet())) {
      if(!Jet1)
      {

        continue;
      }

      // Jet is above threshold?

      //cout << "Running A Jet" << endl;

      // Get the jet constituents

      //AliParticleContainer *fMCContainer = GetParticleContainer(1);
      //UInt_t nMCConstituents = fMCContainer->GetNParticles();

      //AliParticleContainer * MCTracks  = dynamic_cast<AliParticleContainer *> (GetParticleContainer("embeddedTracks"));
      //UInt_t nMCTracks = MCTracks->GetNParticles();

      /*
      TIter nextPartColl(&fParticleCollArray);
      AliParticleContainer* tracks = 0;
      while ((tracks = static_cast<AliParticleContainer*>(nextPartColl())))
      {
      AliParticleContainer* fMCContainer = tracks;
      }
      */

      AliParticleContainer *fTrackCont = JetCont->GetParticleContainer();
      UInt_t nJetConstituents = Jet1->GetNumberOfTracks();



      //cout << nTest << "::::" << nMCConstituents << "::::" << nJetConstituents << endl;


      // Must have at least two constituents
      if( nJetConstituents < 2 )
      {
        continue;
      }


      JetPt_ForThreshold = Jet1->Pt();

      //Check for leading track Pt to reduce cominatorial jets.

      if(Jet1->GetLeadingTrack()->Pt() < 5)
      {

          continue;
          //std::cout << "LEADING TRACK TO SMALL!!!" << std::endl;
      }

      if(JetPt_ForThreshold<fPtThreshold)
      {

        continue;
      }
      else {

        // Initialising the Tagged PYTHIA Jet.
        Bool_t kHasTruthJet = kFALSE;
        AliEmcalJet* TruthJet;
        UInt_t nTruthConstituents;

        if(Jet1->ClosestJet() != NULL)
        {
          kHasTruthJet = kTRUE;
          TruthJet = Jet1->ClosestJet();
          nTruthConstituents = TruthJet->GetNumberOfTracks();

        }


      	// Add other branches below this line ...
        // Initialise Jet shapes
        Int_t fPdgCodes[100] = {};
        Int_t fCurrentPdg = 0;
        Int_t nMothers = 0;
        Double_t jetCharge = 0;
        Double_t jetChargeParticle = 0;

        //Finding Pdg Code using TRUTH Jet.


        if(kHasTruthJet)
        {

          if(TruthJet->Pt()<fPtThreshold)
          {

            continue;
          }

          // Filling the histograms here
          fhJetPt->Fill(Jet1->Pt());
          fhParticleJetPt->Fill(TruthJet->Pt());

          JetPhi=Jet1->Phi();
          if(JetPhi < -1*TMath::Pi())
            JetPhi += (2*TMath::Pi());
          else if (JetPhi > TMath::Pi())
            JetPhi -= (2*TMath::Pi());
          fhJetPhi->Fill(JetPhi);

          JetParticlePhi=TruthJet->Phi();
          if(JetParticlePhi < -1*TMath::Pi())
            JetParticlePhi += (2*TMath::Pi());
          else if (JetParticlePhi > TMath::Pi())
            JetParticlePhi -= (2*TMath::Pi());
          fhParticleJetPhi->Fill(JetParticlePhi);

          fhJetEta->Fill(Jet1->Eta());
          fhParticleJetEta->Fill(TruthJet->Eta());
          // Filling the TTree branch(es) here
          Double_t JetPt;

          JetPt = Jet1->Pt();
          fTreeBranch[0]=Jet1->Pt();
          fTreeBranch[1]=JetPhi;
          fTreeBranch[2]=Jet1->Eta();

          fTreeBranch[9]=TruthJet->Pt();
          fTreeBranch[10]=TruthJet->Phi();
          fTreeBranch[11]=TruthJet->Eta();

          Double_t PtDiff = (Jet1->Pt() - TruthJet->Pt())/TruthJet->Pt();

          PtComparison->Fill(PtDiff);
          fTreeBranch[18]= PtDiff;

          Pt2DCompare->Fill(Jet1->Pt(),TruthJet->Pt());


          //cout << "Has Matched Jet " << endl;
/*
          std::cout << "Jet Phi: "  << Jet1->Phi() << std::endl;
          std::cout << "Jet Eta: "  << Jet1->Eta() << std::endl;
          std::cout << "Jet Energy: "  << Jet1->E() << std::endl;
          std::cout << "Jet Pt:  "  << Jet1->Pt() << std::endl << std::endl;
*/
          for (UInt_t iTruthConst = 0; iTruthConst < nTruthConstituents; iTruthConst++ )
          {
            AliMCParticle* TruthParticle = (AliMCParticle*) TruthJet->Track(iTruthConst);

            AliMCParticle *MotherParticle = (AliMCParticle*)  MCParticleContainer->GetParticle(TruthParticle->GetMother());

            //cout << "Diff Phi: " << abs(MotherParticle->Phi() - Jet1->Phi()) << endl;
            //cout << "Diff Eta: " << abs(MotherParticle->Eta() - Jet1->Eta()) << endl;
            //cout << "Diff R: " << pow(pow((TruthParticle->Phi() - Jet1->Phi()),2) + pow((MotherParticle->Eta() - Jet1->Eta()),2),0.5)  << endl;
            //cout << "Diff Pt: " << abs(TruthParticle->Pt() - Jet1->Pt()) << endl;
            //cout << "Particle Pdg: " << TruthParticle->PdgCode() << endl;
            //cout << "Particle Label: " << TruthParticle->GetLabel() << endl;
            //cout << "Loop Counter : " << iTruthConst << endl  << endl;

            while(MotherParticle->GetMother() > 0)
            {
              MotherParticle = (AliMCParticle*) MCParticleContainer->GetParticle(MotherParticle->GetMother());
              //cout <<"Mother Get Result: " << MotherParticle->GetMother() << endl;
              //cout << "Particle Label: " << MotherParticle->Label() << endl;
              //cout << "Particle Pdg: " << MotherParticle->PdgCode() << endl;
/*
              cout << "Particle Phi: " << MotherParticle->Phi() << endl;
              cout << "Particle Eta: " << MotherParticle->Eta() << endl;
              cout << "Particle Pt: " << MotherParticle->Pt() << endl;
              cout << "Particle Pdg: " << MotherParticle->PdgCode() << endl << endl;
*/
            }

            if(MotherParticle->E() < 3400)
            {
            fCurrentPdg = MotherParticle->PdgCode();
            fPdgCodes[nMothers] = fCurrentPdg;
            nMothers++;
            /*
            cout << "Particle Phi: " << MotherParticle->Phi() << endl;
            cout << "Particle Eta: " << MotherParticle->Eta() << endl;
            cout << "Particle Pt: " << MotherParticle->Pt() << endl;
            cout << "Particle Energy: " << MotherParticle->E() << endl;
            cout << "Particle Pdg: " << MotherParticle->PdgCode() << endl << endl;
            */
            }
          }


            // New Techniqe for determing Jet PDG Code

            Int_t UniquePdgCodes[20] = {};            //To be filled, maximium is that there are  20 uniques
            Double_t UniquePdgFrequency[20] = {};
            Int_t nUniques = 0;

            //Loop of PDG code found
            for(int i = 0; i < nMothers; i++)
            {
              // consider one pdgCode at the time.
              fCurrentPdg = fPdgCodes[i];
              // Loop over unique arry to be filled
              for(int j = 0; j < nMothers; j++)
              {
                // If it hasnt matched and the current unique value is empty list the new value and increment frequncy by 1 and then break
                if(UniquePdgCodes[j] != fCurrentPdg && UniquePdgCodes[j] == 0)
                {
                  UniquePdgCodes[j] = fCurrentPdg;
                  UniquePdgFrequency[j] = UniquePdgFrequency[j] + 1.;
                  nUniques ++;
                  break;
                }
                //Check if the PDG is already lsited if it matched increase the frequency counter by 1
                else if(UniquePdgCodes[j] == fCurrentPdg)
                {
                  UniquePdgFrequency[j] = UniquePdgFrequency[j] + 1.;
                  break;
                }

                // Otherwise the PDG hasnt matched and the value isnt zero check the next value

              }
            }

            // Setting Final Pdg Code
            // normalising frequency array

            for(unsigned int i = 0; i < nUniques; i++)
            {
              UniquePdgFrequency[i] = UniquePdgFrequency[i]/nMothers;
            }

            //Find the index of the max
            int IndexOfMaximum = -1;

            // assigens index of maximum if the over limit factation of mother particles agree

            for(unsigned int i = 0; i < nUniques; i++)
            {
              if(UniquePdgFrequency[i] > MotherFraction)
              {
                IndexOfMaximum = i;
              }
            }

            // Set final PDG Code to be used
            fCurrentPdg = fPdgCodes[IndexOfMaximum];

            if(IndexOfMaximum < 0)
            {
              fCurrentPdg = 0;
            }

                  //Outputs for Checking
/*
                  if(nUniques > 1)
                  {
                    // output PDG Codes for testin
                    cout << "PdgCodes : [" ;
                    for(unsigned int i = 0; i < nMothers; i++)
                    {
                      cout << fPdgCodes[i] << ",";
                    }
                    cout << "]" << endl;

                    // output Unique PDG Codes for testing

                    cout << "Unique PdgCodes : [" ;
                    for(unsigned int i = 0; i < nMothers; i++)
                    {
                      cout << UniquePdgCodes[i] << ",";
                    }
                    cout << "]" << endl;

                    // ouput PDG Codes fort testing

                    cout << "Frequency PdgCodes : [" ;
                    for(unsigned int i = 0; i < nMothers; i++)
                    {
                      cout << UniquePdgFrequency[i] << ",";
                    }
                    cout << "]" << endl;

                    //cout << "Number of Mothers: " << nMothers << endl;
                    //cout << "Number of Uniques: " << nUniques << endl;
                    cout << "Final PDG Choice: " << fCurrentPdg << endl << endl;
                  }
*/

          // Loop over the consituents
          for (UInt_t iJetConst = 0; iJetConst < nJetConstituents; iJetConst++ )
          {
            AliVParticle *JetParticle = Jet1->Track(iJetConst);
            jetCharge += JetParticle->Charge()*pow(JetParticle->Pt(),JetChargeK);
          }

          for (UInt_t iTruthConst = 0; iTruthConst < nTruthConstituents; iTruthConst++ )
          {
            AliMCParticle* TruthParticle = (AliMCParticle*) TruthJet->Track(iTruthConst);
            jetChargeParticle += TruthParticle->Charge()/3*pow(TruthParticle->Pt(),JetChargeK);
          }



        // Normalise the Non Flavoured Jet CHarge

        jetCharge/=pow(Jet1->Pt(),0.5);

        // Normalise Particle level jet charge

        jetChargeParticle/=pow(TruthJet->Pt(),0.5);


        //Put The Jet Charge in the right place
        fTreeBranch[3] = jetCharge;
        JC->Fill(jetCharge);
        fTreeBranch[12] = jetChargeParticle;
        JCParticle->Fill(jetChargeParticle);

        Double_t JetChargeDiff = (jetCharge - jetChargeParticle)/jetChargeParticle;

        JCComparison->Fill(JetChargeDiff);
        fTreeBranch[19] = JetChargeDiff;

        JC2DCompare->Fill(jetCharge,jetChargeParticle);

        //cout << JetChargeDiff << endl;

        PtComparisonVsJCDiff->Fill(Jet1->Pt(),TruthJet->Pt(),JetChargeDiff);

        JCComparisonVsPtDiff->Fill(jetCharge,jetChargeParticle,PtDiff);

        //Split the Jet in to apporpate momentum bin.

        if(JetPt < JetMidPt)
        {
          JCLow->Fill(jetCharge);
          JCParticleLow->Fill(jetChargeParticle);
        }
        else if( JetPt > JetMidPt && JetPt < JetHighPt)
        {
          JCMid->Fill(jetCharge);
          JCParticleMid->Fill(jetChargeParticle);
        }
        else
        {
          JCHigh->Fill(jetCharge);
          JCParticleHigh->Fill(jetChargeParticle);
        }




        //Add Up JetCharge
        if(fCurrentPdg == 2)
        {
          fTreeBranch[4] = jetCharge;
          JCUp->Fill(jetCharge);
          fTreeBranch[13] = jetChargeParticle;
          JCParticleUp->Fill(jetChargeParticle);

          JCComparisonUp->Fill(JetChargeDiff);
          fTreeBranch[20] = JetChargeDiff;

          PtComparisonVsJCDiffUp->Fill(Jet1->Pt(),TruthJet->Pt(),JetChargeDiff);

          JCComparisonVsPtDiffUp->Fill(jetCharge,jetChargeParticle,PtDiff);

            if(JetPt < JetMidPt)
            {
              JCUpLow->Fill(jetCharge);
              JCParticleUpLow->Fill(jetChargeParticle);
            }
            else if( JetPt > JetMidPt && JetPt < JetHighPt)
            {
              JCUpMid->Fill(jetCharge);
              JCParticleUpMid->Fill(jetChargeParticle);
            }
            else
            {
              JCUpHigh->Fill(jetCharge);
              JCParticleUpHigh->Fill(jetChargeParticle);
            }




        }
        //Add Down JetCharge
        else if(fCurrentPdg == 1)
        {
          fTreeBranch[5] = jetCharge;
          JCDown->Fill(jetCharge);
          fTreeBranch[14] = jetChargeParticle;
          JCParticleDown->Fill(jetChargeParticle);

          JCComparisonDown->Fill(JetChargeDiff);
          fTreeBranch[21] = JetChargeDiff;

          PtComparisonVsJCDiffDown->Fill(Jet1->Pt(),TruthJet->Pt(),JetChargeDiff);

          JCComparisonVsPtDiffDown->Fill(jetCharge,jetChargeParticle,PtDiff);

            if(JetPt < JetMidPt)
            {
              JCDownLow->Fill(jetCharge);
              JCParticleDownLow->Fill(jetChargeParticle);
            }
            else if( JetPt > JetMidPt && JetPt < JetHighPt)
            {
              JCDownMid->Fill(jetCharge);
              JCParticleDownMid->Fill(jetChargeParticle);
            }
            else
            {
              JCDownHigh->Fill(jetCharge);
              JCParticleDownHigh->Fill(jetChargeParticle);
            }
        }




        //Add Gluon JetCharge
        else if(fCurrentPdg == 21)
        {
          fTreeBranch[6] = jetCharge;
          JCGluon->Fill(jetCharge);
          fTreeBranch[15] = jetChargeParticle;
          JCParticleGluon->Fill(jetChargeParticle);

          JCComparisonGluon->Fill(JetChargeDiff);
          fTreeBranch[22] = JetChargeDiff;

          PtComparisonVsJCDiffGluon->Fill(Jet1->Pt(),TruthJet->Pt(),JetChargeDiff);

          JCComparisonVsPtDiffGluon->Fill(jetCharge,jetChargeParticle,PtDiff);

            if(JetPt < JetMidPt)
            {
              JCGluonLow->Fill(jetCharge);
              JCParticleGluonLow->Fill(jetChargeParticle);
            }
            else if( JetPt > JetMidPt && JetPt < JetHighPt)
            {
              JCGluonMid->Fill(jetCharge);
              JCParticleGluonMid->Fill(jetChargeParticle);
            }
            else
            {
              JCGluonHigh->Fill(jetCharge);
              JCParticleGluonHigh->Fill(jetChargeParticle);
            }


        }



        //Add Unmatched JetCharge Catagory
        else if(IndexOfMaximum == -1)
        {

          fTreeBranch[7] = jetCharge;
          JCUnmatched->Fill(jetCharge);
          fTreeBranch[17] = jetChargeParticle;
          JCParticleUnmatched->Fill(jetChargeParticle);

          JCComparisonUnmatched->Fill(JetChargeDiff);
          fTreeBranch[24] = JetChargeDiff;

          PtComparisonVsJCDiffUnmatched->Fill(Jet1->Pt(),TruthJet->Pt(),JetChargeDiff);

          JCComparisonVsPtDiffUnmatched->Fill(jetCharge,jetChargeParticle,PtDiff);


            if(JetPt < JetMidPt)
            {
              JCUnmatchedLow->Fill(jetCharge);
              JCParticleUnmatchedLow->Fill(jetChargeParticle);
            }
            else if( JetPt > JetMidPt && JetPt < JetHighPt)
            {
              JCUnmatchedMid->Fill(jetCharge);
              JCParticleUnmatchedMid->Fill(jetChargeParticle);
            }
            else
            {
              JCUnmatchedHigh->Fill(jetCharge);
              JCParticleUnmatchedHigh->Fill(jetChargeParticle);
            }


        }



        //Adding Other Flavour JetCharge
        else
        {
          fTreeBranch[8] = jetCharge;
          JCOther->Fill(jetCharge);
          fTreeBranch[16] = jetChargeParticle;
          JCParticleOther->Fill(jetChargeParticle);

          JCComparisonOther->Fill(JetChargeDiff);
          fTreeBranch[23] = JetChargeDiff;

          PtComparisonVsJCDiffOther->Fill(Jet1->Pt(),TruthJet->Pt(),JetChargeDiff);

          JCComparisonVsPtDiffOther->Fill(jetCharge,jetChargeParticle,PtDiff);


            if(JetPt < JetMidPt)
            {
              JCOtherLow->Fill(jetCharge);
              JCParticleOtherLow->Fill(jetChargeParticle);
            }
            else if( JetPt > JetMidPt && JetPt < JetHighPt)
            {
              JCOtherMid->Fill(jetCharge);
              JCParticleOtherMid->Fill(jetChargeParticle);
            }
            else
            {
              JCOtherHigh->Fill(jetCharge);
              JCParticleOtherHigh->Fill(jetChargeParticle);
            }

        }






        fTreeJets->Fill();




      }

      //cout << "End of Jet" << endl;
      }
    }
  }
  return kTRUE;
}



//________________________________________________________________________
Bool_t AliAnalysisTaskJetChargeFlavourTemplates::RetrieveEventObjects() {
  //
  // retrieve event objects
  //
  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  return kTRUE;
}


//_______________________________________________________________________
void AliAnalysisTaskJetChargeFlavourTemplates::Terminate(Option_t *)
{
  // Called once at the end of the analysis.


  // Normalise historgrams over number of Jets considered

}
