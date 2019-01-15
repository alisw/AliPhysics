//
// Basic analysis task.
//
// Basic analysis task template for analysis jets storing information in both tree
// branches and histograms

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
#include <vector>
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
#include "AliAnalysisTaskRecursiveSoftDrop.h"

#include "FJ_includes.h"

//Globals
using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskRecursiveSoftDrop)

//________________________________________________________________________
AliAnalysisTaskRecursiveSoftDrop::AliAnalysisTaskRecursiveSoftDrop() :
  AliAnalysisTaskEmcalJet("AliAnalysisTaskRecursiveSoftDrop", kTRUE),
  fContainer(0),
  fJetShapeSub(kNoSub),
  fPtThreshold(-9999.),
  fCentSelectOn(kTRUE),
  fCentMin(0),
  fCentMax(10),
  fJetRadius(0.4),
  fSharedFractionPtMin(0.5),
  fReclusteringAlgo(0),
  fhJetPt(0x0),
  fhJetPhi(0x0),
  fhJetEta(0x0),
  fhDetJetPt_Matched(0x0),
  fTreeRecursive_Det(0),
  fTreeRecursive_True(0),
  fAddMedScat(kFALSE),
  fAddMedScatPtFrac(1),
  fAddMedScatN(100),
  fDoSubJetAreaSub(kFALSE)

{
  for(Int_t i=0;i<5;i++){
    fShapesVar_Det[i]=0;
    fShapesVar_True[i]=0;
  }
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskRecursiveSoftDrop::AliAnalysisTaskRecursiveSoftDrop(const char *name) :
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fContainer(0),
  fJetShapeSub(kNoSub),
  fPtThreshold(-9999.),
  fCentSelectOn(kTRUE),
  fCentMin(0),
  fCentMax(10),
  fJetRadius(0.4),
  fSharedFractionPtMin(0.5),
  fReclusteringAlgo(0),
  fhJetPt(0x0),
  fhJetPhi(0x0),
  fhJetEta(0x0),
  fhDetJetPt_Matched(0x0),
  fTreeRecursive_Det(0),
  fTreeRecursive_True(0),
  fAddMedScat(kFALSE),
  fAddMedScatPtFrac(1),
  fAddMedScatN(100),
  fDoSubJetAreaSub(kFALSE)

{
  // Standard constructor.
  for(Int_t i=0;i<5;i++){
    fShapesVar_Det[i]=0;
    fShapesVar_True[i]=0;
  }
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskRecursiveSoftDrop::~AliAnalysisTaskRecursiveSoftDrop()
{
  // Destructor.
}

//________________________________________________________________________
 void AliAnalysisTaskRecursiveSoftDrop::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  TH1::AddDirectory(oldStatus);
  //create a tree used for the MC data and making a 4D response matrix
  const char* nameoutput = GetOutputSlot(2)->GetContainer()->GetName();
  fTreeRecursive_Det = new TTree(nameoutput, nameoutput);
  const char* nameoutput2 = GetOutputSlot(3)->GetContainer()->GetName();
  fTreeRecursive_True = new TTree(nameoutput2, nameoutput2);

  const Int_t intBranches = 5;

  std::vector<TString> fShapesVarNames_Det(intBranches), fShapesVarNames_True(intBranches);


  fShapesVarNames_Det[0] = "Pt";
  fShapesVarNames_Det[1] = "Z";
  fShapesVarNames_Det[2] = "Theta";
  fShapesVarNames_Det[3] = "N";
  fShapesVarNames_Det[4] = "ParentPt";
  fShapesVarNames_True[0] = "Pt_Truth";
  fShapesVarNames_True[1] = "Z_Truth";
  fShapesVarNames_True[2] = "Theta_Truth";
  fShapesVarNames_True[3] = "N_Truth";
  fShapesVarNames_True[4] = "ParentPt_Truth";

  for(Int_t ivar=0; ivar < intBranches; ivar++){
    cout<<"looping over variables"<<endl;
    fTreeRecursive_Det->Branch(fShapesVarNames_Det[ivar].Data(), &fShapesVar_Det[ivar], Form("%s/D", fShapesVarNames_Det[ivar].Data()));
    fTreeRecursive_True->Branch(fShapesVarNames_True[ivar].Data(), &fShapesVar_True[ivar], Form("%s/D", fShapesVarNames_True[ivar].Data()));
  }

  fhJetPt= new TH1F("fhJetPt", "Jet Pt",1500,-0.5,149.5 );
  fOutput->Add(fhJetPt);
  fhJetPhi= new TH1F("fhJetPhi", "Jet Phi",360 , -1.5*(TMath::Pi()), 1.5*(TMath::Pi()));
  fOutput->Add(fhJetPhi);
  fhJetEta= new TH1F("fhJetEta", "Jet Eta",100,-2,2);
  fOutput->Add(fhJetEta);
  fhDetJetPt_Matched= new TH1F("fhDetJetPt_Matched", "Jet Pt",200,-0.5,199.5 );
  fOutput->Add(fhDetJetPt_Matched);

  PostData(1,fOutput);
  PostData(2,fTreeRecursive_Det);
  PostData(3,fTreeRecursive_True);
  // delete [] fShapesVarNames;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskRecursiveSoftDrop::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().


  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskRecursiveSoftDrop::FillHistograms()
{
  if (fCentSelectOn){
    if ((fCent>fCentMax) || (fCent<fCentMin)) return 0;
  }

  if(fJetType == kData){
    AliEmcalJet *Jet1 = NULL; //Original Jet in the event
    AliJetContainer *JetCont= GetJetContainer(0); //Jet Container for event
    Double_t JetPhi=0;
    Double_t JetPt_ForThreshold=0;
    if(JetCont) {
      JetCont->ResetCurrentID();
      while((Jet1=JetCont->GetNextAcceptJet())) {
	if(!Jet1) continue;
	if (fJetShapeSub==kNoSub) JetPt_ForThreshold = Jet1->Pt()-(GetRhoVal(0)*Jet1->Area());
	else JetPt_ForThreshold = Jet1->Pt();
	if(JetPt_ForThreshold<fPtThreshold) {
	  continue;
	}
	else {
	  fhJetPt->Fill(Jet1->Pt());
	  JetPhi=Jet1->Phi();
	  if(JetPhi < -1*TMath::Pi()) JetPhi += (2*TMath::Pi());
	  else if (JetPhi > TMath::Pi()) JetPhi -= (2*TMath::Pi());
	  fhJetPhi->Fill(JetPhi);
	  fhJetEta->Fill(Jet1->Eta());
	  RecursiveParents(Jet1,JetCont,kFALSE);
	}
      }
    }
  }

  if(fJetType == kEmb){
    AliEmcalJet *JetHybridS = NULL; //Subtracted hybrid Jet
    AliEmcalJet *JetHybridUS = NULL; //Unsubtracted Hybrid Jet     //For matching SubtractedHybrid->DetPythia this jet container is also Subtracted Hybrid
    AliEmcalJet *JetPythDet = NULL; //Detector Level Pythia Jet
    AliEmcalJet *JetPythTrue = NULL; //Particle Level Pyhtia Jet
    AliJetContainer *JetContHybridS= GetJetContainer(0); //Jet Container for Subtracted Hybrid Jets
    AliJetContainer *JetContHybridUS= GetJetContainer(1); //Jet Container for Unsubtracted Hybrid Jets
    AliJetContainer *JetContPythDet= GetJetContainer(2); //Jet Container for Detector Level Pyhtia Jets
    AliJetContainer *JetContPythTrue= GetJetContainer(3); //Jet Container for Particle Level Pythia Jets



    Bool_t JetsMatched = kFALSE;
    Double_t JetPtThreshold;
    JetContHybridS->ResetCurrentID();
    JetContHybridUS->ResetCurrentID();
    JetContPythDet->ResetCurrentID();
    JetContPythTrue->ResetCurrentID();

    while((JetHybridS = JetContHybridS->GetNextAcceptJet())){ //Get next Subtracted hybrid jet
      if (fJetShapeSub==kConstSub) JetPtThreshold=JetHybridS->Pt();
      else JetPtThreshold=JetHybridS->Pt()-(GetRhoVal(0)*JetHybridS->Area());
      if ( (!JetHybridS) || (JetPtThreshold<fPtThreshold)) continue; //check pT is above threshold
      Int_t JetNumber=-1;
      for(Int_t i = 0; i<JetContHybridUS->GetNJets(); i++) {
	JetHybridUS = JetContHybridUS->GetJet(i);            //Get unsubtracted jets in order
	if (!JetHybridUS) continue;

	if(JetHybridUS->GetLabel()==JetHybridS->GetLabel()) { //check if it mataches with current subtracted hybrid jet
	  JetNumber=i;
	}
      }
      if(JetNumber==-1) continue;
      JetHybridUS=JetContHybridUS->GetJet(JetNumber);        //Get the matched Unsubtracted jet
      if (JetContHybridUS->AliJetContainer::GetFractionSharedPt(JetHybridUS)<fSharedFractionPtMin) {  //Check that the US closest jet shares a minimum
	continue;                                                                                     // pT with Detector level pythia jet
      }
      JetPythDet=JetHybridUS->ClosestJet();                                                          //get the closest jet that has passed the shared pT cut
      if (!JetHybridUS) {
	Printf("Unsubtracted embedded jet does not exist, returning");
	continue;
      }
      if (!JetPythDet) continue;
      UInt_t rejectionReason = 0;
      if (!(JetContPythDet->AcceptJet(JetPythDet,rejectionReason))) continue;     //Check the detector level just is accepted
      fhDetJetPt_Matched->Fill(JetPythDet->Pt()); //Fill only matched detector level jets for tagging efficiency comparison
      JetPythTrue=JetPythDet->ClosestJet();       //Get the corresponding pythia true jet
      if(!JetPythTrue) continue;
      JetsMatched=kTRUE; //jets have been matched


      fhJetPt->Fill(JetHybridS->Pt());
      Double_t JetPhi=JetHybridS->Phi();
      if(JetPhi < -1*TMath::Pi()) JetPhi += (2*TMath::Pi());
      else if (JetPhi > TMath::Pi()) JetPhi -= (2*TMath::Pi());
      fhJetPhi->Fill(JetPhi);
      fhJetEta->Fill(JetHybridS->Eta());
      RecursiveParents(JetHybridS,JetContHybridS,kFALSE);
      RecursiveParents(JetPythTrue,JetContPythTrue,kTRUE);



     }
  }


  if(fJetType == kTrueDet){
    AliEmcalJet *JetPythDet = NULL; //Detector Level Pythia Jet
    AliEmcalJet *JetPythTrue = NULL; //Particle Level Pyhtia Jet
    AliJetContainer *JetContPythDet= GetJetContainer(0); //Jet Container for Detector Level Pyhtia Jets
    AliJetContainer *JetContPythTrue= GetJetContainer(1); //Jet Container for Particle Level Pythia Jets



    Bool_t JetsMatched = kFALSE;
    Double_t JetPtThreshold;
    JetContPythDet->ResetCurrentID();
    JetContPythTrue->ResetCurrentID();

    while((JetPythDet = JetContPythDet->GetNextAcceptJet())){ //Get next detector level jet
      if (fJetShapeSub==kConstSub) JetPtThreshold=JetPythDet->Pt();
      else JetPtThreshold=JetPythDet->Pt()-(GetRhoVal(0)*JetPythDet->Area());
      if ( (!JetPythDet) || (JetPtThreshold<fPtThreshold)) continue; //check pT is above threshold
      Int_t JetNumber=-1;
      if((JetPythTrue = JetPythDet->ClosestJet())){
	JetsMatched=kTRUE;
      }
      else continue;


      fhJetPt->Fill(JetPythDet->Pt());
      Double_t JetPhi=JetPythDet->Phi();
      if(JetPhi < -1*TMath::Pi()) JetPhi += (2*TMath::Pi());
      else if (JetPhi > TMath::Pi()) JetPhi -= (2*TMath::Pi());
      fhJetPhi->Fill(JetPhi);
      fhJetEta->Fill(JetPythDet->Eta());
      RecursiveParents(JetPythDet,JetContPythDet,kFALSE);
      RecursiveParents(JetPythTrue,JetContPythTrue,kTRUE);



     }
  }


  return kTRUE;
}

//_________________________________________________________________________
void AliAnalysisTaskRecursiveSoftDrop::RecursiveParents(AliEmcalJet *fJet,AliJetContainer *fJetCont,Bool_t bTruth){
  std::vector<fastjet::PseudoJet>  fInputVectors;
  fInputVectors.clear();
  fastjet::PseudoJet  PseudoTracks;
  double xflagalgo=0;
  AliParticleContainer *fTrackCont = fJetCont->GetParticleContainer();

  if (fTrackCont) for (Int_t i=0; i<fJet->GetNumberOfTracks(); i++) {
      AliVParticle *fTrk = fJet->TrackAt(i, fTrackCont->GetArray());
      if (!fTrk) continue;
      PseudoTracks.reset(fTrk->Px(), fTrk->Py(), fTrk->Pz(),fTrk->E());
      PseudoTracks.set_user_index(fJet->TrackAt(i)+100);
      fInputVectors.push_back(PseudoTracks);

    }
  if(fAddMedScat){
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
  }




  fastjet::JetAlgorithm jetalgo(fastjet::antikt_algorithm);
  if(fReclusteringAlgo==0){ xflagalgo=0.5;
    jetalgo=fastjet::kt_algorithm ;}

  if(fReclusteringAlgo==1){ xflagalgo=1.5;
    jetalgo=fastjet::cambridge_algorithm;
  }
  if(fReclusteringAlgo==2){ xflagalgo=2.5;
    jetalgo=fastjet::antikt_algorithm;
  }

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
    Int_t n = 0;
    Double_t jet_pT=jj.perp();
    while(jj.has_parents(j1,j2)){
      n++;
      double area1 = j1.area();
      double area2 = j2.area();
      if((j1.perp()-area1*GetRhoVal(0)) < (j2.perp()-area2*GetRhoVal(0))) swap(j1,j2);
      area1 = j1.area();
      area2 = j2.area();
      double delta_R=j1.delta_R(j2);
      double z = 0;
      if(fJetShapeSub==kNoSub && fDoSubJetAreaSub == kTRUE) z = (j2.perp()-area2*GetRhoVal(0))/((j1.perp()-area1*GetRhoVal(0))+(j2.perp()-area2*GetRhoVal(0)));
      else z=j2.perp()/(j1.perp()+j2.perp());
      if(bTruth) {
	fShapesVar_True[0]=jet_pT;
	fShapesVar_True[1]=z;
	fShapesVar_True[2]=delta_R;
	fShapesVar_True[3]=n;
	fShapesVar_True[4]=jj.perp();
	fTreeRecursive_True->Fill();
      }
      else {
	fShapesVar_Det[0]=jet_pT;
	fShapesVar_Det[1]=z;
	fShapesVar_Det[2]=delta_R;
	fShapesVar_Det[3]=n;
	fShapesVar_Det[4]=jj.perp();
	fTreeRecursive_Det->Fill();
      }
      jj=j1;
    }

  } catch (fastjet::Error) {
    AliError(" [w] FJ Exception caught.");
    //return -1;
  }
  return;
}


//________________________________________________________________________
Bool_t AliAnalysisTaskRecursiveSoftDrop::RetrieveEventObjects() {
  //
  // retrieve event objects
  //
  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  return kTRUE;
}


//_______________________________________________________________________
void AliAnalysisTaskRecursiveSoftDrop::Terminate(Option_t *)
{
  // Called once at the end of the analysis.

}


AliAnalysisTaskRecursiveSoftDrop* AliAnalysisTaskRecursiveSoftDrop::AddTaskRecursiveSoftDrop(
											     const char * njetsData, //data jets
											     const char * njetsTrue, //Pyhthia Particle Level
											     const char * njetsDet,
											     const char * njetsHybridUs,
											     const char * njetsHybridS,
											     const Double_t R,
											     const char * nrhoBase,
											     const char * ntracksData,
											     const char * ntracksTrue,
											     const char * ntracksDet,
											     const char * ntracksHybridUs,
											     const char * ntracksHybridS,
											     const char *type,
											     const char *CentEst,
											     Int_t       pSel,
											     TString     trigClass,
											     TString     kEmcalTriggers,
											     TString     tag,
											     AliAnalysisTaskRecursiveSoftDrop::JetShapeSub jetShapeSub,
											     AliAnalysisTaskRecursiveSoftDrop::JetType fjetType
											     ) {



  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    {
      ::Error("AddTaskRecursiveSoftDrop","No analysis manager found.");
      return 0;
    }
  Bool_t ismc=kFALSE;
  ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
    {
      ::Error("AliAnalysisTaskRecursiveSoftDrop", "This task requires an input event handler");
      return NULL;
    }
  TString wagonName1;
  TString wagonName2;
  TString wagonName3;

  if(fjetType==AliAnalysisTaskRecursiveSoftDrop::kData){
    wagonName1 = Form("AliAnalysisTaskRecursiveSoftDrop_%s_TC%s%s",njetsData,trigClass.Data(),tag.Data());
    wagonName2 = Form("AliAnalysisTaskRecursiveSoftDrop_%s_TC%s%sTree_Det",njetsData,trigClass.Data(),tag.Data());
    wagonName3 = Form("AliAnalysisTaskRecursiveSoftDrop_%s_TC%s%sTree_True",njetsData,trigClass.Data(),tag.Data());
  }
  if(fjetType==AliAnalysisTaskRecursiveSoftDrop::kEmb){
    wagonName1 = Form("AliAnalysisTaskRecursiveSoftDrop_%s_TC%s%s",njetsHybridS,trigClass.Data(),tag.Data());
    wagonName2 = Form("AliAnalysisTaskRecursiveSoftDrop_%s_TC%s%sTree_Det",njetsHybridS,trigClass.Data(),tag.Data());
    wagonName3 = Form("AliAnalysisTaskRecursiveSoftDrop_%s_TC%s%sTree_True",njetsHybridS,trigClass.Data(),tag.Data());
  }
  if(fjetType==AliAnalysisTaskRecursiveSoftDrop::kTrueDet){
    wagonName1 = Form("AliAnalysisTaskRecursiveSoftDrop_%s_TC%s%s",njetsDet,trigClass.Data(),tag.Data());
    wagonName2 = Form("AliAnalysisTaskRecursiveSoftDrop_%s_TC%s%sTree_Det",njetsDet,trigClass.Data(),tag.Data());
    wagonName3 = Form("AliAnalysisTaskRecursiveSoftDrop_%s_TC%s%sTree_True",njetsDet,trigClass.Data(),tag.Data());
  }
  //Configure jet tagger task
  AliAnalysisTaskRecursiveSoftDrop *task = new AliAnalysisTaskRecursiveSoftDrop(wagonName1.Data());

  task->SetJetShapeSub(jetShapeSub);
  task->SetJetType(fjetType);
  task->SetJetRadius(R);

  AliParticleContainer *trackContData=0x0;
  AliParticleContainer *trackContDet=0x0;
  AliParticleContainer *trackContTrue=0x0;
  AliParticleContainer *trackContHybridS=0x0;
  AliParticleContainer *trackContHybridUs=0x0;

  if(fjetType!=AliAnalysisTaskRecursiveSoftDrop::kTrueDet) trackContData = task->AddParticleContainer(ntracksData);
  if(fjetType==AliAnalysisTaskRecursiveSoftDrop::kEmb){
    trackContDet = task->AddParticleContainer(ntracksDet);
    trackContTrue = task->AddMCParticleContainer(ntracksTrue);
    trackContHybridS = task->AddParticleContainer(ntracksHybridS);
    trackContHybridUs = task->AddParticleContainer(ntracksHybridUs);
  }
  if(fjetType==AliAnalysisTaskRecursiveSoftDrop::kTrueDet){
    trackContDet = task->AddParticleContainer(ntracksDet);
    trackContTrue = task->AddMCParticleContainer(ntracksTrue);
  }
  AliJetContainer *JetContData=0x0;
  AliJetContainer *JetContTrue=0x0;
  AliJetContainer *JetContDet=0x0;
  AliJetContainer *JetContHybridUs=0x0;
  AliJetContainer *JetContHybridS=0x0;
  TString strType(type);

  if(fjetType==AliAnalysisTaskRecursiveSoftDrop::kData){
    JetContData = task->AddJetContainer(njetsData,strType,R); //Data
    if(JetContData) {
      JetContData->SetRhoName(nrhoBase);
      JetContData->ConnectParticleContainer(trackContData);
      JetContData->SetPercAreaCut(0.6);
      JetContData->SetJetRadius(R);
      JetContData->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
      if(jetShapeSub==AliAnalysisTaskRecursiveSoftDrop::kConstSub) JetContData->SetAreaEmcCut(-2);
    }
  }
  if(fjetType==AliAnalysisTaskRecursiveSoftDrop::kEmb){
    JetContHybridS = task->AddJetContainer(njetsHybridS,strType,R); //HybridS
    if(JetContHybridS) {
      JetContHybridS->SetRhoName(nrhoBase);
      JetContHybridS->ConnectParticleContainer(trackContHybridS);
      JetContHybridS->SetPercAreaCut(0.6);
      JetContHybridS->SetJetRadius(R);
      JetContHybridS->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
      if(jetShapeSub==AliAnalysisTaskRecursiveSoftDrop::kConstSub) JetContHybridS->SetAreaEmcCut(-2);
    }
    JetContHybridUs = task->AddJetContainer(njetsHybridUs,strType,R); //HybridUs
    if(JetContHybridUs) {
      JetContHybridUs->SetRhoName(nrhoBase);
      JetContHybridUs->ConnectParticleContainer(trackContHybridUs);
      JetContHybridUs->SetPercAreaCut(0.6);
      JetContHybridUs->SetJetRadius(R);
      JetContHybridUs->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
      if(jetShapeSub==AliAnalysisTaskRecursiveSoftDrop::kConstSub) JetContHybridUs->SetAreaEmcCut(-2);
    }
    JetContDet = task->AddJetContainer(njetsDet,strType,R); //Det
    if(JetContDet) {
      JetContDet->SetRhoName(nrhoBase);
      JetContDet->ConnectParticleContainer(trackContDet);
      JetContDet->SetPercAreaCut(0.6);
      JetContDet->SetJetRadius(R);
      JetContDet->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
      if(jetShapeSub==AliAnalysisTaskRecursiveSoftDrop::kConstSub) JetContDet->SetAreaEmcCut(-2);
    }
    JetContTrue = task->AddJetContainer(njetsTrue,strType,R); //True
    if(JetContTrue) {
      JetContTrue->SetRhoName(nrhoBase);
      JetContTrue->ConnectParticleContainer(trackContTrue);
      JetContTrue->SetPercAreaCut(0.6);
      JetContTrue->SetJetRadius(R);
      JetContTrue->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
      if(jetShapeSub==AliAnalysisTaskRecursiveSoftDrop::kConstSub) JetContTrue->SetAreaEmcCut(-2);
    }
  }

  if(fjetType==AliAnalysisTaskRecursiveSoftDrop::kTrueDet){
    JetContDet = task->AddJetContainer(njetsDet,strType,R); //Det
    if(JetContDet) {
      JetContDet->SetRhoName(nrhoBase);
      JetContDet->ConnectParticleContainer(trackContDet);
      JetContDet->SetPercAreaCut(0.6);
      JetContDet->SetJetRadius(R);
      JetContDet->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
      if(jetShapeSub==AliAnalysisTaskRecursiveSoftDrop::kConstSub) JetContDet->SetAreaEmcCut(-2);
    }
    JetContTrue = task->AddJetContainer(njetsTrue,strType,R); //True
    if(JetContTrue) {
      JetContTrue->SetRhoName(nrhoBase);
      JetContTrue->ConnectParticleContainer(trackContTrue);
      JetContTrue->SetPercAreaCut(0.6);
      JetContTrue->SetJetRadius(R);
      JetContTrue->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
      if(jetShapeSub==AliAnalysisTaskRecursiveSoftDrop::kConstSub) JetContTrue->SetAreaEmcCut(-2);
    }
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
  TString contName3(wagonName3);

  if (fjetType == AliAnalysisTaskRecursiveSoftDrop::kEmb){
    contName1 += "_Embedded";
    contName2 += "_Embedded";
    contName2 += "_Embedded";

  }


  TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName1.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,1,coutput1);
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(contName2.Data(), TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,2,coutput2);
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(contName3.Data(), TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,3,coutput3);

  return task;

}
