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
  fhJetPt(0x0),
  fhJetPhi(0x0),
  fhJetEta(0x0),
  fTreeRecursive_Det(0),
  fTreeRecursive_True(0)

{
  for(Int_t i=0;i<4;i++){
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
  fhJetPt(0x0),
  fhJetPhi(0x0),
  fhJetEta(0x0),
  fTreeRecursive_Det(0),
  fTreeRecursive_True(0)
{
  // Standard constructor.
  for(Int_t i=0;i<4;i++){
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

  const Int_t intBranches = 4;

  std::vector<TString> fShapesVarNames_Det(intBranches), fShapesVarNames_True(intBranches);

  
  fShapesVarNames_Det[0] = "Pt";
  fShapesVarNames_Det[1] = "Z";
  fShapesVarNames_Det[2] = "Theta";
  fShapesVarNames_Det[3] = "N";
  fShapesVarNames_True[0] = "Pt_Truth";
  fShapesVarNames_True[1] = "Z_Truth";
  fShapesVarNames_True[2] = "Theta_Truth";
  fShapesVarNames_True[3] = "N_Truth";

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
	RecursiveParents(Jet1,JetCont,1,kFALSE); //Third argument = reclustering algorithm (0=Antikt,1=CA,2=kt)
      }
    }
  }

  return kTRUE;
}

//_________________________________________________________________________
void AliAnalysisTaskRecursiveSoftDrop::RecursiveParents(AliEmcalJet *fJet,AliJetContainer *fJetCont, Int_t ReclusterAlgo,Bool_t bTruth){
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



  fastjet::JetAlgorithm jetalgo(fastjet::antikt_algorithm);
  if(ReclusterAlgo==0){ xflagalgo=0.5;
    jetalgo=fastjet::kt_algorithm ;}
      
  if(ReclusterAlgo==1){ xflagalgo=1.5;
    jetalgo=fastjet::cambridge_algorithm;
  }
  if(ReclusterAlgo==2){ xflagalgo=2.5;
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
      if(j1.perp() < j2.perp()) swap(j1,j2);
      double delta_R=j1.delta_R(j2);
      double z=j2.perp()/(j1.perp()+j2.perp());
      if(bTruth) {
	fShapesVar_True[0]=jet_pT;
	fShapesVar_True[1]=z;
	fShapesVar_True[2]=delta_R;
	fShapesVar_True[3]=n;
	fTreeRecursive_True->Fill();
      }
      else {
	fShapesVar_Det[0]=jet_pT;
	fShapesVar_Det[1]=z;
	fShapesVar_Det[2]=delta_R;
	fShapesVar_Det[3]=n;
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
