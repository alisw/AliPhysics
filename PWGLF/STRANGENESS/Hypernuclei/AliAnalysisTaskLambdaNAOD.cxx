/**************************************************************************
 * Author : Nicole Alice Martin (nicole.alice.martin@cern.ch)                  *
 *                                                                        *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//-----------------------------------------------------------------
//               AliAnalysisTaskLambdaNAOD  class
//  task for the investigation of (anti-)lambda-n bound state
//          uses the V0 finders, based on AODs or ESDS
//-----------------------------------------------------------------


class TTree;
class TParticle;
class TVector3;

class AliESDVertex;
class AliESDv0;

class AliAODVertex;
class AliAODv0;

#include "TChain.h"
#include <Riostream.h>
#include <iostream>
#include <fstream>
#include <string>
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TFile.h"
#include <TRandom3.h>
#include "TPDGCode.h"
#include <TDatabasePDG.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliCentrality.h"
#include "AliV0vertexer.h"
#include "AliPIDResponse.h"
#include "AliMultiplicity.h"
#include "AliVertexerTracks.h"

#include "AliVEvent.h"
#include "AliVTrack.h"

#include "AliESDInputHandler.h" 
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliESDpid.h"
#include "AliESDv0.h"

#include "AliAODInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"

#include "AliCFContainer.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"

#include "AliAnalysisTaskLambdaNAOD.h"
using namespace std;


ClassImp(AliAnalysisTaskLambdaNAOD)


//________________________________________________________________________
AliAnalysisTaskLambdaNAOD::AliAnalysisTaskLambdaNAOD()
: AliAnalysisTaskSE(), 
  fAnalysisType("ESD"),
  fEventHandler(0),
  fESDtrackCutsV0(0),
  fESDpid(0),
  fPIDResponse(0),
  fTreeV0(0),
  fHistNumberOfEvents(0),
  fHistTriggerStat(0),
  fHistLambdaNeutronPtGen(0),
  fHistAntiLambdaNeutronPtGen(0),
  fHistLambdaNeutronInvaMassGen(0),
  fHistAntiLambdaNeutronInvaMassGen(0),
  fHistLambdaNeutronDecayLengthGen(0),
  fHistAntiLambdaNeutronDecayLengthGen(0),
  fHistLambdaNeutronPtAso(0),
  fHistLambdaNeutronPtAsoCuts(0),
  fHistAntiLambdaNeutronPtAso(0),
  fHistAntiLambdaNeutronPtAsoCuts(0),
  fHistLambdaNeutronInvaMassAso(0),
  fHistAntiLambdaNeutronInvaMassAso(0), 
  fHistLambdaNeutronDecayLengthAso(0),
  fHistAntiLambdaNeutronDecayLengthAso(0),
  fof(0),
  fHistArmenterosPodolanskiDeuteronPion(0),
  fHistArmenterosPodolanskiAntiDeuteronPion(0),
  fHistDeDxQA(0),
//fHistdEdxMC(0),
  fNTriggers(5),
  fMCtrue(0),
  fOnlyQA(0),
  fTriggerFired(),
  fV0object(NULL),
  fItrk(0),
  fOutputContainer(NULL)
{
  // default Constructor

  // Define input and output slots here
}

//________________________________________________________________________
AliAnalysisTaskLambdaNAOD::AliAnalysisTaskLambdaNAOD(const char *name)
  : AliAnalysisTaskSE(name), 
    fAnalysisType("ESD"),
    fEventHandler(0),
    fESDtrackCutsV0(0),
    fESDpid(0),
    fPIDResponse(0),
    fTreeV0(0),
    fHistNumberOfEvents(0),
    fHistTriggerStat(0),
    fHistLambdaNeutronPtGen(0),
    fHistAntiLambdaNeutronPtGen(0),
    fHistLambdaNeutronInvaMassGen(0),
    fHistAntiLambdaNeutronInvaMassGen(0),
    fHistLambdaNeutronDecayLengthGen(0),
    fHistAntiLambdaNeutronDecayLengthGen(0),
    fHistLambdaNeutronPtAso(0),
    fHistLambdaNeutronPtAsoCuts(0),
    fHistAntiLambdaNeutronPtAso(0),
    fHistAntiLambdaNeutronPtAsoCuts(0),
    fHistLambdaNeutronInvaMassAso(0),
    fHistAntiLambdaNeutronInvaMassAso(0),
    fHistLambdaNeutronDecayLengthAso(0),
    fHistAntiLambdaNeutronDecayLengthAso(0),
    fof(0),
    fHistArmenterosPodolanskiDeuteronPion(0),
    fHistArmenterosPodolanskiAntiDeuteronPion(0),
    fHistDeDxQA(0),
    //fHistdEdxMC(0),
    fNTriggers(5),
    fMCtrue(0),
    fOnlyQA(0),
    fTriggerFired(),
    fV0object(NULL),
    fItrk(0),
    fOutputContainer(NULL)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(1, TObjArray::Class());
  DefineOutput(2, TTree::Class());

  //ESD Track Cuts for v0's
  fESDtrackCutsV0 = new AliESDtrackCuts("AliESDtrackCutsV0","AliESDtrackCutsV0");
  fESDtrackCutsV0->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCutsV0->SetMinNClustersTPC(60);
  fESDtrackCutsV0->SetMaxChi2PerClusterTPC(5);
  fESDtrackCutsV0->SetRequireTPCRefit(kTRUE);
  //  fESDtrackCutsV0->SetMinNClustersITS(1); // to be tested if this cut is not too strong
  fESDtrackCutsV0->SetEtaRange(-0.9,0.9);

  
  fMCtrue = kTRUE; 
  fOnlyQA = kTRUE;

}

//____________________________________________________________
AliAnalysisTaskLambdaNAOD::~AliAnalysisTaskLambdaNAOD(){
  //
  // Destructor
  //
  //if (fESD) delete fESD;
  if (fESDtrackCutsV0) delete fESDtrackCutsV0;

}

//____________________________________________________________
const Int_t AliAnalysisTaskLambdaNAOD::fgkPdgCode[] = {
  211,                //PionPlus
  -211,               //PionMinus
  2212,               //Proton
  -2212,              //Anti-Proton
  1000010020,         //Deuteron
  -1000010020,        //Anti-Deuteron
  1000020030,         //Helium3
  -1000020030,        //Anti-Helium3
  3122,               //Lambda
  -3122,              //Anti-Lambda
  1010000020,         //LambdaNeutron
  -1010000020,        //Anti-Lambda-Neutron
  1010010030,         //Hypertriton
  -1010010030         //Anti-Hypertriton
};

//____________________________________________________________
const Double_t AliAnalysisTaskLambdaNAOD::fgkMass[] = {
  0.13957,           //Pion
  0.93827,           //Proton
  1.875612,          //Deuteron
  2.808921,          //Triton
  2.80941            //Helium3 Quelle: Wolfram Alpha
};

//________________________________________________________________________
void AliAnalysisTaskLambdaNAOD::UserCreateOutputObjects(){
  
  // New PID object
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
  
  // Create histograms for MC
  //generated histogramms
  fHistLambdaNeutronPtGen = new TH1F("fHistLambdaNeutronPtGen", "Generated  #Lambdan", 201, 0., 10.1);
  fHistLambdaNeutronPtGen->GetYaxis()->SetTitle("Counts");
  fHistLambdaNeutronPtGen->GetXaxis()->SetTitle("#it{p}_{T}  (GeV/#it{c})");

  fHistAntiLambdaNeutronPtGen = new TH1F("fHistAntiLambdaNeutronPtGen", "Generated  #bar{#Lambdan}", 201, 0., 10.1);
  fHistAntiLambdaNeutronPtGen->GetYaxis()->SetTitle("Counts");
  fHistAntiLambdaNeutronPtGen->GetXaxis()->SetTitle("#it{p}_{T}   (GeV/#it{c})");

  fHistLambdaNeutronInvaMassGen = new TH1F("fHistLambdaNeutronInvaMassGen", "Generated #Lambdan ", 100, 2.0, 2.1); 
  fHistLambdaNeutronInvaMassGen->GetYaxis()->SetTitle("Counts");
  fHistLambdaNeutronInvaMassGen->GetXaxis()->SetTitle("Invariant mass (d #pi^{-}) (GeV/#it{c}^{2})");

  fHistAntiLambdaNeutronInvaMassGen = new TH1F("fHistAntiLambdaNeutronInvaMassGen", "Generated  #bar{#Lambdan}", 100, 2.0, 2.1); 
  fHistAntiLambdaNeutronInvaMassGen->GetYaxis()->SetTitle("Counts");
  fHistAntiLambdaNeutronInvaMassGen->GetXaxis()->SetTitle("Invariant mass (#bar{d} #pi^{+}) (GeV/#it{c}^{2})");

  fHistLambdaNeutronDecayLengthGen = new TH1F("fHistLambdaNeutronDecayLengthGen", "Generated  #Lambdan", 401, 0., 400.1);
  fHistLambdaNeutronDecayLengthGen->GetYaxis()->SetTitle("Counts");
  fHistLambdaNeutronDecayLengthGen->GetXaxis()->SetTitle("#it{decay length}  (cm)");

  fHistAntiLambdaNeutronDecayLengthGen = new TH1F("fHistAntiLambdaNeutronDecayLengthGen", "Generated #bar{#Lambdan}", 401, 0., 400.1);
  fHistAntiLambdaNeutronDecayLengthGen->GetYaxis()->SetTitle("Counts");
  fHistAntiLambdaNeutronDecayLengthGen->GetXaxis()->SetTitle("#it{decay length}  (cm)");

  //assoziated (reconstracted) histogramms
  fHistLambdaNeutronPtAso = new TH1F("fHistLambdaNeutronPtAso", "Associated   #Lambdan", 201, 0., 10.1);
  fHistLambdaNeutronPtAso->GetYaxis()->SetTitle("Counts");
  fHistLambdaNeutronPtAso->GetXaxis()->SetTitle("#it{p}_{T}  (GeV/#it{c})");

  fHistLambdaNeutronPtAsoCuts = new TH1F("fHistLambdaNeutronPtAsoCuts", "Associated   #Lambdan Cuts", 201, 0., 10.1);
  fHistLambdaNeutronPtAsoCuts->GetYaxis()->SetTitle("Counts");
  fHistLambdaNeutronPtAsoCuts->GetXaxis()->SetTitle("#it{p}_{T}  (GeV/#it{c})");

  fHistAntiLambdaNeutronPtAsoCuts = new TH1F("fHistAntiLambdaNeutronPtAsoCuts", " associated  #bar{#Lambdan} Cuts", 201, 0., 10.1);
  fHistAntiLambdaNeutronPtAsoCuts->GetYaxis()->SetTitle("Counts");
  fHistAntiLambdaNeutronPtAsoCuts->GetXaxis()->SetTitle("#it{p}_{T}  (GeV/#it{c})");

  fHistAntiLambdaNeutronPtAso = new TH1F("fHistAntiLambdaNeutronPtAso", " associated  #bar{#Lambdan}", 201, 0., 10.1);
  fHistAntiLambdaNeutronPtAso->GetYaxis()->SetTitle("Counts");
  fHistAntiLambdaNeutronPtAso->GetXaxis()->SetTitle("#it{p}_{T}  (GeV/#it{c})");

  fHistLambdaNeutronInvaMassAso = new TH1F("fHistLambdaNeutronInvaMassAso", "Associated  #Lambdan", 100, 2.0, 2.1); 
  fHistLambdaNeutronInvaMassAso->GetYaxis()->SetTitle("Counts");
  fHistLambdaNeutronInvaMassAso->GetXaxis()->SetTitle("Invariant mass (d #pi^{-}) (GeV/#it{c}^{2})");

  fHistAntiLambdaNeutronInvaMassAso = new TH1F("fHistAntiLambdaNeutronInvaMassAso", " Associated  #bar{#Lambdan}", 100, 2.0, 2.1); 
  fHistAntiLambdaNeutronInvaMassAso->GetYaxis()->SetTitle("Counts");
  fHistAntiLambdaNeutronInvaMassAso->GetXaxis()->SetTitle("Invariant mass (#bar{d} #pi^{+})  (GeV/#it{c}^{2})");

  fHistLambdaNeutronDecayLengthAso = new TH1F("fHistLambdaNeutronDecayLengthAso", "Associated  #Lambdan", 401, 0., 400.1);
  fHistLambdaNeutronDecayLengthAso->GetYaxis()->SetTitle("Counts");
  fHistLambdaNeutronDecayLengthAso->GetXaxis()->SetTitle("#it{decay length}  (cm)");

  fHistAntiLambdaNeutronDecayLengthAso = new TH1F("fHistAntiLambdaNeutronDecayLengthAso", "Associated #bar{#Lambdan}", 401, 0., 400.1);
  fHistAntiLambdaNeutronDecayLengthAso->GetYaxis()->SetTitle("Counts");
  fHistAntiLambdaNeutronDecayLengthAso->GetXaxis()->SetTitle("#it{decay length}  (cm)");

  fHistLambdaNeutronPtAso = new TH1F("fHistLambdaNeutronPtAso", "Associated   #Lambdan", 201, 0., 10.1);
  fHistLambdaNeutronPtAso->GetYaxis()->SetTitle("Counts");
  fHistLambdaNeutronPtAso->GetXaxis()->SetTitle("#it{p}_{T}  (GeV/#it{c})");

  //Tree
  //------------ Tree and branch definitions ----------------//
  fTreeV0 = new TTree("tree", "fTreeV0");  
  fTreeV0->Branch("fItrk", &fItrk, "fItrk/I");

  fTreeV0->Branch("fV0object",&fV0object,"fV0object[fItrk]");
  fTreeV0->Branch("fV0finder",fV0finder,"fV0finder[fItrk]/I");
  fTreeV0->Branch("fkMB",fkMB,"fkMB[fItrk]/I");
  fTreeV0->Branch("fkCentral",fkCentral,"fkCentral[fItrk]/I");
  fTreeV0->Branch("fkSemiCentral",fkSemiCentral,"fkSemiCentral[fItrk]/I");
  fTreeV0->Branch("fkEMCEJE",fkEMCEJE,"fkEMCEJE[fItrk]/I");
  fTreeV0->Branch("fkEMCEGA",fkEMCEGA,"fkEMCEGA[fItrk]/I");
 
  fTreeV0->Branch("fPtotN",fPtotN,"fPtotN[fItrk]/D");
  fTreeV0->Branch("fPtotP",fPtotP,"fPtotP[fItrk]/D");
  fTreeV0->Branch("fMotherPt",fMotherPt,"fMotherPt[fItrk]/D");
  fTreeV0->Branch("fdEdxN",fdEdxN,"fdEdxN[fItrk]/D");
  fTreeV0->Branch("fdEdxP",fdEdxP,"fdEdxP[fItrk]/D");
  fTreeV0->Branch("fSignN",fSignN,"fSignN[fItrk]/D");
  fTreeV0->Branch("fSignP",fSignP,"fSignP[fItrk]/D");

  fTreeV0->Branch("fDCAv0",fDCAv0,"fDCAv0[fItrk]/F"); //Dca v0 Daughters
  fTreeV0->Branch("fCosinePAv0",fCosinePAv0,"fCosinePAv0[fItrk]/F"); //Cosine of Pionting Angle
  fTreeV0->Branch("fDecayRadiusTree",fDecayRadiusTree,"fDecayRadiusTree[fItrk]/F"); //decay radius
  
  fTreeV0->Branch("fChargeComboDeuteronPionTree",fChargeComboDeuteronPionTree,"fChargeComboDeuteronPionTree[fItrk]/I"); 
  fTreeV0->Branch("fInvaMassDeuteronPionTree",fInvaMassDeuteronPionTree,"fInvaMassDeuteronPionTree[fItrk]/D"); //invariant mass

  fTreeV0->Branch("fIsCorrectlyAssociated",fIsCorrectlyAssociated,"fIsCorrectlyAssociated[fItrk]/O"); //associated hypertriton

  fTreeV0->Branch("fAmenterosAlphaTree",fAmenterosAlphaTree,"fAmenterosAlphaTree[fItrk]/D");
  fTreeV0->Branch("fAmenterosQtTree",fAmenterosQtTree,"fAmenterosQtTree[fItrk]/D");
  fTreeV0->Branch("fRotationTree",fRotationTree,"fRotationTree[fItrk]/I");
   
  //Armenteros-Podolanski
  fHistArmenterosPodolanskiDeuteronPion= new TH2F("fHistArmenterosPodolanskiDeuteronPion", "Armenteros-Podolanski d #pi^{-}", 200,-1.0,1.0, 500,0,1);
  fHistArmenterosPodolanskiDeuteronPion->GetXaxis()->SetTitle("#alpha");
  fHistArmenterosPodolanskiDeuteronPion->GetYaxis()->SetTitle("q_{t}");
  fHistArmenterosPodolanskiDeuteronPion->SetMarkerStyle(kFullCircle);

  fHistArmenterosPodolanskiAntiDeuteronPion= new TH2F("fHistArmenterosPodolanskiAntiDeuteronPion", "Armenteros-Podolanski #bar{d} #pi^{+}", 200,-1.0,1.0, 500,0,1);
  fHistArmenterosPodolanskiAntiDeuteronPion->GetXaxis()->SetTitle("#alpha");
  fHistArmenterosPodolanskiAntiDeuteronPion->GetYaxis()->SetTitle("q_{t}");
  fHistArmenterosPodolanskiAntiDeuteronPion->SetMarkerStyle(kFullCircle);

  // control histograms
  fof = new TH2F("fof", "OnTheFlyStatus ",5,0.5,5.5,2,-0.5,1.5);
  fof->GetYaxis()->SetBinLabel(1,"offline");
  fof->GetYaxis()->SetBinLabel(2,"onTheFly");
  fof->GetXaxis()->SetBinLabel(1,"total");
  fof->GetXaxis()->SetBinLabel(2,"dcaCut");
  fof->GetXaxis()->SetBinLabel(3,"cosCut");
  fof->GetXaxis()->SetBinLabel(4,"nucleonPID");
  fof->GetXaxis()->SetBinLabel(5,"pionPID");
  //fof->GetXaxis()->SetBinLabel(6,"decayRadiusCut");
  //fof->SetMarkerStyle(kFullCircle);
  
  //histogram to count number of events
  fHistNumberOfEvents = new TH1F("fHistNumberOfEvents", "Number of events", 11, -0.5, 10.5);
  fHistNumberOfEvents ->GetXaxis()->SetTitle("Centrality");
  fHistNumberOfEvents ->GetYaxis()->SetTitle("Entries");
  
  //trigger statitics histogram
  fHistTriggerStat = new TH1F("fHistTriggerStat","Trigger statistics", fNTriggers,-0.5,fNTriggers-0.5);
  
  const Char_t* aTriggerNames[] = { "kMB", "kCentral", "kSemiCentral", "kEMCEJE", "kEMCEGA" };
  for ( Int_t ii=0; ii < fNTriggers; ii++ )
    fHistTriggerStat->GetXaxis()->SetBinLabel(ii+1, aTriggerNames[ii]);

  //QA dE/dx
  fHistDeDxQA = new TH3F("fHistDeDxQA", "QA dE/dx", 400, -6.0, 6.0, 500, 0.0, 2000, 15, -0.5, 14.5);
  fHistDeDxQA->GetYaxis()->SetTitle("TPC Signal");
  fHistDeDxQA->GetXaxis()->SetTitle("p (GeV/c)");
  fHistDeDxQA->GetZaxis()->SetBinLabel(0,"all pos v0 tracks");
  fHistDeDxQA->GetZaxis()->SetBinLabel(1,"all neg v0 tracks");
  fHistDeDxQA->GetZaxis()->SetBinLabel(2,"all neg deuteron");
  fHistDeDxQA->GetZaxis()->SetBinLabel(3,"all pos deuteron");
  fHistDeDxQA->GetZaxis()->SetBinLabel(6,"all selected tracks");
  fHistDeDxQA->GetZaxis()->SetBinLabel(7,"neg deuteron for deuteron+pion");
  fHistDeDxQA->GetZaxis()->SetBinLabel(8,"pos pion for deuteron+pion");
  fHistDeDxQA->GetZaxis()->SetBinLabel(9,"pos deuteron for deuteron+pion");
  fHistDeDxQA->GetZaxis()->SetBinLabel(10,"neg pion for deuteron+pion");
  
  //Define and fill the OutputContainer
  fOutputContainer = new TObjArray(1);
  fOutputContainer->SetOwner(kTRUE);
  fOutputContainer->SetName(GetName()) ;
  fOutputContainer->Add(fHistArmenterosPodolanskiDeuteronPion);
  fOutputContainer->Add(fHistArmenterosPodolanskiAntiDeuteronPion);
  fOutputContainer->Add(fof);
  fOutputContainer->Add(fHistDeDxQA);
  fOutputContainer->Add(fHistNumberOfEvents);
  fOutputContainer->Add(fHistTriggerStat);
  fOutputContainer->Add(fHistLambdaNeutronPtGen);
  fOutputContainer->Add(fHistAntiLambdaNeutronPtGen);
  fOutputContainer->Add(fHistLambdaNeutronInvaMassGen);
  fOutputContainer->Add(fHistAntiLambdaNeutronInvaMassGen); 
  fOutputContainer->Add(fHistLambdaNeutronDecayLengthGen);
  fOutputContainer->Add(fHistAntiLambdaNeutronDecayLengthGen);
  fOutputContainer->Add(fHistLambdaNeutronPtAso);
  fOutputContainer->Add(fHistLambdaNeutronPtAsoCuts);
  fOutputContainer->Add(fHistAntiLambdaNeutronPtAso);
  fOutputContainer->Add(fHistAntiLambdaNeutronPtAsoCuts);
  fOutputContainer->Add(fHistLambdaNeutronInvaMassAso);
  fOutputContainer->Add(fHistAntiLambdaNeutronInvaMassAso);
  fOutputContainer->Add(fHistLambdaNeutronDecayLengthAso);
  fOutputContainer->Add(fHistAntiLambdaNeutronDecayLengthAso);
}

//________________________________________________________________________
void AliAnalysisTaskLambdaNAOD::UserExec(Option_t *){

  // Main loop
  // Called for each event

  //Data
  //get Event-Handler for the trigger information
  fEventHandler= dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!fEventHandler) {
    AliError("Could not get InputHandler");
    //return -1;
    return;
  }
  
  // Monte Carlo
  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!eventHandler){ 
    //Printf("ERROR: Could not retrieve MC event handler");
    fMCtrue = kFALSE;
  }

  AliMCEvent* mcEvent = 0x0;
  AliStack* stack = 0x0;
  if (eventHandler) mcEvent = eventHandler->MCEvent();
  if (!mcEvent){
    //Printf("ERROR: Could not retrieve MC event");
    if (fMCtrue) return;
  }
  
  if (fMCtrue){
    stack = mcEvent->Stack();
    if (!stack) return;
  }
  
  //look for the generated particles
  MCGenerated(stack);
  
  if (SetupEvent() < 0) {
    PostData(1, fOutputContainer);
    return;
  }
  
  //DATA
  AliESDEvent *fESDevent = 0x0;
  AliAODEvent *fAODevent = 0x0;
  
  
  if(!fMCtrue){
    if(!fPIDResponse) {
      AliError("Cannot get pid response");
      return;
    }
  }
 

  Int_t centrality = -1;
  Double_t vertex[3]          = {-100.0, -100.0, -100.0};

  //Initialisation of the event and basic event cuts:
  //1.) ESDs: 
  if (fAnalysisType == "ESD") {
        
    fESDevent=dynamic_cast<AliESDEvent*>(InputEvent());
    if (!fESDevent) {
      AliWarning("ERROR: fESDevent not available \n");
      return;
    }    
    
    //Basic event cuts: 
    //1.) vertex existence
    const AliESDVertex *vertexESD = fESDevent->GetPrimaryVertexTracks();
    if (vertexESD->GetNContributors()<1) 
      {
	// SPD vertex
	vertexESD = fESDevent->GetPrimaryVertexSPD();
	if(vertexESD->GetNContributors()<1) {
	  PostData(1, fOutputContainer);
	  return;
	}  
      }

    vertexESD->GetXYZ(vertex);

    //2. vertex position within 10 cm
    if (TMath::Abs(vertexESD->GetZv()) > 10) return;
    
    //centrality selection in PbPb
    if (fESDevent->GetEventSpecie() == 4) //species == 4 == PbPb
      { // PbPb
	AliCentrality *esdCentrality = fESDevent->GetCentrality();
	centrality = esdCentrality->GetCentralityClass10("V0M"); // centrality percentile determined with V0
	if(!fMCtrue){
	  if (centrality < 0. || centrality > 8. ) return; //select only events with centralities between 0 and 80 %
	}
      }

    //cout << "centrality "<< centrality << endl;
    // count number of events
    fHistNumberOfEvents->Fill(centrality);

 }

//2.) AODs: 
 else if (fAnalysisType == "AOD") {
   
   fAODevent = dynamic_cast<AliAODEvent*>( InputEvent() );
   if (!fAODevent) {
     AliWarning("ERROR: lAODevent not available \n");
     return;
   }
      
   const AliAODVertex *vertexAOD = fAODevent->GetPrimaryVertex();
    if (vertexAOD->GetNContributors()<1) 
      {
	// SPD vertex
	vertexAOD = fAODevent->GetPrimaryVertexSPD();
	if(vertexAOD->GetNContributors()<1) {
	  PostData(1, fOutputContainer);
	  return;
	}  
      }
    vertexAOD->GetXYZ(vertex);

    //2. vertex position within 10 cm
    if (TMath::Abs(vertex[2]) > 10) return;
   
    //centrality selection
   AliCentrality *aodCentrality = fAODevent->GetCentrality();
   centrality = aodCentrality->GetCentralityClass10("V0M"); // centrality percentile determined with V0
   if (centrality < 0 || centrality > 8) return; //select only events with centralities between 0 and 80 %
    
   // count number of events
   fHistNumberOfEvents->Fill(centrality);

 } else {
   
   Printf("Analysis type (ESD or AOD) not specified \n");
   return;
   
 }
 
  
//v0 loop


 fItrk = 0;

 Int_t runNumber = 0;
 runNumber = (InputEvent())->GetRunNumber();
 
 Int_t    nTrackMultiplicity             = -1;
 nTrackMultiplicity = (InputEvent())->GetNumberOfTracks();
  

 for (Int_t ivertex=0; ivertex<(InputEvent()->GetNumberOfV0s()); ivertex++) { //beginn v0 loop
         
   AliESDv0 * v0ESD = 0x0;
   AliAODv0 * v0AOD = 0x0;

   Bool_t v0ChargesCorrect = kTRUE;
   Bool_t testTrackCuts = kFALSE;
   //Bool_t testFilterBit = kFALSE;

   Int_t of = 7;
   Int_t onFlyStatus = 5;
   
   Float_t dcaV0 = -1;
   Float_t cosPointing= 2;
   Float_t decayRadius= -1;
   
   AliVTrack * trackN = 0x0;
   AliVTrack * trackP = 0x0;
   
   Double_t ptotN = -1000;
   Double_t ptotP = -1000;
      
   Int_t chargeComboDeuteronPion = -1; 
   
   Double_t momentumPion[3]={0,0,0};
   Double_t momentumPionRot[3]={0,0,0};
   Double_t momentumDeuteron[3]={0,0,0};
   Double_t momentumDeuteronRot[3]={0,0,0};
   Double_t momentumMother[3] = {0,0,0};
   
   Double_t transversMomentumPion = 0;
   Double_t transversMomentumDeuteron = 0;
   
   Double_t transversMomentumMother = 0;
   Double_t longitudinalMomentumMother = 0;
   
   Double_t totalEnergyMother = 0;
   Double_t energyPion = 0;
   Double_t energyDeuteron = 0;
   
   Double_t rapidity = 2;
   
   TVector3 vecPion(0,0,0);
   TVector3 vecDeuteron(0,0,0);
   TVector3 vecMother(0,0,0);
   
   Double_t alpha = 2;
   Double_t qt = -1;
   
   Double_t thetaPion = 0;
   Double_t thetaDeuteron = 0;
   
   Double_t phi=0;
   Double_t invaMassDeuteronPion = 0;
   
   Double_t e12 = 0;
   Double_t r12   = 0;
   Double_t d22 = 0;
   Double_t dr22   = 0;

   //Tree variables
   fV0finder[fItrk]         = -1;
   fkMB[fItrk]              = -1;
   fkCentral[fItrk]         = -1;
   fkSemiCentral[fItrk]     = -1;
   fkEMCEJE[fItrk]          = -1;
   fkEMCEGA[fItrk]          = -1;
   
   fPtotN[fItrk]            = -1000;
   fPtotP[fItrk]            = -1000;
   fMotherPt[fItrk]         = -1000;
   
   fdEdxN[fItrk]            = -1;
   fdEdxP[fItrk]            = -1;
   fSignN[fItrk]            = 0;
   fSignP[fItrk]            = 0;
   
   fDCAv0[fItrk]            = -1;
   fCosinePAv0[fItrk]       = -2;
   fDecayRadiusTree[fItrk]  = -1;
   
   fInvaMassDeuteronPionTree[fItrk] = 0;
   fChargeComboDeuteronPionTree[fItrk] = -1;
   
   fAmenterosAlphaTree[fItrk] = 2;
   fAmenterosQtTree[fItrk] = -1;
   
   //Get v0 object
   if(fAnalysisType == "ESD")v0ESD = fESDevent->GetV0(ivertex);
   if(fAnalysisType == "AOD")v0AOD = fAODevent->GetV0(ivertex);
   
   
   //distinguish between the two V0 finders: 0 offline V0 finder / 1 online V0 finder
   
   if(fAnalysisType == "ESD") of = v0ESD->GetOnFlyStatus();
   if(fAnalysisType == "AOD") of = v0AOD->GetOnFlyStatus();
   
   if(of)  onFlyStatus= 1;
   if(!of) onFlyStatus= 0;
   
   if(onFlyStatus==0)fof->Fill(1,0);
   if(onFlyStatus==1)fof->Fill(1,1);
   
   //Get dca, cos of pointing angle and decay radius
   
   
   if(fAnalysisType == "ESD")
     { 
       dcaV0 = v0ESD->GetDcaV0Daughters();
       cosPointing = v0ESD->GetV0CosineOfPointingAngle();
       decayRadius = v0ESD->GetRr();
     }
   
   if(fAnalysisType == "AOD")
     { 
       dcaV0 = v0AOD->DcaV0Daughters();
       cosPointing = v0AOD->CosPointingAngle(vertex);
       decayRadius = v0AOD->DecayLengthV0(vertex);
     }


      // select coresponding tracks
      if(fAnalysisType == "ESD")
	{ 
	  trackN = fESDevent->GetTrack(v0ESD->GetIndex(0));	  
	  trackP = fESDevent->GetTrack(v0ESD->GetIndex(1));

	  if (trackN->Charge() > 0) // switch because of bug in V0 interface
	    { 
	      trackN = fESDevent->GetTrack(v0ESD->GetIndex(1));
	      trackP = fESDevent->GetTrack(v0ESD->GetIndex(0));
	      v0ChargesCorrect = kFALSE;
	    }
	  //Track-Cuts
	  testTrackCuts = TrackCuts(trackN,testTrackCuts);
	  if(testTrackCuts == kFALSE) continue;
	  testTrackCuts = TrackCuts(trackP,testTrackCuts);
	  if(testTrackCuts == kFALSE) continue;
	}
      
      if(fAnalysisType == "AOD")
	{ 
	  trackN = dynamic_cast<AliVTrack*>(v0AOD->GetDaughter(0));
	  trackP = dynamic_cast<AliVTrack*>(v0AOD->GetDaughter(1));

	  if (trackN->Charge() > 0) // switch because of bug in V0 interface
	    { 
	      trackN = dynamic_cast<AliVTrack*>(v0AOD->GetDaughter(1));
	      trackP = dynamic_cast<AliVTrack*>(v0AOD->GetDaughter(0));
	      v0ChargesCorrect = kFALSE;
	    }
	  //Test-Filterbit - only use for track base analysis -- not for V0 candidates -- the AOD V0 candidates sould be a copy of the ESD candidates, BUT NOT for AOD0086 -> here there wqas a wider cos-cutto have more candidates in the AODs  
	  //testFilterBit = FilterBit(trackN,testFilterBit);
	  //if(testFilterBit == kFALSE) continue;
	  //testFilterBit = FilterBit(trackP,testFilterBit);
	  //if(testFilterBit == kFALSE) continue;	
	  //Track-Cuts
	  testTrackCuts = TrackCuts(trackN,testTrackCuts);
	  if(testTrackCuts == kFALSE) continue;
	  testTrackCuts = TrackCuts(trackP,testTrackCuts);
	  if(testTrackCuts == kFALSE) continue;
	}
      
      //Get the total momentum for each track ---> at the inner readout of the TPC???? // momenta a always positive
            if(fAnalysisType == "AOD") {
	ptotN = trackN->P(); //GetInnerParam()->GetP(); ******************FIXME: InnerParam
	ptotP = trackP->P(); //GetInnerParam()->GetP(); ******************FIXME: InnerParam
      }

      if(fAnalysisType == "ESD") {
	ptotN =  MomentumInnerParam(trackN,ptotN); 
	ptotP =  MomentumInnerParam(trackP,ptotP); 
      }
      
      //fill QA dEdx with all V0 candidates      
      if(trackP) fHistDeDxQA->Fill(ptotP*trackP->Charge(), trackP->GetTPCsignal(),0);
      if(trackN) fHistDeDxQA->Fill(ptotN*trackN->Charge(), trackN->GetTPCsignal(),1);

      if (dcaV0 > 3 || dcaV0 < 1e-20) continue;
      
      //Check how much the dca cut reduces the statistic (background) for the different VO-finder
      if(onFlyStatus==0)fof->Fill(2,0);
      if(onFlyStatus==1)fof->Fill(2,1);
    
      if (cosPointing < 0.9) continue;
        
      //Check how much the cos-of-the-pointing-angle-cut reduces the statistic (background) for the different VO-finder
      if(onFlyStatus==0)fof->Fill(3,0);
      if(onFlyStatus==1)fof->Fill(3,1);

      //deuteron PID via specific energy loss in the TPC
      Bool_t isDeuteron[3] = {kFALSE,kFALSE,kFALSE}; //0 = posDeuteron, 1 = negDeuteron, 2 = trackN is deuteron

      DeuteronPID(trackP,trackN,ptotP,ptotN,runNumber,isDeuteron);
      
      //if(!isDeuteron[0] && !isDeuteron[1]) continue;
      if((isDeuteron[0]==kFALSE) && (isDeuteron[1]==kFALSE)) continue;


      //Check how much the nuclei PID cuts reduce the statistics (background) for the two V0-finders
      if(onFlyStatus==0)fof->Fill(4,0);
      if(onFlyStatus==1)fof->Fill(4,1);
      
      //Fill the QA dEdx with deuterons and helium3 after the nuclei PID cut
      if(isDeuteron[1]) fHistDeDxQA->Fill(ptotN*trackN->Charge(), trackN->GetTPCsignal(),2);
      if(isDeuteron[0]) fHistDeDxQA->Fill(ptotP*trackP->Charge(), trackP->GetTPCsignal(),3);

      //deuteron PID via specific energy loss in the TPC
      Bool_t isPion[2] = {kFALSE,kFALSE}; //0 = posPion, 1 = negPion
      
      PionPID(trackP,trackN,ptotP,ptotN,runNumber,isPion);

      //if(isDeuteron[0] && !isPion[1]) continue; //pos deuteron and neg Pion
      //if(isDeuteron[1] && !isPion[0]) continue; //neg deuteron and pos Pion

      //Check how much the pion PID cuts reduce the statistics (background) for the two V0-finders
      if(onFlyStatus==0)fof->Fill(5,0);
      if(onFlyStatus==1)fof->Fill(5,1);


      //Save the different charge combinations to differentiat between particles and anit-particles:
      // -/+ = Anti-deuteron + positive Pion
      // -/- = Anti-deuteron + negative Pion
      // +/- = deuteron + negative Pion
      // +/+ = deuteron + positive Pion

      //Like-sign
      if (trackN->Charge()<0 && trackP->Charge()<0 && isDeuteron[1]==kTRUE && isPion[1]==kTRUE) chargeComboDeuteronPion = 1; // -/-
      if (trackN->Charge()>0 && trackP->Charge()>0 && isDeuteron[0]==kTRUE && isPion[0]==kTRUE ) chargeComboDeuteronPion = 3; // +/+

      //unlinke-sign
      if (trackN->Charge()<0 && trackP->Charge()>0 && isDeuteron[1]==kTRUE && isPion[0]==kTRUE ) chargeComboDeuteronPion = 0; // -/+
      //if (trackN->Charge()>0 && trackP->Charge()<0 && isDeuteron[1]==kTRUE && isPion[0]==kTRUE ) chargeComboDeuteronPion = 0; // -/+ //dürfte wegen charge correctur am anfang nicht existiere
      //if (trackN->Charge()>0 && trackP->Charge()<0 && isDeuteron[0]==kTRUE && isPion[1]==kTRUE) chargeComboDeuteronPion = 2; // +/- //dürfte wegen charge correctur am anfang nicht existiere
      if (trackN->Charge()<0 && trackP->Charge()>0 && isDeuteron[0]==kTRUE && isPion[1]==kTRUE) chargeComboDeuteronPion = 2; // +/-

      //if(chargeComboDeuteronPion==0) cout << "chargeN: " << trackN->Charge() << " chargeP: " << trackP->Charge() << endl;


      //Fill the QA dEdx with all selected tracks after the PID cuts
      fHistDeDxQA->Fill(ptotP*trackP->Charge(), trackP->GetTPCsignal(),6);
      fHistDeDxQA->Fill(ptotN*trackN->Charge(), trackN->GetTPCsignal(),6);

 
            
      //get the momenta of the daughters

      if(fAnalysisType == "ESD")
	{
	  if(chargeComboDeuteronPion==0){ //anti-deuteron (isDeuteron[1]==kTRUE), positives pion (isPion[0]==kTRUE), trackN gehört zu anti-deuteron, tackP gehört zu pion
	    
	    v0ESD->GetPPxPyPz(momentumPion[0], momentumPion[1], momentumPion[2]);
	    v0ESD->GetNPxPyPz(momentumDeuteron[0], momentumDeuteron[1], momentumDeuteron[2]);
	    if (!v0ChargesCorrect) {

	      v0ESD->GetNPxPyPz(momentumPion[0], momentumPion[1], momentumPion[2]);
	      v0ESD->GetPPxPyPz(momentumDeuteron[0], momentumDeuteron[1], momentumDeuteron[2]);
	    }
	  }

	  if(chargeComboDeuteronPion==2){ //deuteron (isDeuteron[0]==kTRUE), negative pion (isPion[1]==kTRUE), trackP gehört zu deuteron, tackN gehört zu pion
	    
	    v0ESD->GetNPxPyPz(momentumPion[0], momentumPion[1], momentumPion[2]);
	    v0ESD->GetPPxPyPz(momentumDeuteron[0], momentumDeuteron[1], momentumDeuteron[2]);
	    if (!v0ChargesCorrect){ 

	      v0ESD->GetPPxPyPz(momentumPion[0], momentumPion[1], momentumPion[2]);
	      v0ESD->GetNPxPyPz(momentumDeuteron[0], momentumDeuteron[1], momentumDeuteron[2]);
	    }
	  }

	  if(chargeComboDeuteronPion==1 || chargeComboDeuteronPion==3){ //trackN gehört zu deuteron oder pion, trackP gehört zu deuteron oder pion
	    if(isDeuteron[2]==kTRUE){ //trackN gehört zu deuteron, trackP gehört zum pion

	      v0ESD->GetPPxPyPz(momentumPion[0], momentumPion[1], momentumPion[2]);
	      v0ESD->GetNPxPyPz(momentumDeuteron[0], momentumDeuteron[1], momentumDeuteron[2]);
	      if (!v0ChargesCorrect) {
		
		v0ESD->GetNPxPyPz(momentumPion[0], momentumPion[1], momentumPion[2]);
		v0ESD->GetPPxPyPz(momentumDeuteron[0], momentumDeuteron[1], momentumDeuteron[2]);
	      }
	    }
	    if(isDeuteron[2]==kFALSE){ //trackP gehört zum deuteron, trackN gehört zum pion

	      v0ESD->GetNPxPyPz(momentumPion[0], momentumPion[1], momentumPion[2]);
	      v0ESD->GetPPxPyPz(momentumDeuteron[0], momentumDeuteron[1], momentumDeuteron[2]);
	      if (!v0ChargesCorrect) {
		
		v0ESD->GetPPxPyPz(momentumPion[0], momentumPion[1], momentumPion[2]);
		v0ESD->GetNPxPyPz(momentumDeuteron[0], momentumDeuteron[1], momentumDeuteron[2]);
	      }
	    } 
	  }	 

 
	  //get the momenta of the mother
	  v0ESD->GetPxPyPz(momentumMother[0],momentumMother[1],momentumMother[2]);
	}
    

      if(fAnalysisType == "AOD")
	{
	  if(chargeComboDeuteronPion==0){ //anti-deuteron (isDeuteron[1]==kTRUE), positives pion (isPion[0]==kTRUE), trackN gehört zu anti-deuteron, tackP gehört zu pion 
	    
	    momentumPion[0] = v0AOD->MomPosX();
	    momentumPion[1] = v0AOD->MomPosY();
	    momentumPion[2] = v0AOD->MomPosZ();
	    
	    momentumDeuteron[0] = v0AOD->MomNegX();
	    momentumDeuteron[1] = v0AOD->MomNegY();
	    momentumDeuteron[2] = v0AOD->MomNegZ();
	    
	    if (!v0ChargesCorrect){ 

	      momentumPion[0] = v0AOD->MomNegX();
	      momentumPion[1] = v0AOD->MomNegY();
	      momentumPion[2] = v0AOD->MomNegZ();
	      
	      momentumDeuteron[0] = v0AOD->MomPosX();
	      momentumDeuteron[1] = v0AOD->MomPosY();
	      momentumDeuteron[2] = v0AOD->MomPosZ();
	    }
	  }
	  
	  if (chargeComboDeuteronPion==2){ //deuteron (isDeuteron[0]==kTRUE), negative pion (isPion[1]==kTRUE), trackP gehört zu deuteron, tackN gehört zu pion

	    momentumPion[0] = v0AOD->MomNegX();
	    momentumPion[1] = v0AOD->MomNegY();
	    momentumPion[2] = v0AOD->MomNegZ();
	    
	    momentumDeuteron[0] = v0AOD->MomPosX();
	    momentumDeuteron[1] = v0AOD->MomPosY();
	    momentumDeuteron[2] = v0AOD->MomPosZ();
	    
	    if (!v0ChargesCorrect){

	      momentumPion[0] = v0AOD->MomPosX();
	      momentumPion[1] = v0AOD->MomPosY();
	      momentumPion[2] = v0AOD->MomPosZ();
	      
	      momentumDeuteron[0] = v0AOD->MomNegX();
	      momentumDeuteron[1] = v0AOD->MomNegY();
	      momentumDeuteron[2] = v0AOD->MomNegZ();
	    }
	  }

	  if(chargeComboDeuteronPion==1 || chargeComboDeuteronPion==3){ //trackN gehört zu deuteron oder pion, trackP gehört zu deuteron oder pion
	    if(isDeuteron[2]==kTRUE){ //trackN gehört zu deuteron, trackP gehört zum pion
	  
	      momentumPion[0] = v0AOD->MomPosX();
	      momentumPion[1] = v0AOD->MomPosY();
	      momentumPion[2] = v0AOD->MomPosZ();
	      
	      momentumDeuteron[0] = v0AOD->MomNegX();
	      momentumDeuteron[1] = v0AOD->MomNegY();
	      momentumDeuteron[2] = v0AOD->MomNegZ();
	      
	      if (!v0ChargesCorrect){ 
		
		momentumPion[0] = v0AOD->MomNegX();
		momentumPion[1] = v0AOD->MomNegY();
		momentumPion[2] = v0AOD->MomNegZ();
		
		momentumDeuteron[0] = v0AOD->MomPosX();
		momentumDeuteron[1] = v0AOD->MomPosY();
		momentumDeuteron[2] = v0AOD->MomPosZ();
	      }
	    }
	   
	    if(isDeuteron[2]==kFALSE){ //trackP gehört zum deuteron, trackN gehört zum pion
	      
	      momentumPion[0] = v0AOD->MomNegX();
	      momentumPion[1] = v0AOD->MomNegY();
	      momentumPion[2] = v0AOD->MomNegZ();
	      
	      momentumDeuteron[0] = v0AOD->MomPosX();
	      momentumDeuteron[1] = v0AOD->MomPosY();
	      momentumDeuteron[2] = v0AOD->MomPosZ();
	      
	      if (!v0ChargesCorrect){
		
		momentumPion[0] = v0AOD->MomPosX();
		momentumPion[1] = v0AOD->MomPosY();
		momentumPion[2] = v0AOD->MomPosZ();
		
		momentumDeuteron[0] = v0AOD->MomNegX();
		momentumDeuteron[1] = v0AOD->MomNegY();
		momentumDeuteron[2] = v0AOD->MomNegZ();
	      }
	    }

	    //get the momenta of the mother
	    momentumMother[0] = v0AOD->MomV0X();
	    momentumMother[1] = v0AOD->MomV0Y();
	    momentumMother[2] = v0AOD->MomV0Z();
	  }
	}
      
      //Rapidity - cut    
      transversMomentumPion      = TMath::Sqrt(momentumPion[0]*momentumPion[0]+momentumPion[1]*momentumPion[1]);
      transversMomentumDeuteron  = TMath::Sqrt(momentumDeuteron[0]*momentumDeuteron[0]+momentumDeuteron[1]*momentumDeuteron[1]);
  
      transversMomentumMother=TMath::Sqrt((momentumDeuteron[0]+momentumPion[0])*(momentumDeuteron[0]+momentumPion[0])+(momentumDeuteron[1]+momentumPion[1])*(momentumDeuteron[1]+momentumPion[1]));
      
      longitudinalMomentumMother = TMath::Sqrt((momentumDeuteron[2]+momentumPion[2])*(momentumDeuteron[2]+momentumPion[2]));

      energyDeuteron =  TMath::Sqrt(momentumDeuteron[0]*momentumDeuteron[0]+momentumDeuteron[1]*momentumDeuteron[1]+momentumDeuteron[2]*momentumDeuteron[2]+fgkMass[kMassDeuteron]*fgkMass[kMassDeuteron]);
    
      energyPion = TMath::Sqrt(momentumPion[0]*momentumPion[0]+momentumPion[1]*momentumPion[1]+momentumPion[2]*momentumPion[2]+fgkMass[kMassPion]*fgkMass[kMassPion]);
    
      totalEnergyMother = energyPion + energyDeuteron;
		 		    
      //rapidity cut +- 1
      rapidity = 0.5 * TMath::Log((totalEnergyMother+longitudinalMomentumMother)/( totalEnergyMother-longitudinalMomentumMother));
    
      if(rapidity > 1.0 || rapidity < -1) continue;

    
      //Armanteros-Podolanski
      vecPion.SetXYZ(momentumPion[0],momentumPion[1],momentumPion[2]);
      vecDeuteron.SetXYZ(momentumDeuteron[0],momentumDeuteron[1],momentumDeuteron[2]);
      vecMother.SetXYZ(momentumMother[0],momentumMother[1],momentumMother[2]);
  
      thetaPion = acos((vecPion * vecMother)/(vecPion.Mag() * vecMother.Mag()));
      thetaDeuteron = acos((vecDeuteron * vecMother)/(vecDeuteron.Mag() * vecMother.Mag()));

      alpha = ((vecPion.Mag())*cos(thetaPion)-(vecDeuteron.Mag())*cos(thetaDeuteron))/
	((vecDeuteron.Mag())*cos(thetaDeuteron)+(vecPion.Mag())*cos(thetaPion)) ;
      qt = vecDeuteron.Mag()*sin(thetaDeuteron);


      //Rotation for background calculation
      //Int_t rotation=1; // =1 signal, =2 Rotation of the pion , =3 Rotation of the deuteron

      Double_t fStartAnglePhi=TMath::Pi();
      Double_t fConeAnglePhi=TMath::Pi(); //-0.174;
      phi  = fStartAnglePhi+(2*gRandom->Rndm()-1)*fConeAnglePhi;
     
      for(Int_t rotation=1;rotation<4;rotation++){ //loop for rotation
        
	//Double_t fStartAnglePhi=TMath::Pi();
	//Double_t fConeAnglePhi=TMath::Pi(); //-0.174;
	//phi  = fStartAnglePhi+(2*gRandom->Rndm()-1)*fConeAnglePhi;
          
	//calculate new rotated momenta
	momentumPionRot[0]=TMath::Cos(phi)*momentumPion[0]-TMath::Sin(phi)*momentumPion[1];
	momentumPionRot[1]=TMath::Sin(phi)*momentumPion[0]+TMath::Cos(phi)*momentumPion[1];
	
	momentumDeuteronRot[0]=TMath::Cos(phi)*momentumDeuteron[0]-TMath::Sin(phi)*momentumDeuteron[1];
	momentumDeuteronRot[1]=TMath::Sin(phi)*momentumDeuteron[0]+TMath::Cos(phi)*momentumDeuteron[1];
          
	//invariant mass calculations  
	fIsCorrectlyAssociated[fItrk] = kFALSE;
	
	//factor for the invariant mass calculation, which only include the pion
	e12   = fgkMass[kMassPion]*fgkMass[kMassPion]
	  +momentumPion[0]*momentumPion[0]
	  +momentumPion[1]*momentumPion[1]
	  +momentumPion[2]*momentumPion[2];
	
	r12   = fgkMass[kMassPion]*fgkMass[kMassPion]
	  +momentumPionRot[0]*momentumPionRot[0]
	  +momentumPionRot[1]*momentumPionRot[1]
	  +momentumPion[2]*momentumPion[2];
	
	//factor for the invariant mass calculation, which only include the deuterons
	d22   = fgkMass[kMassDeuteron]*fgkMass[kMassDeuteron]
	  +momentumDeuteron[0]*momentumDeuteron[0]
	  +momentumDeuteron[1]*momentumDeuteron[1]
	  +momentumDeuteron[2]*momentumDeuteron[2];
	dr22   = fgkMass[kMassDeuteron]*fgkMass[kMassDeuteron]
	  +momentumDeuteronRot[0]*momentumDeuteronRot[0]
	  +momentumDeuteronRot[1]*momentumDeuteronRot[1]
	  +momentumDeuteron[2]*momentumDeuteron[2];

	if(rotation == 1){ //signal
	  invaMassDeuteronPion = TMath::Sqrt(TMath::Max(fgkMass[kMassPion]*fgkMass[kMassPion]
							+fgkMass[kMassDeuteron]*fgkMass[kMassDeuteron] 
							+ 2.*(TMath::Sqrt(e12*d22)
							      -momentumPion[0]*momentumDeuteron[0]
							      -momentumPion[1]*momentumDeuteron[1]
							      -momentumPion[2]*momentumDeuteron[2]), 0.));

	  if (fMCtrue){
	    
	    Int_t labelN = 0;
	    Int_t labelP = 0;
	    labelN = trackN->GetLabel();
	    labelP = trackP->GetLabel();
	    
	    TParticle *tparticleDaughterN = stack->Particle(TMath::Abs(labelN));
	    TParticle *tparticleDaughterP = stack->Particle(TMath::Abs(labelP));
	    
	    Int_t labelMotherN = tparticleDaughterN->GetFirstMother();
	    Int_t labelMotherP = tparticleDaughterP->GetFirstMother();
	    
	    TParticle *tparticleMotherN = stack->Particle(TMath::Abs(labelMotherN));
	    TParticle *tparticleMotherP = stack->Particle(TMath::Abs(labelMotherP));
	    
	    //LambdaNeuteron
	    if(tparticleMotherN->GetPdgCode() == fgkPdgCode[kPDGLambdaNeutron] && tparticleMotherP->GetPdgCode() == fgkPdgCode[kPDGLambdaNeutron] && onFlyStatus ==1 && labelMotherN == labelMotherP ){//check mother PDG  and fill the histogramms only for the only V0 finder
	      
	      //cout << "number of daughters: " << tparticleMotherN->GetNDaughters() << endl;;
	      if(tparticleMotherN->GetNDaughters() > 2.) continue;
	      
	      Int_t labelSecondDaughter = tparticleMotherP->GetDaughter(1);
	      Int_t labelFirstDaughter = labelSecondDaughter-1;
	  
	      TParticle *tparticleFirstDaughter = stack->Particle(TMath::Abs(labelFirstDaughter));
	      TParticle *tparticleSecondDaughter = stack->Particle(TMath::Abs(labelSecondDaughter));
	      
	      if(tparticleFirstDaughter->GetPdgCode() == fgkPdgCode[kPDGDeuteron]){//check first daughter PDG
		
		if(tparticleSecondDaughter->GetPdgCode() == fgkPdgCode[kPDGPionMinus]){//check second daughter PDG
		  
		  if(invaMassDeuteronPion < 2.02) continue;
		  
		  fHistLambdaNeutronPtAso->Fill(transversMomentumMother);
		  fHistLambdaNeutronInvaMassAso->Fill(invaMassDeuteronPion);
		  fIsCorrectlyAssociated[fItrk] = kTRUE;
		  fHistLambdaNeutronDecayLengthAso->Fill((decayRadius*2.054)/(tparticleMotherP->P()));
		  //CUTS
		  
		  if((dcaV0 < 0.5) && 
		     (cosPointing > 0.999) && 
		     (decayRadius < 50.0) && 
		     (ptotP > 0.2)){//trackP->GetTPCsignal()>180 &&&&  decayRadius > 1.5 && decayRadius< 50 
		    
		    fHistLambdaNeutronPtAsoCuts->Fill(transversMomentumMother);
		  }
		}//end check second daughter PDG
	      }//end check first daughter PDG
	    }//end LambdaNeutron
	    
	    //Anti-LambdaNeutron
	    if(tparticleMotherN->GetPdgCode() == fgkPdgCode[kPDGAntiLambdaNeutron] && tparticleMotherP->GetPdgCode() == fgkPdgCode[kPDGAntiLambdaNeutron] && onFlyStatus ==1 && labelMotherN == labelMotherP){//check mother PDG  and fill the histogramms only for the only V0 finder
	      
	      Int_t labelSecondDaughter = tparticleMotherN->GetDaughter(1);
	      Int_t labelFirstDaughter = labelSecondDaughter-1;
	      
	      TParticle *tparticleFirstDaughter = stack->Particle(TMath::Abs(labelFirstDaughter));
	      TParticle *tparticleSecondDaughter = stack->Particle(TMath::Abs(labelSecondDaughter));
	      
	      if(tparticleFirstDaughter->GetPdgCode() == fgkPdgCode[kPDGAntiDeuteron]) {//check first daughter PDG
		
		if(tparticleSecondDaughter->GetPdgCode() == fgkPdgCode[kPDGPionPlus]) {//check second daughter PDG
		  
		  if(invaMassDeuteronPion < 2.02) continue;
		  
		  fIsCorrectlyAssociated[fItrk] = kTRUE;
		  fHistAntiLambdaNeutronPtAso->Fill(transversMomentumMother);
		  fHistAntiLambdaNeutronInvaMassAso->Fill(invaMassDeuteronPion);
		  fHistAntiLambdaNeutronDecayLengthAso->Fill((decayRadius*2.054)/(tparticleMotherP->P()));
		  
		  //CUTS
		  if(dcaV0 < 1. && cosPointing > 0.99 && ptotN > 0.2)//&& trackN->GetTPCsignal()>180 && decayRadius > 1.5 && decayRadius< 50
		    {
		      fHistAntiLambdaNeutronPtAsoCuts->Fill(transversMomentumMother);
		    }
		}//end check second daughter PDG
	      }//end check first daughter PDG
	    }//end Anti-LambdaNeutron
	  }//end MC
	}//end rotation == 1, signal
	
	if(rotation == 2){ // rotation of the pion
	  invaMassDeuteronPion = TMath::Sqrt(TMath::Max(fgkMass[kMassPion]*fgkMass[kMassPion]
							+fgkMass[kMassDeuteron]*fgkMass[kMassDeuteron]
							+ 2.*(TMath::Sqrt(r12*d22)
							      -momentumPionRot[0]*momentumDeuteron[0]
							      -momentumPionRot[1]*momentumDeuteron[1]
							      -momentumPion[2]*momentumDeuteron[2]), 0.));
	}//end rotation == 2, rotation of the pion
	
	if(rotation == 3){// Rotation of the deuteron
	  invaMassDeuteronPion = TMath::Sqrt(TMath::Max(fgkMass[kMassPion]*fgkMass[kMassPion]
							+fgkMass[kMassDeuteron]*fgkMass[kMassDeuteron] 
							+ 2.*(TMath::Sqrt(e12*dr22)
							      -momentumPion[0]*momentumDeuteronRot[0]
							      -momentumPion[1]*momentumDeuteronRot[1]
							      -momentumPion[2]*momentumDeuteron[2]), 0.));
	}//end rotation == 3, rotation of the deuteron
	
	//fill the THnSparse and the tree variables
	
	//tree variables which are independent of the particle-species
	fV0finder[fItrk]         = onFlyStatus;
	fkMB[fItrk]              = fTriggerFired[0];
	fkCentral[fItrk]         = fTriggerFired[1];
	fkSemiCentral[fItrk]     = fTriggerFired[2];
	fkEMCEJE[fItrk]          = fTriggerFired[3];
	fkEMCEGA[fItrk]          = fTriggerFired[4];
	
	fPtotN[fItrk]            = trackN->P(); //InnerParam???
	fPtotP[fItrk]            = trackP->P(); //InnerParam???
	fMotherPt[fItrk]         = transversMomentumMother;
	
	fdEdxN[fItrk]            = trackN->GetTPCsignal();
	fdEdxP[fItrk]            = trackP->GetTPCsignal();
	fSignN[fItrk]            = trackN->Charge();
	fSignP[fItrk]            = trackP->Charge();
	
	fDCAv0[fItrk]            = dcaV0;
	fCosinePAv0[fItrk]       = cosPointing;
	fDecayRadiusTree[fItrk]  = decayRadius;
	
	fAmenterosAlphaTree[fItrk] = alpha;
	fAmenterosQtTree[fItrk] = qt;
	fRotationTree[fItrk] = rotation;
	
            
	if (isDeuteron[0] == kTRUE)  //pos deuteron
	  {
	    fInvaMassDeuteronPionTree[fItrk] = invaMassDeuteronPion;
	    fChargeComboDeuteronPionTree[fItrk] = chargeComboDeuteronPion;
	    
	    fItrk++;
	    
	    if(invaMassDeuteronPion < 2.1) 
	      {
		if(chargeComboDeuteronPion == 2)fHistDeDxQA->Fill(ptotP*(trackP->Charge()), trackP->GetTPCsignal(),9);
		if(chargeComboDeuteronPion == 2)fHistDeDxQA->Fill(ptotN*(trackN->Charge()), trackN->GetTPCsignal(),10);
	      }
	    //fHistArmenterosPodolanskiDeuteronPion->Fill(alpha*(-1),qt);
	  }

	if (isDeuteron[1] == kTRUE) 
	  {
	    fInvaMassDeuteronPionTree[fItrk] = invaMassDeuteronPion;
	    fChargeComboDeuteronPionTree[fItrk] = chargeComboDeuteronPion;
	    
	    fItrk++;
	    
	    if(invaMassDeuteronPion < 2.1) 
	      {
		if(chargeComboDeuteronPion == 0)fHistDeDxQA->Fill(ptotN*trackN->Charge(), trackN->GetTPCsignal(),7);
		if(chargeComboDeuteronPion == 0)fHistDeDxQA->Fill(ptotP*trackP->Charge(), trackP->GetTPCsignal(),8);
	      }
	    //fHistArmenterosPodolanskiAntiDeuteronPion->Fill(alpha,qt);
	  }
	
      }//end rotation loop
      
 } //end of v0 loop

  //fill the tree
  fTreeV0->Fill();
 
  // Post output data.
  PostData(1, fOutputContainer);
  PostData(2, fTreeV0);
}      

//________________________________________________________________________
void AliAnalysisTaskLambdaNAOD::Terminate(const Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query

  //get output data and darw 'fHistPt'
  if (!GetOutputData(0)) return;
  //TH1F *hist=(TH1F*)(((TObjArray*)GetOutputData(0))->FindObject("fHistPt"));
  //if (hist) hist->Draw();
}

//_____________________________________________________________________________
Int_t AliAnalysisTaskLambdaNAOD::Initialize() {

 
  // -- Reset Event
  // ----------------
  ResetEvent();

  return 0;
}
//________________________________________________________________________
Int_t AliAnalysisTaskLambdaNAOD::SetupEvent() {
  // Setup Reading of event

  ResetEvent();

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   
  // -- Get Trigger 
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   

  //Bool_t isTriggered = IsTriggered();
  IsTriggered();
  return 0;
}
//_____________________________________________________________________________
void AliAnalysisTaskLambdaNAOD::ResetEvent() {

  // Reset event

  // -- Reset QA
  for (Int_t ii = 0; ii < fNTriggers; ++ii)
    fTriggerFired[ii] = kFALSE;
  
  return;
}
//_______________________________________________________________________
void AliAnalysisTaskLambdaNAOD::BinLogAxis(const THnSparse *h, Int_t axisNumber) {
  //
  // Method for the correct logarithmic binning of histograms
  //
  TAxis *axis = h->GetAxis(axisNumber);
  int bins = axis->GetNbins();

  Double_t from = axis->GetXmin();
  Double_t to = axis->GetXmax();
  Double_t *newBins = new Double_t[bins + 1];
   
  newBins[0] = from;
  Double_t factor = pow(to/from, 1./bins);
  
  for (int i = 1; i <= bins; i++) {
   newBins[i] = factor * newBins[i-1];
  }
  axis->Set(bins, newBins);
  delete [] newBins;
  
}
//________________________________________________________________________
Bool_t AliAnalysisTaskLambdaNAOD::IsTriggered() {
  // Check if Event is triggered and fill Trigger Histogram

  if ((fEventHandler->IsEventSelected() & AliVEvent::kMB))          fTriggerFired[0] = kTRUE;
  if ((fEventHandler->IsEventSelected() & AliVEvent::kCentral))     fTriggerFired[1] = kTRUE;
  if ((fEventHandler->IsEventSelected() & AliVEvent::kSemiCentral)) fTriggerFired[2] = kTRUE;
  if ((fEventHandler->IsEventSelected() & AliVEvent::kEMCEJE))      fTriggerFired[3] = kTRUE;
  if ((fEventHandler->IsEventSelected() & AliVEvent::kEMCEGA))      fTriggerFired[4] = kTRUE;

  Bool_t isTriggered = kFALSE;

  for (Int_t ii=0; ii<fNTriggers; ++ii) {
    if(fTriggerFired[ii]) {
      isTriggered = kTRUE;
      fHistTriggerStat->Fill(ii);
    }
  }

  return isTriggered;
  }
//______________________________________________________________________________
Bool_t AliAnalysisTaskLambdaNAOD::DeuteronPID(AliVTrack *trackP, AliVTrack *trackN, Double_t ptotP, Double_t ptotN, Int_t runNumber, Bool_t isDeuteron[3]) {
 
     //define the arrays for the Bethe-Bloch-Parameters
      Double_t parDeuteron[5] = {0,0,0,0,0};

      if(runNumber < 166500) //LHC10h
	{
	  parDeuteron[0] = 1.45802; // ALEPH parameters for deuterons (pass2)
	  parDeuteron[1] = 27.4992;
	  parDeuteron[2] = 4.00313e-15;
	  parDeuteron[3] = 2.48485;
	  parDeuteron[4] = 8.31768; 
	}
      
      if(runNumber > 166500) //LHC11h
	{ 
	  parDeuteron[0] = 1.11243; // ALEPH parameters for deuterons (pass2)
	  parDeuteron[1] = 26.1144;
	  parDeuteron[2] = 4.00313e-15;
	  parDeuteron[3] = 2.72969 ;
	  parDeuteron[4] = 9.15038; 
	}
 

      //define expected signals for the various species
      Double_t expSignalDeuteronN = 0;
      Double_t expSignalDeuteronP = 0;
  
      isDeuteron[0] = kFALSE;
      isDeuteron[1] = kFALSE;
      isDeuteron[2] = kFALSE;
      
      //for data
      if(!fMCtrue){
       
	expSignalDeuteronN = AliExternalTrackParam::BetheBlochAleph(ptotN/(fgkMass[kMassDeuteron]),parDeuteron[0],parDeuteron[1],parDeuteron[2],parDeuteron[3],parDeuteron[4]);
	expSignalDeuteronP = AliExternalTrackParam::BetheBlochAleph(ptotP/(fgkMass[kMassDeuteron]),parDeuteron[0],parDeuteron[1],parDeuteron[2],parDeuteron[3],parDeuteron[4]);
	
	if(trackP->GetTPCsignal() >= 110 && ///??????????????????????????????????
	   trackP->GetTPCsignal() < 1200 && 
	   (TMath::Abs(trackP->GetTPCsignal() - expSignalDeuteronP)/expSignalDeuteronP) < 0.2 &&
	   ptotP > 0.2 ){
	  
	  if(trackP->Charge() >0)	 	isDeuteron[0] = kTRUE; //pos deuteron
	  if(trackP->Charge() <0)	 	isDeuteron[1] = kTRUE; //neg deuteron
	}
	
	if(trackN->GetTPCsignal() >= 110 && ///??????????????????????????????????
	   trackN->GetTPCsignal() < 1200 && 
	   (TMath::Abs(trackN->GetTPCsignal() - expSignalDeuteronN)/expSignalDeuteronN) < 0.2 &&
	   ptotN > 0.2){ 
	  
	  isDeuteron[2] = kTRUE;
	  
	  if(trackN->Charge() >0)	 isDeuteron[0] = kTRUE; //pos deuteron
	  if(trackN->Charge() <0)	 isDeuteron[1] = kTRUE; //neg deuteron
	}
      }
      
      //for MC
      if(fMCtrue)
	{
	  if(runNumber < 166500) //2010
	    { 
	      expSignalDeuteronN = 0.65*AliExternalTrackParam::BetheBlochAleph(ptotN/(fgkMass[kMassDeuteron]),parDeuteron[0],parDeuteron[1],parDeuteron[2],parDeuteron[3],parDeuteron[4]);
	      expSignalDeuteronP = 0.65*AliExternalTrackParam::BetheBlochAleph(ptotP/(fgkMass[kMassDeuteron]),parDeuteron[0],parDeuteron[1],parDeuteron[2],parDeuteron[3],parDeuteron[4]); 
	    }
	  if(runNumber > 166500) //2011
	    { 
	      expSignalDeuteronN = 0.8*AliExternalTrackParam::BetheBlochAleph(ptotN/(fgkMass[kMassDeuteron]),parDeuteron[0],parDeuteron[1],parDeuteron[2],parDeuteron[3],parDeuteron[4]);
	      expSignalDeuteronP = 0.8*AliExternalTrackParam::BetheBlochAleph(ptotP/(fgkMass[kMassDeuteron]),parDeuteron[0],parDeuteron[1],parDeuteron[2],parDeuteron[3],parDeuteron[4]); 
	    }
	 
	  if(trackP->GetTPCsignal() < 1200 && 
	     (TMath::Abs(trackP->GetTPCsignal() - expSignalDeuteronP)/expSignalDeuteronP) < 0.2 &&
	     ptotP > 0.2 )
	    {
	      if(trackP->Charge() >0)	 	isDeuteron[0] = kTRUE; //pos deuteron
	      if(trackP->Charge() <0)	 	isDeuteron[1] = kTRUE; //neg deuteron	      
	    }

	  if(trackN->GetTPCsignal() < 1200 && 
	     (TMath::Abs(trackN->GetTPCsignal() - expSignalDeuteronN)/expSignalDeuteronN) < 0.2 &&
	     ptotN > 0.2 )
	    {
	      isDeuteron[2] = kTRUE;

	      if(trackN->Charge() >0)	 	isDeuteron[0] = kTRUE; //pos deuteron
	      if(trackN->Charge() <0)	 	isDeuteron[1] = kTRUE; //neg deuteron	      
	    }

	}



     return isDeuteron;
  
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskLambdaNAOD::PionPID(AliVTrack *trackP, AliVTrack *trackN, Double_t ptotP, Double_t ptotN, Int_t runNumber, Bool_t isPion[2]) {
    
  //Pion PID via specific energy loss in the TPC
  //define the array for the Bethe-Bloch-Parameters (only necessary for local MC)
  Double_t parPion[5]  = {0,0,0,0,0};
  
    if(runNumber < 166500) { //LHC10h
        parPion[0]   = 1.45802; // ALEPH parameters for pions (pass2)
        parPion[1]   = 27.4992;
        parPion[2]   = 4.00313e-15;
        parPion[3]   = 2.48485;
        parPion[4]   = 8.31768; 
    }
      
    if(runNumber > 166500) {  //LHC11h
        parPion[0]   = 1.11243; // ALEPH parameters for pions (pass2)
        parPion[1]   = 26.1144;
        parPion[2]   = 4.00313e-15;
        parPion[3]   = 2.72969;
        parPion[4]   = 9.15038; 
    }


    Double_t expSignalPionP = 0;
    Double_t expSignalPionN = 0;

    //for MC
    if(fMCtrue){
      if(runNumber < 166500){ //2010
	expSignalPionP = 0.7*AliExternalTrackParam::BetheBlochAleph(ptotP/(fgkMass[kMassPion]),parPion[0],parPion[1],parPion[2],parPion[3],parPion[4]);
	expSignalPionN = 0.7*AliExternalTrackParam::BetheBlochAleph(ptotN/(fgkMass[kMassPion]),parPion[0],parPion[1],parPion[2],parPion[3],parPion[4]);
      }
      if(runNumber > 166500){ //2011
	expSignalPionP = AliExternalTrackParam::BetheBlochAleph(ptotP/(fgkMass[kMassPion]),parPion[0],parPion[1],parPion[2],parPion[3],parPion[4]);
	expSignalPionN = AliExternalTrackParam::BetheBlochAleph(ptotN/(fgkMass[kMassPion]),parPion[0],parPion[1],parPion[2],parPion[3],parPion[4]);
      }
      
    }
     
    isPion[0] = kFALSE;
    isPion[1] = kFALSE;
    //data
    if(!fMCtrue){
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackP, AliPID::kPion))<3){

	if(trackP->Charge()>0 )	isPion[0] = kTRUE; //pos pion
	if(trackP->Charge()<0 )	isPion[1] = kTRUE; //neg pion
     }
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackN, AliPID::kPion))<3){

	if(trackN->Charge()>0 )	isPion[0] = kTRUE; //pos pion
	if(trackN->Charge()<0 )	isPion[1] = kTRUE; //neg pion
      }
    }

    //MC
    if(fMCtrue){
      if(TMath::Abs(trackP->GetTPCsignal() - expSignalPionP)/expSignalPionP < 0.2
	 && ptotP>0.00001){
 	if(trackP->Charge()>0) isPion[0] = kTRUE; //pos pion
	if(trackP->Charge()<0) isPion[1] = kTRUE; //neg pion
      }
     if(TMath::Abs(trackN->GetTPCsignal() - expSignalPionN)/expSignalPionN < 0.2
        && ptotN>0.00001){
 	if(trackN->Charge()>0) isPion[0] = kTRUE; //pos pion
	if(trackN->Charge()<0) isPion[1] = kTRUE; //neg pion
      }
    }

    return isPion;

}
//______________________________________________________________________________
Bool_t AliAnalysisTaskLambdaNAOD::TrackCuts(AliVTrack *track, Bool_t testTrackCuts) {

  testTrackCuts = kFALSE;

  AliESDtrack *esdtrack = static_cast<AliESDtrack *>(track);
  //if(!esdtrack) testTrackCuts = kFALSE;
  if (fESDtrackCutsV0->AcceptTrack(esdtrack)) testTrackCuts = kTRUE;
//if( testTrackCuts == kTRUE) cout <<   "testTrackCuts im test: " << testTrackCuts << endl;


  return testTrackCuts;
}
//______________________________________________________________________________
/*Bool_t AliAnalysisTaskLambdaNAOD::FilterBit(AliVTrack *track, Bool_t testFilterBit) {

  testFilterBit = kFALSE;

  AliAODTrack *aodtrack = static_cast<AliAODTrack *>(track);
  //if(!aodtrack) testFilterBit = kFALSE;
  if (aodtrack->TestFilterBit(7)) testFilterBit = kTRUE;

  //if(testFilterBit == kTRUE) cout <<   "testFilterBit im test: " << testFilterBit<< endl;

  return testFilterBit;
  }*/
//______________________________________________________________________
void AliAnalysisTaskLambdaNAOD::MCGenerated(AliStack* stack) 
{ 

  // Monte Carlo for genenerated particles
  if (fMCtrue) //MC loop  
    {
 
      Int_t stackN = 0;

      for(stackN = 0; stackN < stack->GetNtrack(); stackN++) //loop over stack
	{

	  const TParticle *tparticleMother = stack->Particle(stackN);
	  Long_t pdgCodeMother = tparticleMother->GetPdgCode();

	  //check which particle the mother is 

	  //LambdaNeutron
	  if(pdgCodeMother == fgkPdgCode[kPDGLambdaNeutron]) //check mother PDG 
	    {
	      MCTwoBodyDecay(stack,tparticleMother,pdgCodeMother,fgkPdgCode[kPDGDeuteron],fgkPdgCode[kPDGPionMinus],fgkMass[kMassDeuteron],fgkMass[kMassPion]);
	    }

	  //Anti-LambdaNeutron
	  if(pdgCodeMother == fgkPdgCode[kPDGAntiLambdaNeutron]) //check mother PDG 
	    {
	      MCTwoBodyDecay(stack,tparticleMother,pdgCodeMother,fgkPdgCode[kPDGAntiDeuteron],fgkPdgCode[kPDGPionPlus],fgkMass[kMassDeuteron],fgkMass[kMassPion]);
	    }

  	      
	}//end loop over stack
      

    }//end MC
}
//_____________________________________________
void AliAnalysisTaskLambdaNAOD::MCTwoBodyDecay(AliStack* stack, const TParticle *tparticleMother, Long_t PDGMother, Long_t PDGFirstDaughter, Long_t PDGSecondDaughter, Double_t massFirstDaughter, Double_t massSecondDaughter){ //function that calculates the invariant mass and the transverse momentum for MC

  Double_t momentumFirstDaughterGen[3]={0,0,0};
  Double_t momentumSecondDaughterGen[3]={0,0,0};
  
  Double_t energyFirstDaughterGen = 0;
  Double_t energySecondDaughterGen = 0;
  
  Double_t transversMomentumMotherGen = 0;
  Double_t longitudinalMomentumMotherGen = 0;
  Double_t totalEnergyMotherGen = 0;
  
  Double_t rapidityGen = 2;
  
  //Int_t labelFirstDaughter = tparticleMother->GetDaughter(0);
  Int_t labelSecondDaughter = tparticleMother->GetDaughter(1);
  Int_t labelFirstDaughter = labelSecondDaughter -1;

  TParticle *tparticleFirstDaughter = stack->Particle(TMath::Abs(labelFirstDaughter));
  TParticle *tparticleSecondDaughter = stack->Particle(TMath::Abs(labelSecondDaughter));

  if(tparticleFirstDaughter->GetPdgCode() == PDGFirstDaughter) //check first daughter PDG
    {
      if(tparticleSecondDaughter->GetPdgCode() == PDGSecondDaughter) //check second daughter PDG
	{
	  
	  momentumFirstDaughterGen[0] = tparticleFirstDaughter->Px();
	  momentumFirstDaughterGen[1] = tparticleFirstDaughter->Py();
	  momentumFirstDaughterGen[2] = tparticleFirstDaughter->Pz();
	  
	  momentumSecondDaughterGen[0] = tparticleSecondDaughter->Px();
	  momentumSecondDaughterGen[1] = tparticleSecondDaughter->Py();
	  momentumSecondDaughterGen[2] = tparticleSecondDaughter->Pz();
	  
	  energyFirstDaughterGen  = tparticleFirstDaughter->Energy();
	  energySecondDaughterGen = tparticleSecondDaughter->Energy();
	  
	  transversMomentumMotherGen=TMath::Sqrt((momentumFirstDaughterGen[0]+momentumSecondDaughterGen[0])*(momentumFirstDaughterGen[0]+momentumSecondDaughterGen[0])+(momentumFirstDaughterGen[1]+momentumSecondDaughterGen[1])*(momentumFirstDaughterGen[1]+momentumSecondDaughterGen[1]));
	  
	  //rapidity cut +- 1
	  //longitudinal momentum of mother
	  longitudinalMomentumMotherGen = TMath::Sqrt((momentumFirstDaughterGen[2]+momentumSecondDaughterGen[2])*(momentumFirstDaughterGen[2]+momentumSecondDaughterGen[2]));
	  
	  //total energy of mother
	  totalEnergyMotherGen = energyFirstDaughterGen + energySecondDaughterGen;
	  
	  //rapidity
	  rapidityGen = 0.5 * TMath::Log( (totalEnergyMotherGen+longitudinalMomentumMotherGen)/(totalEnergyMotherGen-longitudinalMomentumMotherGen));
	  
	  if(rapidityGen > 1.0 || rapidityGen < -1 ) return;

	  //cout << "decay lenght first daughter :" << tparticleFirstDaughter->Rho()  << " decay lenght second daughter :" << tparticleSecondDaughter->Rho()<< endl;
	  //fHistLambdaNeutronDecayLengthGen->Fill(tparticleFirstDaughter->Rho());

	  //calculate the invariant mass
	  Double_t firstDaughterPart = 0;
	  Double_t secondDaughterPart = 0;
	  Double_t invaMass = 0;
		      
	  firstDaughterPart = massFirstDaughter*massFirstDaughter
	    +momentumFirstDaughterGen[0]*momentumFirstDaughterGen[0]
	    +momentumFirstDaughterGen[1]*momentumFirstDaughterGen[1]
	    +momentumFirstDaughterGen[2]*momentumFirstDaughterGen[2];
	  secondDaughterPart = massSecondDaughter*massSecondDaughter
	    +momentumSecondDaughterGen[0]*momentumSecondDaughterGen[0]
	    +momentumSecondDaughterGen[1]*momentumSecondDaughterGen[1]
	    +momentumSecondDaughterGen[2]*momentumSecondDaughterGen[2];
	  
	  invaMass = TMath::Sqrt(TMath::Max(massFirstDaughter*massFirstDaughter
							+massSecondDaughter*massSecondDaughter 
							+ 2.*(TMath::Sqrt(firstDaughterPart*secondDaughterPart)
							      -momentumFirstDaughterGen[0]*momentumSecondDaughterGen[0]
							      -momentumFirstDaughterGen[1]*momentumSecondDaughterGen[1]
							      -momentumFirstDaughterGen[2]*momentumSecondDaughterGen[2]), 0.));

#if 0

							      switch(PDGMother) {
							      case fgkPdgCode[kPDGLambdaNeutron] : 
  							       fHistLambdaNeutronPtGen->Fill(transversMomentumMotherGen);
							       fHistLambdaNeutronInvaMassGen->Fill(invaMass);
							       break;
							      case fgkPdgCode[kPDGAntiLambdaNeutron] :
							       fHistAntiLambdaNeutronPtGen->Fill(transversMomentumMotherGen);
							       fHistAntiLambdaNeutronInvaMassGen->Fill(invaMass);
							       break;
							      default :
							      printf("should not happen!!!! \n");
							      }


#else
	  //LambdaNeutron
	  if(PDGMother == fgkPdgCode[kPDGLambdaNeutron])
	    {  
	      fHistLambdaNeutronPtGen->Fill(transversMomentumMotherGen);
	      fHistLambdaNeutronInvaMassGen->Fill(invaMass);
	      fHistLambdaNeutronDecayLengthGen->Fill(((tparticleFirstDaughter->Rho())*2.054)/(tparticleMother->P()));
	    }
	  //Anti-LambdaNeutron
	  if(PDGMother == fgkPdgCode[kPDGAntiLambdaNeutron])
	    {
	      fHistAntiLambdaNeutronPtGen->Fill(transversMomentumMotherGen);
	      fHistAntiLambdaNeutronInvaMassGen->Fill(invaMass);
	      fHistAntiLambdaNeutronDecayLengthGen->Fill(((tparticleFirstDaughter->Rho())*2.054)/(tparticleMother->P()));
	    }	      

#endif	  
	}//end of check second daughter PDG
    }//end of check first daughter PDG

}
//_____________________________________________
Double_t AliAnalysisTaskLambdaNAOD::MomentumInnerParam(AliVTrack *track, Double_t ptot){ //function to get the momentum at the inner wall of the TPC
 
  AliESDtrack *esdtrack = static_cast<AliESDtrack *>(track);
  ptot = esdtrack->GetInnerParam()->GetP();

  return ptot;
} 
