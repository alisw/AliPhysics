/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 ************************************v**************************************/

#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TList.h"
#include "TTree.h"
#include "TMath.h"
#include "Riostream.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliMultSelection.h"
#include "AliPIDResponse.h"
#include "AliAODMCParticle.h"
#include "AliAODInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "AliVTrack.h"
#include "AliAODv0.h"
#include "TDatabasePDG.h"
#include "AliExternalTrackParam.h"

#include "AliAnalysisTaskK0SPFemto.h"
#include "AliAnalysisK0SPEventCollection.h"

class AliAnalysisTaskK0SPFemto;    // your analysis class

//class AliAnalysisK0SPEventCollection;

class AliPIDResponse;
class AliMultSelection;

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskK0SPFemto) // classimp: necessary for root


AliAnalysisTaskK0SPFemto::AliAnalysisTaskK0SPFemto() :AliAnalysisTaskSE(), 
  fEventCuts(0),
  fCollidingSystem("pp"),  
  Neventi(0),
  fAOD(NULL),
  fIsMC(kFALSE),
  fPIDResponse(0),
  fCentrLowLim(0),
  fCentrUpLim(100),
  fOutputContainer(NULL),
  fHistSparseSignal(0),  
  fHistSparseBkg(0),      
  farrGT(0),
  fTrackBufferSize(20200), // was 18k
  fEventColl(0x0),
  fEvt(0x0),
  fMaxFirstMult(20),//(3000), // 1000 for protons 
  fMaxSecondMult(1000),//(20),  // was 100
  fzVertexBins(10),
  fnCentBins(20),
  fnEventsToMix(7),
  fFilterBit(4),
  fHMtrigger(kFALSE),
  fPDGMassFirst(0.),  
  fPDGcodeFirst(310),
  fPDGMassSecond(0.),  
  fPDGcodeSecond(2212),
  tCentrality(0),
  tSphericity(0),
  tSpherocity(0),
  tKtpair(0),
  tkStar(0),
  tIsCommonParton(0),
  tPtV0(0),
  tPTotV0(0),
  tThetaV0(0),
  tPhiV0(0),
  tDcaPosV0(0),
  tDcaNegV0(0),
  tInvMassK0s(0),
  tInvMassLambda(0),
  tInvMassAntiLambda(0),
  tCosPointingAngleV0(0),

  tPtP(0),
  tPTotP(0),
  tThetaP(0),
  tPhiP(0),
  tSignP(0),
  tDCAxyP(0),
  tDCAzP(0),
  tMassTOFP(0),

  tMCtruepair(0),
  tMCSameMother(0),
  tMCMotherV0(0),
  tMCMotherP(0),
  tMCptcTypeV0(0),
  tMCptcTypeP(0),
  tMCSameGM(0),
  tMotherPDG(0),
  tpdgcodeV0(0),
  tpdgcodeP(0),
  tKstarGen(0),
  
  
  fHistEventMultiplicity(0),
  fHistCentrality(0),
  fHistVertexDistribution(0),
  fHistSphericity(0),
  fHistSpherocity(0),
  fHistMassK0S(0),  
  fHistFirstNPionTPCdEdx(0),
  fHistFirstPPionTPCdEdx(0),
  fHistSecondTPCdEdx(0),
  fHistSecondMassTOFvsPt3sTPC(0)
{
  
  cout<<"I'm taking this dummy constructor!!!!"<<endl;

  fPDGMassSecond = TDatabasePDG::Instance()->GetParticle(fPDGcodeSecond)->Mass(); 
  fPDGMassFirst = TDatabasePDG::Instance()->GetParticle(fPDGcodeFirst)->Mass(); 
  
  // default constructor, don't allocate memory here!
  // this is used by root for IO purpos, it needs to remain empty
}

//_____________________________________________________________________________
AliAnalysisTaskK0SPFemto::AliAnalysisTaskK0SPFemto(const char* name) : 
  AliAnalysisTaskSE(name),
  fEventCuts(0),
  fCollidingSystem("pp"),
  Neventi(0),
  fAOD(NULL),
  fIsMC(kFALSE),
  fPIDResponse(0),
  fCentrLowLim(0),
  fCentrUpLim(100),
  fOutputContainer(NULL),
  fHistSparseSignal(0),  
  fHistSparseBkg(0),      
  farrGT(0),
  fTrackBufferSize(20200), // was 18k
  fEventColl(0x0),
  fEvt(0x0),
  fMaxFirstMult(20),//(3000), // 1000 for protons 
  fMaxSecondMult(1000),//(20),  // was 100
  fzVertexBins(10),
  fnCentBins(20),
  fnEventsToMix(7),
  fFilterBit(4),
  fHMtrigger(kFALSE),
  fPDGMassFirst(0.),  
  fPDGcodeFirst(310),
  fPDGMassSecond(0.),  
  fPDGcodeSecond(2212),
  tCentrality(0),
  tSphericity(0),
  tSpherocity(0),
  tKtpair(0),
  tkStar(0),
  tIsCommonParton(0),
  tPtV0(0),
  tPTotV0(0),
  tThetaV0(0),
  tPhiV0(0),
  tDcaPosV0(0),
  tDcaNegV0(0),
  tInvMassK0s(0),
  tInvMassLambda(0),
  tInvMassAntiLambda(0),
  tCosPointingAngleV0(0),

  tPtP(0),
  tPTotP(0),
  tThetaP(0),
  tPhiP(0),
  tSignP(0),
  tDCAxyP(0),
  tDCAzP(0),
  tMassTOFP(0),

  tMCtruepair(0),
  tMCSameMother(0),
  tMCMotherV0(0),
  tMCMotherP(0),
  tMCptcTypeV0(0),
  tMCptcTypeP(0),
  tMCSameGM(0),
  tMotherPDG(0),
  tpdgcodeV0(0),
  tpdgcodeP(0),
  tKstarGen(0),
  
  
  fHistEventMultiplicity(0),
  fHistCentrality(0),
  fHistVertexDistribution(0),
  fHistSphericity(0),
  fHistSpherocity(0),
  fHistMassK0S(0),  
  fHistFirstNPionTPCdEdx(0),
  fHistFirstPPionTPCdEdx(0),
  fHistSecondTPCdEdx(0),
  fHistSecondMassTOFvsPt3sTPC(0)
{
  cout<<"real constructor"<<endl;

  fPDGMassSecond = TDatabasePDG::Instance()->GetParticle(fPDGcodeSecond)->Mass(); 
  fPDGMassFirst = TDatabasePDG::Instance()->GetParticle(fPDGcodeFirst)->Mass(); 

  // constructor
  DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
  // this chain is created by the analysis manager, so no need to worry about it, 
  // it does its work automatically
  DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
  // you can add more output objects by calling DefineOutput(2, classname::Class())
  // if you add more output objects, make sure to call PostData for all of them, and to
  // make changes to your AddTask macro!
  DefineOutput(2, TTree::Class());     
  DefineOutput(3, TTree::Class()); 

}
//_____________________________________________________________________________
AliAnalysisTaskK0SPFemto::~AliAnalysisTaskK0SPFemto()
{
  // destructor
  if (fOutputContainer){
    delete fOutputContainer;     
    fOutputContainer = 0x0;    
  }
  if(fHistSparseSignal){
    delete fHistSparseSignal;
    fHistSparseSignal = 0x0;
  }
  if(fHistSparseBkg) {
    delete fHistSparseBkg;
    fHistSparseBkg = 0x0;
  }
  if (farrGT)
    delete[] farrGT;
  farrGT=0;

  for (unsigned short i=0; i<fzVertexBins; i++) {
    for (unsigned short j=0; j<fnCentBins; j++) {
      delete fEventColl[i][j];
    }
    delete[] fEventColl[i];
  }

  delete[] fEventColl;
}
  

//_____________________________________________________________________________
void AliAnalysisTaskK0SPFemto::UserCreateOutputObjects()
{
  // create output objects
  //
  // this function is called ONCE at the start of your analysis (RUNTIME)
  // here you ceate the histograms that you want to use 
  //
  // the histograms are in this case added to a tlist, this list is in the end saved
  // to an output event
  //

 
  //file collection
  OpenFile(1);

  fEventColl = new AliAnalysisK0SPEventCollection **[fzVertexBins];   
  for (unsigned short i=0; i<fzVertexBins; i++) {
    fEventColl[i] = new AliAnalysisK0SPEventCollection *[fnCentBins];
    for (unsigned short j=0; j<fnCentBins; j++) {
      fEventColl[i][j] = new AliAnalysisK0SPEventCollection(fnEventsToMix+1, fMaxFirstMult, fMaxSecondMult);
    }
  }

  // Store pointer to global tracks
  farrGT = new Int_t[fTrackBufferSize];
  
  //Define and fill the OutputContainer
  fOutputContainer = new TList();
  fOutputContainer->SetOwner(kTRUE);

  // Create histograms 
  
  fHistEventMultiplicity   = new TH1F( "fHistEventMultiplicity" , "Nb of Events" , 13 , 0.5,13.5);
  fHistEventMultiplicity->GetXaxis()->SetBinLabel(1,"All Events");
  fHistEventMultiplicity->GetXaxis()->SetBinLabel(2,"Events w/PV and PID response");
  fHistEventMultiplicity->GetXaxis()->SetBinLabel(3,"Events w/|Vz|<10cm");
  fHistEventMultiplicity->GetXaxis()->SetBinLabel(4,"Centrality acc");
  fHistEventMultiplicity->GetXaxis()->SetBinLabel(5,"w/o PileUp");
  fHistEventMultiplicity->GetXaxis()->SetBinLabel(6,"Any Events");
  fHistEventMultiplicity->GetXaxis()->SetBinLabel(7,"Central Events");
  fHistEventMultiplicity->GetXaxis()->SetBinLabel(8,"Semi-Central Events");
  fHistEventMultiplicity->GetXaxis()->SetBinLabel(9,"MB Events");
  fHistEventMultiplicity->GetXaxis()->SetBinLabel(10,"kInt7 Events");
  fHistEventMultiplicity->GetXaxis()->SetBinLabel(11,"HM Events");
  fHistEventMultiplicity->GetXaxis()->SetBinLabel(12,"Is Selected Events");
  fOutputContainer->Add(fHistEventMultiplicity);

  // hmult = new TH1I("hmult","Multiplicity distribution (after cuts on event)",30000,-0.5,2999.5);
  // hmult->GetXaxis()->SetTitle("Number of tracklets");
  // fOutputContainer->Add(hmult);

  fHistCentrality = new TH1F("fHistCentrality", "Number of events", 10001, -0.5, 100.5);
  fHistCentrality ->GetXaxis()->SetTitle("Centrality");
  fHistCentrality ->GetYaxis()->SetTitle("Entries");
  fOutputContainer->Add(fHistCentrality);

  fHistVertexDistribution = new TH1F("fHistVertexDistribution", "Primary vertex distribution", 40, -20., 20.);
  fHistVertexDistribution ->GetXaxis()->SetTitle("z_{v} (cm)");
  fHistVertexDistribution ->GetYaxis()->SetTitle("Entries");
  fOutputContainer->Add(fHistVertexDistribution); 

  fHistSphericity = new TH1F("fHistSphericity", "Sphericity Distribution", 40, 0., 1.);
  fHistSphericity ->GetXaxis()->SetTitle("Sphericity");
  fHistSphericity ->GetYaxis()->SetTitle("Entries");
  fOutputContainer->Add(fHistSphericity); 

  fHistSpherocity = new TH1F("fHistSpherocity", "Spherocity Distribution", 40, 0., 1.);
  fHistSpherocity ->GetXaxis()->SetTitle("Spherocity");
  fHistSpherocity ->GetYaxis()->SetTitle("Entries");
  fOutputContainer->Add(fHistSpherocity); 

  fHistMassK0S = new TH1F("fHistMassK0S", "MassK0S Distribution", 800, 0.3, 0.7);
  fHistMassK0S ->GetXaxis()->SetTitle("MassK0S");
  fHistMassK0S ->GetYaxis()->SetTitle("Entries");
  fOutputContainer->Add(fHistMassK0S); 


  fHistFirstNPionTPCdEdx = new TH2F("fHistFirstNPionTPCdEdx", "fHistFirstNPionTPCdEdx", 400, -6.0, 6.0, 500, 0.0, 2000);
  fOutputContainer->Add(fHistFirstNPionTPCdEdx);

  fHistFirstPPionTPCdEdx = new TH2F("fHistFirstPPionTPCdEdx", "fHistFirstPPionTPCdEdx", 400, -6.0, 6.0, 500, 0.0, 2000);
  fOutputContainer->Add(fHistFirstPPionTPCdEdx);

  fHistSecondTPCdEdx = new TH2F("fHistSecondTPCdEdx", "fHistSecondTPCdEdx", 400, -6.0, 6.0, 500, 0.0, 2000);
  fOutputContainer->Add(fHistSecondTPCdEdx);
  fHistSecondMassTOFvsPt3sTPC      = new TH2F("fHistSecondMassTOFvsPt3sTPC", "fHistSecondMassTOFvsPt3sTPC", 200, -2.5, 2.5, 200, -10.0, 10);  
  fOutputContainer->Add(fHistSecondMassTOFvsPt3sTPC);

  fEventCuts.AddQAplotsToList(fOutputContainer);
  PostData(1, fOutputContainer);

  OpenFile(2);
  fHistSparseSignal = new TTree("fHistSparseSignal","fHistSparseSignal");
  /*1 */   fHistSparseSignal->Branch("tSignP",&tSignP,"tSignP/I");
  /*2 */   fHistSparseSignal->Branch("tCentrality",&tCentrality,"tCentrality/F");
  /*3 */   fHistSparseSignal->Branch("tDcaPosV0",&tDcaPosV0,"tDcaPosV0/F");
  /*4 */   fHistSparseSignal->Branch("tDcaNegV0",&tDcaNegV0,"tDcaNegV0/F");
  /*5 */   fHistSparseSignal->Branch("tDCAxyP",&tDCAxyP,"tDCAxyP/F"); 
  /*6 */   fHistSparseSignal->Branch("tDCAzP",&tDCAzP,"tDCAzP/F"); 
  /*7 */   fHistSparseSignal->Branch("tKtpair",&tKtpair,"tKtpair/F"); 
  /*8 */   fHistSparseSignal->Branch("tkStar",&tkStar,"tkStar/F");
  /*9 */   fHistSparseSignal->Branch("tPtV0",&tPtV0,"tPtV0/F");
  /*10*/   fHistSparseSignal->Branch("tPtP",&tPtP,"tPtP/F");
  /*11*/   fHistSparseSignal->Branch("tSphericity",&tSphericity,"tSphericity/F");
  /*12*/   fHistSparseSignal->Branch("tSpherocity",&tSpherocity,"tSpherocity/F");
  /*13*/   fHistSparseSignal->Branch("tInvMassK0s",&tInvMassK0s,"tInvMassK0s/F");
  /*14*/   fHistSparseSignal->Branch("tInvMassLambda",&tInvMassLambda,"tInvMassLambda/F");
  /*15*/   fHistSparseSignal->Branch("tInvMassAntiLambda",&tInvMassAntiLambda,"tInvMassAntiLambda/F");
  /*16*/   fHistSparseSignal->Branch("tCosPointingAngleV0",&tCosPointingAngleV0,"tCosPointingAngleV0/F");
  /*17*/   fHistSparseSignal->Branch("tThetaV0",&tThetaV0,"tThetaV0/F");
  /*18*/   fHistSparseSignal->Branch("tThetaP",&tThetaP,"tThetaP/F");
  /*19*/   fHistSparseSignal->Branch("tPhiV0",&tPhiV0,"tPhiV0/F");
  /*20*/   fHistSparseSignal->Branch("tPhiP",&tPhiP,"tPhiP/F");
  /*21*/   fHistSparseSignal->Branch("tMassTOFP",&tMassTOFP,"tMassTOFP/F");
  if (fIsMC) 
    {
  /*22*/   fHistSparseSignal->Branch("tMCtruepair",&tMCtruepair,"tMCtruepair/I");
  /*23*/   fHistSparseSignal->Branch("tMCSameMother",&tMCSameMother,"tMCSameMother/I");
  /*24*/   fHistSparseSignal->Branch("tMCMotherV0",&tMCMotherV0,"tMCMotherV0/I");
  /*25*/   fHistSparseSignal->Branch("tMCMotherP",&tMCMotherP,"tMCMotherP/I");
  /*26*/   fHistSparseSignal->Branch("tMCptcTypeV0",&tMCptcTypeV0,"tMCptcTypeV0/I");
  /*27*/   fHistSparseSignal->Branch("tMCptcTypeP",&tMCptcTypeP,"tMCptcTypeP/I");
  /*28*/   fHistSparseSignal->Branch("tIsCommonParton",&tIsCommonParton,"tIsCommonParton/O");
  /*29*/   fHistSparseSignal->Branch("tMCSameGM",&tMCSameGM,"tMCSameGM/I");
  /*30*/   fHistSparseSignal->Branch("tMotherPDG",&tMotherPDG,"tMotherPDG/I");
  /*31*/   fHistSparseSignal->Branch("tpdgcodeV0",&tpdgcodeV0,"tpdgcodeV0/I");
  /*32*/   fHistSparseSignal->Branch("tpdgcodeP",&tpdgcodeP,"tpdgcodeP/I");
  /*33*/   fHistSparseSignal->Branch("tKstarGen",&tKstarGen,"tKstarGen/F");
    }
  fHistSparseSignal->SetAutoSave(100000000);
  PostData(2,fHistSparseSignal);


  OpenFile(3);
  fHistSparseBkg = new TTree("fHistSparseBkg","fHistSparseBkg");
  /*1 */   fHistSparseBkg->Branch("tSignP",&tSignP,"tSignP/I");
  /*2 */   fHistSparseBkg->Branch("tCentrality",&tCentrality,"tCentrality/F");
  /*3 */   fHistSparseBkg->Branch("tDcaPosV0",&tDcaPosV0,"tDcaPosV0/F");
  /*4 */   fHistSparseBkg->Branch("tDcaNegV0",&tDcaNegV0,"tDcaNegV0/F");
  /*5 */   fHistSparseBkg->Branch("tDCAxyP",&tDCAxyP,"tDCAxyP/F"); 
  /*6 */   fHistSparseBkg->Branch("tDCAzP",&tDCAzP,"tDCAzP/F"); 
  /*7 */   fHistSparseBkg->Branch("tKtpair",&tKtpair,"tKtpair/F"); 
  /*8 */   fHistSparseBkg->Branch("tkStar",&tkStar,"tkStar/F");
  /*9 */   fHistSparseBkg->Branch("tPtV0",&tPtV0,"tPtV0/F");
  /*10*/   fHistSparseBkg->Branch("tPtP",&tPtP,"tPtP/F");
  /*11*/   fHistSparseBkg->Branch("tSphericity",&tSphericity,"tSphericity/F");
  /*12*/   fHistSparseBkg->Branch("tSpherocity",&tSpherocity,"tSpherocity/F");
  /*13*/   fHistSparseBkg->Branch("tInvMassK0s",&tInvMassK0s,"tInvMassK0s/F");
  /*14*/   fHistSparseBkg->Branch("tInvMassLambda",&tInvMassLambda,"tInvMassLambda/F");
  /*15*/   fHistSparseBkg->Branch("tInvMassAntiLambda",&tInvMassAntiLambda,"tInvMassAntiLambda/F");
  /*16*/   fHistSparseBkg->Branch("tCosPointingAngleV0",&tCosPointingAngleV0,"tCosPointingAngleV0/F");
  /*17*/   fHistSparseBkg->Branch("tThetaV0",&tThetaV0,"tThetaV0/F");
  /*18*/   fHistSparseBkg->Branch("tThetaP",&tThetaP,"tThetaP/F");
  /*19*/   fHistSparseBkg->Branch("tPhiV0",&tPhiV0,"tPhiV0/F");
  /*20*/   fHistSparseBkg->Branch("tPhiP",&tPhiP,"tPhiP/F");
  /*21*/   fHistSparseBkg->Branch("tMassTOFP",&tMassTOFP,"tMassTOFP/F");
  if (fIsMC) 
    {
  /*22*/   fHistSparseBkg->Branch("tMCtruepair",&tMCtruepair,"tMCtruepair/I");
  /*23*/   fHistSparseBkg->Branch("tMCSameMother",&tMCSameMother,"tMCSameMother/I");
  /*24*/   fHistSparseBkg->Branch("tMCMotherV0",&tMCMotherV0,"tMCMotherV0/I");
  /*25*/   fHistSparseBkg->Branch("tMCMotherP",&tMCMotherP,"tMCMotherP/I");
  /*26*/   fHistSparseBkg->Branch("tMCptcTypeV0",&tMCptcTypeV0,"tMCptcTypeV0/I");
  /*27*/   fHistSparseBkg->Branch("tMCptcTypeP",&tMCptcTypeP,"tMCptcTypeP/I");
  /*28*/   fHistSparseBkg->Branch("tIsCommonParton",&tIsCommonParton,"tIsCommonParton/O");
  /*29*/   fHistSparseBkg->Branch("tMCSameGM",&tMCSameGM,"tMCSameGM/I");
  /*30*/   fHistSparseBkg->Branch("tMotherPDG",&tMotherPDG,"tMotherPDG/I");
  /*31*/   fHistSparseBkg->Branch("tpdgcodeV0",&tpdgcodeV0,"tpdgcodeV0/I");
  /*32*/   fHistSparseBkg->Branch("tpdgcodeP",&tpdgcodeP,"tpdgcodeP/I");
  /*33*/   fHistSparseBkg->Branch("tKstarGen",&tKstarGen,"tKstarGen/F");
    }
  fHistSparseBkg->SetAutoSave(100000000);
  PostData(3, fHistSparseBkg );

}
//____________________________________________________________________________
void AliAnalysisTaskK0SPFemto::UserExec(Option_t *) {
  
  // Main loop
  // Called for each event

  // if(Neventi>=1000)
  //   {
  //     PostData(1, fOutputContainer);
  //     PostData(2, fHistSparseSignal );
  //     PostData(3, fHistSparseBkg );
  //     return;
  //   }
  // Neventi+=1;
  // cout<<"Evento numero:"<<Neventi<<endl;

  AliVVertex *vertexmain =0x0;
  //RA// AliCentrality* centrality = 0x0;
  AliMultSelection* centrality = 0x0;
  Double_t lBestPrimaryVtxPos[3] = {-100.0, -100.0, -100.0};

  fHistEventMultiplicity->Fill(1);
  
  //  AliMCEvent   *lMCevent  = 0x0;
  //  AliStack     *lMCstack  = 0x0;
  TClonesArray *arrayMC = 0x0;

  Int_t ntracks = 0;

    
  fAOD = dynamic_cast<AliAODEvent*>( InputEvent() );
  //RA//    cout<<"fAOD: "<<fAOD <<endl;
    
  if (!fAOD) {
    AliWarning("ERROR: AODevent not available \n"); 
    PostData(1, fOutputContainer);
    PostData(2, fHistSparseSignal );
    PostData(3, fHistSparseBkg );
    return;
  }

  if(fHMtrigger == kTRUE){ //modify trigger in event selection    
    fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kHighMultV0);
    //fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kAnyINT);
  }
   
  /// Use the event cut class to apply the required selections
  if (!fEventCuts.AcceptEvent(fAOD)) {   
    PostData(1, fOutputContainer);
    PostData(2, fHistSparseSignal );
    PostData(3, fHistSparseBkg );
    return;
  }
  ntracks = fAOD->GetNumberOfTracks(); 
  
  const AliAODVertex *lPrimaryBestAODVtx = fAOD->GetPrimaryVertex();
  if (!lPrimaryBestAODVtx){
    AliWarning("No prim. vertex in AOD... return!");
    PostData(1, fOutputContainer);
    PostData(2, fHistSparseSignal );
    PostData(3, fHistSparseBkg );
    return;
  }
  vertexmain = (AliVVertex*) lPrimaryBestAODVtx;
  lPrimaryBestAODVtx->GetXYZ( lBestPrimaryVtxPos );

  if (fIsMC) {
    //RA//      Printf("Reading MC truth!!! \n");
    arrayMC = (TClonesArray*) fAOD->GetList()->FindObject(AliAODMCParticle::StdBranchName());

    if (!arrayMC) AliFatal("Error: MC particles branch not found!\n");

  }
    
  // PID object
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  
  UInt_t mask = inputHandler->IsEventSelected();
  /*RA// to see how many events are rejected
     if (!(mask & 0xffffffff)) {
     PostData(1,fOutputContainer );
     return;
     }
  */
  
  fPIDResponse = inputHandler->GetPIDResponse();

  
  if(!fPIDResponse) {
    PostData(1,fOutputContainer );
    PostData(2, fHistSparseSignal );
    PostData(3, fHistSparseBkg );
    AliError("Cannot get pid response");
    return;
  }
  fHistEventMultiplicity->Fill(2);

  if((TMath::Abs(lBestPrimaryVtxPos[2])) > 10.) {
    PostData(1, fOutputContainer);
    PostData(2, fHistSparseSignal );
    PostData(3, fHistSparseBkg );
    return;
  }
  fHistEventMultiplicity->Fill(3);

  //event must not be tagged as pileup
  Bool_t isPileUpSpd=kFALSE;
  isPileUpSpd=fAOD->IsPileupFromSPD();
  if(isPileUpSpd){ 
    PostData(1,fOutputContainer );
    PostData(2, fHistSparseSignal );
    PostData(3, fHistSparseBkg );
    return;
  }
  
  fHistEventMultiplicity->Fill(4);


  Float_t lcentrality = -99.;
  //  AliMultSelection* centrality = 0x0;
  centrality = (AliMultSelection *) fAOD->FindListObject("MultSelection");
  //  cout<<"centrality: "<<centrality<<endl;

 if(fCollidingSystem == "PbPb" ){
    //    lcentrality = centrality->GetCentralityPercentile("V0M"); //FIXME : Centrality in pp and pPb works the same???
    lcentrality = centrality->GetMultiplicityPercentile("V0M"); //FIXME : Centrality in pp and pPb works the same???
    //hmult->Fill(lcentrality);
  }
  else if(fCollidingSystem == "pPb") {//FIXME : I think up AOD have only refmult as mult estimation
    // lcentrality = ((AliAODHeader * )fAODevent->GetHeader())->GetRefMultiplicityComb08(); //-->RIcambiare
    //    lcentrality = ((AliAODHeader * )fAODevent->GetHeader())->GetRefMultiplicity(); //run on phojet
    lcentrality = centrality->GetMultiplicityPercentile("V0A"); //FIXME : Also for pp? Test on kd
    //    cout<<"Centrality: "<<lcentrality<<endl;
    //    hmult->Fill(lcentrality);
    
  }
  else if(fCollidingSystem == "pp") {//FIXME : I think up AOD have only refmult as mult estimation
    // lcentrality = ((AliAODHeader * )fAODevent->GetHeader())->GetRefMultiplicityComb08(); //-->RIcambiare
    //    lcentrality = ((AliAODHeader * )fAODevent->GetHeader())->GetRefMultiplicity(); //run on phojet
    lcentrality = centrality->GetMultiplicityPercentile("V0M"); //FIXME : Also for pp? Test on kd
    //    cout<<"Centrality: "<<lcentrality<<endl;
    //    hmult->Fill(lcentrality);
    
  }

  if ( lcentrality > 199 ){
    //Event didn't pass Event Selections
    PostData(1,fOutputContainer );
    PostData(2, fHistSparseSignal );
    PostData(3, fHistSparseBkg );
    return;
  }
  
  if (lcentrality<fCentrLowLim||lcentrality>=fCentrUpLim){   
    PostData(1,fOutputContainer );
    PostData(2, fHistSparseSignal );
    PostData(3, fHistSparseBkg );
    return;
  }
  
  fHistEventMultiplicity->Fill(5);
  

  // Bool_t isSelectedCentral     = kFALSE;
  // Bool_t isSelectedSemiCentral = kFALSE;
  // Bool_t isSelectedMB          = kFALSE;
  Bool_t isSelectedInt7        = kFALSE;
  Bool_t isSelectedHM          = kFALSE;
  // Bool_t isSelectedAny         = kFALSE;
  Bool_t isSelected            = kFALSE;
   
  
  // isSelectedCentral     = (mask & AliVEvent::kCentral);
  // isSelectedSemiCentral = (mask & AliVEvent::kSemiCentral);
  // isSelectedMB          = (mask & AliVEvent::kMB);
  isSelectedInt7        = (mask & AliVEvent::kINT7);
  isSelectedHM          = (mask & AliVEvent::kHighMultV0);
  // isSelectedAny         = (mask & AliVEvent::kAnyINT);
    
      
  // if(fYear == 2010 && isSelectedMB )
  //   isSelected = kTRUE;
   if( fHMtrigger == kFALSE && isSelectedInt7)
     isSelected = kTRUE;
   else if( fHMtrigger == kTRUE && isSelectedHM)
     isSelected = kTRUE;
   else 
     isSelected = kFALSE;
  
    
  //cout<<isSelectedAny<<" "<<isSelectedInt7<<" "<<isSelectedHM<<endl;

  // if(isSelectedAny)
  //   fHistEventMultiplicity->Fill(6);
  // if(isSelectedCentral)
  //   fHistEventMultiplicity->Fill(7);
  // if(isSelectedSemiCentral)
  //   fHistEventMultiplicity->Fill(8);
  // if(isSelectedMB)
  //   fHistEventMultiplicity->Fill(9);
  if(isSelectedInt7)
    fHistEventMultiplicity->Fill(10);
  if(isSelectedHM)
    fHistEventMultiplicity->Fill(11);
  
  //RA//  cout<<"Trigger mask: "<<fAOD->GetTriggerMask()<<" "<<AliVEvent::kHighMultV0<<endl;
  //RA//  cout<<"Trigger mask: "<<mask<<" "<<AliVEvent::kHighMultV0<<endl;
  //RA// cout<<"Event type  : "<<fAOD->GetEventType()<<endl;
  //RA// FIXME : event selection to be added.. DONE

    
  if(!isSelected){   
    PostData(1, fOutputContainer);
    PostData(2, fHistSparseSignal );
    PostData(3, fHistSparseBkg );
    return;
  }
  

  fHistEventMultiplicity->Fill(12); // is event selected for the analysis
  
  //hmult->Fill(((AliAODHeader * )fAOD->GetHeader())->GetRefMultiplicityComb08()); 

  // cout<<"nTracks: "<<ntracks<<" centrality "<<lcentrality<<endl;
  Double_t fSphericityvalue = CalculateSphericityofEvent(fAOD); 
  Double_t fSpherocityvalue = CalculateSpherocityEvent(fAOD); 
  fHistSphericity->Fill(fSphericityvalue);
  fHistSpherocity->Fill(fSpherocityvalue);
  fHistCentrality->Fill(lcentrality);
  fHistVertexDistribution->Fill(lBestPrimaryVtxPos[2]); 

  const Float_t bfield = (InputEvent())->GetMagneticField();
  int fieldsign;
  if (bfield >=0.) fieldsign = 1;
  else fieldsign = -1;
  // Store the event in the buffer to do mixing
  // ... find vertex... 
  int zBin=0;

  double zStep=2*10/double(fzVertexBins), zStart=-10.;
  
  for (int i=0; i<fzVertexBins; i++) {
    if ((lBestPrimaryVtxPos[2] > zStart+i*zStep) && (lBestPrimaryVtxPos[2] < zStart+(i+1)*zStep)) {
      zBin=i;
      break;
    }
  }
 

  //CENTRALITY!!
  // ... and centrality    //FIXME : find out how centrality wokrs in AOD in pp and pPb
  int centralityBin=0;

  if(fCollidingSystem == "PbPb" ){
    if(lcentrality < 5.) centralityBin=19;      // changed <= with < to be consistent with histogram binning, except last bin 
    else if(lcentrality < 10.) centralityBin=18;
    else if(lcentrality < 15.) centralityBin=17;
    else if(lcentrality < 20.) centralityBin=16;
    else if(lcentrality < 25.) centralityBin=15;
    else if(lcentrality < 30.) centralityBin=14;
    else if(lcentrality < 35.) centralityBin=13;
    else if(lcentrality < 40.) centralityBin=12;
    else if(lcentrality < 45.) centralityBin=11;
    else if(lcentrality < 50.) centralityBin=10;
    else if(lcentrality < 55.) centralityBin=9; 
    else if(lcentrality < 60.) centralityBin=8;
    else if(lcentrality < 65.) centralityBin=7;
    else if(lcentrality < 70.) centralityBin=6;
    else if(lcentrality < 75.) centralityBin=5;
    else if(lcentrality < 80.) centralityBin=4;
    else if(lcentrality < 85.) centralityBin=3;
    else if(lcentrality < 90.) centralityBin=2;   // from here on, wont be filled because the range selected in AddTask is 0-90
    else if(lcentrality < 95.) centralityBin=1;
    else if(lcentrality <= 100.) centralityBin=0;
  }

  else if(fCollidingSystem == "pp"  || fCollidingSystem == "pPb"){ 
    // this should be valid for centrality...
    if(lcentrality < 0.01)       centralityBin=19;  
    else if(lcentrality < 0.1)   centralityBin=18;
    else if(lcentrality < 0.5)   centralityBin=17;
    else if(lcentrality < 1.0)   centralityBin=16;
    else if(lcentrality < 5.0)   centralityBin=15;
    else if(lcentrality < 10.)   centralityBin=14;
    else if(lcentrality < 20.)   centralityBin=13;
    else if(lcentrality < 30.)   centralityBin=12;
    else if(lcentrality < 40.)   centralityBin=11;
    else if(lcentrality < 50.)   centralityBin=10;
    else if(lcentrality < 70.)   centralityBin=9; 
    else if(lcentrality <= 100.) centralityBin=8;
  }
 
  fEventColl[zBin][centralityBin]->FifoShift();
  fEvt = fEventColl[zBin][centralityBin]->fEvt;
 
  //RA//  printf("buffer size: %d\n",fTrackBufferSize);
  //RA//  printf("ntracks: %d\n",ntracks);
  
  for (Int_t igt = 0; igt < fTrackBufferSize; igt++) farrGT[igt] = -1;
  
  AliAODTrack *globaltrack = 0x0;
 
  // Read and store global tracks to retrieve PID information for TPC only tracks
  for (Int_t igt = 0; igt < ntracks; igt++) {
    globaltrack = (AliAODTrack*) fAOD->GetTrack(igt);

    if (!globaltrack) continue; 
    if (globaltrack->GetID()<0 ) continue;
    if (!globaltrack->IsOn(AliAODTrack::kTPCrefit)) continue; // there are such tracks with no TPC clusters

    // Check id is not too big for buffer
    if (globaltrack->GetID()>=fTrackBufferSize) {
      printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n",globaltrack->GetID(),fTrackBufferSize);
      //fHistTrackBufferOverflow->Fill(1);
      continue;
    }

    if ( globaltrack->GetTPCNcls()<=0   ) { // such tracks do not have the TPC refit either, filter map is 2 --> ITS constrained
      //      cout<<" No TPC cl for this global track!!  "<<igt<<endl;
      //      if (!globaltrack->IsOn(AliAODTrack::kTPCrefit)) cout<<" ... and also no tpc refit  "<<globaltrack->GetFilterMap()<<endl;
      continue;
    }
    //cout<<" The array will contain "<<igt<<" , it contains "<<farrGT[globaltrack->GetID()]<<endl; 

    // Warn if we overwrite a track
    if (farrGT[globaltrack->GetID()]>=0) { // Two tracks same id  --> checked that it never happens
      //      cout<<" Array already filled "<<farrGT[globaltrack->GetID()]<<endl;
    } else { 
      farrGT[globaltrack->GetID()] = igt;           // solution adopted in the femto framework
      //    cout<<" Set array, now it contains  "<<farrGT[globaltrack->GetID()]<<endl;
    }
  }
  
  globaltrack = 0x0; 


  //**DICHIARAZIONE VARIABILI
  AliAODTrack *track = 0x0;
  AliVTrack *vtrackg = 0x0;
  AliVTrack *vtrack = 0x0;

  Float_t TPCNCrossedRows=0.;
  Float_t TPCNclsF=0.;

  Bool_t isTOFPIDok = kFALSE;
 
  //  Float_t nsigmaTOFs = 10.;
  Float_t nsigmaTPCs = 10.;
  
  //  Float_t probMis = 0.;
  Double32_t tTOF = 0.;

  Float_t beta  = 0.;
  Float_t gamma = 0.;
  Float_t mass  = 0.;
  Float_t length = 0;
  Float_t ptot = 0;

  Int_t label = 0;
  Int_t PDGcode=0;

  AliPIDResponse::EDetPidStatus statusTOF;
 
  Double_t expectedTimes[AliPID::kSPECIES];

  Float_t dz[2] = {-99.,-99.}; 
  Double_t dzg[2]= {-999.,-999.}; Double_t covarg[3]={-999.,-999.,-999.};
  AliExternalTrackParam etp1;

  Float_t rapiditySecond = 0.;
  Short_t charge = -2;

  int sCount = 0;

  AliReconstructedSecondK0SP::MCSecondOrigin_t mcSecondOrigin = AliReconstructedSecondK0SP::kUnassigned;
  
  Bool_t isP = kFALSE;  // particle
  Bool_t isaP = kFALSE; // anti-particle 

 
  for (Int_t ip = 0; ip < ntracks; ip++)
    {   
      vtrack = fAOD->GetTrack(ip);
      if (!vtrack) continue;
      track = dynamic_cast<AliAODTrack*>(vtrack);
      if(!track) AliFatal("Not a standard AOD");

      if(!track->TestFilterBit(fFilterBit)) continue;
  
      TPCNCrossedRows = track->GetTPCNCrossedRows();
      TPCNclsF = (Float_t)track->GetTPCNclsF();
      if(TPCNCrossedRows<70) continue;
      if(TPCNclsF==0) continue;
      if(TPCNCrossedRows/TPCNclsF<0.8) continue;


      // Get the corresponding global track to use PID --> stored only for global tracks
      // 
      // Check that the array fGTI isn't too small
      // for the track id
      if (-track->GetID()-1 >= fTrackBufferSize) 
  	{
  	  printf ("Exceeding buffer size!!");
  	  continue;
  	}
      if(fFilterBit == 128)
  	vtrackg = fAOD->GetTrack(farrGT[-track->GetID()-1]);
      else
  	vtrackg = track;

      if (!vtrackg) {
  	printf ("No global info! iTrack %d, ID %d\n",ip,track->GetID());
  	continue;
      }
      if (farrGT[-track->GetID()-1]>=ntracks || farrGT[-track->GetID()-1]<0) { /*cout<<"This index is out of range!!"<<farrGT[-track->GetID()-1]<<endl;*/ continue;}
      globaltrack = dynamic_cast<AliAODTrack*>(vtrackg); 
      if(!globaltrack) AliFatal("Not a standard AOD");

      //    cout<<" Filter map for the global track "<<globaltrack->GetFilterMap()<<" "<<globaltrack<<endl;
     
      // IP to PV of tracks
      dz[0] = -99.; dz[1]  = -99.;

      track->GetImpactParameters(&dz[0], &dz[1]);

      //dz[0] = globaltrack->DCA();    // the TPC one should be applied the other biases the CF --> from Maciejs note --> FIXME to be checked 
      //dz[1] = globaltrack->ZAtDCA(); // for those two lines check AliAODTrack.h // FIXME these two lines produce shifted distributions, known problem, asked Mac and Marian. 
      // Btw dont propagate TPC constrained! 
      //    Double_t p[3]; track->GetXYZ(p); // first two elements of p give the same as above, ignore the third, those are original DCA of TPC only tracks (with or w/o constraint?)
      //    cout<<"d xy "<<dz[0]<<"while value is for other "<<p[0]<<endl;
      //    cout<<"d z "<<dz[1]<<"while value is for other "<<p[1]<<endl;
    
      etp1.CopyFromVTrack(vtrackg); 
      etp1.PropagateToDCA(vertexmain,(InputEvent())->GetMagneticField(),100.,dzg,covarg);    

      isP  = kFALSE;
      isaP = kFALSE;

      charge = globaltrack->Charge();

      if (charge>0) {
      	isP  = kTRUE;
      	isaP = kFALSE;
      }
      else{
      	if(charge<0){
      	  isP  = kFALSE;
      	  isaP = kTRUE;
      	}
      }
  
      nsigmaTPCs = fPIDResponse->NumberOfSigmasTPC(globaltrack, (AliPID::EParticleType)AliPID::kProton);     
      if(std::abs(nsigmaTPCs) > 3 ) continue;
      if (TMath::Abs(globaltrack->Eta())> 0.8) continue;
      // nsigmaTOFs = 10.; // be careful with those initialization    
      // probMis = 10.; 
    
      rapiditySecond = 0.5*TMath::Log( (track->E(fPDGMassSecond) + track->Pz()) / (track->E(fPDGMassSecond) - track->Pz() +1.e-13));
    
     
      //TOF
    
      statusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,globaltrack); // this checks kTOFout and kTIMEi https://twiki.cern.ch/twiki/bin/viewauth/ALICE/TOF

      tTOF = 0.0;
      isTOFPIDok = kFALSE;

      if ((statusTOF ==  AliPIDResponse::kDetPidOk)) 
	{  
	  //	  nsigmaTOFs = fPIDResponse->NumberOfSigmasTOF(globaltrack,  (AliPID::EParticleType)AliPID::kProton); 
           
	  isTOFPIDok = kTRUE;
	  //	  probMis = fPIDResponse->GetTOFMismatchProbability(globaltrack);
	  tTOF = globaltrack->GetTOFsignal()-fPIDResponse->GetTOFResponse().GetStartTime(globaltrack->P());
	  globaltrack->GetIntegratedTimes(expectedTimes);
	}
      // else
      // 	{ 
      // 	  probMis = 1.; nsigmaTOFs = 10.; //cout<<"The corresponding global track has no tof pid!"<<endl;
      // 	} 
    
   
      //
      // HERE IS THE PID!
      //
      // 
      Bool_t isMCsecond = kFALSE; 
      Int_t MCptcCodeP2 = -999;
      Int_t MCmumIDP2   = -999;
      Int_t MCmumPDGP2  = -999;
      Int_t MCgrammaPDGP2  = -999;
      Int_t MCgrammaIDP2   = -999;


      AliAODMCParticle *tparticle = 0x0;
      AliAODMCParticle *AncParticle[50]={0};
      Int_t AncPdg[50]={0};
      Int_t AncParticleLabel[50]={0};
      //      Int_t *AncPrticleLabelnew;

      if(fIsMC == kTRUE)
	{
	  label = track->GetLabel();
	  tparticle = (AliAODMCParticle*)arrayMC->At(TMath::Abs(label));
	  PDGcode = tparticle->GetPdgCode();

	  if(TMath::Abs(PDGcode)==fPDGcodeSecond)
	    {
	      isMCsecond= kTRUE;
	      //cout<<"Label: "<<label<<" "<<PDGcode<<endl;
	      /*
	      Int_t mcMotherLabel = tparticle->GetMother();
	      Int_t mcMotherPdg = 0;
	      AliAODMCParticle *mcMother = (AliAODMCParticle*)arrayMC->At(mcMotherLabel);
	    
	      Int_t mcGrandMotherLabel = mcMother->GetMother();
	      Int_t mcGrandMotherPdg = 0; 
	      AliAODMCParticle *mcGrandMother = (AliAODMCParticle*)arrayMC->At(mcGrandMotherLabel);
	    
	      if(mcMotherLabel < 0) {mcMotherPdg = 0;} else {mcMotherPdg = mcMother->GetPdgCode();} 
	      if(mcGrandMotherLabel < 0){mcGrandMotherPdg=0;}else{mcGrandMotherPdg = mcGrandMother->GetPdgCode();}

	      //cout<<"mcMotherlabel: "<<mcMotherLabel<<endl;
	    
	      //  Mum id
	      MCmumIDP2     = mcMotherLabel;
	      MCmumPDGP2    = mcMotherPdg;
	      MCgrammaIDP2  = mcGrandMotherLabel;
	      MCgrammaPDGP2 = mcGrandMotherPdg;
	      */
	      //cout<<"1:------------------------------MCmumIDP2: "<<MCmumIDP2<<"   ----------------------------------MCmumPDGP2: "<<MCmumPDGP2<<endl;
	      
	      Int_t mcMotherLabel = tparticle->GetMother();
	      Int_t mcMotherPdg = 0;
	      Int_t mcGrandMotherLabel = 0;
	      Int_t mcGrandMotherPdg = 0; 
	      

	      AliAODMCParticle *mcMother = (AliAODMCParticle*)arrayMC->At(mcMotherLabel);
	      if(mcMother){
		AliAODMCParticle *mcGrandMother = (AliAODMCParticle*)arrayMC->At(mcGrandMotherLabel);
		mcGrandMotherLabel = mcMother->GetMother();
		mcGrandMotherPdg = 0; 
		    
		//  if (mcMotherLabel < -1) {mcMotherPdg = 0;} else {mcMotherPdg = mcMother->GetPdgCode();} //RAMONA : era questo 02/03/16
		if(mcMotherLabel < 0) {mcMotherPdg = 0;} else {mcMotherPdg = mcMother->GetPdgCode();} 
		if(mcGrandMotherLabel < 0){mcGrandMotherPdg=0;}else{mcGrandMotherPdg = mcGrandMother->GetPdgCode();}
	      }
	      // cout<<"mcMotherlabel: "<<mcMotherLabel<<endl;

	      //Mum id
	      MCmumIDP2     = mcMotherLabel;
	      MCmumPDGP2    = mcMotherPdg;
	      MCgrammaIDP2  = mcGrandMotherLabel;
	      MCgrammaPDGP2 = mcGrandMotherPdg;


	      if (tparticle->IsPhysicalPrimary()) MCptcCodeP2 = 1;
	      else if (tparticle->IsSecondaryFromMaterial()) MCptcCodeP2 = 2;
	      else if (tparticle->IsSecondaryFromWeakDecay()) MCptcCodeP2 = 3;
	      else MCptcCodeP2 = 4;

	      /* try to get info about origin of particles (beginning)*/    
	      AncParticle[0]=tparticle; //assegna solo la particella 0 che è quella di partenza
	      AncPdg[0]= AncParticle[0]->GetPdgCode();
	      AncParticleLabel[0]= AncParticle[0]->GetLabel();

	      for (Int_t i=0; i<49; i++)
		{
		  AncParticleLabel[i+1]=AncParticle[i]->GetMother();//va salvato per ogni k0s selezionata
		  AncParticle[i+1] = static_cast<AliAODMCParticle*>(arrayMC-> At(TMath::Abs(AncParticleLabel[i+1])));  
		  //per i=0 la particella vparticle[1] diventa la madre della vparticle[0]  e avanti così
		  //così si riempie un array di antenate fino alla prima particella generatrice
		  AncPdg[i+1] = AncParticle[i+1]->GetPdgCode(); //salvato per ogni k0s selezionata (quindi aggiungi a feventcoll)

		  if ((AncParticleLabel[i] ==1 || AncParticleLabel[i] ==-1) && (AncPdg[i]==2212) && (AncParticleLabel[i] == AncParticleLabel[i+1] )) 
		    //**boh
		    break;
		}
	      /* try to get info about origin of particles (end)*/
	    }
	}//end of MC loop
      	    
      if(isTOFPIDok && tTOF > 0.) 
	{
	  length = globaltrack->GetIntegratedLength();  // // FIXME length is zero!! from a mail february 2014: this info is not available for AODs 115, use AODs 145
	  ptot = globaltrack->GetTPCmomentum();
      
	  //if (probMis > 0.01) continue; 
	  //if (TMath::Sqrt(nsigmaTOFs*nsigmaTOFs+nsigmaTPCs*nsigmaTPCs)> fnSigmaTPCTOFPIDsecondParticle) continue;   // this cleans the TOF corrected time plot vs p      //FIXME : Original line from mariella

	  beta = length/(tTOF*2.99792457999999984e-02); 
	  //cout<<" rack length  "<<length<<" beta "<<beta<<endl; 
	  gamma = 1/TMath::Sqrt(1 - beta*beta);
	  mass = ptot/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
	  //cout<<"ptc1: "<<length<<" "<<ptot<<" "<<mass<<endl;
	
	}
	else
	{
	  mass=-99.;
	}	
	 
      if(mass-fPDGMassSecond<-0.2 || mass-fPDGMassSecond>0.2) continue; //taglio su massa tof
      
      if(fIsMC == kTRUE)
	fHistSecondTPCdEdx->Fill(globaltrack->Pt()*charge, globaltrack->GetTPCsignal());
      else if(fIsMC == kFALSE)
	fHistSecondTPCdEdx->Fill(globaltrack->GetTPCmomentum()*charge, globaltrack->GetTPCsignal());

      fHistSecondMassTOFvsPt3sTPC->Fill(mass-fPDGMassSecond,charge*track->Pt());    

      //fHistSecondMassTOFvsPt3sTPC3sTOF->Fill(mass-fPDGMassSecond,charge*track->Pt());
	    
      //------------------ Save second particle information
	    
      fEvt->fReconstructedSecond[sCount].sCharge = charge;

      if(fIsMC == kTRUE){
	fEvt->fReconstructedSecond[sCount].sMomentumTruth[0]  = tparticle->Px();
	fEvt->fReconstructedSecond[sCount].sMomentumTruth[1]  = tparticle->Py();
	fEvt->fReconstructedSecond[sCount].sMomentumTruth[2]  = tparticle->Pz();
	fEvt->fReconstructedSecond[sCount].sAncestorParticleLabel = AncParticleLabel;
	fEvt->fReconstructedSecond[sCount].sAncestorPdg = AncPdg;
      }else{
	fEvt->fReconstructedSecond[sCount].sMomentumTruth[0]  = 0.;
	fEvt->fReconstructedSecond[sCount].sMomentumTruth[1]  = 0.;
	fEvt->fReconstructedSecond[sCount].sMomentumTruth[2]  = 0.;
      }

      fEvt->fReconstructedSecond[sCount].sMomentum[0]  = track->Px();
      fEvt->fReconstructedSecond[sCount].sMomentum[1]  = track->Py();
      fEvt->fReconstructedSecond[sCount].sMomentum[2]  = track->Pz();
      fEvt->fReconstructedSecond[sCount].sPt           = track->Pt();
      fEvt->fReconstructedSecond[sCount].sEta          = track->Eta();
      fEvt->fReconstructedSecond[sCount].sPhi          = track->Phi();
      fEvt->fReconstructedSecond[sCount].sTheta        = track->Theta();
	    
      fEvt->fReconstructedSecond[sCount].sRap    = rapiditySecond; 
      fEvt->fReconstructedSecond[sCount].mcSecondOriginType = mcSecondOrigin;
      fEvt->fReconstructedSecond[sCount].isMCptc = isMCsecond;
      fEvt->fReconstructedSecond[sCount].sMCcode = MCptcCodeP2;
      fEvt->fReconstructedSecond[sCount].sPDGcode = PDGcode;
	    
      //  cout<<"2b:------------------------------MCmumIDP2: "<<MCmumIDP2<<"   ----------------------------------MCmumPDGP2: "<<MCmumPDGP2<<endl;
	    
      fEvt->fReconstructedSecond[sCount].sDCAxy = dzg[0];
      fEvt->fReconstructedSecond[sCount].sDCAz  = dzg[1];
      fEvt->fReconstructedSecond[sCount].sMassTOF  = mass;

      fEvt->fReconstructedSecond[sCount].sMCmumIdx = MCmumIDP2;
      fEvt->fReconstructedSecond[sCount].sMCmumPDG  = MCmumPDGP2;
      fEvt->fReconstructedSecond[sCount].sMCgrandmumIdx  = MCgrammaIDP2;
      fEvt->fReconstructedSecond[sCount].sMCgrandmumPDG  = MCgrammaPDGP2;

      if (isP){
	fEvt->fReconstructedSecond[sCount].isP  = kTRUE; 
	fEvt->fReconstructedSecond[sCount].isaP = kFALSE;
      }
      else 
	if (isaP){
	  fEvt->fReconstructedSecond[sCount].isP  = kFALSE;
	  fEvt->fReconstructedSecond[sCount].isaP = kTRUE; 
	}
	    
      fEvt->fReconstructedSecond[sCount].index = TMath::Abs(globaltrack->GetID());

      // Int_t *AncParticleLabelP=0x0;
      // Int_t *AncPdgP=0x0;
      // if(fIsMC && isMCsecond)
      // 	{
      // 	  AncParticleLabelP = fEvt->fReconstructedSecond[sCount].sAncestorParticleLabel;
      // 	  AncPdgP = fEvt->fReconstructedSecond[sCount].sAncestorPdg;
      // 	  for(int j=0;j<50;j++) 
      // 	    cout<<AncParticleLabelP[j]<<" "<<AncPdgP[j]<<endl;
      // 	}



      sCount++;
    
      if (fMaxSecondMult <= sCount){
      	cerr<<"Proton counts exceeded "<<fMaxSecondMult<<"!"<<endl;
      	break;
      }
    }//end track loop

  fEvt->fNumberCandidateSecond = sCount;

  Float_t  ptrackTPCNCrossedRows;
  Float_t  ntrackTPCNCrossedRows;
  Float_t  ptrackTPCNclsF;
  Float_t  ntrackTPCNclsF;
  Float_t rapidityFirst = 0.;
  int fCount = 0;
 
  AliReconstructedFirstK0SP::MCFirstOrigin_t mcFirstOrigin =AliReconstructedFirstK0SP::kUnassigned;
  Int_t fchargeN;
  Int_t fchargeP;

  Int_t iV0s = fAOD->GetNumberOfV0s();
  for(Int_t i=0; i< iV0s; i++)
    {
      AliAODv0* V0 = static_cast<AliAODv0*>(fAOD->GetV0(i));
      if(V0->GetOnFlyStatus())continue;

      AliAODTrack *pTrack=(AliAODTrack *)V0->GetDaughter(0); //0->Positive Daughter
      AliAODTrack *nTrack=(AliAODTrack *)V0->GetDaughter(1); //1->Negative Daughter

      fchargeP = pTrack->Charge();
      fchargeN = nTrack->Charge();

      if (std::abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion)) > 3 ) continue;
      if (std::abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion)) > 3 ) continue;
      if (std::abs(pTrack->Eta())>0.8) continue;
      if (std::abs(nTrack->Eta())>0.8) continue;
      ptrackTPCNCrossedRows = pTrack->GetTPCNCrossedRows();
      ntrackTPCNCrossedRows = nTrack->GetTPCNCrossedRows();
      ptrackTPCNclsF = (Double_t)pTrack->GetTPCNclsF();
      ntrackTPCNclsF = (Double_t)nTrack->GetTPCNclsF();
      if(ptrackTPCNCrossedRows<70 || ntrackTPCNCrossedRows<70) continue;
      if(ptrackTPCNclsF==0 || ntrackTPCNclsF==0) continue;
      if(ptrackTPCNCrossedRows/ptrackTPCNclsF<0.8 || ntrackTPCNCrossedRows/ntrackTPCNclsF<0.8) continue;
      if(V0->DcaV0Daughters()>0.8) continue;
      rapidityFirst=V0->RapK0Short();
      if(TMath::Abs(rapidityFirst)>0.5) continue;
      if(0.4976*(V0->DecayLengthV0(lBestPrimaryVtxPos))/TMath::Sqrt(V0->Ptot2V0())>7*2.6844) continue;

      if(V0->DcaNegToPrimVertex()<0.1) continue;
      if(V0->DcaPosToPrimVertex()<0.1) continue;
      if(V0->CosPointingAngle(lBestPrimaryVtxPos)<0.95) continue;
      if(TMath::Abs(V0->MassLambda()-1.115683)<0.00125) continue;
      if(TMath::Abs(V0->MassAntiLambda()-1.115683)<0.00125) continue;

      fHistMassK0S->Fill(V0->MassK0Short());

      Bool_t isMCfirst  = kFALSE;   
      Int_t MCptcCodeP1 = -999; //1 : primary 2: from weak decay 3: from material
      Int_t MCmumIDP1   = -999;
      Int_t MCmumPDGP1  = -999;
      Int_t MCgrammaPDGP1  = -999;
      Int_t MCgrammaIDP1   = -999;

      AliAODMCParticle *tparticle = 0x0;
      AliAODMCParticle *AncParticle[50]={0};
      Int_t AncPdg[50]={0};
      Int_t AncParticleLabel[50]={0};

      if (fIsMC)   
      	{     
      	  label = pTrack->GetLabel();
      	  tparticle = (AliAODMCParticle*)arrayMC->At(TMath::Abs(label));
      	  Int_t mcMotherLabel1 = tparticle->GetMother();
	 
      	  label = nTrack->GetLabel();
      	  tparticle = (AliAODMCParticle*)arrayMC->At(TMath::Abs(label));
      	  Int_t mcMotherLabel2 = tparticle->GetMother();

      	  tparticle = (AliAODMCParticle*)arrayMC->At(TMath::Abs(mcMotherLabel2));

      	  if(mcMotherLabel1!=mcMotherLabel2) PDGcode=0;
      	  else {
	    PDGcode = tparticle->GetPdgCode();
      	  }


	  if(TMath::Abs(PDGcode)==fPDGcodeFirst)
	    {
	      isMCfirst= kTRUE;
	      //cout<<"Label: "<<label<<" "<<PDGcode<<endl;
	   
	      Int_t mcMotherLabel = tparticle->GetMother();
	      Int_t mcMotherPdg = 0;
	      AliAODMCParticle *mcMother = (AliAODMCParticle*)arrayMC->At(mcMotherLabel);
	    
	      Int_t mcGrandMotherLabel = mcMother->GetMother();
	      Int_t mcGrandMotherPdg = 0; 
	      AliAODMCParticle *mcGrandMother = (AliAODMCParticle*)arrayMC->At(mcGrandMotherLabel);
	    
	      if(mcMotherLabel < 0) {mcMotherPdg = 0;} else {mcMotherPdg = mcMother->GetPdgCode();} 
	      if(mcGrandMotherLabel < 0){mcGrandMotherPdg=0;}else{mcGrandMotherPdg = mcGrandMother->GetPdgCode();}

	      //cout<<"mcMotherlabel: "<<mcMotherLabel<<endl;
	    
	      //  Mum id
	      MCmumIDP1     = mcMotherLabel;
	      MCmumPDGP1    = mcMotherPdg;
	      MCgrammaIDP1  = mcGrandMotherLabel;
	      MCgrammaPDGP1 = mcGrandMotherPdg;

	      //cout<<"1:------------------------------MCmumIDP1: "<<MCmumIDP1<<"   ----------------------------------MCmumPDGP1: "<<MCmumPDGP1<<endl;
	      
	      if (tparticle->IsPhysicalPrimary()){
		MCptcCodeP1 = 1;
	      }
	      else if (tparticle->IsSecondaryFromMaterial()){
		MCptcCodeP1 = 2;
	      }
	      else if (tparticle->IsSecondaryFromWeakDecay()){
		MCptcCodeP1 = 3;
	      }
	      else {
		MCptcCodeP1 = 4;
		//cout<<"-------------------Inside 4 loop !!!!!"<<endl;
	      }

	      /* try to get info about origin of particles (beginning)*/    
	      AncParticle[0]=tparticle; //assegna solo la particella 0 che è quella di partenza
	      AncPdg[0]= AncParticle[0]->GetPdgCode();
	      AncParticleLabel[0]= AncParticle[0]->GetLabel();

	      for (Int_t i=0; i<49; i++)
		{
		  AncParticleLabel[i+1]=AncParticle[i]->GetMother();//va salvato per ogni k0s selezionata
		  AncParticle[i+1] = static_cast<AliAODMCParticle*>(arrayMC-> At(TMath::Abs(AncParticleLabel[i+1])));  
		  //per i=0 la particella particle[1] diventa la madre della particle[0]  e avanti così
		  //così si riempie un array di antenate fino alla prima particella generatrice
		  AncPdg[i+1] = AncParticle[i+1]->GetPdgCode(); //salvato per ogni k0s selezionata (quindi aggiungi a feventcoll)

		  if ((AncParticleLabel[i] ==1 || AncParticleLabel[i] ==-1) && (AncPdg[i]==2212) && (AncParticleLabel[i] == AncParticleLabel[i+1] )) 
		    //**boh
		    break;
		}
	      /* try to get info about origin of particles (end)*/
	    }
      	}//end MC loop


      if(fIsMC == kTRUE)
	{
	  fHistFirstNPionTPCdEdx->Fill(fchargeN*pTrack->Pt(), pTrack->GetTPCsignal());
	  fHistFirstPPionTPCdEdx->Fill(fchargeP*nTrack->Pt(), nTrack->GetTPCsignal());
	}
      else if(fIsMC == kFALSE)
	{
	  fHistFirstNPionTPCdEdx->Fill(fchargeN*pTrack->GetTPCmomentum(), pTrack->GetTPCsignal());
	  fHistFirstPPionTPCdEdx->Fill(fchargeP*nTrack->GetTPCmomentum(), nTrack->GetTPCsignal());
	}

      if(fIsMC == kTRUE){
      	fEvt->fReconstructedFirst[fCount].fMomentumTruth[0]  = tparticle->Px();
      	fEvt->fReconstructedFirst[fCount].fMomentumTruth[1]  = tparticle->Py();
      	fEvt->fReconstructedFirst[fCount].fMomentumTruth[2]  = tparticle->Pz();
	fEvt->fReconstructedFirst[fCount].fAncestorParticleLabel = AncParticleLabel;
        fEvt->fReconstructedFirst[fCount].fAncestorPdg = AncPdg;
      }else{
      	fEvt->fReconstructedFirst[fCount].fMomentumTruth[0]  = 0.;
      	fEvt->fReconstructedFirst[fCount].fMomentumTruth[1]  = 0.;
      	fEvt->fReconstructedFirst[fCount].fMomentumTruth[2]  = 0.;
      }

      fEvt->fReconstructedFirst[fCount].fMomentum[0]  = V0->Px();
      fEvt->fReconstructedFirst[fCount].fMomentum[1]  = V0->Py();
      fEvt->fReconstructedFirst[fCount].fMomentum[2]  = V0->Pz();
      fEvt->fReconstructedFirst[fCount].fPt           = V0->Pt();
      fEvt->fReconstructedFirst[fCount].fEta          = V0->Eta();
      fEvt->fReconstructedFirst[fCount].fPhi          = V0->Phi();
      fEvt->fReconstructedFirst[fCount].fTheta        = V0->Theta();
	    
      fEvt->fReconstructedFirst[fCount].fRap    = rapidityFirst; 
      fEvt->fReconstructedFirst[fCount].mcFirstOriginType = mcFirstOrigin;
      fEvt->fReconstructedFirst[fCount].isMCptc = isMCfirst;
      fEvt->fReconstructedFirst[fCount].fMCcode = MCptcCodeP1;
      fEvt->fReconstructedFirst[fCount].fPDGcode = PDGcode;
	    
      //  cout<<"2b:------------------------------MCmumIDP1: "<<MCmumIDP1<<"   ----------------------------------MCmumPDGP1: "<<MCmumPDGP1<<endl;
	    
      fEvt->fReconstructedFirst[fCount].fDcaPosV0 = V0->DcaPosToPrimVertex();
      fEvt->fReconstructedFirst[fCount].fDcaNegV0 = V0->DcaNegToPrimVertex();
      fEvt->fReconstructedFirst[fCount].fInvMassK0s = V0->MassK0Short();
      fEvt->fReconstructedFirst[fCount].fInvMassLambda = V0->MassLambda();
      fEvt->fReconstructedFirst[fCount].fInvMassAntiLambda = V0->MassAntiLambda();
      fEvt->fReconstructedFirst[fCount].fCosPointingAngle = V0->CosPointingAngle(lBestPrimaryVtxPos);

      fEvt->fReconstructedFirst[fCount].fMCmumIdx = MCmumIDP1;
      fEvt->fReconstructedFirst[fCount].fMCmumPDG  = MCmumPDGP1;
      fEvt->fReconstructedFirst[fCount].fMCgrandmumIdx  = MCgrammaIDP1;
      fEvt->fReconstructedFirst[fCount].fMCgrandmumPDG  = MCgrammaPDGP1;

      fEvt->fReconstructedFirst[fCount].index = TMath::Abs(V0->GetID());
      fEvt->fReconstructedFirst[fCount].indexPosdaughter = TMath::Abs(pTrack->GetID());
      fEvt->fReconstructedFirst[fCount].indexNegdaughter = TMath::Abs(nTrack->GetID());

      fCount++;

      if (fMaxFirstMult <= fCount){
      	cerr<<"K0S counts exceeded "<<fMaxFirstMult<<"!"<<endl;
      	break;
      }
    }//end V0 loop

  fEvt->fNumberCandidateFirst = fCount;

  for (int i=0; i < fEvt->fNumberCandidateFirst; i++) {
    for (int j=0; j<fEvt->fNumberCandidateSecond; j++) {
      if (fEvt->fReconstructedFirst[i].indexPosdaughter == fEvt->fReconstructedSecond[j].index || fEvt->fReconstructedFirst[i].indexNegdaughter == fEvt->fReconstructedSecond[j].index) {
	//cout<<"the track can be both tracks!"<<endl;
	fEvt->fReconstructedFirst[i].doSkipOver = kTRUE;
	fEvt->fReconstructedSecond[j].doSkipOver = kTRUE;
      }
    }
  }
  //--------------------------------------------------------------
  DoPairsh1h2(lcentrality, fieldsign,fSphericityvalue,fSpherocityvalue);  

  // Post output data
  PostData(1, fOutputContainer);
  PostData(2, fHistSparseSignal );
  PostData(3, fHistSparseBkg );
        
}

//----------------------------------------------------------------------------------------------------
  
void AliAnalysisTaskK0SPFemto::DoPairsh1h2 ( const Float_t lcentrality, int fieldsign, const Double_t fSphericityvalue, Double_t fSpherocityvalue ) 
{

  //-----------
  double DcaPosV0  = -999. ;
  double DcaNegV0   = -999. ;  
  double DCAxyP  = -999. ; 
  double DCAzP   = -999. ;  

  double  ptV0 = -999.;
  double  ptP = -999.;

  // Short_t chargeV0 = -999.;
  // Short_t chargeP = -999.;

  //  bool isV0  = kFALSE;
  //  bool isaV0 = kFALSE;
  bool isP  = kFALSE;
  bool isaP = kFALSE;

  //  Int_t  SignV0 = -999;
  Int_t  SignP = -999;

  double phiP  = -999.;
  double phiV0 = -999.;


  double thetaV0 = -999.;
  double thetaP = -999.;
  
  double MassTOFP = -999.; 
  double InvMassK0s = -999.;
  double InvMassLambda = -999.;
  double InvMassAntiLambda = -999.;
  double CosPointingAngleV0 = -999.;


  bool isMC1 = kFALSE;
  bool isMC2 = kFALSE;

  bool isMCvector = kFALSE;
  bool sameMother = kFALSE;
  bool sameGrandMother = kFALSE;
  
  Int_t mcMotherLabelV0 = -999;
  Int_t mcMotherLabelP = -999;

  Int_t mcGrandMotherLabelV0 = -999;
  Int_t mcGrandMotherLabelP = -999;
 
  Int_t typeV0 = -999;
  Int_t typeP = -999;
 
  Int_t mcPDGMotherV0 = 0;
  Int_t mcPDGMotherP = 0;
  //  Int_t mcMotherBin = 0;
 
  Int_t mcPDGGrandMother = 0;
  Int_t mcGrandMotherBin = 0;

  Int_t mcPDGcodeV0 = 0;
  Int_t mcPDGcodeP = 0;
  // Int_t mcPDG1Bin = 0;
  // Int_t mcPDG2Bin = 0;

  int evmultmixed = 0;
  bool multmixedcounted = kFALSE;
 
  double pairKstar   = 0.;
  double pairKstarMC = 0.;
  //  double pairMass  = 0.;
  //  double pairMassE = 0.;
  double pairKt    = 0.;


  for (int i=0; i<fEvt->fNumberCandidateFirst; i++) 
    {
      if (fEvt->fReconstructedFirst[i].doSkipOver) continue;
    
      DcaPosV0         = fEvt->fReconstructedFirst[i].fDcaPosV0;
      DcaNegV0          = fEvt->fReconstructedFirst[i].fDcaNegV0;
      ptV0            = fEvt->fReconstructedFirst[i].fPt;

      //    chargeV0        = fEvt->fReconstructedFirst[i].fCharge;
      isMC1           = fEvt->fReconstructedFirst[i].isMCptc;
      mcMotherLabelV0 = fEvt->fReconstructedFirst[i].fMCmumIdx;
      typeV0          = fEvt->fReconstructedFirst[i].fMCcode;

      mcGrandMotherLabelV0 = fEvt->fReconstructedFirst[i].fMCgrandmumIdx;
      mcPDGcodeV0          = fEvt->fReconstructedFirst[i].fPDGcode;

      mcPDGMotherV0    = fEvt->fReconstructedFirst[i].fMCmumPDG;
      mcPDGGrandMother = fEvt->fReconstructedFirst[i].fMCgrandmumPDG;

      Int_t *AncParticleLabelV0=0x0;
      Int_t *AncPdgV0=0x0;
      if(fIsMC && isMC1)
	{
	  AncParticleLabelV0 = fEvt->fReconstructedFirst[i].fAncestorParticleLabel;
	  AncPdgV0 = fEvt->fReconstructedFirst[i].fAncestorPdg;
	  // for(int mm=0;mm<50;mm++) 
	  //   cout<<AncParticleLabelV0[mm]<<" "<<AncPdgV0[mm]<<endl;
		
	}
      InvMassK0s = fEvt->fReconstructedFirst[i].fInvMassK0s;
      InvMassLambda = fEvt->fReconstructedFirst[i].fInvMassLambda;
      InvMassAntiLambda = fEvt->fReconstructedFirst[i].fInvMassAntiLambda;
      CosPointingAngleV0 =  fEvt->fReconstructedFirst[i].fCosPointingAngle;
      thetaV0 = fEvt->fReconstructedFirst[i].fTheta;
      phiV0 = fEvt->fReconstructedFirst[i].fPhi;

      for (int eventNumber=0; eventNumber<fnEventsToMix+1; eventNumber++) 
	{ 
	  if (!multmixedcounted && eventNumber!=0 && ((fEvt+eventNumber)->fNumberCandidateSecond)!=0) evmultmixed++; 
      
	  for (int j=0; j<(fEvt+eventNumber)->fNumberCandidateSecond; j++) 
	    {
	      if ((fEvt+eventNumber)->fReconstructedSecond[j].doSkipOver) continue;
         
	      DCAxyP         = (fEvt+eventNumber)->fReconstructedSecond[j].sDCAxy;
	      DCAzP          = (fEvt+eventNumber)->fReconstructedSecond[j].sDCAz;
	      ptP            = (fEvt+eventNumber)->fReconstructedSecond[j].sPt;

	      thetaP  = (fEvt+eventNumber)->fReconstructedSecond[j].sTheta;
	      phiP  = (fEvt+eventNumber)->fReconstructedSecond[j].sPhi;
	      MassTOFP  = (fEvt+eventNumber)->fReconstructedSecond[j].sMassTOF;

	      //chargeP        = (fEvt+eventNumber)->fReconstructedSecond[j].sCharge;
	      isMC2           = (fEvt+eventNumber)->fReconstructedSecond[j].isMCptc;
	      mcMotherLabelP = (fEvt+eventNumber)->fReconstructedSecond[j].sMCmumIdx;
	      typeP          = (fEvt+eventNumber)->fReconstructedSecond[j].sMCcode;

	      mcGrandMotherLabelP = (fEvt+eventNumber)->fReconstructedSecond[j].sMCgrandmumIdx;
	      mcPDGcodeP     = (fEvt+eventNumber)->fReconstructedSecond[j].sPDGcode;

	      isP  = (fEvt+eventNumber)->fReconstructedSecond[j].isP;
	      isaP = (fEvt+eventNumber)->fReconstructedSecond[j].isaP;

	      mcPDGMotherP = (fEvt+eventNumber)->fReconstructedSecond[j].sMCmumPDG;


	      Int_t *AncParticleLabelP=0x0;
	      //	      Int_t *AncPdgP=0x0;
	      if(fIsMC==kTRUE && isMC2==kTRUE)
		{
		  AncParticleLabelP = (fEvt+eventNumber)->fReconstructedSecond[j].sAncestorParticleLabel;
		  //		  AncPdgP = (fEvt+eventNumber)->fReconstructedSecond[j].sAncestorPdg;
		  // for(int mm=0;mm<50;mm++) 
		  //   cout<<AncParticleLabelP[mm]<<" "<<AncPdgP[mm]<<endl;
		}

	      if(isP) SignP = 1;
	      else if (isaP) SignP = -1;

	      if(isMC1 && isMC2) isMCvector = kTRUE;
	      else isMCvector = kFALSE;
	
	      if(mcMotherLabelV0 == mcMotherLabelP && mcMotherLabelV0!=-999)sameMother = kTRUE;
	      else sameMother = kFALSE;
		
	      
	      if(mcGrandMotherLabelV0 == mcGrandMotherLabelP){sameGrandMother = kTRUE;}
	    
	      //  cout<<"GM-------------------------------- chargeV0: "<<chargeV0<<" ----------- chargeP: "<<chargeP<<" ---------------------> "<<mcPDGGrandMother<<endl;
	    
	      if(TMath::Abs(mcPDGGrandMother)>= 1 && TMath::Abs(mcPDGGrandMother)<= 6)  mcGrandMotherBin = 1;  //quark
	      else if(TMath::Abs(mcPDGGrandMother)==2212)  mcGrandMotherBin = 2;  //p
	      else if(TMath::Abs(mcPDGGrandMother)==21)  mcGrandMotherBin = 3;    //g
	      else if(TMath::Abs(mcPDGGrandMother)> 400 && TMath::Abs(mcPDGGrandMother)< 500 )  mcGrandMotherBin = 5;  //D meson
	      else if(mcPDGGrandMother!=0) {
		// cout<<"--------------------------------------------------------------------------------> "<<mcPDGMother<<endl;
		mcGrandMotherBin  = 6; 
	      }//
	      
	      // else if(!sameMother)
	      //   mcMotherBin = 15; 
	
	      //cout<<"----------------- outside: "<<mcGrandMotherBin<<endl;


	      //**QUI
	      Bool_t IsCommonParton=kFALSE;
	      if(fIsMC && isMCvector)
		{
		  for (Int_t ii=1; ii<50; ii++)
		    { //I start from one since last element cannot be a parton but is a hadron
		      if (IsCommonParton==kTRUE) break;
		      for (Int_t jj=1; jj<50; jj++)
		      	{//boh
		      	  if ((AncParticleLabelV0[ii] == AncParticleLabelP[jj] ) &&  AncParticleLabelP[jj]!=0 && ( TMath::Abs(AncPdgV0[ii]) <=8 ||  TMath::Abs(AncPdgV0[ii]) ==21)) 
		      	    { //both Xi and Trigger particle have a common ancestor which has to be a quark or a gluon-> therefore te cascade comes form the jet defined by the trigger particle
		      	      IsCommonParton =kTRUE;
		      	      break;
		      	    }
		      	}
		    }
		}
	      //cout<<IsCommonParton<<endl;

	      //Calculate k* for the pair
	      pairKstar = CalculateKstar(fEvt->fReconstructedFirst[i].fMomentum, (fEvt+eventNumber)->fReconstructedSecond[j].sMomentum,fPDGMassFirst, fPDGMassSecond);
	      //mc kstar
	      pairKstarMC = CalculateKstar(fEvt->fReconstructedFirst[i].fMomentumTruth, (fEvt+eventNumber)->fReconstructedSecond[j].sMomentumTruth,fPDGMassFirst, fPDGMassSecond);
		
	      // //Invariant Mass of the pair
	      //pairMass  = CalculateMass(fEvt->fReconstructedFirst[i].fMomentum, (fEvt+eventNumber)->fReconstructedSecond[j].sMomentum,fPDGMassFirst, fPDGMassSecond);

	      //pairMassE  = CalculateMass(fEvt->fReconstructedFirst[i].fMomentum, (fEvt+eventNumber)->fReconstructedSecond[j].sMomentum,5.11e-4, 5.11e-4);
	      // //Kt
	      pairKt = pow(fEvt->fReconstructedFirst[i].fMomentum[0] + (fEvt+eventNumber)->fReconstructedSecond[j].sMomentum[0],2.);
	      pairKt+= pow(fEvt->fReconstructedFirst[i].fMomentum[1] + (fEvt+eventNumber)->fReconstructedSecond[j].sMomentum[1],2.);
	      pairKt = sqrt(pairKt)/2.;

	      if (eventNumber==0) {//Same event pair histogramming

		tSignP             = SignP;                                 
		tCentrality        = lcentrality;                      
		tDcaPosV0          = DcaPosV0;
		tDcaNegV0          = DcaNegV0;
		tDCAxyP            = DCAxyP;                     
		tDCAzP             = DCAzP;                     
		tKtpair             = pairKt;                  
		tkStar             = pairKstar;                 
		tPtV0              = ptV0;                
		tPtP               = ptP;                

		tInvMassK0s = InvMassK0s;
		tInvMassLambda = InvMassLambda;
		tInvMassAntiLambda = InvMassAntiLambda;
		tCosPointingAngleV0 = CosPointingAngleV0;
		tThetaV0 = thetaV0;
		tThetaP = thetaP;
		tPhiV0 = phiV0;
		tPhiP = phiP;
		tMassTOFP = MassTOFP;
		// tDEta              = deta; //**SERVONO?                 
		// tDPhiStar          = dphis;                    
		// tDPhi              = dphi;                
		//tMassPair          = pairMass;                    
		tSphericity        = fSphericityvalue;
		tSpherocity        = fSpherocityvalue;
		// tGammaCoversionMass= pairMassE;                              
		// tDTheta            = dtheta;                                     
  
		if(fIsMC == kTRUE){
		  tMCtruepair    =  isMCvector;
		  tMCSameMother  =  sameMother;
		  tMCMotherV0    =  mcPDGMotherV0;//mcMotherBin;
		  tMCMotherP    =  mcPDGMotherP;//mcMotherBin;
		  tMCptcTypeV0   =  typeV0 ;
		  tMCptcTypeP   =  typeP ;
		  tIsCommonParton = IsCommonParton;
		  tMCSameGM   =  sameGrandMother;
		  tMotherPDG   =  mcGrandMotherBin;
		  tpdgcodeV0    =  mcPDGcodeV0;//mcPDG1Bin;
		  tpdgcodeP   =  mcPDGcodeP;//mcPDG2Bin; 
		  tKstarGen      =  pairKstarMC;            
		}
	    
		fHistSparseSignal->Fill();  
	      }

	      else {//Mixed-event pair histogramming
	    
		tSignP             = SignP;
		tCentrality         = lcentrality;                      
		tDcaPosV0            = DcaPosV0;                   
		tDcaNegV0             = DcaNegV0;                     
		tDCAxyP            = DCAxyP;                     
		tDCAzP             = DCAzP;                     
		tKtpair             = pairKt;                  
		tkStar              = pairKstar;                 
		tPtV0               = ptV0;                
		tPtP               = ptP;                

		tInvMassK0s = InvMassK0s;
		tInvMassLambda = InvMassLambda;
		tInvMassAntiLambda = InvMassAntiLambda;
		tCosPointingAngleV0 = CosPointingAngleV0;
		tThetaV0 = thetaV0;
		tThetaP = thetaP;
		tPhiV0 = phiV0;
		tPhiP = phiP;
		tMassTOFP = MassTOFP;

		// tDEta               = deta;                 
		// tDPhiStar           = dphis;                    
		// tDPhi               = dphi;                
		//tMassPair           = pairMass;                    
		tSphericity         = fSphericityvalue;                      
		tSpherocity         = fSpherocityvalue;                      
		//tGammaCoversionMass = pairMassE;                              
		//tDTheta             = dtheta;                                     

		if(fIsMC == kTRUE){
		  tMCtruepair    =  isMCvector;
		  tMCSameMother  =  sameMother;
		  tMCMotherV0    =  mcPDGMotherV0;//mcMotherBin;
		  tMCMotherP    =  mcPDGMotherP;//mcMotherBin;
		  tMCptcTypeV0   =  typeV0 ;
		  tMCptcTypeP   =  typeP ;
		  tIsCommonParton = IsCommonParton;
		  tMCSameGM   =  sameGrandMother;
		  tMotherPDG   =  mcGrandMotherBin;
		  tpdgcodeV0    =  mcPDGcodeV0;//mcPDG1Bin;
		  tpdgcodeP   =  mcPDGcodeP;//mcPDG2Bin; 
		  tKstarGen      =  pairKstarMC;            
		}
	    
		fHistSparseBkg->Fill();  
	    
	      } //mixed
	
	    } // second part

	}//end event loop

      if (evmultmixed!=0) multmixedcounted = kTRUE;
    
    } // first part
  
  //if(multmixedcounted)  fHistMultiplicityOfMixedEvent->Fill(evmultmixed);
  
}

// //----------------------------------------------------------------------------------------------
// //void AliAnalysisTaskK0SPFemto::DoPairshh (const Float_t lcentrality, int fieldsign) {
// void AliAnalysisTaskK0SPFemto::DoPairshh (const Int_t lcentrality, int fieldsign, const Double_t fSphericityvalue) {
//   return;
// }

// //-----------------------------------------------------------------------------------------------

double AliAnalysisTaskK0SPFemto::CalculateKstar(double momentum1[3], double momentum2[3], double mass1, double mass2) { // Jai S

  // Calculate k* for any pair of particles, regardless of whether the

  // particles have the same mass.

  double kstar = 0.;

  double e1 = 0.;

  double e2 = 0.;

  for(int i = 0; i < 3; i++){

    kstar -= pow(momentum1[i]-momentum2[i],2);

    e1 += pow(momentum1[i],2);

    e2 += pow(momentum2[i],2);

  }

  e1 += pow(mass1,2);

  e1 = sqrt(e1);

  e2 += pow(mass2,2);

  e2 = sqrt(e2);

   

  kstar += pow(e1-e2,2);

 

  double totalMomentumSquared = 0;

  for(int i = 0; i < 3; i++){

    totalMomentumSquared -= pow(momentum1[i]+momentum2[i],2);

  }

  totalMomentumSquared += pow(e1+e2,2);

  kstar -= pow((pow(mass1,2)-pow(mass2,2)),2)/totalMomentumSquared;

  kstar *= -1.;

  kstar = sqrt(kstar); //At this point, we've actually calculated Qinv

  kstar *= 0.5; // kstar is 0.5*Qinv

  return kstar;
}

//-----------------------------------------------------------------------------------------------
// double AliAnalysisTaskK0SPFemto::CalculateMass(double momentum1[3], double momentum2[3], double mass1, double mass2) { // Jai S

//   // Calculate Invariant Mass 
//   TLorentzVector  vP1,vP2,vSum;
  
//   vP1.SetXYZM(momentum1[0],momentum1[1],momentum1[2],mass1); 
//   vP2.SetXYZM(momentum2[0],momentum2[1],momentum2[2],mass2);       
//   vSum=vP1+vP2;

//   double mass = vSum.M();
 
//   return mass;
// }

// //-----------------------------------------------------------------------------------------------
// double AliAnalysisTaskK0SPFemto::CalculateDphiSatR12m(Short_t chg1, Short_t chg2, Int_t magSign, Double_t ptv1, Double_t ptv2, Double_t phi1, Double_t phi2) { // AliFemto framework AliFemtoUser/AliFemtoPairCutRadialDistance.cxx + Dhevan not consistent?
//   /*
//   double rad = 1.2;
//   double afsi0b = 0.07510020733*chg1*magSign*rad/ptv1; // 0.075 = 0.3=e in H-L units*0.5=B/2 calculation on notebook
//   double afsi1b = 0.07510020733*chg2*magSign*rad/ptv2;
//   if (fabs(afsi0b) >=1.) return 9999.; // angle is pi/2 or not defined --> dont cut
//   if (fabs(afsi1b) >=1.) return 9999.; // MN modified these two lines returning 9999 and not kTRUE
//   //  Double_t dps = phi2 - phi1 -TMath::ASin(afsi1b) + TMath::ASin(afsi0b);
//   //  dps = TVector2::Phi_mpi_pi(dps);
  
//   double phi1bis =0.;
//   double phi2bis =0.;
//   phi1bis = phi1-TMath::ASin(afsi0b); 
//   if(phi1bis > 2*PI) phi1bis -= 2*PI;
//   if(phi1bis < 0) phi1bis += 2*PI;
//   phi2bis = phi2 - TMath::ASin(afsi1b);
//   if(phi2bis > 2*PI) phi2bis -= 2*PI;
//   if(phi2bis < 0) phi2bis += 2*PI;
//   double deltaphi = phi2bis - phi1bis;
//   if(deltaphi > PI) deltaphi -= PI;
//   if(deltaphi < -PI) deltaphi += PI;
//   return deltaphi;//dps;
//   */
//   //  cout<<" Dphi "<<dps<<" Dhevan "<<deltaphi<<endl;


//   //from mariella
//   //analitical funcion
  
//   double rad = 1.2;
  
//   double afsi1b = 0.075*chg1*magSign*rad/ptv1; // 0.07510020733 = - 0.3 (= e in Heaviside-Lorentz units) *0.5 (= B in T) /2 (see later for the -), pT in GeV/c
//   double afsi2b = 0.075*chg2*magSign*rad/ptv2;
  
//   if (fabs(afsi1b) >=1.) return 9999.; // angle is pi/2 or not defined --> dont cut
//   if (fabs(afsi2b) >=1.) return 9999.; // MN modified these two lines returning 9999 and not kTRUE
//   double dps = phi2 - phi1  + TMath::ASin(afsi1b) -TMath::ASin(afsi2b); // - sign of e is outside Mariella
  
//   dps = TVector2::Phi_mpi_pi(dps);
  
//   return dps;

// }


// //-----------------------------------------------------------------------------------------------
// double AliAnalysisTaskK0SPFemto::CalculateDPhiStar(Short_t chg1, Short_t chg2, Int_t magSign, Double_t ptv1, Double_t ptv2, Double_t phi1, Double_t phi2,Double_t rad) { //AliFemtoUser/AliFemtoPairCutDetaDphi.h

//   const Double_t unit_factor = 0.299792458 / 2.0;
//   const Double_t b_field = 0.5006670488586 * magSign;

//   Double_t  shift1 = TMath::ASin(unit_factor * chg1 * b_field * rad / ptv1);
//   Double_t  shift2 = TMath::ASin(unit_factor * chg2 * b_field * rad / ptv2);

//   double dps = (phi1 + shift1) - (phi2 + shift2);
  
//   //  dps = TVector2::Phi_mpi_pi(dps); //to be checked

//   return dps; //deltaphi;
  
// }
// //_______________________________________________________________

// Double_t AliAnalysisTaskK0SPFemto::CalculateDeltaEta( Double_t eta1, Double_t eta2 )  {   //AliFemtoUser/AliFemtoPairCutDetaDphi.h
//   const double deta = eta2 - eta1;
//   return deta;
// }
// //_______________________________________________________________
// Double_t AliAnalysisTaskK0SPFemto::CalculateDeltaTheta( Double_t theta1, Double_t theta2 )  {  
//   const double dtheta = theta2 - theta1;
//   return dtheta;
// }

// //-----------------------------------------------------------------------------------------------
// double AliAnalysisTaskK0SPFemto::CalculateDphiSatR12m(Double_t pos1SftR125[3], Double_t pos2SftR125[3]) { // Hans B

//   // Returns delta phi star at R = 1.2 m

//   const Float_t distSft = TMath::Sqrt(TMath::Power(pos1SftR125[0] - pos2SftR125[0],2)
// 				      + TMath::Power(pos1SftR125[1] - pos2SftR125[1],2));
//   return 2.0 * TMath::ATan(distSft/2./(125.)); 

// }

// //-----------------------------------------------------------------------------------------------

// void AliAnalysisTaskK0SPFemto::SetSftPosR125(AliVTrack *track, const Float_t bfield, Double_t priVtx[3], Double_t posSftR125[3] ) {  // Hans B
//   // Sets the spatial position of the track at the radius R=1.25m in the shifted coordinate system
  
  
//   // Initialize the array to something indicating there was no propagation
//   posSftR125[0]=-9999.;
//   posSftR125[1]=-9999.;
//   posSftR125[2]=-9999.;
//   // Make a copy of the track to not change parameters of the track
//   AliExternalTrackParam etp;
//   etp.CopyFromVTrack(track);
  
//   // The global position of the track
//   Double_t xyz[3]={-9999.,-9999.,-9999.};  

//   // The radius we want to propagate to, squared, for faster code
//   const Float_t rSquared = 125.*125.;


//   // Propagation is done in local x of the track
//   for (Float_t x = 58.;x<247.;x+=1.){
//     // Starts at 83 / Sqrt(2) and goes outwards. 85/Sqrt(2) is the smallest local x
//     // for global radius 85 cm. x = 245 is the outer radial limit of the TPC when
//     // the track is straight, i.e. has inifinite pt and doesn't get bent. 
//     // If the track's momentum is smaller than infinite, it will develop a y-component,
//     // which adds to the global radius
//     // We don't change the propagation steps to not mess up things!

//     // Stop if the propagation was not succesful. This can happen for low pt tracks
//     // that don't reach outer radii
//     if (!etp.PropagateTo(x,bfield)) { //cout<<"propagation failed!! and etss is "<<EtaS(posSftR125)<<endl; 
//       break;
//     }
//     etp.GetXYZ(xyz); // GetXYZ returns global coordinates

//     // Calculate the shifted radius we are at, squared. 
//     // Compare squared radii for faster code
//     Float_t shiftedRadiusSquared = (xyz[0]-priVtx[0])*(xyz[0]-priVtx[0])
//       + (xyz[1]-priVtx[1])*(xyz[1]-priVtx[1]);

//     // Roughly reached the radius we want
//     if(shiftedRadiusSquared > rSquared){
      
//       // Bigger loop has bad precision, we're nearly one centimeter too far, 
//       // go back in small steps.
//       while (shiftedRadiusSquared>rSquared) {
//         // Propagate a mm inwards
//         x-=.1;
//         if (!etp.PropagateTo(x,bfield)){
//           // Propagation failed but we're already with a
//           // cm precision at R=1.25m so we only break the 
//           // inner loop
//           //cout<<"propagation failed!! and etss is "<<EtaS(posSftR125)<<endl; 

//           break;
//         }
//         // Get the global position
//         etp.GetXYZ(xyz);
//         // Calculate shifted radius, squared
//         shiftedRadiusSquared = (xyz[0]-priVtx[0])*(xyz[0]-priVtx[0])
// 	  + (xyz[1]-priVtx[1])*(xyz[1]-priVtx[1]);
//       }
//       // We reached R=1.25m with a precission of a cm to a mm,
//       // set the spatial position
//       posSftR125[0]=xyz[0]-priVtx[0];
//       posSftR125[1]=xyz[1]-priVtx[1];
//       posSftR125[2]=xyz[2]-priVtx[2];
//       //cout<<" Pos 125 cm in function end "<<posSftR125[0]<<" "<<posSftR125[1]<<" "<<posSftR125[2]<<endl;
//       // Done
//       return;
//     } // End of if roughly reached radius
//   } // End of coarse propagation loop
// }

// //----------------------------------------------------------------------------------------------
// Double_t AliAnalysisTaskK0SPFemto::ThetaS( Double_t posSftR125[3] ) const { // Hans B
//   // Returns the longitudinal angle of the particles propagated
//   // position at R=1.25m. See
//   // https://edms.cern.ch/file/406391/2/ALICE-INT-2003-038.pdf
//   // for the ALICE coordinate system. Theta is zero at positive z,
//   // pi/2 at z = 0 aka the xy plane and pi at negative z 

//   // R^    ^  
//   //  |   /
//   //  |?'/
//   //  | / ?
//   //  |/____>z
//   // 
//   // Let's compute ?' and ? = pi/2 - ?'
//   // where ?' can even be and should 
//   // sometimes be negative
//   // tan(?') = z/R
//   // ?' = arctan(z/R)
//   // ? = pi/2 - ?'
//   //   = pi/2 - arctan(z/R)
//   // Note that in the doc above theta
//   // is calculated as arccos(z/sqrt(x^2+y^2+z^2))

//   // Array of positions is 85,105,125,..cm,
//   // we take the z position at R=1.25m
//   // return TMath::Pi()/2. - TMath::ATan(fXshifted[2][2]/125.);
//   return TMath::Pi()/2. - TMath::ATan(posSftR125[2]/125.); // ok here R is really there --> transverse plane 
// }
// //_______________________________________________________________
// Double_t AliAnalysisTaskK0SPFemto::EtaS( Double_t posSftR125[3] ) const {  // Hans B
//   // Returns the corresponding eta of a pri. part. 
//   // with this particles pos at R=1.25m

//   // http://en.wikipedia.org/wiki/Pseudorapidity
//   // ? = -ln[ tan(?/2)]
//   // printf("z: %+04.0f, thetaS %+03.2f etaS %+1.2f\n"
//   // ,fXshifted[2][2],ThetaS(),-TMath::Log( TMath::Tan(ThetaS()/2.) ));
//   return -TMath::Log( TMath::Tan(ThetaS(posSftR125 )/2.) );
// }



// //_________________________________________________________________
Double_t AliAnalysisTaskK0SPFemto::CalculateSphericityofEvent(AliAODEvent *aodEvent)
{ //from Oliver
  Double_t Pt_tot = 0.; //total Pt of all protons and v0s in the event

  Double_t S00 = 0.; //Elements of the sphericity matrix
  Double_t S11 = 0.;
  Double_t S10 = 0.;
  
  Int_t NumOfTracks = aodEvent->GetNumberOfTracks();
  if(NumOfTracks<3) return -9999.;//if already at this point not enough tracks are in the event -> return
  

  Int_t NTracks = 0;
  for(Int_t iTrack=0;iTrack<NumOfTracks;iTrack++)
    {
      AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(iTrack));
      
      if(!aodtrack->TestFilterBit(128)) continue;
      
      Double_t Pt = aodtrack->Pt();
      //Double_t Phi = aodtrack->Phi();
      Double_t Px = aodtrack->Px();
      Double_t Py = aodtrack->Py();
      Double_t eta = aodtrack->Eta();
      

      if(!(eta>-0.8 && eta<0.8)) continue;
      if(Pt<0.5) continue;

      Pt_tot += Pt;
      
      S00 += Px*Px/Pt;
      S11 += Py*Py/Pt;
      S10 += Px*Py/Pt;
      NTracks++;
    }
  
  if(NTracks<3) return -9999.;//new flag: check
  
  //normalize to total Pt to obtain a linear form:
  if(Pt_tot == 0.) return -9999.;
  S00 /= Pt_tot;
  S11 /= Pt_tot;
  S10 /= Pt_tot;
  
  //Calculate the trace of the sphericity matrix:
  Double_t T = S00+S11;
  
  //Calculate the determinant of the sphericity matrix:
  Double_t D = S00*S11 - S10*S10;//S10 = S01

  //Calculate the eigenvalues of the sphericity matrix:
  Double_t lambda1 = 0.5*(T + TMath::Sqrt(T*T - 4.*D));
  Double_t lambda2 = 0.5*(T - TMath::Sqrt(T*T - 4.*D));
  
  if((lambda1 + lambda2) == 0.) return -9999.;
  
  Double_t ST = -1.;
  
  if(lambda2>lambda1)
    {
      ST = 2.*lambda1/(lambda1+lambda2);
    }
  else
    {
      ST = 2.*lambda2/(lambda1+lambda2);
    }

  return ST;
}

double AliAnalysisTaskK0SPFemto::CalculateSpherocityEvent(AliAODEvent *evt) {
  float pFull = 0.f;
  float Spherocity = 2.f;

  const float pi = TMath::Pi();

  float pTtot = 0.f;

  std::vector<float> pXVec;
  std::vector<float> pYVec;

  int numOfTracks = evt->GetNumberOfTracks();
  if (numOfTracks < 3)
    return -9999.;

  for (int iTrack = 0; iTrack < numOfTracks; iTrack++) {
    AliAODTrack *track = dynamic_cast<AliAODTrack *>(evt->GetTrack(iTrack));
    if (!track->TestFilterBit(96))
      continue;
    double pt = track->Pt();
    if (TMath::Abs(pt) < 0.5 || TMath::Abs(track->Eta()) > 0.8) {
      continue;
    }
    pTtot += pt;
    pXVec.push_back(track->Px());
    pYVec.push_back(track->Py());
  }

  if (pTtot == 0.f)
    return -9999.;

  const float OneOverPtTotal = 1.f / pTtot;

  float numerator = 0.f;
  float phiparam = 0.f;
  float nx = 0.f;
  float ny = 0.f;

  for (int i = 0; i < 360 / 0.1; ++i) {
    numerator = 0.f;
    phiparam = (pi * i * 0.1 / 180);  // parametrization of the angle
    nx = TMath::Cos(phiparam);  // x component of an unitary vector n
    ny = TMath::Sin(phiparam);  // y component of an unitary vector n

    for (size_t itTrack = 0; itTrack < pXVec.size(); ++itTrack) {
      numerator += TMath::Abs(ny * pXVec[itTrack] - nx * pYVec[itTrack]);  // product between p
      // proyection in XY plane and
      // the unitary vector
    }
    pFull = std::pow((numerator * OneOverPtTotal), 2);

    if (pFull < Spherocity)  // maximization of pFull
      {
	Spherocity = pFull;
      }
  }
  return ((Spherocity) * pi * pi) / 4.0;
}

// //--------------------------------------------------- Methods From AliFemtoESDTrackCut.cxx
// bool AliAnalysisTaskK0SPFemto::IsElectron(float nsigmaTPCE, float nsigmaTPCPi,float nsigmaTPCK, float nsigmaTPCP)
// {
//   //  if(TMath::Abs(nsigmaTPCE)<3 && TMath::Abs(nsigmaTPCPi)>3 && TMath::Abs(nsigmaTPCK)>3 && TMath::Abs(nsigmaTPCP)>3)
//   if(TMath::Abs(nsigmaTPCE)<3)
//     return false;
//   else
//     return true;
// }

// //----------------------------------------------------------

// bool AliAnalysisTaskK0SPFemto::IsPionNSigma(double mom, float nsigmaTPCPi, float nsigmaTOFPi)
// {
//   //sligly changed w.r.t. the original
  
//   return false;

//   // if(mom<0.65){   
//   //   //      if(nsigmaTOFPi<-999.)
//   //   if(nsigmaTOFPi==10)
//   //     {
//   // //use TPC only
//   // if(mom<0.35 && TMath::Abs(nsigmaTPCPi)<3.0) return true;
//   // else if(mom<0.5 && mom>=0.35 && TMath::Abs(nsigmaTPCPi)<3.0) return true;
//   // else if(mom>=0.5 && TMath::Abs(nsigmaTPCPi)<2.0) return true;
//   // else return false;
//   //     } 
 
//   //   else if(TMath::Abs(nsigmaTOFPi)<3.0 && TMath::Abs(nsigmaTPCPi)<3.0) return true; //TPC+TOF
//   // }
//   // //else if(nsigmaTOFPi<-10.) //p > 0.65 + no tof == kfalse
//   // else if(mom>0.65 && nsigmaTOFPi>3) //p > 0.65 + no tof == kfalse
//   //   {
//   //     return false;
//   //   }
//   // else if(mom<1.5 && TMath::Abs(nsigmaTOFPi)<3.0 && TMath::Abs(nsigmaTPCPi)<5.0) return true;
//   // else if(mom>=1.5 && TMath::Abs(nsigmaTOFPi)<2.0 && TMath::Abs(nsigmaTPCPi)<5.0) return true;
  
//   // else 
//   //   return false;
// }
// /*
// //----------------------------------------------------------
// bool AliAnalysisTaskK0SPFemto::IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
// {
//   if (fNsigmaTPCTOF) {
//     if (mom > 0.5) {
//       //        if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP )/TMath::Sqrt(2) < 3.0)
//       if (TMath::Hypot( nsigmaTOFK, nsigmaTPCK ) < fNsigma)
//       return true;
//     }
//     else {
//       if (TMath::Abs(nsigmaTPCK) < fNsigma)
//       return true;
//     }
//   }
//   else {
//     if(mom<0.4)
//       {
//       if(nsigmaTOFK<-999.)
//         {
// 	    if(TMath::Abs(nsigmaTPCK)<2.0) return true;
// 	      }
// 	      else if(TMath::Abs(nsigmaTOFK)<3.0 && TMath::Abs(nsigmaTPCK)<3.0) return true;
//       }
//     else if(mom>=0.4 && mom<=0.6)
//       {
//       if(nsigmaTOFK<-999.)
//         {
// 	    if(TMath::Abs(nsigmaTPCK)<2.0) return true;
// 	      }
// 	      else if(TMath::Abs(nsigmaTOFK)<3.0 && TMath::Abs(nsigmaTPCK)<3.0) return true;
//       }
//     else if(nsigmaTOFK<-999.)
//       {
//       return false;
//       }
//     else if(TMath::Abs(nsigmaTOFK)<3.0 && TMath::Abs(nsigmaTPCK)<3.0) return true;
//   }
//   return false;
// }
// //----------------------------------------------------------
// bool AliAnalysisTaskK0SPFemto::IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP)
// {
//   if (fNsigmaTPCTOF) {
//     if (mom > 0.5) {
// //        if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP )/TMath::Sqrt(2) < 3.0)
//         if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP ) < fNsigma)
//             return true;
//     } else if (TMath::Abs(nsigmaTPCP) < fNsigma) {
//       return true;
//     }
//   }
//   else if (fNsigmaTPConly) {
//     if (TMath::Abs(nsigmaTPCP) < fNsigma)
//       return true;
//   }
//   else {
//     if (mom > 0.8 && mom < 2.5) {
//       if ( TMath::Abs(nsigmaTPCP) < 3.0 && TMath::Abs(nsigmaTOFP) < 3.0)
//         return true;
//     }
//     else if (mom > 2.5) {
//       if ( TMath::Abs(nsigmaTPCP) < 3.0 && TMath::Abs(nsigmaTOFP) < 2.0)
//         return true;
//     }
//     else {
//       if (TMath::Abs(nsigmaTPCP) < 3.0)
//         return true;
//     }
//   }
//   return false;
// }
// */

//-----------------------------------------------------------------------------------------------
void AliAnalysisTaskK0SPFemto::Terminate(const Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query

  if (!GetOutputData(0)) return;
}
