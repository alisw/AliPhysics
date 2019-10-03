
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
 **************************************************************************/
/* $Id$ */

//
//-----------------------------------------------------------------------------
/// \class AliAnalysisTaskNorm
/// Analysis task for muon trigger normalization
/// The output is a list of histograms and collections of counters.
/// The macro class can run on AODs or ESDs.
///
/// \author Cynthia Hadjidakis
//-----------------------------------------------------------------------------

#include "Riostream.h"

//ROOT includes
#include "TChain.h"
#include "TGraph.h"

//STEERBase includes
#include "AliVEvent.h"

//STEER includes
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliInputEventHandler.h"
#include "AliCentrality.h"
#include "AliAODTracklets.h"
#include "AliCounterCollection.h"
#include "AliTriggerConfiguration.h"

//ANALYSIS includes
#include "AliAnalysisManager.h"
#include "AliAnalysisUtils.h"

//PWGmuon include
#include "AliMuonEventCuts.h"

//PWGmuondep include
//#include "AliAnalysisTriggerScalers.h"

#include "AliAnalysisTaskNorm.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskNorm)

//________________________________________________________________________
AliAnalysisTaskNorm::AliAnalysisTaskNorm(const char *name) // All data members should be initialised here
  : AliAnalysisTaskSE(name),
    fRunCounters(0x0),
    fEventCounters(0x0),
    fListVertex(0x0),
    fListV0A(0x0),
    fListZN(0x0),
    fIsESD(kFALSE),
    fIsMC(kFALSE),
    fBeamConf(0),
    fMuonEventCuts(new AliMuonEventCuts("stdEventCuts","stdEventCuts")),
    fSCentEst(0x0), 
    fCentBin(0),
    fSCentBin(0x0),
    fTrackletsBin(0),
    fSTrackletsBin(0x0),
    fV0AMultBin(0),
    fSV0AMultBin(0x0)
{
   // Constructor
   // Define input and output slots here (never in the dummy constructor)
   // Input slot #0 works with a TChain - it is connected to the default input container
  DefineInput(0, TChain::Class());
  // Output slot #1 writes into a AliCounterCollection
  DefineOutput(1, AliCounterCollection::Class()); 
  // Output slot #2 writes into a AliCounterCollection
  DefineOutput(2, AliCounterCollection::Class()); 
  // Output slot #3 writes into a TObjArray
  DefineOutput(3, TObjArray::Class()); 
  // Output slot #4 writes into a TObjArray
  DefineOutput(4, TObjArray::Class()); 
  // Output slot #5 writes into a TObjArray
  DefineOutput(5, TObjArray::Class()); 

}

//_____________________________________________________________________________
AliAnalysisTaskNorm::~AliAnalysisTaskNorm()
{
  // Destructor. Clean-up the output list, but not the histograms that are put inside
  // (the list is owner and will clean-up these histograms). Protect in PROOF case.
  if ( !AliAnalysisManager::GetAnalysisManager()->IsProofMode() ) {
    if (fEventCounters) delete fEventCounters;   
    if (fRunCounters) delete fRunCounters;   
    if (fListVertex) delete fListVertex;
    if (fListV0A) delete fListV0A;
    if (fListZN) delete fListZN;
  }

  if (fSCentEst) delete fSCentEst;
  if (fSCentBin) delete fSCentBin;
  if (fSTrackletsBin) delete fSTrackletsBin;
  if (fSV0AMultBin) delete fSV0AMultBin;

  delete fMuonEventCuts;    
}

//_____________________________________________________________________________
void AliAnalysisTaskNorm::UserCreateOutputObjects()
{

  //initialize centrality bins and estimators
  TString centEst[] = {"CL1","V0A","ZN","any"};//,"V0C","ZNA","ZNC"};//,"ZNA","ZNC","ZPA","ZPC"};
  Int_t nCentEst = sizeof(centEst)/sizeof(centEst[0]);
  TString listOfCentEst = "";

  fSCentEst = new TObjArray(nCentEst);
  fSCentEst->SetOwner(kTRUE);
  for (Int_t iCentEst = 0; iCentEst < nCentEst; iCentEst++) {
    TObjString *obj = new TObjString(centEst[iCentEst].Data());
    fSCentEst->AddAt(obj,iCentEst);
    listOfCentEst += Form("%s",centEst[iCentEst].Data());
    if ( iCentEst < nCentEst - 1 ) listOfCentEst += Form("/");
  }

  Int_t centTabBin[] = {0,2,5,10,20,40,60,80,100};
  Int_t nCentBin = sizeof(centTabBin)/sizeof(centTabBin[0]);
  nCentBin -= 1;
  fCentBin.Set(nCentBin+1, centTabBin);

  fSCentBin = new TObjArray(nCentBin);
  fSCentBin->SetOwner(kTRUE);

  Int_t min = 0 , max = 0;
  for ( Int_t iBin = 0; iBin < nCentBin+1; iBin++ ) {
    max = (Int_t)centTabBin[iBin];
    if ( iBin >= 1 ) {
      TObjString *sobj = new TObjString(Form("%d-%d",min,max));
      fSCentBin->AddAt(sobj,iBin-1);
    }
    min = (Int_t) centTabBin[iBin];
  }

  TString listOfCentBin = "";
  for (Int_t iBin = 0; iBin < nCentBin; iBin++) {
    TObjString *obj = (TObjString*) fSCentBin->At(iBin);
    if(obj) {
      listOfCentBin += Form("%s/",(obj->GetString()).Data());
    }
  }
  listOfCentBin += "other/";
  listOfCentBin += "any";

 
  //  Int_t trTabBin[] = {0,1,30,42,56,72,300};
  Int_t trTabBin[] = {0,1,15,30,35,40,45,50,60,70,100,300};
  Int_t nTrBin = sizeof(trTabBin)/sizeof(trTabBin[0]);
  nTrBin -= 1;
  fTrackletsBin.Set(nTrBin+1,trTabBin);

  fSTrackletsBin = new TObjArray(nTrBin);
  fSTrackletsBin->SetOwner(kTRUE);

  min = max = 0;
  for ( Int_t iBin = 0; iBin < nTrBin+1; iBin++ ) {
    max = trTabBin[iBin];
    if ( iBin >= 1 ) {
      TObjString *sobj = new TObjString(Form("%d-%d",min,max));
      fSTrackletsBin->AddAt(sobj,iBin-1);
    }
    min = trTabBin[iBin];
  }
  TString listOfTrackletsBin = "";
  for (Int_t iBin = 0; iBin < nTrBin; iBin++) {
    TObjString *obj = (TObjString*) fSTrackletsBin->At(iBin);
    if(obj) {
      listOfTrackletsBin += Form("%s/",(obj->GetString()).Data());
    }
  }
  listOfTrackletsBin += "other";

  // printf("%s\n",listOfTrackletsBin.Data());
  Int_t v0AMultTabBin[] = {0,15,40,80,140,200,280,400,500,700};
  Int_t nV0AMultBin = sizeof(v0AMultTabBin)/sizeof(v0AMultTabBin[0]);
  nV0AMultBin -= 1;
  fV0AMultBin.Set(nV0AMultBin+1,v0AMultTabBin);

  fSV0AMultBin = new TObjArray(nV0AMultBin);
  fSV0AMultBin->SetOwner(kTRUE);

  min = max = 0;
  for ( Int_t iBin = 0; iBin < nV0AMultBin+1; iBin++ ) {
    max = v0AMultTabBin[iBin];
    if ( iBin >= 1 ) {
      TObjString *sobj = new TObjString(Form("%d-%d",min,max));
      fSV0AMultBin->AddAt(sobj,iBin-1);
    }
    min = v0AMultTabBin[iBin];
  }
  TString listOfV0AMultBin = "";
  for (Int_t iBin = 0; iBin < nV0AMultBin; iBin++) {
    TObjString *obj = (TObjString*) fSV0AMultBin->At(iBin);
    if(obj) {
      listOfV0AMultBin += Form("%s/",(obj->GetString()).Data());
    }
  }
  listOfV0AMultBin += "other";

  //  printf("%s\n",listOfV0AMultBin.Data());

  //create event counters
  fEventCounters = new AliCounterCollection("eventCounters");
  fEventCounters->AddRubric("run",1000000);
  fEventCounters->AddRubric("trigger",1000000);
  fEventCounters->AddRubric("norm", "MBinMB/MSLinMSL/MULinMUL/MSLinMB/MULinMSL/MULinMB/other/any");
  fEventCounters->AddRubric("physsel", "yes/no");
  //  fEventCounters->AddRubric("isspdpu", "n3d6/n4d6/n5d6/mv/any");
  fEventCounters->AddRubric("isspdpu", "n3d6/n4d6/n5d6/mv/notn4d6/notmv/any");
  //  fEventCounters->AddRubric("ntracklets",listOfTrackletsBin.Data());
  //fEventCounters->AddRubric("v0amult",listOfV0AMultBin.Data());
  //  fEventCounters->AddRubric("outofbunchpileup", "yes/no");
  //  fEventCounters->AddRubric("spdpu", "yes/no");
  fEventCounters->AddRubric("centest",listOfCentEst.Data());
  fEventCounters->AddRubric("centbin", listOfCentBin.Data());
  fEventCounters->AddRubric("vertex", "yes/no");
  fEventCounters->AddRubric("vertexcut", "abszb10/zbm10/zap10");

  fEventCounters->Init();

  //create run counters
  fRunCounters = new AliCounterCollection("runCounters");
  fRunCounters->AddRubric("run",1000000);
  //  fRunCounters->AddRubric("L0b",1000000);
  fRunCounters->AddRubric("l0brate",1000000);
  fRunCounters->AddRubric("l2a",1000000);
  fRunCounters->AddRubric("ncoll",1000000);
  fRunCounters->AddRubric("trigger",1000000);

  fRunCounters->Init();
  
  //histograms for vertex
  fListVertex = new TObjArray(2000);
  fListVertex->SetOwner();

  Float_t zrange = 50;
  Int_t nzbin = 1000, ncbin=150;
  Int_t npileup=20;

  TH1F* hVZMB = new TH1F("hVZMB"," primary vertex z (cm) for MB trigger",nzbin,-zrange,zrange);
  fListVertex->AddAtAndExpand(hVZMB, kVZMB);
  TH1F* hPileupVZMB = new TH1F("hPileupVZMB"," pile-up vertices z (cm) for MB trigger",nzbin,-zrange,zrange);
  fListVertex->AddAtAndExpand(hPileupVZMB, kPileupVZMB);
  TH1F* hVnMB = new TH1F("hVnMB"," number of vertex contributors for MB trigger",ncbin,0,ncbin);
  fListVertex->AddAtAndExpand(hVnMB, kVnCMB);
  TH1F* hNPileupMB = new TH1F("hNPileupMB"," number of pileup vertices for MB trigger",npileup,0,npileup);
  fListVertex->AddAtAndExpand(hNPileupMB, kNPileupMB);
  TH1F* hPileupVnMB = new TH1F("hPileupVnMB"," number of pileup vertices contributors for MB trigger",ncbin,0,ncbin);
  fListVertex->AddAtAndExpand(hPileupVnMB, kPileupnCMB);
  TH1F* hDeltaZMB = new TH1F("hDeltaZMB"," distance between primary vertex z and secondary vertices (cm) for MB trigger",nzbin,-zrange,zrange);
  fListVertex->AddAtAndExpand(hDeltaZMB, kDeltaZMB);

  TH1F* hVZMUL = new TH1F("hVZMUL"," primary vertex z (cm) for MUL trigger",nzbin,-zrange,zrange);
  fListVertex->AddAtAndExpand(hVZMUL, kVZMUL);
  TH1F* hPileupVZMUL = new TH1F("hPileupVZMUL"," pile-up vertices z (cm) for MUL trigger",nzbin,-zrange,zrange);
  fListVertex->AddAtAndExpand(hPileupVZMUL, kPileupVZMUL);
  TH1F* hVnMUL = new TH1F("hVnMUL"," number of contributors for MUL trigger",ncbin,0,ncbin);
  fListVertex->AddAtAndExpand(hVnMUL, kVnCMUL);
  TH1F* hNPileupMUL = new TH1F("hNPileupMUL"," number of pileup vertices for MUL trigger",npileup,0,npileup);
  fListVertex->AddAtAndExpand(hNPileupMUL, kNPileupMUL);
  TH1F* hPileupVnMUL = new TH1F("hPileupVnMUL"," number of pileup vertices contributors for MUL trigger",ncbin,0,ncbin);
  fListVertex->AddAtAndExpand(hPileupVnMUL, kPileupnCMUL);
  TH1F* hDeltaZMUL = new TH1F("hDeltaZMUL"," distance between primary vertex z and secondary vertices (cm) for MUL trigger",nzbin,-zrange,zrange);
  fListVertex->AddAtAndExpand(hDeltaZMUL, kDeltaZMUL);

  //  histograms for V0A multiplicty
  fListV0A = new TObjArray(2000);
  fListV0A->SetOwner();

  Int_t nV0Abin = 350;
  Float_t v0AMin = 0 , v0AMax = 700;
  Float_t centTabBin2[] = {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105};
  Int_t nCentBin2 = sizeof(centTabBin2)/sizeof(centTabBin2[0]);
  nCentBin2--;
  Int_t nV0ACentBin = 105;
  Float_t v0ACentMin = 0, v0ACentMax = 105;
  Int_t nTrackletsBin = 300;
  Int_t nTrackletsMin = 0, nTrackletsMax = 300;


  TH1F* hV0AMultMB = new TH1F("hV0AMultMB"," V0A mutliplicity for CINT7 trigger",nV0Abin,v0AMin,v0AMax);
  fListV0A->AddAtAndExpand(hV0AMultMB, kV0AMB);
  TH1F* hV0AMultMUL = new TH1F("hV0AMultMUL"," V0A mutliplicity for CMUL7 trigger",nV0Abin,v0AMin,v0AMax);
  fListV0A->AddAtAndExpand(hV0AMultMUL, kV0AMUL);
  TH1F* hV0ACentMB = new TH1F("hV0ACentMB"," V0A centrality for CINT7 trigger",nV0ACentBin,v0ACentMin,v0ACentMax);
  fListV0A->AddAtAndExpand(hV0ACentMB, kV0ACentMB);
  TH1F* hV0ACentMBPuCut = new TH1F("hV0ACentMBPuCut"," V0A centrality for CINT7 trigger with SPD pu cut",nV0ACentBin,v0ACentMin,v0ACentMax);
  fListV0A->AddAtAndExpand(hV0ACentMBPuCut, kV0ACentMBPuCut);
  TH1F* hV0ACentMBPuCut2 = new TH1F("hV0ACentMBPuCut2"," V0A centrality for CINT7 trigger with MV pu cut",nV0ACentBin,v0ACentMin,v0ACentMax);
  fListV0A->AddAtAndExpand(hV0ACentMBPuCut2, kV0ACentMBPuCut2);
  TH1F* hV0ACentMUL = new TH1F("hV0ACentMUL"," V0A centrality for CMUL7 trigger",nV0ACentBin,v0ACentMin,v0ACentMax);
  fListV0A->AddAtAndExpand(hV0ACentMUL, kV0ACentMUL);  
  TH1F* hV0ACentMULPuCut = new TH1F("hV0ACentMULPuCut"," V0A centrality for CMUL7 trigger with SPD pu cut",nV0ACentBin,v0ACentMin,v0ACentMax);
  fListV0A->AddAtAndExpand(hV0ACentMULPuCut, kV0ACentMULPuCut);
  TH1F* hV0ACentMULPuCut2 = new TH1F("hV0ACentMULPuCut2"," V0A centrality for CMUL7 trigger with MV pu cut",nV0ACentBin,v0ACentMin,v0ACentMax);
  fListV0A->AddAtAndExpand(hV0ACentMULPuCut2, kV0ACentMULPuCut2);
  TH2F* hV0AMultvsCentMB = new TH2F("hV0AMultvsCentMB"," V0A centrality vs multiplicity for CINT7 trigger",nV0Abin,v0AMin,v0AMax,200,0,200);
  fListV0A->AddAtAndExpand(hV0AMultvsCentMB, kV0AMultvsCentMB);
  TH2F* hV0ACentvsV0CCentMUL = new TH2F("hV0ACentvsV0CCentMUL"," V0A vs V0C centrality for CMUL7 trigger",nCentBin2,centTabBin2,nCentBin2,centTabBin2);//,centMin,centMax,nCent,centMin,centMax);
  hV0ACentvsV0CCentMUL->SetXTitle("V0A Event Multiplicity");
  hV0ACentvsV0CCentMUL->SetYTitle("V0C Event Multiplicity");
  fListV0A->AddAtAndExpand(hV0ACentvsV0CCentMUL, kV0ACentvsV0CCentMUL);
 TH2F* hV0ACentvsV0CCentMB = new TH2F("hV0ACentvsV0CCentMB"," V0A vs V0C centrality for CINT7 trigger",nCentBin2,centTabBin2,nCentBin2,centTabBin2);//,centMin,centMax,nCent,centMin,centMax);
  hV0ACentvsV0CCentMB->SetXTitle("V0A Event Multiplicity");
  hV0ACentvsV0CCentMB->SetYTitle("V0C Event Multiplicity");
  fListV0A->AddAtAndExpand(hV0ACentvsV0CCentMB, kV0ACentvsV0CCentMB);
  TH2F* hV0ACentvsCL1CentMUL = new TH2F("hV0ACentvsCL1CentMUL"," V0A vs CL1 centrality for CMUL7 trigger",nCentBin2,centTabBin2,nCentBin2,centTabBin2);//nCent,centMin,centMax,nCent,centMin,centMax);
  hV0ACentvsCL1CentMUL->SetXTitle("V0A Event Multiplicity");
  hV0ACentvsCL1CentMUL->SetYTitle("CL1 Event Multiplicity");
  fListV0A->AddAtAndExpand(hV0ACentvsCL1CentMUL, kV0ACentvsCL1CentMUL);
  TH2F* hV0ACentvsCL1CentMB = new TH2F("hV0ACentvsCL1CentMB"," V0A vs CL1 centrality for CINT7 trigger",nCentBin2,centTabBin2,nCentBin2,centTabBin2);//nCent,centMin,centMax,nCent,centMin,centMax);
  hV0ACentvsCL1CentMB->SetXTitle("V0A Event Multiplicity");
  hV0ACentvsCL1CentMB->SetYTitle("CL1 Event Multiplicity");
  fListV0A->AddAtAndExpand(hV0ACentvsCL1CentMB, kV0ACentvsCL1CentMB);
 TH2F* hV0CCentvsCL1CentMB = new TH2F("hV0CCentvsCL1CentMB"," V0C vs CL1 centrality for CINT7 trigger",nCentBin2,centTabBin2,nCentBin2,centTabBin2);//nCent,centMin,centMax,nCent,centMin,centMax);
  hV0CCentvsCL1CentMB->SetXTitle("V0C Event Multiplicity");
  hV0CCentvsCL1CentMB->SetYTitle("CL1 Event Multiplicity");
  fListV0A->AddAtAndExpand(hV0CCentvsCL1CentMB, kV0CCentvsCL1CentMB);
  TH2F* hV0AvsSPDTrackletsMB = new TH2F("hV0AvsSPDTrackletsMB"," V0A amplitude vs SPD tracklets for CINT7 trigger",nV0Abin,v0AMin,v0AMax,nTrackletsBin,nTrackletsMin,nTrackletsMax);
  hV0AvsSPDTrackletsMB->SetXTitle("V0A amplitude");
  hV0AvsSPDTrackletsMB->SetYTitle("SPD tracklets nr");
  fListV0A->AddAtAndExpand(hV0AvsSPDTrackletsMB, kV0AvsSPDTrackletsMB);
  TH2F* hV0AvsSPDTrackletsMUL = new TH2F("hV0AvsSPDTrackletsMUL"," V0A amplitude vs SPD tracklets for CMUL7 trigger",nV0Abin,v0AMin,v0AMax,nTrackletsBin,nTrackletsMin,nTrackletsMax);
  hV0AvsSPDTrackletsMUL->SetXTitle("V0A amplitude");
  hV0AvsSPDTrackletsMUL->SetYTitle("SPD tracklets nr");
  fListV0A->AddAtAndExpand(hV0AvsSPDTrackletsMUL, kV0AvsSPDTrackletsMUL);

TH2F* hV0AvsSPDTrackletsMBPuCut = new TH2F("hV0AvsSPDTrackletsMBPuCut"," V0A amplitude vs SPD tracklets for CINT7 trigger with SPD pu cut",nV0Abin,v0AMin,v0AMax,nTrackletsBin,nTrackletsMin,nTrackletsMax);
  hV0AvsSPDTrackletsMBPuCut->SetXTitle("V0A amplitude");
  hV0AvsSPDTrackletsMBPuCut->SetYTitle("SPD tracklets nr");
  fListV0A->AddAtAndExpand(hV0AvsSPDTrackletsMBPuCut, kV0AvsSPDTrackletsMBPuCut);
  TH2F* hV0AvsSPDTrackletsMULPuCut = new TH2F("hV0AvsSPDTrackletsMULPuCut"," V0A amplitude vs SPD tracklets for CMUL7 trigger with SPD pu cut",nV0Abin,v0AMin,v0AMax,nTrackletsBin,nTrackletsMin,nTrackletsMax);
  hV0AvsSPDTrackletsMULPuCut->SetXTitle("V0A amplitude");
  hV0AvsSPDTrackletsMULPuCut->SetYTitle("SPD tracklets nr");
  fListV0A->AddAtAndExpand(hV0AvsSPDTrackletsMULPuCut, kV0AvsSPDTrackletsMULPuCut);

TH2F* hV0AvsSPDTrackletsMBPuCut2 = new TH2F("hV0AvsSPDTrackletsMBPuCut2"," V0A amplitude vs SPD tracklets for CINT7 trigger with MV pu cut",nV0Abin,v0AMin,v0AMax,nTrackletsBin,nTrackletsMin,nTrackletsMax);
  hV0AvsSPDTrackletsMBPuCut2->SetXTitle("V0A amplitude");
  hV0AvsSPDTrackletsMBPuCut2->SetYTitle("SPD tracklets nr");
  fListV0A->AddAtAndExpand(hV0AvsSPDTrackletsMBPuCut2, kV0AvsSPDTrackletsMBPuCut2);
  TH2F* hV0AvsSPDTrackletsMULPuCut2 = new TH2F("hV0AvsSPDTrackletsMULPuCut2"," V0A amplitude vs SPD tracklets for CMUL7 trigger with MV pu cut",nV0Abin,v0AMin,v0AMax,nTrackletsBin,nTrackletsMin,nTrackletsMax);
  hV0AvsSPDTrackletsMULPuCut2->SetXTitle("V0A amplitude");
  hV0AvsSPDTrackletsMUL->SetYTitle("SPD tracklets nr");
  fListV0A->AddAtAndExpand(hV0AvsSPDTrackletsMULPuCut2, kV0AvsSPDTrackletsMULPuCut2);

  //  histograms for ZN multiplicity
  fListZN = new TObjArray(2000);
  fListZN->SetOwner();

  Int_t nZNbin = 2000;
  Float_t zNMin = 0 , zNMax = 2000;
  Int_t nZNCentBin = 110;
  Float_t zNCentMin = 0, zNCentMax = 110;

  TH1F* hZNMultMB = new TH1F("hZNMultMB"," ZN mutliplicity for CINT7 trigger",nZNbin,zNMin,zNMax);
  fListZN->AddAtAndExpand(hZNMultMB, kZNMB);
  TH1F* hZNMultMUL = new TH1F("hZNMultMUL"," ZN mutliplicity for CMUL7 trigger",nZNbin,zNMin,zNMax);
  fListZN->AddAtAndExpand(hZNMultMUL, kZNMUL);
  TH1F* hZNCentMB = new TH1F("hZNCentMB"," ZN centrality distribution for CINT7 trigger",nZNCentBin,zNCentMin,zNCentMax);
  fListZN->AddAtAndExpand(hZNCentMB, kZNCentMB);
  TH1F* hZNCentMUL = new TH1F("hZNCentMUL"," ZN centrality distribution for CMUL7 trigger",nZNCentBin,zNCentMin,zNCentMax);
  fListZN->AddAtAndExpand(hZNCentMUL, kZNCentMUL);
  TH1F* hZNCentMBPuCut = new TH1F("hZNCentMBPuCut"," ZN centrality distribution for CINT7 trigger with SPD pile-up cut",nZNCentBin,zNCentMin,zNCentMax);
  fListZN->AddAtAndExpand(hZNCentMBPuCut, kZNCentMBPuCut);
  TH1F* hZNCentMULPuCut = new TH1F("hZNCentMULPuCut"," ZN centrality distribution for CMUL7 trigger with SPD pile-up cut",nZNCentBin,zNCentMin,zNCentMax);
  fListZN->AddAtAndExpand(hZNCentMULPuCut, kZNCentMULPuCut);
 TH1F* hZNCentMBPuCut2 = new TH1F("hZNCentMBPuCut2"," ZN centrality distribution for CINT7 trigger with MV pile-up cut",nZNCentBin,zNCentMin,zNCentMax);
  fListZN->AddAtAndExpand(hZNCentMBPuCut2, kZNCentMBPuCut2);
  TH1F* hZNCentMULPuCut2 = new TH1F("hZNCentMULPuCut2"," ZN centrality distribution for CMUL7 trigger with MV pile-up cut",nZNCentBin,zNCentMin,zNCentMax);
  fListZN->AddAtAndExpand(hZNCentMULPuCut2, kZNCentMULPuCut2);

  TH2F* hZNMultvsCentMB = new TH2F("hZNMultvsCentMB"," ZN centrality vs multiplicity for CINT7 trigger",nZNbin,zNMin,zNMax,250,-50,200);
  fListZN->AddAtAndExpand(hZNMultvsCentMB, kZNMultvsCentMB);
  TH2F* hZNvsSPDTrackletsMB = new TH2F("hZNvsSPDTrackletsMB"," ZN multiplicity vs SPD tracklets for CINT7 trigger",nZNbin,zNMin,zNMax,nTrackletsBin,nTrackletsMin,nTrackletsMax);
  hZNvsSPDTrackletsMB->SetXTitle("ZN multiplicity");
  hZNvsSPDTrackletsMB->SetYTitle("SPD tracklets nr");
  fListZN->AddAtAndExpand(hZNvsSPDTrackletsMB, kZNvsSPDTrackletsMB);
  TH2F* hZNvsSPDTrackletsMUL = new TH2F("hZNvsSPDTrackletsMUL"," ZN multiplicity vs SPD tracklets for CMUL7 trigger",nZNbin,zNMin,zNMax,nTrackletsBin,nTrackletsMin,nTrackletsMax);
  hZNvsSPDTrackletsMUL->SetXTitle("ZN multiplicity");
  hZNvsSPDTrackletsMUL->SetYTitle("SPD tracklets nr");
  fListZN->AddAtAndExpand(hZNvsSPDTrackletsMUL, kZNvsSPDTrackletsMUL);
  TH2F* hZNvsSPDTrackletsMBPuCut = new TH2F("hZNvsSPDTrackletsMBPuCut"," ZN multiplicity vs SPD tracklets for CINT7 trigger with SPD pile-up cut",nZNbin,zNMin,zNMax,nTrackletsBin,nTrackletsMin,nTrackletsMax);
  hZNvsSPDTrackletsMBPuCut->SetXTitle("ZN multiplicity");
  hZNvsSPDTrackletsMBPuCut->SetYTitle("SPD tracklets nr");
  fListZN->AddAtAndExpand(hZNvsSPDTrackletsMBPuCut, kZNvsSPDTrackletsMBPuCut);
  TH2F* hZNvsSPDTrackletsMULPuCut = new TH2F("hZNvsSPDTrackletsMULPuCut"," ZN multiplicity vs SPD tracklets for CMUL7 trigger with SPD pile-up cut",nZNbin,zNMin,zNMax,nTrackletsBin,nTrackletsMin,nTrackletsMax);
  hZNvsSPDTrackletsMULPuCut->SetXTitle("ZN multiplicity");
  hZNvsSPDTrackletsMULPuCut->SetYTitle("SPD tracklets nr");
  fListZN->AddAtAndExpand(hZNvsSPDTrackletsMULPuCut, kZNvsSPDTrackletsMULPuCut);
TH2F* hZNvsSPDTrackletsMBPuCut2 = new TH2F("hZNvsSPDTrackletsMBPuCut2"," ZN multiplicity vs SPD tracklets for CINT7 trigger with MV pile-up cut",nZNbin,zNMin,zNMax,nTrackletsBin,nTrackletsMin,nTrackletsMax);
  hZNvsSPDTrackletsMBPuCut2->SetXTitle("ZN multiplicity");
  hZNvsSPDTrackletsMBPuCut2->SetYTitle("SPD tracklets nr");
  fListZN->AddAtAndExpand(hZNvsSPDTrackletsMBPuCut2, kZNvsSPDTrackletsMBPuCut2);
  TH2F* hZNvsSPDTrackletsMULPuCut2 = new TH2F("hZNvsSPDTrackletsMULPuCut2"," ZN multiplicity vs SPD tracklets for CMUL7 trigger with MV pile-up cut",nZNbin,zNMin,zNMax,nTrackletsBin,nTrackletsMin,nTrackletsMax);
  hZNvsSPDTrackletsMULPuCut2->SetXTitle("ZN multiplicity");
  hZNvsSPDTrackletsMULPuCut2->SetYTitle("SPD tracklets nr");
  fListZN->AddAtAndExpand(hZNvsSPDTrackletsMULPuCut2, kZNvsSPDTrackletsMULPuCut2);

  PostData(1,fEventCounters);
  PostData(2,fRunCounters);
  PostData(3,fListVertex);
  PostData(4,fListV0A);
  PostData(5,fListZN);

}

//________________________________________________________________________
Bool_t AliAnalysisTaskNorm::IsAODEvent ( const AliVEvent* event )
{
  /// Check if event is from ESD or AOD
  return ( event->IsA() == AliAODEvent::Class() );
}

//________________________________________________________________________
TString AliAnalysisTaskNorm::GetFiredTriggerClasses ( const AliVEvent* event )
{
  /// Check if track is from ESD or AOD
  return ( IsAODEvent(event) ) ? static_cast<const AliAODEvent*>(event)->GetFiredTriggerClasses() : static_cast<const AliESDEvent*>(event)->GetFiredTriggerClasses();
}

//_____________________________________________________________________________
void AliAnalysisTaskNorm::UserExec(Option_t *)
{
  
  Bool_t keepEvent = fMuonEventCuts->IsSelected(fInputHandler);
  TString sListOfTrig = GetFiredTriggerClasses(InputEvent());

  if ( DebugLevel() > 0 ) {
    AliInfo(Form("keepEvent %d - isEsd=%d - isMC=%d - beam conf = %s List of trigger=%s",keepEvent,fIsESD,fIsMC,fBeamConf.Data(),sListOfTrig.Data()));
    fMuonEventCuts->Print();
 }

  if(!keepEvent) return;

  const TObjArray* selectedTrigClasses = fMuonEventCuts->GetSelectedTrigClassesInEvent(InputEvent());
 
  AliAODEvent *aod = 0;
  AliESDEvent *esd = 0;
  if ( !fIsESD ){
    aod = static_cast<AliAODEvent *>(InputEvent());
  }
  else {
    esd = dynamic_cast<AliESDEvent *>(InputEvent());
  }
  if ( !esd && !aod ){
    AliFatal(Form("ERROR: Could not retrieve ESD (fIsESD=%d) or AOD event. Return!",fIsESD));
    return;
  }

  AliCentrality *centrality =  ( fIsESD ) ? esd->GetCentrality() : aod->GetCentrality();

  //Get list of normalization factor
  TList* normFactorList = BuildListOfNormFactor(selectedTrigClasses);
  //Get list of all triggers
  TList* triggerList = BuildListOfTrigger(selectedTrigClasses);
  //Get list of all centrality estimators + bins
  TList* centralityList = BuildListOfCentrality(centrality);
  //Get list of spd pile-up cases
  TList* spdPileUpList = BuildListOfPileUp((AliVEvent*)InputEvent());
  //Get list of ntracklets 
  TList* trackletsList = BuildListOfTracklets((AliVEvent*)InputEvent());
  //Get list of V0A multiplicity bin
  TList* v0AMultList = BuildListOfV0AMult((AliVEvent*)InputEvent());

  if ( DebugLevel() > 0 ){
    Print();
    cout<<"Trigger List"<<endl;
    triggerList->Print();
    cout<<"Norm factor List"<<endl;
    normFactorList->Print();
    cout<<"Centrality List"<<endl;
    centralityList->Print();
    cout<<"Spd Pile-Up List"<<endl;
    spdPileUpList->Print();
    cout<<"Tracklet List"<<endl;
    trackletsList->Print();
    cout<<"V0AMult List"<<endl;
    v0AMultList->Print();
    cout<<"Trigger Mask "<<(InputEvent())->GetTriggerMask()<<endl;
  }

  //Get run number
  Int_t runNr = ( (AliVEvent*) InputEvent() )->GetRunNumber();

  //Get Phys Sel.
  Bool_t physSel = kFALSE;
  UInt_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected());
  if (isSelected) physSel = kTRUE;

  //Get out-of-bunch pile-up
  Bool_t outofbcpileup = kFALSE;
  AliAnalysisUtils analysisUtils;
  outofbcpileup =  analysisUtils.IsOutOfBunchPileUp((AliVEvent*)InputEvent());

  //vertex
  Bool_t isVertex = kFALSE;
  TString sVertexCut = "";
  Int_t nVertexContributors = (((AliVEvent*) InputEvent() )->GetPrimaryVertex())->GetNContributors();
  if (nVertexContributors > 0) isVertex = kTRUE;
  Double_t zVertex = (((AliVEvent*) InputEvent() )->GetPrimaryVertex())->GetZ();
  if (TMath::Abs(zVertex) <= 10) sVertexCut = "vertexcut:abszb10";
  else if (zVertex < -10) sVertexCut = "vertexcut:zbm10";
  else if (zVertex > +10) sVertexCut = "vertexcut:zap10";


  //Fill some histograms for the vertex (primary and secondary)
  if (physSel && isVertex && sVertexCut.Contains("abszb10") ) FillHistoPileUpVertices((AliVEvent*)InputEvent(),selectedTrigClasses);

  //Fill some control histograms for the multiplicity
  if (physSel) FillHistoMult((AliVEvent*)InputEvent(),selectedTrigClasses);

  //Fill event counters
  //  FillEventCounters(runNr,triggerList,normFactorList,physSel,t0pileup,bgID);
  FillEventCounters(runNr,triggerList,normFactorList,centralityList,spdPileUpList,trackletsList,v0AMultList,physSel,isVertex,sVertexCut,outofbcpileup);

  //clean memory
  //DO NOT clean selectedTrigClasses, it is a constant
  delete normFactorList;
  delete triggerList;
  delete centralityList;
  delete spdPileUpList;
  delete trackletsList;
  delete v0AMultList;

  PostData(1, fEventCounters);
  PostData(2, fRunCounters);
  PostData(3, fListVertex);
  PostData(4, fListV0A);
  PostData(5, fListZN);

  return;

}


//_____________________________________________________________________________
void AliAnalysisTaskNorm::Terminate(Option_t *)
{

  fEventCounters = dynamic_cast<AliCounterCollection *>(GetOutputData(1));
  if (!fEventCounters) { AliError("Could not retrieve AliCounterCollection* fEventCounters"); return; }

  fRunCounters = dynamic_cast<AliCounterCollection *>(GetOutputData(2));
  if (!fRunCounters) { AliError("Could not retrieve AliCounterCollection* fRunCounters"); return; }

  fListVertex = dynamic_cast<TObjArray *>(GetOutputData(3));
  if (!fListVertex) { AliError("Could not retrieve TObjArray* fListVertex"); return; }

  fListV0A = dynamic_cast<TObjArray *>(GetOutputData(4));
  if (!fListV0A) { AliError("Could not retrieve TObjArray* fListV0A"); return; }

  fListZN = dynamic_cast<TObjArray *>(GetOutputData(5));
  if (!fListZN) { AliError("Could not retrieve TObjArray* fListZN"); return; }
  
  
}


//_____________________________________________________________________________
void AliAnalysisTaskNorm::Print(Option_t *) const
{

  cout << ClassName() << " - " << GetName() << " - " << fIsESD <<" - "<< fIsMC << endl;

}


//_____________________________________________________________________________
void AliAnalysisTaskNorm::Print(TObjArray *obj) const
{
  
  if (!obj) return;

  for ( Int_t itrig=0; itrig<obj->GetEntries(); ++itrig ) {
    TString sName = ((TObjString*)obj->At(itrig))->GetString();
    cout<<sName<<" - ";
  }
  cout<<endl;
  
}

//________________________________________________________________________
void AliAnalysisTaskNorm::FillHistoPileUpVertices(const AliVEvent *event, const TObjArray *trigClass)
{

  Bool_t isMB = kFALSE;
  Bool_t isMUL = kFALSE;
  
  for ( Int_t itrig = 0; itrig < trigClass->GetEntries(); itrig++ ) {
    TString sTrigClassName = ((TObjString*)trigClass->At(itrig))->GetString();
    if ( sTrigClassName.Contains("CINT7-B-NOPF-ALLNOTRD") || sTrigClassName.Contains("CINT7-B-NOPF-MUFAST") ) 
      isMB = kTRUE;
    if ( sTrigClassName.Contains("CMUL7-B-NOPF-MUON") || sTrigClassName.Contains("CMUL7-B-NOPF-ALLNOTRD") || sTrigClassName.Contains("CMUL7-B-NOPF-MUFAST") ) 
      isMUL = kTRUE; 
  }

  // fill histo related to pileup vertices
  const AliAODEvent *aod = 0;
  const AliESDEvent *esd = 0;
  
  if ( !fIsESD ){
    aod = static_cast<const AliAODEvent *>(event);
  }
  else {
    esd = dynamic_cast<const AliESDEvent *>(event);
  }
  if ( !esd && !aod ){
    AliFatal(Form("ERROR: Could not retrieve ESD (fIsESD=%d) or AOD event. Return!",fIsESD));
    return;
  }

  Double_t zPV = 0;
  Int_t nContrPV = -1;

  //primary vertex  
  AliVVertex *primaryVertex = (fIsESD) ? (AliVVertex*) esd->GetPrimaryVertexSPD() : (AliVVertex*) aod->GetPrimaryVertexSPD();  
  
  if (primaryVertex) {
    zPV = primaryVertex->GetZ();
    nContrPV = primaryVertex->GetNContributors();
    //printf("primary vertex: z= %.2f ncontr = %d\n",zPV,nContrPV);
  }
  if(isMB) {
    ((TH1F*)fListVertex->UncheckedAt(kVZMB))->Fill(zPV);
    ((TH1F*)fListVertex->UncheckedAt(kVnCMB))->Fill(nContrPV);
  }
  if (isMUL){
    ((TH1F*)fListVertex->UncheckedAt(kVZMUL))->Fill(zPV);
    ((TH1F*)fListVertex->UncheckedAt(kVnCMUL))->Fill(nContrPV);
  }
  
  //pile-up vertices
  Int_t nPileUpVertices = (fIsESD) ? esd->GetNumberOfPileupVerticesSPD() : aod->GetNumberOfPileupVerticesSPD();
  AliVVertex *pileupVertex = 0;

  if(isMB)  ((TH1F*)fListVertex->UncheckedAt(kNPileupMB))->Fill(nPileUpVertices);
  if(isMUL)  ((TH1F*)fListVertex->UncheckedAt(kNPileupMUL))->Fill(nPileUpVertices);

  Double_t zV = 0;
  Int_t nContr = -1;
  
  for (Int_t iV = 0; iV < nPileUpVertices; iV++) {
    zV = 0;
    nContr = -1;
    pileupVertex = (fIsESD) ? (AliVVertex*) esd->GetPileupVertexSPD(iV) : (AliVVertex*) aod->GetPileupVertexSPD(iV);
    if (!pileupVertex) continue;
    zV = pileupVertex->GetZ();
    nContr = pileupVertex->GetNContributors();
    //printf("pile-up vertex nr %d/%d: z= %.2f ncontr = %d\n",iV,nPileUpVertices,zV,nContr);
    if(isMB){
      ((TH1F*)fListVertex->UncheckedAt(kPileupVZMB))->Fill(zV);
      ((TH1F*)fListVertex->UncheckedAt(kPileupnCMB))->Fill(nContr);
      ((TH1F*)fListVertex->UncheckedAt(kDeltaZMB))->Fill(zPV-zV);
    }
    if(isMUL){
      ((TH1F*)fListVertex->UncheckedAt(kPileupVZMUL))->Fill(zV);
      ((TH1F*)fListVertex->UncheckedAt(kPileupnCMUL))->Fill(nContr);
      ((TH1F*)fListVertex->UncheckedAt(kDeltaZMUL))->Fill(zPV-zV);
    }
  }

  return;
}

//________________________________________________________________________
void AliAnalysisTaskNorm::FillHistoMult(const AliVEvent *event, const TObjArray *trigClass)
{

  Bool_t isMB = kFALSE;
  Bool_t isMUL = kFALSE;
  
  for ( Int_t itrig = 0; itrig < trigClass->GetEntries(); itrig++ ) {
    TString sTrigClassName = ((TObjString*)trigClass->At(itrig))->GetString();
    if ( sTrigClassName.Contains("CINT7-B-NOPF-ALLNOTRD") ||  sTrigClassName.Contains("CINT7-B-NOPF-MUFAST") ) 
      isMB = kTRUE;
    if ( sTrigClassName.Contains("CMUL7-B-NOPF-MUON") || sTrigClassName.Contains("CMUL7-B-NOPF-ALLNOTRD") || sTrigClassName.Contains("CMUL7-B-NOPF-MUFAST") ) 
      isMUL = kTRUE;
  }

  // fill histo related to V0A multiplicity
  AliAODEvent *aod = 0;
  AliESDEvent *esd = 0;
  
  Bool_t spdPu = kFALSE;
  Bool_t mvPu = kFALSE;

  if ( !fIsESD ){
    aod = static_cast<AliAODEvent *> (const_cast<AliVEvent*> (event) );
    spdPu = aod->IsPileupFromSPD(4,0.6,3.,2.,5.);
    AliAnalysisUtils analysis;
    mvPu = analysis.IsPileUpMV( const_cast<AliVEvent*> (event) );
  }
  else {
    esd = dynamic_cast<AliESDEvent *> (const_cast<AliVEvent*> (event) );
    spdPu = esd->IsPileupFromSPD(4,0.6,3.,2.,5.);
    mvPu = kFALSE;
  }
  if ( !esd && !aod ){
    AliFatal(Form("ERROR: Could not retrieve ESD (fIsESD=%d) or AOD event. Return!",fIsESD));
    return;
  }

  AliVVZERO *v0Data = (fIsESD) ? (AliVVZERO*) esd->GetVZEROData() : (AliVVZERO*) aod->GetVZEROData();
  if ( !v0Data ) return;

  Float_t v0AMult = 0;
  for (Int_t i = 0 ; i < 32 ; i++) v0AMult += v0Data->GetMultiplicityV0A(i);

  AliCentrality *centrality =  dynamic_cast<AliCentrality*> ( ( fIsESD ) ? esd->GetCentrality() : aod->GetCentrality() );
  Float_t centVal = centrality->GetCentralityPercentileUnchecked("V0A");
  Float_t centValV0C = centrality->GetCentralityPercentileUnchecked("V0C");
  Float_t centValCL1 = centrality->GetCentralityPercentileUnchecked("CL1");

  Int_t nSPDTracklets = 0;
  AliAODTracklets *trackletsData = (fIsESD) ? 0 : (AliAODTracklets*) aod->GetTracklets();
  if (trackletsData) nSPDTracklets = trackletsData->GetNumberOfTracklets();

  if(isMB) {
    ((TH1F*)fListV0A->UncheckedAt(kV0AMB))->Fill(v0AMult);
    ((TH1F*)fListV0A->UncheckedAt(kV0ACentMB))->Fill(centVal);
    if (!spdPu)  ((TH1F*)fListV0A->UncheckedAt(kV0ACentMBPuCut))->Fill(centVal);
    if (!mvPu)  ((TH1F*)fListV0A->UncheckedAt(kV0ACentMBPuCut2))->Fill(centVal);
    ((TH2F*)fListV0A->UncheckedAt(kV0AMultvsCentMB))->Fill(v0AMult,centVal);
    ((TH2F*)fListV0A->UncheckedAt(kV0ACentvsV0CCentMB))->Fill(centVal,centValV0C);
    ((TH2F*)fListV0A->UncheckedAt(kV0ACentvsCL1CentMB))->Fill(centVal,centValCL1);
    ((TH2F*)fListV0A->UncheckedAt(kV0CCentvsCL1CentMB))->Fill(centValV0C,centValCL1);
    ((TH2F*)fListV0A->UncheckedAt(kV0AvsSPDTrackletsMB))->Fill(v0AMult,nSPDTracklets);
    if (!spdPu) ((TH2F*)fListV0A->UncheckedAt(kV0AvsSPDTrackletsMBPuCut))->Fill(v0AMult,nSPDTracklets);
    if (!mvPu) ((TH2F*)fListV0A->UncheckedAt(kV0AvsSPDTrackletsMBPuCut2))->Fill(v0AMult,nSPDTracklets);

  }
  if (isMUL){
    ((TH1F*)fListV0A->UncheckedAt(kV0AMUL))->Fill(v0AMult);
    ((TH1F*)fListV0A->UncheckedAt(kV0ACentMUL))->Fill(centVal);
    if (!spdPu)  ((TH1F*)fListV0A->UncheckedAt(kV0ACentMULPuCut))->Fill(centVal);
    if (!mvPu)  ((TH1F*)fListV0A->UncheckedAt(kV0ACentMULPuCut2))->Fill(centVal);
    ((TH2F*)fListV0A->UncheckedAt(kV0ACentvsV0CCentMUL))->Fill(centVal,centValV0C);
    ((TH2F*)fListV0A->UncheckedAt(kV0ACentvsCL1CentMUL))->Fill(centVal,centValCL1);
    ((TH2F*)fListV0A->UncheckedAt(kV0AvsSPDTrackletsMUL))->Fill(v0AMult,nSPDTracklets);
    if (!spdPu) ((TH2F*)fListV0A->UncheckedAt(kV0AvsSPDTrackletsMULPuCut))->Fill(v0AMult,nSPDTracklets);
    if(!mvPu)  ((TH2F*)fListV0A->UncheckedAt(kV0AvsSPDTrackletsMULPuCut2))->Fill(v0AMult,nSPDTracklets);

  }
  
  //Fill histo related to ZN multiplicity
  Double_t znaTower = 0.;           // common PMT of ZNA 
  Double_t zncTower = 0.;           // common PMT of ZNC 
  Bool_t   znaFired = kFALSE;
  Bool_t   zncFired = kFALSE;

  if (fIsESD) {
    AliESDZDC *esdZDC = esd->GetESDZDC();
    if (esdZDC) {

      Int_t detChZNA  = esdZDC->GetZNATDCChannel(); //should be 12 for Run1 and 18 for Run2
      Int_t detChZNC  = esdZDC->GetZNCTDCChannel(); //should be 10 for Run1 and 16 for Run2

      for (Int_t j = 0; j < 4; ++j) 
	if (esdZDC->GetZDCTDCData(detChZNA,j) != 0) 
	  znaFired = kTRUE;
      
      for (Int_t j = 0; j < 4; ++j) 
	if (esdZDC->GetZDCTDCData(detChZNC,j) != 0) 
	  zncFired = kTRUE;   

      const Double_t *ZNAtower = esdZDC->GetZN2TowerEnergy(); 
      const Double_t *ZNCtower = esdZDC->GetZN1TowerEnergy();
      if (znaFired) znaTower = ZNAtower[0];
      if (zncFired) zncTower = ZNCtower[0];
    }
  }
  else {
    AliAODZDC *aodZDC = aod->GetZDCData();
    if (aodZDC) {
      //this works only on full AOD
      const Double_t *ZNAtower = aodZDC->GetZNATowerEnergy(); 
      const Double_t *ZNCtower = aodZDC->GetZNCTowerEnergy();
      znaTower = ZNAtower[0];
      zncTower = ZNCtower[0];
    }
  }


  Double_t zNE = 0;

  if ( !fBeamConf.CompareTo("p-Pb") ) zNE = znaTower;
  else zNE = zncTower;
  
  centVal = 0;
  if ( !fBeamConf.CompareTo("p-Pb") ) centVal = centrality->GetCentralityPercentileUnchecked("ZNA");
  else centVal = centrality->GetCentralityPercentileUnchecked("ZNC");

  //  AliInfo(Form("ZN energy %.2f centrality %.2f ",zNE,centVal));

  if(isMB) {
    ((TH1F*)fListZN->UncheckedAt(kZNMB))->Fill(zNE);
    ((TH1F*)fListZN->UncheckedAt(kZNCentMB))->Fill(centVal);
    //test Cynthia
    //    printf("test2 zNBins = %d zNCentMax= %f\n",((TH1F*)fListZN->UncheckedAt(kZNCentMB))->GetNbinsX(),((TH1F*)fListZN->UncheckedAt(kZNCentMB))->GetXaxis()->GetXmax());
    if (!spdPu) ((TH1F*)fListZN->UncheckedAt(kZNCentMBPuCut))->Fill(centVal);
    if (!mvPu) ((TH1F*)fListZN->UncheckedAt(kZNCentMBPuCut2))->Fill(centVal);
    ((TH2F*)fListZN->UncheckedAt(kZNMultvsCentMB))->Fill(zNE,centVal);
    ((TH2F*)fListZN->UncheckedAt(kZNvsSPDTrackletsMB))->Fill(zNE,nSPDTracklets);
    if (!spdPu) ((TH2F*)fListZN->UncheckedAt(kZNvsSPDTrackletsMBPuCut))->Fill(zNE,nSPDTracklets);
    if (!mvPu) ((TH2F*)fListZN->UncheckedAt(kZNvsSPDTrackletsMBPuCut2))->Fill(zNE,nSPDTracklets);


  }
  if (isMUL){
    ((TH1F*)fListZN->UncheckedAt(kZNMUL))->Fill(zNE);
    ((TH1F*)fListZN->UncheckedAt(kZNCentMUL))->Fill(centVal);
    if (!spdPu) ((TH1F*)fListZN->UncheckedAt(kZNCentMULPuCut))->Fill(centVal);
    if (!mvPu) ((TH1F*)fListZN->UncheckedAt(kZNCentMULPuCut2))->Fill(centVal);
    ((TH2F*)fListZN->UncheckedAt(kZNvsSPDTrackletsMUL))->Fill(zNE,nSPDTracklets);
    if (!spdPu) ((TH2F*)fListZN->UncheckedAt(kZNvsSPDTrackletsMULPuCut))->Fill(zNE,nSPDTracklets);
    if (!mvPu) ((TH2F*)fListZN->UncheckedAt(kZNvsSPDTrackletsMULPuCut2))->Fill(zNE,nSPDTracklets);
  }

  return;
}

//________________________________________________________________________
TList* AliAnalysisTaskNorm::BuildListOfPileUp(const AliVEvent *event)
{
 /// build the list of different cases of spd pileup (different options)
  /// returned TList must be deleted by user
    
  TList* list = new TList();
  list->SetOwner();
  
  //add any
  list->AddLast(new TObjString("isspdpu:any"));

  //Get SPD pile-up (same bunch pile-up)
  const AliAODEvent *aod = 0;
  const AliESDEvent *esd = 0;
  
  if ( !fIsESD ){
    aod = static_cast<const AliAODEvent *>(event);
  }
  else {
    esd = dynamic_cast<const AliESDEvent *>(event);
  }
  if ( !esd && !aod ){
    AliFatal(Form("ERROR: Could not retrieve ESD (fIsESD=%d) or AOD event. Return!",fIsESD));
    return 0;
  }
  
  const Int_t nSpdPileUp = 6;//4 or 6
  //  Bool_t spdpileup[nSpdPileUp] = {kFALSE,kFALSE,kFALSE,kFALSE};
  Bool_t spdpileup[nSpdPileUp] = {kFALSE,kFALSE,kFALSE,kFALSE,kTRUE,kTRUE};
  //  TString sSpdPileUp[nSpdPileUp] = {"n3d6","n3d8","n4d6","n4d8","n5d6","n5d8","n6d6","n6d8","multbinstd","mv"};
  //TString sSpdPileUp[nSpdPileUp] = {"n3d6","n4d6","n5d6","mv"};
  TString sSpdPileUp[nSpdPileUp] = {"n3d6","n4d6","n5d6","mv","notn4d6","notmv"};
  spdpileup[0] =  ( fIsESD ) ? esd->IsPileupFromSPD(3,0.6,3.,2.,5.) : aod->IsPileupFromSPD(3,0.6,3.,2.,5.);
  //spdpileup[1] =  ( fIsESD ) ? esd->IsPileupFromSPD(3,0.8,3.,2.,5.) : aod->IsPileupFromSPD(3,0.8,3.,2.,5.);
  Bool_t pileUp = ( fIsESD ) ? esd->IsPileupFromSPD(4,0.6,3.,2.,5.) : aod->IsPileupFromSPD(4,0.6,3.,2.,5.);
  if (pileUp) {
    spdpileup[1] = kTRUE;
    spdpileup[4] = kFALSE;
  }
  //spdpileup[3] =  ( fIsESD ) ? esd->IsPileupFromSPD(4,0.8,3.,2.,5.) : aod->IsPileupFromSPD(4,0.8,3.,2.,5.);
  spdpileup[2] =  ( fIsESD ) ? esd->IsPileupFromSPD(5,0.6,3.,2.,5.) : aod->IsPileupFromSPD(5,0.6,3.,2.,5.);
  //spdpileup[5] =  ( fIsESD ) ? esd->IsPileupFromSPD(5,0.8,3.,2.,5.) : aod->IsPileupFromSPD(5,0.8,3.,2.,5.);
  //spdpileup[6] =  ( fIsESD ) ? esd->IsPileupFromSPD(6,0.6,3.,2.,5.) : aod->IsPileupFromSPD(6,0.6,3.,2.,5.);
  //spdpileup[7] =  ( fIsESD ) ? esd->IsPileupFromSPD(6,0.8,3.,2.,5.) : aod->IsPileupFromSPD(6,0.8,3.,2.,5.);
  //spdpileup[8] =  kFALSE;
  //AliAODTracklets *tracklets = dynamic_cast<AliAODTracklets*> ((fIsESD) ? 0 : aod->GetTracklets());
  //if (tracklets) {
  //  spdpileup[8] = ( fIsESD ) ? esd->IsPileupFromSPDInMultBins() : aod->IsPileupFromSPDInMultBins();
  //}
  
  AliAnalysisUtils analysis;
  pileUp = analysis.IsPileUpMV( const_cast<AliVEvent*> (event) );
  if (pileUp) {
    spdpileup[3] = kTRUE;
    spdpileup[5] = kFALSE;
  }
  
  //  local SPDInMultBins function
  //  spdpileup[9] =  IsPileupFromSPDInMultBins(event);
  
  for ( Int_t iSpdPileUp = 0; iSpdPileUp < nSpdPileUp; iSpdPileUp++ ) {
    if (spdpileup[iSpdPileUp]) list->AddLast(new TObjString(Form("isspdpu:%s",sSpdPileUp[iSpdPileUp].Data())));
  }
  
  return list;
  
}

//________________________________________________________________________
Bool_t AliAnalysisTaskNorm::IsPileupFromSPDInMultBins(const AliVEvent *event) const {
  
  const AliAODEvent *aod = 0;
  const AliESDEvent *esd = 0;
  
  if ( !fIsESD ){
    aod = static_cast<const AliAODEvent *>(event);
  }
  else {
    esd = dynamic_cast<const AliESDEvent *>(event);
  }
  if ( !esd && !aod ){
    AliFatal(Form("ERROR: Could not retrieve ESD (fIsESD=%d) or AOD event. Return!",fIsESD));
    return 0;
  }

  AliAODTracklets *tracklets = dynamic_cast<AliAODTracklets*> ((fIsESD) ? 0 : aod->GetTracklets());
  if ( !tracklets) return kFALSE;

  Int_t nTracklets = (fIsESD) ? 0 : tracklets->GetNumberOfTracklets();
  if(nTracklets<40) return (fIsESD) ? esd->IsPileupFromSPD(3,0.8) : aod->IsPileupFromSPD(3,0.8);
  else return (fIsESD) ? esd->IsPileupFromSPD(5,0.8) : aod->IsPileupFromSPD(5,0.8);

}

//________________________________________________________________________
TList* AliAnalysisTaskNorm::BuildListOfCentrality(AliCentrality *centrality)
{
  /// build the list of centrality bins for each estimator
  /// returned TList must be deleted by user
    
  TList* list = new TList();
  list->SetOwner();
  
  if ( centrality ){
    for ( Int_t iMethod = 0; iMethod < fSCentEst->GetEntriesFast(); iMethod++ ) {
      
      TString sMethod = ((TObjString*)fSCentEst->At(iMethod))->GetString();
      TString sMethodCentrality;

      if ( !sMethod.CompareTo("any") ) {
	list->AddLast(new TObjString(Form("centest:any/centbin:any")));
	continue;
      }

      //pA or Ap?
      if ( !sMethod.CompareTo("ZN") ) {
	if ( !fBeamConf.CompareTo("p-Pb") ) sMethodCentrality = "ZNA";
	else sMethodCentrality = "ZNC";
      }
      else sMethodCentrality = sMethod;

      Float_t centVal = centrality->GetCentralityPercentileUnchecked( sMethodCentrality.Data() );
      //      if (centVal == 0) centVal = -1; Now centval==0 goes to bin 0-5
      if ( !sMethod.CompareTo("ZN") && centVal == 101 ) centVal = 100;
      Int_t centBin = -1;
      Int_t index = 0;
      while (centBin == -1 && index < fCentBin.GetSize()-1 ) {
	index++; 
	if ( centVal >= 0 && centVal <= fCentBin.At(index) ) centBin = index - 1;
      }

      if(DebugLevel() > 0) {
	if (centBin !=-1) {
	  TString sCentBin = ((TObjString*)fSCentBin->At(centBin))->GetString();
	  printf("method %s-%s centrality %f bin found %d bin from centrality table %s\n",sMethod.Data(),sMethodCentrality.Data(),centVal,centBin, sCentBin.Data());
	}
      }
       
      if (centBin == -1) {
	printf("method %s centrality %f centbin:other \n",sMethod.Data(),centVal);
	list->AddLast(new TObjString(Form("centest:%s/centbin:other",sMethod.Data() )));
	list->AddLast(new TObjString(Form("centest:%s/centbin:any",sMethod.Data())));
      }
      else {
	//Add specific event centrality estimator and bin
	TString sCentBin = ((TObjString*)fSCentBin->At(centBin))->GetString();
	list->AddLast(new TObjString(Form("centest:%s/centbin:%s",sMethod.Data(),sCentBin.Data())));
	list->AddLast(new TObjString(Form("centest:%s/centbin:any",sMethod.Data())));
      }
    }
  }
  
  return list;
}



//________________________________________________________________________
TList* AliAnalysisTaskNorm::BuildListOfV0AMult(const AliVEvent *event) {
  /// build the list of V0A multiplicity bins 
  /// returned TList must be deleted by user
    
  TList* list = new TList();
  list->SetOwner();

  const AliAODEvent *aod = 0;
  const AliESDEvent *esd = 0;
  
  if ( !fIsESD ){
    aod = static_cast<const AliAODEvent *>(event);
  }
  else {
    esd = dynamic_cast<const AliESDEvent *>(event);
  }
  if ( !esd && !aod ){
    AliFatal(Form("ERROR: Could not retrieve ESD (fIsESD=%d) or AOD event. Return!",fIsESD));
    list->AddLast(new TObjString(Form("v0amult:other")));
    return list;
  }

  AliVVZERO *v0Data = (fIsESD) ? (AliVVZERO*) esd->GetVZEROData() : (AliVVZERO*) aod->GetVZEROData();
  if ( !v0Data ) {
    list->AddLast(new TObjString(Form("v0amult:other")));
    return list;
  }

  Float_t v0AMult = 0;
  for (Int_t i = 0 ; i < 32 ; i++) v0AMult += v0Data->GetMultiplicityV0A(i);

  Int_t bin = -1, index = -1;

  while (bin == -1 && index < fV0AMultBin.GetSize()-1 ) {
    index++; 
    if ( v0AMult >= 0 && v0AMult < fV0AMultBin.At(index) ) bin = index - 1;
  }
  
  if(DebugLevel() > 0) {
    if (bin !=-1) {
      TString sBin = ((TObjString*)fSV0AMultBin->At(bin))->GetString();
      printf("test v0Amult %.2f bin found %d bin from v0AMult table %s\n",v0AMult,bin, sBin.Data());
    }
    else 
      printf("test v0AMult %.2f no bin\n",v0AMult);
  }
  
  if (bin == -1) {
    list->AddLast(new TObjString(Form("v0amult:other")));
  }
  else {
    //Add specific event centrality estimator and bin
    TString sBin = ((TObjString*)fSV0AMultBin->At(bin))->GetString();
    list->AddLast(new TObjString(Form("v0amult:%s",sBin.Data())));
  }

  return list;
}



//________________________________________________________________________
TList* AliAnalysisTaskNorm::BuildListOfTracklets(const AliVEvent *event) {
  /// build the list of tracklets multiplicity bins 
  /// returned TList must be deleted by user
    
  TList* list = new TList();
  list->SetOwner();

  const AliAODEvent *aod = 0;
  const AliESDEvent *esd = 0;
  
  if ( !fIsESD ){
    aod = static_cast<const AliAODEvent *>(event);
  }
  else {
    esd = dynamic_cast<const AliESDEvent *>(event);
  }
  if ( !esd && !aod ){
    AliFatal(Form("ERROR: Could not retrieve ESD (fIsESD=%d) or AOD event. Return!",fIsESD));
    list->AddLast(new TObjString(Form("ntracklets:other")));
    return list;
  }

  AliAODTracklets *tracklets = dynamic_cast<AliAODTracklets*> ((fIsESD) ? 0 : aod->GetTracklets());
  if ( !tracklets) {
    list->AddLast(new TObjString(Form("ntracklets:other")));
    return list;
  }

  Int_t nTracklets = (fIsESD) ? 0 : tracklets->GetNumberOfTracklets();
  Int_t trBin = -1, index = -1;

  while (trBin == -1 && index < fTrackletsBin.GetSize()-1 ) {
    index++; 
    if ( nTracklets >= 0 && nTracklets < fTrackletsBin.At(index) ) trBin = index - 1;
  }
  
  if(DebugLevel() > 0) {
    if (trBin !=-1) {
      TString sTrBin = ((TObjString*)fSTrackletsBin->At(trBin))->GetString();
      printf("test ntracklets %d bin found %d bin from ntracklets table %s\n",nTracklets,trBin, sTrBin.Data());
    }
    else 
      printf("test ntracklets %d no bin\n",nTracklets);
  }
  
  if (trBin == -1) {
    list->AddLast(new TObjString(Form("ntracklets:other")));
  }
  else {
    //Add specific event centrality estimator and bin
    TString sTrBin = ((TObjString*)fSTrackletsBin->At(trBin))->GetString();
    list->AddLast(new TObjString(Form("ntracklets:%s",sTrBin.Data())));
  }

  return list;
}



//________________________________________________________________________
TList* AliAnalysisTaskNorm::BuildListOfTrigger(const TObjArray *obj)
{
  /// build the list of trigger for the counters from the selected trigger objarray
  /// returned TList must be deleted by user
    
  TList* list = new TList();
  list->SetOwner();
  
  // add case any
  list->AddLast(new TObjString("trigger:any"));

  if ( obj ){
    for ( Int_t itrig = 0; itrig<obj->GetEntries(); itrig++ ) {
      TString sTrigClassName = ((TObjString*)obj->At(itrig))->GetString();
      //Add specific trigger
      list->AddLast(new TObjString(Form("trigger:%s",sTrigClassName.Data())));
    }
  }
  
  // add case other if no specific trigger was found
  if (list->GetSize() == 1) list->AddLast(new TObjString("trigger:other"));

  return list;
}


//_____________________________________________________________________________
Bool_t AliAnalysisTaskNorm::CheckPattern ( TString trigName, TObjArray* keepArray, TObjArray* rejectArray )
{
  for ( Int_t ipat=0; ipat<rejectArray->GetEntries(); ++ipat ) {
    if ( trigName.Contains(rejectArray->At(ipat)->GetName() ) ) return kFALSE;
  } // loop on reject pattern
  
  for ( Int_t ipat=0; ipat<keepArray->GetEntries(); ++ipat ) {
    if ( trigName.Contains(keepArray->At(ipat)->GetName() ) ) return kTRUE;
  } // loop on keep pattern
  
  return ( keepArray->GetEntries() == 0 ) ? kTRUE : kFALSE;
}


//_____________________________________________________________________________
TObjArray* AliAnalysisTaskNorm::BuildArrayOfTrigger ( const TObjArray* triggerArray, TString keepPattern, TString rejectPattern )
{

  /// build the selected list of trigger from the trigger objarray
  /// returned TObjArray must be deleted by user

  TObjArray* selectedList = new TObjArray();
  selectedList->SetOwner();
  TObjArray* rejectArray = rejectPattern.Tokenize(",");
  TObjArray* keepArray = keepPattern.Tokenize(",");
  
  for ( Int_t iTrig = 0; iTrig < triggerArray->GetEntries(); iTrig++ ){
    TString currTrigName = ((TObjString*)triggerArray->At(iTrig))->GetName();
    if ( CheckPattern(currTrigName, keepArray, rejectArray) ) selectedList->AddLast(new TObjString(currTrigName.Data()));
  }

  delete rejectArray;
  delete keepArray;

  return selectedList;

}

//________________________________________________________________________
TList* AliAnalysisTaskNorm::BuildListOfNormFactor(const TObjArray *trig ) {
  // returned TList must be deleted by user
  TList *list = new TList();
  list->SetOwner();

  //fEventCounters->AddRubric("norm", "MBinMB/MSLinMSL/MULinMUL/MSLinMB/MULinMSL/MULinMB/any");
  
  //any event
  list->AddLast(new TObjString("norm:any"));
  
  if ( trig ) {

    Bool_t isMUL = kFALSE, isMSL = kFALSE, isMB = kFALSE;
    Bool_t isMULInMSL = kFALSE, isMSLInMB = kFALSE, isMULInMB = kFALSE;
    TString sNorm;
    
    for ( Int_t itrig=0; itrig< trig->GetEntries(); itrig++ ) {
      TString sTrigClassName = ((TObjString*) trig->At(itrig))->GetString();
      if ( sTrigClassName.Contains("CINT7-B-NOPF-ALLNOTRD") || sTrigClassName.Contains("CINT7-B-NOPF-MUFAST")) 
	isMB = kTRUE;
      if ( sTrigClassName.Contains("CMUL7-B-NOPF-MUON") || sTrigClassName.Contains("CMUL7-B-NOPF-MUFAST") ) 
	isMUL = kTRUE;
      if ( sTrigClassName.Contains("CMSL7-B-NOPF-MUON") || sTrigClassName.Contains("CMSL7-B-NOPF-MUFAST") ) 
	isMSL = kTRUE;
      if ( sTrigClassName.Contains("CINT7-B-NOPF-ALLNOTRD&0MSL") || sTrigClassName.Contains("CINT7-B-NOPF-MUFAST&0MSL") ) 
	isMSLInMB = kTRUE;
      if ( sTrigClassName.Contains("CMSL7-B-NOPF-MUON&0MUL") ||  sTrigClassName.Contains("CMSL7-B-NOPF-MUFAST&0MUL") ) 
	isMULInMSL = kTRUE;
      if ( sTrigClassName.Contains("CINT7-B-NOPF-ALLNOTRD&0MUL") || sTrigClassName.Contains("CINT7-B-NOPF-MUFAST&0MUL") ) 
	isMULInMB = kTRUE;
    }
    
    if (isMB){
      sNorm = "MBinMB"; 
      list->AddLast(new TObjString(Form("norm:%s",sNorm.Data())));
    }
    if (isMSL){
      sNorm = "MSLinMSL"; 
      list->AddLast(new TObjString(Form("norm:%s",sNorm.Data())));
    }
    if (isMUL){
      sNorm = "MULinMUL"; 
      list->AddLast(new TObjString(Form("norm:%s",sNorm.Data())));
    }
    if (isMSLInMB){
      sNorm = "MSLinMB"; 
      list->AddLast(new TObjString(Form("norm:%s",sNorm.Data())));
    }
    if (isMULInMSL){
      sNorm = "MULinMSL"; 
      list->AddLast(new TObjString(Form("norm:%s",sNorm.Data())));
    }
    if (isMULInMB){
      sNorm = "MULinMB"; 
      list->AddLast(new TObjString(Form("norm:%s",sNorm.Data())));
    }

    // add case other if no specific trigger was found
    if ( list->GetSize() == 1 ) list->AddLast(new TObjString("norm:other"));

  }  

  return list;
  
}

//_____________________________________________________________________________
void AliAnalysisTaskNorm::FillEventCounters( Int_t runNr, TList *triggerList, TList *normFactorList, TList *centralityList, TList *spdPileUpList, TList * /*trackletsList*/, TList * /*v0AMultList*/, Bool_t physsel, Bool_t isVertex, TString sVertexCut, Bool_t /*outofbunchpileup*/) {

  TIter nextTriggerKey(triggerList);
  TObjString *triggerKey = 0x0;

  TIter nextNormFactorKey(normFactorList);
  TObjString *normFactorKey = 0x0;

  TIter nextCentralityKey(centralityList);
  TObjString *centralityKey = 0x0;

  TIter nextSpdPileUpKey(spdPileUpList);
  TObjString *spdPileUpKey = 0x0;

  //  TIter nextTrackletsKey(trackletsList);
  // TObjString *trackletsKey = 0x0;

  //TIter nextV0AMultKey(v0AMultList);
  //TObjString *v0AMultKey = 0x0;

  TString selected;
  selected = (physsel) ? "physsel:yes" : "physsel:no";
  selected += (isVertex) ? "/vertex:yes" : "/vertex:no";
  selected += Form("/%s",sVertexCut.Data());

  //  selected += (outofbunchpileup) ? "/outofbunchpileup:yes" : "/outofbunchpileup:no";
  

  //Loop over triggerList
  while ( ( triggerKey = (TObjString*) nextTriggerKey() )  ) {
    
    //Loop over normFactor List
    nextNormFactorKey.Reset();
    while ( ( normFactorKey = (TObjString*) nextNormFactorKey() ) ) {

      //Loop over centrality List
      nextCentralityKey.Reset();
      while ( ( centralityKey = (TObjString*) nextCentralityKey() ) ) {

	//Loop over spd pileup list
	nextSpdPileUpKey.Reset();
	while ( ( spdPileUpKey = (TObjString*) nextSpdPileUpKey() ) ) {

	  //Loop over ntracklets list
	  //nextTrackletsKey.Reset();
	  //while ( ( trackletsKey = (TObjString*) nextTrackletsKey() ) ) {

	    //Loop over V0Amult list
	    //nextV0AMultKey.Reset();
	    //while ( ( v0AMultKey = (TObjString*) nextV0AMultKey() ) ) {
	      
	      TString sEventCounters;
	      //sEventCounters = Form("run:%d/%s/%s/%s/%s/%s/%s/%s",runNr,selected.Data(),triggerKey->GetName(),normFactorKey->GetName(),centralityKey->GetName(),spdPileUpKey->GetName(),trackletsKey->GetName(),v0AMultKey->GetName());
	      sEventCounters = Form("run:%d/%s/%s/%s/%s/%s",runNr,selected.Data(),triggerKey->GetName(),normFactorKey->GetName(),centralityKey->GetName(),spdPileUpKey->GetName());
	      fEventCounters->Count(sEventCounters);
	      //	printf("event counters = %s\n",sEventCounters.Data());

	      // }//end loop on V0A multiplcity list
	      // }//end loop on ntracklets list
	}// end loop on spd pile up list
      }//end loop on centrality list
    }//end loop on normfactor list
  }//end loop on triggerList
  
  return;
}
 
