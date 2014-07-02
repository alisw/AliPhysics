/**************************************************************************
 * Contributors are not mentioned at all.                                 *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright noticxse appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//-----------------------------------------------------------------
//                 AliAnalysisTaskNucleiv2 class
//-----------------------------------------------------------------

class TTree;
class TParticle;
class TVector3;

#include "AliAnalysisManager.h"
#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <AliStack.h>

class AliESDVertex;
class AliAODVertex;
class AliESDv0;
class AliAODv0; 
class AliCascadeVertexer;

#include <iostream>

#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TNtuple.h"
#include "TGraph.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TChain.h"
#include "Riostream.h"
#include "AliLog.h"
#include "AliCascadeVertexer.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"
#include "AliAODEvent.h"
#include "AliInputEventHandler.h"
#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliAnalysisTaskNucleiv2.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"
#include "AliFlowCandidateTrack.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowEventSimple.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowEvent.h"
#include "AliFlowTrack.h"
#include "TProfile.h"

ClassImp(AliAnalysisTaskNucleiv2)

using std::cout;
using std::endl;
    
//________________________________________________________________________
AliAnalysisTaskNucleiv2::AliAnalysisTaskNucleiv2() 
: AliAnalysisTaskSE(), 
  fAnalysisType(0), 
  fCollidingSystems(0), 
  fDataType(0),
  fFillNtuple(0),
  fCentralityMin(0),
  fCentralityMax(0), 
  fCutsRP(0),
  fNullCuts(0),
  fFlowEvent(0),
  fListHist(0), 
  fHistEventMultiplicity(0), 
  fHistTrackMultiplicity(0),
  fHistTrackMultiplicityCentral(0),    
  fHistTrackMultiplicitySemiCentral(0),
  fHistTrackMultiplicityMB(0),
  fhBB(0),
  fhBBDeu(0),
  fhPtDeu(0),
  fhTOF(0),
  fhMassTOF(0),
  EPVzAvsCentrality(0), 
  EPVzCvsCentrality(0), 
  EPTPCvsCentrality(0), 
  EPVzvsCentrality(0), 
  EPTPCpvsCentrality(0), 
  EPTPCnvsCentrality(0), 
  fSubEventDPhiv205(0), 
  fSubEventDPhiv2new05(0),
  fSubEventDPhiv22040(0), 
  fSubEventDPhiv2new2040(0),
  fSubEventDPhiv24060(0), 
  fSubEventDPhiv2new4060(0),
  hCos2DeltaPhiVzAvsCentrality(0),
  hCos2DeltaPhiVzCvsCentrality(0),
  hCos2DeltaPhiVzMvsCentrality(0),
  hCos2DeltaPhiTPCfvsCentrality(0),
  hCos2DeltaPhiTPCpvsCentrality(0),
  hCos2DeltaPhiTPCnvsCentrality(0),
  hEvPlaneTPCvsEvPVz05(0),                      
  hEvPlaneTPCvsEvPVz075(0), 
  hEvPlaneTPCvsEvPVz1530(0),
  hEvPlaneTPCvsEvPVz3050(0),                      
  hEvPlaneTPCvsEvPVz2040(0),                      
  hEvPlaneTPCvsEvPVz4060(0),       
  hCos2DeltaPhivsPt075(0),                      
  hCos2DeltaPhiVZEROvsPt075(0),                 
  hCos2DeltaPhivsPt1530(0),                      
  hCos2DeltaPhiVZEROvsPt1530(0),                 
  hCos2DeltaPhivsPt3050(0),                      
  hCos2DeltaPhiVZEROvsPt3050(0),                 
  hCos2DeltaPhivsPt05(0),                      
  hCos2DeltaPhiVZEROvsPt05(0),                 
  hCos2DeltaPhivsPt2040(0),                      
  hCos2DeltaPhiVZEROvsPt2040(0),                 
  hCos2DeltaPhivsPt4060(0),                      
  hCos2DeltaPhiVZEROvsPt4060(0),           
  fESDtrackCuts(0),
  fESDtrackCutsEP(0),
  fPIDResponse(0),
  fNtuple1(0),
  tCharge(0),
  tPhi(0),              
  trpangleTPC(0),       
  tPDGCode(0),          
  tPDGCodeMum(0),       
  tIsPrimaryTr(0),
  fNtuple2(0) ,
  tCentralityMC(0),
  tPDGCodeMC(0),
  tPDGCodeMumMC(0),
  tIsPrimary(0),
  tEtaMC(0),
  tPtMC(0),
  tYMC(0)
{
  // Dummy Constructor 
  fESDtrackCuts = new AliESDtrackCuts("AliESDtrackCuts","AliESDtrackCuts");
  fESDtrackCutsEP = new AliESDtrackCuts("AliESDtrackCutsEP","AliESDtrackCutsEP");
  //
  Initialize();
}

//________________________________________________________________________
AliAnalysisTaskNucleiv2::AliAnalysisTaskNucleiv2(const char *name) 
  :  AliAnalysisTaskSE(name), 
     fAnalysisType(0), 
     fCollidingSystems(0), 
     fDataType(0),
     fFillNtuple(0),
     fCentralityMin(0),
     fCentralityMax(0),
     fCutsRP(0),  
     fNullCuts(0),
     fFlowEvent(0),
     fListHist(0), 
     fHistEventMultiplicity(0), 
     fHistTrackMultiplicity(0),
     fHistTrackMultiplicityCentral(0),    
     fHistTrackMultiplicitySemiCentral(0),
     fHistTrackMultiplicityMB(0),
     fhBB(0),
     fhBBDeu(0),
     fhPtDeu(0),
     fhTOF(0),
     fhMassTOF(0),
     EPVzAvsCentrality(0), 
     EPVzCvsCentrality(0), 
     EPTPCvsCentrality(0), 
     EPVzvsCentrality(0), 
     EPTPCpvsCentrality(0), 
     EPTPCnvsCentrality(0), 
     fSubEventDPhiv205(0), 
     fSubEventDPhiv2new05(0),
     fSubEventDPhiv22040(0), 
     fSubEventDPhiv2new2040(0),
     fSubEventDPhiv24060(0), 
     fSubEventDPhiv2new4060(0),
     hCos2DeltaPhiVzAvsCentrality(0),
     hCos2DeltaPhiVzCvsCentrality(0),
     hCos2DeltaPhiVzMvsCentrality(0),
     hCos2DeltaPhiTPCfvsCentrality(0),
     hCos2DeltaPhiTPCpvsCentrality(0),
     hCos2DeltaPhiTPCnvsCentrality(0),
     hEvPlaneTPCvsEvPVz05(0),                      
     hEvPlaneTPCvsEvPVz075(0), 
     hEvPlaneTPCvsEvPVz1530(0),
     hEvPlaneTPCvsEvPVz3050(0),                      
     hEvPlaneTPCvsEvPVz2040(0),                      
     hEvPlaneTPCvsEvPVz4060(0),   
     hCos2DeltaPhivsPt075(0),                      
     hCos2DeltaPhiVZEROvsPt075(0),                 
     hCos2DeltaPhivsPt1530(0),                      
     hCos2DeltaPhiVZEROvsPt1530(0),                 
     hCos2DeltaPhivsPt3050(0),                      
     hCos2DeltaPhiVZEROvsPt3050(0),                 
     hCos2DeltaPhivsPt05(0),                      
     hCos2DeltaPhiVZEROvsPt05(0),                 
     hCos2DeltaPhivsPt2040(0),                      
     hCos2DeltaPhiVZEROvsPt2040(0),                 
     hCos2DeltaPhivsPt4060(0),                      
     hCos2DeltaPhiVZEROvsPt4060(0),           
     fESDtrackCuts(0),
     fESDtrackCutsEP(0),
     fPIDResponse(0),
     fNtuple1(0),
     tCharge(0),
     tPhi(0),              
     trpangleTPC(0),       
     tPDGCode(0),          
     tPDGCodeMum(0),       
     tIsPrimaryTr(0),
     fNtuple2(0) ,
     tCentralityMC(0),
     tPDGCodeMC(0),
     tPDGCodeMumMC(0),
     tIsPrimary(0),
     tEtaMC(0),
     tPtMC(0),
     tYMC(0)
     
{
  // Define input and output slots here
  // Input slot #0 works with a TChain
  //DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList container ()

  //
  // create track cuts
  //
  fESDtrackCuts = new AliESDtrackCuts("AliESDtrackCuts","AliESDtrackCuts");
  //
  Initialize();

  DefineInput(0, TChain::Class());
  
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
  DefineOutput(4, AliFlowEventSimple::Class());
}

void AliAnalysisTaskNucleiv2::Initialize()
{
  //
  // updating parameters in case of changes
  //
  fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE,kTRUE);
  fESDtrackCuts->SetMaxDCAToVertexXY(3);
  fESDtrackCuts->SetMaxDCAToVertexZ(2);
  fESDtrackCuts->SetEtaRange(-0.8,0.8);
  //
  //
 
  fESDtrackCutsEP = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();  
  //  Printf("Initizialize\n");
 
}

//________________________________________________________________________
Float_t AliAnalysisTaskNucleiv2::GetEventPlaneForCandidate(AliESDtrack* track0, const TVector2* q,AliEventplane *pl, const TVector2* qsub1, const TVector2* qsub2){
 
  // remove autocorrelations 
 
  TArrayF* qx = 0x0;
  TArrayF* qy = 0x0;
  TVector2 qcopy; 
  // if(!fEtaGap){
  qx = pl->GetQContributionXArray();
  qy = pl->GetQContributionYArray();
  qcopy = *q;
  // }else {
  // if(d->Eta()<0.){
  //   qx = pl->GetQContributionXArraysub1();
  //   qy = pl->GetQContributionYArraysub1();
  //   qcopy = *qsub1;
  // }else{
  //   qx = pl->GetQContributionXArraysub2();
  //   qy = pl->GetQContributionYArraysub2();
  //   qcopy = *qsub2;
  // }
  //}
  
  TVector2 q0;
  if((track0->GetID()) < qx->fN){
    q0.Set(qx->At(track0->GetID()),qy->At(track0->GetID()));
  }
  
  qcopy = qcopy - q0;
  
  return qcopy.Phi()/2.;
  
}
//_____________________________________________________________________________
template <typename T> void AliAnalysisTaskNucleiv2::SetNullCuts(T* event)
{
    //Set null cuts
  if (fDebug) cout << " fCutsRP " << fCutsRP << endl;
  fCutsRP->SetEvent(event, MCEvent());
  fNullCuts->SetParamType(AliFlowTrackCuts::kGlobal);
  fNullCuts->SetPtRange(+1, -1); // select nothing QUICK
  fNullCuts->SetEtaRange(+1, -1); // select nothing VZERO
  fNullCuts->SetEvent(event, MCEvent());
}
//______________________________________________________________________
void AliAnalysisTaskNucleiv2::PrepareFlowEvent(Int_t iMulti, AliFlowEvent *FlowEv) const
{
    //Prepare flow events
    FlowEv->ClearFast();
    FlowEv->Fill(fCutsRP, fNullCuts);
    FlowEv->SetReferenceMultiplicity(iMulti);
    FlowEv->DefineDeadZone(0, 0, 0, 0);
    //  FlowEv->TagSubeventsInEta(-0.7, 0, 0, 0.7);
}

//________________________________________________________________________
Float_t AliAnalysisTaskNucleiv2::GetPhi0Pi(Float_t phi){
  // Sets the phi angle in the range 0-pi
  Float_t result=phi;
  while(result<0){
    result=result+TMath::Pi();
  }
  while(result>TMath::Pi()){
    result=result-TMath::Pi();
  }
   return result;
}

//_____________________________________________________________________________
void AliAnalysisTaskNucleiv2::SetCentralityParameters(Double_t CentralityMin, Double_t CentralityMax)
{
    // Set a centrality range ]min, max] 
  fCentralityMin = CentralityMin;
  fCentralityMax = CentralityMax;
  
}
//==================DEFINITION OF OUTPUT OBJECTS==============================

void AliAnalysisTaskNucleiv2::UserCreateOutputObjects()
{
  //-------------------------------------------------------
  fNullCuts = new AliFlowTrackCuts("null_cuts");

  AliFlowCommonConstants* cc = AliFlowCommonConstants::GetMaster();
  cc->SetNbinsMult(10000);
  cc->SetMultMin(0);
  cc->SetMultMax(10000);
  
  cc->SetNbinsPt(100);
  cc->SetPtMin(0);
  cc->SetPtMax(50);
    
  cc->SetNbinsPhi(180);
  cc->SetPhiMin(0.0);
  cc->SetPhiMax(TMath::TwoPi());
  
  cc->SetNbinsEta(30);
  cc->SetEtaMin(-7.0);
  cc->SetEtaMax(+7.0);
  
  cc->SetNbinsQ(500);
  cc->SetQMin(0.0);
  cc->SetQMax(3.0);
  
  //------------------------------------------------------

  fListHist = new TList();
  fListHist->SetOwner();  // IMPORTANT!

  if(! fHistEventMultiplicity ){

    fHistEventMultiplicity   = new TH1F( "fHistEventMultiplicity" , "Nb of Events" , 12 , -0.5,11.5);

    fHistEventMultiplicity->GetXaxis()->SetBinLabel(1,"All Events");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(2,"Events w/PV");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(3,"Events w/|Vz|<10cm");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(4,"Central Events");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(5,"Semi-Central Events");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(6,"MB Events");
    //from HF
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(7,"nEventsAnal");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(8,"nEvSelected");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(9,"nCandidatesSelected");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(10,"out of pt bounds");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(11,"mismatch lab");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(12,"non valid TPC EP");
    fListHist->Add(fHistEventMultiplicity);
  }

  if(! fHistTrackMultiplicity ){
    fHistTrackMultiplicity  = new TH2F( "fHistTrackMultiplicity", "Nb of Tracks MB Events |Vz| < 10", 25000,0, 25000,105,-0.5,104.5);
    fHistTrackMultiplicity->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicity->GetYaxis()->SetTitle("Percentile");
    fListHist->Add(fHistTrackMultiplicity);
  } 

  if(! fHistTrackMultiplicityCentral ){
    fHistTrackMultiplicityCentral  = new TH2F( "fHistTrackMultiplicityCentral", "Nb of Tracks MB Events |Vz| < 10", 25000,0, 25000,105,-0.5,104.5);
    fHistTrackMultiplicityCentral->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicityCentral->GetYaxis()->SetTitle("Percentile");
    fListHist->Add(fHistTrackMultiplicityCentral);
  } 
  if(! fHistTrackMultiplicitySemiCentral ){
    fHistTrackMultiplicitySemiCentral  = new TH2F( "fHistTrackMultiplicitySemiCentral", "Nb of Tracks MB Events |Vz| < 10", 25000,0, 25000,105,-0.5,104.5);
    fHistTrackMultiplicitySemiCentral->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicitySemiCentral->GetYaxis()->SetTitle("Percentile");
    fListHist->Add(fHistTrackMultiplicitySemiCentral);
  } 
  if(! fHistTrackMultiplicityMB ){
    fHistTrackMultiplicityMB  = new TH2F( "fHistTrackMultiplicityMB", "Nb of Tracks MB Events |Vz| < 10", 25000,0, 25000,105,-0.5,104.5);
    fHistTrackMultiplicityMB->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicityMB->GetYaxis()->SetTitle("Percentile");
    fListHist->Add(fHistTrackMultiplicityMB);
  } 
 
  if(! fhBB ){
    fhBB = new TH2F( "fhBB" , "BetheBlochTPC" , 240,-6,6,1000,0,1000);
    fListHist->Add(fhBB);
  }
  
  if(! fhBBDeu ){
    fhBBDeu = new TH2F( "fhBBDeu" , "BetheBlochTPC - Deuteron" , 240,-6,6,1000,0,1000);
    fListHist->Add(fhBBDeu);
  }
 
  if(!fhPtDeu  ){
    fhPtDeu = new TH2F( "fhPtDeu" , "pt corretto vs pt track - Deuteron" , 120,0,6,120,0,6);
    fListHist->Add(fhPtDeu);
  }

  if(! fhTOF ){
    fhTOF = new TH2F( "fhTOF" , "Scatter Plot TOF" , 240,-6,6,1000,0,1.2);
    fListHist->Add(fhTOF);
  }
  if(! fhMassTOF){
    fhMassTOF=new TH1F ("fhMassTOF","Particle Mass - TOF", 300,-5 ,5);
    fListHist->Add(fhMassTOF);
  }
  
  EPVzAvsCentrality = new TH2D("EPVzAvsCentrality", "EPVzAvsCentrality", 80, -2, 2,105,-0.5,105.5);
  fListHist->Add(EPVzAvsCentrality);
  EPVzCvsCentrality = new TH2D("EPVzCvsCentrality", "EPVzCvsCentrality", 80, -2, 2,105,-0.5,105.5);
  fListHist->Add(EPVzCvsCentrality);
  EPTPCvsCentrality = new TH2D("EPTPCvsCentrality", "EPTPCvsCentrality", 80, -2, 2,105,-0.5,105.5);
  fListHist->Add(EPTPCvsCentrality);
    
    
  EPVzvsCentrality = new TH2D("EPVzvsCentrality", "EPVzvsCentrality", 80, -2, 2,105,-0.5,105.5);
  fListHist->Add(EPVzvsCentrality);
  EPTPCpvsCentrality = new TH2D("EPTPCpvsCentrality", "EPTPCpvsCentrality", 80, -2, 2,105,-0.5,105.5);
  fListHist->Add(EPTPCpvsCentrality);
  EPTPCnvsCentrality = new TH2D("EPTPCnvsCentrality", "EPTPCnvsCentrality", 80, -2, 2,105,-0.5,105.5);
  fListHist->Add(EPTPCnvsCentrality);

  //1
 
  fSubEventDPhiv205 = new TProfile("fSubEventDPhiv205", "fSubEventDPhiv205", 3, 0, 3);
  fSubEventDPhiv205->GetXaxis()->SetBinLabel(1, "<cos(2(#Psi_{a} - #Psi_{b}))>");
  fSubEventDPhiv205->GetXaxis()->SetBinLabel(2, "<cos(2(#Psi_{a} - #Psi_{c}>))");
  fSubEventDPhiv205->GetXaxis()->SetBinLabel(3, "<cos(2(#Psi_{b} - #Psi_{c}>))");
  fListHist->Add(fSubEventDPhiv205);
      
  fSubEventDPhiv2new05 = new TProfile("fSubEventDPhiv2new05", "fSubEventDPhiv2new05", 3, 0, 3);
  fSubEventDPhiv2new05->GetXaxis()->SetBinLabel(1, "<cos(2(#Psi_{a} - #Psi_{b}))>");
  fSubEventDPhiv2new05->GetXaxis()->SetBinLabel(2, "<cos(2(#Psi_{a} - #Psi_{c}>))");
  fSubEventDPhiv2new05->GetXaxis()->SetBinLabel(3, "<cos(2(#Psi_{b} - #Psi_{c}>))");
  fListHist->Add(fSubEventDPhiv2new05);
  
  //2

  fSubEventDPhiv22040 = new TProfile("fSubEventDPhiv22040", "fSubEventDPhiv22040", 3, 0, 3);
  fSubEventDPhiv22040->GetXaxis()->SetBinLabel(1, "<cos(2(#Psi_{a} - #Psi_{b}))>");
  fSubEventDPhiv22040->GetXaxis()->SetBinLabel(2, "<cos(2(#Psi_{a} - #Psi_{c}>))");
  fSubEventDPhiv22040->GetXaxis()->SetBinLabel(3, "<cos(2(#Psi_{b} - #Psi_{c}>))");
  fListHist->Add(fSubEventDPhiv22040);
      
  fSubEventDPhiv2new2040 = new TProfile("fSubEventDPhiv2new2040", "fSubEventDPhiv2new2040", 3, 0, 3);
  fSubEventDPhiv2new2040->GetXaxis()->SetBinLabel(1, "<cos(2(#Psi_{a} - #Psi_{b}))>");
  fSubEventDPhiv2new2040->GetXaxis()->SetBinLabel(2, "<cos(2(#Psi_{a} - #Psi_{c}>))");
  fSubEventDPhiv2new2040->GetXaxis()->SetBinLabel(3, "<cos(2(#Psi_{b} - #Psi_{c}>))");
  fListHist->Add(fSubEventDPhiv2new2040);

  //3
 
  fSubEventDPhiv24060 = new TProfile("fSubEventDPhiv24060", "fSubEventDPhiv24060", 3, 0, 3);
  fSubEventDPhiv24060->GetXaxis()->SetBinLabel(1, "<cos(2(#Psi_{a} - #Psi_{b}))>");
  fSubEventDPhiv24060->GetXaxis()->SetBinLabel(2, "<cos(2(#Psi_{a} - #Psi_{c}>))");
  fSubEventDPhiv24060->GetXaxis()->SetBinLabel(3, "<cos(2(#Psi_{b} - #Psi_{c}>))");
  fListHist->Add(fSubEventDPhiv24060);
      
  fSubEventDPhiv2new4060 = new TProfile("fSubEventDPhiv2new4060", "fSubEventDPhiv2new4060", 3, 0, 3);
  fSubEventDPhiv2new4060->GetXaxis()->SetBinLabel(1, "<cos(2(#Psi_{a} - #Psi_{b}))>");
  fSubEventDPhiv2new4060->GetXaxis()->SetBinLabel(2, "<cos(2(#Psi_{a} - #Psi_{c}>))");
  fSubEventDPhiv2new4060->GetXaxis()->SetBinLabel(3, "<cos(2(#Psi_{b} - #Psi_{c}>))");
  fListHist->Add(fSubEventDPhiv2new4060);
  
  //------------------

  hCos2DeltaPhiVzAvsCentrality  = new TH2F("hCos2DeltaPhiVzAvsCentrality" ,"hCos2DeltaPhiVzAvsCentrality" ,100,-1.1,1.1,105,-0.5,105.5);
  hCos2DeltaPhiVzCvsCentrality  = new TH2F("hCos2DeltaPhiVzCvsCentrality" ,"hCos2DeltaPhiVzCvsCentrality" ,100,-1.1,1.1,105,-0.5,105.5);
  hCos2DeltaPhiVzMvsCentrality  = new TH2F("hCos2DeltaPhiVzMvsCentrality" ,"hCos2DeltaPhiVzMvsCentrality" ,100,-1.1,1.1,105,-0.5,105.5);
  hCos2DeltaPhiTPCfvsCentrality = new TH2F("hCos2DeltaPhiTPCfvsCentrality","hCos2DeltaPhiTPCfvsCentrality",100,-1.1,1.1,105,-0.5,105.5);
  hCos2DeltaPhiTPCpvsCentrality = new TH2F("hCos2DeltaPhiTPCpvsCentrality","hCos2DeltaPhiTPCpvsCentrality",100,-1.1,1.1,105,-0.5,105.5);
  hCos2DeltaPhiTPCnvsCentrality = new TH2F("hCos2DeltaPhiTPCnvsCentrality","hCos2DeltaPhiTPCnvsCentrality",100,-1.1,1.1,105,-0.5,105.5);

  fListHist->Add(hCos2DeltaPhiVzAvsCentrality);
  fListHist->Add(hCos2DeltaPhiVzCvsCentrality);
  fListHist->Add(hCos2DeltaPhiVzMvsCentrality);
  fListHist->Add(hCos2DeltaPhiTPCfvsCentrality);
  fListHist->Add(hCos2DeltaPhiTPCpvsCentrality);
  fListHist->Add(hCos2DeltaPhiTPCnvsCentrality);

  hEvPlaneTPCvsEvPVz05   = new TH2F("hEvPlaneTPCvsEvPVz05"  ,"hEvPlaneTPCvsEvPVz05"  ,100,-2*TMath::Pi(),2*TMath::Pi(),100,-2*TMath::Pi(),2*TMath::Pi());                      
  hEvPlaneTPCvsEvPVz075  = new TH2F("hEvPlaneTPCvsEvPVz075" ,"hEvPlaneTPCvsEvPVz075" ,100,-2*TMath::Pi(),2*TMath::Pi(),100,-2*TMath::Pi(),2*TMath::Pi()); 
  hEvPlaneTPCvsEvPVz1530 = new TH2F("hEvPlaneTPCvsEvPVz1530","hEvPlaneTPCvsEvPVz1530",100,-2*TMath::Pi(),2*TMath::Pi(),100,-2*TMath::Pi(),2*TMath::Pi());
  hEvPlaneTPCvsEvPVz3050 = new TH2F("hEvPlaneTPCvsEvPVz3050","hEvPlaneTPCvsEvPVz3050",100,-2*TMath::Pi(),2*TMath::Pi(),100,-2*TMath::Pi(),2*TMath::Pi());                      
  hEvPlaneTPCvsEvPVz2040 = new TH2F("hEvPlaneTPCvsEvPVz2040","hEvPlaneTPCvsEvPVz2040",100,-2*TMath::Pi(),2*TMath::Pi(),100,-2*TMath::Pi(),2*TMath::Pi());                      
  hEvPlaneTPCvsEvPVz4060 = new TH2F("hEvPlaneTPCvsEvPVz4060","hEvPlaneTPCvsEvPVz4060",100,-2*TMath::Pi(),2*TMath::Pi(),100,-2*TMath::Pi(),2*TMath::Pi());   
 
 
  fListHist->Add(hEvPlaneTPCvsEvPVz05);                      
  fListHist->Add(hEvPlaneTPCvsEvPVz075); 
  fListHist->Add(hEvPlaneTPCvsEvPVz1530);
  fListHist->Add(hEvPlaneTPCvsEvPVz3050);                      
  fListHist->Add(hEvPlaneTPCvsEvPVz2040);                      
  fListHist->Add(hEvPlaneTPCvsEvPVz4060);   

  hCos2DeltaPhivsPt075          = new TH2F("hCos2DeltaPhivsPt075",      "hCos2DeltaPhivsPt075"       ,75, -1.1,1.1,210,-5.5,5.5);
  hCos2DeltaPhiVZEROvsPt075     = new TH2F("hCos2DeltaPhiVZEROvsPt075", "hCos2DeltaPhiVZEROvsPt075"  ,75, -1.1,1.1,210,-5.5,5.5);
  hCos2DeltaPhivsPt1530         = new TH2F("hCos2DeltaPhivsPt1530",     "hCos2DeltaPhivsPt1530"      ,75, -1.1,1.1,210,-5.5,5.5);
  hCos2DeltaPhiVZEROvsPt1530    = new TH2F("hCos2DeltaPhiVZEROvsPt1530","hCos2DeltaPhiVZEROvsPt1530" ,75, -1.1,1.1,210,-5.5,5.5);
  hCos2DeltaPhivsPt3050         = new TH2F("hCos2DeltaPhivsPt3050",     "hCos2DeltaPhivsPt3050"      ,75, -1.1,1.1,210,-5.5,5.5);
  hCos2DeltaPhiVZEROvsPt3050    = new TH2F("hCos2DeltaPhiVZEROvsPt3050","hCos2DeltaPhiVZEROvsPt3050" ,75, -1.1,1.1,210,-5.5,5.5);
  hCos2DeltaPhivsPt05           = new TH2F("hCos2DeltaPhivsPt05",       "hCos2DeltaPhivsPt05"        ,75, -1.1,1.1,210,-5.5,5.5);
  hCos2DeltaPhiVZEROvsPt05      = new TH2F("hCos2DeltaPhiVZEROvsPt05",  "hCos2DeltaPhiVZEROvsPt05"   ,75, -1.1,1.1,210,-5.5,5.5);
  hCos2DeltaPhivsPt2040         = new TH2F("hCos2DeltaPhivsPt2040",     "hCos2DeltaPhivsPt2040"      ,75, -1.1,1.1,210,-5.5,5.5);
  hCos2DeltaPhiVZEROvsPt2040    = new TH2F("hCos2DeltaPhiVZEROvsPt2040","hCos2DeltaPhiVZEROvsPt2040" ,75, -1.1,1.1,210,-5.5,5.5);
  hCos2DeltaPhivsPt4060         = new TH2F("hCos2DeltaPhivsPt4060",     "hCos2DeltaPhivsPt4060"      ,75, -1.1,1.1,210,-5.5,5.5);
  hCos2DeltaPhiVZEROvsPt4060    = new TH2F("hCos2DeltaPhiVZEROvsPt4060","hCos2DeltaPhiVZEROvsPt4060" ,75, -1.1,1.1,210,-5.5,5.5);

  //--------------
  fListHist->Add(hCos2DeltaPhivsPt075);
  fListHist->Add(hCos2DeltaPhiVZEROvsPt075);
  fListHist->Add(hCos2DeltaPhivsPt1530);
  fListHist->Add(hCos2DeltaPhiVZEROvsPt1530);
  fListHist->Add(hCos2DeltaPhivsPt3050);
  fListHist->Add(hCos2DeltaPhiVZEROvsPt3050);
  fListHist->Add(hCos2DeltaPhivsPt05);    
  fListHist->Add(hCos2DeltaPhiVZEROvsPt05);    
  fListHist->Add(hCos2DeltaPhivsPt2040); 
  fListHist->Add(hCos2DeltaPhiVZEROvsPt2040);    	      
  fListHist->Add(hCos2DeltaPhivsPt4060);    	      
  fListHist->Add(hCos2DeltaPhiVZEROvsPt4060); 


  if(! fNtuple1 ) {

    fNtuple1 = new TTree("fNtuple1","fNtuple1");
 
    fNtuple1->Branch("tCentrality"   ,&tCentrality   ,"tCentrality[2]/D"    );
    fNtuple1->Branch("tPulls"        ,&tPulls        ,"tPulls[3]/D"   );
    fNtuple1->Branch("tMomentum"     ,&tMomentum     ,"tMomentum[4]/D"   );
    fNtuple1->Branch("tDCA" 	     ,&tDCA 	     ,"tDCA[2]/D" 	 );
    fNtuple1->Branch("tisTOF"        ,&tisTOF	     ,"tisTOF[2]/I"	 );
    fNtuple1->Branch("tTOFtrack"     ,&tTOFtrack     ,"tTOFtrack[3]/D"   );
    fNtuple1->Branch("tCharge"       ,&tCharge       ,"tCharge/I"        );
    fNtuple1->Branch("tPhi"          ,&tPhi          ,"tPhi/D"           );
    fNtuple1->Branch("trpangleTPC"   ,&trpangleTPC   ,"trpangleTPC/D"    );
    fNtuple1->Branch("trpangleVZERO" ,&trpangleVZERO ,"trpangleVZERO[3]/D"  );
 
    if(fDataType == "SIM"){
      fNtuple1->Branch("tPDGCode"      ,&tPDGCode      ,"tPDGCode/D"       );
      fNtuple1->Branch("tPDGCodeMum"   ,&tPDGCode      ,"tPDGCodeMum/D"    );
      fNtuple1->Branch("tIsPrimaryTr"  ,&tIsPrimaryTr  ,"tIsPrimaryTr/D"   );
      fNtuple1->Branch("tIsSecondaryTr",&tIsSecondaryTr,"tIsSecondaryTr[2]/D" );
    }
     
  }
  
  
  if(fDataType == "SIM"){
  
    if(! fNtuple2 ) {
      
      fNtuple2 = new TTree("fNtuple2","fNtuple2");
      
      fNtuple2->Branch("tCentralityMC"   ,&tCentralityMC   ,"tCentralityMC/D"    );
      fNtuple2->Branch("tVertexCoordMC"  ,&tVertexCoordMC  ,"tVertexCoordMC[3]/D");
      fNtuple2->Branch("tMomentumMC"     ,&tMomentumMC     ,"tMomentumMC[3]/D"   );
      fNtuple2->Branch("tPDGCodeMC"      ,&tPDGCodeMC      ,"tPDGCodeMC/D"       );
      fNtuple2->Branch("tPDGCodeMumMC"   ,&tPDGCodeMC      ,"tPDGCodeMumMC/D"    );
      fNtuple2->Branch("tIsPrimary"      ,&tIsPrimary      ,"tIsPrimary/D"       );
      fNtuple2->Branch("tIsSecondary"    ,&tIsSecondary    ,"tIsSecondary[2]/D"     );
      fNtuple2->Branch("tEtaMC"          ,&tEtaMC	   ,"tEtaMC/D"	       );
      fNtuple2->Branch("tPtMC"           ,&tPtMC	   ,"tPtMC/D"	       );
      fNtuple2->Branch("tYMC"            ,&tYMC	           ,"tYMC/D"	       );

    } 
    
  }
 
  PostData(1,  fListHist);
  PostData(2,  fNtuple1);
  PostData(3,  fNtuple2);

  fFlowEvent = new AliFlowEvent(10000);
  PostData(4, fFlowEvent);

}// end UserCreateOutputObjects


//====================== USER EXEC ========================

void AliAnalysisTaskNucleiv2::UserExec(Option_t *) 
{
  // Main loop
  // Called for EACH event
  //  cout<<"AliAnalysisTaskNucleiv2 Starting UserExec"<<endl;

  Info("AliAnalysisTaskNucleiv2","Starting UserExec");  
  
  AliVEvent *event = InputEvent();
  if (!event) { Printf("ERROR: Could not retrieve event"); return; }
  
  AliStack *stack = 0;
  
  if(fDataType == "SIM"){
    AliMCEvent *mcEvent = MCEvent();
    if (!mcEvent) { 
      Printf("ERROR: Could not retrieve MC event"); 
      return; 
    }
    Printf("MC particles: %d", mcEvent->GetNumberOfTracks());
    stack = mcEvent->Stack();
    if( !stack ) { 
      Printf( "Stack not available"); 
      return; 
    }
  }
  
  // create pointer to event
 
  if(fAnalysisType == "ESD"){
    
    AliESDEvent* lESDevent = dynamic_cast<AliESDEvent*>(event);
    if (!lESDevent) {
      AliError("Cannot get the ESD event");
      return;
    }  
    
    fHistEventMultiplicity->Fill(1);
    fHistEventMultiplicity->Fill(7);
  
    //_____________________________________________________
    //   Centrality  
    
    AliCentrality *centrality = lESDevent->GetCentrality();
    Float_t percentile=centrality->GetCentralityPercentile("V0M");
    
    Int_t TrackNumber = lESDevent->GetNumberOfTracks();
    fHistTrackMultiplicity->Fill(TrackNumber,percentile); //tracce per evento
    
    //______________________________________________________
    // PID
    
    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse=inputHandler->GetPIDResponse(); 
    
    //=================================================================
        
    Float_t impactXY=-999., impactZ=-999.;
    Double_t pinTPC=0.,poutTPC=0.,TPCSignal=0.;
    
    ULong_t  status=0;
    Bool_t   isTPC=kFALSE;
    
    
    // Primary vertex cut
  
    const AliESDVertex *vtx = lESDevent->GetPrimaryVertexTracks();
    
    if(vtx->GetNContributors()<1) {
      
      // SPD vertex cut
      vtx = lESDevent->GetPrimaryVertexSPD();
      
      if(vtx->GetNContributors()<1) {
	Info("AliAnalysisTaskHelium3Pi","No good vertex, skip event");
	return; // NO GOOD VERTEX, SKIP EVENT 
      }
    }
    
    fHistEventMultiplicity->Fill(2); // analyzed events with PV
    
    if(TMath::Abs(vtx->GetZv())>10) return;
    fHistEventMultiplicity->Fill(3);

    Bool_t isSelectedCentral     = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kCentral);
    Bool_t isSelectedSemiCentral = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kSemiCentral);
    Bool_t isSelectedMB          = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
    
    fHistTrackMultiplicity->Fill(TrackNumber,percentile); 
    
    Int_t eventtype = -999;
    
    //  cout<<"ET 1: "<<eventtype<<endl;
    
    if(fDataType == "REAL"){
   
      if(isSelectedCentral){
	// if(percentile<0)return;
	// if(percentile>=7.5)return;
	fHistEventMultiplicity->Fill(4);
	fHistTrackMultiplicityCentral->Fill(TrackNumber,percentile); 
	eventtype =1;
      }
      
      if(isSelectedSemiCentral){
	// if(percentile<15)return;
	// if(percentile>=50)return;
	fHistEventMultiplicity->Fill(5);
	fHistTrackMultiplicitySemiCentral->Fill(TrackNumber,percentile); 
	eventtype =2;
      }
      
      if(isSelectedMB){
	if(percentile<0)return;
	if(percentile>=80)return;
	fHistEventMultiplicity->Fill(6);
	fHistTrackMultiplicityMB->Fill(TrackNumber,percentile); 
	eventtype =3;
      }
    
      //    cout<<"ET 2: "<<eventtype<<endl;
      
      if(eventtype!=1 && eventtype!=2 && eventtype!=3 )return;
    }
    
    if(fDataType == "SIM"){
      cout<<"Take SIM event"<<endl;
      eventtype = -999;
      if(eventtype!=-999)return;
    }
    //----------------
    SetNullCuts(lESDevent);
    PrepareFlowEvent(TrackNumber,fFlowEvent);    //Calculate event plane Qvector and EP resolution for inclusive
  
    AliEventplane *pl=lESDevent->GetEventplane();
    
    if(fDataType == "REAL"){
      
      if(!pl ){
	AliError("AliAnalysisTaskSENucleiv2::UserExec:no eventplane! v2 analysis without eventplane not possible!\n");
	fHistEventMultiplicity->Fill(12);
      }
    }
     
   
    //Event plane from FLOW
    //=============================================V0EP from Alex======================================================================
    
    Double_t qxEPa = 0, qyEPa = 0;
    Double_t qxEPc = 0, qyEPc = 0;
    Double_t qxEP = 0, qyEP = 0;
    
    Double_t evPlAngV0A = pl->CalculateVZEROEventPlane(lESDevent, 8, 2, qxEPa, qyEPa);
    Double_t evPlAngV0C = pl->CalculateVZEROEventPlane(lESDevent, 9, 2, qxEPc, qyEPc);
    Double_t evPlAngV0  = pl->CalculateVZEROEventPlane(lESDevent,10, 2, qxEP,  qyEP);
    
    Double_t Qx2  = 0, Qy2  = 0;
    Double_t Qx2p = 0, Qy2p = 0;
    Double_t Qx2n = 0, Qy2n = 0;
    
    for (Int_t iT = 0; iT < TrackNumber; iT++){
      AliESDtrack* track = lESDevent->GetTrack(iT);
    
      if (!track)
	continue;
      
      if ((TMath::Abs(track->Eta()) > 0.8) || (track->Pt() < 0.2) || (track->GetTPCNcls() < 70) || (track->Pt() >= 20.0))
	continue;
      
      if(!fESDtrackCutsEP->AcceptTrack(track))
	continue;
      
      if(track->Eta()>0 && track->Eta()<0.8){
        
	Qx2p += TMath::Cos(2*track->Phi());
	Qy2p += TMath::Sin(2*track->Phi());
      }
      if(track->Eta()<0 && track->Eta()> -0.8){
	
	Qx2n += TMath::Cos(2*track->Phi());
	Qy2n += TMath::Sin(2*track->Phi());
      }
         
      Qx2 += TMath::Cos(2*track->Phi());
      Qy2 += TMath::Sin(2*track->Phi());
      
    }
    
    Double_t evPlAngTPC  = TMath::ATan2(Qy2, Qx2)/2.;
    Double_t evPlAngTPCn = TMath::ATan2(Qy2n, Qx2n)/2.;
    Double_t evPlAngTPCp = TMath::ATan2(Qy2p, Qx2p)/2.;
    
    EPVzAvsCentrality->Fill(evPlAngV0A,percentile);
    EPVzCvsCentrality->Fill(evPlAngV0C,percentile);
    EPTPCvsCentrality->Fill(evPlAngTPC,percentile);
    
    EPTPCnvsCentrality->Fill(evPlAngTPCn,percentile);
    EPTPCpvsCentrality->Fill(evPlAngTPCp,percentile);
    EPVzvsCentrality->Fill(evPlAngV0,percentile);
    
    // For TPC resolution
    // Cos(2*(tpc-v0A)) vs Centrality
    // Cos(2*(tpc-v0C)) vs Centrality
    // Cos(2*(v0A-v0C)) vs Centrality
    // For VZEROM resolution :
    // Cos(2*(V0-TPCp)) vs Centrality
    // Cos(2*(V0-TPCn)) vs Centrality
    // Cos(2*(TPCp-TPCn)) vs Centrality

    hCos2DeltaPhiVzAvsCentrality ->Fill(TMath::Cos(2.*(evPlAngV0-evPlAngTPCp))  ,percentile);
    hCos2DeltaPhiVzCvsCentrality ->Fill(TMath::Cos(2.*(evPlAngV0-evPlAngTPCn))  ,percentile);
    hCos2DeltaPhiVzMvsCentrality ->Fill(TMath::Cos(2.*(evPlAngTPCp-evPlAngTPCn)),percentile);
    hCos2DeltaPhiTPCfvsCentrality->Fill(TMath::Cos(2.*(evPlAngTPC-evPlAngV0A))  ,percentile);
    hCos2DeltaPhiTPCpvsCentrality->Fill(TMath::Cos(2.*(evPlAngTPC-evPlAngV0C))  ,percentile);
    hCos2DeltaPhiTPCnvsCentrality->Fill(TMath::Cos(2.*(evPlAngV0C-evPlAngV0A))  ,percentile);
    

    // e volendo 3 per il vzero come sotto

    if(percentile>=0 || percentile<=5){
      //This is v0 resolution 
      fSubEventDPhiv205->Fill(0.5, TMath::Cos(2.*(evPlAngV0A-evPlAngTPC))); // vzeroa - tpc
      fSubEventDPhiv205->Fill(1.5, TMath::Cos(2.*(evPlAngV0A-evPlAngV0C))); // vzeroa - vzeroc
      fSubEventDPhiv205->Fill(2.5, TMath::Cos(2.*(evPlAngV0C-evPlAngTPC))); // tpc - vzeroc
      
      
      fSubEventDPhiv2new05->Fill(0.5, TMath::Cos(2.*(evPlAngV0-evPlAngTPCp))  ); // vzero - tpcp
      fSubEventDPhiv2new05->Fill(1.5, TMath::Cos(2.*(evPlAngV0-evPlAngTPCn))  ); // vzero - tpcn
      fSubEventDPhiv2new05->Fill(2.5, TMath::Cos(2.*(evPlAngTPCp-evPlAngTPCn))); // tpcp - tpcn
      
      hEvPlaneTPCvsEvPVz05  ->Fill(  evPlAngTPC,evPlAngV0);    
    }
    
    if(percentile>=20 || percentile<=40){
      fSubEventDPhiv22040->Fill(0.5, TMath::Cos(2.*(evPlAngV0A-evPlAngTPC))); // vzeroa - tpc
      fSubEventDPhiv22040->Fill(1.5, TMath::Cos(2.*(evPlAngV0A-evPlAngV0C))); // vzeroa - vzeroc
      fSubEventDPhiv22040->Fill(2.5, TMath::Cos(2.*(evPlAngV0C-evPlAngTPC))); // tpc - vzeroc
      
      
      fSubEventDPhiv2new2040->Fill(0.5, TMath::Cos(2.*(evPlAngV0-evPlAngTPCp))); // vzero - tpcp
      fSubEventDPhiv2new2040->Fill(1.5, TMath::Cos(2.*(evPlAngV0-evPlAngTPCn))); // vzero - tpcn
      fSubEventDPhiv2new2040->Fill(2.5, TMath::Cos(2.*(evPlAngTPCp-evPlAngTPCn))); // tpcp - tpcn
      
      hEvPlaneTPCvsEvPVz2040   ->Fill(  evPlAngTPC,evPlAngV0); 
      
    }
    
    if(percentile>=40 || percentile<=60){
      fSubEventDPhiv24060->Fill(0.5, TMath::Cos(2.*(evPlAngV0A-evPlAngTPC))); // vzeroa - tpc
      fSubEventDPhiv24060->Fill(1.5, TMath::Cos(2.*(evPlAngV0A-evPlAngV0C))); // vzeroa - vzeroc
      fSubEventDPhiv24060->Fill(2.5, TMath::Cos(2.*(evPlAngV0C-evPlAngTPC))); // tpc - vzeroc
            
      fSubEventDPhiv2new4060->Fill(0.5, TMath::Cos(2.*(evPlAngV0-evPlAngTPCp))); // vzero - tpcp
      fSubEventDPhiv2new4060->Fill(1.5, TMath::Cos(2.*(evPlAngV0-evPlAngTPCn))); // vzero - tpcn
      fSubEventDPhiv2new4060->Fill(2.5, TMath::Cos(2.*(evPlAngTPCp-evPlAngTPCn))); // tpcp - tpcn

      hEvPlaneTPCvsEvPVz4060 ->Fill(  evPlAngTPC,evPlAngV0);
    }
  
    if(percentile>=0 || percentile<=7.5)               
      hEvPlaneTPCvsEvPVz075->Fill(  evPlAngTPC,evPlAngV0);
    if(percentile>=15 || percentile<=30)
      hEvPlaneTPCvsEvPVz1530->Fill(  evPlAngTPC,evPlAngV0);
    if(percentile>=30 || percentile<=50)
      hEvPlaneTPCvsEvPVz3050->Fill(  evPlAngTPC,evPlAngV0);
    
    //====================================================================================================================
    
    // To remove auto-correlation
    TVector2 *q = 0x0;
    TVector2 *qsub1=0x0;
    TVector2 *qsub2=0x0;
    qsub1 = pl->GetQsub1();
    qsub2 = pl->GetQsub2();
            
    q = pl->GetQVector();

    trpangleTPC      = evPlAngTPC;
    trpangleVZERO[0] = evPlAngV0;
    trpangleVZERO[1] = evPlAngV0A;
    trpangleVZERO[2] = evPlAngV0C;
    
    // cout<<"rpangle TPC: "<<rpangleTPC<<endl;
  
    Int_t isTOF=0;
    Int_t isoutTPC=0;
    //  cout<<"TRack number MC "<<TrackNumber<<endl;
    
    for (Int_t j=0; j<TrackNumber; j++) { //loop on tracks
      
      AliESDtrack *esdtrack=lESDevent->GetTrack(j);
      if (!fESDtrackCuts->AcceptTrack(esdtrack)) continue;
      
      status  = (ULong_t)esdtrack->GetStatus();
      
      isTPC    = (((status) & AliESDtrack::kTPCin)   != 0);
      isTOF    = ((((status) & AliESDtrack::kTOFout) != 0) && (((status) & AliESDtrack::kTIME) != 0));
      isoutTPC = (((status) & AliESDtrack::kTPCout)  != 0);
      
      TPCSignal=esdtrack->GetTPCsignal(); 
      
      if (TPCSignal<10)continue;
      if (TPCSignal>1000)continue;
      if (!isTPC)continue;
      
      if(!esdtrack->GetTPCInnerParam())continue;
      AliExternalTrackParam trackIn(*esdtrack->GetInnerParam()); 
      pinTPC= trackIn.GetP(); 

      fhBB->Fill(pinTPC*esdtrack->GetSign(),TPCSignal);

      if(isTOF){
	if(!esdtrack->GetOuterParam())continue;    
	AliExternalTrackParam trackOut(*esdtrack->GetOuterParam()); 
	poutTPC = trackOut.GetP();  
	fhTOF->Fill(poutTPC*esdtrack->GetSign(),(esdtrack->GetIntegratedLength()/esdtrack->GetTOFsignal())/2.99792458e-2);
    
      }

      Int_t   fIdxInt[200]; //dummy array
      Int_t   nClustersTPC = esdtrack->GetTPCclusters(fIdxInt);
      Float_t chi2PerClusterTPC = esdtrack->GetTPCchi2()/(Float_t)(nClustersTPC);

      esdtrack->GetImpactParameters(impactXY, impactZ);

      Float_t deutExp  = -999;
      Float_t tritExp  = -999;
      Float_t hel3Exp  = -999;
  
      if(fDataType == "REAL"){
	deutExp  = AliExternalTrackParam::BetheBlochAleph(pinTPC/(0.938*2),1.45802,27.4992,4.00313e-15,2.48485,8.31768);
	tritExp  = AliExternalTrackParam::BetheBlochAleph(pinTPC/(0.938*3),1.45802,27.4992,4.00313e-15,2.48485,8.31768);
	hel3Exp  = 4*AliExternalTrackParam::BetheBlochAleph(2*pinTPC/(0.938*3),1.74962,27.4992,4.00313e-15,2.42485,8.31768);
      }

      if(fDataType == "SIM"){
	Double_t parMC[5] = {1.17329, 27.4992, 4.00313e-15, 2.1204316, 4.1373729}; // NEW!!!
	deutExp = AliExternalTrackParam::BetheBlochAleph(pinTPC/(0.938*2),parMC[0],parMC[1],parMC[2],parMC[3],parMC[4]);
	tritExp = AliExternalTrackParam::BetheBlochAleph(pinTPC/(0.938*3),parMC[0],parMC[1],parMC[2],parMC[3],parMC[4]);
	hel3Exp = 4*AliExternalTrackParam::BetheBlochAleph(2*pinTPC/(0.938*3),parMC[0],parMC[1],parMC[2],parMC[3],parMC[4]);
      }
    
      Double_t pullTPC     = (TPCSignal - deutExp)/(0.07*deutExp);
      Double_t pullTPChel3 = (TPCSignal - hel3Exp)/(0.07*hel3Exp);
      Double_t pullTPCtrit = (TPCSignal - tritExp)/(0.07*tritExp);
    
      tPulls[0] = pullTPC;
      tPulls[1] = pullTPChel3;
      tPulls[2] = pullTPCtrit;

      //Fill the tree

      tCentrality[0] = percentile;
      tCentrality[1] = eventtype;

    
      tMomentum[0] = esdtrack->Px();
      tMomentum[1] = esdtrack->Py();
      tMomentum[2] = esdtrack->Pz();
    
      //Corrected momentum from Alexander
      Double_t pT = esdtrack->Pt()/(1 - 0.333303/TMath::Power(esdtrack->Pt() + 0.651111, 5.27268));
      tMomentum[3] = pT;
    
      tDCA[0]      = impactXY;
      tDCA[1]      = impactZ;

      tisTOF[0] = isTOF;
      tisTOF[1] = isoutTPC;
 
      Double_t p   = esdtrack->P();
      Double_t tof = esdtrack->GetTOFsignal()-fPIDResponse->GetTOFResponse().GetStartTime(p);
    
      tTOFtrack[0] = poutTPC;
      tTOFtrack[1] = tof;                             //ps = Time
      tTOFtrack[2] = esdtrack->GetIntegratedLength(); //cm
    
      tCharge =  esdtrack->Charge();
      tPhi    =  esdtrack->Phi();
    
      Float_t beta = 0;
      Float_t gamma = 0;
      Float_t deltaMass = 0;

      if(fDataType == "REAL"){
	if(TMath::Abs(pinTPC) < 6 && TMath::Abs(pullTPC) < 3){

	  fhBBDeu->Fill(pinTPC*esdtrack->GetSign(),TPCSignal);
	
	  if(pinTPC < 1.0 && fFillNtuple == kTRUE)
	    fNtuple1->Fill();
	
	  if(tTOFtrack[1] > 0 ){
	    beta = tTOFtrack[2]/(tTOFtrack[1] * 2.99792457999999984e-02);
	    gamma = 1/TMath::Sqrt(1 - beta*beta);
	    deltaMass = poutTPC/TMath::Sqrt(gamma*gamma - 1) - 1.8756;
	    fhMassTOF->Fill(deltaMass);
	  }
	
	  if(pinTPC >= 1.0 && tTOFtrack[1] > 0 && TMath::Abs(deltaMass)<0.5 && fFillNtuple == kTRUE )
	    fNtuple1->Fill();
	}
      
	if(pinTPC < 10. && TMath::Abs(pullTPChel3) < 3){
	  if(fFillNtuple == kTRUE)
	    fNtuple1->Fill();
	  fhBBDeu->Fill(pinTPC*esdtrack->GetSign(),TPCSignal);
	}
      
      
	//Fill also things for flow package

	if(TMath::Abs(pullTPC)<2){
	
	  fhPtDeu->Fill(esdtrack->Pt(),pT);
	
	  //Remove AutoCorrelation
	  trpangleTPC = GetEventPlaneForCandidate(esdtrack,q,pl,qsub1,qsub2);

	  Float_t deltaphiTPC=2*GetPhi0Pi(tPhi-trpangleTPC);
	  Float_t deltaphiV0 =2*GetPhi0Pi(tPhi-trpangleVZERO[0]);
	
	
	  if(pinTPC < 1.0){
	  
	    //Here

	    //==================================================================
	    //----------------------Flow of Inclusive Particles --------------------------------------------------------
	 
	    AliFlowTrack *sTrack = new AliFlowTrack();
	    sTrack->Set(esdtrack);
	    sTrack->SetID(esdtrack->GetID());
	    sTrack->SetForRPSelection(kTRUE);
	    sTrack->SetForPOISelection(kTRUE);
	  
	    for(int iRPs=0; iRPs!=fFlowEvent->NumberOfTracks(); ++iRPs)
	      {
		//   cout << " no of rps " << iRPs << endl;
		AliFlowTrack *iRP = dynamic_cast<AliFlowTrack*>(fFlowEvent->GetTrack( iRPs ));
		if (!iRP) continue;
		if (!iRP->InRPSelection()) continue;
		if( sTrack->GetID() == iRP->GetID())
		  {
		    if(fDebug) printf(" was in RP set");
		    //       cout << sTrack->GetID() <<"   ==  " << iRP->GetID() << " was in RP set" <<endl;
		    iRP->SetForRPSelection(kFALSE);
		    // fFlowEvent->SetNumberOfRPs(fFlowEvent->GetNumberOfRPs() - 1);
		  }
	      } //end of for loop on RPs
	    fFlowEvent->InsertTrack(((AliFlowTrack*) sTrack));
	    fFlowEvent->SetNumberOfPOIs(fFlowEvent->GetNumberOfPOIs()+1);
	    

	    //=================================================================

	    //  if(tCharge > 0){
	  
	    if(tCentrality[0]>0 && tCentrality[0]<7.5){
	      hCos2DeltaPhivsPt075           ->Fill(TMath::Cos(deltaphiTPC),tMomentum[3]*tCharge);	       
	      hCos2DeltaPhiVZEROvsPt075      ->Fill(TMath::Cos(deltaphiV0) ,tMomentum[3]*tCharge);       
	    
	    }

	    if(tCentrality[0]>0 && tCentrality[0]<5){
	      hCos2DeltaPhivsPt05           ->Fill(TMath::Cos(deltaphiTPC),tMomentum[3]*tCharge);	       
	      hCos2DeltaPhiVZEROvsPt05      ->Fill(TMath::Cos(deltaphiV0) ,tMomentum[3]*tCharge);       
	    }
	    
	    if(tCentrality[0]>15 && tCentrality[0]<30){
	      hCos2DeltaPhivsPt1530           ->Fill(TMath::Cos(deltaphiTPC),tMomentum[3]*tCharge);
	      hCos2DeltaPhiVZEROvsPt1530      ->Fill(TMath::Cos(deltaphiV0) ,tMomentum[3]*tCharge);
	    }
	
	    if(tCentrality[0]>20 && tCentrality[0]<40){
	      hCos2DeltaPhivsPt2040           ->Fill(TMath::Cos(deltaphiTPC),tMomentum[3]*tCharge);
	      hCos2DeltaPhiVZEROvsPt2040      ->Fill(TMath::Cos(deltaphiV0) ,tMomentum[3]*tCharge);
	    }
	  
	    if(tCentrality[0]>30 && tCentrality[0]<50){
	      hCos2DeltaPhivsPt3050           ->Fill(TMath::Cos(deltaphiTPC),tMomentum[3]*tCharge);     
	      hCos2DeltaPhiVZEROvsPt3050      ->Fill(TMath::Cos(deltaphiV0) ,tMomentum[3]*tCharge);        
	    }
	    
	    if(tCentrality[0]>40 && tCentrality[0]<60){
	      hCos2DeltaPhivsPt4060           ->Fill(TMath::Cos(deltaphiTPC),tMomentum[3]*tCharge);
	      hCos2DeltaPhiVZEROvsPt4060      ->Fill(TMath::Cos(deltaphiV0) ,tMomentum[3]*tCharge);
	    } 
	  }
	  
	  if(pinTPC > 1.0 && pinTPC  < 6.0 && tTOFtrack[1] > 0 && TMath::Abs(deltaMass)<0.5){
	  
		  
	    if(tCentrality[0]>0 && tCentrality[0]<7.5){
	      hCos2DeltaPhivsPt075           ->Fill(TMath::Cos(deltaphiTPC),tMomentum[3]*tCharge);	       
	      hCos2DeltaPhiVZEROvsPt075      ->Fill(TMath::Cos(deltaphiV0) ,tMomentum[3]*tCharge);       
	    
	    }

	    if(tCentrality[0]>0 && tCentrality[0]<5){
	      hCos2DeltaPhivsPt05           ->Fill(TMath::Cos(deltaphiTPC),tMomentum[3]*tCharge);	       
	      hCos2DeltaPhiVZEROvsPt05      ->Fill(TMath::Cos(deltaphiV0) ,tMomentum[3]*tCharge);       
	    }
	    
	    if(tCentrality[0]>15 && tCentrality[0]<30){
	      hCos2DeltaPhivsPt1530           ->Fill(TMath::Cos(deltaphiTPC),tMomentum[3]*tCharge);
	      hCos2DeltaPhiVZEROvsPt1530      ->Fill(TMath::Cos(deltaphiV0) ,tMomentum[3]*tCharge);
	    }
	
	    if(tCentrality[0]>20 && tCentrality[0]<40){
	      hCos2DeltaPhivsPt2040           ->Fill(TMath::Cos(deltaphiTPC),tMomentum[3]*tCharge);
	      hCos2DeltaPhiVZEROvsPt2040      ->Fill(TMath::Cos(deltaphiV0) ,tMomentum[3]*tCharge);
	    }
	  
	    if(tCentrality[0]>30 && tCentrality[0]<50){
	      hCos2DeltaPhivsPt3050           ->Fill(TMath::Cos(deltaphiTPC),tMomentum[3]*tCharge);     
	      hCos2DeltaPhiVZEROvsPt3050      ->Fill(TMath::Cos(deltaphiV0) ,tMomentum[3]*tCharge);        
	    }
	  
	    if(tCentrality[0]>40 && tCentrality[0]<60){
	      hCos2DeltaPhivsPt4060           ->Fill(TMath::Cos(deltaphiTPC),tMomentum[3]*tCharge);
	      hCos2DeltaPhiVZEROvsPt4060      ->Fill(TMath::Cos(deltaphiV0) ,tMomentum[3]*tCharge);
	    } 
	  }
	}
      }
        
      if(fDataType == "SIM"){
      
	Int_t  label = TMath::Abs(esdtrack->GetLabel());
	TParticle * part = stack->Particle(label);
	Int_t PDGCode=part->GetPdgCode();
	
	Int_t motherPDG=0;
      
	Int_t mumid = part->GetFirstMother();
	if(mumid>-1){
	  TParticle *mother=(TParticle*)stack->Particle(mumid);
	  motherPDG = mother->GetPdgCode();
	}
      
	if( PDGCode == 1000010020 || PDGCode == -1000010020 || PDGCode == 1000020030 || PDGCode == -1000020030){
	
	  tPDGCode=PDGCode;
	  tPDGCodeMum = motherPDG;
	
	  tIsPrimaryTr      = stack->IsPhysicalPrimary(label);
	  tIsSecondaryTr[0] = stack->IsSecondaryFromMaterial(label);
	  tIsSecondaryTr[1] = stack->IsSecondaryFromWeakDecay(label);
	
	  fNtuple1->Fill();
	  fhPtDeu->Fill(esdtrack->Pt(),pT);
	
	  if(tTOFtrack[1] > 0){
	    beta = tTOFtrack[2]/(tTOFtrack[1] * 2.99792457999999984e-02);
	    gamma = 1/TMath::Sqrt(1 - beta*beta);
	    fhMassTOF->Fill(poutTPC/TMath::Sqrt(gamma*gamma - 1) - 1.8756);
	  }
	}
      }
    
    }   //track
  
    //==END RECONSTRUCTION==
    
    
    // MC truth

    if(fDataType == "SIM"){
      
      for (Int_t iMC=0; iMC<stack->GetNtrack(); iMC++){
	
	const TParticle *tparticle = stack->Particle(iMC);
	Long_t PDGCode = tparticle->GetPdgCode();
	
	Double_t eta = tparticle->Eta();
	Double_t pt  = tparticle->Pt();
	Double_t rap = tparticle->Y();
      
	//check which particle it is 
	
	Float_t codemoth = 0;
	
	Int_t indexMoth=tparticle->GetFirstMother();
	
	if(indexMoth>=0){
	  TParticle* moth = stack->Particle(indexMoth);
	  codemoth = TMath::Abs(moth->GetPdgCode());
	}
	
	//d, 3He
	
	if( PDGCode == 1000010020 || PDGCode == -1000010020 || PDGCode == 1000020030 || PDGCode == -1000020030){
	  
	
	  
	  tCentralityMC = percentile;
	  
	  tVertexCoordMC[0] = vtx->GetXv();
	  tVertexCoordMC[1] = vtx->GetYv();
	  tVertexCoordMC[2] = vtx->GetZv();
	  
	  tVertexCoordMC[0] = tparticle->Px();
	  tVertexCoordMC[1] = tparticle->Py();
	  tVertexCoordMC[2] = tparticle->Pz();
	  
	  tPDGCodeMC = PDGCode;
	  tPDGCodeMumMC = codemoth;
	  tEtaMC = eta;
	  
	  tPtMC = pt;
	  tYMC  = rap;
	  
	  tIsPrimary      = stack->IsPhysicalPrimary(iMC);
	  tIsSecondary[0] = stack->IsSecondaryFromMaterial(iMC);
	  tIsSecondary[1] = stack->IsSecondaryFromWeakDecay(iMC);
	  
	  
	  fNtuple2->Fill();
	}
      }
      
    }
  }
  
  PostData(1, fListHist);
  PostData(2,  fNtuple1);
  PostData(3,  fNtuple2);
  
  PostData(4, fFlowEvent);

} //end userexec


//________________________________________________________________________

void AliAnalysisTaskNucleiv2::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
}

