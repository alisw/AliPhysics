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

ClassImp(AliAnalysisTaskNucleiv2)

using std::cout;
using std::endl;
    
//________________________________________________________________________
AliAnalysisTaskNucleiv2::AliAnalysisTaskNucleiv2() 
: AliAnalysisTaskSE(), 
  fAnalysisType("ESD"), 
  fCollidingSystems(0), 
  fDataType("REAL"),
  fFillNtuple(kFALSE),
  fListHistCascade(0), 
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
  hRPangleTPCvsCentrality(0),           //RESOLUTION Histrograms
  hPlaneResoTPCvsCentrality(0),
  hRPangleVZEROvsCentrality(0),
  hRPangleVZEROAvsCentrality(0),
  hRPangleVZEROCvsCentrality(0),
  hPlaneResoVZEROvsCentrality(0),
  hPlaneResoVZEROAvsCentrality(0),
  hPlaneResoVZEROCvsCentrality(0),
  hCosPhivsPt(0),                         
  hSinPhivsPt(0),                         
  hPhivsPt(0),                         
  hAntiCosPhivsPt(0),                     
  hAntiSinPhivsPt(0),                     
  hAntiPhivsPt(0),                     
  hCosDeltaPhivsPt075(0),                      
  hSinDeltaPhivsPt075(0),                      
  hDeltaPhivsPt075(0),                      
  hCosDeltaPhiVZEROvsPt075(0),                 	      
  hSinDeltaPhiVZEROvsPt075(0),                 	      
  hDeltaPhiVZEROvsPt075(0),                 
  hCosDeltaPhivsPt1530(0),                     
  hSinDeltaPhivsPt1530(0),                     
  hDeltaPhivsPt1530(0),                     
  hCosDeltaPhiVZEROvsPt1530(0),                
  hSinDeltaPhiVZEROvsPt1530(0),                
  hDeltaPhiVZEROvsPt1530(0),                
  hCosDeltaPhivsPt3050(0),                     
  hSinDeltaPhivsPt3050(0),                     
  hDeltaPhivsPt3050(0),                     
  hCosDeltaPhiVZEROvsPt3050(0),                
  hSinDeltaPhiVZEROvsPt3050(0),                
  hDeltaPhiVZEROvsPt3050(0),                
  hCosDeltaPhivsPt1550(0),                     
  hSinDeltaPhivsPt1550(0),                     
  hDeltaPhivsPt1550(0),                     
  hCosDeltaPhiVZEROvsPt1550(0),                
  hSinDeltaPhiVZEROvsPt1550(0),                
  hDeltaPhiVZEROvsPt1550(0),                
  hAntiCosDeltaPhivsPt075(0),                  
  hAntiSinDeltaPhivsPt075(0),                  
  hAntiDeltaPhivsPt075(0),                  
  hAntiCosDeltaPhiVZEROvsPt075(0),             
  hAntiSinDeltaPhiVZEROvsPt075(0),             
  hAntiDeltaPhiVZEROvsPt075(0),             
  hAntiCosDeltaPhivsPt1530(0),                 
  hAntiSinDeltaPhivsPt1530(0),                 
  hAntiDeltaPhivsPt1530(0),                 
  hAntiCosDeltaPhiVZEROvsPt1530(0),            
  hAntiSinDeltaPhiVZEROvsPt1530(0),            
  hAntiDeltaPhiVZEROvsPt1530(0),            
  hAntiCosDeltaPhivsPt3050(0),                 
  hAntiSinDeltaPhivsPt3050(0),                 
  hAntiDeltaPhivsPt3050(0),                 
  hAntiCosDeltaPhiVZEROvsPt3050(0),            
  hAntiSinDeltaPhiVZEROvsPt3050(0),            
  hAntiDeltaPhiVZEROvsPt3050(0),            
  hAntiCosDeltaPhivsPt1550(0),                 
  hAntiSinDeltaPhivsPt1550(0),                 
  hAntiDeltaPhivsPt1550(0),                 
  hAntiCosDeltaPhiVZEROvsPt1550(0),            
  hAntiSinDeltaPhiVZEROvsPt1550(0),            
  hAntiDeltaPhiVZEROvsPt1550(0),     
  fESDtrackCuts(0),
  fPIDResponse(0),
  fNtuple1(0),
  tCentrality(0),
  tTPCMomentum(0),      
  tdEdx(0),             
  tEta(0),              
  tITSclustermap(0),    
  tCharge(0),           
  tPtCorr(0),           
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
  //
  Initialize();
}

//________________________________________________________________________
AliAnalysisTaskNucleiv2::AliAnalysisTaskNucleiv2(const char *name,const char *datatype,Bool_t filltree) 
  : AliAnalysisTaskSE(name), 
    fAnalysisType("ESD"), 
    fCollidingSystems(0), 
    fDataType(datatype),
    fFillNtuple(filltree),
    fListHistCascade(0), 
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
    hRPangleTPCvsCentrality(0),           //RESOLUTION Histrograms
    hPlaneResoTPCvsCentrality(0),
    hRPangleVZEROvsCentrality(0),
    hRPangleVZEROAvsCentrality(0),
    hRPangleVZEROCvsCentrality(0),
    hPlaneResoVZEROvsCentrality(0),
    hPlaneResoVZEROAvsCentrality(0),
    hPlaneResoVZEROCvsCentrality(0),
    hCosPhivsPt(0),                         
    hSinPhivsPt(0),                         
    hPhivsPt(0),                         
    hAntiCosPhivsPt(0),                     
    hAntiSinPhivsPt(0),                     
    hAntiPhivsPt(0),                     
    hCosDeltaPhivsPt075(0),                      
    hSinDeltaPhivsPt075(0),                      
    hDeltaPhivsPt075(0),                      
    hCosDeltaPhiVZEROvsPt075(0),                 	      
    hSinDeltaPhiVZEROvsPt075(0),                 	      
    hDeltaPhiVZEROvsPt075(0),                 
    hCosDeltaPhivsPt1530(0),                     
    hSinDeltaPhivsPt1530(0),                     
    hDeltaPhivsPt1530(0),                     
    hCosDeltaPhiVZEROvsPt1530(0),                
    hSinDeltaPhiVZEROvsPt1530(0),                
    hDeltaPhiVZEROvsPt1530(0),                
    hCosDeltaPhivsPt3050(0),                     
    hSinDeltaPhivsPt3050(0),                     
    hDeltaPhivsPt3050(0),                     
    hCosDeltaPhiVZEROvsPt3050(0),                
    hSinDeltaPhiVZEROvsPt3050(0),                
    hDeltaPhiVZEROvsPt3050(0),                
    hCosDeltaPhivsPt1550(0),                     
    hSinDeltaPhivsPt1550(0),                     
    hDeltaPhivsPt1550(0),                     
    hCosDeltaPhiVZEROvsPt1550(0),                
    hSinDeltaPhiVZEROvsPt1550(0),                
    hDeltaPhiVZEROvsPt1550(0),                
    hAntiCosDeltaPhivsPt075(0),                  
    hAntiSinDeltaPhivsPt075(0),                  
    hAntiDeltaPhivsPt075(0),                  
    hAntiCosDeltaPhiVZEROvsPt075(0),             
    hAntiSinDeltaPhiVZEROvsPt075(0),             
    hAntiDeltaPhiVZEROvsPt075(0),             
    hAntiCosDeltaPhivsPt1530(0),                 
    hAntiSinDeltaPhivsPt1530(0),                 
    hAntiDeltaPhivsPt1530(0),                 
    hAntiCosDeltaPhiVZEROvsPt1530(0),            
    hAntiSinDeltaPhiVZEROvsPt1530(0),            
    hAntiDeltaPhiVZEROvsPt1530(0),            
    hAntiCosDeltaPhivsPt3050(0),                 
    hAntiSinDeltaPhivsPt3050(0),                 
    hAntiDeltaPhivsPt3050(0),                 
    hAntiCosDeltaPhiVZEROvsPt3050(0),            
    hAntiSinDeltaPhiVZEROvsPt3050(0),            
    hAntiDeltaPhiVZEROvsPt3050(0),            
    hAntiCosDeltaPhivsPt1550(0),                 
    hAntiSinDeltaPhivsPt1550(0),                 
    hAntiDeltaPhivsPt1550(0),                 
    hAntiCosDeltaPhiVZEROvsPt1550(0),            
    hAntiSinDeltaPhiVZEROvsPt1550(0),            
    hAntiDeltaPhiVZEROvsPt1550(0),     
    fESDtrackCuts(0),
    fPIDResponse(0),
    fNtuple1(0),
    tCentrality(0),
    tTPCMomentum(0),      
    tdEdx(0),             
    tEta(0),              
    tITSclustermap(0),    
    tCharge(0),           
    tPtCorr(0),           
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
  // Output slot #0 writes into a TList container (Cascade)

  //
  // create track cuts
  //
  fESDtrackCuts = new AliESDtrackCuts("AliESDtrackCuts","AliESDtrackCuts");
  //
  Initialize();
  SetDataType(datatype);
  SetFillTree(filltree);

  DefineInput(0, TChain::Class());
  
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());

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
 
  //  Printf("Initizialize\n");
 
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


//==================DEFINITION OF OUTPUT OBJECTS==============================

void AliAnalysisTaskNucleiv2::UserCreateOutputObjects()
{
  fListHistCascade = new TList();
  fListHistCascade->SetOwner();  // IMPORTANT!

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
    fListHistCascade->Add(fHistEventMultiplicity);
  }

  if(! fHistTrackMultiplicity ){
    fHistTrackMultiplicity  = new TH2F( "fHistTrackMultiplicity", "Nb of Tracks MB Events |Vz| < 10", 25000,0, 25000,105,-0.5,104.5);
    fHistTrackMultiplicity->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicity->GetYaxis()->SetTitle("Percentile");
    fListHistCascade->Add(fHistTrackMultiplicity);
  } 

  if(! fHistTrackMultiplicityCentral ){
    fHistTrackMultiplicityCentral  = new TH2F( "fHistTrackMultiplicityCentral", "Nb of Tracks MB Events |Vz| < 10", 25000,0, 25000,105,-0.5,104.5);
    fHistTrackMultiplicityCentral->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicityCentral->GetYaxis()->SetTitle("Percentile");
    fListHistCascade->Add(fHistTrackMultiplicityCentral);
  } 
  if(! fHistTrackMultiplicitySemiCentral ){
    fHistTrackMultiplicitySemiCentral  = new TH2F( "fHistTrackMultiplicitySemiCentral", "Nb of Tracks MB Events |Vz| < 10", 25000,0, 25000,105,-0.5,104.5);
    fHistTrackMultiplicitySemiCentral->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicitySemiCentral->GetYaxis()->SetTitle("Percentile");
    fListHistCascade->Add(fHistTrackMultiplicitySemiCentral);
  } 
  if(! fHistTrackMultiplicityMB ){
    fHistTrackMultiplicityMB  = new TH2F( "fHistTrackMultiplicityMB", "Nb of Tracks MB Events |Vz| < 10", 25000,0, 25000,105,-0.5,104.5);
    fHistTrackMultiplicityMB->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicityMB->GetYaxis()->SetTitle("Percentile");
    fListHistCascade->Add(fHistTrackMultiplicityMB);
  } 
 
  if(! fhBB ){
    fhBB = new TH2F( "fhBB" , "BetheBlochTPC" , 240,-6,6,1000,0,1000);
    fListHistCascade->Add(fhBB);
  }
  
  if(! fhBBDeu ){
    fhBBDeu = new TH2F( "fhBBDeu" , "BetheBlochTPC - Deuteron" , 240,-6,6,1000,0,1000);
    fListHistCascade->Add(fhBBDeu);
  }
 
  if(!fhPtDeu  ){
    fhPtDeu = new TH2F( "fhPtDeu" , "pt corretto vs pt track - Deuteron" , 120,0,6,120,0,6);
    fListHistCascade->Add(fhPtDeu);
  }

  if(! fhTOF ){
    fhTOF = new TH2F( "fhTOF" , "Scatter Plot TOF" , 240,-6,6,1000,0,1.2);
    fListHistCascade->Add(fhTOF);
  }
  if(! fhMassTOF){
    fhMassTOF=new TH1F ("fhMassTOF","Particle Mass - TOF", 300,-5 ,5);
    fListHistCascade->Add(fhMassTOF);
  }
  
  if(!hRPangleTPCvsCentrality){
    hRPangleTPCvsCentrality     = new TH2F("hRPangleTPCvsCentrality"    ,"hRPangleTPCvsCentrality"    ,100,0,TMath::Pi(),210,-0.5,105.5);
    fListHistCascade->Add(hRPangleTPCvsCentrality);
  }
  if(!hPlaneResoTPCvsCentrality){
    hPlaneResoTPCvsCentrality   = new TH2F("hPlaneResoTPCvsCentrality"  ,"hPlaneResoTPCvsCentrality"  ,500,-1,1,210,-0.5,105.5);
    fListHistCascade->Add(hPlaneResoTPCvsCentrality);
  }
  if(!hRPangleVZEROvsCentrality){
    hRPangleVZEROvsCentrality   = new TH2F("hRPangleVZEROvsCentrality"  ,"hRPangleVZEROvsCentrality"  ,100,0,TMath::Pi(),210,-0.5,105.5);
    fListHistCascade->Add(hRPangleVZEROvsCentrality);
  }
  
  if(!hRPangleVZEROAvsCentrality){
    hRPangleVZEROAvsCentrality   = new TH2F("hRPangleVZEROAvsCentrality"  ,"hRPangleVZEROAvsCentrality"  ,100,0,TMath::Pi(),210,-0.5,105.5);
    fListHistCascade->Add(hRPangleVZEROAvsCentrality);
  }

  if(!hRPangleVZEROCvsCentrality){
    hRPangleVZEROCvsCentrality   = new TH2F("hRPangleVZEROCvsCentrality"  ,"hRPangleVZEROCvsCentrality"  ,100,0,TMath::Pi(),210,-0.5,105.5);
    fListHistCascade->Add(hRPangleVZEROCvsCentrality);
  }

  if(!hPlaneResoVZEROvsCentrality){
    hPlaneResoVZEROvsCentrality = new TH2F("hPlaneResoVZEROvsCentrality","hPlaneResoVZEROvsCentrality",500,-1,1,210,-0.5,105.5);
    fListHistCascade->Add(hPlaneResoVZEROvsCentrality);
  }

  if(!hPlaneResoVZEROAvsCentrality){
    hPlaneResoVZEROAvsCentrality = new TH2F("hPlaneResoVZEROAvsCentrality","hPlaneResoVZEROAvsCentrality",500,-1,1,210,-0.5,105.5);
    fListHistCascade->Add(hPlaneResoVZEROAvsCentrality);
  }

  if(!hPlaneResoVZEROCvsCentrality){
    hPlaneResoVZEROCvsCentrality = new TH2F("hPlaneResoVZEROCvsCentrality","hPlaneResoVZEROCvsCentrality",500,-1,1,210,-0.5,105.5);
    fListHistCascade->Add(hPlaneResoVZEROCvsCentrality);
  }
 
  hCosPhivsPt      = new TH2F("hCosPhivsPt","hCosPhivsPt",200,-1,1,60,0,3.0);	     
  hSinPhivsPt      = new TH2F("hSinPhivsPt","hSinPhivsPt",100, 0,1,60,0,3.0);	     
  hPhivsPt         = new TH2F("hPhivsPt","hPhivsPt"      ,36, 0,2*TMath::Pi(),60,0,3.0);
  
  hAntiCosPhivsPt      = new TH2F("hAntiCosPhivsPt","hAntiCosPhivsPt",200,-1,1,60,0,3.0);	     
  hAntiSinPhivsPt      = new TH2F("hAntiSinPhivsPt","hAntiSinPhivsPt",100, 0,1,60,0,3.0);	     
  hAntiPhivsPt         = new TH2F("hAntiPhivsPt","hAntiPhivsPt"      ,36, 0,2*TMath::Pi(),60,0,3.0);
  

  hCosDeltaPhivsPt075 = new TH2F("hCosDeltaPhivsPt075","hCosDeltaPhivsPt075",200,-1,1,60,0,3.0);
  hSinDeltaPhivsPt075 = new TH2F("hSinDeltaPhivsPt075","hSinDeltaPhivsPt075",100, 0,1,60,0,3.0);
  hDeltaPhivsPt075    = new TH2F("hDeltaPhivsPt075","hDeltaPhivsPt075"      ,36, 0,2*TMath::Pi(),60,0,3.0);
  
  hCosDeltaPhiVZEROvsPt075 = new TH2F("hCosDeltaPhiVZEROvsPt075","hCosDeltaPhiVZEROvsPt075",200,-1,1,60,0,3.0);	      
  hSinDeltaPhiVZEROvsPt075 = new TH2F("hSinDeltaPhiVZEROvsPt075","hSinDeltaPhiVZEROvsPt075",100, 0,1,60,0,3.0);	      
  hDeltaPhiVZEROvsPt075    = new TH2F("hDeltaPhiVZEROvsPt075","hDeltaPhiVZEROvsPt075"      ,36, 0,2*TMath::Pi(),60,0,3.0);

  hCosDeltaPhivsPt1530 = new TH2F("hCosDeltaPhivsPt1530","hCosDeltaPhivsPt1530",200,-1,1,60,0,3.0);
  hSinDeltaPhivsPt1530 = new TH2F("hSinDeltaPhivsPt1530","hSinDeltaPhivsPt1530",100, 0,1,60,0,3.0);
  hDeltaPhivsPt1530    = new TH2F("hDeltaPhivsPt1530","hDeltaPhivsPt1530"      ,36, 0,2*TMath::Pi(),60,0,3.0);
  
  hCosDeltaPhiVZEROvsPt1530 = new TH2F("hCosDeltaPhiVZEROvsPt1530","hCosDeltaPhiVZEROvsPt1530",200,-1,1,60,0,3.0);	      
  hSinDeltaPhiVZEROvsPt1530 = new TH2F("hSinDeltaPhiVZEROvsPt1530","hSinDeltaPhiVZEROvsPt1530",100, 0,1,60,0,3.0);	      
  hDeltaPhiVZEROvsPt1530    = new TH2F("hDeltaPhiVZEROvsPt1530","hDeltaPhiVZEROvsPt1530"      ,36, 0,2*TMath::Pi(),60,0,3.0);

  hCosDeltaPhivsPt3050 = new TH2F("hCosDeltaPhivsPt3050","hCosDeltaPhivsPt3050",200,-1,1,60,0,3.0);
  hSinDeltaPhivsPt3050 = new TH2F("hSinDeltaPhivsPt3050","hSinDeltaPhivsPt3050",100, 0,1,60,0,3.0);
  hDeltaPhivsPt3050    = new TH2F("hDeltaPhivsPt3050","hDeltaPhivsPt3050"      ,36, 0,2*TMath::Pi(),60,0,3.0);
  
  hCosDeltaPhiVZEROvsPt3050 = new TH2F("hCosDeltaPhiVZEROvsPt3050","hCosDeltaPhiVZEROvsPt3050",200,-1,1,60,0,3.0);	      
  hSinDeltaPhiVZEROvsPt3050 = new TH2F("hSinDeltaPhiVZEROvsPt3050","hSinDeltaPhiVZEROvsPt3050",100, 0,1,60,0,3.0);	      
  hDeltaPhiVZEROvsPt3050    = new TH2F("hDeltaPhiVZEROvsPt3050","hDeltaPhiVZEROvsPt3050"      ,36, 0,2*TMath::Pi(),60,0,3.0);

  hCosDeltaPhivsPt1550 = new TH2F("hCosDeltaPhivsPt1550","hCosDeltaPhivsPt1550",200,-1,1,60,0,3.0);
  hSinDeltaPhivsPt1550 = new TH2F("hSinDeltaPhivsPt1550","hSinDeltaPhivsPt1550",100, 0,1,60,0,3.0);
  hDeltaPhivsPt1550    = new TH2F("hDeltaPhivsPt1550","hDeltaPhivsPt1550"      ,36, 0,2*TMath::Pi(),60,0,3.0);
  
  hCosDeltaPhiVZEROvsPt1550 = new TH2F("hCosDeltaPhiVZEROvsPt1550","hCosDeltaPhiVZEROvsPt1550",200,-1,1,60,0,3.0);	      
  hSinDeltaPhiVZEROvsPt1550 = new TH2F("hSinDeltaPhiVZEROvsPt1550","hSinDeltaPhiVZEROvsPt1550",100, 0,1,60,0,3.0);	      
  hDeltaPhiVZEROvsPt1550    = new TH2F("hDeltaPhiVZEROvsPt1550","hDeltaPhiVZEROvsPt1550"      ,36, 0,2*TMath::Pi(),60,0,3.0);

  
  hAntiCosDeltaPhivsPt075 = new TH2F("hAntiCosDeltaPhivsPt075","hAntiCosDeltaPhivsPt075",200,-1,1,60,0,3.0);
  hAntiSinDeltaPhivsPt075 = new TH2F("hAntiSinDeltaPhivsPt075","hAntiSinDeltaPhivsPt075",100, 0,1,60,0,3.0);
  hAntiDeltaPhivsPt075    = new TH2F("hAntiDeltaPhivsPt075","hAntiDeltaPhivsPt075"      ,36, 0,2*TMath::Pi(),60,0,3.0);
  
  hAntiCosDeltaPhiVZEROvsPt075 = new TH2F("hAntiCosDeltaPhiVZEROvsPt075","hAntiCosDeltaPhiVZEROvsPt075",200,-1,1,60,0,3.0);	      
  hAntiSinDeltaPhiVZEROvsPt075 = new TH2F("hAntiSinDeltaPhiVZEROvsPt075","hAntiSinDeltaPhiVZEROvsPt075",100, 0,1,60,0,3.0);	      
  hAntiDeltaPhiVZEROvsPt075    = new TH2F("hAntiDeltaPhiVZEROvsPt075","hAntiDeltaPhiVZEROvsPt075"      ,36, 0,2*TMath::Pi(),60,0,3.0);

  hAntiCosDeltaPhivsPt1530 = new TH2F("hAntiCosDeltaPhivsPt1530","hAntiCosDeltaPhivsPt1530",200,-1,1,60,0,3.0);
  hAntiSinDeltaPhivsPt1530 = new TH2F("hAntiSinDeltaPhivsPt1530","hAntiSinDeltaPhivsPt1530",100, 0,1,60,0,3.0);
  hAntiDeltaPhivsPt1530    = new TH2F("hAntiDeltaPhivsPt1530","hAntiDeltaPhivsPt1530"      ,36, 0,2*TMath::Pi(),60,0,3.0);
  
  hAntiCosDeltaPhiVZEROvsPt1530 = new TH2F("hAntiCosDeltaPhiVZEROvsPt1530","hAntiCosDeltaPhiVZEROvsPt1530",200,-1,1,60,0,3.0);	      
  hAntiSinDeltaPhiVZEROvsPt1530 = new TH2F("hAntiSinDeltaPhiVZEROvsPt1530","hAntiSinDeltaPhiVZEROvsPt1530",100, 0,1,60,0,3.0);	      
  hAntiDeltaPhiVZEROvsPt1530    = new TH2F("hAntiDeltaPhiVZEROvsPt1530","hAntiDeltaPhiVZEROvsPt1530"      ,36, 0,2*TMath::Pi(),60,0,3.0);

  hAntiCosDeltaPhivsPt3050 = new TH2F("hAntiCosDeltaPhivsPt3050","hAntiCosDeltaPhivsPt3050",200,-1,1,60,0,3.0);
  hAntiSinDeltaPhivsPt3050 = new TH2F("hAntiSinDeltaPhivsPt3050","hAntiSinDeltaPhivsPt3050",100, 0,1,60,0,3.0);
  hAntiDeltaPhivsPt3050    = new TH2F("hAntiDeltaPhivsPt3050","hAntiDeltaPhivsPt3050"      ,36, 0,2*TMath::Pi(),60,0,3.0);
  
  hAntiCosDeltaPhiVZEROvsPt3050 = new TH2F("hAntiCosDeltaPhiVZEROvsPt3050","hAntiCosDeltaPhiVZEROvsPt3050",200,-1,1,60,0,3.0);	      
  hAntiSinDeltaPhiVZEROvsPt3050 = new TH2F("hAntiSinDeltaPhiVZEROvsPt3050","hAntiSinDeltaPhiVZEROvsPt3050",100, 0,1,60,0,3.0);	      
  hAntiDeltaPhiVZEROvsPt3050    = new TH2F("hAntiDeltaPhiVZEROvsPt3050","hAntiDeltaPhiVZEROvsPt3050"      ,36, 0,2*TMath::Pi(),60,0,3.0);

  hAntiCosDeltaPhivsPt1550 = new TH2F("hAntiCosDeltaPhivsPt1550","hAntiCosDeltaPhivsPt1550",200,-1,1,60,0,3.0);
  hAntiSinDeltaPhivsPt1550 = new TH2F("hAntiSinDeltaPhivsPt1550","hAntiSinDeltaPhivsPt1550",100, 0,1,60,0,3.0);
  hAntiDeltaPhivsPt1550    = new TH2F("hAntiDeltaPhivsPt1550","hAntiDeltaPhivsPt1550"      ,36, 0,2*TMath::Pi(),60,0,3.0);
  
  hAntiCosDeltaPhiVZEROvsPt1550 = new TH2F("hAntiCosDeltaPhiVZEROvsPt1550","hAntiCosDeltaPhiVZEROvsPt1550",200,-1,1,60,0,3.0);	      
  hAntiSinDeltaPhiVZEROvsPt1550 = new TH2F("hAntiSinDeltaPhiVZEROvsPt1550","hAntiSinDeltaPhiVZEROvsPt1550",100, 0,1,60,0,3.0);	      
  hAntiDeltaPhiVZEROvsPt1550    = new TH2F("hAntiDeltaPhiVZEROvsPt1550","hAntiDeltaPhiVZEROvsPt1550"      ,36, 0,2*TMath::Pi(),60,0,3.0);

  //--------------
  fListHistCascade->Add(hCosPhivsPt);                         
  fListHistCascade->Add(hSinPhivsPt);                         
  fListHistCascade->Add(hPhivsPt);                         
  fListHistCascade->Add(hAntiCosPhivsPt);                     
  fListHistCascade->Add(hAntiSinPhivsPt);                     
  fListHistCascade->Add(hAntiPhivsPt);                     
  fListHistCascade->Add(hCosDeltaPhivsPt075);                      
  fListHistCascade->Add(hSinDeltaPhivsPt075);                      
  fListHistCascade->Add(hDeltaPhivsPt075);                      
  fListHistCascade->Add(hCosDeltaPhiVZEROvsPt075);                 	      
  fListHistCascade->Add(hSinDeltaPhiVZEROvsPt075);                 	      
  fListHistCascade->Add(hDeltaPhiVZEROvsPt075);                 
  fListHistCascade->Add(hCosDeltaPhivsPt1530);                     
  fListHistCascade->Add(hSinDeltaPhivsPt1530);                     
  fListHistCascade->Add(hDeltaPhivsPt1530);                     
  fListHistCascade->Add(hCosDeltaPhiVZEROvsPt1530);                
  fListHistCascade->Add(hSinDeltaPhiVZEROvsPt1530);                
  fListHistCascade->Add(hDeltaPhiVZEROvsPt1530);                
  fListHistCascade->Add(hCosDeltaPhivsPt3050);                     
  fListHistCascade->Add(hSinDeltaPhivsPt3050);                     
  fListHistCascade->Add(hDeltaPhivsPt3050);                     
  fListHistCascade->Add(hCosDeltaPhiVZEROvsPt3050);                
  fListHistCascade->Add(hSinDeltaPhiVZEROvsPt3050);                
  fListHistCascade->Add(hDeltaPhiVZEROvsPt3050);                
  fListHistCascade->Add(hCosDeltaPhivsPt1550);                     
  fListHistCascade->Add(hSinDeltaPhivsPt1550);                     
  fListHistCascade->Add(hDeltaPhivsPt1550);                     
  fListHistCascade->Add(hCosDeltaPhiVZEROvsPt1550);                
  fListHistCascade->Add(hSinDeltaPhiVZEROvsPt1550);                
  fListHistCascade->Add(hDeltaPhiVZEROvsPt1550);                
  fListHistCascade->Add(hAntiCosDeltaPhivsPt075);                  
  fListHistCascade->Add(hAntiSinDeltaPhivsPt075);                  
  fListHistCascade->Add(hAntiDeltaPhivsPt075);                  
  fListHistCascade->Add(hAntiCosDeltaPhiVZEROvsPt075);             
  fListHistCascade->Add(hAntiSinDeltaPhiVZEROvsPt075);             
  fListHistCascade->Add(hAntiDeltaPhiVZEROvsPt075);             
  fListHistCascade->Add(hAntiCosDeltaPhivsPt1530);                 
  fListHistCascade->Add(hAntiSinDeltaPhivsPt1530);                 
  fListHistCascade->Add(hAntiDeltaPhivsPt1530);                 
  fListHistCascade->Add(hAntiCosDeltaPhiVZEROvsPt1530);            
  fListHistCascade->Add(hAntiSinDeltaPhiVZEROvsPt1530);            
  fListHistCascade->Add(hAntiDeltaPhiVZEROvsPt1530);            
  fListHistCascade->Add(hAntiCosDeltaPhivsPt3050);                 
  fListHistCascade->Add(hAntiSinDeltaPhivsPt3050);                 
  fListHistCascade->Add(hAntiDeltaPhivsPt3050);                 
  fListHistCascade->Add(hAntiCosDeltaPhiVZEROvsPt3050);            
  fListHistCascade->Add(hAntiSinDeltaPhiVZEROvsPt3050);            
  fListHistCascade->Add(hAntiDeltaPhiVZEROvsPt3050);            
  fListHistCascade->Add(hAntiCosDeltaPhivsPt1550);                 
  fListHistCascade->Add(hAntiSinDeltaPhivsPt1550);                 
  fListHistCascade->Add(hAntiDeltaPhivsPt1550);                 
  fListHistCascade->Add(hAntiCosDeltaPhiVZEROvsPt1550);            
  fListHistCascade->Add(hAntiSinDeltaPhiVZEROvsPt1550);            
  fListHistCascade->Add(hAntiDeltaPhiVZEROvsPt1550);            

  if(! fNtuple1 ) {

    fNtuple1 = new TTree("fNtuple1","fNtuple1");
 
    fNtuple1->Branch("tEventNumber"  ,&tEventNumber  ,"tEventNumber[7]/D");
    fNtuple1->Branch("tCentrality"   ,&tCentrality   ,"tCentrality/D"    );
    fNtuple1->Branch("tVertexCoord"  ,&tVertexCoord  ,"tVertexCoord[3]/D");
    fNtuple1->Branch("tPIDITS"       ,&tPIDITS       ,"tPIDITS[9]/D"     );
    fNtuple1->Branch("tPIDTPC"       ,&tPIDTPC       ,"tPIDTPC[9]/D"     );
    fNtuple1->Branch("tPIDTOF"       ,&tPIDTOF       ,"tPIDTOF[9]/D"     );
    fNtuple1->Branch("tPulls"        ,&tPulls        ,"tPulls[2]/D"   );
    fNtuple1->Branch("tMomentum"     ,&tMomentum     ,"tMomentum[3]/D"   );
    fNtuple1->Branch("tTPCMomentum"  ,&tTPCMomentum  ,"tTPCMomentum/D"   );
    fNtuple1->Branch("tdEdx"	     ,&tdEdx	     ,"tdEdx/D"	         );
    fNtuple1->Branch("tEta"	     ,&tEta	     ,"tEta/D"	         );
    fNtuple1->Branch("tDCA" 	     ,&tDCA 	     ,"tDCA[2]/D" 	 );
    fNtuple1->Branch("tTracksTPC"    ,&tTracksTPC    ,"tTracksTPC[2]/D"  );
    fNtuple1->Branch("tITSclustermap",&tITSclustermap,"tITSclustermap"   );
    fNtuple1->Branch("tITSsample"    ,&tITSsample    ,"tITSsample[4]/D"  );
    fNtuple1->Branch("tisTOF[2]"     ,&tisTOF	     ,"tisTOF[2]/I"	 );
    fNtuple1->Branch("tTOFtrack"     ,&tTOFtrack     ,"tTOFtrack[3]/D"   );
    fNtuple1->Branch("tCharge"       ,&tCharge       ,"tCharge/I"        );
    fNtuple1->Branch("tPtCorr"       ,&tPtCorr       ,"tPtCorr/D"        );
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
      
      fNtuple2->Branch("tEventNumberMC"  ,&tEventNumberMC  ,"tEventNumberMC[6]/D");
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

  PostData(1,  fListHistCascade);
  PostData(2,  fNtuple1);
  PostData(3,  fNtuple2);
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
  
  //!*********************!//
  //!  Define variables   !//
  //!*********************!//
  
  Double_t evNumber     = 0.;
  Double_t runNumber    = 0.;
  Double_t BCNumber     = 0.;
  Double_t OrbitNumber  = 0.;
  Double_t PeriodNumber = 0.;

  
  Double_t ITSsample[4];
  for(Int_t i=0;i<4;i++)ITSsample[i]=0;
  
  Double_t xPrimaryVertex=0.,yPrimaryVertex=0.,zPrimaryVertex=0.;
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
  
  xPrimaryVertex=vtx->GetXv();
  yPrimaryVertex=vtx->GetYv();
  zPrimaryVertex=vtx->GetZv();  

  if(TMath::Abs(zPrimaryVertex)>10) return;
  fHistEventMultiplicity->Fill(3);

  Bool_t isSelectedCentral     = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kCentral);
  Bool_t isSelectedSemiCentral = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kSemiCentral);
  Bool_t isSelectedMB          = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
  
  fHistTrackMultiplicity->Fill(TrackNumber,percentile); 
  
  Int_t eventtype = -999;
  
  //  cout<<"ET 1: "<<eventtype<<endl;

  if(fDataType == "REAL"){
   
    if(isSelectedCentral){
      if(percentile<0)return;
      if(percentile>=7.5)return;
      fHistEventMultiplicity->Fill(4);
      fHistTrackMultiplicityCentral->Fill(TrackNumber,percentile); 
      eventtype =1;
    }
    
    if(isSelectedSemiCentral){
      if(percentile<15)return;
      if(percentile>=50)return;
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
    //cout<<"ET MC: "<<eventtype<<endl;
    if(eventtype!=-999)return;
      
  }
  
  evNumber    = lESDevent->GetEventNumberInFile();
  runNumber   = lESDevent->GetRunNumber(); 
  BCNumber    = lESDevent->GetBunchCrossNumber();
  OrbitNumber = lESDevent->GetOrbitNumber();
  PeriodNumber= lESDevent->GetPeriodNumber();
  
  AliEventplane *pl=lESDevent->GetEventplane();
  
  if(fDataType == "REAL"){
    
    if(!pl ){
      AliError("AliAnalysisTaskSENucleiv2::UserExec:no eventplane! v2 analysis without eventplane not possible!\n");
      fHistEventMultiplicity->Fill(12);
    }
  }
  
  //Event plane
  
  Double_t rpangleTPC    =0;
  Double_t rpangleVZERO  =0;

  Double_t rpangleeventATPC =0;
  Double_t rpangleeventBTPC =0;
  Double_t deltaPsiTPC      =0;
  Double_t planeresoTPC     =0;

  Double_t rpangleeventAVZERO =0;
  Double_t rpangleeventBVZERO =0;
  Double_t rpangleeventCVZERO =0;
  Double_t deltaPsiVZERO      =0;
  Double_t planeresoVZERO     =0;

  //For candidate removal from TPC EP

  TVector2 *qsub1=0x0;
  TVector2 *qsub2=0x0;

  rpangleTPC   = pl->GetEventplane("Q");
  
  if(fDataType == "REAL"){
    if(rpangleTPC<0){
      fHistEventMultiplicity->Fill(11);
      return;
    }
  }
  
  //TPC resolution 
   
  hRPangleTPCvsCentrality->Fill(rpangleTPC,percentile);

  qsub1 = pl->GetQsub1();
  qsub2 = pl->GetQsub2();

  if(fDataType == "REAL"){  
    if(!qsub1 || !qsub2){
      AliError("AliAnalysisTaskSENucleiv2::UserExec:no qsub1 or qsub2!\n");
      return;
    }
 
    //TPC event Plane
    
    rpangleeventATPC = qsub1->Phi()/2.;
    rpangleeventBTPC = qsub2->Phi()/2.;
    
    deltaPsiTPC =rpangleeventATPC-rpangleeventBTPC;
    
    if(TMath::Abs(deltaPsiTPC)>TMath::Pi()/2.){
      if(deltaPsiTPC>0.) deltaPsiTPC-=TMath::Pi();
      else deltaPsiTPC +=TMath::Pi();
    } // difference of subevents reaction plane angle cannot be bigger than phi/2
    
    planeresoTPC = TMath::Cos(2.*deltaPsiTPC); // reaction plane resolution
    
    hPlaneResoTPCvsCentrality->Fill(planeresoTPC,percentile);
  
    //VZERO event plane
  
    rpangleVZERO = GetPhi0Pi(pl->GetEventplane("V0",lESDevent,2));
    rpangleeventBVZERO=GetPhi0Pi(pl->GetEventplane("V0A",lESDevent,2));
    rpangleeventCVZERO=GetPhi0Pi(pl->GetEventplane("V0C",lESDevent,2));

    hRPangleVZEROvsCentrality->Fill(rpangleVZERO,percentile);
    hRPangleVZEROAvsCentrality->Fill(rpangleeventBVZERO,percentile);
    hRPangleVZEROCvsCentrality->Fill(rpangleeventCVZERO,percentile);

    //Resolution V0 : V0M - V0A
    rpangleeventAVZERO = rpangleVZERO;
    deltaPsiVZERO =rpangleeventAVZERO-rpangleeventBVZERO;

    if(TMath::Abs(deltaPsiVZERO)>TMath::Pi()/2.){
      if(deltaPsiVZERO>0.) deltaPsiVZERO-=TMath::Pi();
      else deltaPsiVZERO +=TMath::Pi();
    } // difference of subevents reaction plane angle cannot be bigger than phi/2
  
    planeresoVZERO = TMath::Cos(2.*deltaPsiVZERO);

    hPlaneResoVZEROvsCentrality->Fill(planeresoVZERO,percentile);

    //Resolution V0 : V0M - V0C

    deltaPsiVZERO =rpangleeventAVZERO-rpangleeventCVZERO;

    if(TMath::Abs(deltaPsiVZERO)>TMath::Pi()/2.){
      if(deltaPsiVZERO>0.) deltaPsiVZERO-=TMath::Pi();
      else deltaPsiVZERO +=TMath::Pi();
    } // difference of subevents reaction plane angle cannot be bigger than phi/2
  
    planeresoVZERO = TMath::Cos(2.*deltaPsiVZERO);

    hPlaneResoVZEROAvsCentrality->Fill(planeresoVZERO,percentile);
    
    //Resolution V0 : V0A - V0C

    deltaPsiVZERO =rpangleeventBVZERO-rpangleeventCVZERO;

    if(TMath::Abs(deltaPsiVZERO)>TMath::Pi()/2.){
      if(deltaPsiVZERO>0.) deltaPsiVZERO-=TMath::Pi();
      else deltaPsiVZERO +=TMath::Pi();
    } // difference of subevents reaction plane angle cannot be bigger than phi/2
  
    planeresoVZERO = TMath::Cos(2.*deltaPsiVZERO);

    hPlaneResoVZEROCvsCentrality->Fill(planeresoVZERO,percentile);
  
    //rpangleeventBVZERO

    if(TMath::Abs(rpangleTPC-rpangleVZERO)>10)return;
  
    trpangleTPC      = rpangleTPC;
    trpangleVZERO[0] = rpangleVZERO;
    trpangleVZERO[1] = rpangleeventBVZERO;
    trpangleVZERO[2] = rpangleeventCVZERO;
  }
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
    Float_t hel3Exp  = -999;
  
    if(fDataType == "REAL"){
      deutExp  = AliExternalTrackParam::BetheBlochAleph(pinTPC/(0.938*2),1.45802,27.4992,4.00313e-15,2.48485,8.31768);
      hel3Exp = 4*AliExternalTrackParam::BetheBlochAleph(2*pinTPC/(0.938*3),1.74962,27.4992,4.00313e-15,2.42485,8.31768);
    }

    if(fDataType == "SIM"){
      Double_t parMC[5] = {1.17329, 27.4992, 4.00313e-15, 2.1204316, 4.1373729}; // NEW!!!
      deutExp = AliExternalTrackParam::BetheBlochAleph(pinTPC/(0.938*2),parMC[0],parMC[1],parMC[2],parMC[3],parMC[4]);
      hel3Exp = 4*AliExternalTrackParam::BetheBlochAleph(2*pinTPC/(0.938*3),parMC[0],parMC[1],parMC[2],parMC[3],parMC[4]);
    }
    
    Double_t pullTPC = (TPCSignal - deutExp)/(0.07*deutExp);
    Double_t pullTPChel3 = (TPCSignal - hel3Exp)/(0.07*hel3Exp);
  
    tPulls[0] = pullTPC;
    tPulls[1] = pullTPChel3;

    //Fill the tree

    tEventNumber[0] = evNumber    ;
    tEventNumber[1] = runNumber   ;
    tEventNumber[2] = BCNumber    ;
    tEventNumber[3] = OrbitNumber ;
    tEventNumber[4] = PeriodNumber;
    tEventNumber[5] = TrackNumber;
    tEventNumber[6] = eventtype ;

    tCentrality = percentile;
        
    tVertexCoord[0] = xPrimaryVertex;
    tVertexCoord[1] = yPrimaryVertex;
    tVertexCoord[2] = zPrimaryVertex;

    if(TMath::Abs(zPrimaryVertex)>10)continue;

    // bbtheo = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)0,esdtrack,(AliPID::EParticleType) 7);

    tPIDITS[0] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)0,esdtrack,(AliPID::EParticleType) 0);
    tPIDITS[1] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)0,esdtrack,(AliPID::EParticleType) 1);
    tPIDITS[2] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)0,esdtrack,(AliPID::EParticleType) 2);
    tPIDITS[3] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)0,esdtrack,(AliPID::EParticleType) 3);
    tPIDITS[4] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)0,esdtrack,(AliPID::EParticleType) 4);
    tPIDITS[5] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)0,esdtrack,(AliPID::EParticleType) 5);
    tPIDITS[6] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)0,esdtrack,(AliPID::EParticleType) 6);
    tPIDITS[7] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)0,esdtrack,(AliPID::EParticleType) 7);
    tPIDITS[8] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)0,esdtrack,(AliPID::EParticleType) 8);

    tPIDTPC[0] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)1,esdtrack,(AliPID::EParticleType) 0);
    tPIDTPC[1] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)1,esdtrack,(AliPID::EParticleType) 1);
    tPIDTPC[2] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)1,esdtrack,(AliPID::EParticleType) 2);
    tPIDTPC[3] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)1,esdtrack,(AliPID::EParticleType) 3);
    tPIDTPC[4] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)1,esdtrack,(AliPID::EParticleType) 4);
    tPIDTPC[5] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)1,esdtrack,(AliPID::EParticleType) 5);
    tPIDTPC[6] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)1,esdtrack,(AliPID::EParticleType) 6);
    tPIDTPC[7] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)1,esdtrack,(AliPID::EParticleType) 7);
    tPIDTPC[8] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)1,esdtrack,(AliPID::EParticleType) 8);

    tPIDTOF[0] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)3,esdtrack,(AliPID::EParticleType) 0);
    tPIDTOF[1] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)3,esdtrack,(AliPID::EParticleType) 1);
    tPIDTOF[2] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)3,esdtrack,(AliPID::EParticleType) 2);
    tPIDTOF[3] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)3,esdtrack,(AliPID::EParticleType) 3);
    tPIDTOF[4] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)3,esdtrack,(AliPID::EParticleType) 4);
    tPIDTOF[5] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)3,esdtrack,(AliPID::EParticleType) 5);
    tPIDTOF[6] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)3,esdtrack,(AliPID::EParticleType) 6);
    tPIDTOF[7] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)3,esdtrack,(AliPID::EParticleType) 7);
    tPIDTOF[8] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)3,esdtrack,(AliPID::EParticleType) 8);

    tMomentum[0] = esdtrack->Px();
    tMomentum[1] = esdtrack->Py();
    tMomentum[2] = esdtrack->Pz();
    
    tTPCMomentum = pinTPC;
    tEta         = esdtrack->Eta();
    tdEdx        = TPCSignal;
    tDCA[0]      = impactXY;
    tDCA[1]      = impactZ;

    // if(nClustersTPC<80)
    //   cout<<"!!!!! TPC cls: "<<nClustersTPC<<endl;

    tTracksTPC[0] = chi2PerClusterTPC;
    tTracksTPC[1] = nClustersTPC;

    tITSclustermap = esdtrack->GetITSClusterMap();
    esdtrack->GetITSdEdxSamples(ITSsample);
    tITSsample[0] = ITSsample[0];
    tITSsample[1] = ITSsample[1];
    tITSsample[2] = ITSsample[2];
    tITSsample[3] = ITSsample[3];

    tisTOF[0] = isTOF;
    tisTOF[1] = isoutTPC;
 
    Double_t p = esdtrack->P();
    Double_t tof = esdtrack->GetTOFsignal()-fPIDResponse->GetTOFResponse().GetStartTime(p);

    tTOFtrack[0] = poutTPC;
    tTOFtrack[1] = tof;    //esdtrack->GetTOFsignal();         //ps = Time
    tTOFtrack[2] = esdtrack->GetIntegratedLength(); //cm
    
    tCharge =  esdtrack->Charge();
    tPhi    =  esdtrack->Phi();

    //Corrected momentum from Alexander

    Double_t pT = esdtrack->Pt()/(1 - 0.333303/TMath::Power(esdtrack->Pt() + 0.651111, 5.27268));
    tPtCorr = pT;

    if(fDataType == "REAL"){
      if(pinTPC < 3. && TMath::Abs(pullTPC) < 3){
	if(fFillNtuple == kTRUE)
	  fNtuple1->Fill();
	fhBBDeu->Fill(pinTPC*esdtrack->GetSign(),TPCSignal);
      }
      
      if(pinTPC < 10. && TMath::Abs(pullTPChel3) < 3){
	if(fFillNtuple == kTRUE)
	  fNtuple1->Fill();
	fhBBDeu->Fill(pinTPC*esdtrack->GetSign(),TPCSignal);
      }
      
      if(TMath::Abs(pullTPC)<2){
	
	fhPtDeu->Fill(esdtrack->Pt(),pT);
	
	Float_t deltaphiTPC=2*GetPhi0Pi(tPhi-trpangleTPC);
	Float_t deltaphiV0 =2*GetPhi0Pi(tPhi-trpangleVZERO[0]);
	
	
	if(tTPCMomentum < 1.0){
	  
	  if(tCharge > 0){
	    
	    hCosPhivsPt           ->Fill(TMath::Cos(tPhi),tPtCorr);
	    hSinPhivsPt           ->Fill(TMath::Sin(tPhi),tPtCorr);
	    hPhivsPt              ->Fill(tPhi,tPtCorr);
	    
	    if(tCentrality>0 && tCentrality<7.5){
	      hCosDeltaPhivsPt075        ->Fill(TMath::Cos(deltaphiTPC),tPtCorr);   
	      hSinDeltaPhivsPt075        ->Fill(TMath::Sin(deltaphiTPC),tPtCorr);   
	      hDeltaPhivsPt075           ->Fill(deltaphiTPC,tPtCorr);	       
	      hCosDeltaPhiVZEROvsPt075   ->Fill(TMath::Cos(deltaphiV0),tPtCorr);    
	      hSinDeltaPhiVZEROvsPt075   ->Fill(TMath::Sin(deltaphiV0),tPtCorr);    
	      hDeltaPhiVZEROvsPt075      ->Fill(deltaphiV0,tPtCorr);       
	      
	    }
	    
	    if(tCentrality>15 && tCentrality<30){
	      hCosDeltaPhivsPt1530        ->Fill(TMath::Cos(deltaphiTPC),tPtCorr);   
	      hSinDeltaPhivsPt1530        ->Fill(TMath::Sin(deltaphiTPC),tPtCorr);   
	      hDeltaPhivsPt1530           ->Fill(deltaphiTPC,tPtCorr);	       
	      hCosDeltaPhiVZEROvsPt1530   ->Fill(TMath::Cos(deltaphiV0),tPtCorr);    
	      hSinDeltaPhiVZEROvsPt1530   ->Fill(TMath::Sin(deltaphiV0),tPtCorr);    
	      hDeltaPhiVZEROvsPt1530      ->Fill(deltaphiV0,tPtCorr);         
	    }
	
	    if(tCentrality>30 && tCentrality<50){
	      hCosDeltaPhivsPt3050        ->Fill(TMath::Cos(deltaphiTPC),tPtCorr);   
	      hSinDeltaPhivsPt3050        ->Fill(TMath::Sin(deltaphiTPC),tPtCorr);   
	      hDeltaPhivsPt3050           ->Fill(deltaphiTPC,tPtCorr);	       
	      hCosDeltaPhiVZEROvsPt3050   ->Fill(TMath::Cos(deltaphiV0),tPtCorr);    
	      hSinDeltaPhiVZEROvsPt3050   ->Fill(TMath::Sin(deltaphiV0),tPtCorr);    
	      hDeltaPhiVZEROvsPt3050      ->Fill(deltaphiV0,tPtCorr);             
	    }
	    
	    if(tCentrality>15 && tCentrality<50){
	      hCosDeltaPhivsPt1550        ->Fill(TMath::Cos(deltaphiTPC),tPtCorr);   
	      hSinDeltaPhivsPt1550        ->Fill(TMath::Sin(deltaphiTPC),tPtCorr);   
	      hDeltaPhivsPt1550           ->Fill(deltaphiTPC,tPtCorr);	       
	      hCosDeltaPhiVZEROvsPt1550   ->Fill(TMath::Cos(deltaphiV0),tPtCorr);    
	      hSinDeltaPhiVZEROvsPt1550   ->Fill(TMath::Sin(deltaphiV0),tPtCorr);    
	      hDeltaPhiVZEROvsPt1550      ->Fill(deltaphiV0,tPtCorr); 
	    }
	  }
	  
	  if(tCharge < 0){
	    
	    hAntiCosPhivsPt           ->Fill(TMath::Cos(tPhi),tPtCorr);
	    hAntiSinPhivsPt           ->Fill(TMath::Sin(tPhi),tPtCorr);
	    hAntiPhivsPt              ->Fill(tPhi,tPtCorr);
	    if(tCentrality>0 && tCentrality<7.5){
	      hAntiCosDeltaPhivsPt075        ->Fill(TMath::Cos(deltaphiTPC),tPtCorr);   
	      hAntiSinDeltaPhivsPt075        ->Fill(TMath::Sin(deltaphiTPC),tPtCorr);   
	      hAntiDeltaPhivsPt075           ->Fill(deltaphiTPC,tPtCorr);	       
	      hAntiCosDeltaPhiVZEROvsPt075   ->Fill(TMath::Cos(deltaphiV0),tPtCorr);    
	      hAntiSinDeltaPhiVZEROvsPt075   ->Fill(TMath::Sin(deltaphiV0),tPtCorr);    
	      hAntiDeltaPhiVZEROvsPt075      ->Fill(deltaphiV0,tPtCorr); 
	    }

	    if(tCentrality>15 && tCentrality<30){
	      hAntiCosDeltaPhivsPt1530        ->Fill(TMath::Cos(deltaphiTPC),tPtCorr);   
	      hAntiSinDeltaPhivsPt1530        ->Fill(TMath::Sin(deltaphiTPC),tPtCorr);   
	      hAntiDeltaPhivsPt1530           ->Fill(deltaphiTPC,tPtCorr);	       
	      hAntiCosDeltaPhiVZEROvsPt1530   ->Fill(TMath::Cos(deltaphiV0),tPtCorr);    
	      hAntiSinDeltaPhiVZEROvsPt1530   ->Fill(TMath::Sin(deltaphiV0),tPtCorr);    
	      hAntiDeltaPhiVZEROvsPt1530      ->Fill(deltaphiV0,tPtCorr);         
	    }
	  
	    if(tCentrality>30 && tCentrality<50){
	      hAntiCosDeltaPhivsPt3050        ->Fill(TMath::Cos(deltaphiTPC),tPtCorr);   
	      hAntiSinDeltaPhivsPt3050        ->Fill(TMath::Sin(deltaphiTPC),tPtCorr);   
	      hAntiDeltaPhivsPt3050           ->Fill(deltaphiTPC,tPtCorr);	       
	      hAntiCosDeltaPhiVZEROvsPt3050   ->Fill(TMath::Cos(deltaphiV0),tPtCorr);    
	      hAntiSinDeltaPhiVZEROvsPt3050   ->Fill(TMath::Sin(deltaphiV0),tPtCorr);    
	      hAntiDeltaPhiVZEROvsPt3050      ->Fill(deltaphiV0,tPtCorr);                
	    }
	  
	    if(tCentrality>15 && tCentrality<50){
	      hAntiCosDeltaPhivsPt1550        ->Fill(TMath::Cos(deltaphiTPC),tPtCorr);   
	      hAntiSinDeltaPhivsPt1550        ->Fill(TMath::Sin(deltaphiTPC),tPtCorr);   
	      hAntiDeltaPhivsPt1550           ->Fill(deltaphiTPC,tPtCorr);	       
	      hAntiCosDeltaPhiVZEROvsPt1550   ->Fill(TMath::Cos(deltaphiV0),tPtCorr);    
	      hAntiSinDeltaPhiVZEROvsPt1550   ->Fill(TMath::Sin(deltaphiV0),tPtCorr);    
	      hAntiDeltaPhiVZEROvsPt1550      ->Fill(deltaphiV0,tPtCorr);                
	    }
	  
	  }
	  
	}

	if(tTPCMomentum > 1.0 && tTPCMomentum < 3.0){
	  if(tPIDTOF[5]>2 || tPIDTOF[5]< -2)continue;
	
	  if(tCharge > 0){
	    
	    hCosPhivsPt           ->Fill(TMath::Cos(tPhi),tPtCorr);
	    hSinPhivsPt           ->Fill(TMath::Sin(tPhi),tPtCorr);
	    hPhivsPt              ->Fill(tPhi,tPtCorr);
	    
	    
	    if(tCentrality>0 && tCentrality<7.5){
	      hCosDeltaPhivsPt075        ->Fill(TMath::Cos(deltaphiTPC),tPtCorr);   
	      hSinDeltaPhivsPt075        ->Fill(TMath::Sin(deltaphiTPC),tPtCorr);   
	      hDeltaPhivsPt075           ->Fill(deltaphiTPC,tPtCorr);	       
	      hCosDeltaPhiVZEROvsPt075   ->Fill(TMath::Cos(deltaphiV0),tPtCorr);    
	      hSinDeltaPhiVZEROvsPt075   ->Fill(TMath::Sin(deltaphiV0),tPtCorr);    
	      hDeltaPhiVZEROvsPt075      ->Fill(deltaphiV0,tPtCorr);                
	      
	    }
	    
	    if(tCentrality>15 && tCentrality<30){
	      hCosDeltaPhivsPt1530        ->Fill(TMath::Cos(deltaphiTPC),tPtCorr);   
	      hSinDeltaPhivsPt1530        ->Fill(TMath::Sin(deltaphiTPC),tPtCorr);   
	      hDeltaPhivsPt1530           ->Fill(deltaphiTPC,tPtCorr);	       
	      hCosDeltaPhiVZEROvsPt1530   ->Fill(TMath::Cos(deltaphiV0),tPtCorr);    
	      hSinDeltaPhiVZEROvsPt1530   ->Fill(TMath::Sin(deltaphiV0),tPtCorr);    
	      hDeltaPhiVZEROvsPt1530      ->Fill(deltaphiV0,tPtCorr);          
	    }
	
	    if(tCentrality>30 && tCentrality<50){
	      hCosDeltaPhivsPt3050        ->Fill(TMath::Cos(deltaphiTPC),tPtCorr);   
	      hSinDeltaPhivsPt3050        ->Fill(TMath::Sin(deltaphiTPC),tPtCorr);   
	      hDeltaPhivsPt3050           ->Fill(deltaphiTPC,tPtCorr);	       
	      hCosDeltaPhiVZEROvsPt3050   ->Fill(TMath::Cos(deltaphiV0),tPtCorr);    
	      hSinDeltaPhiVZEROvsPt3050   ->Fill(TMath::Sin(deltaphiV0),tPtCorr);    
	      hDeltaPhiVZEROvsPt3050      ->Fill(deltaphiV0,tPtCorr);      
	    
	    }
	  
	    if(tCentrality>15 && tCentrality<50){
	      hCosDeltaPhivsPt1550        ->Fill(TMath::Cos(deltaphiTPC),tPtCorr);   
	      hSinDeltaPhivsPt1550        ->Fill(TMath::Sin(deltaphiTPC),tPtCorr);   
	      hDeltaPhivsPt1550           ->Fill(deltaphiTPC,tPtCorr);	       
	      hCosDeltaPhiVZEROvsPt1550   ->Fill(TMath::Cos(deltaphiV0),tPtCorr);    
	      hSinDeltaPhiVZEROvsPt1550   ->Fill(TMath::Sin(deltaphiV0),tPtCorr);    
	      hDeltaPhiVZEROvsPt1550      ->Fill(deltaphiV0,tPtCorr);       
	    }
	  }
	
	  if(tCharge < 0){
	  
	    hAntiCosPhivsPt           ->Fill(TMath::Cos(tPhi),tPtCorr);
	    hAntiSinPhivsPt           ->Fill(TMath::Sin(tPhi),tPtCorr);
	    hAntiPhivsPt              ->Fill(tPhi,tPtCorr);

	    if(tCentrality>0 && tCentrality<7.5){
	      hAntiCosDeltaPhivsPt075        ->Fill(TMath::Cos(deltaphiTPC),tPtCorr);   
	      hAntiSinDeltaPhivsPt075        ->Fill(TMath::Sin(deltaphiTPC),tPtCorr);   
	      hAntiDeltaPhivsPt075           ->Fill(deltaphiTPC,tPtCorr);	       
	      hAntiCosDeltaPhiVZEROvsPt075   ->Fill(TMath::Cos(deltaphiV0),tPtCorr);    
	      hAntiSinDeltaPhiVZEROvsPt075   ->Fill(TMath::Sin(deltaphiV0),tPtCorr);    
	      hAntiDeltaPhiVZEROvsPt075      ->Fill(deltaphiV0,tPtCorr);                

	  
	    }

	    if(tCentrality>15 && tCentrality<30){
	      hAntiCosDeltaPhivsPt1530        ->Fill(TMath::Cos(deltaphiTPC),tPtCorr);   
	      hAntiSinDeltaPhivsPt1530        ->Fill(TMath::Sin(deltaphiTPC),tPtCorr);   
	      hAntiDeltaPhivsPt1530           ->Fill(deltaphiTPC,tPtCorr);	       
	      hAntiCosDeltaPhiVZEROvsPt1530   ->Fill(TMath::Cos(deltaphiV0),tPtCorr);    
	      hAntiSinDeltaPhiVZEROvsPt1530   ->Fill(TMath::Sin(deltaphiV0),tPtCorr);    
	      hAntiDeltaPhiVZEROvsPt1530      ->Fill(deltaphiV0,tPtCorr);  
	    
	      
	    }
	
	    if(tCentrality>30 && tCentrality<50){
	      hAntiCosDeltaPhivsPt3050        ->Fill(TMath::Cos(deltaphiTPC),tPtCorr);   
	      hAntiSinDeltaPhivsPt3050        ->Fill(TMath::Sin(deltaphiTPC),tPtCorr);   
	      hAntiDeltaPhivsPt3050           ->Fill(deltaphiTPC,tPtCorr);	       
	      hAntiCosDeltaPhiVZEROvsPt3050   ->Fill(TMath::Cos(deltaphiV0),tPtCorr);    
	      hAntiSinDeltaPhiVZEROvsPt3050   ->Fill(TMath::Sin(deltaphiV0),tPtCorr);    
	      hAntiDeltaPhiVZEROvsPt3050      ->Fill(deltaphiV0,tPtCorr);   
	    
             
	    }
	  
	    if(tCentrality>15 && tCentrality<50){
	      hAntiCosDeltaPhivsPt1550        ->Fill(TMath::Cos(deltaphiTPC),tPtCorr);   
	      hAntiSinDeltaPhivsPt1550        ->Fill(TMath::Sin(deltaphiTPC),tPtCorr);   
	      hAntiDeltaPhivsPt1550           ->Fill(deltaphiTPC,tPtCorr);	       
	      hAntiCosDeltaPhiVZEROvsPt1550   ->Fill(TMath::Cos(deltaphiV0),tPtCorr);    
	      hAntiSinDeltaPhiVZEROvsPt1550   ->Fill(TMath::Sin(deltaphiV0),tPtCorr);    
	      hAntiDeltaPhiVZEROvsPt1550      ->Fill(deltaphiV0,tPtCorr);    
	    
            
	    }

	  }
	
	}


	if(tTOFtrack[1] > 0){
	  Double_t beta = tTOFtrack[2]/(tTOFtrack[1] * 2.99792457999999984e-02);
	  Float_t gamma = 1/TMath::Sqrt(1 - beta*beta);
	  fhMassTOF->Fill(poutTPC/TMath::Sqrt(gamma*gamma - 1) - 1.8756);
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
	  Double_t beta = tTOFtrack[2]/(tTOFtrack[1] * 2.99792457999999984e-02);
	  Float_t gamma = 1/TMath::Sqrt(1 - beta*beta);
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
	
	tEventNumberMC[0] = evNumber    ;
	tEventNumberMC[1] = runNumber   ;
	tEventNumberMC[2] = BCNumber    ;
	tEventNumberMC[3] = OrbitNumber ;
	tEventNumberMC[4] = PeriodNumber;
	tEventNumberMC[5] = TrackNumber;
	  
	tCentralityMC = percentile;
	
	tVertexCoordMC[0] = xPrimaryVertex;
	tVertexCoordMC[1] = yPrimaryVertex;
	tVertexCoordMC[2] = zPrimaryVertex;
	
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
  PostData(1, fListHistCascade);
  PostData(2,  fNtuple1);
  PostData(3,  fNtuple2);
} //end userexec


//________________________________________________________________________

void AliAnalysisTaskNucleiv2::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
}

