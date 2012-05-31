/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

//  *************************************************************
//  Checks the quality assurance 
//  by comparing with reference data
//  contained in a DB
//  -------------------------------------------------------------
//  W. Ferrarese + P. Cerello Feb 2008
//  INFN Torino
//  Melinda Siciliano Aug 2008


// --- ROOT system ---
#include <TH2.h>
// --- Standard library ---

// --- AliRoot header files ---
#include "AliITSQADataMakerRec.h"
#include "AliITSQASPDDataMakerRec.h"
#include "AliITSQASDDDataMakerRec.h"
#include "AliITSQASSDDataMakerRec.h"
#include "AliQAv1.h"
#include "AliQAChecker.h"
#include "AliITSQAChecker.h"
#include "AliITSRecPoint.h"
#include "AliITSRecPointContainer.h"
#include "AliRawReader.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliMultiplicity.h"
#include "AliITSgeomTGeo.h"

//class TH2;
//class TH2F;
class AliESDVertex;
class AliLog;
class TTree;

ClassImp(AliITSQADataMakerRec)

//____________________________________________________________________________ 
AliITSQADataMakerRec::AliITSQADataMakerRec(Bool_t kMode, Short_t subDet, Short_t ldc) :
AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kITS), "ITS Quality Assurance Data Maker"),
  fkOnline(kMode),
  fSubDetector(subDet),
  fLDC(ldc),
  fRunNumber(0),
  fEventNumber(0),
  fSelectedTaskIndex(AliQAv1::kNULLTASKINDEX),
  fSPDDataMaker(NULL),
  fSDDDataMaker(NULL),
  fSSDDataMaker(NULL)

{
  //ctor used to discriminate OnLine-Offline analysis
  if(fSubDetector < 0 || fSubDetector > 3) {
    AliError("Error: fSubDetector number out of range; return\n");
  }

  // Initialization for RAW data 
  if(fSubDetector == 0 || fSubDetector == 1) {
    AliDebug(AliQAv1::GetQADebugLevel(),"AliITSQADM::Create SPD DataMakerRec\n");
    fSPDDataMaker = new AliITSQASPDDataMakerRec(this,fkOnline);
  }
  if(fSubDetector == 0 || fSubDetector == 2) {
    AliDebug(AliQAv1::GetQADebugLevel(),"AliITSQADM::Create SDD DataMakerRec\n");
    fSDDDataMaker = new AliITSQASDDDataMakerRec(this,fkOnline);
  }
  if(fSubDetector == 0 || fSubDetector == 3) {
    AliDebug(AliQAv1::GetQADebugLevel(),"AliITSQADM::Create SSD DataMakerRec\n");
    fSSDDataMaker = new AliITSQASSDDataMakerRec(this,fkOnline);
  }
}

//____________________________________________________________________________ 
AliITSQADataMakerRec::~AliITSQADataMakerRec(){
  // destructor
  if(fSPDDataMaker)delete fSPDDataMaker;
  if(fSDDDataMaker)delete fSDDDataMaker;
  if(fSSDDataMaker)delete fSSDDataMaker;
}

//____________________________________________________________________________ 
AliITSQADataMakerRec::AliITSQADataMakerRec(const AliITSQADataMakerRec& qadm) :
  AliQADataMakerRec(),
  fkOnline(qadm.fkOnline),
  fSubDetector(qadm.fSubDetector),
  fLDC(qadm.fLDC),
  fRunNumber(qadm.fRunNumber),
  fEventNumber(qadm.fEventNumber),
  fSelectedTaskIndex(qadm.fSelectedTaskIndex),
  fSPDDataMaker(NULL),
  fSDDDataMaker(NULL),
  fSSDDataMaker(NULL)

{
  //copy ctor 
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle());
}

//__________________________________________________________________
AliITSQADataMakerRec& AliITSQADataMakerRec::operator = (const AliITSQADataMakerRec& qac )
{
  // Equal operator.
  this->~AliITSQADataMakerRec();
  new(this) AliITSQADataMakerRec(qac);
  return *this;
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  AliDebug(AliQAv1::GetQADebugLevel(),"AliITSQADM::Start of ITS Cycle\n");
  ResetEventTrigClasses(); // reset triggers list to select all histos
  ResetEvCountCycle();
  //  
  if(fSubDetector == 0 || fSubDetector == 1) fSPDDataMaker->StartOfDetectorCycle();
  if(fSubDetector == 0 || fSubDetector == 2) fSDDDataMaker->StartOfDetectorCycle();
  if(fSubDetector == 0 || fSubDetector == 3) fSSDDataMaker->StartOfDetectorCycle();
}

//____________________________________________________________________________
void AliITSQADataMakerRec::StartOfCycle(AliQAv1::TASKINDEX_t task, Int_t run, const Bool_t sameCycle) 
{ 
  // Start a cycle of QA data acquistion
  fSelectedTaskIndex=task;
  AliQADataMakerRec::StartOfCycle(task,run,sameCycle);
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray** list)
{
  // launch the QA checking
  ResetEventTrigClasses();
  //
  AliInfo(Form("End of Dedetctor Cycle called for %s\n",AliQAv1::GetTaskName(task).Data() ));
  //
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    //
    if(!AliQAv1::Instance()->IsEventSpecieSet(specie)) continue;
    SetEventSpecie(AliRecoParam::ConvertIndex(specie));
    Int_t idnumber=list[specie]->GetUniqueID();
    //printf("specie %s \t id number == %d\n",AliRecoParam::GetEventSpecieName(specie),idnumber);
    if(idnumber==40||idnumber==0){
      //AliInfo(Form("No check for %s\n",AliQAv1::GetTaskName(task).Data() ))
      continue;
    } //skip kDigitsR and not filled TobjArray specie
    else{
      AliDebug(AliQAv1::GetQADebugLevel(),"AliITSDM instantiates checker with Run(AliQAv1::kITS, task, list[specie])\n"); 
      if(fSubDetector == 0 || fSubDetector == 1) fSPDDataMaker->EndOfDetectorCycle(task, list);//[/*GetEventSpecie()*/specie]);
      if(fSubDetector == 0 || fSubDetector == 2) fSDDDataMaker->EndOfDetectorCycle(task, list);//[/*GetEventSpecie()*/specie]);
      if(fSubDetector == 0 || fSubDetector == 3) fSSDDataMaker->EndOfDetectorCycle(task, list);//[/*GetEventSpecie()*/specie]);
      
      
      AliQAChecker *qac = AliQAChecker::Instance();
      AliITSQAChecker *qacb = (AliITSQAChecker *) qac->GetDetQAChecker(0);
      Int_t subdet=GetSubDet();
      qacb->SetSubDet(subdet);
      
      if(subdet== 0 ){
	qacb->SetTaskOffset(fSPDDataMaker->GetOffset(task,specie), fSDDDataMaker->GetOffset(task,specie), fSSDDataMaker->GetOffset(task,specie)); //Setting the offset for the QAChecker list
	qacb->SetHisto(fSPDDataMaker->GetTaskHisto(task), fSDDDataMaker->GetTaskHisto(task), fSSDDataMaker->GetTaskHisto(task));
      }
      else
	if(subdet!=0){
	  Int_t offset=GetDetTaskOffset(subdet, task,specie);
	  qacb->SetDetTaskOffset(subdet,offset);
	  Int_t histo=GetDetTaskHisto(subdet, task);
	  qacb->SetDetHisto(subdet,histo);
	}
      
      qac->Run( AliQAv1::kITS , task, list);
      
    }//end else unique id
    
  }//end for
}

//____________________________________________________________________________ 
//void AliITSQADataMakerRec::EndOfDetectorCycle(const char * /*fgDataName*/)
//{
//eventually used for different  AliQAChecker::Instance()->Run
//}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::InitRaws() {
  // Initialization of RAW data histograms  

  //if(fRawsQAList[AliRecoParam::AConvert(fEventSpecie)]->GetEntries()) return;
	
  if(fSubDetector == 0 || fSubDetector == 1) {
    AliDebug(AliQAv1::GetQADebugLevel(),"AliITSQADM:: SPD InitRaws\n");
    fSPDDataMaker->InitRaws();
  }
  if(fSubDetector == 0 || fSubDetector == 2) {
    AliDebug(AliQAv1::GetQADebugLevel(),"AliITSQADM:: SDD InitRaws\n");
    
    fSDDDataMaker->SetOffset(AliQAv1::kRAWS, fRawsQAList[AliRecoParam::AConvert(fEventSpecie)]->GetEntries(),AliRecoParam::AConvert(fEventSpecie));
    fSDDDataMaker->InitRaws();
  }
  if(fSubDetector == 0 || fSubDetector == 3) {
    AliDebug(AliQAv1::GetQADebugLevel(),"AliITSQADM:: SSD InitRaws\n");
    
    fSSDDataMaker->SetOffset(AliQAv1::kRAWS, fRawsQAList[AliRecoParam::AConvert(fEventSpecie)]->GetEntries(),AliRecoParam::AConvert(fEventSpecie));
    fSSDDataMaker->InitRaws();
  }
  fRawsQAList[AliRecoParam::AConvert(fEventSpecie)]->SetUniqueID(10);
  //
  ClonePerTrigClass(AliQAv1::kRAWS); // this should be the last line
}

//____________________________________________________________________________
void AliITSQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{ 
  // Fill QA for RAW   
  //return ; 

  SetRunNumber(rawReader->GetRunNumber());

  if(fSubDetector == 0 || fSubDetector == 1)  {
    fSPDDataMaker->MakeRaws(rawReader) ; 
  }
  
  if(fSubDetector == 0 || fSubDetector == 2) {
    fSDDDataMaker->MakeRaws(rawReader) ; 
  }

  if(fSubDetector == 0 || fSubDetector == 3) fSSDDataMaker->MakeRaws(rawReader);
  //
  IncEvCountCycleRaws();
  IncEvCountTotalRaws();
  //
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::InitDigits()
{

  // Initialization for DIGITS
  if(fSubDetector == 0 || fSubDetector == 1) {
    AliDebug(AliQAv1::GetQADebugLevel(),"AliITSQADM:: SPD InitDigits\n");

    fSPDDataMaker->InitDigits();
  }
  if(fSubDetector == 0 || fSubDetector == 2) {
    AliDebug(AliQAv1::GetQADebugLevel(),"AliITSQADM:: SDD InitDigits\n");
    fSDDDataMaker->SetOffset(AliQAv1::kDIGITSR, fDigitsQAList[AliRecoParam::AConvert(fEventSpecie)]->GetEntries(),AliRecoParam::AConvert(fEventSpecie));

    fSDDDataMaker->InitDigits();
  }
  if(fSubDetector == 0 || fSubDetector == 3) {
    AliDebug(AliQAv1::GetQADebugLevel(),"AliITSQADM:: SSD InitDigits\n");
    fSSDDataMaker->SetOffset(AliQAv1::kDIGITSR, fDigitsQAList[AliRecoParam::AConvert(fEventSpecie)]->GetEntries(),AliRecoParam::AConvert(fEventSpecie));

    fSSDDataMaker->InitDigits();
  }
  fDigitsQAList[AliRecoParam::AConvert(fEventSpecie)]->SetUniqueID(40);
  //
  ClonePerTrigClass(AliQAv1::kDIGITS); // this should be the last line
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::MakeDigits(TTree * digitsTree)
{

  
  // Fill QA for recpoints
  if(fSubDetector == 0 || fSubDetector == 1) {
    fSPDDataMaker->MakeDigits(digitsTree) ; 
  }
  
  if(fSubDetector == 0 || fSubDetector == 2) {
    fSDDDataMaker->MakeDigits(digitsTree) ; 
    
  }
  
  if(fSubDetector == 0 || fSubDetector == 3) fSSDDataMaker->MakeDigits(digitsTree);
  //
  IncEvCountCycleDigits();
  IncEvCountTotalDigits();
  //
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::InitRecPoints()
{

  // Initialization for RECPOINTS


  //if(fRecPointsQAList[AliRecoParam::AConvert(fEventSpecie)]->GetEntries()) return;
  if(fSubDetector == 0 || fSubDetector == 1) {
    AliDebug(AliQAv1::GetQADebugLevel(),"AliITSQADM:: SPD InitRecPoints\n");
    fSPDDataMaker->InitRecPoints();
  }
  if(fSubDetector == 0 || fSubDetector == 2) {
    AliDebug(AliQAv1::GetQADebugLevel(),"AliITSQADM:: SDD InitRecPoints\n");
    fSDDDataMaker->SetOffset(AliQAv1::kRECPOINTS, fRecPointsQAList[AliRecoParam::AConvert(fEventSpecie)]->GetEntries(), AliRecoParam::AConvert(fEventSpecie));
    fSDDDataMaker->InitRecPoints();
  }
  if(fSubDetector == 0 || fSubDetector == 3) {
    AliDebug(AliQAv1::GetQADebugLevel(),"AliITSQADM:: SSD InitRecPoints\n");
    fSSDDataMaker->SetOffset(AliQAv1::kRECPOINTS, fRecPointsQAList[AliRecoParam::AConvert(fEventSpecie)]->GetEntries(),AliRecoParam::AConvert(fEventSpecie));
    fSSDDataMaker->InitRecPoints();
  }

  fRecPointsQAList[AliRecoParam::AConvert(fEventSpecie)]->SetUniqueID(20);
  if(fSubDetector == 0){
    Int_t offset = fRecPointsQAList [AliRecoParam::AConvert(fEventSpecie)]->GetEntries();
    const Bool_t expert   = kTRUE ; 
    const Bool_t image    = kTRUE ; 
    TH2F* hPhiEta[6];
    for (Int_t iLay=0;iLay<6;iLay++) {
      hPhiEta[iLay]=new TH2F(Form("Phi_vs_Eta_ITS_Layer%d",iLay+1),Form("Phi_vs_Eta_ITS_Layer%d",iLay+1),30,-1.5,1.5,200,0.,2*TMath::Pi());
      hPhiEta[iLay]->GetXaxis()->SetTitle("Pseudorapidity");
      hPhiEta[iLay]->GetYaxis()->SetTitle("#varphi [rad]");
      Add2RecPointsList(hPhiEta[iLay], iLay + offset, !expert, image);	    
      //delete hPhiEta[iLay];
    }
	  
  }
  //
  ClonePerTrigClass(AliQAv1::kRECPOINTS); // this should be the last line	
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::MakeRecPoints(TTree * clustersTree)
{ 
  // Fill QA for recpoints

  if(fSubDetector == 0 || fSubDetector == 1) {
    fSPDDataMaker->MakeRecPoints(clustersTree) ; 
  }
    
  if(fSubDetector == 0 || fSubDetector == 2) {
    fSDDDataMaker->MakeRecPoints(clustersTree) ; 
  }
  
  if(fSubDetector == 0 || fSubDetector == 3) fSSDDataMaker->MakeRecPoints(clustersTree);


  
  if(fSubDetector == 0){

    // Check id histograms already created for this Event Specie
    AliITSRecPointContainer* rpcont=AliITSRecPointContainer::Instance();
    TClonesArray *recpoints =NULL;
    if(fkOnline){
      rpcont->FetchClusters(0,clustersTree,GetEventNumber());
    } 
    else{
      rpcont->FetchClusters(0,clustersTree);
    }
    if(!rpcont->GetStatusOK()){
      AliError("cannot access to ITS recpoints");
      return;
    }
  
    Int_t offset = fRecPointsQAList [AliRecoParam::AConvert(fEventSpecie)]->GetEntries();
    Float_t cluGlo[3] = {0.,0.,0.};
    Int_t lay, lad, det; 
    // Fill QA for recpoints
    for(Int_t module=0; module<rpcont->GetNumberOfModules();module++){
      //  AliInfo(Form("Module %d\n",module));
      recpoints = rpcont->UncheckedGetClusters(module);
      AliITSgeomTGeo::GetModuleId(module, lay, lad, det);
      for(Int_t j=0;j<recpoints->GetEntries();j++){
	AliITSRecPoint *rcp = (AliITSRecPoint*)recpoints->At(j);    
	//Check id histograms already created for this Event Specie
	rcp->GetGlobalXYZ(cluGlo);
	Double_t rad=TMath::Sqrt(cluGlo[0]*cluGlo[0]+cluGlo[1]*cluGlo[1]+cluGlo[2]*cluGlo[2]);
	Double_t phi= TMath::Pi() + TMath::ATan2(-cluGlo[1],-cluGlo[0]);
	Double_t theta = TMath::ACos(cluGlo[2]/rad);
	Double_t eta = 100.;
	if(AreEqual(rad,0.) == kFALSE) {
	  if(theta<=1.e-14){ eta=30.; }
	  else { eta = -TMath::Log(TMath::Tan(theta/2.));}
	}
	//	printf("=========================>hlt   rcp->GetLayer() = %d \n",rcp->GetLayer());
	FillRecPointsData(rcp->GetLayer() + offset - 6,eta,phi);
      }
    }
  }
  //
  IncEvCountCycleRecPoints();
  IncEvCountTotalRecPoints();
  //
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::FillRecPoint(AliITSRecPoint rcp)
{

  // Fill QA for recpoints
  Float_t cluGlo[3] = {0.,0.,0.};
  Int_t offset = fRecPointsQAList [AliRecoParam::AConvert(fEventSpecie)]->GetEntries();
  // Check id histograms already created for this Event Specie
  rcp.GetGlobalXYZ(cluGlo);
  Double_t rad=TMath::Sqrt(cluGlo[0]*cluGlo[0]+cluGlo[1]*cluGlo[1]+cluGlo[2]*cluGlo[2]);
  Double_t phi= TMath::Pi() + TMath::ATan2(-cluGlo[1],-cluGlo[0]);
  Double_t theta = TMath::ACos(cluGlo[2]/rad);
  Double_t eta = 100.;
  if(AreEqual(rad,0.)==kFALSE) {
    if(theta<=1.e-14){eta=30.;}
    else    {eta = -TMath::Log(TMath::Tan(theta/2.));}
  }
  FillRecPointsData( rcp.GetLayer() + offset - 6, eta,phi);	

}

//____________________________________________________________________________ 
TH2F *AliITSQADataMakerRec::GetITSGlobalHisto(Int_t layer)
{

  Int_t offset = fRecPointsQAList [AliRecoParam::AConvert(fEventSpecie)]->GetEntries();
  return ((TH2F *) GetRecPointsData( layer + offset - 6));//local distribution
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::InitESDs()
{

  // Create ESDs histograms in ESDs subdir

  Bool_t expertHistogram = kTRUE;

  TH1F *hESDClustersMI = 
    new TH1F("hESDClustersMI", "N ITS clusters per track (MI); N clusters; Counts",
	     7, -0.5, 6.5);
  hESDClustersMI->Sumw2();
  hESDClustersMI->SetMinimum(0);
  Add2ESDsList(hESDClustersMI, 0);

  TH1F *hESDClusterMapMI =
    new TH1F("hESDClusterMapMI", "N tracks with point Layer (MI); Layer; N tracks",
	     6, -0.5, 5.5);
  hESDClusterMapMI->Sumw2();
  hESDClusterMapMI->SetMinimum(0);
  Add2ESDsList(hESDClusterMapMI, 1, expertHistogram);

  TH1F *hESDClustersSA = 
    new TH1F("hESDClustersSA", "N ITS clusters per track (SA); N clusters; Counts",
	     7, -0.5, 6.5);
  hESDClustersSA->Sumw2();
  hESDClustersSA->SetMinimum(0);
  Add2ESDsList(hESDClustersSA, 2);

  TH1F *hESDClusterMapSA =
    new TH1F("hESDClusterMapSA", "N tracks with point Layer (SA); Layer; N tracks",
	     6, -0.5, 5.5);
  hESDClusterMapSA->Sumw2();
  hESDClusterMapSA->SetMinimum(0);
  Add2ESDsList(hESDClusterMapSA, 3, expertHistogram);

  TH1F *hSPDVertexX = 
    new TH1F("hSPDVertexX","SPD Vertex x; x [cm]; N events",
	     10000,-2,2);
  hSPDVertexX->Sumw2();
  Add2ESDsList(hSPDVertexX, 4);

  TH1F *hSPDVertexY = 
    new TH1F("hSPDVertexY","SPD Vertex y; y [cm]; N events",
	     10000,-2,2);
  hSPDVertexY->Sumw2();
  Add2ESDsList(hSPDVertexY, 5);

  TH1F *hSPDVertexZ = 
    new TH1F("hSPDVertexZ","SPD Vertex Z; z [cm]; N events",
	     10000,-20,20);
  hSPDVertexZ->Sumw2();
  Add2ESDsList(hSPDVertexZ, 6);

  TH1F *hSPDVertexContrOverMult =
    new TH1F("hSPDVertexContrOverMult","SPD Vertex: contributors / multiplicity; N contributors / SPD multiplicity; N events",
	     100,-4,20);
  hSPDVertexContrOverMult->Sumw2();
  Add2ESDsList(hSPDVertexContrOverMult, 7, expertHistogram);

  TH1F *hTrkVertexX = 
    new TH1F("hTrkVertexX","ITS+TPC Trk Vertex x; x [cm]; N events",
	     10000,-2,2);
  hTrkVertexX->Sumw2();
  Add2ESDsList(hTrkVertexX, 8, expertHistogram);

  TH1F *hTrkVertexY = 
    new TH1F("hTrkVertexY","ITS+TPC Trk Vertex y; y [cm]; N events",
	     10000,-2,2);
  hTrkVertexY->Sumw2();
  Add2ESDsList(hTrkVertexY, 9, expertHistogram);

  TH1F *hTrkVertexZ = 
    new TH1F("hTrkVertexZ","ITS+TPC Trk Vertex Z; z [cm]; N events",
	     10000,-20,20);
  hTrkVertexZ->Sumw2();
  Add2ESDsList(hTrkVertexZ, 10, expertHistogram);

  TH1F *hTrkVertexContrOverITSrefit5 =
    new TH1F("hTrkVertexContrOverITSrefit5","ITS+TPC Trk Vertex: contributors / tracks; N contributors / N trks kITSrefit with 5 or 6 clusters; N events",
	     100,-4,2);
  hTrkVertexContrOverITSrefit5->Sumw2();
  Add2ESDsList(hTrkVertexContrOverITSrefit5, 11, expertHistogram);

  TH1F *hSPDTrkVertexDeltaX =
    new TH1F("hSPDTrkVertexDeltaX","Comparison of SPD and Trk vertices: x; xSPD-xTrk [cm]; N events",
	     1000,-1,1);
  hSPDTrkVertexDeltaX->Sumw2();
  Add2ESDsList(hSPDTrkVertexDeltaX, 12, expertHistogram);
    
  TH1F *hSPDTrkVertexDeltaY =
    new TH1F("hSPDTrkVertexDeltaY","Comparison of SPD and Trk vertices: y; ySPD-yTrk [cm]; N events",
	     1000,-1,1);
  hSPDTrkVertexDeltaY->Sumw2();
  Add2ESDsList(hSPDTrkVertexDeltaY, 13, expertHistogram);
    
  TH1F *hSPDTrkVertexDeltaZ =
    new TH1F("hSPDTrkVertexDeltaZ","Comparison of SPD and Trk vertices: z; zSPD-zTrk [cm]; N events",
	     1000,-1,1);
  hSPDTrkVertexDeltaZ->Sumw2();
  Add2ESDsList(hSPDTrkVertexDeltaZ, 14);
    
  // SPD Tracklets

  TH1F* hSPDTracklets = 
    new TH1F("hSPDTracklets","N SPD Tracklets; N tracklets; Counts",300,0.,300.);
  hSPDTracklets->Sumw2();
  Add2ESDsList(hSPDTracklets, 15); 

  TH2F* hSPDTrackletsvsFiredChips0 = 
    new TH2F("hSPDTrackletsvsFiredChips0","N SPD Tracklets vs N FiredChips Layer0",
	     300,0.,300.,300,0.,300.);
  hSPDTrackletsvsFiredChips0->GetXaxis()->SetTitle("N FiredChips Layer0"); 
  hSPDTrackletsvsFiredChips0->GetYaxis()->SetTitle("N SPD Tracklets"); 
  hSPDTrackletsvsFiredChips0->Sumw2();
  Add2ESDsList(hSPDTrackletsvsFiredChips0, 16, expertHistogram ); 

  TH2F* hSPDTrackletsvsFiredChips1 = 
    new TH2F("hSPDTrackletsvsFiredChips1","N SPD Tracklets vs N FiredChips Layer1",
	     300,0.,300.,300,0.,300.);
  hSPDTrackletsvsFiredChips1->GetXaxis()->SetTitle("N FiredChips Layer1"); 
  hSPDTrackletsvsFiredChips1->GetYaxis()->SetTitle("N SPD Tracklets"); 
  hSPDTrackletsvsFiredChips1->Sumw2();
  Add2ESDsList(hSPDTrackletsvsFiredChips1, 17, expertHistogram); 

  TH2F* hSPDFiredChips1vsFiredChips0 = 
    new TH2F("hSPDFiredChips1vsFiredChips0","N FiredChips Layer1 vs N FiredChips Layer0",
	     300,0.,300.,300,0.,300.);
  hSPDFiredChips1vsFiredChips0->GetXaxis()->SetTitle("N FiredChips Layer0"); 
  hSPDFiredChips1vsFiredChips0->GetYaxis()->SetTitle("N FiredChips Layer1"); 
  hSPDFiredChips1vsFiredChips0->Sumw2();
  Add2ESDsList(hSPDFiredChips1vsFiredChips0, 18, expertHistogram ); 
    
  TH1F* hSPDTrackletsDePhi = 
    new TH1F("hSPDTrackletsDePhi","DeltaPhi SPD Tracklets; DeltaPhi [rad]; N events",200,-0.2,0.2);
  hSPDTrackletsDePhi->Sumw2();
  Add2ESDsList(hSPDTrackletsDePhi, 19); 
    
  TH1F* hSPDTrackletsPhi = 
    new TH1F("hSPDTrackletsPhi","Phi SPD Tracklets; Phi [rad]; N events",1000,0.,2*TMath::Pi());
  hSPDTrackletsPhi->Sumw2();
  Add2ESDsList(hSPDTrackletsPhi, 20); 
    
  TH1F* hSPDTrackletsDeTheta = 
    new TH1F("hSPDTrackletsDeTheta","DeltaTheta SPD Tracklets; DeltaTheta [rad]; N events",200,-0.2,0.2);
  hSPDTrackletsDeTheta->Sumw2();
  Add2ESDsList(hSPDTrackletsDeTheta, 21); 

  TH1F* hSPDTrackletsTheta = 
    new TH1F("hSPDTrackletsTheta","Theta SPD Tracklets; Theta [rad]; N events",500,0.,TMath::Pi());
  hSPDTrackletsTheta->Sumw2();
  Add2ESDsList(hSPDTrackletsTheta, 22); 

  // map of layers skipped by tracking (set in AliITSRecoParam)
  TH1F *hESDSkippedLayers = 
    new TH1F("hESDSkippedLayers", "Map of layers skipped by tracking; Layer; Skipped",
	     6, -0.5, 5.5);
  hESDSkippedLayers->Sumw2();
  hESDSkippedLayers->SetMinimum(0);
  Add2ESDsList(hESDSkippedLayers, 23, expertHistogram);

  fESDsQAList[AliRecoParam::AConvert(fEventSpecie)]->SetUniqueID(30);
  //
  ClonePerTrigClass(AliQAv1::kESDS); // this should be the last line
}

//____________________________________________________________________________
void AliITSQADataMakerRec::MakeESDs(AliESDEvent *esd)
{
  // Make QA data from ESDs

  // Check id histograms already created for this Event Specie
  //  if ( ! GetESDsData(0) )
  //    InitESDs() ;
 
  const Int_t nESDTracks = esd->GetNumberOfTracks();
  Int_t nITSrefit5 = 0; 

  Int_t idet,status;
  Float_t xloc,zloc;

  // loop on tracks
  AliInfo(Form("Filling histograms for ESD. Number of tracks %d",nESDTracks)); 
  for(Int_t i = 0; i < nESDTracks; i++) {
    
    AliESDtrack *track = esd->GetTrack(i);
    
    Int_t nclsITS = track->GetNcls(0);

    Bool_t itsrefit=kFALSE,tpcin=kFALSE,itsin=kFALSE;
    if ((track->GetStatus() & AliESDtrack::kITSrefit)) itsrefit=kTRUE;
    if ((track->GetStatus() & AliESDtrack::kTPCin)) tpcin=kTRUE;     
    if ((track->GetStatus() & AliESDtrack::kITSin)) itsin=kTRUE;     
    if(nclsITS>=5 && itsrefit) nITSrefit5++;

    if(tpcin) {
      FillESDsData(0,nclsITS);
    }
    if(itsin && !tpcin){
      FillESDsData(2,nclsITS);
    }

    for(Int_t layer=0; layer<6; layer++) {

      if(TESTBIT(track->GetITSClusterMap(),layer)) {
	if(tpcin) {
	  FillESDsData(1,layer);
	} else {
	  FillESDsData(3,layer);
	}
      }
      track->GetITSModuleIndexInfo(layer,idet,status,xloc,zloc);
      if(status==3) SetESDsDataBinContent(23,layer,1);
    }     

  } // end loop on tracks

  // vertices
  const AliESDVertex *vtxSPD = esd->GetPrimaryVertexSPD();
  const AliESDVertex *vtxTrk = esd->GetPrimaryVertexTracks();

  Int_t mult = ((AliMultiplicity*)(esd->GetMultiplicity()))->GetNumberOfTracklets();
  AliInfo(Form("Multiplicity %d ; Number of SPD vert contributors %d",mult,vtxSPD->GetNContributors()));
  if(mult>0)
    FillESDsData(7,(Float_t)(vtxSPD->GetNContributors())/(Float_t)mult);

  if(nITSrefit5>0)
    FillESDsData(11,(Float_t)(vtxTrk->GetNIndices())/(Float_t)nITSrefit5);

  if(vtxSPD->GetNContributors()>0) {
    FillESDsData(4,vtxSPD->GetXv());
    FillESDsData(5,vtxSPD->GetYv());
    FillESDsData(6,vtxSPD->GetZv());
  }

  if(vtxTrk->GetNContributors()>0) {
    FillESDsData(8,vtxTrk->GetXv());
    FillESDsData(9,vtxTrk->GetYv());
    FillESDsData(10,vtxTrk->GetZv());
  }

  if(vtxSPD->GetNContributors()>0 && 
     vtxTrk->GetNContributors()>0) {
    FillESDsData(12,vtxSPD->GetXv()-vtxTrk->GetXv());
    FillESDsData(13,vtxSPD->GetYv()-vtxTrk->GetYv());
    FillESDsData(14,vtxSPD->GetZv()-vtxTrk->GetZv());
  }

  // SPD Tracklets
  FillESDsData(15,mult);

  Short_t nFiredChips0 = ((AliMultiplicity*)(esd->GetMultiplicity()))->GetNumberOfFiredChips(0);
  Short_t nFiredChips1 = ((AliMultiplicity*)(esd->GetMultiplicity()))->GetNumberOfFiredChips(1);
  FillESDsData(16,nFiredChips0,mult);
  FillESDsData(17,nFiredChips1,mult);
  FillESDsData(18,nFiredChips0,nFiredChips1);

  // Loop over tracklets
  for (Int_t itr=0; itr<mult; ++itr) {
    Float_t dePhiTr   = ((AliMultiplicity*)(esd->GetMultiplicity()))->GetDeltaPhi(itr);
    Float_t deThetaTr = ((AliMultiplicity*)(esd->GetMultiplicity()))->GetDeltaTheta(itr);
    Float_t phiTr   = ((AliMultiplicity*)(esd->GetMultiplicity()))->GetPhi(itr);
    Float_t thetaTr = ((AliMultiplicity*)(esd->GetMultiplicity()))->GetTheta(itr);
    FillESDsData(19,dePhiTr);
    FillESDsData(20,phiTr);
    FillESDsData(21,deThetaTr);
    FillESDsData(22,thetaTr);
  } // end loop on tracklets
  //
  IncEvCountCycleESDs();
  IncEvCountTotalESDs();
  //
}

//_________________________________________________________________
Int_t AliITSQADataMakerRec::GetDetTaskOffset(Int_t subdet,AliQAv1::TASKINDEX_t task, Int_t specie)
{
  //number of booked histos for the QAchecking Raws offset
  Int_t offset=0;
  switch(subdet)
    {
    case 1:
      offset=fSPDDataMaker->GetOffset(task,specie);
      //return offset;
      break;
    case 2:
      offset=fSDDDataMaker->GetOffset(task,specie);
      //return offset;
      break;
    case 3:
      offset=fSSDDataMaker->GetOffset(task,specie);
      //return offset;
      break;
    default:
      AliWarning("No specific subdetector (SPD, SDD, SSD) selected!! Offset set to zero \n");
      offset=0;
      //return offset;
      break;
    }
  return offset;
}

//____________________________________________________________________

Bool_t AliITSQADataMakerRec::AreEqual(Double_t a1,Double_t a2)
{
  const Double_t kEpsilon= 1.e-14;
  return TMath::Abs(a1-a2)<=kEpsilon*TMath::Abs(a1);      
}

//_________________________________________________________________
Int_t AliITSQADataMakerRec::GetDetTaskHisto(Int_t subdet,AliQAv1::TASKINDEX_t task)
{
  //return the number of histo booked for each the Raws Task 

  Int_t histo=0;
  switch(subdet)
    {
    case 1:
      histo=fSPDDataMaker->GetTaskHisto(task);
      //return histo;
      break;
    case 2:
      histo=fSDDDataMaker->GetTaskHisto(task);
      //return histo;
      break;
    case 3:
      histo=fSSDDataMaker->GetTaskHisto(task);
      //return histo;
      break;
    default:
      AliWarning("No specific subdetector (SPD, SDD, SSD) selected!! Offset set to zero \n");
      histo=0;
      //return histo;
      break;
    }
  //return offset;
  return histo;
}


//____________________________________________________________________

void AliITSQADataMakerRec::ResetDetector(AliQAv1::TASKINDEX_t task)
{
  //reset the detector histograms for a given task
  AliQADataMakerRec::ResetDetector(task);

  if(fSubDetector==0||fSubDetector==1)fSPDDataMaker->ResetDetector(task);
  
  if(fSubDetector==0||fSubDetector==2)fSDDDataMaker->ResetDetector(task);

  if(fSubDetector==0||fSubDetector==3)fSSDDataMaker->ResetDetector(task);
  
}


//____________________________________________________________________

AliITSDDLModuleMapSDD *AliITSQADataMakerRec::GetDDLSDDModuleMap()
{
  //return the SDD module map
  if(fSubDetector==2){return fSDDDataMaker->GetDDLSDDModuleMap();}
  else {return NULL;}
}

//____________________________________________________________________
Bool_t AliITSQADataMakerRec::ListExists(AliQAv1::TASKINDEX_t task) const
{
  //Check the existence of a list for a given task
  Bool_t havethelist=kFALSE;
  if( ( task == AliQAv1::kRAWS && fRawsQAList ) ||
      ( task == AliQAv1::kRECPOINTS && fRecPointsQAList ) ||
      ( task == AliQAv1::kDIGITSR && fDigitsQAList ) ||
      ( task == AliQAv1::kESDS && fESDsQAList ) ) havethelist=kTRUE;
  return havethelist;

}
