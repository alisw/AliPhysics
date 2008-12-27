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

// --- ROOT system ---
#include <TH2.h>
#include <TTree.h>
// --- Standard library ---

// --- AliRoot header files ---
#include "AliITSQADataMakerRec.h"
#include "AliITSQASPDDataMakerRec.h"
#include "AliITSQASDDDataMakerRec.h"
#include "AliITSQASSDDataMakerRec.h"
#include "AliLog.h"
#include "AliQA.h"
#include "AliQAChecker.h"
#include "AliITSQAChecker.h"
#include "AliRawReader.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliMultiplicity.h"

ClassImp(AliITSQADataMakerRec)

//____________________________________________________________________________ 
AliITSQADataMakerRec::AliITSQADataMakerRec(Bool_t kMode, Short_t subDet, Short_t ldc) :
AliQADataMakerRec(AliQA::GetDetName(AliQA::kITS), "ITS Quality Assurance Data Maker"),
fkOnline(kMode),
fHLTMode(0),
fSubDetector(subDet),
fLDC(ldc),
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
    AliDebug(1,"AliITSQADM::Create SPD DataMakerRec\n");
	fSPDDataMaker = new AliITSQASPDDataMakerRec(this,fkOnline);
  }
  if(fSubDetector == 0 || fSubDetector == 2) {
	AliDebug(1,"AliITSQADM::Create SDD DataMakerRec\n");
	fSDDDataMaker = new AliITSQASDDDataMakerRec(this,fkOnline);
	if(fkOnline){SetHLTMode(fSDDDataMaker->GetHLTMode()); }
  }
  if(fSubDetector == 0 || fSubDetector == 3) {
	AliDebug(1,"AliITSQADM::Create SSD DataMakerRec\n");
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
fHLTMode(qadm.fHLTMode),
fSubDetector(qadm.fSubDetector),
fLDC(qadm.fLDC),
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
  AliDebug(1,"AliITSQADM::Start of ITS Cycle\n");
  if(fSubDetector == 0 || fSubDetector == 1) fSPDDataMaker->StartOfDetectorCycle();
  if(fSubDetector == 0 || fSubDetector == 2) fSDDDataMaker->StartOfDetectorCycle();
  if(fSubDetector == 0 || fSubDetector == 3) fSSDDataMaker->StartOfDetectorCycle();
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray** list)
{
  // launch the QA checking

  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    SetEventSpecie(specie) ; 
    AliDebug(1,"AliITSDM instantiates checker with Run(AliQA::kITS, task, list[specie])\n"); 
    if(fSubDetector == 0 || fSubDetector == 1) fSPDDataMaker->EndOfDetectorCycle(task, list[specie]);
    if(fSubDetector == 0 || fSubDetector == 2) fSDDDataMaker->EndOfDetectorCycle(task, list[specie]);
    if(fSubDetector == 0 || fSubDetector == 3) fSSDDataMaker->EndOfDetectorCycle(task, list[specie]);
  }
  
  AliQAChecker *qac = AliQAChecker::Instance();
  AliITSQAChecker *qacb = (AliITSQAChecker *) qac->GetDetQAChecker(0);
  Int_t subdet=GetSubDet();
  qacb->SetSubDet(subdet);
 
  if(subdet== 0 ){
    qacb->SetTaskOffset(fSPDDataMaker->GetOffset(task), fSDDDataMaker->GetOffset(task), fSSDDataMaker->GetOffset(task)); //Setting the offset for the QAChecker list
  }
 else
    if(subdet!=0){
      Int_t offset=GetDetTaskOffset(subdet, task);
      qacb->SetDetTaskOffset(subdet,offset);
    }

  qac->Run( AliQA::kITS , task, list);  //temporary skipping the checking
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::EndOfDetectorCycle(const char * /*fgDataName*/)
{
  //eventually used for different  AliQAChecker::Instance()->Run
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::InitRaws()
{  
  // Initialization for RAW data 
	if(fSubDetector == 0 || fSubDetector == 1) {
	  AliDebug(1,"AliITSQADM:: SPD InitRaws\n");
	  fSPDDataMaker->InitRaws();
	}
	if(fSubDetector == 0 || fSubDetector == 2) {
 	  AliDebug(1,"AliITSQADM:: SDD InitRaws\n");
	  fSDDDataMaker->InitRaws();
	}
	if(fSubDetector == 0 || fSubDetector == 3) {
	  AliDebug(1,"AliITSQADM:: SSD InitRaws\n");
	  fSSDDataMaker->InitRaws();
	}
}

//____________________________________________________________________________
void AliITSQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{ 
  // Fill QA for RAW   
  if(fSubDetector == 0 || fSubDetector == 1) fSPDDataMaker->MakeRaws(rawReader);
  if(fSubDetector == 0 || fSubDetector == 2) fSDDDataMaker->MakeRaws(rawReader);
  if(fSubDetector == 0 || fSubDetector == 3) fSSDDataMaker->MakeRaws(rawReader);
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::InitRecPoints()
{
  // Initialization for RECPOINTS
  if(fSubDetector == 0 || fSubDetector == 1) {
	AliDebug(1,"AliITSQADM:: SPD InitRecPoints\n");
    fSPDDataMaker->InitRecPoints();
  }
  if(fSubDetector == 0 || fSubDetector == 2) {
	AliDebug(1,"AliITSQADM:: SDD InitRecPoints\n");
	fSDDDataMaker->InitRecPoints();
  }
  if(fSubDetector == 0 || fSubDetector == 3) {
	AliDebug(1,"AliITSQADM:: SSD InitRecPoints\n");
	fSSDDataMaker->InitRecPoints();
  }
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::MakeRecPoints(TTree * clustersTree)
{
  // Fill QA for recpoints
  if(fSubDetector == 0 || fSubDetector == 1) fSPDDataMaker->MakeRecPoints(clustersTree);
  if(fSubDetector == 0 || fSubDetector == 2) fSDDDataMaker->MakeRecPoints(clustersTree);
  if(fSubDetector == 0 || fSubDetector == 3) fSSDDataMaker->MakeRecPoints(clustersTree);
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
  hSPDTrackletsvsFiredChips0->GetXaxis()->SetTitle("N SPD Tracklets"); 
  hSPDTrackletsvsFiredChips0->GetYaxis()->SetTitle("N FiredChips Layer0"); 
  hSPDTrackletsvsFiredChips0->Sumw2();
  Add2ESDsList(hSPDTrackletsvsFiredChips0, 16, expertHistogram ); 

  TH2F* hSPDTrackletsvsFiredChips1 = 
    new TH2F("hSPDTrackletsvsFiredChips1","N SPD Tracklets vs N FiredChips Layer1",
              300,0.,300.,300,0.,300.);
  hSPDTrackletsvsFiredChips1->GetXaxis()->SetTitle("N SPD Tracklets"); 
  hSPDTrackletsvsFiredChips1->GetYaxis()->SetTitle("N FiredChips Layer1"); 
  hSPDTrackletsvsFiredChips1->Sumw2();
  Add2ESDsList(hSPDTrackletsvsFiredChips1, 17, expertHistogram); 
    
  TH1F* hSPDTrackletsDePhi = 
    new TH1F("hSPDTrackletsDePhi","DeltaPhi SPD Tracklets; DeltaPhi [rad]; N events",200,-0.2,0.2);
  hSPDTrackletsDePhi->Sumw2();
  Add2ESDsList(hSPDTrackletsDePhi, 18); 
    
  TH1F* hSPDTrackletsPhi = 
    new TH1F("hSPDTrackletsPhi","Phi SPD Tracklets; Phi [rad]; N events",1000,0.,2*TMath::Pi());
  hSPDTrackletsPhi->Sumw2();
  Add2ESDsList(hSPDTrackletsPhi, 19); 
    
  TH1F* hSPDTrackletsTheta = 
    new TH1F("hSPDTrackletsTheta","Theta SPD Tracklets; Theta [rad]; N events",500,0.,TMath::Pi());
  hSPDTrackletsTheta->Sumw2();
  Add2ESDsList(hSPDTrackletsTheta, 20); 

  // map of layers skipped by tracking (set in AliITSRecoParam)
  TH1F *hESDSkippedLayers = 
    new TH1F("hESDSkippedLayers", "Map of layers skipped by tracking; Layer; Skipped",
	     6, -0.5, 5.5);
  hESDSkippedLayers->Sumw2();
  hESDSkippedLayers->SetMinimum(0);
  Add2ESDsList(hESDSkippedLayers, 21, expertHistogram);


  return;
}

//____________________________________________________________________________
void AliITSQADataMakerRec::MakeESDs(AliESDEvent *esd)
{
  // Make QA data from ESDs
  
  const Int_t nESDTracks = esd->GetNumberOfTracks();
  Int_t nITSrefit5 = 0; 

  Int_t idet,status;
  Float_t xloc,zloc;

  // loop on tracks
  for(Int_t i = 0; i < nESDTracks; i++) {
    
    AliESDtrack *track = esd->GetTrack(i);
    
    Int_t nclsITS = track->GetNcls(0);

    Bool_t itsrefit=kFALSE,tpcin=kFALSE;
    if ((track->GetStatus() & AliESDtrack::kITSrefit)) itsrefit=kTRUE;
    if ((track->GetStatus() & AliESDtrack::kTPCin)) tpcin=kTRUE;     
    if(nclsITS>=5 && itsrefit) nITSrefit5++;

    if(tpcin) {
      GetESDsData(0)->Fill(nclsITS);
    } else {
      GetESDsData(2)->Fill(nclsITS);
    }

    for(Int_t layer=0; layer<6; layer++) {

      if(TESTBIT(track->GetITSClusterMap(),layer)) {
	if(tpcin) {
	  GetESDsData(1)->Fill(layer);
	} else {
	  GetESDsData(3)->Fill(layer);
	}
      }
      track->GetITSModuleIndexInfo(layer,idet,status,xloc,zloc);
      if(status==3) GetESDsData(21)->SetBinContent(layer,1);
    }     

  } // end loop on tracks

  // vertices
  const AliESDVertex *vtxSPD = esd->GetPrimaryVertexSPD();
  const AliESDVertex *vtxTrk = esd->GetPrimaryVertex();

  Int_t mult = ((AliMultiplicity*)(esd->GetMultiplicity()))->GetNumberOfTracklets();
  if(mult>0)
    GetESDsData(7)->Fill((Float_t)(vtxSPD->GetNContributors())/(Float_t)mult);

  if(nITSrefit5>0)
    GetESDsData(11)->Fill((Float_t)(vtxTrk->GetNIndices())/(Float_t)nITSrefit5);

  if(vtxSPD->GetNContributors()>0) {
    GetESDsData(4)->Fill(vtxSPD->GetXv());
    GetESDsData(5)->Fill(vtxSPD->GetYv());
    GetESDsData(6)->Fill(vtxSPD->GetZv());
  }

  if(vtxTrk->GetNContributors()>0) {
    GetESDsData(8)->Fill(vtxTrk->GetXv());
    GetESDsData(9)->Fill(vtxTrk->GetYv());
    GetESDsData(10)->Fill(vtxTrk->GetZv());
  }

  if(vtxSPD->GetNContributors()>0 && 
     vtxTrk->GetNContributors()>0) {
    GetESDsData(12)->Fill(vtxSPD->GetXv()-vtxTrk->GetXv());
    GetESDsData(13)->Fill(vtxSPD->GetYv()-vtxTrk->GetYv());
    GetESDsData(14)->Fill(vtxSPD->GetZv()-vtxTrk->GetZv());
  }

  // SPD Tracklets
  GetESDsData(15)->Fill(mult);

  Short_t nFiredChips0 = ((AliMultiplicity*)(esd->GetMultiplicity()))->GetNumberOfFiredChips(0);
  Short_t nFiredChips1 = ((AliMultiplicity*)(esd->GetMultiplicity()))->GetNumberOfFiredChips(1);
  GetESDsData(16)->Fill(nFiredChips0,mult);
  GetESDsData(17)->Fill(nFiredChips1,mult);

  // Loop over tracklets
  for (Int_t itr=0; itr<mult; ++itr) {
    Float_t dePhiTr = ((AliMultiplicity*)(esd->GetMultiplicity()))->GetDeltaPhi(itr);
    Float_t phiTr   = ((AliMultiplicity*)(esd->GetMultiplicity()))->GetPhi(itr);
    Float_t thetaTr = ((AliMultiplicity*)(esd->GetMultiplicity()))->GetTheta(itr);
    GetESDsData(18)->Fill(dePhiTr);
    GetESDsData(19)->Fill(phiTr);
    GetESDsData(20)->Fill(thetaTr);
  } // end loop on tracklets

  return;
}

//_________________________________________________________________
Int_t AliITSQADataMakerRec::GetDetTaskOffset(Int_t subdet,AliQA::TASKINDEX_t task)
{
  switch(subdet)
    {

      Int_t offset;
    case 1:
      offset=fSPDDataMaker->GetOffset(task);
      return offset;
      break;
    case 2:
      offset=fSDDDataMaker->GetOffset(task);
      return offset;
      break;
    case 3:
      offset=fSSDDataMaker->GetOffset(task);
      return offset;
      break;
    default:
      AliWarning("No specific subdetector (SPD, SDD, SSD) selected!! Offset set to zero \n");
      offset=0;
      return offset;
      break;
    }
  //return offset;
}
