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


//-----------------------------------------------------------------------
// Author : A. Mastroserio
//-----------------------------------------------------------------------


#ifndef ALIANALYSISTASKSPD_CXX
#define ALIANALYSISTASKSPD_CXX

#include "AliAnalysisTaskSPD.h"
#include "TClonesArray.h"
#include "TBranch.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TFile.h"
#include "AliSPDUtils.h"
#include "AliESDEvent.h"
#include "TChain.h"
#include "AliLog.h"
#include "AliESDInputHandlerRP.h"
#include "AliAnalysisManager.h"
#include "AliMultiplicity.h"
#include "AliCDBPath.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBStorage.h"
#include "AliGeomManager.h"
#include "AliITSRecPoint.h"
#include "AliITSsegmentationSPD.h"
ClassImp(AliAnalysisTaskSPD)
  //__________________________________________________________________________
  AliAnalysisTaskSPD::AliAnalysisTaskSPD() :
    fSegSPD(0x0),
    fOutput(0x0),
    fRunNb(999),
    fOCDBLocation("local://$ALICE_ROOT/OCDB"),
    fHI(kFALSE),
    fTest(kFALSE)
{
  //
  //Default ctor
  //
}
//___________________________________________________________________________
AliAnalysisTaskSPD::AliAnalysisTaskSPD(const Char_t* name) :
  AliAnalysisTaskSE(name),
  fSegSPD(0x0),
  fOutput(0x0),
  fRunNb(999),
  fOCDBLocation("local://$ALICE_ROOT/OCDB"),
  fHI(kFALSE),
  fTest(kFALSE)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskSPD","Calling Constructor");

  DefineOutput(1,TList::Class());
  // 
}

//___________________________________________________________________________
AliAnalysisTaskSPD& AliAnalysisTaskSPD::operator=(const AliAnalysisTaskSPD& c) 
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c) ;
    fSegSPD = c.fSegSPD ;
    fOutput = c.fOutput ;
    fRunNb = c.fRunNb;
    fOCDBLocation = c.fOCDBLocation;
    fHI = c.fHI;
    fTest = c.fTest;
  }
  return *this;
}

//___________________________________________________________________________
AliAnalysisTaskSPD::AliAnalysisTaskSPD(const AliAnalysisTaskSPD& c) :
  AliAnalysisTaskSE(c),
  fSegSPD(c.fSegSPD),
  fOutput(c.fOutput),
  fRunNb(c.fRunNb),
  fOCDBLocation(c.fOCDBLocation),
  fHI(c.fHI),
  fTest(c.fTest)
{
  //
  // Copy Constructor
  //

}

//___________________________________________________________________________
AliAnalysisTaskSPD::~AliAnalysisTaskSPD() {
  //
  //destructor
  //
 
  Info("~AliAnalysisTaskSPD","Calling Destructor");
 if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return; 
 if (fSegSPD) delete fSegSPD ;
   if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
}
//___________________________________________________________________________
void AliAnalysisTaskSPD::UserCreateOutputObjects() {

  Info("CreateOutputObjects","CreateOutputObjects of task %s", GetName());
  if(fRunNb!=999){
    // Geometry is loaded (including ITS alignment) by AliTaskCDBconnect (A.G. 14/10/2011) 
    if (fTest) LoadGeometryFromOCDB();
  }
  fSegSPD = new AliITSsegmentationSPD();

  //slot #1
  OpenFile(1);
  fOutput = new TList();
  fOutput->SetOwner();
  
  //
  // Booking rec points related histograms
  //0
  TH2D *hLocalMapL1 = new TH2D("hLocalMapL1"," Local coordinates  - Layer 1",330,-16.5,16.5,205,-0.5,40.5); // safe limits for local coordinates in a module : z = -4,4, x = -1,1;
  hLocalMapL1->SetXTitle("direction Z [cm]");
  hLocalMapL1->SetYTitle("direction X [cm]");
  fOutput->AddLast(hLocalMapL1);
  //1
  TH2D *hLocalMapL2 = new TH2D("hLocalMapL2"," Local coordinates  - Layer 2",330,-16.5,16.5,405,-0.5,80.5);
  hLocalMapL2->SetXTitle("direction Z [cm]");
  hLocalMapL2->SetYTitle("direction X [cm]");
  fOutput->AddLast(hLocalMapL2);
  //2
  TH1F *hClusterModYield = new TH1F("hClusterModYield","Cluster yield in modules",241,-0.5,240.5);
  hClusterModYield->SetXTitle("module number");
  fOutput->AddLast(hClusterModYield);
  // 3
  TH1F *hClusterYield = new TH1F("hClusterYield","Cluster yield per chip",1200,-0.5,1199.5) ;
  hClusterYield->SetXTitle("chip key");
  fOutput->AddLast(hClusterYield); 
  //4
  TH1F *hClusterYieldOnline = new TH1F("hClusterYieldOnline","Cluster yield per chip (online coord eq*60+hs*10+chip)",1200,-0.5,1199.5) ;
  hClusterYieldOnline->SetXTitle("chip");
  fOutput->AddLast(hClusterYieldOnline);
  //5
  TH1F *hFiredChip = new TH1F("hFiredChip","Fired chip (at least one cluster)",1200,-0.5,1199.5);
  hFiredChip->SetXTitle("chip key");
  fOutput->AddLast(hFiredChip);
  //6
  TH1F *hFOFiredChip = new TH1F("hFOFiredChip","FO Fired chip ",1200,-0.5,1199.5);
  hFOFiredChip->SetXTitle("chip key");
  fOutput->AddLast(hFOFiredChip);
  //7
  TH1F *hFOgood = new TH1F("hFOgood"," FO-good (at least one cluster) ",1200,-0.5,1199.5);
  hFOgood->SetXTitle("chip key");
  fOutput->AddLast(hFOgood);
  //8
  TH2F *hFOgoodPerBC = new TH2F("hFOgoodPerBCmod4"," FO-good signals in BCmod4 ",1200,-0.5,1199.5,4,-0.5,3.5);
  hFOgoodPerBC->SetXTitle("chip key");
  fOutput->AddLast(hFOgoodPerBC);
  //9
  TH2F *hFiredChipsPerBC = new TH2F("hFiredChipsPerBCmod4"," fired chips in BCmod4 ",1200,-0.5,1199.5,4,-0.5,3.5);
  hFiredChipsPerBC->SetXTitle("chip key");
  fOutput->AddLast(hFiredChipsPerBC);
  //10
  TH1F *hFOnoisy = new TH1F("hFOnoisy","FO-noisy (no cluster)",1200,-0.5,1199.5);
  hFOnoisy->SetXTitle("chip key");
  fOutput->AddLast(hFOnoisy);

  //
  // Booking ESD related histograms
  //

  Int_t nbin, nTrMax;
  if(fHI) {
   nbin   =  800;
   nTrMax = 8000;
  } else {
  nbin =  500;
  nTrMax = 500; 
  }

  // 11
  TH1I *hTracklets = new TH1I("hNtracklets","Tracklet distribution",nbin,0,nTrMax);
  hTracklets->SetXTitle("# Tracklets");
  fOutput->AddLast(hTracklets);
  //12
  TH2F *hSPDphivsSPDeta= new TH2F("hSPDphivsSPDeta", "Tracklets - #varphi vs #eta",120,-3.,3,360,0.,2*TMath::Pi());
  hSPDphivsSPDeta->GetXaxis()->SetTitle("Pseudorapidity #eta");
  hSPDphivsSPDeta->GetYaxis()->SetTitle("#varphi [rad]");
  fOutput->AddLast(hSPDphivsSPDeta);
  //13
  TH1F *hSPDphiZpos= new TH1F("hSPDphiZpos", "Tracklets - #varphi (Z>0)",360,0.,2*TMath::Pi());
  hSPDphiZpos->SetXTitle("#varphi [rad]");
  fOutput->AddLast(hSPDphiZpos);
  //14
  TH1F *hSPDphiZneg= new TH1F("hSPDphiZneg", "Tracklets - #varphi (Z<0)",360,0.,2*TMath::Pi());
  hSPDphiZneg->SetXTitle("#varphi [rad]");
  fOutput->AddLast(hSPDphiZneg);
  //15
  TH1F *hVertexZ = new TH1F("hVertexZ","Vertex Z distribution",300,-15,15);
  hVertexZ->SetXTitle("Z Vertex [cm]");
  fOutput->AddLast(hVertexZ); 
  //16
  TH2F *hTracklVsClu1 = new TH2F("hTrackVsClu1","Tracklet Vs Clusters Layer 1",nbin/2,0,1.5*nTrMax,nbin/2,0,nTrMax);
  hTracklVsClu1->SetXTitle("clusters SPD Layer 1");
  hTracklVsClu1->SetYTitle("tracklets");
  fOutput->AddLast(hTracklVsClu1);
  //17
  TH2F *hTracklVsClu2 = new TH2F("hTrackVsClu2","Tracklet Vs Clusters Layer 2",nbin/2,0,1.5*nTrMax,nbin/2,0,nTrMax);
  hTracklVsClu2->SetXTitle("clusters SPD Layer 2");
  hTracklVsClu2->SetYTitle("tracklets");
  fOutput->AddLast(hTracklVsClu2);
  //18
  TH1I *hEventsProcessed = new TH1I("hEventsProcessed","Number of processed events",1,0,1) ;
  fOutput->AddLast(hEventsProcessed);

  PostData(1,fOutput);
}
//_________________________________________________
void AliAnalysisTaskSPD::UserExec(Option_t *)
{
  //
  // Main loop function
  //
  AliESDInputHandlerRP *hand = dynamic_cast<AliESDInputHandlerRP*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(!hand) {
    printf("No AliESDInputHandlerRP \n");  
    return;
  }

  AliESDEvent *ESD = hand->GetEvent();
  if(!ESD) {
    printf("No AliESDEvent \n");
    return;
  }

  Bool_t recP = kTRUE;
  TTree * treeRP = hand->GetTreeR("ITS");
  if(!treeRP) {
    //AliWarning("No ITS RecPoints tree ");
    recP=kFALSE;
  }

  

  // ESD related histograms
 
   const AliESDVertex *vertex = ESD->GetPrimaryVertexSPD();
   const AliMultiplicity *mult = ESD->GetMultiplicity();
 
   // Event selection
   if(!vertex) return;
   if(!vertex->GetStatus()) return;
   if(vertex->GetNContributors() < 1) return;

  ((TH1I*)fOutput->At(18))->Fill(0);
   
  ((TH1I*)fOutput->At(11))->Fill(mult->GetNumberOfTracklets());
  UInt_t bc = (UInt_t)ESD->GetBunchCrossNumber();
  for(Int_t iChipKey=0; iChipKey < 1200; iChipKey++){
    if(mult->TestFiredChipMap(iChipKey)) {
     ((TH1F*)fOutput->At(5))->Fill(iChipKey);
     if(bc>0)((TH2F*)fOutput->At(9))->Fill(iChipKey,bc%4);   
     }
    if(mult->TestFastOrFiredChips(iChipKey)) ((TH1F*)fOutput->At(6))->Fill(iChipKey);
    if(mult->TestFastOrFiredChips(iChipKey) && mult->TestFiredChipMap(iChipKey)) {
      ((TH1F*)fOutput->At(7))->Fill(iChipKey);
      if(bc>0) ((TH2F*)fOutput->At(8))->Fill(iChipKey,bc%4);
     
    }
    if(mult->TestFastOrFiredChips(iChipKey) && !mult->TestFiredChipMap(iChipKey)) ((TH1F*)fOutput->At(10))->Fill(iChipKey);
  }
  
  
  Double_t vtxPos[3] = {999, 999, 999};
  if(vertex){
    vertex->GetXYZ(vtxPos);
 
    ((TH1F*)fOutput->At(15))->Fill(vtxPos[2]);

    ((TH2F*)fOutput->At(16))->Fill(mult->GetNumberOfITSClusters(0),mult->GetNumberOfTracklets());
    ((TH2F*)fOutput->At(17))->Fill(mult->GetNumberOfITSClusters(1),mult->GetNumberOfTracklets());
 
    for(Int_t iTracklet =0; iTracklet < mult->GetNumberOfTracklets(); iTracklet++){

      Float_t phiTr= mult->GetPhi(iTracklet);
      Float_t etaTr =mult->GetEta(iTracklet);

      ((TH2F*)fOutput->At(12))->Fill(etaTr,phiTr);

      // Z pos or Z neg
      Float_t z = vtxPos[2] + 3.9 / TMath::Tan(2 * TMath::ATan(TMath::Exp(- etaTr)));
      if(z>0) ((TH1F*)fOutput->At(13))->Fill(phiTr);
      else ((TH1F*)fOutput->At(14))->Fill(phiTr);
    }
  } 
  
  
  if(recP){
  // RecPoints info

  TClonesArray statITSrec("AliITSRecPoint");
  TClonesArray *ITSCluster = &statITSrec;
  TBranch* branch=treeRP->GetBranch("ITSRecPoints");
  if(!branch) {
    printf("NO treeRP branch available. Exiting...\n");
    return;
  }

  branch->SetAddress(&ITSCluster);

  for(Int_t iMod=0;iMod<240;iMod++){
    branch->GetEvent(iMod);
    Int_t nrecp = statITSrec.GetEntries();
    if(nrecp>0) ((TH1F*)fOutput->At(2))->Fill(iMod,nrecp);
  
    for(Int_t irec=0;irec<nrecp;irec++) {
      AliITSRecPoint *recp = (AliITSRecPoint*)statITSrec.At(irec);
      Int_t lay=recp->GetLayer();
      if(lay>1) continue;

      // ----  Filling maps (local coordinates rearranged) -----
      Float_t local[3]={-1,-1};
      local[1]=recp->GetDetLocalX();
      local[0]=recp->GetDetLocalZ();
      Int_t eq = AliSPDUtils::GetOnlineEqIdFromOffline(iMod);
      Int_t hs = AliSPDUtils::GetOnlineHSFromOffline(iMod);
      Int_t row, col;
      fSegSPD->LocalToDet(0.5,local[0],row,col);
      Int_t chip = AliSPDUtils::GetOnlineChipFromOffline(iMod,col);
     
      Double_t locx, locz, equip;
      Double_t corrlocz =0;
      if(lay==0) corrlocz=local[0];
      else corrlocz = -local[0];
      // rearranging local geometry
      if(eq<10) equip=eq;
      else equip=eq-10;
      if(eq<10){
	if(chip<5) locz =corrlocz +8 +4;
	else locz = corrlocz+4;
      } else {
	if(chip<5) locz = corrlocz-8-4;
	else locz = corrlocz-4;
      }
      // filling maps
      if(lay==0){
	locx = (local[1]+1) + hs*2 +equip*4;
	((TH2D*)fOutput->At(0))->Fill(locz,locx);
	
      } else {
	locx = (local[1]+1) + (hs-2)*2 +equip*8;
	((TH2D*)fOutput->At(1))->Fill(locz,locx);
      }
      // ---- End Filling maps (local coordinates rearranged) -----
    
      ((TH1F*)fOutput->At(3))->Fill(AliSPDUtils::GetOfflineChipKeyFromOnline(eq,hs,chip));
      ((TH1F*)fOutput->At(4))->Fill(eq*60+hs*10+chip);
   
    }
  }
 }// end if rec points are available
 
  
  
  
  /* PostData(0) is taken care of by AliAnalysisTaskSE */
  PostData(1,fOutput) ;

}


//___________________________________________________________________________
void AliAnalysisTaskSPD::Terminate(Option_t*)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliAnalysisTaskSE::Terminate();
}



//___________________________________________________________________________
void AliAnalysisTaskSPD::LoadGeometryFromOCDB(){
  //method to get the gGeomanager
  // it is called at the CreatedOutputObject stage
  // to comply with the CAF environment

  if(fTest){ 
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage(fOCDBLocation.Data());
  man->SetRun(fRunNb);  
  }

  if(!AliCDBManager::Instance()){
  AliWarning("No CDB MANAGER, geometry can not be loaded");
  return;
  }

  AliCDBEntry* obj =(AliCDBManager::Instance())->Get(AliCDBPath("GRP", "Geometry", "Data"));
  if (obj)
      AliGeomManager::SetGeometry((TGeoManager*)obj->GetObject());
  AliGeomManager::GetNalignable("ITS");
  AliGeomManager::ApplyAlignObjsFromCDB("ITS"); 
}


#endif
