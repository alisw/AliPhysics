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
#include <TProfile2D.h>
#include <TH2D.h>
#include <TBranch.h>
#include <TTree.h>
#include <TGaxis.h>
#include <TMath.h>
// --- Standard library ---

// --- AliRoot header files ---
#include "AliITSQASDDDataMakerRec.h"
#include "AliLog.h"
#include "AliQA.h"
#include "AliQAChecker.h"
#include "AliRawReader.h"
#include "AliITSRawStreamSDD.h"
#include "AliITSRecPoint.h"
#include "AliITSgeomTGeo.h"

#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"


ClassImp(AliITSQASDDDataMakerRec)

//____________________________________________________________________________ 
AliITSQASDDDataMakerRec::AliITSQASDDDataMakerRec(AliITSQADataMakerRec *aliITSQADataMakerRec, Bool_t kMode, Short_t ldc) :
TObject(),
fAliITSQADataMakerRec(aliITSQADataMakerRec),
fkOnline(kMode),
fLDC(ldc),
fSDDhTask(0),
fGenOffset(0),
fTimeBinSize(1),
fDDLModuleMap(0)
{
  //ctor used to discriminate OnLine-Offline analysis
  if(fLDC < 0 || fLDC > 4) {
	AliError("Error: LDC number out of range; return\n");
  }
  //fDDLModuleMap=NULL;
}

//____________________________________________________________________________ 
AliITSQASDDDataMakerRec::AliITSQASDDDataMakerRec(const AliITSQASDDDataMakerRec& qadm) :
TObject(),
fAliITSQADataMakerRec(qadm.fAliITSQADataMakerRec),
fkOnline(qadm.fkOnline),
fLDC(qadm.fLDC),
fSDDhTask(qadm.fSDDhTask),
fGenOffset(qadm.fGenOffset),
fTimeBinSize(1),
fDDLModuleMap(0)
{
  //copy ctor 
  fAliITSQADataMakerRec->SetName((const char*)qadm.fAliITSQADataMakerRec->GetName()) ; 
  fAliITSQADataMakerRec->SetTitle((const char*)qadm.fAliITSQADataMakerRec->GetTitle());
  fDDLModuleMap=NULL;
}

//____________________________________________________________________________ 
AliITSQASDDDataMakerRec::~AliITSQASDDDataMakerRec(){
  // destructor
  // 	if(fDDLModuleMap) delete fDDLModuleMap;
}
//__________________________________________________________________
AliITSQASDDDataMakerRec& AliITSQASDDDataMakerRec::operator = (const AliITSQASDDDataMakerRec& qac )
{
  // Equal operator.
  this->~AliITSQASDDDataMakerRec();
  new(this) AliITSQASDDDataMakerRec(qac);
  return *this;
}

//____________________________________________________________________________ 
void AliITSQASDDDataMakerRec::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  AliDebug(1,"AliITSQADM::Start of SDD Cycle\n");
}

//____________________________________________________________________________ 
void AliITSQASDDDataMakerRec::EndOfDetectorCycle(AliQA::TASKINDEX_t /*task*/, TObjArray* /*list*/)
{
  // launch the QA checking
  AliDebug(1,"AliITSDM instantiates checker with Run(AliQA::kITS, task, list)\n"); 
}

//____________________________________________________________________________ 
void AliITSQASDDDataMakerRec::InitRaws()
{ 
  // Initialization for RAW data - SDD -
  fGenOffset = (fAliITSQADataMakerRec->fRawsQAList)->GetEntries();

  AliCDBEntry *ddlMapSDD = AliCDBManager::Instance()->Get("ITS/Calib/DDLMapSDD");
  Bool_t cacheStatus = AliCDBManager::Instance()->GetCacheFlag();

  if( !ddlMapSDD){
    AliError("Calibration object retrieval failed! SDD will not be processed");
    fDDLModuleMap = NULL;
    return;
  }  
  fDDLModuleMap = (AliITSDDLModuleMapSDD*)ddlMapSDD->GetObject();
  if(!cacheStatus)ddlMapSDD->SetObject(NULL);
  ddlMapSDD->SetOwner(kTRUE);
  if(!cacheStatus)delete ddlMapSDD;
 
  Int_t lay, lad, det;
  Int_t LAY = -1;  //, LAD = -1;
  char hname0[50];
  Int_t indexlast = 0;
  Int_t index1 = 0;

  if(fLDC == 1 || fLDC == 2) LAY = 2;
  if(fLDC == 3 || fLDC == 4) LAY = 3;
  
  if(fkOnline) {
    AliInfo("Book Online Histograms for SDD\n");
  }
  else {
    AliInfo("Book Offline Histograms for SDD\n ");
  }
  TH1D *h0 = new TH1D("ModPattern","HW Modules pattern",fgknSDDmodules,-0.5,259.5);
  h0->GetXaxis()->SetTitle("Module Number");
  h0->GetYaxis()->SetTitle("Counts");
  fAliITSQADataMakerRec->Add2RawsList((new TH1D(*h0)),0+fGenOffset, kTRUE);
  delete h0;
  fSDDhTask++;
  if(fLDC==0 || fLDC==1 || fLDC==2){
    TH1D *h1 = new TH1D("LadPatternL3","Ladder pattern L3",14,0.5,14.5);  
    h1->GetXaxis()->SetTitle("Ladder Number on Lay3");
    h1->GetYaxis()->SetTitle("Counts");
    fAliITSQADataMakerRec->Add2RawsList((new TH1D(*h1)),1+fGenOffset, kTRUE);
	delete h1;
    fSDDhTask++;
  }	
  if(fLDC==0 || fLDC==3 || fLDC==4){
    TH1D *h2 = new TH1D("LadPatternL4","Ladder pattern L4",22,0.5,22.5);  
    h2->GetXaxis()->SetTitle("Ladder Number on Lay4");
    h2->GetYaxis()->SetTitle("Counts");
    fAliITSQADataMakerRec->Add2RawsList((new TH1D(*h2)),2+fGenOffset, kTRUE);
	delete h2;
    fSDDhTask++;
  }
  if(fLDC==0 || fLDC==1 || fLDC==2){
	for(Int_t i=1; i<=fgkLADDonLAY3; i++) {
      sprintf(hname0,"ModPattern_L3_%d",i);
      TH1D *h3 = new TH1D(hname0,hname0,6,0.5,6.5);
      h3->GetXaxis()->SetTitle("Module Number");
      h3->GetYaxis()->SetTitle("Counts");
      fAliITSQADataMakerRec->Add2RawsList((new TH1D(*h3)),i-1+3+fGenOffset, kTRUE);
	  delete h3;
      fSDDhTask++;
    }
  }
  if(fLDC==0 || fLDC==3 || fLDC==4){
    for(Int_t i=1; i<=fgkLADDonLAY4; i++) {
      sprintf(hname0,"ModPattern_L4_%d",i);
	  TH1D *h4 = new TH1D(hname0,hname0,8,0.5,8.5);
      h4->GetXaxis()->SetTitle("Module Number");
      h4->GetYaxis()->SetTitle("Counts");
      fAliITSQADataMakerRec->Add2RawsList((new TH1D(*h4)),i-1+17+fGenOffset, kTRUE);
	  delete h4;
      fSDDhTask++;
    }
  }

  Int_t indexlast1 = 0;
  Int_t indexlast2 = 0;

  if(fkOnline) {
	fTimeBinSize = 4;
    indexlast = 0;
    index1 = 0;
	indexlast1 = fSDDhTask;
    indexlast2 = 0;
    char *hname[3];
    for(Int_t i=0; i<3; i++) hname[i]= new char[50];
    for(Int_t moduleSDD =0; moduleSDD<fgknSDDmodules; moduleSDD++){
		for(Int_t iside=0;iside<fgknSide;iside++){
			AliITSgeomTGeo::GetModuleId(moduleSDD+fgkmodoffset, lay, lad, det);
			sprintf(hname[0],"chargeMapFSE_L%d_%d_%d_%d",lay,lad,det,iside);
			sprintf(hname[1],"ChargeMapForSingleEvent_L%d_%d_%d_%d",lay,lad,det,iside);
			sprintf(hname[2],"hmonoDMap_L%d_%d_%d_%d",lay,lad,det,iside);
			TProfile2D *fModuleChargeMapFSE = new TProfile2D(hname[0],hname[1],256/fTimeBinSize,-0.5,255.5,256,-0.5,255.5);
			fModuleChargeMapFSE->GetXaxis()->SetTitle("Time Bin");
			fModuleChargeMapFSE->GetYaxis()->SetTitle("Anode");
			fAliITSQADataMakerRec->Add2RawsList((new TProfile2D(*fModuleChargeMapFSE)),indexlast1 + index1 + fGenOffset);
			delete fModuleChargeMapFSE;
			
			fSDDhTask++;
			index1++;	 
			indexlast2 = indexlast1 + index1;
		}
	}

    for(Int_t moduleSDD =0; moduleSDD<fgknSDDmodules; moduleSDD++){
		for(Int_t iside=0;iside<fgknSide;iside++){
			AliITSgeomTGeo::GetModuleId(moduleSDD+fgkmodoffset, lay, lad, det);
			sprintf(hname[0],"chargeMap_L%d_%d_%d_%d",lay,lad,det,iside);
			sprintf(hname[1],"ChargeMap_L%d_%d_%d_%d",lay,lad,det,iside);
			TProfile2D *fModuleChargeMap = new TProfile2D(hname[0],hname[1],256/fTimeBinSize,-0.5,255.5,256,-0.5,255.5);
			fModuleChargeMap->GetXaxis()->SetTitle("Time Bin");
			fModuleChargeMap->GetYaxis()->SetTitle("Anode");
			fAliITSQADataMakerRec->Add2RawsList((new TProfile2D(*fModuleChargeMap)),indexlast1 + index1 + fGenOffset);
			delete fModuleChargeMap;

			fSDDhTask++;
			index1++;	 
			indexlast2 = indexlast1 + index1;
		}
	}
  
}  // kONLINE

  AliDebug(1,Form("%d SDD Raws histograms booked\n",fSDDhTask));
}


//____________________________________________________________________________
void AliITSQASDDDataMakerRec::MakeRaws(AliRawReader* rawReader)
{ 
  // Fill QA for RAW - SDD -
  if(!fDDLModuleMap){
    AliError("SDD DDL module map not available - skipping SDD QA");
    return;
  }
  if(rawReader->GetType() != 7) return;  // skips non physical triggers
  AliDebug(1,"entering MakeRaws\n");                 
  rawReader->SelectEquipment(17,fgkeqOffset,fgkeqOffset + fgkDDLidRange); 

  rawReader->Reset();                         
  AliITSRawStreamSDD s(rawReader); 
  s.SetDDLModuleMap(fDDLModuleMap);
  Int_t lay, lad, det; 

  Int_t index=0;
  if(fkOnline) {
    for(Int_t moduleSDD =0; moduleSDD<fgknSDDmodules; moduleSDD++){
        for(Int_t iside=0;iside<fgknSide;iside++) {
          if(fSDDhTask > 39 + index) fAliITSQADataMakerRec->GetRawsData(39 + index +fGenOffset)->Reset();
          index++;
        }
    }
  }
  Int_t cnt = 0;
  Int_t ildcID = -1;
  Int_t iddl = -1;
  Int_t isddmod = -1;
  Int_t coord1, coord2, signal, moduleSDD, ioffset, iorder, activeModule, index1;
  while(s.Next()) {
    ildcID = rawReader->GetLDCId();
    iddl = rawReader->GetDDLID() - fgkDDLIDshift;
    isddmod = s.GetModuleNumber(iddl,s.GetCarlosId());
    if(isddmod==-1){
      AliDebug(1,Form("Found module with iddl: %d, s.GetCarlosId: %d \n",iddl,s.GetCarlosId() ));
      continue;
    }
    if(s.IsCompletedModule()) {
      AliDebug(1,Form("IsCompletedModule == KTRUE\n"));
      continue;
    } 
    
    coord1 = s.GetCoord1();
    coord2 = s.GetCoord2();
    signal = s.GetSignal();
    
    moduleSDD = isddmod - fgkmodoffset;
    if(moduleSDD < 0 || moduleSDD>fgknSDDmodules) {
      AliDebug(1,Form( "Module SDD = %d, resetting it to 1 \n",moduleSDD));
      moduleSDD = 1;
    }
    fAliITSQADataMakerRec->GetRawsData(0 +fGenOffset)->Fill(moduleSDD); 
    
    AliITSgeomTGeo::GetModuleId(isddmod, lay, lad, det);
    ioffset = 3;
    iorder = 1;
    if(lay==4) { 
      ioffset += 14;
      iorder = 2;   
    } 
    fAliITSQADataMakerRec->GetRawsData(iorder +fGenOffset)->Fill(lad);
    fAliITSQADataMakerRec->GetRawsData(ioffset+lad-1 +fGenOffset)->Fill(det); //-1 because ladder# starts from 1    
    
    Short_t iside = s.GetChannel();
    activeModule = moduleSDD;
    index1 = activeModule * 2 + iside;
    
    if(index1<0){
      AliDebug(1,Form("Wrong index number %d - patched to 0\n",index1));
      index1 = 0;
    }
    
    if(fkOnline) {
      if(fSDDhTask > 39 + index1) {
        ((TProfile2D *)(fAliITSQADataMakerRec->GetRawsData(39 + index1 +fGenOffset)))->Fill(coord2, coord1, signal);
        ((TProfile2D *)(fAliITSQADataMakerRec->GetRawsData(39 + index1 + 260*2 +fGenOffset)))->Fill(coord2, coord1, signal);
      }
    }
    cnt++;
    if(!(cnt%10000)) AliDebug(1,Form(" %d raw digits read",cnt));
  }
  AliDebug(1,Form("Event completed, %d raw digits read",cnt)); 
 }

//____________________________________________________________________________ 
void AliITSQASDDDataMakerRec::InitRecPoints()
{
  // Initialization for RECPOINTS - SDD -
  fGenOffset = (fAliITSQADataMakerRec->fRecPointsQAList)->GetEntries();
  
  TH1F *h0 = new TH1F("Lay3TotCh","Layer 3 total charge",1000,-0.5, 499.5);
  h0->GetXaxis()->SetTitle("ADC value");
  h0->GetYaxis()->SetTitle("Entries");
  fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h0)), 0 +fGenOffset);
  delete h0;
  fSDDhTask++;
 
  TH1F *h1 = new TH1F("Lay4TotCh","Layer 4 total charge",1000,-0.5, 499.5);
  h1->GetXaxis()->SetTitle("ADC value");
  h1->GetYaxis()->SetTitle("Entries");
  fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h1)), 1 +fGenOffset);
  delete h1;
  fSDDhTask++;

    
  char hisnam[50];
  for(Int_t i=1; i<=3; i++){
    sprintf(hisnam,"Charge_L3_Strip%d",i);
    TH1F *h2 = new TH1F(hisnam,hisnam,1000,-0.5, 499.5);
    fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h2)),i+1 +fGenOffset);
	delete h2;
    fSDDhTask++;
  }
  
  for(Int_t i=1; i<=4; i++){
    sprintf(hisnam,"Charge_L4_Strip%d",i);
    TH1F *h3 = new TH1F(hisnam,hisnam,1000,-0.5, 499.5);
    fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h3)),i+4 +fGenOffset);
	delete h3;
    fSDDhTask++;
  }
  
  TH1F *h4 = new TH1F("ModPatternRP","Modules pattern RP",fgknSDDmodules,239.5,499.5); 
  h4->GetXaxis()->SetTitle("Module number");
  h4->GetYaxis()->SetTitle("Entries");
  fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h4)),9 +fGenOffset);
  delete h4;
  fSDDhTask++;
  TH1F *h5 = new TH1F("ModPatternL3 RP","Ladder pattern L3 RP",14,0.5,14.5);  
  h5->GetXaxis()->SetTitle("Ladder #, Layer 3");
  h5->GetYaxis()->SetTitle("Entries");
  fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h5)),10 +fGenOffset);
  delete h5;
  fSDDhTask++;
  TH1F *h6 = new TH1F("ModPatternL4 RP","Ladder pattern L4 RP",22,0.5,22.5); 
  h6->GetXaxis()->SetTitle("Ladder #, Layer 4");
  h6->GetYaxis()->SetTitle("Entries");
  fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h6)),11 +fGenOffset);
  delete h6;
  fSDDhTask++;
  TH2F *h7 = new TH2F("Local Coord Distrib","Local Coord Distrib",1000,-4,4,1000,-4,4);
  h7->GetXaxis()->SetTitle("X local coord, drift, cm");
  h7->GetYaxis()->SetTitle("Z local coord, anode, cm");
  fAliITSQADataMakerRec->Add2RecPointsList((new TH2F(*h7)),12 +fGenOffset);
  delete h7;
  fSDDhTask++;
  TH2F *h8 = new TH2F("Global Coord Distrib","Global Coord Distrib",6000,-30,30,6000,-30,30);
  h8->GetYaxis()->SetTitle("Y glob coord, cm");
  h8->GetXaxis()->SetTitle("X glob coord, cm");
  fAliITSQADataMakerRec->Add2RecPointsList((new TH2F(*h8)),13 +fGenOffset);
  delete h8;
  fSDDhTask++;
   
  for(Int_t iLay=0; iLay<=1; iLay++){
    sprintf(hisnam,"hr_Layer%d",iLay+3);
    TH1F *h9 = new TH1F(hisnam,hisnam,100,10.,30.);
    h9->GetXaxis()->SetTitle("r (cm)");
    h9->GetXaxis()->CenterTitle();
    h9->GetYaxis()->SetTitle("Entries");
    fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h9)),iLay+14 +fGenOffset);
	delete h9;
    fSDDhTask++;
  }

  for(Int_t iLay=0; iLay<=1; iLay++){
    sprintf(hisnam,"hphi_Layer%d",iLay+3);
    TH1F *h10 = new TH1F(hisnam,hisnam,100,-TMath::Pi(),TMath::Pi());
    h10->GetXaxis()->SetTitle("#varphi (rad)");
    h10->GetXaxis()->CenterTitle();
    h10->GetYaxis()->SetTitle("Entries");
    fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h10)),iLay+16 +fGenOffset);
	delete h10;
    fSDDhTask++;
  }

  AliDebug(1,Form("%d SDD Recs histograms booked\n",fSDDhTask));
}

//____________________________________________________________________________ 
void AliITSQASDDDataMakerRec::MakeRecPoints(TTree * clustersTree)
{
  // Fill QA for RecPoints - SDD -
  Int_t lay, lad, det; 
  TBranch *branchRecP = clustersTree->GetBranch("ITSRecPoints");
  if (!branchRecP) { 
    AliError("can't get the branch with the ITS clusters !");
    return;
  }
  static TClonesArray statRecpoints("AliITSRecPoint") ;
  TClonesArray *recpoints = &statRecpoints;
  branchRecP->SetAddress(&recpoints);
  Int_t npoints = 0;      
  Float_t cluglo[3]={0.,0.,0.}; 
  for(Int_t module=0; module<clustersTree->GetEntries();module++){
    branchRecP->GetEvent(module);
    npoints += recpoints->GetEntries();
    AliITSgeomTGeo::GetModuleId(module, lay, lad, det);
    //printf("modnumb %d, lay %d, lad %d, det %d \n",module, lay, lad, det);
    
    for(Int_t j=0;j<recpoints->GetEntries();j++){
      AliITSRecPoint *recp = (AliITSRecPoint*)recpoints->At(j);    
      fAliITSQADataMakerRec->GetRecPointsData(9 +fGenOffset)->Fill(module);
      recp->GetGlobalXYZ(cluglo);
      Float_t rad=TMath::Sqrt(cluglo[0]*cluglo[0]+cluglo[1]*cluglo[1]); 
      Float_t phi=TMath::ATan2(cluglo[1],cluglo[0]);
      if(recp->GetLayer() ==2) {
	fAliITSQADataMakerRec->GetRecPointsData(0 +fGenOffset)->Fill(recp->GetQ()) ;
	fAliITSQADataMakerRec->GetRecPointsData(10 +fGenOffset)->Fill(lad);
	fAliITSQADataMakerRec->GetRecPointsData(14 +fGenOffset)->Fill(rad);
	fAliITSQADataMakerRec->GetRecPointsData(16 +fGenOffset)->Fill(phi);
	fAliITSQADataMakerRec->GetRecPointsData(9 +fGenOffset)->Fill(module);
	fAliITSQADataMakerRec->GetRecPointsData(12 +fGenOffset)->Fill(recp->GetDetLocalX(),recp->GetDetLocalZ());
	fAliITSQADataMakerRec->GetRecPointsData(13 +fGenOffset)->Fill(cluglo[0],cluglo[1]);
      }
      else if(recp->GetLayer() ==3) {
	fAliITSQADataMakerRec->GetRecPointsData(1 +fGenOffset)->Fill(recp->GetQ()) ;
	fAliITSQADataMakerRec->GetRecPointsData(11 +fGenOffset)->Fill(lad);
	fAliITSQADataMakerRec->GetRecPointsData(15 +fGenOffset)->Fill(rad);
	fAliITSQADataMakerRec->GetRecPointsData(17 +fGenOffset)->Fill(phi);
	fAliITSQADataMakerRec->GetRecPointsData(9 +fGenOffset)->Fill(module);
	fAliITSQADataMakerRec->GetRecPointsData(12 +fGenOffset)->Fill(recp->GetDetLocalX(),recp->GetDetLocalZ());
	fAliITSQADataMakerRec->GetRecPointsData(13 +fGenOffset)->Fill(cluglo[0],cluglo[1]);
      }
    }
  }
  statRecpoints.Clear();
}

