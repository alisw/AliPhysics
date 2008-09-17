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
//  M.Siciliano Aug 2008 QA RecPoints and HLT mode
//  INFN Torino

// --- ROOT system ---
#include <TProfile2D.h>
#include <TH2D.h>
#include <TBranch.h>
#include <TTree.h>
#include <TGaxis.h>
#include <TMath.h>
#include <TDirectory.h>
#include <TSystem.h>
// --- Standard library ---

// --- AliRoot header files ---
#include "AliITSQASDDDataMakerRec.h"
#include "AliLog.h"
#include "AliQA.h"
#include "AliQAChecker.h"
#include "AliRawReader.h"
#include "AliITSRawStream.h"
#include "AliITSRawStreamSDD.h"
#include "AliITSRawStreamSDDCompressed.h"
#include "AliITSDetTypeRec.h"
#include "AliITSRecPoint.h"
#include "AliITSgeomTGeo.h"
#include "AliITSHLTforSDD.h"
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
fDDLModuleMap(0),
fHLTMode(0),
fHLTSDD(0)
{
  //ctor used to discriminate OnLine-Offline analysis
  if(fLDC < 0 || fLDC > 4) {
	AliError("Error: LDC number out of range; return\n");
  }
  if(!fkOnline){AliInfo("Offline mode: HLT set from AliITSDetTypeRec for SDD\n");}
  else
    if(fkOnline){
      AliInfo("Online mode: HLT set from environment for SDD\n");
      SetHLTModeFromEnvironment();
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
fDDLModuleMap(0),
fHLTMode(qadm.fHLTMode),
fHLTSDD( qadm.fHLTSDD)
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
  if(!ddlMapSDD)
    {
      AliError("Calibration object retrieval failed! SDD will not be processed");
      fDDLModuleMap = NULL;
      return;
    }
  fDDLModuleMap = (AliITSDDLModuleMapSDD*)ddlMapSDD->GetObject();
  if(!cacheStatus)ddlMapSDD->SetObject(NULL);
  ddlMapSDD->SetOwner(kTRUE);
  if(!cacheStatus)
    {
      delete ddlMapSDD;
    }

  if(fkOnline==kFALSE){
    AliInfo("Offline mode: HLTforSDDobject used \n");
    AliCDBEntry *hltforSDD = AliCDBManager::Instance()->Get("ITS/Calib/HLTforSDD");
    if(!hltforSDD){
      AliError("Calibration object retrieval failed! SDD will not be processed");    
      fHLTSDD=NULL;
      return;
    }  
    fHLTSDD = (AliITSHLTforSDD*)hltforSDD->GetObject();
    if(!cacheStatus)hltforSDD->SetObject(NULL);
    hltforSDD->SetOwner(kTRUE);
    if(!cacheStatus)
      {
	delete hltforSDD;
      }
  }
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
  TH1D *h0 = new TH1D("SDDModPattern","HW Modules pattern",fgknSDDmodules,239.5,499.5);
  h0->GetXaxis()->SetTitle("Module Number");
  h0->GetYaxis()->SetTitle("Counts");
  fAliITSQADataMakerRec->Add2RawsList((new TH1D(*h0)),0+fGenOffset,kTRUE,kFALSE);
  delete h0;
  fSDDhTask++;
  if(fLDC==0 || fLDC==1 || fLDC==2){
    TH1D *h1 = new TH1D("SDDLadPatternL3","Ladder pattern L3",14,0.5,14.5);  
    h1->GetXaxis()->SetTitle("Ladder Number on Lay3");
    h1->GetYaxis()->SetTitle("Counts");
    fAliITSQADataMakerRec->Add2RawsList((new TH1D(*h1)),1+fGenOffset, kTRUE,kFALSE);
	delete h1;
    fSDDhTask++;
  }	
  if(fLDC==0 || fLDC==3 || fLDC==4){
    TH1D *h2 = new TH1D("SDDLadPatternL4","Ladder pattern L4",22,0.5,22.5);  
    h2->GetXaxis()->SetTitle("Ladder Number on Lay4");
    h2->GetYaxis()->SetTitle("Counts");
    fAliITSQADataMakerRec->Add2RawsList((new TH1D(*h2)),2+fGenOffset, kTRUE,kFALSE);
	delete h2;
    fSDDhTask++;
  }
  if(fLDC==0 || fLDC==1 || fLDC==2){
	for(Int_t i=1; i<=fgkLADDonLAY3; i++) {
      sprintf(hname0,"SDDModPattern_L3_%d",i);
      TH1D *h3 = new TH1D(hname0,hname0,6,0.5,6.5);
      h3->GetXaxis()->SetTitle("Module Number");
      h3->GetYaxis()->SetTitle("Counts");
      fAliITSQADataMakerRec->Add2RawsList((new TH1D(*h3)),i-1+3+fGenOffset, kTRUE, kFALSE);
	  delete h3;
      fSDDhTask++;
    }
  }
  if(fLDC==0 || fLDC==3 || fLDC==4){
    for(Int_t i=1; i<=fgkLADDonLAY4; i++) {
      sprintf(hname0,"SDDModPattern_L4_%d",i);
      TH1D *h4 = new TH1D(hname0,hname0,8,0.5,8.5);
      h4->GetXaxis()->SetTitle("Module Number");
      h4->GetYaxis()->SetTitle("Counts");
      fAliITSQADataMakerRec->Add2RawsList((new TH1D(*h4)),i-1+17+fGenOffset, kTRUE,kFALSE);
	  delete h4;
      fSDDhTask++;
    }
  }

  //zPhi distribution using ladder and modules numbers
      TH2D *hphil3 = new TH2D("SDDphizL3","SDD #varphiz Layer3 ",6,0.5,6.5,14,0.5,14.5);
      hphil3->GetXaxis()->SetTitle("z[#Module L3 ]");
      hphil3->GetYaxis()->SetTitle("#varphi[#Ladder L3]");
      fAliITSQADataMakerRec->Add2RawsList((new TH2D(*hphil3)),39+fGenOffset, kFALSE,kTRUE);
      delete hphil3;
      fSDDhTask++;

      TH2D *hphil4 = new TH2D("SDDphizL4","SDD #varphiz Layer4 ",8,0.5,8.5,22,0.5,22.5);
      hphil4->GetXaxis()->SetTitle("z[#Module L4]");
      hphil4->GetYaxis()->SetTitle("#varphi[#Ladder L4]");
      fAliITSQADataMakerRec->Add2RawsList((new TH2D(*hphil4)),40+fGenOffset, kFALSE,kTRUE);
      delete hphil4;
      fSDDhTask++;


  Int_t indexlast1 = 0;
  Int_t indexlast2 = 0;

  if(fkOnline) {
	fTimeBinSize = 4;
    indexlast = 0;
    index1 = 0;
    indexlast1 = fSDDhTask;
    //   cout<<"Last of the offline "<<fSDDhTask<<endl;
    indexlast2 = 0;
    char *hname[3];
    for(Int_t i=0; i<3; i++) hname[i]= new char[50];
    for(Int_t moduleSDD =0; moduleSDD<fgknSDDmodules; moduleSDD++){
		for(Int_t iside=0;iside<fgknSide;iside++){
			AliITSgeomTGeo::GetModuleId(moduleSDD+fgkmodoffset, lay, lad, det);
			sprintf(hname[0],"SDDchargeMapFSE_L%d_%d_%d_%d",lay,lad,det,iside);
			sprintf(hname[1],"SDDChargeMapForSingleEvent_L%d_%d_%d_%d",lay,lad,det,iside);
			sprintf(hname[2],"SDDhmonoDMap_L%d_%d_%d_%d",lay,lad,det,iside);
			TProfile2D *fModuleChargeMapFSE = new TProfile2D(hname[0],hname[1],256/fTimeBinSize,-0.5,255.5,256,-0.5,255.5);
			fModuleChargeMapFSE->GetXaxis()->SetTitle("Time Bin");
			fModuleChargeMapFSE->GetYaxis()->SetTitle("Anode");
			fAliITSQADataMakerRec->Add2RawsList((new TProfile2D(*fModuleChargeMapFSE)),indexlast1 + index1 + fGenOffset,kTRUE,kFALSE);
			delete fModuleChargeMapFSE;
			
			fSDDhTask++;
			index1++;	 
			indexlast2 = indexlast1 + index1;
		}
	}

    for(Int_t moduleSDD =0; moduleSDD<fgknSDDmodules; moduleSDD++){
		for(Int_t iside=0;iside<fgknSide;iside++){
			AliITSgeomTGeo::GetModuleId(moduleSDD+fgkmodoffset, lay, lad, det);
			sprintf(hname[0],"SDDchargeMap_L%d_%d_%d_%d",lay,lad,det,iside);
			sprintf(hname[1],"SDDChargeMap_L%d_%d_%d_%d",lay,lad,det,iside);
			TProfile2D *fModuleChargeMap = new TProfile2D(hname[0],hname[1],256/fTimeBinSize,-0.5,255.5,256,-0.5,255.5);
			fModuleChargeMap->GetXaxis()->SetTitle("Time Bin");
			fModuleChargeMap->GetYaxis()->SetTitle("Anode");
			fAliITSQADataMakerRec->Add2RawsList((new TProfile2D(*fModuleChargeMap)),indexlast1 + index1 + fGenOffset,kTRUE,kFALSE);
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

  rawReader->Reset();       
  AliITSRawStream *stream;
  
  if(fkOnline==kTRUE)
    {
      if(GetHLTMode()==kTRUE)
	{
	  AliInfo("Online  mode: HLT C compressed mode used for SDD\n");
	  stream = new AliITSRawStreamSDDCompressed(rawReader); }
      else{ 
	AliInfo("Online  mode: HLT A mode used for SDD\n");
	stream = new AliITSRawStreamSDD(rawReader);}     
    }
  else 
    {
      if(fHLTSDD->IsHLTmodeC()==kTRUE){
	  AliInfo("Offline  mode: HLT C compressed mode used for SDD\n");
	stream = new AliITSRawStreamSDDCompressed(rawReader);
      }else 
	{
	AliInfo("Offline  mode: HLT A mode used for SDD\n");
	stream = new AliITSRawStreamSDD(rawReader);
      }
    }
  
  //ckeck on HLT mode
                  

  
  //  AliITSRawStreamSDD s(rawReader); 
  stream->SetDDLModuleMap(fDDLModuleMap);
  
  Int_t lay, lad, det; 

  Int_t index=0;
  if(fkOnline) {
    for(Int_t moduleSDD =0; moduleSDD<fgknSDDmodules; moduleSDD++){
        for(Int_t iside=0;iside<fgknSide;iside++) {
          if(fSDDhTask > 41 + index) fAliITSQADataMakerRec->GetRawsData(41 + index +fGenOffset)->Reset();
          index++;
        }
    }
  }
  
  Int_t cnt = 0;
  Int_t ildcID = -1;
  Int_t iddl = -1;
  Int_t isddmod = -1;
  Int_t coord1, coord2, signal, moduleSDD, ioffset, iorder, activeModule, index1;
  
  while(stream->Next()) {
    ildcID = rawReader->GetLDCId();
    iddl = rawReader->GetDDLID() - fgkDDLIDshift;
    isddmod = fDDLModuleMap->GetModuleNumber(iddl,stream->GetCarlosId());
    //  cout<<"isddmod "<<isddmod<<endl;
    if(isddmod==-1){
      AliDebug(1,Form("Found module with iddl: %d, stream->GetCarlosId: %d \n",iddl,stream->GetCarlosId() ));
      continue;
    }
    if(stream->IsCompletedModule()) {
      AliDebug(1,Form("IsCompletedModule == KTRUE\n"));
      continue;
    } 
    
    coord1 = stream->GetCoord1();
    coord2 = stream->GetCoord2();
    signal = stream->GetSignal();
    
     moduleSDD = isddmod - fgkmodoffset;

    if(isddmod <fgkmodoffset|| isddmod>fgknSDDmodules+fgkmodoffset-1) {
      AliDebug(1,Form( "Module SDD = %d, resetting it to 1 \n",isddmod));
      isddmod = 1;
    }

    fAliITSQADataMakerRec->GetRawsData(0 +fGenOffset)->Fill(isddmod); 
    
    AliITSgeomTGeo::GetModuleId(isddmod, lay, lad, det);
    ioffset = 3;
    iorder = 1;
    if(lay==4) { 
      ioffset += 14;
      iorder = 2;   
    } 
    fAliITSQADataMakerRec->GetRawsData(iorder +fGenOffset)->Fill(lad);
    fAliITSQADataMakerRec->GetRawsData(ioffset+lad-1 +fGenOffset)->Fill(det); //-1 because ladder# starts from 1    
    
    if(lay==3)    fAliITSQADataMakerRec->GetRawsData(39+fGenOffset)->Fill(det,lad);
    if(lay==4) { 
      fAliITSQADataMakerRec->GetRawsData(40+fGenOffset)->Fill(det,lad);}

    Short_t iside = stream->GetChannel();
    activeModule = moduleSDD;
    index1 = activeModule * 2 + iside;
    
    if(index1<0){
      AliDebug(1,Form("Wrong index number %d - patched to 0\n",index1));
      index1 = 0;
    }
    
    //    cout<< fSDDhTask<<endl;
    if(fkOnline) {
      if(fSDDhTask > 41 + index1) {
	//	cout<<index1<<endl;
	//cout<<"Fill"<< fAliITSQADataMakerRec->GetRawsData(41 + index1 +fGenOffset)->GetName()<<endl;
        ((TProfile2D *)(fAliITSQADataMakerRec->GetRawsData(41 + index1 +fGenOffset)))->Fill(coord2, coord1, signal);
        ((TProfile2D *)(fAliITSQADataMakerRec->GetRawsData(41 + index1 + 260*2 +fGenOffset)))->Fill(coord2, coord1, signal);
      }
    }
    cnt++;
    if(!(cnt%10000)) AliDebug(1,Form(" %d raw digits read",cnt));
  }
  AliDebug(1,Form("Event completed, %d raw digits read",cnt)); 
  delete stream;
  stream = NULL; 
 }

//____________________________________________________________________________ 
void AliITSQASDDDataMakerRec::InitRecPoints()
{
  // Initialization for RECPOINTS - SDD -

 fGenOffset = (fAliITSQADataMakerRec->fRecPointsQAList)->GetEntries();
  Int_t nOnline=1;
  Int_t  nOnline2=1;
  Int_t  nOnline3=1; 
  Int_t  nOnline4=1;
  if(fkOnline)
    {
      //      cout<<"konline"<<endl;
      nOnline=4;
      nOnline2=28;
      nOnline3=64;
      nOnline4=14;
    }

  
  TH1F *h0 = new TH1F("SDDLay3TotCh","Layer 3 total charge",1000/nOnline,-0.5, 499.5); //position number 0
  h0->GetXaxis()->SetTitle("ADC value");
  h0->GetYaxis()->SetTitle("Entries");
  fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h0)), 0 +fGenOffset,kFALSE);
  delete h0;
  fSDDhTask++;
 
  TH1F *h1 = new TH1F("SDDLay4TotCh","Layer 4 total charge",1000/nOnline,-0.5, 499.5);//position number 1
  h1->GetXaxis()->SetTitle("ADC value");
  h1->GetYaxis()->SetTitle("Entries");
  fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h1)), 1 +fGenOffset,kFALSE);
  delete h1;
  fSDDhTask++;

  char hisnam[50];
  TH2F *h2 = new TH2F("SDDGlobalCoordDistribYX","YX Global Coord Distrib",5600/nOnline2,-28,28,5600/nOnline2,-28,28);//position number 2
  h2->GetYaxis()->SetTitle("Y[cm]");
  h2->GetXaxis()->SetTitle("X[cm]");
  fAliITSQADataMakerRec->Add2RecPointsList((new TH2F(*h2)),2+fGenOffset,kTRUE);
  delete h2;
  fSDDhTask++;

  TH2F *h3 = new TH2F("SDDGlobalCoordDistribRZ","RZ Global Coord Distrib",6400/nOnline3,-32,32,1400/nOnline4,12,26);//position number 3
  h3->GetYaxis()->SetTitle("R[cm]");
  h3->GetXaxis()->SetTitle("Z[cm]");
  fAliITSQADataMakerRec->Add2RecPointsList((new TH2F(*h3)),3+fGenOffset,kTRUE);
  delete h3;
  fSDDhTask++;
  
  TH2F *h4 = new TH2F("SDDGlobalCoordDistribL3PHIZ","#varphi Z Global Coord Distrib L3",6400/nOnline3,-32,32,360/nOnline,-TMath::Pi(),TMath::Pi());//position number 4
  h4->GetYaxis()->SetTitle("#phi[rad]");
  h4->GetXaxis()->SetTitle("Z[cm]");
  fAliITSQADataMakerRec->Add2RecPointsList((new TH2F(*h4)),4+fGenOffset,kFALSE);
  delete h4;
  fSDDhTask++;

  TH2F *h5 = new TH2F("SDDGlobalCoordDistribL4PHIZ","#varphi Z Global Coord Distrib L4",6400/nOnline3,-32,32,360/nOnline,-TMath::Pi(),TMath::Pi());//position number 5
  h5->GetYaxis()->SetTitle("#phi[rad]");
  h5->GetXaxis()->SetTitle("Z[cm]");
  fAliITSQADataMakerRec->Add2RecPointsList((new TH2F(*h5)),5+fGenOffset,kFALSE);
  delete h5;
  fSDDhTask++;
  
  TH1F *h6 = new TH1F("SDDModPatternRP","Modules pattern RP",fgknSDDmodules,239.5,499.5); //position number 6
  h6->GetXaxis()->SetTitle("Module number"); //spd offset = 240
  h6->GetYaxis()->SetTitle("Entries");
  fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h6)),6 +fGenOffset,kTRUE);
  delete h6;
  fSDDhTask++;
  TH1F *h7 = new TH1F("SDDLadPatternL3RP","Ladder pattern L3 RP",14,0.5,14.5);  //position number 7
  h7->GetXaxis()->SetTitle("Ladder #, Layer 3");
  h7->GetYaxis()->SetTitle("Entries");
  fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h7)),7 +fGenOffset,kTRUE);
  delete h7;
  fSDDhTask++;
  TH1F *h8 = new TH1F("SDDLadPatternL4RP","Ladder pattern L4 RP",22,0.5,22.5); //position number 8
  h8->GetXaxis()->SetTitle("Ladder #, Layer 4");
  h8->GetYaxis()->SetTitle("Entries");
  fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h8)),8 +fGenOffset,kTRUE);
  delete h8;
  fSDDhTask++;
  TH2F *h9 = new TH2F("SDDLocalCoordDistrib","Local Coord Distrib",1000/nOnline,-4,4,1000/nOnline,-4,4);//position number 9
  h9->GetXaxis()->SetTitle("X local coord, drift, cm");
  h9->GetYaxis()->SetTitle("Z local coord, anode, cm");
  fAliITSQADataMakerRec->Add2RecPointsList((new TH2F(*h9)),9 +fGenOffset,kTRUE);
  delete h9;
  fSDDhTask++;


    TH1F *h10 = new TH1F("SDDrdistrib_Layer3" ,"SDD r distribution Layer3" ,100,14.,18.);//position number 10 (L3)
    h10->GetXaxis()->SetTitle("r[cm]");
    h10->GetXaxis()->CenterTitle();
    h10->GetYaxis()->SetTitle("Entries");
    fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h10)),10 +fGenOffset,kTRUE);
    delete h10;
    fSDDhTask++;

    TH1F *h11 = new TH1F("SDDrdistrib_Layer4" ,"SDD r distribution Layer4" ,100,22.,26.);// and position number 11 (L4)
    h11->GetXaxis()->SetTitle("r[cm]");
    h11->GetXaxis()->CenterTitle();
    h11->GetYaxis()->SetTitle("Entries");
    fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h11)),11 +fGenOffset,kTRUE);
    delete h11;
    fSDDhTask++;

  for(Int_t iLay=0; iLay<=1; iLay++){
    sprintf(hisnam,"SDDphidistrib_Layer%d",iLay+3);
    TH1F *h12 = new TH1F(hisnam,hisnam,180,-TMath::Pi(),TMath::Pi());//position number 12 (L3) and position number 13 (L4)
    h12->GetXaxis()->SetTitle("#varphi[rad]");
    h12->GetXaxis()->CenterTitle();
    h12->GetYaxis()->SetTitle("Entries");
    fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h12)),iLay+12+fGenOffset,kTRUE);
    delete h12;
    fSDDhTask++;
  }

  if(fkOnline)
    {
      TH2F *h14 = new TH2F("SDDGlobalCoordDistribYXFSE","YX Global Coord Distrib FSE",5600/nOnline2,-28,28,5600/nOnline2,-28,28);//position number 14
      h14->GetYaxis()->SetTitle("Y[cm]");
      h14->GetXaxis()->SetTitle("X[cm]");
      fAliITSQADataMakerRec->Add2RecPointsList((new TH2F(*h14)),14+fGenOffset,kTRUE);
      delete h14;
      fSDDhTask++;
      
      TH2F *h15 = new TH2F("SDDGlobalCoordDistribRZFSE","RZ Global Coord Distrib FSE",Int_t(6400/nOnline3),-32,32,1400/nOnline4,12,26);//position number 15
      h15->GetYaxis()->SetTitle("R[cm]");
      h15->GetXaxis()->SetTitle("Z[cm]");
      fAliITSQADataMakerRec->Add2RecPointsList((new TH2F(*h15)),15+fGenOffset,kTRUE);
      delete h15;
      fSDDhTask++;
      
    }//online



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
  if(fkOnline)
    {
      for(Int_t i=14;i<16;i++)
	{
	  fAliITSQADataMakerRec->GetRecPointsData(i+fGenOffset)->Reset();
	}
    }
  for(Int_t module=0; module<clustersTree->GetEntries();module++){
    branchRecP->GetEvent(module);
    npoints += recpoints->GetEntries();
    AliITSgeomTGeo::GetModuleId(module, lay, lad, det);
    //printf("modnumb %d, lay %d, lad %d, det %d \n",module, lay, lad, det);
    for(Int_t j=0;j<recpoints->GetEntries();j++){
      AliITSRecPoint *recp = (AliITSRecPoint*)recpoints->At(j);    
      fAliITSQADataMakerRec->GetRecPointsData(6 +fGenOffset)->Fill(module);//modpatternrp
      recp->GetGlobalXYZ(cluglo);
      Float_t rad=TMath::Sqrt(cluglo[0]*cluglo[0]+cluglo[1]*cluglo[1]); 
      Float_t phi=TMath::ATan2(cluglo[1],cluglo[0]);
      fAliITSQADataMakerRec->GetRecPointsData(9 +fGenOffset)->Fill(recp->GetDetLocalX(),recp->GetDetLocalZ());//local distribution
      fAliITSQADataMakerRec->GetRecPointsData(2 +fGenOffset)->Fill(cluglo[0],cluglo[1]);//global distribution YX
      fAliITSQADataMakerRec->GetRecPointsData(3 +fGenOffset)->Fill(cluglo[2],rad);//global distribution rz
      if(fkOnline)
	{
	  fAliITSQADataMakerRec->GetRecPointsData(14 +fGenOffset)->Fill(cluglo[0],cluglo[1]);//global distribution YX FSE
	  fAliITSQADataMakerRec->GetRecPointsData(15 +fGenOffset)->Fill(cluglo[2],rad);//global distribution rz FSE
	}
      if(recp->GetLayer() ==2) {
	fAliITSQADataMakerRec->GetRecPointsData(0 +fGenOffset)->Fill(recp->GetQ()) ;//total charge of layer 3
	fAliITSQADataMakerRec->GetRecPointsData(7 +fGenOffset)->Fill(lad);//lad pattern layer 3
	fAliITSQADataMakerRec->GetRecPointsData(10 +fGenOffset)->Fill(rad);//r distribution layer 3
	fAliITSQADataMakerRec->GetRecPointsData(12 +fGenOffset)->Fill(phi);// phi distribution layer 3
	fAliITSQADataMakerRec->GetRecPointsData(4 +fGenOffset)->Fill(cluglo[2],phi);// phi distribution layer 3
      }
      else if(recp->GetLayer() ==3) {
	fAliITSQADataMakerRec->GetRecPointsData(1 +fGenOffset)->Fill(recp->GetQ()) ;//total charge layer 4
	fAliITSQADataMakerRec->GetRecPointsData(8 +fGenOffset)->Fill(lad);//ladpatternlayer4
	fAliITSQADataMakerRec->GetRecPointsData(11 +fGenOffset)->Fill(rad);//r distribution
	fAliITSQADataMakerRec->GetRecPointsData(13 +fGenOffset)->Fill(phi);//phi distribution
	fAliITSQADataMakerRec->GetRecPointsData(5 +fGenOffset)->Fill(cluglo[2],phi);// phi distribution layer 4
      }
    }
  }
  statRecpoints.Clear();

}

//_______________________________________________________________

void AliITSQASDDDataMakerRec::SetHLTModeFromEnvironment()
{

   Int_t  hltmode= ::atoi(gSystem->Getenv("HLT_MODE"));

   if(hltmode==1)
     {
       AliInfo("Online mode: HLT mode A selected from environment for SDD\n");
       SetHLTMode(kFALSE);
     }
   else
     if(hltmode==2)
       {
       AliInfo("Online mode: HLT mode C compressed selected from environment for SDD\n");
       SetHLTMode(kTRUE);
       }
}
