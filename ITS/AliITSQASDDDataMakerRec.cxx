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
#include "AliQAv1.h"
#include "AliQAChecker.h"
#include "AliRawReader.h"
#include "AliITSRawStream.h"
#include "AliITSRawStreamSDD.h"
#include "AliITSRawStreamSDDCompressed.h"
#include "AliITSDetTypeRec.h"
#include "AliITSdigit.h"
#include "AliITSRecPoint.h"
#include "AliITSgeomTGeo.h"
#include "AliITSHLTforSDD.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "Riostream.h"
#include "AliITSdigitSDD.h"
#include "AliITS.h"
#include "AliRunLoader.h"
#include "AliITSLoader.h"
#include "AliITSDetTypeRec.h"



ClassImp(AliITSQASDDDataMakerRec)

//____________________________________________________________________________ 
AliITSQASDDDataMakerRec::AliITSQASDDDataMakerRec(AliITSQADataMakerRec *aliITSQADataMakerRec, Bool_t kMode, Short_t ldc) :
TObject(),
fAliITSQADataMakerRec(aliITSQADataMakerRec),
fkOnline(kMode),
fLDC(ldc),
fSDDhRawsTask(0),
fSDDhDigitsTask(0),
fSDDhRecPointsTask(0),
fGenRawsOffset(0),
fGenDigitsOffset(0),
fGenRecPointsOffset(0),
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
fSDDhRawsTask(qadm.fSDDhRawsTask),
fSDDhDigitsTask(qadm.fSDDhDigitsTask),
fSDDhRecPointsTask(qadm.fSDDhRecPointsTask),
fGenRawsOffset(qadm.fGenRawsOffset),
fGenDigitsOffset(qadm.fGenDigitsOffset),
fGenRecPointsOffset(qadm.fGenRecPointsOffset),
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
  AliDebug(AliQAv1::GetQADebugLevel(),"AliITSQADM::Start of SDD Cycle\n");
}

//____________________________________________________________________________ 
void AliITSQASDDDataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t /*task*/, TObjArray* /*list*/)
{
  // launch the QA checking
  AliDebug(AliQAv1::GetQADebugLevel(),"AliITSDM instantiates checker with Run(AliQAv1::kITS, task, list)\n"); 
}

//____________________________________________________________________________ 
void AliITSQASDDDataMakerRec::InitRaws()
{ 
  // Initialization for RAW data - SDD -
  const Bool_t expert   = kTRUE ; 
  const Bool_t saveCorr = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  fGenRawsOffset = (fAliITSQADataMakerRec->fRawsQAList[AliRecoParam::kDefault])->GetEntries();
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
  Int_t indexlast = 0;
  Int_t index1 = 0;

  if(fkOnline) 
    {
      AliInfo("Book Online Histograms for SDD\n");
    }
  else 
    {
      AliInfo("Book Offline Histograms for SDD\n ");
    }
  TH1D *h0 = new TH1D("SDDModPattern","HW Modules pattern",fgknSDDmodules,239.5,499.5); //0
  h0->GetXaxis()->SetTitle("Module Number");
  h0->GetYaxis()->SetTitle("Counts");
  fAliITSQADataMakerRec->Add2RawsList((new TH1D(*h0)),0+fGenRawsOffset, expert, !image, !saveCorr);
  delete h0;
  fSDDhRawsTask++;
  
  //zPhi distribution using ladder and modules numbers
  TH2D *hphil3 = new TH2D("SDDphizL3","SDD #varphiz Layer3 ",6,0.5,6.5,14,0.5,14.5);
  hphil3->GetXaxis()->SetTitle("z[#Module L3 ]");
  hphil3->GetYaxis()->SetTitle("#varphi[#Ladder L3]");
  fAliITSQADataMakerRec->Add2RawsList((new TH2D(*hphil3)),1+fGenRawsOffset, !expert, image, saveCorr); 
  delete hphil3;
  fSDDhRawsTask++;
  
  TH2D *hphil4 = new TH2D("SDDphizL4","SDD #varphiz Layer4 ",8,0.5,8.5,22,0.5,22.5); 
  hphil4->GetXaxis()->SetTitle("z[#Module L4]");
  hphil4->GetYaxis()->SetTitle("#varphi[#Ladder L4]");
  fAliITSQADataMakerRec->Add2RawsList((new TH2D(*hphil4)),2+fGenRawsOffset, !expert, image, saveCorr); 
  delete hphil4;
  fSDDhRawsTask++;
  

  if(fkOnline) 
    {

      //DDL Pattern 
      TH2D *hddl = new TH2D("SDDDDLPattern","SDD DDL Pattern ",24,-0.5,23.5,24,-0.5,23.5); 
      hddl->GetXaxis()->SetTitle("Channel");
      hddl->GetYaxis()->SetTitle("#DDL");
      fAliITSQADataMakerRec->Add2RawsList((new TH2D(*hddl)),3+fGenRawsOffset, expert, !image, !saveCorr);
      delete hddl;
      fSDDhRawsTask++;
      Int_t indexlast1 = 0;
  
      fTimeBinSize = 4;
      indexlast = 0;
      index1 = 0;
      indexlast1 = fSDDhRawsTask;
      char *hname[3];
      for(Int_t i=0; i<3; i++) hname[i]= new char[50];
      for(Int_t moduleSDD =0; moduleSDD<fgknSDDmodules; moduleSDD++){
	for(Int_t iside=0;iside<fgknSide;iside++){
	  AliITSgeomTGeo::GetModuleId(moduleSDD+fgkmodoffset, lay, lad, det);
	  sprintf(hname[0],"SDDchargeMapFSE_L%d_%d_%d_%d",lay,lad,det,iside);
	  sprintf(hname[1],"SDDChargeMapForSingleEvent_L%d_%d_%d_%d",lay,lad,det,iside);
	  //	  sprintf(hname[2],"SDDhmonoDMap_L%d_%d_%d_%d",lay,lad,det,iside);
	  TProfile2D *fModuleChargeMapFSE = new TProfile2D(hname[0],hname[1],256/fTimeBinSize,-0.5,255.5,256,-0.5,255.5);
	  fModuleChargeMapFSE->GetXaxis()->SetTitle("Time Bin");
	  fModuleChargeMapFSE->GetYaxis()->SetTitle("Anode");
	  fAliITSQADataMakerRec->Add2RawsList((new TProfile2D(*fModuleChargeMapFSE)),indexlast1 + index1 + fGenRawsOffset, expert, !image, !saveCorr);
	  delete fModuleChargeMapFSE;
	  
	  fSDDhRawsTask++;
	  index1++;	 
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
	  fAliITSQADataMakerRec->Add2RawsList((new TProfile2D(*fModuleChargeMap)),indexlast1 + index1 + fGenRawsOffset, expert, !image, !saveCorr);
	  delete fModuleChargeMap;
	  
	  fSDDhRawsTask++;
	  index1++;	 
	}
      }
      
    }  // kONLINE
  
  AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SDD Raws histograms booked\n",fSDDhRawsTask));
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
  AliDebug(AliQAv1::GetQADebugLevel(),"entering MakeRaws\n");                 
  
  rawReader->Reset();       
  AliITSRawStream *stream;
  
  if(fkOnline==kTRUE)
    {
      if(GetHLTMode()==kTRUE)
	{
	  //AliInfo("Online  mode: HLT C compressed mode used for SDD\n");
	  stream = new AliITSRawStreamSDDCompressed(rawReader); }
      else{ 
	//AliInfo("Online  mode: HLT A mode used for SDD\n");
	stream = new AliITSRawStreamSDD(rawReader);}     
    }
  else 
    {
      if(fHLTSDD->IsHLTmodeC()==kTRUE){
	//AliInfo("Offline  mode: HLT C compressed mode used for SDD\n");
	stream = new AliITSRawStreamSDDCompressed(rawReader);
      }else 
	{
	  //AliInfo("Offline  mode: HLT A mode used for SDD\n");
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
	if(fSDDhRawsTask > 4 + index) fAliITSQADataMakerRec->GetRawsData(4 + index +fGenRawsOffset)->Reset();   
	// 4  because the 2D histos for single events start after the fourth position
	index++;
      }
    }
  }
  
  Int_t cnt = 0;
  Int_t ildcID = -1;
  Int_t iddl = -1;
  Int_t isddmod = -1;
  Int_t coord1, coord2, signal, moduleSDD, activeModule, index1; 
  
  while(stream->Next()) {
    ildcID = rawReader->GetLDCId();
    iddl = rawReader->GetDDLID() - fgkDDLIDshift;
    
    isddmod = fDDLModuleMap->GetModuleNumber(iddl,stream->GetCarlosId());
    if(isddmod==-1){
      AliDebug(AliQAv1::GetQADebugLevel(),Form("Found module with iddl: %d, stream->GetCarlosId: %d \n",iddl,stream->GetCarlosId()));
      continue;
    }
    if(stream->IsCompletedModule()) {
      AliDebug(AliQAv1::GetQADebugLevel(),Form("IsCompletedModule == KTRUE\n"));
      continue;
    } 
    if(stream->IsCompletedDDL()) {
      AliDebug(AliQAv1::GetQADebugLevel(),Form("IsCompletedDDL == KTRUE\n"));
      continue;
    } 
    
    coord1 = stream->GetCoord1();
    coord2 = stream->GetCoord2();
    signal = stream->GetSignal();
    
    moduleSDD = isddmod - fgkmodoffset;
    
    if(isddmod <fgkmodoffset|| isddmod>fgknSDDmodules+fgkmodoffset-1) {
      AliDebug(AliQAv1::GetQADebugLevel(),Form( "Module SDD = %d, resetting it to 1 \n",isddmod));
      isddmod = 1;
    }
    
    AliITSgeomTGeo::GetModuleId(isddmod, lay, lad, det);

    
    fAliITSQADataMakerRec->GetRawsData( 0 + fGenRawsOffset )->Fill(isddmod);   
    
    if(lay==3)    fAliITSQADataMakerRec->GetRawsData(1+fGenRawsOffset)->Fill(det,lad); 
    if(lay==4) { 
      fAliITSQADataMakerRec->GetRawsData(2+fGenRawsOffset)->Fill(det,lad);}  
    
    Short_t iside = stream->GetChannel();
    

    

    if(fkOnline) {

      fAliITSQADataMakerRec->GetRawsData(3+fGenRawsOffset)->Fill(2*(stream->GetCarlosId())+iside,iddl);

      activeModule = moduleSDD;
      index1 = activeModule * 2 + iside;
      
      if(index1<0){
        AliDebug(AliQAv1::GetQADebugLevel(),Form("Wrong index number %d - patched to 0\n",index1));
	index1 = 0;
      }      
      fAliITSQADataMakerRec->GetRawsData(3+fGenRawsOffset)->Fill(2*(stream->GetCarlosId())+iside,iddl);
      if(fSDDhRawsTask > 4 + index1) {                                  
        ((TProfile2D *)(fAliITSQADataMakerRec->GetRawsData(4 + index1 +fGenRawsOffset)))->Fill(coord2, coord1, signal);     
        ((TProfile2D *)(fAliITSQADataMakerRec->GetRawsData(4 + index1 + 260*2 +fGenRawsOffset)))->Fill(coord2, coord1, signal); 
      }
    }
    cnt++;
    if(!(cnt%10000)) AliDebug(AliQAv1::GetQADebugLevel(),Form(" %d raw digits read",cnt));
  }
  AliDebug(AliQAv1::GetQADebugLevel(),Form("Event completed, %d raw digits read",cnt)); 
  delete stream;
  stream = NULL; 
}

//____________________________________________________________________________ 
void AliITSQASDDDataMakerRec::InitDigits()
{ 

  //printf(" ======================================================> Init digits\n " );
  // Initialization for DIGIT data - SDD -  
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ;
  
  fGenDigitsOffset = (fAliITSQADataMakerRec->fDigitsQAList[AliRecoParam::kDefault])->GetEntries();
  //fSDDhTask must be incremented by one unit every time a histogram is ADDED to the QA List
  TH1F* h0=new TH1F("SDD DIGITS Module Pattern","SDD DIGITS Module Pattern",260,239.5,499.5);       //hmod
  h0->GetXaxis()->SetTitle("SDD Module Number");
  h0->GetYaxis()->SetTitle("# DIGITS");
  fAliITSQADataMakerRec->Add2DigitsList(h0,fGenDigitsOffset, !expert, image);
  fSDDhDigitsTask ++;
  TH1F* h1=new TH1F("SDD Anode Distribution","DIGITS Anode Distribution",512,-0.5,511.5);      //hanocc
  h1->GetXaxis()->SetTitle("Anode Number");
  h1->GetYaxis()->SetTitle("# DIGITS");
  fAliITSQADataMakerRec->Add2DigitsList(h1,1+fGenDigitsOffset, !expert, image);
  fSDDhDigitsTask ++;
  TH1F* h2=new TH1F("SDD Tbin Distribution","DIGITS Tbin Distribution",256,-0.5,255.5);      //htbocc
  h2->GetXaxis()->SetTitle("Tbin Number");
  h2->GetYaxis()->SetTitle("# DIGITS");
  fAliITSQADataMakerRec->Add2DigitsList(h2,2+fGenDigitsOffset, !expert, image);
  fSDDhDigitsTask ++;
  TH1F* h3=new TH1F("SDD ADC Counts Distribution","DIGITS ADC Counts Distribution",200,0.,1024.);          //hsig
  h3->GetXaxis()->SetTitle("ADC Value");
  h3->GetYaxis()->SetTitle("# DIGITS");
  fAliITSQADataMakerRec->Add2DigitsList(h3,3+fGenDigitsOffset, !expert, image);
  fSDDhDigitsTask ++;
  AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SDD Digits histograms booked\n",fSDDhDigitsTask));
}

//____________________________________________________________________________
void AliITSQASDDDataMakerRec::MakeDigits(TTree * digits)
{ 
  //printf(" ======================================================> make digits\n " );
  // Fill QA for DIGIT - SDD -
  //AliITS *fITS  = (AliITS*)gAlice->GetModule("ITS");
  //fITS->SetTreeAddress();
  //TClonesArray *iITSdigits  = fITS->DigitsAddress(1);


  TBranch *branchD = digits->GetBranch("ITSDigitsSDD");

  if (!branchD) {
    AliError("can't get the branch with the ITS SDD digits !");
    return;
  }

  static TClonesArray statDigits("AliITSdigitSDD");
  TClonesArray *iITSdigits = &statDigits;
  branchD->SetAddress(&iITSdigits);

  for(Int_t i=0; i<260; i++){
    Int_t nmod=i+240;
    digits->GetEvent(nmod);
    Int_t ndigits = iITSdigits->GetEntries();
    fAliITSQADataMakerRec->GetDigitsData(fGenDigitsOffset)->Fill(nmod,ndigits);
    //printf(" Filled:  =======================================> %s \t %i \t %i \n",fAliITSQADataMakerRec->GetDigitsData(fGenDigitsOffset)->GetName(), nmod, ndigits );
    for (Int_t idig=0; idig<ndigits; idig++) {
      AliITSdigit *dig=(AliITSdigit*)iITSdigits->UncheckedAt(idig);
      Int_t iz=dig->GetCoord1();  // cell number z
      Int_t ix=dig->GetCoord2();  // cell number x
      Int_t sig=dig->GetSignal();
      fAliITSQADataMakerRec->GetDigitsData(1+fGenDigitsOffset)->Fill(iz);
      fAliITSQADataMakerRec->GetDigitsData(2+fGenDigitsOffset)->Fill(ix);
      fAliITSQADataMakerRec->GetDigitsData(3+fGenDigitsOffset)->Fill(sig);
    }
  }

}

//____________________________________________________________________________ 
void AliITSQASDDDataMakerRec::InitRecPoints()
{
  // Initialization for RECPOINTS - SDD -
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  fGenRecPointsOffset = (fAliITSQADataMakerRec->fRecPointsQAList[AliRecoParam::kDefault])->GetEntries();

  Int_t nOnline=1;
  Int_t  nOnline2=1;
  Int_t  nOnline3=1; 
  Int_t  nOnline4=1;
  if(fkOnline)
    {
      nOnline=4;
      nOnline2=28;
      nOnline3=64;
      nOnline4=14;
    }

  
  TH1F *h0 = new TH1F("SDDLay3TotCh","Layer 3 total charge",1000/nOnline,-0.5, 499.5); //position number 0
  h0->GetXaxis()->SetTitle("ADC value");
  h0->GetYaxis()->SetTitle("Entries");
  fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h0)), 0 +fGenRecPointsOffset, !expert, image);
  delete h0;
  fSDDhRecPointsTask++;
 
  TH1F *h1 = new TH1F("SDDLay4TotCh","Layer 4 total charge",1000/nOnline,-0.5, 499.5);//position number 1
  h1->GetXaxis()->SetTitle("ADC value");
  h1->GetYaxis()->SetTitle("Entries");
  fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h1)), 1 +fGenRecPointsOffset, !expert, image);
  delete h1;
  fSDDhRecPointsTask++;

  char hisnam[50];
  TH2F *h2 = new TH2F("SDDGlobalCoordDistribYX","YX Global Coord Distrib",5600/nOnline2,-28,28,5600/nOnline2,-28,28);//position number 2
  h2->GetYaxis()->SetTitle("Y[cm]");
  h2->GetXaxis()->SetTitle("X[cm]");
  fAliITSQADataMakerRec->Add2RecPointsList((new TH2F(*h2)),2+fGenRecPointsOffset, expert, !image);
  delete h2;
  fSDDhRecPointsTask++;

  TH2F *h3 = new TH2F("SDDGlobalCoordDistribRZ","RZ Global Coord Distrib",6400/nOnline3,-32,32,1400/nOnline4,12,26);//position number 3
  h3->GetYaxis()->SetTitle("R[cm]");
  h3->GetXaxis()->SetTitle("Z[cm]");
  fAliITSQADataMakerRec->Add2RecPointsList((new TH2F(*h3)),3+fGenRecPointsOffset, expert, !image);
  delete h3;
  fSDDhRecPointsTask++;
  
  TH2F *h4 = new TH2F("SDDGlobalCoordDistribL3PHIZ","#varphi Z Global Coord Distrib L3",6400/nOnline3,-32,32,360/nOnline,-TMath::Pi(),TMath::Pi());//position number 4
  h4->GetYaxis()->SetTitle("#phi[rad]");
  h4->GetXaxis()->SetTitle("Z[cm]");
  fAliITSQADataMakerRec->Add2RecPointsList((new TH2F(*h4)),4+fGenRecPointsOffset, !expert, image);
  delete h4;
  fSDDhRecPointsTask++;

  TH2F *h5 = new TH2F("SDDGlobalCoordDistribL4PHIZ","#varphi Z Global Coord Distrib L4",6400/nOnline3,-32,32,360/nOnline,-TMath::Pi(),TMath::Pi());//position number 5
  h5->GetYaxis()->SetTitle("#phi[rad]");
  h5->GetXaxis()->SetTitle("Z[cm]");
  fAliITSQADataMakerRec->Add2RecPointsList((new TH2F(*h5)),5+fGenRecPointsOffset, !expert, image);
  delete h5;
  fSDDhRecPointsTask++;
  
  TH1F *h6 = new TH1F("SDDModPatternRP","Modules pattern RP",fgknSDDmodules,239.5,499.5); //position number 6
  h6->GetXaxis()->SetTitle("Module number"); //spd offset = 240
  h6->GetYaxis()->SetTitle("Entries");
  fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h6)),6 +fGenRecPointsOffset, expert, !image);
  delete h6;
  fSDDhRecPointsTask++;
  TH1F *h7 = new TH1F("SDDLadPatternL3RP","Ladder pattern L3 RP",14,0.5,14.5);  //position number 7
  h7->GetXaxis()->SetTitle("Ladder #, Layer 3");
  h7->GetYaxis()->SetTitle("Entries");
  fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h7)),7 +fGenRecPointsOffset, expert, !image);
  delete h7;
  fSDDhRecPointsTask++;
  TH1F *h8 = new TH1F("SDDLadPatternL4RP","Ladder pattern L4 RP",22,0.5,22.5); //position number 8
  h8->GetXaxis()->SetTitle("Ladder #, Layer 4");
  h8->GetYaxis()->SetTitle("Entries");
  fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h8)),8 +fGenRecPointsOffset, expert, !image);
  delete h8;
  fSDDhRecPointsTask++;
  TH2F *h9 = new TH2F("SDDLocalCoordDistrib","Local Coord Distrib",1000/nOnline,-4,4,1000/nOnline,-4,4);//position number 9
  h9->GetXaxis()->SetTitle("X local coord, drift, cm");
  h9->GetYaxis()->SetTitle("Z local coord, anode, cm");
  fAliITSQADataMakerRec->Add2RecPointsList((new TH2F(*h9)),9 +fGenRecPointsOffset, expert, !image);
  delete h9;
  fSDDhRecPointsTask++;


    TH1F *h10 = new TH1F("SDDrdistrib_Layer3" ,"SDD r distribution Layer3" ,100,14.,18.);//position number 10 (L3)
    h10->GetXaxis()->SetTitle("r[cm]");
    h10->GetXaxis()->CenterTitle();
    h10->GetYaxis()->SetTitle("Entries");
    fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h10)),10 +fGenRecPointsOffset, expert, !image);
    delete h10;
    fSDDhRecPointsTask++;

    TH1F *h11 = new TH1F("SDDrdistrib_Layer4" ,"SDD r distribution Layer4" ,100,22.,26.);// and position number 11 (L4)
    h11->GetXaxis()->SetTitle("r[cm]");
    h11->GetXaxis()->CenterTitle();
    h11->GetYaxis()->SetTitle("Entries");
    fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h11)),11 +fGenRecPointsOffset, expert, !image);
    delete h11;
    fSDDhRecPointsTask++;

  for(Int_t iLay=0; iLay<=1; iLay++){
    sprintf(hisnam,"SDDphidistrib_Layer%d",iLay+3);
    TH1F *h12 = new TH1F(hisnam,hisnam,180,-TMath::Pi(),TMath::Pi());//position number 12 (L3) and position number 13 (L4)
    h12->GetXaxis()->SetTitle("#varphi[rad]");
    h12->GetXaxis()->CenterTitle();
    h12->GetYaxis()->SetTitle("Entries");
    fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h12)),iLay+12+fGenRecPointsOffset, expert, !image);
    delete h12;
    fSDDhRecPointsTask++;
  }

  if(fkOnline)
    {
      TH2F *h14 = new TH2F("SDDGlobalCoordDistribYXFSE","YX Global Coord Distrib FSE",5600/nOnline2,-28,28,5600/nOnline2,-28,28);//position number 14
      h14->GetYaxis()->SetTitle("Y[cm]");
      h14->GetXaxis()->SetTitle("X[cm]");
      fAliITSQADataMakerRec->Add2RecPointsList((new TH2F(*h14)),14+fGenRecPointsOffset, expert, !image);
      delete h14;
      fSDDhRecPointsTask++;
      
      TH2F *h15 = new TH2F("SDDGlobalCoordDistribRZFSE","RZ Global Coord Distrib FSE",Int_t(6400/nOnline3),-32,32,1400/nOnline4,12,26);//position number 15
      h15->GetYaxis()->SetTitle("R[cm]");
      h15->GetXaxis()->SetTitle("Z[cm]");
      fAliITSQADataMakerRec->Add2RecPointsList((new TH2F(*h15)),15+fGenRecPointsOffset, expert, !image);
      delete h15;
      fSDDhRecPointsTask++;
      
    }//online

  //printf("%d SDD Recs histograms booked\n",fSDDhRecPointsTask);


  AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SDD Recs histograms booked\n",fSDDhRecPointsTask));


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
	  fAliITSQADataMakerRec->GetRecPointsData(i+fGenRecPointsOffset)->Reset();
	}
    }
  for(Int_t module=0; module<clustersTree->GetEntries();module++){
    branchRecP->GetEvent(module);
    npoints += recpoints->GetEntries();
    AliITSgeomTGeo::GetModuleId(module, lay, lad, det);
    //printf("modnumb %d, lay %d, lad %d, det %d \n",module, lay, lad, det);
    for(Int_t j=0;j<recpoints->GetEntries();j++){
      AliITSRecPoint *recp = (AliITSRecPoint*)recpoints->At(j);    
      fAliITSQADataMakerRec->GetRecPointsData(6 +fGenRecPointsOffset)->Fill(module);//modpatternrp
      recp->GetGlobalXYZ(cluglo);
      Float_t rad=TMath::Sqrt(cluglo[0]*cluglo[0]+cluglo[1]*cluglo[1]); 
      Float_t phi=TMath::ATan2(cluglo[1],cluglo[0]);
      fAliITSQADataMakerRec->GetRecPointsData(9 +fGenRecPointsOffset)->Fill(recp->GetDetLocalX(),recp->GetDetLocalZ());//local distribution
      fAliITSQADataMakerRec->GetRecPointsData(2 +fGenRecPointsOffset)->Fill(cluglo[0],cluglo[1]);//global distribution YX
      fAliITSQADataMakerRec->GetRecPointsData(3 +fGenRecPointsOffset)->Fill(cluglo[2],rad);//global distribution rz
      if(fkOnline)
	{
	  fAliITSQADataMakerRec->GetRecPointsData(14 +fGenRecPointsOffset)->Fill(cluglo[0],cluglo[1]);//global distribution YX FSE
	  fAliITSQADataMakerRec->GetRecPointsData(15 +fGenRecPointsOffset)->Fill(cluglo[2],rad);//global distribution rz FSE
	}
      if(recp->GetLayer() ==2) {
	fAliITSQADataMakerRec->GetRecPointsData(0  +fGenRecPointsOffset)->Fill(recp->GetQ()) ;//total charge of layer 3
	fAliITSQADataMakerRec->GetRecPointsData(7  +fGenRecPointsOffset)->Fill(lad);//lad pattern layer 3
	fAliITSQADataMakerRec->GetRecPointsData(10 +fGenRecPointsOffset)->Fill(rad);//r distribution layer 3
	fAliITSQADataMakerRec->GetRecPointsData(12 +fGenRecPointsOffset)->Fill(phi);// phi distribution layer 3
	fAliITSQADataMakerRec->GetRecPointsData(4  +fGenRecPointsOffset)->Fill(cluglo[2],phi);// phi distribution layer 3
      }
      else if(recp->GetLayer() ==3) {
	fAliITSQADataMakerRec->GetRecPointsData(1  +fGenRecPointsOffset)->Fill(recp->GetQ()) ;//total charge layer 4
	fAliITSQADataMakerRec->GetRecPointsData(8  +fGenRecPointsOffset)->Fill(lad);//ladpatternlayer4
	fAliITSQADataMakerRec->GetRecPointsData(11 +fGenRecPointsOffset)->Fill(rad);//r distribution
	fAliITSQADataMakerRec->GetRecPointsData(13 +fGenRecPointsOffset)->Fill(phi);//phi distribution
	fAliITSQADataMakerRec->GetRecPointsData(5  +fGenRecPointsOffset)->Fill(cluglo[2],phi);// phi distribution layer 4
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


//_______________________________________________________________

Int_t AliITSQASDDDataMakerRec::GetOffset(AliQAv1::TASKINDEX_t task)
{
  Int_t offset=0;
  if( task == AliQAv1::kRAWS )
    {
      offset=fGenRawsOffset;  
    }
  else
    if( task == AliQAv1::kRECPOINTS )
      {
	offset=fGenRecPointsOffset;   
      }
    else
      if(task == AliQAv1::kDIGITSR )
	{
	  offset=fGenDigitsOffset;
	}
    else AliInfo("No task has been selected. Offset set to zero.\n");
  return offset;
}

//_______________________________________________________________

Int_t AliITSQASDDDataMakerRec::GetTaskHisto(AliQAv1::TASKINDEX_t task)
{

  Int_t histotot=0;

  if( task == AliQAv1::kRAWS )
    {
      histotot=fSDDhRawsTask ;  
    }
  else
    if( task == AliQAv1::kRECPOINTS )
      {
	histotot=fSDDhRecPointsTask;   
      }
    else
      if(task == AliQAv1::kDIGITSR)
	{
	  histotot=fSDDhDigitsTask;
	}
    else AliInfo("No task has been selected. TaskHisto set to zero.\n");
  return histotot;
}
