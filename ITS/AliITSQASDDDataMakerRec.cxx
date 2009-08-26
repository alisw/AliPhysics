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
//  M.Siciliano Aug 2008 QA RecPoints 
//  INFN Torino

// --- ROOT system ---

#include <TProfile2D.h>
#include <TH2D.h>
#include <TBranch.h>
#include <TTree.h>
#include <TGaxis.h>
#include <TMath.h>
#include <TF1.h>
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
fGoodAnodes(0),
fBadAnodes(0),
fGoodAnodesCurrent(0),
fBadAnodesCurrent(0)
{
  //ctor used to discriminate OnLine-Offline analysis
  if(fLDC < 0 || fLDC > 6) {
	AliError("Error: LDC number out of range; return\n");
  }
	fGenRawsOffset = new Int_t[AliRecoParam::kNSpecies];
	fGenRecPointsOffset = new Int_t[AliRecoParam::kNSpecies];
	fGenDigitsOffset = new Int_t[AliRecoParam::kNSpecies];
	for(Int_t i=0; i<AliRecoParam::kNSpecies; i++) {
		fGenRawsOffset[i] = 0;
		fGenRecPointsOffset[i] = 0;
		fGenDigitsOffset[i]=0;
	}
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
fGoodAnodes(qadm.fGoodAnodes),
fBadAnodes(qadm.fBadAnodes),
fGoodAnodesCurrent(qadm.fGoodAnodesCurrent),
fBadAnodesCurrent(qadm.fBadAnodesCurrent){
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
	if(fkOnline) {
		AnalyseBNG(); // Analyse Baseline, Noise, Gain
		AnalyseINJ(); // Analyse Injectors
	}
  
	AliDebug(AliQAv1::GetQADebugLevel(),"AliITSDM instantiates checker with Run(AliQAv1::kITS, task, list)\n"); 
}

//____________________________________________________________________________ 
Int_t AliITSQASDDDataMakerRec::InitRaws()
{ 
  // Initialization for RAW data - SDD -
  const Bool_t expert   = kTRUE ; 
  const Bool_t saveCorr = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  Int_t rv = 0 ; 
  AliCDBEntry *ddlMapSDD = AliCDBManager::Instance()->Get("ITS/Calib/DDLMapSDD");
  Bool_t cacheStatus = AliCDBManager::Instance()->GetCacheFlag();
  if(!ddlMapSDD)
    {
      AliError("Calibration object retrieval failed! SDD will not be processed");
      fDDLModuleMap = NULL;
      return rv;
    }
  fDDLModuleMap = (AliITSDDLModuleMapSDD*)ddlMapSDD->GetObject();
  if(!cacheStatus)ddlMapSDD->SetObject(NULL);
  ddlMapSDD->SetOwner(kTRUE);
  if(!cacheStatus)
    {
      delete ddlMapSDD;
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
   rv = fAliITSQADataMakerRec->Add2RawsList((new TH1D(*h0)),0+fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image, !saveCorr);
  delete h0;
  fSDDhRawsTask++;
  
  //zPhi distribution using ladder and modules numbers
  TH2D *hphil3 = new TH2D("SDDphizL3","SDD #varphiz Layer3 ",6,0.5,6.5,14,0.5,14.5);
  hphil3->GetXaxis()->SetTitle("z[#Module L3 ]");
  hphil3->GetYaxis()->SetTitle("#varphi[#Ladder L3]");
   rv = fAliITSQADataMakerRec->Add2RawsList((new TH2D(*hphil3)),1+fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image, saveCorr); 
  delete hphil3;
  fSDDhRawsTask++;
  
  TH2D *hphil4 = new TH2D("SDDphizL4","SDD #varphiz Layer4 ",8,0.5,8.5,22,0.5,22.5); 
  hphil4->GetXaxis()->SetTitle("z[#Module L4]");
  hphil4->GetYaxis()->SetTitle("#varphi[#Ladder L4]");
   rv = fAliITSQADataMakerRec->Add2RawsList((new TH2D(*hphil4)),2+fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image, saveCorr); 
  delete hphil4;
  fSDDhRawsTask++;
  

  if(fkOnline) 
    {

      //DDL Pattern 
      TH2D *hddl = new TH2D("SDDDDLPattern","SDD DDL Pattern ",24,-0.5,23.5,24,-0.5,23.5); 
      hddl->GetXaxis()->SetTitle("Channel");
      hddl->GetYaxis()->SetTitle("#DDL");
      rv = fAliITSQADataMakerRec->Add2RawsList((new TH2D(*hddl)),3+fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image, !saveCorr);
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
	   rv = fAliITSQADataMakerRec->Add2RawsList((new TProfile2D(*fModuleChargeMapFSE)),indexlast1 + index1 + fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image, !saveCorr);
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
	   rv = fAliITSQADataMakerRec->Add2RawsList((new TProfile2D(*fModuleChargeMap)),indexlast1 + index1 + fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image, !saveCorr);
	  delete fModuleChargeMap;
	  
	  fSDDhRawsTask++;
	  index1++;	 
	}
      }
      
	  
    }  // kONLINE
  
  AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SDD Raws histograms booked\n",fSDDhRawsTask));
  return rv ; 
}


//____________________________________________________________________________
Int_t AliITSQASDDDataMakerRec::MakeRaws(AliRawReader* rawReader)
{ 
  // Fill QA for RAW - SDD -
	Int_t rv = 0;
  // Check id histograms already created for this Event Specie

  if(!fDDLModuleMap){
    AliError("SDD DDL module map not available - skipping SDD QA");
    return rv;
  }
  if(rawReader->GetType() != 7) return rv;  // skips non physical triggers
  AliDebug(AliQAv1::GetQADebugLevel(),"entering MakeRaws\n");                 
  rawReader->Reset();       
  AliITSRawStream *stream=AliITSRawStreamSDD::CreateRawStreamSDD(rawReader);
   stream->SetDDLModuleMap(fDDLModuleMap);
  
  Int_t lay, lad, det; 
  
  Int_t index=0;
  if(fkOnline) {
    for(Int_t moduleSDD =0; moduleSDD<fgknSDDmodules; moduleSDD++){
      for(Int_t iside=0;iside<fgknSide;iside++) {
		if(fSDDhRawsTask > 4 + index) fAliITSQADataMakerRec->GetRawsData(4 + index +fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Reset();   
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
    fAliITSQADataMakerRec->GetRawsData( 0 + fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()] )->Fill(isddmod);   
    if(lay==3)    fAliITSQADataMakerRec->GetRawsData(1+fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(det,lad); 
    if(lay==4) { 
      fAliITSQADataMakerRec->GetRawsData(2+fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(det,lad);}  
    
    Short_t iside = stream->GetChannel();
   

    if(fkOnline) {

      fAliITSQADataMakerRec->GetRawsData(3+fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(2*(stream->GetCarlosId())+iside,iddl);

      activeModule = moduleSDD;
      index1 = activeModule * 2 + iside;
      
      if(index1<0){
        AliDebug(AliQAv1::GetQADebugLevel(),Form("Wrong index number %d - patched to 0\n",index1));
	index1 = 0;
      }      

      if(fSDDhRawsTask > 4 + index1) {                                  
        ((TProfile2D *)(fAliITSQADataMakerRec->GetRawsData(4 + index1 +fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()])))->Fill(coord2, coord1, signal);     
        ((TProfile2D *)(fAliITSQADataMakerRec->GetRawsData(4 + index1 + 260*2 +fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()])))->Fill(coord2, coord1, signal); 
      }
    }
    cnt++;
    if(!(cnt%10000)) AliDebug(AliQAv1::GetQADebugLevel(),Form(" %d raw digits read",cnt));
  }
  AliDebug(AliQAv1::GetQADebugLevel(),Form("Event completed, %d raw digits read",cnt)); 
  delete stream;
  stream = NULL; 

//	if(fkOnline) {
//		AnalyseBNG(); // Analyse Baseline, Noise, Gain
//		AnalyseINJ(); // Analyse Injectors
//	}


  return rv ; 
}

//____________________________________________________________________________ 
Int_t AliITSQASDDDataMakerRec::InitDigits()
{ 


  // Initialization for DIGIT data - SDD -  
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ;
  Int_t rv = 0 ; 
//  fGenDigitsOffset = (fAliITSQADataMakerRec->fDigitsQAList[AliRecoParam::kDefault])->GetEntries();
  //fSDDhTask must be incremented by one unit every time a histogram is ADDED to the QA List
  TH1F* h0=new TH1F("SDD DIGITS Module Pattern","SDD DIGITS Module Pattern",260,239.5,499.5);       //hmod
  h0->GetXaxis()->SetTitle("SDD Module Number");
  h0->GetYaxis()->SetTitle("# DIGITS");
  rv = fAliITSQADataMakerRec->Add2DigitsList(h0,fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image);
  fSDDhDigitsTask ++;
  // printf("Add %s \t the task offset is %i\n",fAliITSQADataMakerRec->GetDigitsData(fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->GetName() , fSDDhDigitsTask );
  TH1F* h1=new TH1F("SDD Anode Distribution","DIGITS Anode Distribution",512,-0.5,511.5);      //hanocc
  h1->GetXaxis()->SetTitle("Anode Number");
  h1->GetYaxis()->SetTitle("# DIGITS");
  rv = fAliITSQADataMakerRec->Add2DigitsList(h1,1+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image);
  fSDDhDigitsTask ++;
  //printf("Add %s \t the task offset is %i\n",fAliITSQADataMakerRec->GetDigitsData(1+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->GetName() , fSDDhDigitsTask );
  TH1F* h2=new TH1F("SDD Tbin Distribution","DIGITS Tbin Distribution",256,-0.5,255.5);      //htbocc
  h2->GetXaxis()->SetTitle("Tbin Number");
  h2->GetYaxis()->SetTitle("# DIGITS");
  rv = fAliITSQADataMakerRec->Add2DigitsList(h2,2+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image);
  fSDDhDigitsTask ++;
  //printf("Add %s \t the task offset is %i\n",fAliITSQADataMakerRec->GetDigitsData(2+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->GetName() , fSDDhDigitsTask );
  TH1F* h3=new TH1F("SDD ADC Counts Distribution","DIGITS ADC Counts Distribution",200,0.,1024.);          //hsig
  h3->GetXaxis()->SetTitle("ADC Value");
  h3->GetYaxis()->SetTitle("# DIGITS");
  rv = fAliITSQADataMakerRec->Add2DigitsList(h3,3+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image);
  fSDDhDigitsTask ++;
  //printf("Add %s \t the task offset is %i\n",fAliITSQADataMakerRec->GetDigitsData(3+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->GetName() , fSDDhDigitsTask );
  AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SDD Digits histograms booked\n",fSDDhDigitsTask));
  return rv ; 
}

//____________________________________________________________________________
Int_t AliITSQASDDDataMakerRec::MakeDigits(TTree * digits)
{ 

  // Fill QA for DIGIT - SDD -
  //AliITS *fITS  = (AliITS*)gAlice->GetModule("ITS");
  //fITS->SetTreeAddress();
  //TClonesArray *iITSdigits  = fITS->DigitsAddress(1);


  Int_t rv = 0 ; 

  TBranch *branchD = digits->GetBranch("ITSDigitsSDD");

  if (!branchD) {
    AliError("can't get the branch with the ITS SDD digits !");
    return rv ;
  }
  // Check id histograms already created for this Event Specie
//  if ( ! fAliITSQADataMakerRec->GetDigitsData(fGenDigitsOffset) )
//    rv = InitDigits() ;
  
  static TClonesArray statDigits("AliITSdigitSDD");
  TClonesArray *iITSdigits = &statDigits;
  branchD->SetAddress(&iITSdigits);

  for(Int_t i=0; i<260; i++){
    Int_t nmod=i+240;
    digits->GetEvent(nmod);
    Int_t ndigits = iITSdigits->GetEntries();
    fAliITSQADataMakerRec->GetDigitsData(fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(nmod,ndigits);

    for (Int_t idig=0; idig<ndigits; idig++) {
      AliITSdigit *dig=(AliITSdigit*)iITSdigits->UncheckedAt(idig);
      Int_t iz=dig->GetCoord1();  // cell number z
      Int_t ix=dig->GetCoord2();  // cell number x
      Int_t sig=dig->GetSignal();
      fAliITSQADataMakerRec->GetDigitsData(1+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(iz);
      fAliITSQADataMakerRec->GetDigitsData(2+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(ix);
      fAliITSQADataMakerRec->GetDigitsData(3+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(sig);
    }
  }
  return rv ; 
}

//____________________________________________________________________________ 
Int_t AliITSQASDDDataMakerRec::InitRecPoints()
{

	//AliInfo("Initialize SDD recpoints histos\n");
  // Initialization for RECPOINTS - SDD -
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  Int_t rv = 0 ; 
//  fGenRecPointsOffset = (fAliITSQADataMakerRec->fRecPointsQAList[AliRecoParam::kDefault])->GetEntries();

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

  //AliInfo(Form("fAliITSQADataMakerRec->GetEventSpecie() %d\n",fAliITSQADataMakerRec->GetEventSpecie()));
  //AliInfo(Form("fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] %d\n",fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]));
  TH1F *h0 = new TH1F("SDDLay3TotCh","Layer 3 total charge",1000/nOnline,-0.5, 499.5); //position number 0
  h0->GetXaxis()->SetTitle("ADC value");
  h0->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h0)), 0 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image);
  //delete h0;
  fSDDhRecPointsTask++;
 
  TH1F *h1 = new TH1F("SDDLay4TotCh","Layer 4 total charge",1000/nOnline,-0.5, 499.5);//position number 1
  h1->GetXaxis()->SetTitle("ADC value");
  h1->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h1)), 1 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image);
  delete h1;
  fSDDhRecPointsTask++;

  char hisnam[50];
  TH2F *h2 = new TH2F("SDDGlobalCoordDistribYX","YX Global Coord Distrib",5600/nOnline2,-28,28,5600/nOnline2,-28,28);//position number 2
  h2->GetYaxis()->SetTitle("Y[cm]");
  h2->GetXaxis()->SetTitle("X[cm]");
  rv = fAliITSQADataMakerRec->Add2RecPointsList((new TH2F(*h2)),2+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, !image);
  delete h2;
  fSDDhRecPointsTask++;

  TH2F *h3 = new TH2F("SDDGlobalCoordDistribRZ","RZ Global Coord Distrib",6400/nOnline3,-32,32,1400/nOnline4,12,26);//position number 3
  h3->GetYaxis()->SetTitle("R[cm]");
  h3->GetXaxis()->SetTitle("Z[cm]");
  rv = fAliITSQADataMakerRec->Add2RecPointsList((new TH2F(*h3)),3+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, !image);
  delete h3;
  fSDDhRecPointsTask++;
  
  TH2F *h4 = new TH2F("SDDGlobalCoordDistribL3PHIZ","#varphi Z Global Coord Distrib L3",6400/nOnline3,-32,32,360/nOnline,-TMath::Pi(),TMath::Pi());//position number 4
  h4->GetYaxis()->SetTitle("#phi[rad]");
  h4->GetXaxis()->SetTitle("Z[cm]");
  rv = fAliITSQADataMakerRec->Add2RecPointsList((new TH2F(*h4)),4+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image);
  delete h4;
  fSDDhRecPointsTask++;

  TH2F *h5 = new TH2F("SDDGlobalCoordDistribL4PHIZ","#varphi Z Global Coord Distrib L4",6400/nOnline3,-32,32,360/nOnline,-TMath::Pi(),TMath::Pi());//position number 5
  h5->GetYaxis()->SetTitle("#phi[rad]");
  h5->GetXaxis()->SetTitle("Z[cm]");
  rv = fAliITSQADataMakerRec->Add2RecPointsList((new TH2F(*h5)),5+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image);
  delete h5;
  fSDDhRecPointsTask++;
  
  TH1F *h6 = new TH1F("SDDModPatternRP","Modules pattern RP",fgknSDDmodules,239.5,499.5); //position number 6
  h6->GetXaxis()->SetTitle("Module number"); //spd offset = 240
  h6->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h6)),6 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);
  delete h6;
  fSDDhRecPointsTask++;
  TH1F *h7 = new TH1F("SDDLadPatternL3RP","Ladder pattern L3 RP",14,0.5,14.5);  //position number 7
  h7->GetXaxis()->SetTitle("Ladder #, Layer 3");
  h7->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h7)),7 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);
  delete h7;
  fSDDhRecPointsTask++;
  TH1F *h8 = new TH1F("SDDLadPatternL4RP","Ladder pattern L4 RP",22,0.5,22.5); //position number 8
  h8->GetXaxis()->SetTitle("Ladder #, Layer 4");
  h8->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h8)),8 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);
  delete h8;
  fSDDhRecPointsTask++;
  TH2F *h9 = new TH2F("SDDLocalCoordDistrib","Local Coord Distrib",1000/nOnline,-4,4,1000/nOnline,-4,4);//position number 9
  h9->GetXaxis()->SetTitle("X local coord, drift, cm");
  h9->GetYaxis()->SetTitle("Z local coord, anode, cm");
  rv = fAliITSQADataMakerRec->Add2RecPointsList((new TH2F(*h9)),9 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);
  delete h9;
  fSDDhRecPointsTask++;

	//AliInfo("Create SDD recpoints histos\n");

    TH1F *h10 = new TH1F("SDDrdistrib_Layer3" ,"SDD r distribution Layer3" ,100,14.,18.);//position number 10 (L3)
    h10->GetXaxis()->SetTitle("r[cm]");
    h10->GetXaxis()->CenterTitle();
    h10->GetYaxis()->SetTitle("Entries");
    rv = fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h10)),10 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);
    delete h10;
    fSDDhRecPointsTask++;

    TH1F *h11 = new TH1F("SDDrdistrib_Layer4" ,"SDD r distribution Layer4" ,100,22.,26.);// and position number 11 (L4)
    h11->GetXaxis()->SetTitle("r[cm]");
    h11->GetXaxis()->CenterTitle();
    h11->GetYaxis()->SetTitle("Entries");
    rv = fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h11)),11 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);
    delete h11;
    fSDDhRecPointsTask++;

  for(Int_t iLay=0; iLay<=1; iLay++){
    sprintf(hisnam,"SDDphidistrib_Layer%d",iLay+3);
    TH1F *h12 = new TH1F(hisnam,hisnam,180,-TMath::Pi(),TMath::Pi());//position number 12 (L3) and position number 13 (L4)
    h12->GetXaxis()->SetTitle("#varphi[rad]");
    h12->GetXaxis()->CenterTitle();
    h12->GetYaxis()->SetTitle("Entries");
    rv = fAliITSQADataMakerRec->Add2RecPointsList((new TH1F(*h12)),iLay+12+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);
    delete h12;
    fSDDhRecPointsTask++;
  }

  if(fkOnline)
    {
      TH2F *h14 = new TH2F("SDDGlobalCoordDistribYXFSE","YX Global Coord Distrib FSE",5600/nOnline2,-28,28,5600/nOnline2,-28,28);//position number 14
      h14->GetYaxis()->SetTitle("Y[cm]");
      h14->GetXaxis()->SetTitle("X[cm]");
      rv = fAliITSQADataMakerRec->Add2RecPointsList((new TH2F(*h14)),14+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);
      delete h14;
      fSDDhRecPointsTask++;
      
      TH2F *h15 = new TH2F("SDDGlobalCoordDistribRZFSE","RZ Global Coord Distrib FSE",Int_t(6400/nOnline3),-32,32,1400/nOnline4,12,26);//position number 15
      h15->GetYaxis()->SetTitle("R[cm]");
      h15->GetXaxis()->SetTitle("Z[cm]");
      rv = fAliITSQADataMakerRec->Add2RecPointsList((new TH2F(*h15)),15+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);
      delete h15;
      fSDDhRecPointsTask++;
      
    }//online

  AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SDD Recs histograms booked\n",fSDDhRecPointsTask));

  return rv ; 
}

//____________________________________________________________________________ 
Int_t AliITSQASDDDataMakerRec::MakeRecPoints(TTree * clustersTree)
{
 // Fill QA for RecPoints - SDD -
  Int_t rv = 0 ; 

  //AliInfo(Form("fAliITSQADataMakerRec->GetEventSpecie() %d\n",fAliITSQADataMakerRec->GetEventSpecie()));
  //AliInfo(Form("fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] %d\n",fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]));
  // Check id histograms already created for this Event Specie
//  if ( ! fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]) )
//    rv = InitRecPoints() ;
  Int_t lay, lad, det; 
  //AliInfo("get the branch with the ITS clusters !\n");
  TBranch *branchRecP = clustersTree->GetBranch("ITSRecPoints");
  if (!branchRecP) {
    AliError("can't get the branch with the ITS clusters !");
    return rv ;
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
	  fAliITSQADataMakerRec->GetRecPointsData(i+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Reset();
	}
    }
	for(Int_t module=0; module<clustersTree->GetEntries();module++){
		//AliInfo(Form("Module %d\n",module));
		branchRecP->GetEvent(module);
		npoints += recpoints->GetEntries();
		//AliInfo(Form("modnumb %d, npoints %d, total points %d\n",module, recpoints->GetEntries(),npoints));
		AliITSgeomTGeo::GetModuleId(module, lay, lad, det);
		//AliInfo(Form("modnumb %d, lay %d, lad %d, det %d \n",module, lay, lad, det));
		//Bool_t kSDD = kFALSE;
		//if(lay == 3 || lay == 4) kSDD = kTRUE;
		//if(!kSDD) continue;
		//AliInfo(Form("modnumb %d, entries %d\n",module, recpoints->GetEntries()));
		for(Int_t j=0;j<recpoints->GetEntries();j++){
			//AliInfo(Form("modnumb %d, entry %d \n",module, j));
			AliITSRecPoint *recp = (AliITSRecPoint*)recpoints->At(j);    
			fAliITSQADataMakerRec->GetRecPointsData(6 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(module);//modpatternrp
			recp->GetGlobalXYZ(cluglo);
			Float_t rad=TMath::Sqrt(cluglo[0]*cluglo[0]+cluglo[1]*cluglo[1]); 
			Float_t phi=TMath::ATan2(cluglo[1],cluglo[0]);
			fAliITSQADataMakerRec->GetRecPointsData(9 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(recp->GetDetLocalX(),recp->GetDetLocalZ());//local distribution
			fAliITSQADataMakerRec->GetRecPointsData(2 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(cluglo[0],cluglo[1]);//global distribution YX
			fAliITSQADataMakerRec->GetRecPointsData(3 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(cluglo[2],rad);//global distribution rz
			if(fkOnline) {
				fAliITSQADataMakerRec->GetRecPointsData(14 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(cluglo[0],cluglo[1]);//global distribution YX FSE
				fAliITSQADataMakerRec->GetRecPointsData(15 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(cluglo[2],rad);//global distribution rz FSE
			}
			if(recp->GetLayer() == 2) {
				fAliITSQADataMakerRec->GetRecPointsData(0  +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(recp->GetQ()) ;//total charge of layer 3
				fAliITSQADataMakerRec->GetRecPointsData(7  +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(lad);//lad pattern layer 3
				fAliITSQADataMakerRec->GetRecPointsData(10 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(rad);//r distribution layer 3
				fAliITSQADataMakerRec->GetRecPointsData(12 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(phi);// phi distribution layer 3
				fAliITSQADataMakerRec->GetRecPointsData(4  +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(cluglo[2],phi);// phi distribution layer 3
			} else if(recp->GetLayer() == 3) {
				fAliITSQADataMakerRec->GetRecPointsData(1  +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(recp->GetQ()) ;//total charge layer 4
				fAliITSQADataMakerRec->GetRecPointsData(8  +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(lad);//ladpatternlayer4
				fAliITSQADataMakerRec->GetRecPointsData(11 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(rad);//r distribution
				fAliITSQADataMakerRec->GetRecPointsData(13 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(phi);//phi distribution
				fAliITSQADataMakerRec->GetRecPointsData(5  +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(cluglo[2],phi);// phi distribution layer 4
			}
		}
	}
	statRecpoints.Clear();
	return rv ; 
}

//_______________________________________________________________

Int_t AliITSQASDDDataMakerRec::GetOffset(AliQAv1::TASKINDEX_t task)
{
  Int_t offset=0;
  if( task == AliQAv1::kRAWS )
    {
      offset=fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()];  
    }
  else if(task == AliQAv1::kDIGITSR )
    {
      offset=fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()];
    }
  else if( task == AliQAv1::kRECPOINTS )
    {
      offset=fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()];   
    }
  return offset;
}

//_______________________________________________________________

void AliITSQASDDDataMakerRec::SetOffset(AliQAv1::TASKINDEX_t task, Int_t offset, Int_t specie) {
  // Returns offset number according to the specified task
  if( task == AliQAv1::kRAWS ) {
    fGenRawsOffset[specie]=offset;
  }
  else if( task == AliQAv1::kDIGITSR ) {
    fGenDigitsOffset[specie]=offset;
  }
  else if( task == AliQAv1::kRECPOINTS ) {
    fGenRecPointsOffset[specie]=offset;
  }
}

//_______________________________________________________________

Int_t AliITSQASDDDataMakerRec::GetTaskHisto(AliQAv1::TASKINDEX_t task)
{

  Int_t histotot=0;

  if( task == AliQAv1::kRAWS )
    {
      histotot=fSDDhRawsTask ;  
    }
  else if(task == AliQAv1::kDIGITSR)
    {
      histotot=fSDDhDigitsTask;
    }
  else if( task == AliQAv1::kRECPOINTS )
    {
      histotot=fSDDhRecPointsTask;   
    }
  else {
    AliInfo("No task has been selected. TaskHisto set to zero.\n");
  }
  return histotot;
}

//_______________________________________________________________

void AliITSQASDDDataMakerRec::AnalyseBNG()
{

// get file time for Previous test
	AliInfo("AnalyseBNG\n");
	Int_t bngtimeBasPrevious; 
	FILE *fpinPreviousBas = fopen( "SDDbasHistos.time", "r" );
	if(fpinPreviousBas) {
	  fscanf(fpinPreviousBas,"%d",&bngtimeBasPrevious);
	  fclose (fpinPreviousBas);
	} else 
          bngtimeBasPrevious = 0;
        Int_t bngtimeBasCurrent = 0; 
        gSystem->Exec("perl -e '@d=localtime ((stat(shift))[9]); printf \"%02d%02d%02d%02d%02d\n\", $d[5]-100,$d[4]+1,$d[3],$d[2],$d[1]'  SDDbasHistos.root > SDDbasHistos.time");
        FILE *fpinBas = fopen( "SDDbasHistos.time", "r" );
        fscanf(fpinBas,"%d",&bngtimeBasCurrent);
     	if(bngtimeBasCurrent>bngtimeBasPrevious )AliInfo(Form("bngtimeBasCurrent %d, bngtimeBasPrevious %d\n",bngtimeBasCurrent,bngtimeBasPrevious));
       
	Bool_t kAnalyseBas = kTRUE;
	if(bngtimeBasCurrent <= bngtimeBasPrevious) kAnalyseBas = kFALSE;
	if(kAnalyseBas) {
	        // new bng file found
		bngtimeBasPrevious = bngtimeBasCurrent;
		Bool_t kFilesExist = kTRUE;
		TFile basFile("SDDbasHistos.root");
		if(basFile.IsZombie()) kFilesExist = kFALSE;
		if(kFilesExist) {
		  AnodeStatus();
		  AnalyseHistos(1); // Baseline
		  AnalyseHistos(2); // Uncorrected Noise
		  AnalyseHistos(3); // Common Mode Noise
		  AnalyseHistos(4); // Corrected Noise
		  gSystem->Exec("cp SDDbasHistos.root SDDbasHistosPrevious.root");
		} else {
			AliInfo("file SDDbasHistos.root not found \n");
		}
	}
	fclose (fpinBas);
//	delete fpinBas;

	Int_t bngtimeGainPrevious; 
	FILE *fpinPrevious = fopen( "SDDgainHistos.time", "r" );
	if(fpinPrevious) {
	  fscanf(fpinPrevious,"%d",&bngtimeGainPrevious);
	  fclose (fpinPrevious);
	} else 
	  bngtimeGainPrevious = 0;
        Int_t bngtimeGainCurrent = 0; 
        gSystem->Exec("perl -e '@d=localtime ((stat(shift))[9]); printf \"%02d%02d%02d%02d%02d\n\", $d[5]-100,$d[4]+1,$d[3],$d[2],$d[1]'  SDDgainHistos.root > SDDgainHistos.time");
        FILE *fpin = fopen( "SDDgainHistos.time", "r" );
        fscanf(fpin,"%d",&bngtimeGainCurrent);
        if(bngtimeGainCurrent>bngtimeGainPrevious )AliInfo(Form("bngtimeGainCurrent %d, bngtimeGainPrevious %d\n",bngtimeGainCurrent,bngtimeGainPrevious));
       
	Bool_t kAnalyse = kTRUE;
	if(bngtimeGainCurrent <= bngtimeGainPrevious) kAnalyse = kFALSE;
	if(kAnalyse) {
	        // new bng file found
		bngtimeGainPrevious = bngtimeGainCurrent;
		Bool_t kFilesExist = kTRUE;
		TFile gainFile("SDDgainHistos.root");
		if(gainFile.IsZombie()) kFilesExist = kFALSE;
		if(kFilesExist) {
		  AnalyseHistos(5); // Gain
		  gSystem->Exec("cp SDDgainHistos.root SDDgainHistosPrevious.root");
		} else {
			AliInfo("file SDDgainHistos.root not found \n");
		}
	}
	fclose (fpin);
//	delete fpin;

}

//_______________________________________________________________

void AliITSQASDDDataMakerRec::AnalyseINJ()
{
// get file time for last test

	AliInfo("AnalyseINJ\n");
	Int_t injtimePrevious; 
	FILE *fpinPrevious = fopen( "SDDinjectHistos.time", "r" );
	if(fpinPrevious) {
          fscanf(fpinPrevious,"%d",&injtimePrevious);
	  fclose (fpinPrevious);
	} else 
	  injtimePrevious = 0;
	Int_t injtimeCurrent = 0; 
        gSystem->Exec("perl -e '@d=localtime ((stat(shift))[9]); printf \"%02d%02d%02d%02d%02d\n\", $d[5]-100,$d[4]+1,$d[3],$d[2],$d[1]'  SDDinjectHistos.root > SDDinjectHistos.time");
        FILE *fpin = fopen( "SDDinjectHistos.time", "r" );
	fscanf(fpin,"%d",&injtimeCurrent);
       	if(injtimeCurrent>injtimePrevious )AliInfo(Form("injtimeCurrent %d, injtimePrevious %d\n",injtimeCurrent,injtimePrevious));
       
	Bool_t kAnalyse = kTRUE;
	if(injtimeCurrent <= injtimePrevious) kAnalyse = kFALSE;
	if(kAnalyse) {
		// new inj file found
		injtimePrevious = injtimeCurrent;
		Bool_t kFilesExist = kTRUE;
      		TFile gainFile("SDDinjectHistos.root");
       		if(gainFile.IsZombie()) kFilesExist = kFALSE;
		
		if(kFilesExist) {
			AnalyseHistos(6); // Drift Speed
			gSystem->Exec("cp SDDinjectHistos.root SDDinjectHistosPrevious.root");
		} else {
			AliInfo("file(s) SDDinjectHistos.root not found \n");
		}
	}
	fclose (fpin);
//	delete fpin;
}

//_______________________________________________________________

void AliITSQASDDDataMakerRec::AnodeStatus()
{
	char *hnamePrevious = new char[50];
	fGoodAnodes = 0;

       	TFile basFilePrevious("SDDbasHistosPrevious.root");
       	if(!basFilePrevious.IsZombie()) {
	  for(Int_t ddl =0; ddl<fDDLModuleMap->GetNDDLs(); ddl++){
       		for(Int_t crx =0; crx<fDDLModuleMap->GetNModPerDDL(); crx++){
       			for(Int_t iside=0;iside<fgknSide;iside++){
       				Int_t moduleSDD = fDDLModuleMap->GetModuleNumber(ddl,crx);
       		//AliITSgeomTGeo::GetModuleId(moduleSDD+fgkmodoffset, lay, lad, det);
       				sprintf(hnamePrevious,"hgood%02dc%02ds%d",ddl,crx,iside);
       			//AliInfo(Form("get histo %s\n",hnamePrevious));
       				TH1F *hgood = (TH1F *) basFilePrevious.Get(hnamePrevious);
       				if(!hgood) continue;
       				for(Int_t i=0; i<hgood->GetNbinsX();i++) {
       					fAnodeMap[moduleSDD-fgkmodoffset][iside][i] = hgood->GetBinContent(i);
       					if(fAnodeMap[moduleSDD-fgkmodoffset][iside][i]) fGoodAnodes++;
       				}
       				delete hgood;
       			}
       		}
	  }
	  basFilePrevious.Close();
	}
	fGoodAnodesCurrent = 0;
	fBadAnodesCurrent = 0;
	char *hname = new char[50];
	Int_t nChangedStatus = 0;
	Bool_t CurrentAnodeMap[fgknSDDmodules][fgknSide][fgknAnode];	
      	TFile basFile("SDDbasHistos.root");
	if(!basFile.IsZombie()) {
	  for(Int_t ddl =0; ddl<fDDLModuleMap->GetNDDLs(); ddl++){
       		for(Int_t crx =0; crx<fDDLModuleMap->GetNModPerDDL(); crx++){
       			for(Int_t iside=0;iside<fgknSide;iside++){
       				Int_t moduleSDD = fDDLModuleMap->GetModuleNumber(ddl,crx);
       			//AliITSgeomTGeo::GetModuleId(moduleSDD+fgkmodoffset, lay, lad, det);
       				sprintf(hname,"hgood%02dc%02ds%d",ddl,crx,iside);
       			//AliInfo(Form("get histo %s\n",hname));
       				TH1F *hgood = (TH1F *) basFile.Get(hname);
       				if(!hgood) continue;
       				for(Int_t i=0; i<hgood->GetNbinsX();i++) {
       					CurrentAnodeMap[moduleSDD-fgkmodoffset][iside][i] = hgood->GetBinContent(i);
       					if(CurrentAnodeMap[moduleSDD-fgkmodoffset][iside][i]) fGoodAnodesCurrent++;
       					else fBadAnodesCurrent++;
       					if(CurrentAnodeMap[moduleSDD-fgkmodoffset][iside][i] != fAnodeMap[moduleSDD-fgkmodoffset][iside][i]) {
       						fAnodeMap[moduleSDD-fgkmodoffset][iside][i] = CurrentAnodeMap[moduleSDD-fgkmodoffset][iside][i];
       						nChangedStatus++;
						//	AliWarning(Form("DDL %d, CRX %d, Side %d, Anode %d changed status to %d \n",ddl,crx,iside,i,fAnodeMap[moduleSDD-fgkmodoffset][iside][i]));
       					}
       				}
       				delete hgood;
       			}
       		}
	  }
	  basFile.Close();
	}

	AliWarning(Form("Number of good anodes changed from %d to %d, that is %f %%\n",fGoodAnodes,fGoodAnodesCurrent,((Float_t) TMath::Abs(fGoodAnodes-fGoodAnodesCurrent))/(fBadAnodesCurrent+fGoodAnodesCurrent)));
	if(fGoodAnodesCurrent != fGoodAnodes) {
		fGoodAnodes = fGoodAnodesCurrent;
	}
	AliWarning(Form("Number of bad anodes changed from %d to %d, that is %f %%\n",fBadAnodes,fBadAnodesCurrent,((Float_t) TMath::Abs(fBadAnodes-fBadAnodesCurrent))/(fBadAnodesCurrent+fGoodAnodesCurrent)));
	if(fBadAnodesCurrent != fBadAnodes) {
		fBadAnodes = fBadAnodesCurrent;
	}
//	delete hname;
}

//_______________________________________________________________

void AliITSQASDDDataMakerRec::AnalyseHistos(Int_t type)
{

	if(type < 1 || type > 6) {
	  AliWarning(Form("Wrong type (%d), must be between 1 and 6\n",type));
	  return;
	}

	Double_t Current[fgknSDDmodules][fgknSide][fgknAnode];	
	char *hnamePrevious = new char[50];
	TString *fnamePrevious=NULL;

		if(type < 5) fnamePrevious = new TString("SDDbasHistosPrevious.root");
		else if(type == 5) fnamePrevious = new TString("SDDgainHistosPrevious.root");
		else if(type == 6) fnamePrevious = new TString("SDDinjectHistosPrevious.root");
		TFile *gainFilePrevious = new TFile(fnamePrevious->Data());
		if(!gainFilePrevious->IsZombie()) {

		  for(Int_t ddl =0; ddl<fDDLModuleMap->GetNDDLs(); ddl++){
		    for(Int_t crx =0; crx<fDDLModuleMap->GetNModPerDDL(); crx++){
				for(Int_t iside=0;iside<fgknSide;iside++){
					Int_t moduleSDD = fDDLModuleMap->GetModuleNumber(ddl,crx);
			//AliITSgeomTGeo::GetModuleId(moduleSDD+fgkmodoffset, lay, lad, det);
					if(type == 1) sprintf(hnamePrevious,"hbase%02dc%02ds%d",ddl,crx,iside);
					else if(type == 2) sprintf(hnamePrevious,"hnois%02dc%02ds%d",ddl,crx,iside);
					else if(type == 3) sprintf(hnamePrevious,"hcmn%02dc%02ds%d",ddl,crx,iside);
					else if(type == 4) sprintf(hnamePrevious,"hcorn%02dc%02ds%d",ddl,crx,iside);
					else if(type == 5) sprintf(hnamePrevious,"hgain%02dc%02ds%d",ddl,crx,iside);
					else if(type == 6) sprintf(hnamePrevious,"hdrsp%02dc%02ds%d",ddl,crx,iside);
				//AliInfo(Form("get histo %s\n",hnamePrevious));
					TH1F *hhist = (TH1F *) gainFilePrevious->Get(hnamePrevious);
					if(!hhist) continue;
					for(Int_t i=0; i<hhist->GetNbinsX();i++) {
						Current[moduleSDD-fgkmodoffset][iside][i] = hhist->GetBinContent(i);
					}
					delete hhist;
				}
			}
		  }
		  gainFilePrevious->Close();
		  delete gainFilePrevious;
		}
		delete fnamePrevious;

	Float_t xmin = 0.;
	Float_t xmax = 0;
	Int_t nbins = 1;
	TH1F *hDist = 0;
	TH1F *hDistDiff = 0;
	if(type == 1) {
		xmin = 0.;
		xmax = 500.;
		nbins = (Int_t)(xmax-xmin);
 		hDist = new TH1F("hBaseline","Baseline",nbins,xmin,xmax);
		hDistDiff = new TH1F("hBaselineDiff","Baseline Difference",200,-100.,100.);
	} else if(type == 2) {
		xmin = 0.;
		xmax = 10.;
		nbins = (Int_t) (100.*(xmax-xmin));
		hDist = new TH1F("hNoiseUnc","Noise (before correction)",nbins,xmin,xmax);
		hDistDiff = new TH1F("hNoiseUncDiff","Noise (before correction) Difference",200,-10.,10.);
	} else if(type == 3) {
		xmin = 0.;
		xmax = 10.;
		nbins = (Int_t)( 100.*(xmax-xmin));
		hDist = new TH1F("hNoiseCMN","Noise (common mode)",nbins,xmin,xmax);
		hDistDiff = new TH1F("hNoiseCMNDiff","Noise (common mode) Difference",200,-10.,10.);
	} else if(type == 4) {
		xmin = 0.;
		xmax = 10.;
		nbins = (Int_t)(100.*(xmax-xmin));
		hDist = new TH1F("hNoiseCor","Noise (after correction)",nbins,xmin,xmax);
		hDistDiff = new TH1F("hNoiseCorDiff","Noise (after correction) Difference",200,-10.,10.);
	} else if(type == 5) {
		xmin = 0.;
		xmax = 5.;
		nbins = (Int_t)(100.*(xmax-xmin));
		hDist = new TH1F("hGain","Gain",nbins,xmin,xmax);
		hDistDiff = new TH1F("hGainDiff","Gain Difference",200,-10.,10.);
	} else if(type == 6) {
		xmin = 0.;
		xmax = 10.;
		nbins = (Int_t)(100.*(xmax-xmin));
		hDist = new TH1F("hDriftSpeed","Drift Speed",nbins,xmin,xmax);
		hDistDiff = new TH1F("hDriftSpeedDiff","Drift Speed Difference",200,-10.,10.);
	}

	Float_t	binw = (xmax-xmin)/nbins;

	TString *fnamePrevious2=NULL;

		if(type < 5) fnamePrevious2 = new TString("SDDbasHistosPrevious.root");
		else if(type == 5) fnamePrevious2 = new TString("SDDgainHistosPrevious.root");
		else if(type == 6) fnamePrevious2 = new TString("SDDinjectHistosPrevious.root");
		TFile *gainFile = new TFile(fnamePrevious2->Data());
		if(!gainFile->IsZombie()) {

		  for(Int_t ddl =0; ddl<fDDLModuleMap->GetNDDLs(); ddl++){
			for(Int_t crx =0; crx<fDDLModuleMap->GetNModPerDDL(); crx++){
				for(Int_t iside=0;iside<fgknSide;iside++){
					Int_t moduleSDD = fDDLModuleMap->GetModuleNumber(ddl,crx);
			//AliITSgeomTGeo::GetModuleId(moduleSDD+fgkmodoffset, lay, lad, det);
					if(type == 1) sprintf(hnamePrevious,"hbase%02dc%02ds%d",ddl,crx,iside);
					else if(type == 2) sprintf(hnamePrevious,"hnois%02dc%02ds%d",ddl,crx,iside);
					else if(type == 3) sprintf(hnamePrevious,"hcmn%02dc%02ds%d",ddl,crx,iside);
					else if(type == 4) sprintf(hnamePrevious,"hcorn%02dc%02ds%d",ddl,crx,iside);
					else if(type == 5) sprintf(hnamePrevious,"hgain%02dc%02ds%d",ddl,crx,iside);
					else if(type == 6) sprintf(hnamePrevious,"hdrsp%02dc%02ds%d",ddl,crx,iside);
				//AliInfo(Form("get histo %s\n",hname));
					TH1F *hhist = (TH1F *) gainFile->Get(hnamePrevious);
					if(!hhist) continue;
					for(Int_t i=0; i<hhist->GetNbinsX();i++) {
						if(!fAnodeMap[moduleSDD-fgkmodoffset][iside][i]) continue;
						hDist->Fill(hhist->GetBinContent(i));
						hDistDiff->Fill(hhist->GetBinContent(i)-Current[moduleSDD-fgkmodoffset][iside][i]);
					}
					delete hhist;
				}
			}
		  }
		  gainFile->Close();
		  delete gainFile;
		}
		delete fnamePrevious2;

	TF1 ff("ff", "gaus", xmin+0.1, xmax-0.1);
	hDist->Fit("ff","NWR");
//	hDist->Fit("gaus","","",xmin+0.1, xmax-0.1);
//	Float_t ChiSquared = (Float_t) ff.GetChisquare();
//	Int_t NDF = ff.GetNumberFitPoints() - ff.GetNpar();
	Float_t average = (Float_t) ff.GetParameter(1);
	Float_t sigma = (Float_t) ff.GetParameter(2);
//	Float_t mean = hDist->GetMean();
//	Float_t rms = hDist->GetRMS();
	Int_t badB = 0;
	for(Int_t i=0; i<hDist->GetNbinsX();i++) {
//		if(type < 6) 
	  if(TMath::Abs(i*binw-average) > 4.*sigma) badB += (Int_t)hDist->GetBinContent(i);
//		else
//			if(TMath::Abs(i-mean) > 4*rms) badB += hDist->GetBinContent(i);
	}
        Double_t denomi = hDist->GetEntries();
        if(denomi == 0) {
	  denomi = 1;
	  badB = 0; 
        } 
	if(type == 1) {
		AliInfo(Form("Number of anodes with baseline out of 4*sigma from average: %d, %f%%\n",badB,100.*((Float_t) badB)/denomi));
	} else if(type == 2) {
		AliInfo(Form("Number of anodes with uncorrected noise out of 4*sigma from average: %d, %f%%\n",badB,100.*((Float_t) badB)/denomi));
	} else if(type == 3) {
		AliInfo(Form("Number of anodes with common mode noise out of 4*sigma from average: %d, %f%%\n",badB,100.*((Float_t) badB)/denomi));
	} else if(type == 4) {
		AliInfo(Form("Number of anodes with corrected noise out of 4*sigma from average: %d, %f%%\n",badB,100.*((Float_t) badB)/denomi));
	} else if(type == 5) {
		AliInfo(Form("Number of anodes with gain out of 4*sigma from average: %d, %f%%\n",badB,100.*((Float_t) badB)/denomi));
	} else if(type == 6) {
		Int_t badspeed = (Int_t)hDist->GetBinContent(1);
		AliInfo(Form("Number of anodes with drift speed equal to 0: %d\n",badspeed));
		AliInfo(Form("Number of anodes with drift speed out of 4*sigma from average: %d, %f%%\n",badB-badspeed,100.*((Float_t) (badB-badspeed))/(denomi-badspeed)));
	}
	
	TH1F *hDistHistoryCurrent = NULL;
	TH1F *hDistHistoryPrevious = NULL;

	TFile *gainHistoryFile=NULL;
	if(type < 5) 
		gainHistoryFile = new TFile("SDDbasHistory.root","UPDATE");
	else if(type ==5) 
		gainHistoryFile = new TFile("SDDgainHistory.root","UPDATE");
	else if(type == 6)
		gainHistoryFile = new TFile("SDDinjectHistory.root","UPDATE");
	hDist->Write();
	hDistDiff->Write();
	//AliInfo("SDDgainHistory.root file opened\n");
	if(!gainHistoryFile->IsZombie()) {
		if(type == 1) hDistHistoryPrevious = (TH1F *) gainHistoryFile->Get("hBaselineHistory");
		else if(type == 2) hDistHistoryPrevious = (TH1F *) gainHistoryFile->Get("hNoiseUncHistory");
		else if(type == 3) hDistHistoryPrevious = (TH1F *) gainHistoryFile->Get("hNoiseCMNHistory");
		else if(type == 4) hDistHistoryPrevious = (TH1F *) gainHistoryFile->Get("hNoiseCorHistory");
		else if(type == 5) hDistHistoryPrevious = (TH1F *) gainHistoryFile->Get("hGainHistory");
		else if(type == 6) hDistHistoryPrevious = (TH1F *) gainHistoryFile->Get("hDriftSpeedHistory");
		//AliInfo(Form("hDistHistoryPrevious %x\n",hDistHistoryPrevious));
	
		if(!hDistHistoryPrevious) {
			if(type == 1) hDistHistoryCurrent = new TH1F("hBaselineHistory","Average Baseline History",2,0,2);
			else if(type == 2) hDistHistoryCurrent = new TH1F("hNoiseUncHistory","Average Uncorrected Noise History",2,0,2);
			else if(type == 3) hDistHistoryCurrent = new TH1F("hNoiseCMNHistory","Average Common Mode Noise History",2,0,2);
			else if(type == 4) hDistHistoryCurrent = new TH1F("hNoiseCorHistory","Average Corrected Noise History",2,0,2);
			else if(type == 5) hDistHistoryCurrent = new TH1F("hGainHistory","Average Gain History",2,0,2);
			else if(type == 6) hDistHistoryCurrent = new TH1F("hDriftSpeedHistory","Average Drift Speed History",2,0,2);
			//AliInfo(Form("hDistHistoryCurrent 1 %x\n",hDistHistoryCurrent));
//			if(type < 6) {
				hDistHistoryCurrent->SetBinContent(1,average);
				hDistHistoryCurrent->SetBinError(1,sigma);
/*
			} else {
				hDistHistoryCurrent->SetBinContent(1,mean);
				hDistHistoryCurrent->SetBinError(1,rms);
			}
*/
		} else {
			if(type == 1) hDistHistoryCurrent = new TH1F("hBaselineHistory","Average Baseline History",hDistHistoryPrevious->GetNbinsX()+1,0,hDistHistoryPrevious->GetNbinsX()+1);
			else if(type == 2) hDistHistoryCurrent = new TH1F("hNoiseUncHistory","Average Uncorrected Noise History",hDistHistoryPrevious->GetNbinsX()+1,0,hDistHistoryPrevious->GetNbinsX()+1);
			else if(type == 3) hDistHistoryCurrent = new TH1F("hNoiseCMNHistory","Average Common Mode Noise History",hDistHistoryPrevious->GetNbinsX()+1,0,hDistHistoryPrevious->GetNbinsX()+1);
			else if(type == 4) hDistHistoryCurrent = new TH1F("hNoiseCorHistory","Average Corrected Noise History",hDistHistoryPrevious->GetNbinsX()+1,0,hDistHistoryPrevious->GetNbinsX()+1);
			else if(type == 5) hDistHistoryCurrent = new TH1F("hGainHistory","Average Gain History",hDistHistoryPrevious->GetNbinsX()+1,0,hDistHistoryPrevious->GetNbinsX()+1);
			else if(type == 6) hDistHistoryCurrent = new TH1F("hDriftSpeedHistory","Average Drift Speed History",hDistHistoryPrevious->GetNbinsX()+1,0,hDistHistoryPrevious->GetNbinsX()+1);
			//AliInfo(Form("hBaselineHistory 2 %x\n",hDistHistory));
			for(Int_t i=0;i<hDistHistoryPrevious->GetNbinsX();i++) {
				hDistHistoryCurrent->SetBinContent(i,hDistHistoryPrevious->GetBinContent(i));
				hDistHistoryCurrent->SetBinError(i,hDistHistoryPrevious->GetBinError(i));
			}
			hDistHistoryCurrent->SetBinContent(hDistHistoryPrevious->GetNbinsX(),average);
			hDistHistoryCurrent->SetBinError(hDistHistoryPrevious->GetNbinsX(),sigma);
		}
	}
	hDistHistoryCurrent->Write();
	gainHistoryFile->Close();
	delete gainHistoryFile;
//	delete hname;
	delete hDist;
	delete hDistDiff;
//        if(hDistHistoryCurrent) delete hDistHistoryCurrent;
//	if(hDistHistoryPrevious) delete hDistHistoryPrevious;
}//_______________________________________________________________

