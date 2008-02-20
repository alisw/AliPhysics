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
fSDDhRaws(0),
fSDDhRecs(0),
fRawsOffset(0),
fRecsOffset(0),
fSDDDDLModuleMap(0)
{
  //ctor used to discriminate OnLine-Offline analysis
  if(fLDC < 0 || fLDC > 4) {
	AliError("Error: LDC number out of range; return\n");
  }
  for(Int_t i=0;i<2*fgknSDDmodules;i++){
    fModuleChargeMap[i] = NULL;
  }
}

//____________________________________________________________________________ 
AliITSQASDDDataMakerRec::AliITSQASDDDataMakerRec(const AliITSQASDDDataMakerRec& qadm) :
TObject(),
fAliITSQADataMakerRec(qadm.fAliITSQADataMakerRec),
fkOnline(qadm.fkOnline),
fLDC(qadm.fLDC),
fSDDhRaws(qadm.fSDDhRaws),
fSDDhRecs(qadm.fSDDhRecs),
fRawsOffset(qadm.fRawsOffset),
fRecsOffset(qadm.fRecsOffset),
fSDDDDLModuleMap(0)
{
  //copy ctor 
  fAliITSQADataMakerRec->SetName((const char*)qadm.fAliITSQADataMakerRec->GetName()) ; 
  fAliITSQADataMakerRec->SetTitle((const char*)qadm.fAliITSQADataMakerRec->GetTitle());
}

//____________________________________________________________________________ 
AliITSQASDDDataMakerRec::~AliITSQASDDDataMakerRec(){
  // destructor

  for(Int_t i=0;i<2*fgknSDDmodules;i++){
    if(fModuleChargeMap[i]) delete fModuleChargeMap[i];
  }
  
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
void AliITSQASDDDataMakerRec::EndOfDetectorCycle(AliQA::TASKINDEX /*task*/, TObjArray* /*list*/)
{
  // launch the QA checking
  AliDebug(1,"AliITSDM instantiates checker with Run(AliQA::kITS, task, list)\n"); 
  
  //AliQAChecker::Instance()->Run( AliQA::kITS , task, list);
}

//____________________________________________________________________________ 
void AliITSQASDDDataMakerRec::InitRaws()
{ 
  // Initialization for RAW data - SDD -
  fRawsOffset = (fAliITSQADataMakerRec->fRawsQAList)->GetEntries();

  AliCDBEntry *ddlMapSDD = AliCDBManager::Instance()->Get("ITS/Calib/DDLMapSDD");
  if( !ddlMapSDD){
    AliFatal("Calibration object retrieval failed! ");
  }  
  AliITSDDLModuleMapSDD *ddlsdd=(AliITSDDLModuleMapSDD*)ddlMapSDD->GetObject();
  ddlMapSDD->SetOwner(kTRUE);
  fSDDDDLModuleMap = ddlsdd;
 
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
  TH1F *h0 = new TH1F("ModPattern","Modules pattern",fgknSDDmodules,-0.5,259.5);
  h0->GetXaxis()->SetTitle("Module Number");
  h0->GetYaxis()->SetTitle("Counts");
  fAliITSQADataMakerRec->Add2RawsList(h0,0+fRawsOffset);
  fSDDhRaws++;
  if(fLDC==0 || fLDC==1 || fLDC==2){
    TH1F *h1 = new TH1F("LadPatternL3","Ladder pattern L3",14,0.5,14.5);  
    h1->GetXaxis()->SetTitle("Ladder Number on Lay3");
    h1->GetYaxis()->SetTitle("Counts");
    fAliITSQADataMakerRec->Add2RawsList(h1,1+fRawsOffset);
    fSDDhRaws++;
  }	
  if(fLDC==0 || fLDC==3 || fLDC==4){
    TH1F *h2 = new TH1F("LadPatternL4","Ladder pattern L4",22,0.5,22.5);  
    h2->GetXaxis()->SetTitle("Ladder Number on Lay4");
    h2->GetYaxis()->SetTitle("Counts");
    fAliITSQADataMakerRec->Add2RawsList(h2,2+fRawsOffset);
    fSDDhRaws++;
  }
  if(fLDC==0 || fLDC==1 || fLDC==2){
    TH1D *h3[fgkLADDonLAY3] ; 
    for(Int_t i=1; i<=fgkLADDonLAY3; i++) {
      sprintf(hname0,"ModPattern_L3_%d",i);
      h3[i-1] = new TH1D(hname0,hname0,6,0.5,6.5);
      h3[i-1]->GetXaxis()->SetTitle("Channel Number");
      h3[i-1]->GetYaxis()->SetTitle("Counts");
      fAliITSQADataMakerRec->Add2RawsList(h3[i-1],i-1+3+fRawsOffset);
      fSDDhRaws++;
    }
  }
  if(fLDC==0 || fLDC==3 || fLDC==4){
    TH1D *h4[fgkLADDonLAY4] ; 
    for(Int_t i=1; i<=fgkLADDonLAY4; i++) {
      sprintf(hname0,"ModPattern_L4_%d",i);
      h4[i-1] = new TH1D(hname0,hname0,8,0.5,8.5);
      h4[i-1]->GetXaxis()->SetTitle("Channel Number");
      h4[i-1]->GetYaxis()->SetTitle("Counts");
      fAliITSQADataMakerRec->Add2RawsList(h4[i-1],i-1+17+fRawsOffset);
      fSDDhRaws++;
    }
  }

  Int_t indexlast1 = 0;
  Int_t indexlast2 = 0;

  if(fkOnline) {
    indexlast = 0;
    TH1D *h5[fgknSDDmodules*2] ; 
    index1 = 0;

/*
    if(fLDC == 0 || fLDC == 1 || fLDC == 4){
      
	for(Int_t lad=0; lad<fgkLADDonLAY3; lad++){
	for(Int_t ichan=0; det<3; det++){
	  for(Int_t iside=0; iside<fgknSide; iside++){
	    sprintf(hname0,"ProjYMap_L%d_%d_%d_%d",lay,lad,det+sideflag*3,iside);  //need to define "lay"
	    h5[index1] = new TH1D(hname0,hname0,256,-0.5,255.5);
	    h5[index1]->GetXaxis()->SetTitle("Anode Number");
	    h5[index1]->GetYaxis()->SetTitle("Counts");
	    fAliITSQADataMakerRec->Add2RawsList(h5[index1],index1+39+fRawsOffset);
		fSDDhRaws++;
	    index1++;
	    indexlast = index1+39;
	  }
	}	
      }
    }
*/  
    for(Int_t moduleSDD=0; moduleSDD<fgknSDDmodules; moduleSDD++){
      if((moduleSDD >= 0 && moduleSDD < 36) || (moduleSDD >= 84 && moduleSDD < 180)) {
	for(Int_t iside=0; iside<fgknSide; iside++){
	  AliITSgeomTGeo::GetModuleId(moduleSDD+fgkmodoffset, lay, lad, det);
	  //if(fLDC == 0 || ((fLDC == 1 || fLDC == 2) && lay == 2) || ((fLDC == 3 || fLDC == 4) && lay == 3)
	  sprintf(hname0,"ProjYMap_L%d_%d_%d_%d",lay,lad,det,iside);
	  h5[index1] = new TH1D(hname0,hname0,256,-0.5,255.5);
	  h5[index1]->GetXaxis()->SetTitle("Anode Number");
	  h5[index1]->GetYaxis()->SetTitle("Counts");
	  fAliITSQADataMakerRec->Add2RawsList(h5[index1],index1+39+fRawsOffset);
	  fSDDhRaws++;
	  index1++;
	  indexlast = index1+39;
	}
      }
    }
    
    TH1D *h6[fgknSDDmodules*8] ; 
    Int_t indextot = 0;
    indexlast1 = 0;
    for(Int_t htype=0; htype<4; htype++){
      for(Int_t moduleSDD=0; moduleSDD<fgknSDDmodules; moduleSDD++){
	if((moduleSDD >= 0 && moduleSDD < 36) || (moduleSDD >= 84 && moduleSDD < 180)) {	 
	  for(Int_t iside=0; iside<fgknSide; iside++){
	    AliITSgeomTGeo::GetModuleId(moduleSDD+fgkmodoffset, lay, lad, det);
	    if(htype == 0) sprintf(hname0,"CountsVSAnode_L%d_%d_%d_%d",lay,lad,det,iside);
	    if(htype == 1) sprintf(hname0,"ChargeVSAnode_L%d_%d_%d_%d",lay,lad,det,iside);
	    if(htype == 2) sprintf(hname0,"CountsVSTbin_L%d_%d_%d_%d",lay,lad,det,iside);
	    if(htype == 3) sprintf(hname0,"ChargeVSTbin_L%d_%d_%d_%d",lay,lad,det,iside);
	    h6[indextot] = new TH1D(hname0,hname0,256,-0.5,255.5);
	    if(htype == 0){
	      h6[indextot]->GetXaxis()->SetTitle("Anode Number");
	      h6[indextot]->GetYaxis()->SetTitle("Counts");
	    }
	    if(htype == 1){
	      h6[indextot]->GetXaxis()->SetTitle("Anode Number");
	      h6[indextot]->GetYaxis()->SetTitle("Charge (ADC)");
	    }
	    if(htype == 2){
	      h6[indextot]->GetXaxis()->SetTitle("Time bin");
	      h6[indextot]->GetYaxis()->SetTitle("Counts");
	    }
	    if(htype == 3){
	      h6[indextot]->GetXaxis()->SetTitle("Time bin");
	      h6[indextot]->GetYaxis()->SetTitle("Charge (ADC)");
	    }
	    fAliITSQADataMakerRec->Add2RawsList(h6[indextot],indexlast + indextot + fRawsOffset);
	    fSDDhRaws++;
	    indextot++;
	    indexlast1 = indexlast + indextot;
	  }
	}
      }
    }


    indexlast2 = 0;
    index1 = 0;
    char *hname[3];
    for(Int_t i=0; i<3; i++) hname[i]= new char[50];
    for(Int_t moduleSDD =0; moduleSDD<fgknSDDmodules; moduleSDD++){
      if((moduleSDD >= 0 && moduleSDD < 36) || (moduleSDD >= 84 && moduleSDD < 180)) {
	for(Int_t iside=0;iside<fgknSide;iside++){
	AliITSgeomTGeo::GetModuleId(moduleSDD+fgkmodoffset, lay, lad, det);
	sprintf(hname[0],"chargeMap_L%d_%d_%d_%d",lay,lad,det,iside);
	sprintf(hname[1],"TotalCharge_L%d_%d_%d_%d",lay,lad,det,iside);
	sprintf(hname[2],"hmonoDMap_L%d_%d_%d_%d",lay,lad,det,iside);
	fModuleChargeMap[index1] = new TH2D(hname[0],hname[1],256,-0.5,255.5,256,-0.5,255.5);
	fModuleChargeMap[index1]->GetXaxis()->SetTitle("Time Bin");
	fModuleChargeMap[index1]->GetYaxis()->SetTitle("Anode");
	fAliITSQADataMakerRec->Add2RawsList((new TH2D(*fModuleChargeMap[index1])),indexlast1 + index1 + fRawsOffset);
	fSDDhRaws++;
       	index1++;	 
	indexlast2 = indexlast1 + index1;
      }
    }
  }
  
  
  TH2D *h7[fgknSDDmodules*2] ;
  index1 = 0;
  for(Int_t moduleSDD=0; moduleSDD<fgknSDDmodules; moduleSDD++){
    if((moduleSDD >= 0 && moduleSDD < 36) || (moduleSDD >= 84 && moduleSDD < 180)) {
      for(Int_t iside=0; iside<fgknSide; iside++){
	AliITSgeomTGeo::GetModuleId(moduleSDD+fgkmodoffset, lay, lad, det);
	sprintf(hname0,"Anode_vs_Charge_L%d_%d_%d_%d",lay,lad,det,iside);
	h7[index1] = new TH2D(hname0,hname0,fgknSDDmodules*2,-0.5,-0.5+fgknSDDmodules*2,256,0.5,256.5);
	h7[index1]->GetXaxis()->SetTitle("Charge Value (ADC)");
	h7[index1]->GetYaxis()->SetTitle("Anode");
	fAliITSQADataMakerRec->Add2RawsList(h7[index1],indexlast2 + index1 + fRawsOffset);
	fSDDhRaws++;
	index1++;
      }
    }
  }
  /*
  TH2F *h8[3];
  for(Int_t htype =0; htype<3; htype++){
    if(htype == 0) sprintf(hname0,"BLcompare");
    if(htype == 1) sprintf(hname0,"NoiseCompare");
    if(htype == 2) sprintf(hname0,"VdriftCompare");
    h8[htype]= new TH2F(hname0,hname0,520,0,260,100,-0.5,99.5);
    if(htype == 0){
      h8[htype]->GetXaxis()->SetTitle("Module Number");
      h8[htype]->GetYaxis()->SetTitle("Time #");
    }
    if(htype == 1){
      h8[htype]->GetXaxis()->SetTitle("Module Number");
      h8[htype]->GetYaxis()->SetTitle("Time #)");
    }
    if(htype == 2){
      h8[htype]->GetXaxis()->SetTitle("Module Number");
      h8[htype]->GetYaxis()->SetTitle("Time #");
    }

    fAliITSQADataMakerRec->Add2RawsList(h8[htype],indexlast2 + index1 + htype+1 + fRawsOffset);
    fSDDhRaws++;
    printf("added histos %c, at position %d \n", hname0, indexlast2 + index1 + htype+1 + fRawsOffset);
  }
  */
}  // kONLINE


  AliDebug(1,Form("%d SDD Raws histograms booked\n",fSDDhRaws));
}


//____________________________________________________________________________
void AliITSQASDDDataMakerRec::MakeRaws(AliRawReader* rawReader)
{ 
  // Fill QA for RAW - SDD -
  if(rawReader->GetType() != 7) return;  // skips non physical triggers
  AliDebug(1,"entering MakeRaws\n");                 
  rawReader->SelectEquipment(17,fgkeqOffset,fgkeqOffset + fgkDDLidRange); 

  /*
  if(rawReader->GetEventId()!=fEvtId){
    TFile *DAoutput = new TFile::Open(filename);
    TH1F *BLhisto;
    for(Int_t imod=0; imod<nSDDmodules; imod++){ 
      BLhisto = (TH1F*)DAoutput->Get("BLhistoname[imod]");
      mean = BLhisto->take mean;
      fAliITSQADataMakerRec->GetRawsData(i+1887)->Fill(mean, imod);
      Noisehisto....;
      Vdrifthisto...;
      <Q>histo....;
    }
  }
  fEvtId==rawReader->GetEventId();
  */

  rawReader->Reset();                         
  AliITSRawStreamSDD s(rawReader); 
  s.SetDDLModuleMap(fSDDDDLModuleMap);
  Int_t lay, lad, det; 

  Int_t index=0;
  if(fkOnline) {
    for(Int_t moduleSDD =0; moduleSDD<fgknSDDmodules; moduleSDD++){
      if((moduleSDD >= 0 && moduleSDD < 36) || (moduleSDD >= 84 && moduleSDD < 180)) {
        for(Int_t iside=0;iside<fgknSide;iside++) {
          if(fSDDhRaws > 39+12 * 132 + index) fAliITSQADataMakerRec->GetRawsData(39+12 * 132 + index +fRawsOffset)->Reset();
          index++;
        }
      }
    }
  }
  Int_t cnt = 0;
  Int_t ildcID = -1;
  Int_t iddl = -1;
  Int_t isddmod = -1;
  Int_t coord1, coord2, signal, moduleSDD, ioffset, iorder, activeModule, index1 ,nBins , iside, bin ;
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
    fAliITSQADataMakerRec->GetRawsData(0 +fRawsOffset)->Fill(moduleSDD); 
    
    AliITSgeomTGeo::GetModuleId(isddmod, lay, lad, det);
    ioffset = 3;
    iorder = 1;
    if(lay==4) { 
      ioffset += 14;
      iorder = 2;   
    } 
    fAliITSQADataMakerRec->GetRawsData(iorder +fRawsOffset)->Fill(lad);
    fAliITSQADataMakerRec->GetRawsData(ioffset+lad-1 +fRawsOffset)->Fill(det); //-1 because ladder# starts from 1    
    
    Short_t iside = s.GetChannel();
    activeModule = moduleSDD;
    if(moduleSDD > 35) activeModule -= 48;
    index1 = activeModule * 2 + iside;
    
    if(index1<0){
      AliDebug(1,Form("Wrong index number %d - patched to 0\n",index1));
      index1 = 0;
    }
    
    if(fkOnline) {
      if(fSDDhRaws > 39+12 * 132 + index1) {
	
        fAliITSQADataMakerRec->GetRawsData(39+ 2 * 132 + index1 +fRawsOffset)->Fill(coord1); 
        fAliITSQADataMakerRec->GetRawsData(39+ 4 * 132 + index1 +fRawsOffset)->Fill(coord1,signal); 
        fAliITSQADataMakerRec->GetRawsData(39+ 6 * 132 + index1 +fRawsOffset)->Fill(coord2);	
        fAliITSQADataMakerRec->GetRawsData(39+ 8 * 132 + index1 +fRawsOffset)->Fill(coord2,signal); 
        fAliITSQADataMakerRec->GetRawsData(39+12 * 132 + index1 +fRawsOffset)->Fill(signal,coord1);
        ((TH2D *)(fAliITSQADataMakerRec->GetRawsData(39+10 * 132 + index1 +fRawsOffset)))->Fill(coord2, coord1, signal);

      }
    }
    cnt++;
    if(!(cnt%10000)) AliDebug(1,Form(" %d raw digits read",cnt));
  }
  AliDebug(1,Form("Event completed, %d raw digits read",cnt));  

  if(fkOnline) {
    TH1D *ptr = NULL;
    nBins = 256;
    for(moduleSDD=0; moduleSDD<fgknSDDmodules; moduleSDD++){
      if((moduleSDD >= 0 && moduleSDD < 36) || (moduleSDD >= 84 && moduleSDD < 180)) {
	for(iside=0; iside<fgknSide; iside++){
	  activeModule = moduleSDD;
	  if(moduleSDD > 35) activeModule -= 48;
	  index1 = activeModule * 2 + iside;
	  if(fSDDhRaws > 39 + 2 * 132 + index1) {
	    ptr = ((TH2D *) (fAliITSQADataMakerRec->GetRawsData(39+10 * 132 + index1 +fRawsOffset)))->ProjectionY();
	    for(bin=0; bin<nBins; bin++) fAliITSQADataMakerRec->GetRawsData(index1+39 +fRawsOffset)->Fill(bin,ptr->GetBinContent(bin+1) );
	  }  
	}
      }
    }
    for(Int_t i=0; i<fSDDhRaws; i++){
      Int_t entries = static_cast<Int_t>(fAliITSQADataMakerRec->GetRawsData(i +fRawsOffset)->GetEntries());
      if(entries != 0)
	AliDebug(1,Form("histo %d, name %s , entries %d ",i,fAliITSQADataMakerRec->GetRawsData(i +fRawsOffset)->GetName(),entries));
    }
  }
}

//____________________________________________________________________________ 
void AliITSQASDDDataMakerRec::InitRecPoints()
{
  // Initialization for RECPOINTS - SDD -
  fRecsOffset = (fAliITSQADataMakerRec->fRecPointsQAList)->GetEntries();
  
  TH1F * h0 = new TH1F("Lay3TotCh","Layer 3 total charge",1000,-0.5, 499.5);
  h0->GetXaxis()->SetTitle("ADC value");
  h0->GetYaxis()->SetTitle("Entries");
  fAliITSQADataMakerRec->Add2RecPointsList(h0, 0 +fRecsOffset);
  fSDDhRecs++;
 
  TH1F * h1 = new TH1F("Lay4TotCh","Layer 4 total charge",1000,-0.5, 499.5);
  h1->GetXaxis()->SetTitle("ADC value");
  h1->GetYaxis()->SetTitle("Entries");
  fAliITSQADataMakerRec->Add2RecPointsList(h1, 1 +fRecsOffset);
  fSDDhRecs++;

    
  char hisnam[50];
  TH1F *h2[3]; 
  for(Int_t i=1; i<=3; i++){
    sprintf(hisnam,"Charge_L3_Strip%d",i);
    h2[i] = new TH1F(hisnam,hisnam,1000,-0.5, 499.5);
    fAliITSQADataMakerRec->Add2RecPointsList(h2[i],i+1 +fRecsOffset);
    fSDDhRecs++;
  }
  
  TH1F *h3[4]; 
  for(Int_t i=1; i<=4; i++){
    sprintf(hisnam,"Charge_L4_Strip%d",i);
    h3[i] = new TH1F(hisnam,hisnam,1000,-0.5, 499.5);
    fAliITSQADataMakerRec->Add2RecPointsList(h3[i],i+4 +fRecsOffset);
    fSDDhRecs++;
  }
  
  TH1F *h4 = new TH1F("ModPatternRP","Modules pattern RP",fgknSDDmodules,239.5,499.5); 
  h4->GetXaxis()->SetTitle("Module number");
  h4->GetYaxis()->SetTitle("Entries");
  fAliITSQADataMakerRec->Add2RecPointsList(h4,9 +fRecsOffset);
  fSDDhRecs++;
  TH1F *h5 = new TH1F("ModPatternL3 RP","Ladder pattern L3 RP",14,0.5,14.5);  
  h5->GetXaxis()->SetTitle("Ladder #, Layer 3");
  h5->GetYaxis()->SetTitle("Entries");
  fAliITSQADataMakerRec->Add2RecPointsList(h5,10 +fRecsOffset);
  fSDDhRecs++;
  TH1F *h6 = new TH1F("ModPatternL4 RP","Ladder pattern L4 RP",22,0.5,22.5); 
  h6->GetXaxis()->SetTitle("Ladder #, Layer 4");
  h6->GetYaxis()->SetTitle("Entries");
  fAliITSQADataMakerRec->Add2RecPointsList(h6,11 +fRecsOffset);
  fSDDhRecs++;
  TH2F *h7 = new TH2F("Local Coord Distrib","Local Coord Distrib",1000,-4,4,1000,-4,4);
  h7->GetXaxis()->SetTitle("X local coord, drift, cm");
  h7->GetYaxis()->SetTitle("Z local coord, anode, cm");
  fAliITSQADataMakerRec->Add2RecPointsList(h7,12 +fRecsOffset);
  fSDDhRecs++;
  TH2F *h8 = new TH2F("Global Coord Distrib","Global Coord Distrib",6000,-30,30,6000,-30,30);
  h8->GetYaxis()->SetTitle("Y glob coord, cm");
  h8->GetXaxis()->SetTitle("X glob coord, cm");
  fAliITSQADataMakerRec->Add2RecPointsList(h8,13 +fRecsOffset);
  fSDDhRecs++;
  
  TH1F *h9[2]; 
  for(Int_t iLay=0; iLay<=1; iLay++){
    sprintf(hisnam,"hr_Layer%d",iLay+3);
    h9[iLay] = new TH1F(hisnam,hisnam,100,10.,30.);
    h9[iLay]->GetXaxis()->SetTitle("r (cm)");
    h9[iLay]->GetXaxis()->CenterTitle();
    h9[iLay]->GetYaxis()->SetTitle("Entries");
    fAliITSQADataMakerRec->Add2RecPointsList(h9[iLay],iLay+14 +fRecsOffset);
    fSDDhRecs++;
  }
  TH1F *h10[2]; 
  for(Int_t iLay=0; iLay<=1; iLay++){
    sprintf(hisnam,"hphi_Layer%d",iLay+3);
    h10[iLay] = new TH1F(hisnam,hisnam,100,-TMath::Pi(),TMath::Pi());
    h10[iLay]->GetXaxis()->SetTitle("#varphi (rad)");
    h10[iLay]->GetXaxis()->CenterTitle();
    h10[iLay]->GetYaxis()->SetTitle("Entries");
    fAliITSQADataMakerRec->Add2RecPointsList(h10[iLay],iLay+16 +fRecsOffset);
    fSDDhRecs++;
  }

  AliDebug(1,Form("%d SDD Recs histograms booked\n",fSDDhRecs));
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
  TClonesArray * recpoints = new TClonesArray("AliITSRecPoint") ;
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
      fAliITSQADataMakerRec->GetRecPointsData(9 +fRecsOffset)->Fill(module);
      recp->GetGlobalXYZ(cluglo);
      Float_t rad=TMath::Sqrt(cluglo[0]*cluglo[0]+cluglo[1]*cluglo[1]); 
      Float_t phi=TMath::ATan2(cluglo[1],cluglo[0]);
      if(recp->GetLayer() ==2) {
	fAliITSQADataMakerRec->GetRecPointsData(0 +fRecsOffset)->Fill(recp->GetQ()) ;
	fAliITSQADataMakerRec->GetRecPointsData(10 +fRecsOffset)->Fill(lad);
	fAliITSQADataMakerRec->GetRecPointsData(14 +fRecsOffset)->Fill(rad);
	fAliITSQADataMakerRec->GetRecPointsData(16 +fRecsOffset)->Fill(phi);
	fAliITSQADataMakerRec->GetRecPointsData(9 +fRecsOffset)->Fill(module);
	fAliITSQADataMakerRec->GetRecPointsData(12 +fRecsOffset)->Fill(recp->GetDetLocalX(),recp->GetDetLocalZ());
	fAliITSQADataMakerRec->GetRecPointsData(13 +fRecsOffset)->Fill(cluglo[0],cluglo[1]);
      }
      else if(recp->GetLayer() ==3) {
	fAliITSQADataMakerRec->GetRecPointsData(1 +fRecsOffset)->Fill(recp->GetQ()) ;
	fAliITSQADataMakerRec->GetRecPointsData(11 +fRecsOffset)->Fill(lad);
	fAliITSQADataMakerRec->GetRecPointsData(15 +fRecsOffset)->Fill(rad);
	fAliITSQADataMakerRec->GetRecPointsData(17 +fRecsOffset)->Fill(phi);
	fAliITSQADataMakerRec->GetRecPointsData(9 +fRecsOffset)->Fill(module);
	fAliITSQADataMakerRec->GetRecPointsData(12 +fRecsOffset)->Fill(recp->GetDetLocalX(),recp->GetDetLocalZ());
	fAliITSQADataMakerRec->GetRecPointsData(13 +fRecsOffset)->Fill(cluglo[0],cluglo[1]);
      }
    }
  }
  recpoints->Delete();
  delete recpoints;

}

