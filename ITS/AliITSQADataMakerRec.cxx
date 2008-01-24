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
//  W. Ferrarese Gen 2008
//  INFN Torino

// --- ROOT system ---
#include <TH2D.h>
#include <TBranch.h>
#include <TTree.h>
#include <TGaxis.h>
#include <TMath.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliITSQADataMakerRec.h"
#include "AliLog.h"
#include "AliQA.h"
#include "AliQAChecker.h"
#include "AliRawReader.h"
#include "AliITSRawStreamSDD.h"
#include "AliITSRecPoint.h"
#include "AliITSgeomTGeo.h"
#include "AliRawReader.h"

#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"


ClassImp(AliITSQADataMakerRec)

//____________________________________________________________________________ 
AliITSQADataMakerRec::AliITSQADataMakerRec() : 
AliQADataMakerRec(AliQA::GetDetName(AliQA::kITS), "SDD Quality Assurance Data Maker"),
fkOnline(kFALSE),
fLDC(0),
fnSDDHistos(0),
fSDDDDLModuleMap(0)
{ 
  // Default constructor 
}

//____________________________________________________________________________ 
AliITSQADataMakerRec::AliITSQADataMakerRec(Int_t ldc, AliITSDDLModuleMapSDD *sddDDLModuleMap, Bool_t kMode) :
AliQADataMakerRec(AliQA::GetDetName(AliQA::kITS), "SDD Quality Assurance Data Maker"),
fkOnline(kMode),
fLDC(ldc),
fnSDDHistos(0),
fSDDDDLModuleMap(sddDDLModuleMap)
{
  //ctor no more useful since the SDDModuleMap is taken in the InitRaws 
  //ctor used to discriminate OnLine-Offline analysis
}

//____________________________________________________________________________ 
AliITSQADataMakerRec::AliITSQADataMakerRec(Int_t ldc, Bool_t kMode) :
AliQADataMakerRec(AliQA::GetDetName(AliQA::kITS), "SDD Quality Assurance Data Maker"),
fkOnline(kMode),
fLDC(ldc),
fnSDDHistos(0),
fSDDDDLModuleMap(0)
{
  //ctor used to discriminate OnLine-Offline analysis
}

//____________________________________________________________________________ 
AliITSQADataMakerRec::AliITSQADataMakerRec(const AliITSQADataMakerRec& qadm) :
AliQADataMakerRec(qadm),
fkOnline(qadm.fkOnline),
fLDC(qadm.fLDC),
fnSDDHistos(qadm.fnSDDHistos),
fSDDDDLModuleMap(0)
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
void AliITSQADataMakerRec::StartOfDetectorCycle() const
{
  //Detector specific actions at start of cycle
  AliDebug(1,"AliITSQADM::Start of ITS Cycle\n");
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::EndOfDetectorCycle(AliQA::TASKINDEX /*task*/, TObjArray* /*list*/)
{
  // launch the QA checking
  AliDebug(1,"AliITSDM instantiates checker with Run(AliQA::kITS, task, list)\n"); 
  
  //AliQAChecker::Instance()->Run( AliQA::kITS , task, list);
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::EndOfDetectorCycle(const char * /* fgDataName */)
{
  //possibly used for different  AliQAChecker::Instance()->Run
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::InitRaws()
{  
  // Initialization for RAW data
	InitSPDRaws();
	InitSDDRaws();
	InitSSDRaws();
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::InitSPDRaws()
{  
  // Initialization for RAW data - SPD -
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::InitSDDRaws()
{  
  // create SDD histo of raw
  
  AliCDBEntry *ddlMapSDD = AliCDBManager::Instance()->Get("ITS/Calib/DDLMapSDD");
  if( !ddlMapSDD){
    AliFatal("Calibration object retrieval failed! ");
  }  
  AliITSDDLModuleMapSDD *ddlsdd=(AliITSDDLModuleMapSDD*)ddlMapSDD->GetObject();
  ddlMapSDD->SetOwner(kTRUE);
  fSDDDDLModuleMap = ddlsdd;

  Int_t lay, lad, det;
  
  if(fkOnline) {
    AliInfo("Book Online Histograms\n");
  }
  else {
    AliInfo("Book Offline Histograms\n ");
  }

  TH1F *h0 = new TH1F("ModPattern","Modules pattern",fgknSDDmodules,-0.5,259.5); 
  Add2RawsList(h0,0);
  fnSDDHistos++;
  TH1F *h1 = new TH1F("LadPatternL3","Ladder pattern L3",14,0.5,14.5);  
  Add2RawsList(h1,1);
  fnSDDHistos++;
  TH1F *h2 = new TH1F("LadPatternL4","Ladder pattern L4",22,0.5,22.5);  
  Add2RawsList(h2,2);
  fnSDDHistos++;

  Char_t *hname0[fgkLADDonLAY3] ; 
  TH1D *h3[fgkLADDonLAY3] ; 
  for(Int_t i=1; i<=fgkLADDonLAY3; i++) {
    hname0[i-1] = new Char_t[20];
    sprintf(hname0[i-1],"ModPattern_L3_%d",i);
    h3[i-1] = new TH1D(hname0[i-1],hname0[i-1],6,0.5,6.5);
    Add2RawsList(h3[i-1],i-1+3);
    fnSDDHistos++;
  }

  Char_t *hname1[fgkLADDonLAY4] ;
  TH1D *h4[fgkLADDonLAY4] ; 
  for(Int_t i=1; i<=fgkLADDonLAY4; i++) {
    hname1[i-1] = new Char_t[20];
    sprintf(hname1[i-1],"ModPattern_L4_%d",i);
    h4[i-1] = new TH1D(hname1[i-1],hname1[i-1],8,0.5,8.5);
    Add2RawsList(h4[i-1],i-1+17);  
    fnSDDHistos++;
  }

  if(fkOnline) {
    Int_t indexlast = 0;
    Char_t *hname2[fgknSDDmodules*2] ;
    TH1D *h5[fgknSDDmodules*2] ; 
    Int_t index1 = 0;
    for(Int_t moduleSDD=0; moduleSDD<fgknSDDmodules; moduleSDD++){
      if((moduleSDD >= 0 && moduleSDD < 36) || (moduleSDD >= 84 && moduleSDD < 180)) {
	for(Int_t iside=0; iside<fgknSide; iside++){
	  //Int_t index1 = moduleSDD * 2 + iside;
	  hname2[index1] = new Char_t[50];
	  AliITSgeomTGeo::GetModuleId(moduleSDD, lay, lad, det);
	  sprintf(hname2[index1],"ProjYMap_L%d_%d_%d_%d",lay,lad,det,iside);
	  
	  h5[index1] = new TH1D(hname2[index1],hname2[index1],256,-0.5,255.5);
	  Add2RawsList(h5[index1],index1+39);
	  fnSDDHistos++;
	  index1++;
	  indexlast = index1+39;
	}
      }
    }

    Char_t *hname3[fgknSDDmodules*8] ;
    TH1D *h6[fgknSDDmodules*8] ; 
    Int_t indextot = 0;
    Int_t indexlast1 = 0;
    for(Int_t htype=0; htype<4; htype++){
      for(Int_t moduleSDD=0; moduleSDD<fgknSDDmodules; moduleSDD++){
	if((moduleSDD >= 0 && moduleSDD < 36) || (moduleSDD >= 84 && moduleSDD < 180)) {
	  for(Int_t iside=0; iside<fgknSide; iside++){
	    //Int_t index1 = moduleSDD*2 + iside;
	    hname3[indextot] = new Char_t[50]; 
	    AliITSgeomTGeo::GetModuleId(moduleSDD+fgkmodoffset, lay, lad, det);
	    if(htype == 0) sprintf(hname3[indextot],"CountsVSAnode_L%d_%d_%d_%d",lay,lad,det,iside);
	    if(htype == 1) sprintf(hname3[indextot],"ChargeVSAnode_L%d_%d_%d_%d",lay,lad,det,iside);
	    if(htype == 2) sprintf(hname3[indextot],"CountsVSTbin_L%d_%d_%d_%d",lay,lad,det,iside);
	    if(htype == 3) sprintf(hname3[indextot],"ChargeVSTbin_L%d_%d_%d_%d",lay,lad,det,iside);
	    h6[indextot] = new TH1D(hname3[indextot],hname3[indextot],256,-0.5,255.5);
	    Add2RawsList(h6[indextot],indexlast + indextot);
	    fnSDDHistos++;
	    indextot++;
	    indexlast1 = indexlast + indextot;
	  }
	}
      }
    }
    
    Int_t indexlast2 = 0;
    Char_t *hname4[fgknSDDmodules*2];
    TH2D *h7[fgknSDDmodules*2] ;
    index1 = 0;
    for(Int_t moduleSDD=0; moduleSDD<fgknSDDmodules; moduleSDD++){
      if((moduleSDD >= 0 && moduleSDD < 36) || (moduleSDD >= 84 && moduleSDD < 180)) {
	for(Int_t iside=0; iside<fgknSide; iside++){
	  //Int_t index1 = moduleSDD * 2 + iside;
	  AliITSgeomTGeo::GetModuleId(moduleSDD+fgkmodoffset, lay, lad, det);
	  hname4[index1] = new Char_t[50]; 
	  sprintf(hname4[index1],"Anode_vs_Charge_L%d_%d_%d_%d",lay,lad,det,iside);
	  h7[index1] = new TH2D(hname4[index1],hname4[index1],fgknSDDmodules*2,-0.5,-0.5+fgknSDDmodules*2,256,0.5,256.5);
	  Add2RawsList(h7[index1],indexlast1 + index1);
	  fnSDDHistos++;
	  index1++;
	  indexlast2 = indexlast1 + index1;
	}
      }
    }  

    Char_t *hname[3][2 * fgknSDDmodules] ;
    index1 = 0;
    for(Int_t moduleSDD =0; moduleSDD<fgknSDDmodules; moduleSDD++){
      if((moduleSDD >= 0 && moduleSDD < 36) || (moduleSDD >= 84 && moduleSDD < 180)) {
	for(Int_t iside=0;iside<fgknSide;iside++){
	  //Int_t index1 = moduleSDD * 2 + iside ;
	  for(Int_t i=0; i<3; i++) hname[i][index1]= new Char_t[50];
	  AliITSgeomTGeo::GetModuleId(moduleSDD+fgkmodoffset, lay, lad, det);
	  sprintf(hname[0][index1],"chargeMap_L%d_%d_%d_%d",lay,lad,det,iside);
	  sprintf(hname[1][index1],"TotalCharge_L%d_%d_%d_%d",lay,lad,det,iside);
	  sprintf(hname[2][index1],"hmonoDMap_L%d_%d_%d_%d",lay,lad,det,iside);
	  fModuleChargeMap[index1] = new TH2D(hname[0][index1],hname[1][index1],256,-0.5,255.5,256,-0.5,255.5);
	  Add2RawsList(fModuleChargeMap[index1],indexlast2 + index1);
	  fnSDDHistos++;
	  fmonoD[index1] = new TH1D(hname[2][index1],hname[2][index1],256,-0.5,255.5);
	  index1++;
	}
      }
    }
  }


  AliDebug(1,Form("%d SDD histograms booked\n",fnSDDHistos));
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::InitSSDRaws()
{  
  // Initialization for RAW data - SSD -
}

//____________________________________________________________________________
void AliITSQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{ 
  // Fill QA for RAW
	MakeSPDRaws(rawReader);
	MakeSDDRaws(rawReader);
	MakeSSDRaws(rawReader);
}

//____________________________________________________________________________
void AliITSQADataMakerRec::MakeSPDRaws(AliRawReader* /* rawReader */)
{
  // Fill QA for RAW - SPD - 
}

//____________________________________________________________________________
void AliITSQADataMakerRec::MakeSDDRaws(AliRawReader* rawReader)
{ 

  //Fills Raw QA list of histos
  if(rawReader->GetType() != 7) return;  // skips non physical triggers
  AliDebug(1,"entering MakeRaws\n");             
  rawReader->SelectEquipment(17,fgkeqOffset,fgkeqOffset + fgkDDLidRange); 

  rawReader->Reset();                         
  AliITSRawStreamSDD s(rawReader); 
  s.SetDDLModuleMap(fSDDDDLModuleMap);
  Int_t lay, lad, det; 

  Int_t index=0;
  if(fkOnline) {
    for(Int_t moduleSDD =0; moduleSDD<fgknSDDmodules; moduleSDD++){
      if((moduleSDD >= 0 && moduleSDD < 36) || (moduleSDD >= 84 && moduleSDD < 180)) {
	for(Int_t iside=0;iside<fgknSide;iside++) {
	  if(fnSDDHistos > 39+12 * 132 + index) GetRawsData(39+12 *132 + index)->Reset();
	  index++;
	}
      }
    }
  }
  Int_t cnt = 0;
  while(s.Next()) {
    Int_t iddl = rawReader->GetDDLID() - fgkDDLIDshift;
    Int_t isddmod = s.GetModuleNumber(iddl,s.GetCarlosId());
    if(isddmod==-1){
      AliDebug(1,Form("Found module with iddl: %d, s.GetCarlosId: %d \n",iddl,s.GetCarlosId() ));
    }
    if(s.IsCompletedModule()) {
      AliDebug(1,Form("IsCompletedModule == KTRUE\n"));
      continue;
    } 
    Int_t coord1 = s.GetCoord1();
    Int_t coord2 = s.GetCoord2();
    Int_t signal = s.GetSignal();
    Int_t moduleSDD = isddmod - fgkmodoffset;
    if(moduleSDD < 0 || moduleSDD>fgknSDDmodules+fgkmodoffset) {
      AliDebug(1,Form( "Module SDD = %d, resetting it to 1 \n",moduleSDD));
      moduleSDD = 1;
    }
    GetRawsData(0)->Fill(moduleSDD); 
    
    AliITSgeomTGeo::GetModuleId(isddmod, lay, lad, det);
    //printf("modnumb %d, lay %d, lad %d, det %d \n",isddmod, lay, lad, det);
    Int_t ioffset = 3;
    Int_t iorder = 1;
    if(lay==4) { 
      ioffset += 14;
      iorder = 2;   
    } 
    GetRawsData(iorder)->Fill(lad);
    GetRawsData(ioffset+lad-1)->Fill(det); //-1 because ladder# starts from 1    
    
    Short_t iside = s.GetChannel();
    
    Int_t activeModule = moduleSDD - fgkmodoffset ;
    if(moduleSDD > 35) activeModule -= 48;
    Int_t index1 = activeModule * 2 + iside;
    
    if(index1<0){
      AliDebug(1,Form("Wrong index number %d - patched to 0\n",index1));
      index1 = 0;
    }
    
    if(fkOnline) {
      if(fnSDDHistos > 39+12 * 132 + index1) {
        GetRawsData(39+ 2 * 132 + index1)->Fill(coord1);
        GetRawsData(39+ 4 * 132 + index1)->Fill(coord1,signal);
        GetRawsData(39+ 6 * 132 + index1)->Fill(coord2);
        GetRawsData(39+ 8 * 132 + index1)->Fill(coord2,signal);
        GetRawsData(39+10 * 132 + index1)->Fill(signal,coord1);
        ((TH2D *)(GetRawsData(39+12 * 132 + index1)))->Fill(coord2, coord1, signal);
      }
    }
    cnt++;
    if(!(cnt%10000)) AliDebug(1,Form(" %d raw digits read",cnt));
  }
  AliDebug(1,Form("Event completed, %d raw digits read",cnt));
  
  if(fkOnline) {
    Int_t nBins = 256;
    for(Int_t moduleSDD=0; moduleSDD<fgknSDDmodules; moduleSDD++){
      if((moduleSDD >= 0 && moduleSDD < 36) || (moduleSDD >= 84 && moduleSDD < 180)) {
	for(Int_t iside=0; iside<fgknSide; iside++){
	  Int_t activeModule = moduleSDD - fgkmodoffset;
	  if(moduleSDD > 35) activeModule -= 48;
	  Int_t index1 = activeModule * 2 + iside;
	  if(fnSDDHistos > 39 + 2 * 132 + index1) {
	    fmonoD[index1] = ((TH2D *) (GetRawsData(39+12 * 132 + index1)))->ProjectionY();
	    for(Int_t bin=0; bin<nBins; bin++) GetRawsData(index1+39)->Fill(bin,fmonoD[index1]->GetBinContent(bin+1) );
	  }  
	}
      }
    }
    for(Int_t i=0; i<fnSDDHistos; i++){
      Int_t entries = static_cast<Int_t>(GetRawsData(i)->GetEntries());
      if(entries != 0)
	AliDebug(1,Form("histo %d, name %s , entries %d ",i,GetRawsData(i)->GetName(),entries));
    }
  }

}

//____________________________________________________________________________
void AliITSQADataMakerRec::MakeSSDRaws(AliRawReader* /* rawReader */)
{ 
  // Fill QA for RAW - SSD -
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::InitRecPoints()
{
  // Initialization for RECPOINTS
	InitSPDRecPoints();
	InitSDDRecPoints();
	InitSSDRecPoints();
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::InitSPDRecPoints()
{
  // Initialization for RECPOINTS - SPD -
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::InitSDDRecPoints()
{
  // create SDD histo of RecPoints
  
  TH1F * h0 = new TH1F("Lay3TotCh","Layer 3 total charge",1000,-0.5, 499.5);
  h0->GetXaxis()->SetTitle("ADC Counts");
  Add2RecPointsList(h0, 0);
 
  TH1F * h1 = new TH1F("Lay4TotCh","Layer 4 total charge",1000,-0.5, 499.5);
  h1->GetXaxis()->SetTitle("ADC Counts");
  Add2RecPointsList(h1, 1);

    
  Char_t *hisnam[3];
  TH1F *h2[3]; 
  for(Int_t i=1; i<=3; i++){
    hisnam[i] = new Char_t[50];
    sprintf(hisnam[i],"Charge_L3_Strip%d",i);
    h2[i] = new TH1F(hisnam[i],hisnam[i],1000,-0.5, 499.5);
    Add2RecPointsList(h2[i],i+1);
  }
  
  Char_t *hisnam2[4] ;
  TH1F *h3[4]; 
  for(Int_t i=1; i<=4; i++){
    hisnam2[i] = new Char_t[50];
    sprintf(hisnam2[i],"Charge_L4_Strip%d",i);
    h3[i] = new TH1F(hisnam2[i],hisnam2[i],1000,-0.5, 499.5);
    Add2RecPointsList(h3[i],i+4);
  }
  
  TH1F *h4 = new TH1F("ModPatternRP","Modules pattern RP",fgknSDDmodules,-0.5,259.5); 
  h4->GetXaxis()->SetTitle("Module #");
  Add2RecPointsList(h4,9);
  TH1F *h5 = new TH1F("ModPatternL3 RP","Ladder pattern L3 RP",14,0.5,14.5);  
  h5->GetXaxis()->SetTitle("Ladder #, Layer 3");
  Add2RecPointsList(h5,10);
  TH1F *h6 = new TH1F("ModPatternL4 RP","Ladder pattern L4 RP",22,0.5,22.5); 
  h6->GetXaxis()->SetTitle("Ladder #, Layer 4");
  Add2RecPointsList(h6,11);
  TH2F *h7 = new TH2F("Local Coord Distrib","Local Coord Distrib",110,-5.5,5.5,110,-5.5,5.5);
  h7->GetXaxis()->SetTitle("X local coord, drift, cm");
  h7->GetYaxis()->SetTitle("Z local coord, anode, cm");
  Add2RecPointsList(h7,12);
  TH2F *h8 = new TH2F("Global Coord Distrib","Global Coord Distrib",80,-40,40,80,-40,40);
  h8->GetYaxis()->SetTitle("Y glob coord, cm");
  h8->GetXaxis()->SetTitle("X glob coord, cm");
  Add2RecPointsList(h8,13);
  
  Char_t *name;
  TH1F *h9[2]; 
  for(Int_t iLay=0; iLay<=1; iLay++){
    name = new Char_t[50];
    sprintf(name,"hr_Layer%d",iLay+3);
    h9[iLay] = new TH1F(name,name,100,5.,35.);
    h9[iLay]->GetXaxis()->SetTitle("r (cm)");
    h9[iLay]->GetXaxis()->CenterTitle();
    Add2RecPointsList(h9[iLay],iLay+14);
  }
  TH1F *h10[2]; 
  for(Int_t iLay=0; iLay<=1; iLay++){
    name = new Char_t[50];
    sprintf(name,"hphi_Layer%d",iLay+3);
    h10[iLay] = new TH1F(name,name,100,-TMath::Pi(),TMath::Pi());
    h10[iLay]->GetXaxis()->SetTitle("#varphi (rad)");
    h10[iLay]->GetXaxis()->CenterTitle();
    Add2RecPointsList(h10[iLay],iLay+16);
  }

  
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::InitSSDRecPoints()
{
  // Initialization for RECPOINTS - SSD -
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::MakeRecPoints(TTree * clustersTree)
{
  // Fill QA for recpoints
	MakeSPDRecPoints(clustersTree);
	MakeSDDRecPoints(clustersTree);
	MakeSSDRecPoints(clustersTree);
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::MakeSPDRecPoints(TTree * /* clustersTree */)
{
  // Fill QA for recpoints - SPD -
}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::MakeSDDRecPoints(TTree * clustersTree)
{
  // makes data from RecPoints
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
  for(Int_t module=240; module<500;module++){
    branchRecP->GetEvent(module);
    npoints += recpoints->GetEntries();
    AliITSgeomTGeo::GetModuleId(module, lay, lad, det);
    //printf("modnumb %d, lay %d, lad %d, det %d \n",module, lay, lad, det);
    
    for(Int_t j=0;j<recpoints->GetEntries();j++){
      AliITSRecPoint *recp = (AliITSRecPoint*)recpoints->At(j);    
      GetRecPointsData(9)->Fill(module);
      recp->GetGlobalXYZ(cluglo);
      Float_t rad=TMath::Sqrt(cluglo[0]*cluglo[0]+cluglo[1]*cluglo[1]); 
      Float_t phi=TMath::ATan2(cluglo[1],cluglo[0]);
      if(recp->GetLayer() ==2) {
	GetRecPointsData(0)->Fill( recp->GetQ()) ;
	GetRecPointsData(10)->Fill(lad);
	GetRecPointsData(14)->Fill(rad);
	GetRecPointsData(16)->Fill(phi);
      }
      else if(recp->GetLayer() ==3) {
	GetRecPointsData(1)->Fill( recp->GetQ()) ;
	GetRecPointsData(11)->Fill(lad);
	GetRecPointsData(15)->Fill(rad);
	GetRecPointsData(17)->Fill(phi);
      }
      else AliWarning(Form("Wrong SDD Layer: %d",recp->GetLayer())); 
      GetRecPointsData(12)->Fill(recp->GetDetLocalX(),recp->GetDetLocalZ());
      GetRecPointsData(13)->Fill(cluglo[0],cluglo[1]);

    }
  }
  recpoints->Delete();
  delete recpoints;

}

//____________________________________________________________________________ 
void AliITSQADataMakerRec::MakeSSDRecPoints(TTree * /*clustersTree*/)
{
  // Fill QA for recpoints - SSD -
}


