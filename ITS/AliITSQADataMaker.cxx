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
//  W. Ferrarese Nov 2007
//  INFN Torino

// --- ROOT system ---
#include <TH2D.h>
#include <TBranch.h>
#include <TTree.h>
#include <TGaxis.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliITSQADataMaker.h"
#include "AliLog.h"
#include "AliQA.h"
#include "AliQAChecker.h"
#include "AliRawReader.h"
#include "AliITSRawStreamSDD.h"
#include "AliITSRecPoint.h"
#include "AliITSgeomTGeo.h"
#include "AliRawReader.h"


ClassImp(AliITSQADataMaker)

//____________________________________________________________________________ 
AliITSQADataMaker::AliITSQADataMaker() : 
  AliQADataMaker(AliQA::GetDetName(AliQA::kITS), "SDD Quality Assurance Data Maker")
{ 
  fkOnline = kFALSE;
  // ctor 
}

//____________________________________________________________________________ 
AliITSQADataMaker::AliITSQADataMaker(Int_t ldc, Bool_t kMode) :
  AliQADataMaker(AliQA::GetDetName(AliQA::kITS), "SDD Quality Assurance Data Maker")
{
  //ctor used to discriminate OnLine-Offline analysis
  fkOnline = kMode;
  fLDC = ldc; 
}

//____________________________________________________________________________ 
AliITSQADataMaker::AliITSQADataMaker(const AliITSQADataMaker& qadm) :
  AliQADataMaker()
{
  //copy ctor 
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliITSQADataMaker& AliITSQADataMaker::operator = (const AliITSQADataMaker& qac )
{
  // Equal operator.
  this->~AliITSQADataMaker();
  new(this) AliITSQADataMaker(qac);
  return *this;
}

//____________________________________________________________________________ 
void AliITSQADataMaker::StartOfDetectorCycle() const
{
  //Detector specific actions at start of cycle
  AliDebug(1,"AliITSQADM::Start of ITS Cycle\n");
}

//____________________________________________________________________________ 
void AliITSQADataMaker::EndOfDetectorCycle(AliQA::TASKINDEX task, TList *list)
{
  // launch the QA checking
  AliDebug(1,"AliITSDM instantiates checker with Run(AliQA::kITS, task, list)\n"); 
  
  AliQAChecker::Instance()->Run( AliQA::kITS , task, list);
}

//____________________________________________________________________________ 
void AliITSQADataMaker::EndOfDetectorCycle(const char * fgDataName)
{
  //eventually used for different  AliQAChecker::Instance()->Run
}

//____________________________________________________________________________ 
void AliITSQADataMaker::InitRaws()
{  
  // create SDD histo of raw
  
  Int_t lay, lad, det;
  
  if(fkOnline) {
    AliInfo("Book Online Histograms\n");
  }
  else {
    AliInfo("Book Offline Histograms\n ");
  }
  if(fkOnline) {
    Char_t *hname[3][2 * fgknSDDmodules] ;
    for(Int_t moduleSDD =0; moduleSDD<fgknSDDmodules; moduleSDD++){
      for(Int_t iside=0;iside<fgknSide;iside++){
	Int_t index1 = moduleSDD * 2 + iside ;
	for(Int_t i=0; i<3; i++) hname[i][index1]= new Char_t[50];
	sprintf(hname[0][index1],"chargeMap%dSide%d",moduleSDD,iside);
	sprintf(hname[1][index1],"Total Charge, module number %d , Side%d",moduleSDD,iside);
	sprintf(hname[2][index1],"hmonoDMap%dSide%d",moduleSDD,iside);
	fModuleChargeMap[index1] = new TH2D(hname[0][index1],hname[1][index1],256,-0.5,255.5,256,-0.5,255.5);
	fmonoD[index1] = new TH1D(hname[2][index1],hname[2][index1],256,-0.5,255.5);
      }
    }
  }
  
  TH1F *h0 = new TH1F("ModPattern","Modules pattern",fgknSDDmodules,-0.5,259.5); 
  Add2RawsList(h0,0);
  TH1F *h1 = new TH1F("ModPatternL3","Modules pattern L3",14,0.5,14.5);  
  Add2RawsList(h1,1);
  TH1F *h2 = new TH1F("ModPatternL4","Modules pattern L4",22,0.5,22.5);  
  Add2RawsList(h2,2);

  Char_t *hname0[fgkLADDonLAY3] ; 
  TH1D *h3[fgkLADDonLAY3] ; 
  for(Int_t i=1; i<=fgkLADDonLAY3; i++) {
    hname0[i-1] = new Char_t[20];
    sprintf(hname0[i-1],"ModPattern_L3_%d",i);
    h3[i-1] = new TH1D(hname0[i-1],hname0[i-1],6,0.5,6.5);
    Add2RawsList(h3[i-1],i-1+3);
  }

  Char_t *hname1[fgkLADDonLAY4] ;
  TH1D *h4[fgkLADDonLAY4] ; 
  for(Int_t i=1; i<=fgkLADDonLAY4; i++) {
    hname1[i-1] = new Char_t[20];
    sprintf(hname1[i-1],"ModPattern_L4_%d",i);
    h4[i-1] = new TH1D(hname1[i-1],hname1[i-1],8,0.5,8.5);
    Add2RawsList(h4[i-1],i-1+17);  
  }

  if(fkOnline) {
    Char_t *hname2[fgknSDDmodules*2] ;
    TH1D *h5[fgknSDDmodules*2] ; 
    for(Int_t moduleSDD=0; moduleSDD<fgknSDDmodules; moduleSDD++){
      for(Int_t iside=0; iside<fgknSide; iside++){
	Int_t index1 = moduleSDD * 2 + iside;
	hname2[index1] = new Char_t[50];
	sprintf(hname2[index1],"ProjYMap%dSide%d",moduleSDD+1,iside);
	h5[index1] = new TH1D(hname2[index1],hname2[index1],256,-0.5,255.5);
	Add2RawsList(h5[index1],index1+39);
      }
    }
    
    Char_t *hname3[fgknSDDmodules*8] ;
    TH1D *h6[fgknSDDmodules*8] ; 
    for(Int_t htype=0; htype<4; htype++){
      for(Int_t moduleSDD=0; moduleSDD<fgknSDDmodules; moduleSDD++){
	for(Int_t iside=0; iside<fgknSide; iside++){
	  Int_t index1 = moduleSDD*2 + iside;
	  hname3[index1 + htype*2*fgknSDDmodules] = new Char_t[50]; 
	  AliITSgeomTGeo::GetModuleId(moduleSDD+fgkmodoffset, lay, lad, det);
	  Int_t indextot = index1 + htype*2*fgknSDDmodules;
	  if(htype == 0) sprintf(hname3[indextot],"CountsVSAnode_L%d_%d_%d_%d",lay,lad,det,iside);
	  if(htype == 1) sprintf(hname3[indextot],"ChargeVSAnode_L%d_%d_%d_%d",lay,lad,det,iside);
	  if(htype == 2) sprintf(hname3[indextot],"CountsVSTbin_L%d_%d_%d_%d",lay,lad,det,iside);
	  if(htype == 3) sprintf(hname3[indextot],"ChargeVSTbin_L%d_%d_%d_%d",lay,lad,det,iside);
	  h6[indextot] = new TH1D(hname3[indextot],hname3[indextot],256,-0.5,255.5);
	  //cout << "add at index " << indextot+39+fgknSDDmodules*2 << " with name " << hname3[indextot] << endl;
	  Add2RawsList(h6[indextot],39 + fgknSDDmodules*2 + indextot);
	}			
      }
    }
    
    Char_t *hname4[fgknSDDmodules*2];
    TH2D *h7[fgknSDDmodules*2] ;
    for(Int_t moduleSDD=0; moduleSDD<fgknSDDmodules; moduleSDD++){
      for(Int_t iside=0; iside<fgknSide; iside++){
	Int_t index1 = moduleSDD * 2 + iside;
	AliITSgeomTGeo::GetModuleId(moduleSDD+fgkmodoffset, lay, lad, det);
	hname4[index1] = new Char_t[50]; 
	sprintf(hname4[index1],"Anode_vs_Charge_L%d_%d_%d_%d",lay,lad,det,iside);
	h7[index1] = new TH2D(hname4[index1],hname4[index1],fgknSDDmodules*2,-0.5,-0.5+fgknSDDmodules*2,256,0.5,256.5);
	Add2RawsList(h7[index1],2639 + index1);
      }
    }
  }
 
}

//____________________________________________________________________________
void AliITSQADataMaker::MakeRaws(AliRawReader* rawReader)
{ 
  //Fills Raw QA list of histos

  AliDebug(1,"entering MakeRaws\n");
  rawReader->RequireHeader(kFALSE);               
  rawReader->SelectEvents(7);                    
  rawReader->SelectEquipment(17,fgkeqOffset+1,fgkeqOffset + fgkDDLidRange); 
  rawReader->Reset();                         
  AliITSRawStreamSDD s(rawReader); 
  Int_t lay, lad, det; 

  Int_t cnt = 0;
  while(s.Next()){
    Int_t iddl = rawReader->GetDDLID() - fgkDDLIDshift;
    Int_t isddmod = s.GetModuleNumber(iddl,s.GetCarlosId());
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
    //AliLog("modnumb %d, lay %d, lad %d, det %d \n",isddmod, lay, lad, det);
    Int_t ioffset = 3;
    Int_t iorder = 1;
    if(lay==4) { 
      ioffset += 14;
      iorder = 2;   
    } 
    GetRawsData(iorder)->Fill(lad);
    GetRawsData(ioffset+lad-1)->Fill(det); //-1 because ladder# starts from 1    
    
    Short_t iside = s.GetChannel();  
    Int_t index1 = moduleSDD * 2 + iside;
    
    if(index1<0){
      AliDebug(1,Form("Wrong index number %d - patched to 0\n",index1));
      index1 = 0;
    }
    
    if(fkOnline) {
      fModuleChargeMap[index1]->Fill(coord2, coord1, signal);
    
      GetRawsData(39+ 2 * fgknSDDmodules + index1)->Fill(coord1); 
      GetRawsData(39+ 4 * fgknSDDmodules + index1)->Fill(coord1,signal); 
      GetRawsData(39+ 6 * fgknSDDmodules + index1)->Fill(coord2);	
      GetRawsData(39+ 8 * fgknSDDmodules + index1)->Fill(coord2,signal); 
      GetRawsData(2639+index1)->Fill(signal,coord1);
    }
    cnt++;
    if(!(cnt%10000)) AliDebug(1,Form(" %d raw digits read",cnt));
  }
  AliDebug(1,Form("Event completed, %d raw digits read",cnt));
  
  if(fkOnline) {
    Int_t nBins = 256;
    for(Int_t moduleSDD=0; moduleSDD<fgknSDDmodules; moduleSDD++){
      for(Int_t iside=0; iside<fgknSide; iside++){
	Int_t index1 = moduleSDD * 2 + iside;
	fmonoD[index1] = fModuleChargeMap[index1]->ProjectionY();
	for(Int_t bin=0; bin<nBins; bin++)
	  GetRawsData(index1+39)->Fill(bin,fmonoD[index1]->GetBinContent(bin+1) );
      }
    }  
    for(Int_t i=0; i<3159; i++){
      Int_t entries = static_cast<Int_t>(GetRawsData(i)->GetEntries());
      if(entries != 0)
	AliDebug(1,Form("histo %d, name %s , entries %d ",i,GetRawsData(i)->GetName(),entries));
    }
    rawReader->RequireHeader(kTRUE); 
  }
}



//____________________________________________________________________________ 
void AliITSQADataMaker::InitRecPoints()
{
  // create SDD histo of RecPoints
  
  TH1F * h0 = new TH1F("Lay3TotCh","Layer 3 total charge",1000,-0.5, 499.5);
  Add2RecPointsList(h0, 0);
 
  TH1F * h1 = new TH1F("Lay4TotCh","Layer 4 total charge",1000,-0.5, 499.5);
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
  
}


//____________________________________________________________________________ 
void AliITSQADataMaker::MakeRecPoints(TTree * clustersTree)
{
  // makes data from RecPoints
  /*  
  TBranch *branchRecP = clustersTree->GetBranch("ITSRecPoints");
  if (!branchRecP) { 
    AliError("can't get the branch with the ITS clusters !");
    return;
  }
  TObjArray * recpoints = new TObjArray(100) ;
  branchRecP->SetAddress(&recpoints);
  branchRecP->GetEntry(0);
  
  TIter next(recpoints);
  AliITSRecPoint * rp;
  while ( rp = dynamic_cast<AliITSRecPoint *>(next())  ) {
    if(rp->GetLayer() ==3) GetRecPointsData(0)->Fill( rp->GetQ()) ;
    else if(rp->GetLayer() ==4) GetRecPointsData(1)->Fill( rp->GetQ()) ;
  }
  recpoints->Delete();
  delete recpoints;
  */
}
