/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                         *
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

/////////////////////////////////////////////////////////////////////
//
//       Class performing TOF Trigger
//       Cosmic_Multi_muon: Cosmic Multi-Muonic Event Trigger (L0)
//       ppMB: p-p Minimum Bias Event Trigger (L0)
//       UltraPer_Coll: Ultra-Peripheral Collision Event Trigger (L0)
//       High_Mult: High Multiplicity Event Trigger (L0)
//       Jet: Events with Jet Topology Trigger (L1)
//
//       A.Silenzi: added CTTM map,
//                        method to fill LTM matrix from raw data,
//                        method to retrieve the TOF pre-trigger for TRD detector
//                        
//
/////////////////////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TTree.h>
#include <TMath.h>

#include "AliLoader.h"
#include "AliLog.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliTriggerInput.h"
#include "AliRawReader.h"

#include "AliTOFRawStream.h"
#include "AliTOFrawData.h"
#include "AliTOFdigit.h"
#include "AliTOFGeometry.h"
#include "AliTOFTrigger.h"
#include "AliTOFTriggerMask.h"

#include "AliCDBManager.h"
#include "AliCDBEntry.h"


extern AliRun* gAlice;

AliTOFTriggerMask* AliTOFTrigger:: fTOFTrigMap=NULL;
AliTOFTriggerMask* AliTOFTrigger:: fTOFTrigMask=NULL;

Int_t AliTOFTrigger::fgFromTriggertoDCS[72] = {
  0,1,4,5, 8, 9,12,13,16,17,20,21,24,25,28,29,32,33,36,37,40,41,44,45,48,49,52,53,56,57,60,61,64,65,68,69,
  3,2,7,6,11,10,15,14,19,18,23,22,27,26,31,30,35,34,39,38,43,42,47,46,51,50,55,54,59,58,63,62,67,66,71,70
}; // dcs to trigger mapping

//-------------------------------------------------------------------------
ClassImp(AliTOFTrigger)

//----------------------------------------------------------------------
  AliTOFTrigger::AliTOFTrigger() :
    AliTriggerDetector(),
    fHighMultTh(1000),
    fppMBTh(4),//4:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    fMultiMuonTh(4),
    fUPTh(2),
    fdeltaminpsi(150), //150
    fdeltamaxpsi(170), //170
    fdeltaminro(70),
    fdeltamaxro(110),
    fstripWindow(2),
    fSel1(0),
    fSel2(0),
    fSel3(0),
    fSel4(0),
    fNCrateOn(0),
    fNMaxipadOn(0),
    fNMaxipadOnAll(0),
    fStartTimeHit(0.0),
    fTimeWidthTrigger(25.0)

{
  //main ctor
  for (Int_t i=0;i<kNCTTM;i++) fLTMarray[i] = kFALSE;

  for (Int_t i=0;i<kNLTM;i++){

    for (Int_t j=0;j<kNLTMchannels;j++){
      fLTMmatrix[i][j]=kFALSE;
    }
    if (i<kNCTTM){
      for (Int_t j=0;j<kNCTTMchannels;j++){
	fCTTMmatrixFront[i][j]=kFALSE;
	fCTTMmatrixBack[i][j]=kFALSE;
      }
    }
  }

  fPowerMask[0] = 1;
  for(Int_t i=1;i <= kNCTTMchannels;i++){
      fPowerMask[i] = fPowerMask[i-1]*2;
  }

  SetName("TOF");
  CreateInputs();

  if(!fTOFTrigMap) fTOFTrigMap = new AliTOFTriggerMask();

}

//----------------------------------------------------------------------

AliTOFTrigger::AliTOFTrigger(Int_t HighMultTh, Int_t ppMBTh, Int_t MultiMuonTh, Int_t UPTh, Float_t deltaminpsi, Float_t deltamaxpsi, Float_t deltaminro, Float_t deltamaxro, Int_t stripWindow,Float_t startTimeWindow,Float_t widthTimeWindow) :
  AliTriggerDetector(),
  fHighMultTh(HighMultTh),
  fppMBTh(ppMBTh),
  fMultiMuonTh(MultiMuonTh),
  fUPTh(UPTh),
  fdeltaminpsi(deltaminpsi),
  fdeltamaxpsi(deltamaxpsi),
  fdeltaminro(deltaminro),
  fdeltamaxro(deltamaxro),
  fstripWindow(stripWindow),
  fSel1(0),
  fSel2(0),
  fSel3(0),
  fSel4(0),
  fNCrateOn(0),
  fNMaxipadOn(0),
  fNMaxipadOnAll(0),
  fStartTimeHit(startTimeWindow),
  fTimeWidthTrigger(widthTimeWindow)
{
  //ctor with thresholds for triggers
  for (Int_t i=0;i<kNCTTM;i++) fLTMarray[i] = kFALSE;
  for (Int_t i=0;i<kNLTM;i++){
    for (Int_t j=0;j<kNLTMchannels;j++){
      fLTMmatrix[i][j]=kFALSE;
    }
    if (i<kNCTTM){
      for (Int_t j=0;j<kNCTTMchannels;j++){
	fCTTMmatrixFront[i][j]=kFALSE;
	fCTTMmatrixBack[i][j]=kFALSE;
      }
    }
  }

  fPowerMask[0] = 1;
  for(Int_t i=1;i <= kNCTTMchannels;i++){
      fPowerMask[i] = fPowerMask[i-1]*2;
  }

  SetName("TOF");
  CreateInputs();

  if(!fTOFTrigMap) fTOFTrigMap = new AliTOFTriggerMask();
}


#if 0 /*** COPY CONSTRUCTOR SUPPRESSED **/
//____________________________________________________________________________

AliTOFTrigger::AliTOFTrigger(const AliTOFTrigger & tr):
  AliTriggerDetector(tr),
  fHighMultTh(tr.fHighMultTh),
  fppMBTh(tr.fppMBTh),
  fMultiMuonTh(tr.fMultiMuonTh),
  fUPTh(tr.fUPTh),
  fdeltaminpsi(tr.fdeltaminpsi),
  fdeltamaxpsi(tr.fdeltamaxpsi),
  fdeltaminro(tr.fdeltaminro),
  fdeltamaxro(tr.fdeltamaxro),
  fstripWindow(tr.fstripWindow),
  fSel1(tr.fSel1),
  fSel2(tr.fSel2),
  fSel3(tr.fSel3),
  fSel4(tr.fSel4),
  fNCrateOn(tr.fNCrateOn),
  fNMaxipadOn(tr.fNMaxipadOn),
  fNMaxipadOnAll(tr.fNMaxipadOnAll)
{
  //copy ctor
  for (Int_t i=0;i<kNCTTM;i++) fLTMarray[i] = kFALSE;
  for (Int_t i=0;i<kNLTM;i++){
    for (Int_t j=0;j<kNLTMchannels;j++){
      fLTMmatrix[i][j]=tr.fLTMmatrix[i][j];
    }
    if (i<kNCTTM){
      for (Int_t j=0;j<kNCTTMchannels;j++){
	fCTTMmatrixFront[i][j]=tr.fCTTMmatrixFront[i][j];
	fCTTMmatrixBack[i][j]=tr.fCTTMmatrixBack[i][j];
      }
    }
  }

  fPowerMask[0] = 1;
  for(Int_t i=1;i <= kNCTTMchannels;i++){
      fPowerMask[i] = fPowerMask[i-1]*2;
  }

  SetName(tr.GetName());
  //fInputs=&(tr.GetInputs());
  CreateInputs();

}
#endif /*** COPY CONTRUCTOR SUPPRESSED ***/

//----------------------------------------------------------------------

void AliTOFTrigger::CreateInputs()
{
  // creating inputs
  // Do not create inputs again!!
  if( fInputs.GetEntriesFast() > 0 ) return;

  //LoadActiveMask();

  fInputs.AddLast(new AliTriggerInput("TOF_Cosmic_MultiMuon_L0","TOF",0));
  fInputs.AddLast(new AliTriggerInput("0OIN","TOF",0)); // was "TOF_pp_MB_L0"
  fInputs.AddLast(new AliTriggerInput("0OM2","TOF",0)); // was "TOF_PbPb_MB2_L0"
  fInputs.AddLast(new AliTriggerInput("0OM3","TOF",0)); // was "TOF_PbPb_MB3_L0"
  fInputs.AddLast(new AliTriggerInput("0OUP","TOF",0)); // was "TOF_UltraPer_Coll_L0"
  fInputs.AddLast(new AliTriggerInput("0OMU","TOF",0)); // new trigger (150 < DeltaPhi < 180) and 2 <= N_pad <= 6

  fInputs.AddLast(new AliTriggerInput("0OHM","TOF",0)); // was "TOF_High_Mult_L0"
  fInputs.AddLast(new AliTriggerInput("TOF_Jet_L1","TOF",0));

}

//----------------------------------------------------------------------
void AliTOFTrigger::Trigger() {
  fTOFTrigMap->ResetMask();
  if(!fTOFTrigMask) LoadActiveMask();

  //triggering method
  fSel1=0;
  fSel2=0;
  fSel3=0;
  fSel4=0;

  CreateLTMMatrix();

  Int_t nchonFront = 0;
  Int_t nchonBack = 0;
  Int_t nchonTot = 0;
  Int_t nSectOn = 0; // °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
  Int_t DeSlots = -1;// °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
  Int_t AntiDeSlots = -1;// °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
  Int_t nchonFrontBack = 0;
  Int_t nchonFront1 = 0;
  Int_t nchonBack1 = 0;
  Int_t nchonFrontBack1 = 0;
  Int_t mindeltapsi = (Int_t)fdeltaminpsi/10;
  Int_t maxdeltapsi = (Int_t)fdeltamaxpsi/10;
  Int_t mindeltaro = (Int_t)fdeltaminro/10;
  Int_t maxdeltaro = (Int_t)fdeltamaxro/10;

  for (Int_t i=0;i<kNCTTM;i++){
    for (Int_t j=0;j<kNCTTMchannels;j++){
      if (fCTTMmatrixFront[i][j]){
	nchonFront++;
	fTOFTrigMap->SetON(i,j);
      }
    }
  }

  for (Int_t i=kNCTTM;i<(kNCTTM*2);i++){
    for (Int_t j=0;j<kNCTTMchannels;j++){
      if (fCTTMmatrixBack[i-kNCTTM][j]){
	nchonBack++;
	fTOFTrigMap->SetON(i,j);
      }
    }
  }

  nchonTot = nchonFront + nchonBack;
//  fNMaxipadOn = nchonTot;
  for(Int_t i=0;i<kNCTTM;i++) { if(fLTMarray[i]) nSectOn++; } 

  //pp Minimum Bias Trigger
  if (nchonTot >= fppMBTh) {
    SetInput("0OIN");
    fSel1=1;
    //printf("0OIN - MB\n");
  }

  // PbPb MB
  if (nchonTot >= 2) {
    SetInput("0OM2");
    fSel2=1;
  }
  if (nchonTot >= 3) {
    SetInput("0OM3");
    fSel3=1;
  }

  //High Multiplicity Trigger
  if (nchonTot >= fHighMultTh) {
    SetInput("0OHM");
    //printf("0OHM - High Mult\n");
  }


  //MultiMuon Trigger
  nchonFront = 0;
  nchonBack = 0;
  nchonFrontBack = 0;

  Bool_t boolCTTMor = kFALSE;

  for (Int_t i=0;i<(kNCTTM/2);i++){
    Int_t iopp = i+kNCTTM/2;
    for (Int_t j=0;j<kNCTTMchannels;j++){
      if (fCTTMmatrixFront[i][j]){
	Int_t minj = j-fstripWindow;
	Int_t maxj = j+fstripWindow;
	if (minj<0) minj =0;
	if (maxj>=kNCTTMchannels) maxj = kNCTTMchannels-1;
	boolCTTMor = kFALSE;
	for (Int_t k = minj;k<=maxj;k++){
	  boolCTTMor |= fCTTMmatrixFront[iopp][k];
	}
	if (boolCTTMor) {
	  nchonFront++;
	}
      }

      if (fCTTMmatrixBack[i][j]){
	Int_t minj = j-fstripWindow;
	Int_t maxj = j+fstripWindow;
	if (minj<0) minj =0;
	if (maxj>=kNCTTMchannels) maxj =kNCTTMchannels-1;
	boolCTTMor = kFALSE;
	for (Int_t k = minj;k<=maxj;k++){
	  boolCTTMor |= fCTTMmatrixBack[iopp][k];
	}
	if (boolCTTMor) {
	  nchonBack++;
	}
      }
    }
  }

  nchonFrontBack = nchonFront+nchonBack;

  nchonFront1 = 0;
  nchonBack1 = 0;
  nchonFrontBack1 = 0;

  boolCTTMor = kFALSE;
  for (Int_t i=0;i<(kNCTTM/2);i++){
    Int_t i2max = (kNCTTM-1)-i+1;
    Int_t i2min = (kNCTTM-1)-i-1;
    if (i2max >=kNCTTM) i2max = kNCTTM-1;
    if (i2min==i) i2min = kNCTTM-1-i;
    for (Int_t j=0;j<kNCTTMchannels;j++){
      Int_t j2min = j-fstripWindow;
      Int_t j2max = j+fstripWindow;
      if (j2min<0) j2min =0;
      if (j2max>=kNCTTMchannels) j2max =kNCTTMchannels-1;
      if (fCTTMmatrixFront[i][j]){
	boolCTTMor = kFALSE;
	for (Int_t i2=i2min;i2<=i2max;i2++){
	  for (Int_t j2 = j2min;j2<=j2max;j2++){
	    boolCTTMor |= fCTTMmatrixFront[i2][j2];
	  }
	  if (boolCTTMor) {
	    nchonFront++;
	  }
	}
      }
      if (fCTTMmatrixBack[i][j]){
	  boolCTTMor = kFALSE;
	  for (Int_t i2=i2min;i2<=i2max;i2++){
	      for (Int_t j2 = j2min;j2<=j2max;j2++){
		  boolCTTMor |= fCTTMmatrixBack[i2][j2];
	      }
	  }
	  if (boolCTTMor) {
	      nchonBack++;
	  }
      }
    }
  }
  
  nchonFrontBack1 = nchonFront1+nchonBack1;
  
  if (nchonFrontBack >= fMultiMuonTh || nchonFrontBack1 >= fMultiMuonTh) {
      SetInput("TOF_Cosmic_MultiMuon_L0");
  }
  
  //Ultra-Peripheral collision Trigger
  
  // °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
  // DeSlots = (k+1)_Array Element - k_Array Element
  // AntiDeSlots = kNCTTM - DeSlots
  
  if(nchonTot >= 2 && nchonTot <= 6)// && nchonFront >=2 && nchonBack <=6)
  {  
      // printf("nHitMaxipad CLASSE: %i \n",fNMaxipadOn);
      // printf("Total Number per event of Switched-On sectors : %i \n", nSectOn);
      // printf("mindeltapsi %i \n", mindeltapsi);
      //printf("maxdeltapsi %i \n", maxdeltapsi); 
      for(Int_t i = 0; i < kNCTTM; i++){    
	  if(fLTMarray[i]){
	      // printf(" i-sect On:  %i\n",i);
	      for(Int_t j = i+1; j < kNCTTM; j++){
		  if(fLTMarray[j]) {
		      //  printf(" j-sect On:  %i\n",j);
		      DeSlots = j-i;
		      AntiDeSlots = kNCTTM - DeSlots;
		      //printf("DeSlots = %i \n",DeSlots);
		      //printf("AntiDeSlots = %i \n",AntiDeSlots);
		      if(DeSlots >= mindeltapsi && DeSlots <= maxdeltapsi){
			  fSel4=1;
			  SetInput("0OUP");
			  //printf("trigger On with DeSlot \n");
		      }
		      if(AntiDeSlots >= mindeltapsi && AntiDeSlots <= maxdeltapsi){
			  fSel4=1;
			  SetInput("0OUP");
			  //printf("trigger On with AntiDeSlot \n");
		      }
		      
		      
		      if(DeSlots >= mindeltaro && DeSlots <= maxdeltaro){
			  fSel4=1;
			  SetInput("0OUP");
			  //printf("trigger On with DeSlot \n");
		      }
		      if(AntiDeSlots >= mindeltaro && AntiDeSlots <= maxdeltaro){
			  fSel4=1;
			  SetInput("0OUP");
			  //printf("trigger On with AntiDeSlot \n");
		      }	
		      
		      if(DeSlots >= 15 && DeSlots <= 18){
			SetInput("0OMU");
		      }
		      else if(AntiDeSlots >= 15 && AntiDeSlots <= 18){
			SetInput("0OMU");
		      }	
		  }
	      }    
	  }
      }
  } 
}

//-----------------------------------------------------------------------------
void AliTOFTrigger::CreateLTMMatrix() {
  //creating LTMMatrix
  //initialization
  CreateLTMMatrixFromDigits();
  CreateCTTMMatrix();
}

//-------------------------------------------------------------------------

void AliTOFTrigger::CreateLTMMatrixFromDigits() {
  //
  // Create LTM matrix by TOF digits
  //

  //initialization
  for (Int_t i=0;i<kNCTTM;i++) fLTMarray[i]= kFALSE;
  for (Int_t i=0;i<kNLTM;i++){
    for (Int_t j=0;j<kNLTMchannels;j++){
      fLTMmatrix[i][j]=kFALSE;
    }
  }
  for (Int_t i=0;i<kNCTTM;i++){
      for (Int_t j=0;j<kNCTTMchannels;j++){
	  fCTTMmatrixFront[i][j]=kFALSE;
	  fCTTMmatrixBack[i][j]=kFALSE;
      }
  }


  AliRunLoader *rl;
  rl = AliRunLoader::Instance();

  Int_t ncurrevent = rl->GetEventNumber();
  rl->GetEvent(ncurrevent);

  AliLoader * tofLoader = rl->GetLoader("TOFLoader");

  tofLoader->LoadDigits("read");
  TTree *treeD = tofLoader->TreeD();
  if (treeD == 0x0) {
    AliFatal("AliTOFTrigger: Can not get TreeD");
    return;
  }

  TBranch *branch = treeD->GetBranch("TOF");
  if (!branch) {
    AliError("can't get the branch with the TOF digits !");
    return;
  }
  TClonesArray *tofDigits =new TClonesArray("AliTOFdigit",  1000);
  branch->SetAddress(&tofDigits);
  treeD->GetEvent(0);
  Int_t ndigits = tofDigits->GetEntriesFast();
  Int_t detind[5]; //detector indexes: 0 -> sector
  //                                   1 -> plate(modulo)
  //                                   2 -> strip
  //                                   3 -> padz
  //                                   4 -> padx

  for (Int_t i=0;i<ndigits;i++){
    AliTOFdigit * digit = (AliTOFdigit*)tofDigits->UncheckedAt(i);
    detind[0] = digit->GetSector();
    detind[1] = digit->GetPlate();
    detind[2] = digit->GetStrip();
    detind[3] = digit->GetPadz();
    detind[4] = digit->GetPadx();

    Int_t indexLTM[2] = {-1,-1};
    GetLTMIndex(detind,indexLTM);

    //Float_t timedigit = digit->GetTdc()*AliTOFGeometry::TdcBinWidth()*1E-3; // decalibrated time digit in ns
    Float_t timedigit = digit->GetTdcND()*AliTOFGeometry::TdcBinWidth()*1E-3; // time digit in ns

    Float_t pos[3];
    AliTOFGeometry::GetPosPar(detind, pos);
    Float_t length = 0.;
    for (Int_t ic = 0; ic < 3; ic++) length += pos[ic] * pos[ic];
    length = TMath::Sqrt(length);
    timedigit -= length * 0.0333564095198152043; // subtract the minimal time in ns for the current channel

    if(timedigit > fStartTimeHit - 0.5 && timedigit < fStartTimeHit + fTimeWidthTrigger - 0.5)
      fLTMmatrix[indexLTM[0]][indexLTM[1]] = kTRUE;

//    fLTMarray[indexLTM[0]%36] = kTRUE; //dimensione MAX array 36 = kNCTTM 
  }


  tofLoader->UnloadDigits();
  //   rl->UnloadgAlice();

}

//-----------------------------------------------------------------------------

void AliTOFTrigger::CreateLTMMatrixFromRaw(AliRawReader *fRawReader) {
  //
  // Create LTM matrix by TOF raw data
  //
  fTOFTrigMap->ResetMask();

  //initialization
  for (Int_t i=0;i<kNLTM;i++){
    for (Int_t j=0;j<kNLTMchannels;j++){
      fLTMmatrix[i][j]=kFALSE;
    }
  }
  for (Int_t i=0;i<kNCTTM;i++){
      for (Int_t j=0;j<kNCTTMchannels;j++){
	  fCTTMmatrixFront[i][j]=kFALSE;
	  fCTTMmatrixBack[i][j]=kFALSE;
      }
  }

  if(fRawReader){
    AliTOFRawStream * tofRawStream = new AliTOFRawStream();

    Int_t inholes = 0;

    //if(!GetLoader()->TreeS()) {MakeTree("S");  MakeBranch("S");}
    
    Clear();
    tofRawStream->SetRawReader(fRawReader);
    
    //ofstream ftxt;
    //if (fVerbose==2) ftxt.open("TOFsdigitsRead.txt",ios::app);
    
    TClonesArray staticRawData("AliTOFrawData",10000);
    staticRawData.Clear();
    TClonesArray * clonesRawData = &staticRawData;
    
    Int_t dummy = -1;
    Int_t detectorIndex[5] = {-1, -1, -1, -1, -1};
    Int_t digit[2];
    //Int_t track = -1;
    //Int_t last = -1;
    
    Int_t indexDDL = 0;
    Int_t iRawData = 0;
    AliTOFrawData *tofRawDatum = 0;
    for (indexDDL=0; indexDDL<AliDAQ::NumberOfDdls("TOF"); indexDDL++) {
      
      fRawReader->Reset();
      tofRawStream->LoadRawDataBuffersV2(indexDDL);
      
      clonesRawData = tofRawStream->GetRawData();
      if (clonesRawData->GetEntriesFast()!=0) AliInfo(Form(" TOF raw data number = %3d", clonesRawData->GetEntriesFast()));
      for (iRawData = 0; iRawData<clonesRawData->GetEntriesFast(); iRawData++) {
	
        tofRawDatum = (AliTOFrawData*)clonesRawData->UncheckedAt(iRawData);
	
        //if (tofRawDatum->GetTOT()==-1 || tofRawDatum->GetTOF()==-1) continue;
        if (tofRawDatum->GetTOF()==-1) continue;
	
        SetBit(indexDDL, tofRawDatum->GetTRM(), tofRawDatum->GetTRMchain(),
               tofRawDatum->GetTDC(), tofRawDatum->GetTDCchannel());

        dummy = detectorIndex[3];
        detectorIndex[3] = detectorIndex[4];//padz
        detectorIndex[4] = dummy;//padx

        digit[0] = tofRawDatum->GetTOF();
        digit[1] = tofRawDatum->GetTOT();

        dummy = detectorIndex[3];
        detectorIndex[3] = detectorIndex[4];//padx
        detectorIndex[4] = dummy;//padz

        // Do not reconstruct anything in the holes
        if (detectorIndex[0]==13 || detectorIndex[0]==14 || detectorIndex[0]==15 ) { // sectors with holes
          if (detectorIndex[1]==2) { // plate with holes
            inholes++;
            continue;
          }
        }

        tofRawDatum = 0;
      } // while loop

      clonesRawData->Clear();

    } // DDL Loop

    //if (fVerbose==2) ftxt.close();

    if (inholes) AliWarning(Form("Clusters in the TOF holes: %d",inholes));
    delete tofRawStream;
    tofRawStream = NULL;

  }

}
//-----------------------------------------------------------------------------
void AliTOFTrigger::PrepareTOFMapFromRaw(AliRawReader *fRawReader,Int_t deltaBC) {
  //
  // Create LTM matrix by TOF raw data
  //
  if(!fTOFTrigMap) fTOFTrigMap = new AliTOFTriggerMask();
  fTOFTrigMap->ResetMask();
  LoadActiveMask();

  if(fRawReader){
   AliTOFRawStream * tofRawStream = new AliTOFRawStream();

    tofRawStream->SetRawReader(fRawReader);
        
    TClonesArray staticRawData("AliTOFrawData",10000);
    staticRawData.Clear();
    TClonesArray * clonesRawData = &staticRawData;
    
    Int_t indexDDL = 0;
    Int_t iRawData = 0;
    AliTOFrawData *tofRawDatum = 0;
    for (indexDDL=0; indexDDL<AliDAQ::NumberOfDdls("TOF"); indexDDL++) {
      
      fRawReader->Reset();
      tofRawStream->LoadRawDataBuffersV2(indexDDL);
      
      clonesRawData = tofRawStream->GetRawData();
      for (iRawData = 0; iRawData<clonesRawData->GetEntriesFast(); iRawData++) {
	
        tofRawDatum = (AliTOFrawData*)clonesRawData->UncheckedAt(iRawData);
	
        //if (tofRawDatum->GetTOT()==-1 || tofRawDatum->GetTOF()==-1) continue;
        if (tofRawDatum->GetTOF()==-1) continue;

	Int_t nTRM = tofRawDatum->GetTRM();
	Int_t iChain = tofRawDatum->GetTRMchain();
	Int_t iTDC = tofRawDatum->GetTDC();
	Int_t iCH=tofRawDatum->GetTDCchannel();

	if(nTRM==3 && iTDC>=12 && iTDC<=14 && indexDDL%2==1){ // DDL number to LTM number mapping
	  Int_t iLTMindex=-1;
	  Int_t iChannelIndex=-1;
	  switch(indexDDL%AliTOFGeometry::NDDL()){
	  case 1:
	    iLTMindex=1;
	    break;
	  case 3:
	    iLTMindex=36;
	    break;
	  default:
	    break;
	  }
	  iLTMindex+=2*(Int_t)(indexDDL/AliTOFGeometry::NDDL());
	  if(iChain==0 && indexDDL<36)
	    iLTMindex--;
	  if(iChain==0 && indexDDL>=36)
	    iLTMindex++;
	  iChannelIndex=iCH+iTDC*AliTOFGeometry::NCh()-12*AliTOFGeometry::NCh();
	  Int_t index[2]={iLTMindex,iChannelIndex};
	  
	  UInt_t indexDDL    = fgFromTriggertoDCS[index[0]];


	  if(fTOFTrigMask->IsON(indexDDL,index[1]) && TMath::Abs(tofRawDatum->GetTOF()-deltaBC) < 500) fTOFTrigMap->SetON(index[0],index[1]);


//	  if(!fTOFTrigMask->IsON(indexDDL,index[1])){
//	    printf("TOF problem (%i,%i) - (%i,%i) \n",indexDDL,index[1],index[0],index[1]);
//	  }
//	  else
//	    printf("TOF time = %i\n",tofRawDatum->GetTOF());

	}
	
        tofRawDatum = 0;
      } // while loop

      clonesRawData->Clear();

    } // DDL Loop


    delete tofRawStream;
    tofRawStream = NULL;

  }

}
//-----------------------------------------------------------------------------
void AliTOFTrigger::PrepareTOFMapFromDigit(TTree *treeD, Float_t startTimeHit, Float_t timeWidthTrigger) {
  if(!fTOFTrigMap) fTOFTrigMap = new AliTOFTriggerMask();
  LoadActiveMask();

  fTOFTrigMap->ResetMask();
  if (treeD == 0x0) {
    return;
  }

  TBranch *branch = treeD->GetBranch("TOF");
  if (!branch) {
    return;
  }
  TClonesArray *tofDigits =new TClonesArray("AliTOFdigit",  1000);
  branch->SetAddress(&tofDigits);
  treeD->GetEvent(0);
  Int_t ndigits = tofDigits->GetEntriesFast();
  Int_t detind[5]; //detector indexes: 0 -> sector
  //                                   1 -> plate(modulo)
  //                                   2 -> strip
  //                                   3 -> padz
  //                                   4 -> padx

  for (Int_t i=0;i<ndigits;i++){
    AliTOFdigit * digit = (AliTOFdigit*)tofDigits->UncheckedAt(i);
    detind[0] = digit->GetSector();
    detind[1] = digit->GetPlate();
    detind[2] = digit->GetStrip();
    detind[3] = digit->GetPadz();
    detind[4] = digit->GetPadx();

    Float_t pos[3];
    AliTOFGeometry::GetPosPar(detind, pos);
    Float_t timedigit = digit->GetTdcND()*AliTOFGeometry::TdcBinWidth()*1E-3; // time digit in ns
    timedigit -= TMath::Sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2])*0.0333564095198152043; // correct for minimum time-of-flight
    if (!(timedigit>startTimeHit-0.5 && timedigit<startTimeHit+timeWidthTrigger-0.5)) continue;
    
    Int_t indexLTM[2] = {-1,-1};
    GetLTMIndex(detind,indexLTM);

    UInt_t channelCTTM = indexLTM[1]/2;
    UInt_t indexDDL    = fgFromTriggertoDCS[indexLTM[0]];

    // skip digits for channels not included in the trigger active mask
    if (!fTOFTrigMask->IsON(indexDDL,channelCTTM)) continue;
    if (!fTOFTrigMap->IsON(indexLTM[0],channelCTTM)) continue;
    fTOFTrigMap->SetON(indexLTM[0],channelCTTM);

    if(fTOFTrigMask->IsON(indexLTM[0],indexLTM[1])) fTOFTrigMap->SetON(indexLTM[0],indexLTM[1]);
  }

}
//-----------------------------------------------------------------------------
void AliTOFTrigger::GetLTMIndex(const Int_t * const detind, Int_t *indexLTM) {
  //
  // getting LTMmatrix indexes for current digit
  //

  if (detind[1]==0 || detind[1]==1 || (detind[1]==2 && detind[2]<=7)) {
    if (detind[4]<24){
      indexLTM[0] = detind[0]*2;
    }
    else {
      indexLTM[0] = detind[0]*2+1;
    }
  }
  else {
    if (detind[4]<24){
      indexLTM[0] = detind[0]*2+36;
    }
    else {
      indexLTM[0] = (detind[0]*2+1)+36;
    }
  }

  if (indexLTM[0]<36) {
    if (detind[1] ==0){
      indexLTM[1] = detind[2];
    }
    else if (detind[1] ==1){
      indexLTM[1] = detind[2]+19;
    }
    else if (detind[1] ==2){
      indexLTM[1] = detind[2]+19*2;
    }
    else{
      //      AliError("Smth Wrong!!!");
    }
  }
  else {
    if (detind[1] ==2){
      indexLTM[1] = detind[2]-8;
    }
    else if (detind[1] ==3){
      indexLTM[1] = detind[2]+7;
    }
    else if (detind[1] ==4){
      indexLTM[1] = detind[2]+26;
    }
    else{
      //      AliError("Smth Wrong!!!");
    }
  }

}
//-------------------------------------------------------------------------
/*
// to be checked because of warning problems
void AliTOFTrigger::PrintMap()
{
  //
  //
  //

  for(Int_t i = 0; i<kNLTM;i++) {
    if(i<36) {
      printf("| %d | %d |%d |%d |%d |%d |%d |%d |%d |%d |%d |%d |%d |%d |%d |%d |%d |%d |%d |%d |%d |%d |%d |%d |%d |\n",
	     (fCTTMmatrixFront[i][0])?1:0,(fCTTMmatrixFront[i][1])?1:0,(fCTTMmatrixFront[i][2])?1:0, \
	     (fCTTMmatrixFront[i][3])?1:0,(fCTTMmatrixFront[i][4])?1:0,(fCTTMmatrixFront[i][5])?1:0, \
	     (fCTTMmatrixFront[i][6])?1:0,(fCTTMmatrixFront[i][7])?1:0,(fCTTMmatrixFront[i][8])?1:0, \
	     (fCTTMmatrixFront[i][9])?1:0,(fCTTMmatrixFront[i][10])?1:0,(fCTTMmatrixFront[i][11])?1:0, \
	     (fCTTMmatrixFront[i][12])?1:0,(fCTTMmatrixFront[i][13])?1:0,(fCTTMmatrixFront[i][14])?1:0,	\
	     (fCTTMmatrixFront[i][15])?1:0,(fCTTMmatrixFront[i][16])?1:0,(fCTTMmatrixFront[i][17])?1:0,	\
	     (fCTTMmatrixFront[i][18])?1:0,(fCTTMmatrixFront[i][19])?1:0,(fCTTMmatrixFront[i][20])?1:0,	\
	     (fCTTMmatrixFront[i][21])?1:0,(fCTTMmatrixFront[i][22])?1:0,(fCTTMmatrixFront[i][23])?1:0);
    } else {
      printf("| %d | %d |%d |%d |%d |%d |%d |%d |%d |%d |%d |%d |%d |%d |%d |%d |%d |%d |%d |%d |%d |%d |%d |%d |%d |\n",
	     (fCTTMmatrixBack[i][0])?1:0,(fCTTMmatrixBack[i][1])?1:0,(fCTTMmatrixBack[i][2])?1:0, \
	     (fCTTMmatrixBack[i][3])?1:0,(fCTTMmatrixBack[i][4])?1:0,(fCTTMmatrixBack[i][5])?1:0, \
	     (fCTTMmatrixBack[i][6])?1:0,(fCTTMmatrixBack[i][7])?1:0,(fCTTMmatrixBack[i][8])?1:0, \
	     (fCTTMmatrixBack[i][9])?1:0,(fCTTMmatrixBack[i][10])?1:0,(fCTTMmatrixBack[i][11])?1:0, \
	     (fCTTMmatrixBack[i][12])?1:0,(fCTTMmatrixBack[i][13])?1:0,(fCTTMmatrixBack[i][14])?1:0, \
	     (fCTTMmatrixBack[i][15])?1:0,(fCTTMmatrixBack[i][16])?1:0,(fCTTMmatrixBack[i][17])?1:0, \
	     (fCTTMmatrixBack[i][18])?1:0,(fCTTMmatrixBack[i][19])?1:0,(fCTTMmatrixBack[i][20])?1:0, \
	     (fCTTMmatrixBack[i][21])?1:0,(fCTTMmatrixBack[i][22])?1:0,(fCTTMmatrixBack[i][23])?1:0);
    }
  }

}
*/
//-------------------------------------------------------------------------

void AliTOFTrigger::GetMapMatrix(Bool_t map[][24]) const
{
  //
  // Returns CTTM map
  //

  for(Int_t i = 0; i<kNLTM;i++)
    for(Int_t j = 0; j<kNCTTMchannels;j++)
      map[i][j]=(i<36)?fCTTMmatrixFront[i][j]:fCTTMmatrixBack[i-36][j];

}
//-------------------------------------------------------------------------

void AliTOFTrigger::GetMap(Bool_t **map) const
{
  //
  // Returns CTTM map
  //

  for(Int_t i = 0; i<kNLTM;i++)
    for(Int_t j = 0; j<kNCTTMchannels;j++)
      map[i][j]=(i<36)?fCTTMmatrixFront[i][j]:fCTTMmatrixBack[i-36][j];

}


//-------------------------------------------------------------------------
void AliTOFTrigger::GetTRDmap(Bool_t **map) const
{
  //
  // Retriev the bit map sent to the TRD detector
  //
    
    for(int i = 0; i<kNLTM;i++)
	for(int j = 0; j<kNLTMtoTRDchannels;j++){
	    map[i][j]=kFALSE;
	}

    for(int i = 0; i<kNLTM/2;i++)
    for(int j = 0; j<AliTOFTrigger::kNCTTMchannels;j++){
	UInt_t uTRDbit=j/3;
	if(fCTTMmatrixFront[i][j]) map[i][uTRDbit]=kTRUE;
    }
  for(int i = kNLTM/2; i<kNLTM;i++)
      for(int j = 0; j<AliTOFTrigger::kNCTTMchannels;j++){
	  UInt_t uTRDbit=j/3;
	  if(fCTTMmatrixBack[i-kNLTM/2][j]) map[i][uTRDbit]=kTRUE;
    }
  
}
//-------------------------------------------------------------------------
void AliTOFTrigger::GetTRDmapMatrix(Bool_t map[][8]) const
{
  //
  // Retriev the bit map sent to the TRD detector
  //
    
    for(int i = 0; i<kNLTM;i++)
	for(int j = 0; j<kNLTMtoTRDchannels;j++){
	    map[i][j]=kFALSE;
	}

    for(int i = 0; i<kNLTM/2;i++)
    for(int j = 0; j<AliTOFTrigger::kNCTTMchannels;j++){
	UInt_t uTRDbit=j/3;
	if(fCTTMmatrixFront[i][j]) map[i][uTRDbit]=kTRUE;
    }
  for(int i = kNLTM/2; i<kNLTM;i++)
      for(int j = 0; j<AliTOFTrigger::kNCTTMchannels;j++){
	  UInt_t uTRDbit=j/3;
	  if(fCTTMmatrixBack[i-kNLTM/2][j]) map[i][uTRDbit]=kTRUE;
    }
  
}

//-------------------------------------------------------------------------
void AliTOFTrigger::SetBit(Int_t *detind)
{
  //
  // Sets CTTM map element corresponding to detector element 'detind'
  //

  Int_t index[2];
  GetCTTMIndex(detind,index);
  if(index[0]<36)
    fCTTMmatrixFront[index[0]][index[1]]=kTRUE;
  else
    fCTTMmatrixBack[index[0]-36][index[1]]=kTRUE;

}

//-------------------------------------------------------------------------
void AliTOFTrigger::SetBit(Int_t nDDL, Int_t nTRM, Int_t iChain,
                           Int_t iTDC, Int_t iCH)
{
  //
  // Sets CTTM map element corresponding to equipment ID
  // labelled by number nDDL, nTRM, iChain, iTDC, iCH
  //

    if(nTRM==3 && iTDC>=12 && iTDC<=14 && nDDL%2==1){ // DDL number to LTM number mapping
//       getchar();
    Int_t iLTMindex=-1;
    Int_t iChannelIndex=-1;
    switch(nDDL%AliTOFGeometry::NDDL()){
    case 1:
      iLTMindex=1;
      break;
    case 3:
      iLTMindex=36;
      break;
    default:
      AliError("Call this function only if(nTRM==3 && iTDC>12 && iTDC<14 && nDDL%2==1) ");
      break;
    }
    iLTMindex+=2*(Int_t)(nDDL/AliTOFGeometry::NDDL());
    if(iChain==0 && nDDL<36)
      iLTMindex--;
    if(iChain==0 && nDDL>=36)
      iLTMindex++;
    iChannelIndex=iCH+iTDC*AliTOFGeometry::NCh()-12*AliTOFGeometry::NCh();
    Int_t index[2]={iLTMindex,iChannelIndex};
    if (index[0]<36){
      fCTTMmatrixFront[index[0]][index[1]]=kTRUE;
      fLTMmatrix[index[0]][index[1]*2]=kTRUE;
    }
    else{
	fCTTMmatrixBack[index[0]-36][index[1]]=kTRUE;
	fLTMmatrix[index[0]][index[1]*2]=kTRUE;
    }
    }

}
//-------------------------------------------------------------------------

void AliTOFTrigger::ResetBit(Int_t *detind)
{
  //
  // Sets CTTM map element corresponding to detector element 'detind'
  //

  Int_t index[2];
  GetCTTMIndex(detind,index);
  if(index[0]<36)
    fCTTMmatrixFront[index[0]][index[1]]=kFALSE;
  else
    fCTTMmatrixBack[index[0]-36][index[1]]=kFALSE;

}

//-------------------------------------------------------------------------
void AliTOFTrigger::ResetBit(Int_t nDDL, Int_t nTRM, Int_t iChain,
			     Int_t iTDC, Int_t iCH)
{
  //
  // Sets CTTM map element corresponding to equipment ID
  // labelled by number nDDL, nTRM, iChain, iTDC, iCH
  //

  if(nTRM==3 && iTDC>12 && iTDC<14 && nDDL%2==1){ // DDL number to LTM number mapping
    Int_t iLTMindex=-1;
    Int_t iChannelIndex=-1;
    switch(nDDL%AliTOFGeometry::NDDL()){
    case 1:
      iLTMindex=1;
      break;
    case 3:
      iLTMindex=36;
      break;
    default:
      AliError("Call this function only if(nTRM==3 && iTDC>12 && iTDC<14 && nDDL%2==1) ");
      break;
    }
    iLTMindex+=2*(Int_t)(nDDL/AliTOFGeometry::NDDL());
    if(iChain==0 && nDDL<36)
      iLTMindex--;
    if(iChain==0 && nDDL>=36)
      iLTMindex++;
    iChannelIndex=iCH+iTDC*AliTOFGeometry::NCh()-12*AliTOFGeometry::NCh();
    Int_t index[2]={iLTMindex,iChannelIndex};
    if (index[0]<36){
      fCTTMmatrixFront[index[0]][index[1]]=kFALSE;
    }
    else{
      fCTTMmatrixBack[index[0]-36][index[1]]=kFALSE;
    }
  }
  else
    AliError("Call this function only if(nTRM==3 && iTDC>12 && iTDC<14 && nDDL%2==1) ");

}
//-------------------------------------------------------------------------

Bool_t AliTOFTrigger::GetBit(Int_t *detind)
{
  //
  // Returns CTTM map element corresponding to detector element 'detind'
  //

  Int_t index[2];
  GetCTTMIndex(detind,index);
  return (index[0]<36)?fCTTMmatrixFront[index[0]][index[1]]:fCTTMmatrixBack[index[0]-36][index[1]];

}

//-------------------------------------------------------------------------
Bool_t AliTOFTrigger::GetBit(Int_t nDDL, Int_t nTRM, Int_t iChain,
                             Int_t iTDC, Int_t iCH)
{
  //
  // Returns CTTM map element corresponding to equipment ID
  // labelled by number nDDL, nTRM, iChain, iTDC, iCH
  //

  if ( !(nTRM==3 && iTDC>12 && iTDC<14 && nDDL%2==1) ) {
    AliWarning("Call this function only if(nTRM==3 && iTDC>12 && iTDC<14) ");
    return kFALSE;
  }
  //if (nTRM==3 && iTDC>12 && iTDC<14 && nDDL%2==1) { // DDL number to LTM number mapping

  UInt_t iLTMindex=0;
  UInt_t iChannelindex=0;
  switch(nDDL%AliTOFGeometry::NDDL()) {
  case 1:
    iLTMindex=1;
    break;
  case 3:
    iLTMindex=36;
    break;
  default:
    AliError("something wrong");
    break;
  }
  iLTMindex+=2*(Int_t)(nDDL/AliTOFGeometry::NDDL());

  if (iChain==1) return kFALSE; // AdC

  if (nDDL<36)
    iLTMindex--;
  if (nDDL>=36)
    iLTMindex++;
  iChannelindex=iCH+iTDC*AliTOFGeometry::NCh()-12*AliTOFGeometry::NCh();
  Int_t index[2]={static_cast<Int_t>(iLTMindex),static_cast<Int_t>(iChannelindex)};
  return (index[0]<36)?fCTTMmatrixFront[index[0]][index[1]]:fCTTMmatrixBack[index[0]-36][index[1]];

}

//-------------------------------------------------------------------------

void AliTOFTrigger::CreateCTTMMatrix() {
  //
  // Create CTTM bit map
  //

  LoadActiveMask();

  fNMaxipadOnAll=0;
  fNMaxipadOn=0;

  for(Int_t i = 0; i<kNLTM;i++){
	UInt_t currentMask = fPowerMask[kNCTTMchannels]-1;
	if(fTOFTrigMask) currentMask=fTOFTrigMask->GetTriggerMask(fgFromTriggertoDCS[i]);
	if(i<kNCTTM){
	    for(Int_t j = 0; j<kNCTTMchannels;j++){
		fCTTMmatrixFront[i][j]=fLTMmatrix[i][2*j]||fLTMmatrix[i][2*j+1];
		if(fCTTMmatrixFront[i][j]) fNMaxipadOnAll++;
		if(!(currentMask & fPowerMask[j])) fCTTMmatrixFront[i][j]=0;
		if(fCTTMmatrixFront[i][j]){
		    fNMaxipadOn++;
		    fLTMarray[i] = kTRUE;
		}
	    }
	}
	else{
	    for(Int_t j = 0; j<kNCTTMchannels;j++){
		fCTTMmatrixBack[i-kNCTTM][j]=fLTMmatrix[i][2*j]||fLTMmatrix[i][2*j+1];;
		if(fCTTMmatrixBack[i-kNCTTM][j]) fNMaxipadOnAll++;
		if(!(currentMask & fPowerMask[j])) fCTTMmatrixBack[i-kNCTTM][j]=0;
		if(fCTTMmatrixBack[i-kNCTTM][j]){
		    fNMaxipadOn++;
		    fLTMarray[i-kNCTTM] = kTRUE;
		}
	    }
	}
    }
  
    fNCrateOn = 0; 
    for(Int_t j=0; j < kNCTTM; j++) {if(fLTMarray[j]) fNCrateOn++;}

}     
//-----------------------------------------------------------------------------

void AliTOFTrigger::GetCTTMIndex(Int_t *detind, Int_t *indexCTTM) {
  //
  // Returns CTTM index corresponding to the detector element detind
  //

  GetLTMIndex(detind,indexCTTM);
  indexCTTM[1]/=2;

}
//-----------------------------------------------------------------------------
void AliTOFTrigger::LoadActiveMask(){
//
// Load OCDB current mask
//

    AliCDBManager *cdb = AliCDBManager::Instance();
    if(cdb->GetRun() < 0 || !(cdb->GetDefaultStorage())){
	if(!(cdb->GetDefaultStorage())){
	    cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
	    printf("AliTOFTrigger (WARNING): probably CDB first instance - Default Sorage set to \"local://$ALICE_ROOT/OCDB\"\n");
	}
	if(cdb->GetRun() < 0){
	    cdb->SetRun(0);
         printf("AliTOFTrigger (WARNING): probably CDB first instance - number of run set to 0\n");
	}
    }

    AliCDBEntry *cdbe = cdb->Get("TRIGGER/TOF/TriggerMask");
    if(!cdbe) return;
    fTOFTrigMask= (AliTOFTriggerMask *)cdbe->GetObject();
    
//     UInt_t maskArray[kNLTM];
//     if(fTOFTrigMask == NULL) fTOFTrigMask = new AliTOFTriggerMask();
//     for (Int_t k = 0; k < kNLTM ; k++) maskArray[k] = fPowerMask[kNCTTMchannels]-1;
//     //for (Int_t k = 0; k < kNLTM ; k+=2) maskArray[k] = 0;
    
//     fTOFTrigMask->SetTriggerMaskArray(maskArray);
}


//-----------------------------------------------------------------------------
AliTOFTrigger::~AliTOFTrigger()
{
  // dtor

}

//-----------------------------------------------------------------------------
AliTOFTrigger& AliTOFTrigger::operator=(const AliTOFTrigger &/*source*/)
{
  // ass. op.
  return *this;

}

