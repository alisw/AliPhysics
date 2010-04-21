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

extern AliRun* gAlice;

//-------------------------------------------------------------------------
ClassImp(AliTOFTrigger)

//----------------------------------------------------------------------
  AliTOFTrigger::AliTOFTrigger() :
    AliTriggerDetector(),
    fHighMultTh(1000),
    fppMBTh(4),
    fMultiMuonTh(2),
    fUPTh(2),
    fdeltaminpsi(150),
    fdeltamaxpsi(170),
    fdeltaminro(70),
    fdeltamaxro(110),
    fstripWindow(2)
{
  //main ctor
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
  SetName("TOF");
  CreateInputs();
}

//----------------------------------------------------------------------

AliTOFTrigger::AliTOFTrigger(Int_t HighMultTh, Int_t ppMBTh, Int_t MultiMuonTh, Int_t UPTh, Float_t deltaminpsi, Float_t deltamaxpsi, Float_t deltaminro, Float_t deltamaxro, Int_t stripWindow) :
  AliTriggerDetector(),
  fHighMultTh(HighMultTh),
  fppMBTh(ppMBTh),
  fMultiMuonTh(MultiMuonTh),
  fUPTh(UPTh),
  fdeltaminpsi(deltaminpsi),
  fdeltamaxpsi(deltamaxpsi),
  fdeltaminro(deltaminro),
  fdeltamaxro(deltamaxro),
  fstripWindow(stripWindow)

{
  //ctor with thresholds for triggers
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
  SetName("TOF");
  CreateInputs();
}

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
  fstripWindow(tr.fstripWindow)
{
  //copy ctor
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
  SetName(tr.GetName());
  CreateInputs();
  //fInputs=&(tr.GetInputs());
}

//----------------------------------------------------------------------

void AliTOFTrigger::CreateInputs()
{
  // creating inputs
  // Do not create inputs again!!
  if( fInputs.GetEntriesFast() > 0 ) return;

  fInputs.AddLast(new AliTriggerInput("TOF_Cosmic_MultiMuon_L0","TOF",0));
  fInputs.AddLast(new AliTriggerInput("0OIN","TOF",0)); // was "TOF_pp_MB_L0"
  fInputs.AddLast(new AliTriggerInput("0OX1","TOF",0)); // was "TOF_UltraPer_Coll_L0"

  fInputs.AddLast(new AliTriggerInput("0OHM","TOF",0)); // was "TOF_High_Mult_L0"
  fInputs.AddLast(new AliTriggerInput("TOF_Jet_L1","TOF",0));

}

//----------------------------------------------------------------------
void AliTOFTrigger::Trigger() {
  //triggering method

  CreateLTMMatrix();
  Int_t nchonFront = 0;
  Int_t nchonBack = 0;
  Int_t nchonTot = 0;
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
      if (fCTTMmatrixFront[i][j]) nchonFront++;
    }
  }

  for (Int_t i=kNCTTM;i<(kNCTTM*2);i++){
    for (Int_t j=0;j<kNCTTMchannels;j++){
      if (fCTTMmatrixBack[i-kNCTTM][j]) nchonBack++;
    }
  }

  nchonTot = nchonFront + nchonBack;

  //pp Minimum Bias Trigger
  if (nchonTot >= fppMBTh) {
    SetInput("0OIN");
  }

  //High Multiplicity Trigger
  if (nchonTot >= fHighMultTh) {
    SetInput("0OHM");
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
  Bool_t boolpsi = kFALSE;
  Bool_t boolro = kFALSE;
  if (nchonTot == fUPTh){
    for (Int_t i=0;i<kNCTTM;i++){
      for (Int_t j=0;j<kNCTTMchannels;j++){
	Int_t minipsi = i+mindeltapsi;
	Int_t maxipsi = i+maxdeltapsi;
	if (minipsi>=kNCTTM) minipsi = mindeltapsi-((kNCTTM-1)-i)-1;
	if (maxipsi>=kNCTTM) maxipsi = maxdeltapsi-((kNCTTM-1)-i)-1;
	Int_t miniro = i+mindeltaro;
	Int_t maxiro = i+maxdeltaro;
	if (miniro>=kNCTTM) miniro = mindeltaro-((kNCTTM-1)-i)-1;
	if (maxiro>=kNCTTM) maxiro = maxdeltaro-((kNCTTM-1)-i)-1;
	Int_t j2min = j-fstripWindow;
	Int_t j2max = j+fstripWindow;
	if (j2min<0) j2min =0;
	if (j2max>=kNCTTMchannels) j2max =kNCTTMchannels-1;
	if (fCTTMmatrixFront[i][j]){
	  for (Int_t i2=minipsi;i2<=maxipsi;i2++){
	    for (Int_t j2 = j2min;j2<=j2max;j2++){
	      if (fCTTMmatrixFront[i2][j2]) {
		SetInput("0OX1");
		boolpsi = kTRUE;
		//exiting loops
		j2 = j2max+1;
		i2 = maxipsi+1;
		j=kNCTTMchannels;
		i=kNCTTM;
	      }
	    }
	  }
	  if (!boolpsi){
	    for (Int_t i2=miniro;i2<=maxiro;i2++){
	      for (Int_t j2 = j2min;j2<=j2max;j2++){
		if (fCTTMmatrixFront[i2][j2]) {
		  SetInput("0OX1");
		  boolro = kTRUE;
		  //exiting loops
		  j2 = j2max+1;
		  i2 = maxiro+1;
		  j=kNCTTMchannels;
		  i=kNCTTM;
		}
	      }
	    }
	  }
	}

	else if (fCTTMmatrixBack[i][j]){
	  for (Int_t i2=minipsi;i2<=maxipsi;i2++){
	    for (Int_t j2 = j2min;j2<=j2max;j2++){
	      if (fCTTMmatrixBack[i2][j2]) {
		SetInput("0OX1");
		boolpsi = kTRUE;
		//exiting loops
		j2 = j2max+1;
		i2 = maxipsi+1;
		j=kNCTTMchannels;
		i=kNCTTM;
	      }
	    }
	  }
	  if (!boolpsi){
	    for (Int_t i2=miniro;i2<=maxiro;i2++){
	      for (Int_t j2 = j2min;j2<=j2max;j2++){
		if (fCTTMmatrixBack[i2][j2]) {
		  SetInput("0OX1");
		  boolro = kTRUE;
		  //exiting loops
		  j2 = j2max+1;
		  i2 = maxiro+1;
		  j=kNCTTMchannels;
		  i=kNCTTM;
		}
	      }
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
}

//-------------------------------------------------------------------------

void AliTOFTrigger::CreateLTMMatrixFromDigits() {
  //
  // Create LTM matrix by TOF digits
  //

  //initialization
  for (Int_t i=0;i<kNLTM;i++){
    for (Int_t j=0;j<kNLTMchannels;j++){
      fLTMmatrix[i][j]=kFALSE;
    }
  }
  AliRunLoader *rl;
  rl = AliRunLoader::Instance();

  Int_t ncurrevent = rl->GetEventNumber();
  rl->GetEvent(ncurrevent);

  AliLoader * tofLoader = rl->GetLoader("TOFLoader");

  tofLoader->LoadDigits("read");
  TTree *treeD = tofLoader->TreeD();
  if (treeD == 0x0)
    {
      AliFatal("AliTOFTrigger: Can not get TreeD");
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
  //                                   1 -> plate
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

    fLTMmatrix[indexLTM[0]][indexLTM[1]] = kTRUE;
  }

  tofLoader->UnloadDigits();
  //   rl->UnloadgAlice();
  CreateCTTMMatrix();

}

//-----------------------------------------------------------------------------

void AliTOFTrigger::CreateLTMMatrixFromRaw(AliRawReader *fRawReader) {
  //
  // Create LTM matrix by TOF raw data
  //

  //initialization
  for (Int_t i=0;i<kNLTM;i++){
    for (Int_t j=0;j<kNLTMchannels;j++){
      fLTMmatrix[i][j]=kFALSE;
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
      tofRawStream->LoadRawData(indexDDL);
      
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
      AliError("Smth Wrong!!!");
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
      AliError("Smth Wrong!!!");
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
    for(int j = 0; j<kNLTMtoTRDchannels;j++)
      map[i][j]=kFALSE;

  for(int i = 0; i<kNLTM/2;i++)
    for(int j = 0; j<AliTOFTrigger::kNCTTMchannels;j++){
      UInt_t uTRDbit=j/3;
      if(fCTTMmatrixFront[i][j]) map[i][uTRDbit]=kTRUE;
    }
  for(int i = kNLTM/2; i<kNLTM;i++)
    for(int j = 0; j<AliTOFTrigger::kNCTTMchannels;j++){
      UInt_t uTRDbit=j/3;
      if(fCTTMmatrixBack[i - kNLTM/2][j]) map[i][uTRDbit]=kTRUE;
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
    fCTTMmatrixBack[index[0]][index[1]]=kTRUE;

}

//-------------------------------------------------------------------------
void AliTOFTrigger::SetBit(Int_t nDDL, Int_t nTRM, Int_t iChain,
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
    if (index[0]<36)
      fCTTMmatrixFront[index[0]][index[1]]=kTRUE;
    else
      fCTTMmatrixBack[index[0]][index[1]]=kTRUE;
  }
  else
    AliError("Call this function only if(nTRM==3 && iTDC>12 && iTDC<14 && nDDL%2==1) ");

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
    fCTTMmatrixBack[index[0]][index[1]]=kFALSE;

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
    if (index[0]<36)
      fCTTMmatrixFront[index[0]][index[1]]=kFALSE;
    else
      fCTTMmatrixBack[index[0]][index[1]]=kFALSE;
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
  return (index[0]<36)?fCTTMmatrixFront[index[0]][index[1]]:fCTTMmatrixBack[index[0]][index[1]];

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
  Int_t index[2]={iLTMindex,iChannelindex};
  return (index[0]<36)?fCTTMmatrixFront[index[0]][index[1]]:fCTTMmatrixBack[index[0]][index[1]];

}

//-------------------------------------------------------------------------

void AliTOFTrigger::CreateCTTMMatrix() {
  //
  // Create CTTM bit map
  //

  for(Int_t i = 0; i<kNLTM;i++){
    if(i<kNCTTM){
      for(Int_t j = 0; j<kNCTTMchannels;j++)
        fCTTMmatrixFront[i][j]=fLTMmatrix[i][2*j]||fLTMmatrix[i][2*j+1];
    }
    else{
      for(Int_t j = 0; j<kNCTTMchannels;j++)
        fCTTMmatrixBack[i-kNCTTM][j]=fLTMmatrix[i][2*j]||fLTMmatrix[i][2*j+1];;
    }
  }
}     
//-----------------------------------------------------------------------------

void AliTOFTrigger::GetCTTMIndex(Int_t *detind, Int_t *indexCTTM) {
  //
  // Returns CTTM index corresponding to the detector element detind
  //

  GetLTMIndex(detind,indexCTTM);
  indexCTTM[1]/=2;

}
