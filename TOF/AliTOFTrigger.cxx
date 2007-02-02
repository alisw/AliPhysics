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
/////////////////////////////////////////////////////////////////////

#include <TClonesArray.h>

#include "AliLoader.h"
#include "AliLog.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliTriggerInput.h"

#include "AliTOFdigit.h"
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
  AliTriggerDetector(),
  fHighMultTh(0),
  fppMBTh(0),
  fMultiMuonTh(0),
  fUPTh(0),
  fdeltaminpsi(0),
  fdeltamaxpsi(0),
  fdeltaminro(0),
  fdeltamaxro(0),
  fstripWindow(0)
{
  //copy ctor
  fHighMultTh=tr.fHighMultTh;
  fppMBTh=tr.fppMBTh;
  fMultiMuonTh=tr.fMultiMuonTh;
  fUPTh=tr.fUPTh;
  fdeltaminpsi = tr.fdeltaminpsi;
  fdeltamaxpsi = tr.fdeltamaxpsi;
  fdeltaminro = tr.fdeltaminro;
  fdeltamaxro = tr.fdeltamaxro;
  fstripWindow = tr.fstripWindow;
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
   
   fInputs.AddLast(new AliTriggerInput("TOF_Cosmic_MultiMuon_L0","Cosmic Multimuon Topology",0x01));
   fInputs.AddLast(new AliTriggerInput("TOF_pp_MB_L0","pp Minimum Bias",0x02));
   fInputs.AddLast(new AliTriggerInput("TOF_UltraPer_Coll_L0","Ultra Peripheral Collisions",0x04));

   fInputs.AddLast(new AliTriggerInput("TOF_High_Mult_L0","High Multiplicity",0x08));
   fInputs.AddLast(new AliTriggerInput("TOF_Jet_L1","Jet Search",0x10));
}

//----------------------------------------------------------------------
void AliTOFTrigger::Trigger(){
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
      fCTTMmatrixFront[i][j]=kFALSE;
      fCTTMmatrixBack[i][j]=kFALSE;
    }
  }
  for (Int_t i=0;i<kNCTTM;i++){
    for (Int_t j=0;j<kNCTTMchannels;j++){
      fCTTMmatrixFront[i][j] = (fLTMmatrix[i][j*2] || fLTMmatrix[i][j*2+1]);
      if (fCTTMmatrixFront[i][j]) nchonFront++; 
    }
  }

  for (Int_t i=kNCTTM;i<(kNCTTM*2);i++){
    for (Int_t j=0;j<kNCTTMchannels;j++){
      fCTTMmatrixBack[i-kNCTTM][j] = (fLTMmatrix[i][j*2] || fLTMmatrix[i][j*2+1]);
      if (fCTTMmatrixBack[i-kNCTTM][j]) nchonBack++; 
    }
  }

  nchonTot = nchonFront + nchonBack;

  //pp Minimum Bias Trigger
  if (nchonTot >= fppMBTh) {
    SetInput("TOF_pp_MB_L0");
  }

  //High Multiplicity Trigger
  if (nchonTot >= fHighMultTh) {
    SetInput("TOF_High_Mult_L0");
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
		SetInput("TOF_UltraPer_Coll_L0");
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
		  SetInput("TOF_UltraPer_Coll_L0");
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
		SetInput("TOF_UltraPer_Coll_L0");
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
		  SetInput("TOF_UltraPer_Coll_L0");
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
void AliTOFTrigger::CreateLTMMatrix(){
  //creating LTMMatrix
  //initialization
  for (Int_t i=0;i<kNLTM;i++){
    for (Int_t j=0;j<kNLTMchannels;j++){
      fLTMmatrix[i][j]=kFALSE;
    }
  }
  AliRunLoader *rl;
  rl = gAlice->GetRunLoader();
  
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
}
//-----------------------------------------------------------------------------
void AliTOFTrigger::GetLTMIndex(Int_t *detind, Int_t *indexLTM){
  //getting LTMmatrix indexes for current digit
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

  if (indexLTM[0]<36){ 
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
  
