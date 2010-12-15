/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// Base Class used to find                                                //
// the reconstructed points for ITS Upgrade                               // 
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliITSUpgradeClusterFinder.h"
#include "AliITSRawStreamSPD.h"
#include "AliLog.h"
#include <string.h>
#include <TObjString.h>


AliITSUpgradeClusterFinder::AliITSUpgradeClusterFinder() : 
  fNhitsLeft(0),
  fOldModule(-1),
  fClusterTypeFlag(kTRUE),
  fFindClusterType(kFALSE),
  fClusterTypeOrigCol(0),
  fClusterTypeOrigRow(0),
  fColSum(0),fRowSum(0),
  fCharge(0),
  fClusterWidthMaxCol(0),
  fClusterWidthMinCol(0),
  fClusterWidthMaxRow(0),
  fClusterWidthMinRow(0),
  fChargeArray(0x0)
{ 
  fChargeArray = new TObjArray();
  fTmpLabel[0]=-5;
  fTmpLabel[1]=-5;
  fTmpLabel[2]=-5;
   for(int il=0; il<10;il++) fLabels[il]=-5;
}

//___________________________________________________________________________________
AliITSUpgradeClusterFinder::~AliITSUpgradeClusterFinder() {}
//___________________________________________________________________________________
void AliITSUpgradeClusterFinder::StartEvent() {
  NewEvent();
}
//___________________________________________________________________________________
Int_t AliITSUpgradeClusterFinder::ProcessHitOnline(Int_t layer,  UInt_t col, UInt_t row, UShort_t charge, Int_t label[3]) {
  if(layer < 6 ) { 
    return ProcessHit(layer,col,row, charge, label);
  }
  else {
    printf("ERROR: UpgradeClusterFinder::ProcessHitOnline: Out of bounds: layer,col,row, charge = %d,%d,%d,%d\n",layer ,col,row,charge);
    return 1;
  }
}
//___________________________________________________________________________________
Int_t AliITSUpgradeClusterFinder::ProcessHit(Int_t layer , UInt_t col, UInt_t row, UShort_t charge, Int_t label[3]) {
  
  if (layer>=6) {
   printf("ERROR: UpgradeClusterFinder::ProcessHit: Out of bounds: layer ,col,row, charge, label(0,1,2)= %d,%d,%d,%d,(%d,%d,%d)\n",layer ,col,row,charge,label[0],label[1],label[2]); 
   return 1;
  }

  // do we need to perform clustering on previous module?
  if (fOldModule!=-1 && (Int_t)layer!=fOldModule) {
    fChargeArray->AddLast(new TObjString(Form("%i %i %i %i %i %i",col,row,charge,label[0],label[1],label[2])));
    DoModuleClustering(fOldModule,charge);
    NewModule();
  }
  // add hit
  fChargeArray->AddLast(new TObjString(Form("%i %i %i %i %i %i",col,row,charge,label[0],label[1],label[2])));

  fOldModule=layer;
  fHitCol[fNhitsLeft]=col;
  fHitRow[fNhitsLeft]=row;
  fHits[col][row]=kTRUE;
  fTmpLabel[0]=label[0];
  fTmpLabel[1]=label[1];
  fTmpLabel[2]=label[2];
  fNhitsLeft++;
  return 0;
}
//___________________________________________________________________________________
void AliITSUpgradeClusterFinder::FinishEvent() {
  if (fNhitsLeft>0) {
    DoModuleClustering(fOldModule,fCharge);
  }
}
//___________________________________________________________________________________
UInt_t AliITSUpgradeClusterFinder::GetClusterCount(Int_t layer) {
  if (layer>=6) {
    printf("ERROR: UpgradeClusterFinder::GetClusterCount: Module out of bounds: layer = %d\n",layer);
    return 0;
  }
  return fClusterList[layer].GetNrEntries();
}
//___________________________________________________________________________________
Float_t AliITSUpgradeClusterFinder::GetClusterMeanCol(Int_t layer, UInt_t index) {
  if (layer>=6) {
    printf("ERROR: UpgradeClusterFinder::GetClusterMeanCol: Module out of bounds: layer = %d\n",layer);
    return 0;
  }
  return fClusterList[layer].GetColIndex(index);
}
//___________________________________________________________________________________
Float_t AliITSUpgradeClusterFinder::GetClusterMeanRow(Int_t layer, UInt_t index) {
  if (layer>=6) {
    printf("ERROR: UpgradeClusterFinder::GetClusterMeanRow: Module out of bounds: layer = %d\n",layer);
    return 0;
  }
  return fClusterList[layer].GetRowIndex(index);
}
//___________________________________________________________________________________
UInt_t AliITSUpgradeClusterFinder::GetClusterSize(Int_t layer, UInt_t index) {
  if (layer>=6) {
    printf("ERROR: UpgradeClusterFinder::GetClusterSize: Module out of bounds: layer = %d\n",layer);
    return 0;
  }
  return fClusterList[layer].GetSizeIndex(index);
}
//___________________________________________________________________________________
UInt_t AliITSUpgradeClusterFinder::GetClusterWidthZ(Int_t layer, UInt_t index) {
  if (layer>=6) {
    printf("ERROR: UpgradeClusterFinder::GetClusterWidthZ: Module out of bounds: layer = %d\n",layer);
    return 0;
  }
  return fClusterList[layer].GetWidthZIndex(index);
}
//___________________________________________________________________________________
UInt_t AliITSUpgradeClusterFinder::GetClusterWidthPhi(Int_t layer, UInt_t index) {
  if (layer>=6) {
    printf("ERROR: UpgradeClusterFinder::GetClusterWidthPhi: Module out of bounds: layer = %d\n",layer);
    return 0;
  }
  return fClusterList[layer].GetWidthPhiIndex(index);
}
//___________________________________________________________________________________
UInt_t AliITSUpgradeClusterFinder::GetClusterType(Int_t layer, UInt_t index) {
  if (layer>=6) {
    printf("ERROR: UpgradeClusterFinder::GetClusterType: Module out of bounds: layer = %d\n",layer);
    return 0;
  }
  return fClusterList[layer].GetTypeIndex(index);
}
//___________________________________________________________________________________
void AliITSUpgradeClusterFinder::PrintAllClusters() {
  for (Int_t layer=0; layer<6; layer++) {
    PrintClusters(layer);
  }
}
//___________________________________________________________________________________
void AliITSUpgradeClusterFinder::PrintClusters(Int_t layer) {
  if (layer>=6) {
    printf("ERROR: UpgradeClusterFinder::PrintClusters: Out of bounds: layer = %d\n",layer);
    return;
  }
  if(fClusterList[layer].GetNrEntries()==0) {
    printf("no cluster list entries. Exiting... \n");
    return;
  }
  for (UInt_t c=0; c<fClusterList[layer].GetNrEntries(); c++) {
  printf("layer  %d , z,y=%f,%f , size=%d , type=%d labels=%d %d %d (label printout to be implemented...)\n",layer,fClusterList[layer].GetColIndex(c),fClusterList[layer].GetRowIndex(c),fClusterList[layer].GetSizeIndex(c),fClusterList[layer].GetTypeIndex(c),999,999,999);  
}  
}
//___________________________________________________________________________________
void AliITSUpgradeClusterFinder::NewEvent() {
  for (Int_t i=0; i<6; i++) {
    fClusterList[i].Clear();
  }
  NewModule();
  fOldModule = -1;
}
//___________________________________________________________________________________
void AliITSUpgradeClusterFinder::NewModule() {
  fNhitsLeft=0;
  memset(fHits,0,999*999*sizeof(Bool_t));
}
//___________________________________________________________________________________
Int_t AliITSUpgradeClusterFinder::DoModuleClustering(Int_t Layer, UShort_t charge) {
  UInt_t maxHits=fNhitsLeft;
  for (UInt_t hit=0; hit<maxHits; hit++) {
    if (fClusterTypeFlag) fFindClusterType = kTRUE;

    fClusterWidthMinCol = fHitCol[hit];
    fClusterWidthMinRow = fHitRow[hit];
    fClusterWidthMaxCol = fHitCol[hit];
    fClusterWidthMaxRow = fHitRow[hit];
    
 
    fClusterTypeOrigCol = fHitCol[hit];
    fClusterTypeOrigRow = fHitRow[hit];
    memset(fClusterTypeArea,0,kMAXCLUSTERTYPESIDEZ*kMAXCLUSTERTYPESIDEY*sizeof(Bool_t));
    fColSum=0;
    fRowSum=0;
    UInt_t size = FindClusterRecu(fClusterTypeOrigCol,fClusterTypeOrigRow,charge);

    if(size==1){
      fCharge=GetPixelCharge(fColSum,fRowSum);
      AddLabelIndex(fColSum,fRowSum);
    }
    if (size>0) {
     if(size==2) printf("DoModuleClustering, size 2, labels :  %i  %i  %i \n",fLabels[0],fLabels[1],fLabels[2]);
      fClusterList[Layer].Insert((Float_t)fColSum/size, (Float_t)fRowSum/size, size, GetClusterWidthZ(), GetClusterWidthPhi(), GetClusterType(size), fCharge,fLabels);
      fCharge=0;
      for(Int_t i=0; i<10; i++) fLabels[i]=-5;
    }
    if (fNhitsLeft==0) break;
  }
  return 0;
}
//___________________________________________________________________________________
UInt_t AliITSUpgradeClusterFinder::FindClusterRecu(Int_t col, Int_t row, UShort_t charge) {
  if (col<0 || !fHits[col][row]) {return 0;}
  fHits[col][row]=kFALSE;
  fColSum+=col;
  fRowSum+=row;
  fNhitsLeft--;
  UInt_t size=1;

  if (col>fClusterWidthMaxCol) fClusterWidthMaxCol = col;
  if (col<fClusterWidthMinCol) fClusterWidthMinCol = col;
  if (row>fClusterWidthMaxRow) fClusterWidthMaxRow = row;
  if (row<fClusterWidthMaxRow) fClusterWidthMinRow = row;

  if (fFindClusterType) {
    Short_t diffz = fClusterTypeOrigCol - col;
    Short_t diffy = row - fClusterTypeOrigRow;
    if (diffz>=kMAXCLUSTERTYPESIDEZ || diffy>=kMAXCLUSTERTYPESIDEY) {
      fFindClusterType=kFALSE;
    }
    else {
      if (diffz==-1) {
	ShiftClusterTypeArea(kSHIFTRIGHT);
	diffz=0;
      }
      if (diffy==-1) {
	ShiftClusterTypeArea(kSHIFTDOWN);
	diffy=0;
      }
      fClusterTypeArea[diffz][diffy] = kTRUE;
    }
  }
  // straight:
  size+=FindClusterRecu(col  ,row-1,charge);
  fCharge+=GetPixelCharge(col,row-1);
  AddLabelIndex(col,row-1);

  size+=FindClusterRecu(col-1,row  ,charge);
  fCharge+=GetPixelCharge(col-1,row);
  AddLabelIndex(col-1,row);

  size+=FindClusterRecu(col  ,row+1,charge);
  fCharge+=GetPixelCharge(col,row+1);
  AddLabelIndex(col,row+1);

  size+=FindClusterRecu(col+1,row ,charge );
  fCharge+=GetPixelCharge(col+1,row);
  AddLabelIndex(col+1,row);

  // diagonal:
  size+=FindClusterRecu(col-1,row-1,charge);
  fCharge+=GetPixelCharge(col-1,row-1);
  AddLabelIndex(col-1,row-1);

  size+=FindClusterRecu(col-1,row+1,charge);
  fCharge+=GetPixelCharge(col-1,row+1);
  AddLabelIndex(col-1,row+1);

  size+=FindClusterRecu(col+1,row-1,charge);
  fCharge+=GetPixelCharge(col+1,row-1);
  AddLabelIndex(col+1,row-1);

  size+=FindClusterRecu(col+1,row+1,charge);
  fCharge+=GetPixelCharge(col+1,row+1);
  AddLabelIndex(col+1,row+1);

  return size;
}
//___________________________________________________________________________________
void AliITSUpgradeClusterFinder::ShiftClusterTypeArea(UInt_t direction) {
  if (direction == kSHIFTRIGHT) {
    fClusterTypeOrigCol++;
    Bool_t tmpClusterTypeArea[kMAXCLUSTERTYPESIDEZ][kMAXCLUSTERTYPESIDEY];
    memset(tmpClusterTypeArea,0,kMAXCLUSTERTYPESIDEZ*kMAXCLUSTERTYPESIDEY*sizeof(Bool_t));
    for (UInt_t z=0; z<kMAXCLUSTERTYPESIDEZ; z++) {
      for (UInt_t y=0; y<kMAXCLUSTERTYPESIDEY; y++) {
	if (fClusterTypeArea[z][y]) {
	  if (z==kMAXCLUSTERTYPESIDEZ-1) {
	    fFindClusterType=kFALSE;
	    return;
	  }
	  else {
	    tmpClusterTypeArea[z+1][y] = kTRUE;
	  }
	}
      }
    }
    memcpy(fClusterTypeArea,tmpClusterTypeArea,kMAXCLUSTERTYPESIDEZ*kMAXCLUSTERTYPESIDEY*sizeof(Bool_t));
  }
  else if (direction == kSHIFTDOWN) {
    fClusterTypeOrigRow--;
    Bool_t tmpClusterTypeArea[kMAXCLUSTERTYPESIDEZ][kMAXCLUSTERTYPESIDEY];
    memset(tmpClusterTypeArea,0,kMAXCLUSTERTYPESIDEZ*kMAXCLUSTERTYPESIDEY*sizeof(Bool_t));
    for (UInt_t z=0; z<kMAXCLUSTERTYPESIDEZ; z++) {
      for (UInt_t y=0; y<kMAXCLUSTERTYPESIDEY; y++) {
	if (fClusterTypeArea[z][y]) {
	  if (y==kMAXCLUSTERTYPESIDEY-1) {
	    fFindClusterType=kFALSE;
	    return;
	  }
	  else {
	    tmpClusterTypeArea[z][y+1] = kTRUE;
	  }
	}
      }
    }
    memcpy(fClusterTypeArea,tmpClusterTypeArea,kMAXCLUSTERTYPESIDEZ*kMAXCLUSTERTYPESIDEY*sizeof(Bool_t));
  }
}



//___________________________________________________________________________________
UShort_t AliITSUpgradeClusterFinder::GetCharge(Int_t layer,UInt_t index ) {
  return fClusterList[layer].GetCharge(index);
}
//___________________________________________________________________________________
Int_t * AliITSUpgradeClusterFinder::GetLabels(Int_t layer,UInt_t index) {
       return fClusterList[layer].GetLabels(index);
}

//___________________________________________________________________________________
UInt_t AliITSUpgradeClusterFinder::GetClusterWidthZ() {
  return fClusterWidthMaxCol-fClusterWidthMinCol+1;
}
//___________________________________________________________________________________
UInt_t AliITSUpgradeClusterFinder::GetClusterWidthPhi() {
  return fClusterWidthMaxRow-fClusterWidthMinRow+1;
}
//___________________________________________________________________________________
UInt_t AliITSUpgradeClusterFinder::GetClusterType(UInt_t size) {
  // category 'other':
  if (!fFindClusterType || size>4) 
    return 0;

  // X
  if (size==1) 
    return 1;

  // X
  // X
  if (size==2 && 
      fClusterTypeArea[0][1] && 
      fClusterTypeArea[0][0] ) 
    return 2;

  // XX
  if (size==2 && 
      fClusterTypeArea[1][0] &&
      fClusterTypeArea[0][0] ) 
    return 3;

  // X
  // X
  // X
  if (size==3 && 
      fClusterTypeArea[0][2] && 
      fClusterTypeArea[0][1] && 
      fClusterTypeArea[0][0] ) 
    return 4;
  
  // XX
  // X      and equivalent...
  if (size==3 && 
      (
       (
	fClusterTypeArea[0][0] && 
	fClusterTypeArea[0][1] && 
	(fClusterTypeArea[1][1] || 
	 fClusterTypeArea[1][0])
	)
       ||
       (
	fClusterTypeArea[1][0] && 
	fClusterTypeArea[1][1] && 
	(fClusterTypeArea[0][1] || 
	 fClusterTypeArea[0][0])
	)
       )
      ) 
    return 5;
  
  // XXX
  if (size==3 && 
      fClusterTypeArea[2][0] && 
      fClusterTypeArea[1][0] && 
      fClusterTypeArea[0][0] )
    return 6;
  
  // X
  // X
  // X
  // X
  if (size==4 &&
      fClusterTypeArea[0][3] && 
      fClusterTypeArea[0][2] &&
      fClusterTypeArea[0][1] && 
      fClusterTypeArea[0][0] ) 
    return 7;
  
  // XX
  // XX
  if (size==4 && 
      fClusterTypeArea[1][1] && 
      fClusterTypeArea[1][0] && 
      fClusterTypeArea[0][1] &&
      fClusterTypeArea[0][0] ) 
    return 8;

  // XX
  // X
  // X     and equivalent...
  if (size==4 &&
      (
       (
	fClusterTypeArea[1][0] && 
	fClusterTypeArea[1][1] && 
	fClusterTypeArea[1][2] && 
	fClusterTypeArea[0][0]
	)
       ||
       (
	fClusterTypeArea[1][0] && 
	fClusterTypeArea[1][1] && 
	fClusterTypeArea[1][2] && 
	fClusterTypeArea[0][2]
	)
       ||
       (
	fClusterTypeArea[0][0] && 
	fClusterTypeArea[0][1] && 
	fClusterTypeArea[0][2] && 
	fClusterTypeArea[1][0]
	)
       ||
       (
	fClusterTypeArea[0][0] && 
	fClusterTypeArea[0][1] && 
	fClusterTypeArea[0][2] && 
	fClusterTypeArea[1][2]
	)
       )
      )
    return 9;

  // X
  // XX
  // X     and equivalent...
  if (size==4 &&
      (
       (
	fClusterTypeArea[1][0] && 
	fClusterTypeArea[1][1] && 
	fClusterTypeArea[1][2] && 
	fClusterTypeArea[0][1]
	)
       ||
       (
	fClusterTypeArea[0][0] && 
	fClusterTypeArea[0][1] && 
	fClusterTypeArea[0][2] && 
	fClusterTypeArea[1][1]
	)
       )
      )
    return 10;

  // X
  // XXX    and equivalent...
  if (size==4 &&
      (
       (
	fClusterTypeArea[2][0] && 
	fClusterTypeArea[2][1] && 
	fClusterTypeArea[1][1] && 
	fClusterTypeArea[0][1]
	)
       ||
       (
	fClusterTypeArea[0][0] && 
	fClusterTypeArea[2][1] && 
	fClusterTypeArea[1][1] && 
	fClusterTypeArea[0][1]
	)
       ||
       (
	fClusterTypeArea[2][1] && 
	fClusterTypeArea[2][0] && 
	fClusterTypeArea[1][0] && 
	fClusterTypeArea[0][0]
	)
       ||
       (
	fClusterTypeArea[0][1] && 
	fClusterTypeArea[2][0] && 
	fClusterTypeArea[1][0] && 
	fClusterTypeArea[0][0]
	)
       )
      )
    return 11;

  //  X
  // XXX     and equivalent...
  if (size==4 &&
      (
       (
	fClusterTypeArea[1][0] && 
	fClusterTypeArea[2][1] && 
	fClusterTypeArea[1][1] && 
	fClusterTypeArea[0][1]
	)
       ||
       (
	fClusterTypeArea[1][1] && 
	fClusterTypeArea[2][0] && 
	fClusterTypeArea[1][0] && 
	fClusterTypeArea[0][0]
	)
       )
      )
    return 12;

  //  X
  // X      and equivalent...
  if (size==2 &&
      (
       (
	fClusterTypeArea[0][0] &&
	fClusterTypeArea[1][1] 
	)
       ||
       (
	fClusterTypeArea[1][0] &&
	fClusterTypeArea[0][1] 
	)
       )
      )
    return 13;

  //  X
  //  X
  // X      and equivalent...
  if (size==3 &&
      (
       (
	fClusterTypeArea[0][0] &&
	fClusterTypeArea[0][1] && 
	fClusterTypeArea[1][2] 
	)
       ||
       (
	fClusterTypeArea[1][0] &&
	fClusterTypeArea[1][1] && 
	fClusterTypeArea[0][2] 
	)
       ||
       (
	fClusterTypeArea[1][0] &&
	fClusterTypeArea[0][1] && 
	fClusterTypeArea[0][2] 
	)
       ||
       (
	fClusterTypeArea[0][0] &&
	fClusterTypeArea[1][1] && 
	fClusterTypeArea[1][2] 
	)
       )
      )
    return 14;

  //  XX
  // X      and equivalent...
  if (size==3 &&
      (
       (
	fClusterTypeArea[0][0] &&
	fClusterTypeArea[1][0] && 
	fClusterTypeArea[2][1] 
	)
       ||
       (
	fClusterTypeArea[0][1] &&
	fClusterTypeArea[1][1] && 
	fClusterTypeArea[2][0] 
	)
       ||
       (
	fClusterTypeArea[0][0] &&
	fClusterTypeArea[1][1] && 
	fClusterTypeArea[2][1] 
	)
       ||
       (
	fClusterTypeArea[0][1] &&
	fClusterTypeArea[1][0] && 
	fClusterTypeArea[2][0] 
	)
       )
      )
    return 15;



  return 0;
}

//___________________________________________________________________________________________________

UInt_t AliITSUpgradeClusterFinder::GetPixelCharge(UInt_t col, UInt_t row){
  Int_t q=0;
  for(Int_t entry =0; entry < fChargeArray->GetEntries(); entry++){
    TObjString *s = (TObjString*)fChargeArray->At(entry);
    TString name = s->GetString();
    if(!name.Contains(Form("%i %i",col,row)))
      continue;
    TObjArray *array = name.Tokenize(" ");
    TString charge = ((TObjString*)array->At(2))->String();
    TString rowS, colS;
    rowS = ((TObjString*)array->At(0))->String();
    colS = ((TObjString*)array->At(1))->String();
    delete array;
    q=charge.Atoi();
    return q;

  }
  return q;
}
//____________________________________________________

void AliITSUpgradeClusterFinder::AddLabelIndex(UInt_t col, UInt_t row){

 for(Int_t entry =0; entry < fChargeArray->GetEntries(); entry++){
    TObjString *s = (TObjString*)fChargeArray->At(entry);
    TString name = s->GetString();
 if(!name.Contains(Form("%i %i",col,row)))
  continue;
 TObjArray *array = name.Tokenize(" ");
 TString index[3];
 Int_t label[3];
  for(Int_t i=0; i<3; i++){
  index[i]=((TObjString*)array->At(3+i))->String();
  label[i]=index[i].Atoi();

  }
 SetLabels(label);
 delete array;
 }
}
//____________________________________________________

void AliITSUpgradeClusterFinder::SetLabels(Int_t label[3]){


 Bool_t isAssigned[3]={kFALSE,kFALSE,kFALSE};

 for(Int_t i=0; i<10; i++){
  for(Int_t k=0; k<3; k++){
   //printf(" fLabels[%i]=%i  label[%i]=%i \n",i,fLabels[i],k,label[k]);
    if(fLabels[i]!=label[k] && label[k]>-1 && fLabels[i]<0) {
      if(!isAssigned[k]) {
     //  printf("assignign...\n.");
    printf("Assignign  fLabels[%i]=%i  label[%i]=%i \n",i,fLabels[i],k,label[k]);
        fLabels[i+k]=label[k];
        isAssigned[k]=kTRUE;
      }
     }
   }
 }
}

