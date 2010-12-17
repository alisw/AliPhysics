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
#include "AliITSsegmentationUpgrade.h"
#include "AliITSRecPoint.h"
#include "AliITSDigitUpgrade.h"
#include "AliITSRawStreamSPD.h"
#include "AliLog.h"
#include <string.h>
#include <TObjString.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <TMath.h>

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
  fChargeArray(0x0),
  fRecPoints(0x0)
{ 
  fChargeArray = new TObjArray();
  fRecPoints = new TClonesArray("AliITSRecPoint",3000);
  fTmpLabel[0]=-5;
  fTmpLabel[1]=-5;
  fTmpLabel[2]=-5;
  for(int il=0; il<10;il++) fLabels[il]=-5;
}

//___________________________________________________________________________________
AliITSUpgradeClusterFinder::~AliITSUpgradeClusterFinder() {
  if(fChargeArray) delete fChargeArray;
  if(fRecPoints) delete fRecPoints; 
}
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
  //
  // Cluster shape
  //
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
  //
  //...self explaining
  //
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
  //
  // Adding cluster labels
  //

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
  //
  //
  //

  Bool_t isAssigned[3]={kFALSE,kFALSE,kFALSE};

  for(Int_t i=0; i<10; i++){
    for(Int_t k=0; k<3; k++){
      if(fLabels[i]!=label[k] && label[k]>-1 && fLabels[i]<0) {
	if(!isAssigned[k]) {
	  //  printf("assignign...\n.");
	  //  printf("Assignign  fLabels[%i]=%i  label[%i]=%i \n",i,fLabels[i],k,label[k]);
	  fLabels[i+k]=label[k];
	  isAssigned[k]=kTRUE;
	}
      }
    }
  }
}

//____________________________________________________
void AliITSUpgradeClusterFinder::MakeRecPointBranch(TTree *treeR){
  //
  // Creating the branch (see AliITSUpgradeReconstructor::Reconstruct)
  //

  if(!fRecPoints)fRecPoints = new TClonesArray("AliITSRecPoint",1000);
  if (treeR) {
    TBranch *branch = treeR->GetBranch("ITSRecPoints");
    if (branch) return ;
    else branch = treeR->Branch("ITSRecPoints",&fRecPoints); 
  }
}
//____________________________________________________
void AliITSUpgradeClusterFinder::SetRecPointTreeAddress(TTree *treeR){
  //
  // Addressing the branch (see AliITSUpgradeReconstructor::Reconstruct)
  //
  if(!treeR) return;
  if(!fRecPoints) fRecPoints = new TClonesArray("AliITSRecPoint",1000);

  TBranch *branch;
  branch = treeR->GetBranch("ITSRecPoints");
  if (branch) {
    branch->SetAddress(&fRecPoints);
  } else AliError("No RecPoint branch available");

}
//____________________________________________________
void AliITSUpgradeClusterFinder::DigitsToRecPoints(TObjArray *digList) {
  //
  // the clusterization is performed here
  //
  AliITSsegmentationUpgrade *segmentation2 = 0x0;
  AliITSRecPoint  recpnt;
  if(!segmentation2) segmentation2 = new AliITSsegmentationUpgrade();
  Int_t nClusters =0;
  TClonesArray &lrecp = *fRecPoints;

  for(Int_t ilayer=0; ilayer < 6 ;ilayer ++){
    TClonesArray *pArrDig= (TClonesArray*)digList->At(ilayer);
    StartEvent();
    AliInfo(Form("layer %i with digit entries %i",ilayer,pArrDig->GetEntries()));
    for(Int_t ientr =0; ientr < pArrDig->GetEntries() ; ientr++){
      AliITSDigitUpgrade *dig = (AliITSDigitUpgrade*)pArrDig->At(ientr);
      Int_t colz=dig->GetzPixelNumber();
      Int_t rowx=dig->GetxPixelNumber();
      Double_t hitcharge= (dig->GetNelectrons());
      ProcessHitOnline(ilayer,colz, rowx,(Short_t)hitcharge,dig->GetTracks());
    }//ientr
    FinishEvent();
    
    for(UInt_t nClu = 0; nClu <  GetClusterCount(ilayer); nClu++){
      UShort_t charge = GetCharge(ilayer, nClu);
      recpnt.SetQ(charge);
      recpnt.SetLayer(ilayer);
      Int_t *lab=GetLabels(ilayer,nClu);
      for(Int_t l=0; l<3; l++) {recpnt.SetLabel(lab[l],l);}

      Bool_t check2;
      Double_t xcheck2=0.;
      Double_t ycheck2=0.;
      Double_t zcheck2=0.;
      Double_t xzl2[2]={0.,0.};
      Double_t XpixC2,ZpixC2=0.;

      XpixC2 = GetClusterMeanRow(ilayer, nClu);
      ZpixC2 = GetClusterMeanCol(ilayer, nClu);
      xzl2[0] = XpixC2*(segmentation2->GetCellSizeX(ilayer))+0.5*(segmentation2-> GetCellSizeX(ilayer));
      xzl2[1] = ZpixC2*(segmentation2->GetCellSizeZ(ilayer))+0.5*(segmentation2->GetCellSizeZ(ilayer))-(segmentation2->GetHalfLength(ilayer));
      check2 = segmentation2->DetToGlobal(ilayer,xzl2[0], xzl2[1],xcheck2,ycheck2,zcheck2);
      recpnt.SetType(GetClusterType(ilayer,nClu ));
      // recpnt.SetLocalCoord(xzl2[0],xzl2[1]); //temporary solution (no LocalToTrack Matrix)
      //from global to tracking system coordinate
      // global coordinate -> local coordinate getting alpha angle of the recpoint
      Float_t xclg = xcheck2;//upgrade clusters global coordinate ( ITS official: GetX tracking coordinate)
      Float_t yclg = ycheck2;
      Float_t zclg = zcheck2;
      Double_t phiclu1rad, phiclu1deg;
      phiclu1rad=TMath::ATan2(yclg,xclg);//cluster phi angle (rad)
      if (phiclu1rad<0) phiclu1rad+=TMath::TwoPi();//from 0 to 360
      else if (phiclu1rad>=TMath::TwoPi()) phiclu1rad-=TMath::TwoPi();//

      phiclu1deg=180.*phiclu1rad/TMath::Pi();// in deg
      Int_t ladder;// virtual segmentation starting from the cluster phi

      ladder=(Int_t)phiclu1deg/18;// in which ladder the cluster is
      Double_t alpha= (ladder*18.+9.)*TMath::Pi()/180.;//angle at the center of the ladder (rad)

      //alpha rotation
      Float_t xclu1_tr = xclg*TMath::Cos(alpha)-yclg*TMath::Sin(alpha);
      Float_t yclu1 = yclg*TMath::Cos(alpha)+xclg*TMath::Sin(alpha);
      Float_t xclu1 = TMath::Sqrt(xclu1_tr*xclu1_tr+yclu1*yclu1);
      Float_t zclu1 = zclg;
      Double_t phi_trk= (phiclu1rad-alpha);// cluster angle in the rotated system (rad)

      yclu1=xclu1*phi_trk; // tracking system coordinate: r*phi
      recpnt.SetX(0.);
      recpnt.SetY(yclu1);
      recpnt.SetZ(zclu1);

      Double_t xsize, zsize;
      segmentation2->GetSegmentation(ilayer,xsize, zsize);
      recpnt.SetSigmaY2(xsize/TMath::Sqrt(12)*xsize/TMath::Sqrt(12));
      recpnt.SetSigmaZ2(zsize/TMath::Sqrt(12)*zsize/TMath::Sqrt(12));
      new(lrecp[nClusters++]) AliITSRecPoint(recpnt);
      //Int_t idx = fRecPoints->GetEntries();
      AliInfo(Form("recpoint : Nelectrons %f (entry %i)",recpnt.GetQ(),fRecPoints->GetEntries()));
    }//cluster list entries
  }//ilayer
}

