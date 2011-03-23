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
#include "AliITSRecPointU.h"
#include "AliITSDigitUpgrade.h"
#include "AliITSRawStreamSPD.h"
#include "AliITSUPixelModule.h"
#include "AliLog.h"
#include <string.h>
#include <TObjString.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <TMath.h>

AliITSUpgradeClusterFinder::AliITSUpgradeClusterFinder() : 
  TObject(),
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
  fClusterList(0x0),
  fChargeArray(0x0),
  fRecPoints(0x0),
  fNSectors()
{ 
  //
  // Default constructor
  //
  AliITSsegmentationUpgrade *s = new AliITSsegmentationUpgrade();
  fNSectors = s->GetNSectors();
  delete s;
  fClusterList = new AliITSUpgradeClusterList*[fNSectors];  
 for(Int_t imod =0; imod < fNSectors; imod++){
  fClusterList[imod] = new AliITSUpgradeClusterList();
} 

  fChargeArray = new TObjArray();
  fChargeArray->SetOwner(kTRUE);
  fRecPoints = new TClonesArray("AliITSRecPointU",3000);
  fTmpLabel[0]=-5;
  fTmpLabel[1]=-5;
  fTmpLabel[2]=-5;
  for(int il=0; il<10;il++) fLabels[il]=-5;
  for(int k=0; k<999999; k++){
    fHitCol[k]=999;
    fHitRow[k]=999;
  }
}

//___________________________________________________________________________________
AliITSUpgradeClusterFinder::~AliITSUpgradeClusterFinder() {
  if(fChargeArray) delete fChargeArray;
  if(fRecPoints) delete fRecPoints; 
  if(fClusterList)delete [] fClusterList;
}
//___________________________________________________________________________________
void AliITSUpgradeClusterFinder::StartEvent() {
  NewEvent();
}
//___________________________________________________________________________________
Int_t AliITSUpgradeClusterFinder::ProcessHit(Int_t module , UInt_t col, UInt_t row, UShort_t charge, Int_t label[3]) {
  //
  // Adds one pixel to the cluster 
  //
  AliDebug(2,Form("module,col,row,charge,label(0,1,2) = %d,%d,%d,%d,(%d,%d,%d)\n",module ,col,row,charge,label[0],label[1],label[2])); 
  if (module>=fNSectors) {
    AliError(Form("Out of bounds: module ,col,row, charge, label(0,1,2)= %d,%d,%d,%d,(%d,%d,%d)\n",module ,col,row,charge,label[0],label[1],label[2])); 
    return 1;
  }
  // do we need to perform clustering on previous module?
  if (fOldModule!=-1 && (Int_t)module!=fOldModule) {
    //fChargeArray->AddLast(new TObjString(Form("%i %i %i %i %i %i",col,row,charge,label[0],label[1],label[2])));
    AliITSUPixelModule *pix = new AliITSUPixelModule(module,col, row, charge, label);
    fChargeArray->AddLast(pix);
    
    DoModuleClustering(fOldModule,charge);
    NewModule();
  }
  // add hit
  AliITSUPixelModule *pix = new AliITSUPixelModule(module,col, row, charge, label);
    fChargeArray->AddLast(pix);
//  fChargeArray->AddLast(new TObjString(Form("%i %i %i %i %i %i",col,row,charge,label[0],label[1],label[2])));

  fOldModule=module;
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
UInt_t AliITSUpgradeClusterFinder::GetClusterCount(Int_t module) const {
  //
  // number of clusters in the module
  // 
  if (module>fNSectors ) {
    printf("ERROR: UpgradeClusterFinder::GetClusterCount: Module out of bounds: module = %d\n",module);
    return 0;
  }
  return fClusterList[module]->GetNrEntries();
}
//___________________________________________________________________________________
Float_t AliITSUpgradeClusterFinder::GetClusterMeanCol(Int_t module, UInt_t index) {
  //
  // cluster position in terms of colums : mean column ID
  //
  if (module>=fNSectors) {
    printf("ERROR: UpgradeClusterFinder::GetClusterMeanCol: Module out of bounds: module = %d\n",module);
    return 0;
  }
  return fClusterList[module]->GetColIndex(index);
}
//___________________________________________________________________________________
Float_t AliITSUpgradeClusterFinder::GetClusterMeanRow(Int_t module, UInt_t index) {
  //
  // cluster position in terms of rows : mean row ID
  //
  if (module>=fNSectors) {
    printf("ERROR: UpgradeClusterFinder::GetClusterMeanRow: Module out of bounds: module = %d\n",module);
    return 0;
  }
  return fClusterList[module]->GetRowIndex(index);
}
//___________________________________________________________________________________
UInt_t AliITSUpgradeClusterFinder::GetClusterSize(Int_t module, UInt_t index) {
  //
  // total number of pixels of the cluster 
  //
  if (module>=fNSectors) {
    printf("ERROR: UpgradeClusterFinder::GetClusterSize: Module out of bounds: module = %d\n",module);
    return 0;
  }
  return fClusterList[module]->GetSizeIndex(index);
}
//___________________________________________________________________________________
UInt_t AliITSUpgradeClusterFinder::GetClusterWidthZ(Int_t module, UInt_t index) {
  //
  // # pixels of the cluster in Z direction
  //
  
  if (module>=fNSectors) {
    printf("ERROR: UpgradeClusterFinder::GetClusterWidthZ: Module out of bounds: module = %d\n",module);
    return 0;
  }
  return fClusterList[module]->GetWidthZIndex(index);
}
//___________________________________________________________________________________
UInt_t AliITSUpgradeClusterFinder::GetClusterWidthPhi(Int_t module, UInt_t index) {
  //
  // # pixels of the cluster in phi direction (XY plane)
  //

  if (module>=fNSectors) {
    printf("ERROR: UpgradeClusterFinder::GetClusterWidthPhi: Module out of bounds: module = %d\n",module);
    return 0;
  }
  return fClusterList[module]->GetWidthPhiIndex(index);
}
//___________________________________________________________________________________
UInt_t AliITSUpgradeClusterFinder::GetClusterType(Int_t module, UInt_t index) {
  //
  // cluster shape
  //

  if (module>=fNSectors) {
    printf("ERROR: UpgradeClusterFinder::GetClusterType: Module out of bounds: layer = %d\n",module);
    return 0;
  }
  return fClusterList[module]->GetTypeIndex(index);
}
//___________________________________________________________________________________
void AliITSUpgradeClusterFinder::PrintAllClusters() {
  //
  // printout of the cluster information
  // 
  
  for (Int_t module=0; module<fNSectors; module++) {
    PrintClusters(module);
  }
}
//___________________________________________________________________________________
void AliITSUpgradeClusterFinder::PrintClusters(Int_t module) {
  //
  // printout of the cluster information
  //   
  
  if (module>=fNSectors) {
    printf("ERROR: UpgradeClusterFinder::PrintClusters: Out of bounds: layer = %d\n",module);
    return;
  }
  if(fClusterList[module]->GetNrEntries()==0) {
    printf("no cluster list entries. Exiting... \n");
    return;
  }
 // for (UInt_t c=0; c<fClusterList[layer]->GetNrEntries(); c++) {
    //printf("layer  %d , z,y=%f,%f , size=%d , type=%d labels=%d %d %d (label printout to be implemented...)\n",layer,fClusterList[layer].GetColIndex(c),fClusterList[layer].GetRowIndex(c),fClusterList[layer].GetSizeIndex(c),fClusterList[layer].GetTypeIndex(c),999,999,999);  
 // }  
}
//___________________________________________________________________________________
void AliITSUpgradeClusterFinder::NewEvent() {
  //
  // Cleaning up and preparation for the clustering procedure
  //
  
  for (Int_t i=0; i<fNSectors; i++) {
   fClusterList[i]->Clear();
 }
  //NewModule();
  fOldModule = -1;
}
//___________________________________________________________________________________
void AliITSUpgradeClusterFinder::NewModule() {
  //
  // Initializations
  //
  fNhitsLeft=0;
  memset(fHits,0,999*999*sizeof(Bool_t));
  fChargeArray->Clear();
}
//___________________________________________________________________________________
Int_t AliITSUpgradeClusterFinder::DoModuleClustering(Int_t module, UShort_t charge) {
  //
  // Clustering and cluster-list container filling
  //
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
      if(size>1) AliDebug(2,Form("DoModuleClustering, size %i , labels :  %i  %i  %i \n",size,fLabels[0],fLabels[1],fLabels[2]));
      if(size > kMAXCLUSTERTYPESIDEZ*kMAXCLUSTERTYPESIDEY ) return 0; // such clusters are rejected. Temporary fix. not clear why....`
      fClusterList[module]->Insert((Float_t)fColSum/size, (Float_t)fRowSum/size, size, GetClusterWidthZ(), GetClusterWidthPhi(), GetClusterType(size), fCharge,fLabels);
      fCharge=0;
      for(Int_t i=0; i<10; i++) fLabels[i]=-5;
    }
    if (fNhitsLeft==0) break;
  }
  return 0;
}
//___________________________________________________________________________________
UInt_t AliITSUpgradeClusterFinder::FindClusterRecu(Int_t col, Int_t row, UShort_t charge) {
  //
  // Pixel selection for the clusters (adjiacent pixels or pixels at the corners) 
  //
  if (col<0 || !fHits[col][row]) {
    return 0;
  }
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
  //
  // Cluster checks 
  //
  
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
UShort_t AliITSUpgradeClusterFinder::GetCharge(Int_t module,UInt_t index ) {
  return fClusterList[module]->GetCharge(index);
}
//___________________________________________________________________________________
Int_t * AliITSUpgradeClusterFinder::GetLabels(Int_t module,UInt_t index) {
  return fClusterList[module]->GetLabels(index);
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
  //AliInfo(Form(" entrate charge array %i ", fChargeArray->GetEntries()));
  for(Int_t entry =0; entry < fChargeArray->GetEntries(); entry++){
/*    TObjString *s = (TObjString*)fChargeArray->At(entry);
    TString name = s->GetString();
    if(!name.Contains(Form("%i %i",col,row)))
      continue;
    AliInfo(Form(" 1 entry %i ", entry));
    TObjArray *array = name.Tokenize(" ");
    array->SetOwner(kTRUE);
    AliInfo(Form(" 2 entry %i ", entry));
    TString charge = ((TObjString*)array->At(2))->String();
    
    TString rowS, colS;
    rowS = ((TObjString*)array->At(0))->String();
    colS = ((TObjString*)array->At(1))->String();
    AliInfo(Form(" 3 prima del delete entry %i ", entry));
    array->Clear();
    delete array;
    AliInfo(Form(" 4 dopo il delete  entry %i ", entry));
    q=charge.Atoi();
*/
    AliITSUPixelModule *pixMod = (AliITSUPixelModule*)fChargeArray->At(entry);
  //  pixMod->PrintInfo();
    if(col!=pixMod->GetCol())continue;
    if(row!=pixMod->GetRow())continue;
    q= pixMod->GetCharge();

      }
  return q;
}
//____________________________________________________

void AliITSUpgradeClusterFinder::AddLabelIndex(UInt_t col, UInt_t row){
  //
  // Adding cluster labels
  //

  for(Int_t entry =0; entry < fChargeArray->GetEntries(); entry++){
/*    TObjString *s = (TObjString*)fChargeArray->At(entry);
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
*/ 
   AliITSUPixelModule *pix= (AliITSUPixelModule*)fChargeArray->At(entry);
   if(col!=pix->GetCol())continue;
   if(row!=pix->GetRow())continue;
    Int_t label[3]={-1,-1,-1};
    for(Int_t i=0;i<3;i++){
    label[i] = pix->GetLabels(i);
    }   
    SetLabels(label);
  }
}
//____________________________________________________

void AliITSUpgradeClusterFinder::SetLabels(Int_t label[3]){
  
  //Set the MC lables to the cluster

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

  if(!fRecPoints)fRecPoints = new TClonesArray("AliITSRecPointU",1000);
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
  if(!fRecPoints) fRecPoints = new TClonesArray("AliITSRecPointU",1000);

  TBranch *branch;
  branch = treeR->GetBranch("ITSRecPoints");
  if (branch) {
    branch->SetAddress(&fRecPoints);
  } else AliError("No RecPoint branch available");

}
//____________________________________________________
void AliITSUpgradeClusterFinder::DigitsToRecPoints(const TObjArray *digList) {
  //
  // the clusterization is performed here
  //
  AliITSsegmentationUpgrade *segmentation = new AliITSsegmentationUpgrade(); 
  AliITSRecPointU  recpnt;
  Int_t nClusters =0;
  TClonesArray &lrecp = *fRecPoints;
  for(Int_t ilayer=0; ilayer < 6 ;ilayer ++){
    NewModule();
    TClonesArray *pArrDig= (TClonesArray*)digList->At(ilayer);
    StartEvent(); 
    pArrDig->Sort();
    
    for(Int_t ientr =0; ientr < pArrDig->GetEntries() ; ientr++){
      AliITSDigitUpgrade *dig = (AliITSDigitUpgrade*)pArrDig->At(ientr);
      Int_t colz=dig->GetzPixelNumber();
      Int_t rowx=dig->GetxPixelNumber();
      Double_t hitcharge= (dig->GetNelectrons());
      ProcessHit(dig->GetModule(),colz, rowx,(Short_t)hitcharge,dig->GetTracks());
    }//ientr
    FinishEvent();
    for(Int_t module=0; module<fNSectors; module++){

  //  printf(" module loop %i \n", module); 
    for(UInt_t nClu = 0; nClu <  GetClusterCount(module); nClu++){
      //printf(" nclu in getclustercount %i \n", nClu);
      UShort_t charge = GetCharge(module, nClu);
      recpnt.SetQ(charge);
      recpnt.SetLayer(ilayer);
      recpnt.SetModule(module);
      Int_t *lab=GetLabels(module,nClu);
      for(Int_t l=0; l<3; l++) {recpnt.SetLabel(lab[l],l);}

      //Bool_t check2;
      //Double_t xcheck2=0.;
      //Double_t ycheck2=0.;
      //Double_t zcheck2=0.;
      Double_t xzl2[2]={0.,0.};
      Double_t xPixC2,zPixC2=0.;

      xPixC2 = GetClusterMeanRow(module, nClu);
      zPixC2 = GetClusterMeanCol(module, nClu);
      //cout << "zPixC2 "<< zPixC2 << endl;
      xzl2[0] = xPixC2*(segmentation->GetCellSizeX(ilayer))+0.5*(segmentation-> GetCellSizeX(ilayer));
      xzl2[1] = zPixC2*(segmentation->GetCellSizeZ(ilayer))+0.5*(segmentation->GetCellSizeZ(ilayer))-(segmentation->GetHalfLength(ilayer)); 
      //cout << " vediamo che positione ha il recpoint !!!! zl = "<< xzl2[1] << endl;
      //check2 = segmentation->DetToGlobal(ilayer,xzl2[0], xzl2[1],xcheck2,ycheck2,zcheck2);
      recpnt.SetType(GetClusterType(module,nClu ));
      recpnt.SetLocalCoord(xzl2[0],xzl2[1]); //temporary solution (no LocalToTrack Matrix)
      //from global to tracking system coordinate
      // global coordinate -> local coordinate getting alpha angle of the recpoint
    /////
      Double_t yclu1 = 0.;//upgrade clusters global coordinate ( ITS official: GetX tracking coordinate)
      Double_t zclu1 = 0.;//upgrade clusters global coordinate ( ITS official: GetX tracking coordinate)
     // Float_t xclg = xcheck2;//upgrade clusters global coordinate ( ITS official: GetX tracking coordinate)
     // Float_t yclg = ycheck2;
     // Float_t zclg = zcheck2;
      Bool_t detr=kFALSE;
      detr = segmentation->DetToTrack(ilayer,module, xzl2[0],xzl2[1], yclu1, zclu1);      
//      printf( " det to track in clusterfinder il %i xl %f zl %f y track %f z track %f module %i \n", ilayer, xzl2[0] , xzl2[1] , yclu1, zclu1, module);      

//////////////////////////
      recpnt.SetX(0.);
      recpnt.SetY(yclu1);
      recpnt.SetZ(zclu1);
      
      
      Double_t xsize, zsize;
      segmentation->GetSegmentation(ilayer,xsize, zsize);
      recpnt.SetSigmaY2(xsize/TMath::Sqrt(12)*xsize/TMath::Sqrt(12));
      recpnt.SetSigmaZ2(zsize/TMath::Sqrt(12)*zsize/TMath::Sqrt(12));
      new(lrecp[nClusters++]) AliITSRecPointU(recpnt);
      //Int_t idx = fRecPoints->GetEntries();
      AliDebug(1,Form("recpoint : Nelectrons %f (entry %i)",recpnt.GetQ(),fRecPoints->GetEntries()));
    
  }//cluster list entries
  }//module
  }//ilayer
  if(segmentation) delete segmentation;
}

