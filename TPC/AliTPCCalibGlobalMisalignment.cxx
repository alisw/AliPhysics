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


////////////////////////////////////////////////////////////////////////////
//                                                                        //
// AliTPCCalibGlobalMisalignment class                                        //
// The class calculates the space point distortions due to simple         // 
// misalignments like shifts in caresian coordinates or a rotation        //
// of the TPC read out planes (A and C side)                              //
// Optionaly possible to use it for visualization of the alignemnt form the Alignment OCDB //
// fUseGeoManager has to be set to kTRUE to enable this option
//                                                                        //
// date: 06/05/2010                                                       //
// Authors: Stefan Rossegger, Jim Thomas, Magnus Mager                    //
////////////////////////////////////////////////////////////////////////////

#include "AliTPCCalibGlobalMisalignment.h"
#include "TMath.h"
#include "TGeoMatrix.h"
#include "AliTPCROC.h"
#include "AliTPCcalibDB.h"
#include "AliTPCParam.h"
#include <TGeoPhysicalNode.h>
//
#include "AliAlignObjParams.h"
#include "AliGeomManager.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"

AliTPCCalibGlobalMisalignment::AliTPCCalibGlobalMisalignment()
  : AliTPCCorrection("mialign","Misalignment"),
    fXShift(0.),fYShift(0.),fZShift(0.),
    fRotPhiA(0.),fRotPhiC(0.),
    fdRPhiOffsetA(0.), 
    fdRPhiOffsetC(0.), 
    fQuadrantQ0(0),   //OROC medium pads -delta ly+ - ly - shift
    fQuadrantRQ0(0),  //OROC medium pads -delta ly+ - ly - rotation 
    fQuadrantQ1(0),   //OROC long   pads -delta ly+ - ly - shift
    fQuadrantQ2(0),   //OROC long   pads -shift
    fQuadrantRQ1(0),  //OROC long   pads -delta ly+ - ly - rotation
    fQuadrantRQ2(0),  //OROC long   pads -rotation
    fMatrixGlobal(0), // global Alignment common
    fMatrixGlobalDelta(0), // global Alignment Delta A side-c side
    fArraySector(0)   // fArraySector
{
  //
  // default constructor
  //
}

AliTPCCalibGlobalMisalignment::~AliTPCCalibGlobalMisalignment() {
  //
  // destructor
  //
  delete fQuadrantQ0;   //OROC medium pads -delta ly+ - ly - shift
  delete fQuadrantRQ0;  //OROC medium pads -delta ly+ - ly - rotation 
  delete fQuadrantQ1;   //OROC long   pads -delta ly+ - ly - shift
  delete fQuadrantQ2;   //OROC long   pads -shift
  delete fQuadrantRQ1;  //OROC long   pads -delta ly+ - ly - rotation
  delete fQuadrantRQ2;  //OROC long   pads -rotation
  delete fMatrixGlobal; // global matrix
  delete fMatrixGlobal; // global matrix
  delete fArraySector;  // sector matrices
}
 
void AliTPCCalibGlobalMisalignment::SetQuadranAlign(const TVectorD *quadrantQ0, const TVectorD *quadrantRQ0, const TVectorD *quadrantQ1,const TVectorD *quadrantRQ1,  const TVectorD *quadrantQ2,  const TVectorD *quadrantRQ2){
  //
  // Set quadrant alignment
  // 6 vectors for 36 (super) sectors
  //
  if (quadrantQ0) fQuadrantQ0   = new TVectorD(*quadrantQ0);
  if (quadrantRQ0) fQuadrantRQ0 = new TVectorD(*quadrantRQ0);
  //
  if (quadrantQ1) fQuadrantQ1   = new TVectorD(*quadrantQ1);
  if (quadrantQ1) fQuadrantRQ1  = new TVectorD(*quadrantRQ1);
  if (quadrantQ2) fQuadrantQ2   = new TVectorD(*quadrantQ2);
  if (quadrantQ2) fQuadrantRQ2  = new TVectorD(*quadrantRQ2);
}



void AliTPCCalibGlobalMisalignment::SetAlignGlobal(const TGeoMatrix * matrixGlobal){
  //
  // Set global misalignment
  // Object is OWNER 
  // 
  if (fMatrixGlobal) delete fMatrixGlobal;
  fMatrixGlobal=0;
  if (matrixGlobal) fMatrixGlobal = new TGeoHMatrix(*matrixGlobal);
}

void AliTPCCalibGlobalMisalignment::SetAlignGlobalDelta(const TGeoMatrix * matrixGlobalDelta){
  //
  // Set global misalignment
  // Object is OWNER 
  // 
  if (fMatrixGlobalDelta) delete fMatrixGlobalDelta;
  fMatrixGlobalDelta=0;
  if (matrixGlobalDelta) fMatrixGlobalDelta = new TGeoHMatrix(*matrixGlobalDelta);
}

void AliTPCCalibGlobalMisalignment::SetAlignSectors(const TObjArray *arraySector){
  //
  // Set misalignment TObjArray of TGeoMatrices  - for each sector
  // Object is OWNER
  // 
  if (fArraySector) delete fArraySector;
  fArraySector=0;
  if (arraySector) fArraySector = (TObjArray*)arraySector->Clone();
}


//void AliTPCCalibGlobalMisalignment::Init() {
//  //
// // Initialization funtion
//  //

//  // nothing to be initialized, results of this calibration class will go to the global aligment structure

//}

//void AliTPCCalibGlobalMisalignment::Update(const TTimeStamp &/*timeStamp*/) {
//  //
//  // Update function 
//  //
//
//  // nothing to be updated, results of this calibration class will go to the global aligment structure
//
//}



void AliTPCCalibGlobalMisalignment::GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]) {
  //
  // Calculates the simple correction due to a shift (in x,y,z) or an rotation of the TPC (around z)
  //  
  static AliTPCROC *tpcRoc =AliTPCROC::Instance();  
  Double_t xref  = ( tpcRoc->GetPadRowRadii(0,0)+tpcRoc->GetPadRowRadii(36,tpcRoc->GetNRows(36)-1))*0.5;
  Double_t xquadrant  = tpcRoc->GetPadRowRadii(36,53); // row 53 from uli
  Double_t xIO  = ( tpcRoc->GetPadRowRadii(0,tpcRoc->GetNRows(0)-1)+tpcRoc->GetPadRowRadii(36,0))*0.5;
  Double_t r=0, phi=0;
  r   = TMath::Sqrt( x[0]*x[0] + x[1]*x[1] );
  phi = TMath::ATan2(x[1],x[0]);
  // Getsector number
  Double_t sec=TMath::Nint(-0.5+(phi*9./TMath::Pi()));
  if (sec<0) sec+=18;
  Int_t isec = TMath::Nint(sec);    
  if  (roc%36>=18) isec+=18;
  //
  // Get the point on the local coordiante frame
  //
  Double_t alpha=(sec+0.5)*TMath::Pi()/9;
  Double_t pos[3]={0,0,x[2]};
  pos[0]=  TMath::Cos(alpha)*x[0]+TMath::Sin(alpha)*x[1];
  pos[1]= -TMath::Sin(alpha)*x[0]+TMath::Cos(alpha)*x[1];
  if (pos[0]>tpcRoc->GetPadRowRadiiUp(0))  isec+=36;

  //
  // apply quadrant alignment if available - in local coordinate frame
  //
  //
  Double_t posQG[3]={x[0],x[1],x[2]};
  {
    Double_t dly=0;
    Bool_t isQ0 = (pos[0]<xquadrant)&&(pos[0]>xIO);
    Bool_t isQ1 = (pos[0]>xquadrant);
    Double_t  sign = (pos[1]>0)? 1.: -1.;
    if (isQ0){
      if (fQuadrantQ0)  dly+=sign*(*fQuadrantQ0)[isec%36];  // shift in cm
      if (fQuadrantRQ0) dly+=sign*(*fQuadrantRQ0)[isec%36]*(pos[0]-xref);      
    }
    if (isQ1){
      if (fQuadrantQ1)  dly+=sign*(*fQuadrantQ1)[isec%36];  // shift in cm
      if (fQuadrantRQ1) dly+=sign*(*fQuadrantRQ1)[isec%36]*(pos[0]-xref);      
      if (fQuadrantQ2)  dly+=(*fQuadrantQ2)[isec%36];  // shift in cm
      if (fQuadrantRQ2) dly+=(*fQuadrantRQ2)[isec%36]*(pos[0]-xref);      
    }
    // Tranform the corrected point to the global frame
    posQG[0]=  TMath::Cos(alpha)*pos[0]-TMath::Sin(alpha)*(pos[1]+dly);
    posQG[1]=  TMath::Sin(alpha)*pos[0]+TMath::Cos(alpha)*(pos[1]+dly);
  }
  //
  // rotation of the read-out planes
  if  (roc%36<18) // A side
    phi += fRotPhiA;
  else         // C side
    phi += fRotPhiC;
  
  // Simply adding a constant dRPHi residual. PURELY FOR CALIBRATION PURPOSES
  if  (roc%36<18) // A side
    phi += fdRPhiOffsetA/r;
  else         // C side
    phi += fdRPhiOffsetC/r;
  
  dx[0] = r * TMath::Cos(phi) - x[0];
  dx[1] = r * TMath::Sin(phi) - x[1]; 
  dx[2] = 0.; 

  // Simple shifts
  dx[0] -= fXShift;
  dx[1] -= fYShift;
  dx[2] -= fZShift;
  // quadrant shifts
  dx[0] += (posQG[0]-x[0]);
  dx[1] += (posQG[1]-x[1]);
  //
  //
  if (fMatrixGlobal){
    // apply global alignment matrix
    Double_t ppos[3]={x[0],x[1],x[2]};
    Double_t pposC[3]={x[0],x[1],x[2]};
    fMatrixGlobal->LocalToMaster(ppos,pposC);
    dx[0]+=pposC[0]-ppos[0];
    dx[1]+=pposC[1]-ppos[1];
    dx[2]+=pposC[2]-ppos[2];
  }
  if (fMatrixGlobalDelta){
    // apply global alignment matrix A-C Side side
    Double_t ppos[3]={x[0],x[1],x[2]};
    Double_t pposC[3]={x[0],x[1],x[2]};
    fMatrixGlobalDelta->LocalToMaster(ppos,pposC);
    Double_t ssign=(roc%36<18) ? 1.:-1.;
    dx[0]+=ssign*(pposC[0]-ppos[0]);
    dx[1]+=ssign*(pposC[1]-ppos[1]);
    dx[2]+=ssign*(pposC[2]-ppos[2]);
  }

  if (fArraySector){
    // apply global alignment matrix
    TGeoMatrix  *mat = (TGeoMatrix*)fArraySector->At(isec);
    if (mat){
      Double_t ppos[3]={x[0],x[1],x[2]};
      Double_t pposC[3]={x[0],x[1],x[2]};
      mat->LocalToMaster(ppos,pposC);
      dx[0]+=pposC[0]-ppos[0];
      dx[1]+=pposC[1]-ppos[1];
      dx[2]+=pposC[2]-ppos[2];
    }
  }
}

void AliTPCCalibGlobalMisalignment::Print(Option_t* /*option*/ ) const {
  //
  // Print function to check the settings 
  //
  printf("%s",GetTitle());  
  printf(" - Trivial Misalignments for calibration purposes: \n");
  printf(" - X-Shift: %1.3f cm, Y-Shift: %1.3f cm, Z-Shift: %1.3f cm \n",fXShift,fYShift,fZShift);
  printf(" - Phi-Rotations: A side: %1.5f rad, C side: %1.5f rad\n",fRotPhiA,fRotPhiC);
  printf(" - dRPhi offsets: A side: %1.5f cm, C side: %1.5f cm\n",fdRPhiOffsetA,fdRPhiOffsetC);
 
 
}

void AliTPCCalibGlobalMisalignment::AddAlign(const  AliTPCCalibGlobalMisalignment & add){
  //
  // Add the alignmnet to current object
  //
  fXShift+=add.fXShift;               // Shift in global X [cm]
  fYShift+=add.fYShift;               // Shift in global Y [cm]
  fZShift+=add.fZShift;               // Shift in global Z [cm]

  fRotPhiA+=add.fRotPhiA;      // simple rotation of A side read-out plane around the Z axis [rad]
  fRotPhiC+=add.fRotPhiC;      // simple rotation of C side read-out plane around the Z axis [rad]
  fdRPhiOffsetA+=add.fdRPhiOffsetA;  // add a constant offset of dRPhi (or local Y) in [cm]: purely for calibration purposes!
  fdRPhiOffsetC+=add.fdRPhiOffsetC;  // add a constant offset of dRPhi (or local Y) in [cm]: purely for calibration purposes!
  //
  // Quadrant alignment
  //
  if (add.fQuadrantQ0) {
    if (fQuadrantQ0) fQuadrantQ0->Add(*(add.fQuadrantQ0));
    if (!fQuadrantQ0) fQuadrantQ0 = (TVectorD*)(add.fQuadrantQ0->Clone());
  }
  if (add.fQuadrantRQ0) {
    if (fQuadrantRQ0) fQuadrantRQ0->Add(*(add.fQuadrantRQ0));
    if (!fQuadrantRQ0) fQuadrantRQ0 = (TVectorD*)(add.fQuadrantRQ0->Clone());
  }
  //
  if (add.fQuadrantQ1) {
    if (fQuadrantQ1) fQuadrantQ1->Add(*(add.fQuadrantQ1));
    if (!fQuadrantQ1) fQuadrantQ1 = (TVectorD*)(add.fQuadrantQ1->Clone());
  }
  if (add.fQuadrantRQ1) {
    if (fQuadrantRQ1) fQuadrantRQ1->Add(*(add.fQuadrantRQ1));
    if (!fQuadrantRQ1) fQuadrantRQ1 = (TVectorD*)(add.fQuadrantRQ1->Clone());
  }
  //
  if (add.fQuadrantQ2) {
    if (fQuadrantQ2) fQuadrantQ2->Add(*(add.fQuadrantQ2));
    if (!fQuadrantQ2) fQuadrantQ2 = (TVectorD*)(add.fQuadrantQ2->Clone());
  }
  if (add.fQuadrantRQ2) {
    if (fQuadrantRQ2) fQuadrantRQ2->Add(*(add.fQuadrantRQ2));
    if (!fQuadrantRQ2) fQuadrantRQ2 = (TVectorD*)(add.fQuadrantRQ2->Clone());
  }
  //
  // Global alignment - use native ROOT representation
  //
  if (add.fMatrixGlobal){
    if (!fMatrixGlobal)  fMatrixGlobal = new TGeoHMatrix(*(add.fMatrixGlobal)); 
    if (fMatrixGlobal)   ((TGeoHMatrix*)fMatrixGlobal)->Multiply(add.fMatrixGlobal); 
  }
  if (add.fArraySector){
    if (!fArraySector) {SetAlignSectors(add.fArraySector);
    }else{
      for (Int_t isec=0; isec<72; isec++){
	TGeoHMatrix *mat0= (TGeoHMatrix*)fArraySector->At(isec);
	TGeoHMatrix *mat1= (TGeoHMatrix*)add.fArraySector->At(isec);
	if (mat0&&mat1) mat0->Multiply(mat1);
      }
    }
  }
}


AliTPCCalibGlobalMisalignment *  AliTPCCalibGlobalMisalignment::CreateOCDBAlign(){
  //
  // Create  AliTPCCalibGlobalMisalignment from OCDB Alignment entry
  // OCDB has to be initialized before in user code
  // All storages (defualt and specific)  and run number 
  //
  AliCDBEntry * entry = AliCDBManager::Instance()->Get("TPC/Align/Data");
  if (!entry){
    printf("Missing alignmnet entry. OCDB not initialized?\n");
    return 0;
  }
  TClonesArray * array = (TClonesArray*)entry->GetObject();
  Int_t entries = array->GetEntries();
  TGeoHMatrix matrixGlobal;
  TObjArray *alignArrayOCDB= new TObjArray(73);  // sector misalignment + global misalignment
  //                                            // global is number 72
  //
  { for (Int_t i=0;i<entries; i++){
      //
      //
      TGeoHMatrix matrix;
      AliAlignObjParams *alignP = (AliAlignObjParams*)array->UncheckedAt(i);
      alignP->GetMatrix(matrix);
      Int_t imod;
      AliGeomManager::ELayerID ilayer;
      alignP->GetVolUID(ilayer, imod);
      if (ilayer==AliGeomManager::kInvalidLayer) {
	alignArrayOCDB->AddAt(matrix.Clone(),72);
	alignP->GetMatrix(matrixGlobal);
      }else{
	Int_t sector=imod;
	if (ilayer==AliGeomManager::kTPC2) sector+=36;
	alignArrayOCDB->AddAt(matrix.Clone(),sector);
      }
    }
  }
  AliTPCCalibGlobalMisalignment *align = new  AliTPCCalibGlobalMisalignment;
  align->SetAlignGlobal(&matrixGlobal);
  align->SetAlignSectors(alignArrayOCDB);
  return align;
}


AliTPCCalibGlobalMisalignment *  AliTPCCalibGlobalMisalignment::CreateMeanAlign(const AliTPCCalibGlobalMisalignment *alignIn){
  //
  // Create new object, disantangle common mean alignmnet and sector alignment
  //
  // 1. Try to get mean alignment
  // 2. Remove mean alignment from sector alignment
  // 3. Create new object
  //
  TObjArray * array = alignIn->GetAlignSectors();
  TObjArray * arrayNew = new TObjArray(72);
  //
  //Get mean transformation
  TGeoHMatrix matrix;  
  {for (Int_t isec=0; isec<72; isec++){
      const TGeoMatrix* cmatrix=(TGeoMatrix*)array->At(isec);
      if (!cmatrix) continue;
      matrix.Multiply(cmatrix);
    }}
  TGeoHMatrix matrixMean(matrix);
  matrixMean.SetDx(matrix.GetTranslation()[0]/72.);
  matrixMean.SetDy(matrix.GetTranslation()[1]/72.);
  matrixMean.SetDz(matrix.GetTranslation()[2]/72.);
  Double_t rotation[12];
  {for (Int_t i=0; i<12; i++) {
      rotation[i]=1.0;
      if (TMath::Abs(matrix.GetRotationMatrix()[i]-1.)>0.1){
      rotation[i]=matrix.GetRotationMatrix()[i]/72.;
      }
    }}
  matrixMean.SetRotation(rotation);
  TGeoHMatrix matrixInv = matrixMean.Inverse();
  //
  {for (Int_t isec=0; isec<72; isec++){
      TGeoHMatrix* amatrix=(TGeoHMatrix*)(array->At(isec)->Clone());
      if (!amatrix) continue;
      amatrix->Multiply(&matrixInv);
      arrayNew->AddAt(amatrix,isec);
    }}
  if (alignIn->GetAlignGlobal()) matrixMean.Multiply((alignIn->GetAlignGlobal()));
  AliTPCCalibGlobalMisalignment *alignOut = new  AliTPCCalibGlobalMisalignment;
  alignOut->SetAlignGlobal(&matrixMean);  
  alignOut->SetAlignSectors(arrayNew);  
  /*
    Checks transformation:   
    AliTPCCalibGlobalMisalignment *  alignIn =  AliTPCCalibGlobalMisalignment::CreateOCDBAlign()
    AliTPCCalibGlobalMisalignment * alignOut =  AliTPCCalibGlobalMisalignment::CreateMeanAlign(alignIn)
    alignOutM= (AliTPCCalibGlobalMisalignment*)alignOut->Clone();
    alignOutS= (AliTPCCalibGlobalMisalignment*)alignOut->Clone();
    alignOutS->SetAlignGlobal(0);  
    alignOutM->SetAlignSectors(0);  
    //
    AliTPCCorrection::AddVisualCorrection(alignOut,0);
    AliTPCCorrection::AddVisualCorrection(alignOutM,1);
    AliTPCCorrection::AddVisualCorrection(alignOutS,2);
    AliTPCCorrection::AddVisualCorrection(alignIn,3);
    
    TF1 f0("f0","AliTPCCorrection::GetCorrSector(x,85,0.9,1,0)",0,18);
    TF1 f1("f1","AliTPCCorrection::GetCorrSector(x,85,0.9,1,1)",0,18);
    TF1 f2("f2","AliTPCCorrection::GetCorrSector(x,85,0.9,1,2)",0,18);
    TF1 f3("f3","AliTPCCorrection::GetCorrSector(x,85,0.9,1,3)",0,18);
    f0->SetLineColor(1);
    f1->SetLineColor(2);
    f2->SetLineColor(3);
    f3->SetLineColor(4);
    f0->Draw();
    f1->Draw("same");
    f2->Draw("same");
    f3->Draw("same");

    TF2 f2D("f2D","AliTPCCorrection::GetCorrSector(x,y,0.9,1,0)-AliTPCCorrection::GetCorrSector(x,y,0.9,1,3)",0,18,85,245);
  */
  return alignOut;
}


void AliTPCCalibGlobalMisalignment::DumpAlignment( AliTPCCalibGlobalMisalignment* align, TTreeSRedirector *pcstream, const char *name){
  //
  // Dump alignment per sector into tree
  //
  TObjArray * array = align->GetAlignSectors();
  if (!array) return;
  //
  //Get mean transformation
  TGeoHMatrix matrix;  
  {for (Int_t isec=0; isec<72; isec++){
      TGeoHMatrix* cmatrix=(TGeoHMatrix*)array->At(isec);
      TGeoHMatrix* cmatrixDown=(TGeoHMatrix*)array->At(isec%36);
      TGeoHMatrix* cmatrixUp=(TGeoHMatrix*)array->At(isec%36+36);
      TGeoHMatrix diff(*cmatrixDown);
      diff.Multiply(&(cmatrixUp->Inverse()));
      (*pcstream)<<name<<
	"isec="<<isec<<
	"m0.="<<cmatrix<<
	"diff.="<<&diff<<
	"\n";
    }
  }
  
}
