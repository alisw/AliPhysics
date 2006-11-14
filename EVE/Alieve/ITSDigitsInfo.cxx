// $Header$

//__________________________________________________________________________
// ITSDigitsInfo
//
//

#include <Reve/TTreeTools.h>

#include "ITSDigitsInfo.h"
//#include <AliITSresponseSDD.h>
#include <AliITSCalibrationSDD.h>
#include <AliITSdigit.h>
#include <AliITSdigitSPD.h>

using namespace Reve;
using namespace Alieve;
using namespace std;

ClassImp(ITSDigitsInfo)

/**************************************************************************/

ITSDigitsInfo::ITSDigitsInfo() :
  TObject(),
  ReferenceCount(),
  fSPDmap(), fSDDmap(), fSSDmap(),
  fTree (0),
  fGeom (0),
  fSegSPD(0), fSegSDD(0), fSegSSD(0)
{}

/**************************************************************************/

ITSDigitsInfo:: ~ITSDigitsInfo() 
{
  map<Int_t, TClonesArray*>::iterator j;
  for(j = fSPDmap.begin(); j != fSPDmap.end(); ++j)
    delete j->second;
  for(j = fSDDmap.begin(); j != fSDDmap.end(); ++j)
    delete j->second;
  for(j = fSSDmap.begin(); j != fSSDmap.end(); ++j)
    delete j->second;

  delete fSegSPD; delete fSegSDD; delete fSegSSD; 
  delete fGeom;
  delete fTree;
}

/**************************************************************************/

void ITSDigitsInfo::SetTree(TTree* tree)
{
  static const Exc_t eH("ITSDigitsInfo::SetTree ");

  if(fGeom == 0) {
    fGeom = new AliITSgeom();
    fGeom->ReadNewFile("$REVESYS/alice-data/ITSgeometry.det");
    if(fGeom == 0)
      throw(eH + "can not load ITS geometry \n");
  }

  fTree = tree;

  SetITSSegmentation();
  
  // create tables for scaling
  // lowest scale factor refers to unscaled ITS module
  fSPDScaleX[0] = 1;
  fSPDScaleZ[0] = 1;
  fSDDScaleX[0] = 1;
  fSDDScaleZ[0] = 1;
  fSSDScale [0] = 1;
  // spd lows resolution is in the level of 8x2 readout chips
  Int_t nx = 8; // fSegSPD->Npx()/8; // 32
  Int_t nz = 6; // fSegSPD->Npz()/2; // 128

  fSPDScaleX[1] = Int_t(nx); 
  fSPDScaleZ[1] = Int_t(nz); 
  fSPDScaleX[2] = Int_t(nx*2); 
  fSPDScaleZ[2] = Int_t(nz*2); 
  fSPDScaleX[3] = Int_t(nx*3); 
  fSPDScaleZ[3] = Int_t(nz*3); 
  fSPDScaleX[4] = Int_t(nx*4); 
  fSPDScaleZ[4] = Int_t(nz*4); 


  fSDDScaleX[1] = 2;
  fSDDScaleZ[1] = 2;
  fSDDScaleX[2] = 8;
  fSDDScaleZ[2] = 8;
  fSDDScaleX[3] = 16;
  fSDDScaleZ[3] = 16;
  fSDDScaleX[4] = 25;
  fSDDScaleZ[4] = 25;

  fSSDScale[1] = 3;
  fSSDScale[2] = 9;
  fSSDScale[3] = 20;
  fSSDScale[4] = 30;

  
  // lowest scale factor refers unscaled ITS module
  fSPDScaleX[0] = 1;
  fSPDScaleZ[0] = 1;
  fSDDScaleX[0] = 1;
  fSDDScaleZ[0] = 1;
  fSSDScale [0] = 1;
}

/**************************************************************************/

void ITSDigitsInfo::SetITSSegmentation()
{
  // SPD
  fSegSPD = new AliITSsegmentationSPD(fGeom);
  //SPD geometry  
  Int_t m;
  Float_t fNzSPD=160;
  Float_t fZ1pitchSPD=0.0425; Float_t fZ2pitchSPD=0.0625;
  Float_t fHlSPD=3.48;

  fSPDZCoord[0]=fZ1pitchSPD -fHlSPD;
  for (m=1; m<fNzSPD; m++) {
    Double_t dz=fZ1pitchSPD;
    if (m==31 || m==32 || m==63  || m==64  || m==95 || m==96 || 
        m==127 || m==128) dz=fZ2pitchSPD; 
    fSPDZCoord[m]=fSPDZCoord[m-1]+dz;
  }
  
  for (m=0; m<fNzSPD; m++) {
    Double_t dz=1.*fZ1pitchSPD;
    if (m==31 || m==32 || m==63  || m==64  || m==95 || m==96 || 
	m==127 || m==128) dz=1.*fZ2pitchSPD; 
    fSPDZCoord[m]-=dz;
  }
  
  // end of SPD geometry
  
  // SDD
  // response replaced by Calibration (March 2006).
  AliITSresponseSDD*   resp1 = new AliITSresponseSDD();
  AliITSCalibrationSDD* cal1 = new AliITSCalibrationSDD;
  cal1->SetResponse(resp1);
  fSegSDD = new AliITSsegmentationSDD(fGeom, cal1);

  // SSD
  fSegSSD = new AliITSsegmentationSSD(fGeom);
}

void ITSDigitsInfo::GetSPDLocalZ(Int_t j, Float_t& z)
{
  z = fSPDZCoord[j];
}

/**************************************************************************/

TClonesArray* ITSDigitsInfo::GetDigits(Int_t mod, Int_t subdet)
{
  switch(subdet) {
  case 0: {
    TClonesArray* digitsSPD = 0;
    map<Int_t, TClonesArray*>::iterator i = fSPDmap.find(mod);
    if(i == fSPDmap.end()) {
      fTree->SetBranchAddress("ITSDigitsSPD",&digitsSPD);
      fTree->GetEntry(mod);
      fSPDmap[mod] = digitsSPD;
      return digitsSPD;
    } else {
      return i->second;
    }
    break;
  }
  case 1: {
    TClonesArray* digitsSDD = 0;
    map<Int_t, TClonesArray*>::iterator i = fSDDmap.find(mod);
    if(i == fSDDmap.end()) {
      fTree->SetBranchAddress("ITSDigitsSDD",&digitsSDD);
      fTree->GetEntry(mod);
      fSDDmap[mod] = digitsSDD;
      return digitsSDD;
    } else {
      return i->second;
    }
    break;
  }
  case 2: {
    TClonesArray* digitsSSD = 0;
    map<Int_t, TClonesArray*>::iterator i = fSSDmap.find(mod);
    if(i == fSSDmap.end()) {
      fTree->SetBranchAddress("ITSDigitsSSD",&digitsSSD);
      fTree->GetEntry(mod);
      fSSDmap[mod] = digitsSSD;
      return digitsSSD;
    } else {
      return i->second;
    }
    break;
  }
  default:
    return 0;
  } //end switch
  return 0;
}


/**************************************************************************/

void ITSDigitsInfo::Print(Option_t* ) const
{
  printf("*********************************************************\n");
  printf("SPD module dimension (%f,%f) \n",fSegSPD->Dx()*0.0001, fSegSPD->Dz()*0.0001);
  printf("SPD first,last module:: %d,%d \n", fGeom->GetStartSPD(),fGeom->GetLastSPD() );
  printf("SPD num cells per module (x::%d,z::%d)\n",fSegSPD->Npx(), fSegSPD->Npz());
  Int_t iz=0,ix = 0;
  printf("SPD dimesion of (%d,%d) in pixel(%f,%f) \n", ix,iz, fSegSPD->Dpx(ix), fSegSPD->Dpz(iz));
  iz = 32;
  printf("SPD dimesion of pixel (%d,%d) are (%f,%f) \n", ix,iz, fSegSPD->Dpx(ix)*0.001, fSegSPD->Dpz(iz)*0.001);
 
  printf("*********************************************************\n");
  printf("SDD module dimension (%f,%f) \n",fSegSDD->Dx()*0.0001, fSegSDD->Dz()*0.0001);
  printf("SDD first,last module:: %d,%d \n", fGeom->GetStartSDD(),fGeom->GetLastSDD() );
  printf("SDD num cells per module (x::%d,z::%d)\n",fSegSDD->Npx(), fSegSDD->Npz());
  printf("SDD dimesion of pixel are (%f,%f) \n", fSegSDD->Dpx(1)*0.001,fSegSDD->Dpz(1)*0.001);
  printf("*********************************************************\n");
  printf("SSD module dimension (%f,%f) \n",fSegSSD->Dx()*0.0001, fSegSSD->Dz()*0.0001);
  printf("SSD first,last module:: %d,%d \n", fGeom->GetStartSSD(),fGeom->GetLastSSD() );
  printf("SSD strips in module %d \n",fSegSSD->Npx());
  printf("SSD strip sizes are (%f,%f) \n", fSegSSD->Dpx(1),fSegSSD->Dpz(1));
  fSegSSD->SetLayer(5); Float_t ap,an;  fSegSSD->Angles(ap,an);
  printf("SSD layer 5 stereoP %f stereoN %f angle \n",ap,an); 
  fSegSSD->SetLayer(6);  fSegSSD->Angles(ap,an);
  printf("SSD layer 6 stereoP %f stereoN %f angle \n",ap,an); 
}


/*
  printf("num cells %d,%d scaled %d,%d \n",fSegSPD->Npz(),fSegSPD->Npx(),Nz,Nx);
  printf("%d digits in ITSModule %d\n",ne, module);
  Float_t zn = i*(3.48*2)/Nz - 3.48 ;
  Float_t xo =  -fSegSPD->Dx()*0.00005 + fSegSPD->Dpx(0)*od->GetCoord2()*0.0001;
  Float_t xn =  -fSegSPD->Dx()*0.00005 + j*0.0001*fSegSPD->Dx()/Nx;
  Float_t dpx = 0.0001*fSegSPD->Dx()/Nx;
  Float_t dpz = 3.48*2/Nz;
  printf("Z::original (%3f) scaled (%3f, %3f) \n", zo, zn-dpz/2, zn+dpz/2);
  printf("X::original (%3f) scaled (%3f, %3f) \n", xo, xn-dpx/2, xn+dpx/2);
  printf("%d,%d maped to %d,%d \n", od->GetCoord1(), od->GetCoord2(), i,j );
*/        
