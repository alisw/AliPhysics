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


#include <TClass.h>
#include <TFile.h>
#include <TSystem.h>

#include "AliMagF.h"
#include "AliMagWrapCheb.h"
#include "AliLog.h"

ClassImp(AliMagF)

const Double_t AliMagF::fgkSol2DipZ    =  -700.;  

//_______________________________________________________________________
AliMagF::AliMagF():
  TVirtualMagField(),
  fMeasuredMap(0),
  fMapType(k5kG),
  fSolenoid(0),
  fBeamType(kNoBeamField),
  fBeamEnergy(0),
  //
  fInteg(0),
  fPrecInteg(0),
  fFactorSol(1.),
  fFactorDip(1.),
  fMax(15),
  fDipoleOFF(kFALSE),
  //
  fQuadGradient(0),
  fDipoleField(0),
  fCCorrField(0), 
  fACorr1Field(0),
  fACorr2Field(0),
  fParNames("","")
{
  // Default constructor
  //
}

//_______________________________________________________________________
AliMagF::AliMagF(const char *name, const char* title, Int_t integ, 
		 Double_t factorSol, Double_t factorDip, 
		 Double_t fmax, BMap_t maptype, const char* path,
		 BeamType_t bt, Double_t be):
  TVirtualMagField(name),
  fMeasuredMap(0),
  fMapType(maptype),
  fSolenoid(0),
  fBeamType(bt),
  fBeamEnergy(be),
  //
  fInteg(integ),
  fPrecInteg(1),
  fFactorSol(1.),
  fFactorDip(1.),
  fMax(fmax),
  fDipoleOFF(factorDip==0.),
  //
  fQuadGradient(0),
  fDipoleField(0),
  fCCorrField(0), 
  fACorr1Field(0),
  fACorr2Field(0),
  fParNames("","")
{
  // Initialize the field with Geant integration option "integ" and max field "fmax,
  // Impose scaling of parameterized L3 field by factorSol and of dipole by factorDip.
  // The "be" is the energy of the beam in GeV/nucleon
  //
  SetTitle(title);
  if(integ<0 || integ > 2) {
    AliWarning(Form("Invalid magnetic field flag: %5d; Helix tracking chosen instead",integ));
    fInteg = 2;
  }
  if (fInteg == 0) fPrecInteg = 0;
  //
  const char* parname = 0;
  //  
  if (fMapType == k2kG) {
    fSolenoid = 2.;
    parname = fDipoleOFF ? "Sol12_Dip0_Hole":"Sol12_Dip6_Hole";
  } else if (fMapType == k5kG) {
    fSolenoid = 5.;
    parname = fDipoleOFF ? "Sol30_Dip0_Hole":"Sol30_Dip6_Hole";
  } else if (fMapType == k5kGUniform) {
    fSolenoid = 5.;
    parname = "Sol30_Dip6_Uniform";
  } else {
    AliFatal(Form("Unknown field identifier %d is requested\n",fMapType)); 
  }
  //
  SetDataFileName(path);
  SetParamName(parname);
  //
  SetFactorSol(factorSol);
  SetFactorDip(factorDip);
  LoadParameterization();
  InitMachineField(fBeamType,fBeamEnergy);
}

//_______________________________________________________________________
AliMagF::AliMagF(const AliMagF &src):
  TVirtualMagField(src),
  fMeasuredMap(0),
  fMapType(src.fMapType),
  fSolenoid(src.fSolenoid),
  fBeamType(src.fBeamType),
  fBeamEnergy(src.fBeamEnergy),
  fInteg(src.fInteg),
  fPrecInteg(src.fPrecInteg),
  fFactorSol(src.fFactorSol),
  fFactorDip(src.fFactorDip),
  fMax(src.fMax),
  fDipoleOFF(src.fDipoleOFF),
  fQuadGradient(src.fQuadGradient),
  fDipoleField(src.fDipoleField),
  fCCorrField(src.fCCorrField), 
  fACorr1Field(src.fACorr1Field),
  fACorr2Field(src.fACorr2Field),
  fParNames(src.fParNames)
{
  if (src.fMeasuredMap) fMeasuredMap = new AliMagWrapCheb(*src.fMeasuredMap);
}

//_______________________________________________________________________
AliMagF::~AliMagF()
{
  delete fMeasuredMap;
}

//_______________________________________________________________________
Bool_t AliMagF::LoadParameterization()
{
  if (fMeasuredMap) {
    AliError(Form("Field data %s are already loaded from %s\n",GetParamName(),GetDataFileName()));
    return kTRUE;
  }
  //
  char* fname = gSystem->ExpandPathName(GetDataFileName());
  TFile* file = TFile::Open(fname);
  if (!file) {
    AliError(Form("Failed to open magnetic field data file %s\n",fname)); 
    return kFALSE;
  }
  //
  fMeasuredMap = dynamic_cast<AliMagWrapCheb*>(file->Get(GetParamName()));
  if (!fMeasuredMap) {
    AliError(Form("Did not find field %s in %s\n",GetParamName(),fname)); 
    return kFALSE;
  }
  file->Close();
  delete file;
  return kTRUE;
}


//_______________________________________________________________________
void AliMagF::Field(const Double_t *xyz, Double_t *b)
{
  // Method to calculate the field at point  xyz
  //
  //  b[0]=b[1]=b[2]=0.0;
  if (fMeasuredMap && xyz[2]>fMeasuredMap->GetMinZ() && xyz[2]<fMeasuredMap->GetMaxZ()) {
    fMeasuredMap->Field(xyz,b);
    if (xyz[2]>fgkSol2DipZ || fDipoleOFF) for (int i=3;i--;) b[i] *= fFactorSol;
    else                                  for (int i=3;i--;) b[i] *= fFactorDip;    
  }
  else MachineField(xyz, b);
  //
}

//_______________________________________________________________________
Double_t AliMagF::GetBz(const Double_t *xyz) const
{
  // Method to calculate the field at point  xyz
  //
  if (fMeasuredMap && xyz[2]>fMeasuredMap->GetMinZ() && xyz[2]<fMeasuredMap->GetMaxZ()) {
    double bz = fMeasuredMap->GetBz(xyz);
    return (xyz[2]>fgkSol2DipZ || fDipoleOFF) ? bz*fFactorSol : bz*fFactorDip;    
  }
  else return 0.;
}

//_______________________________________________________________________
AliMagF& AliMagF::operator=(const AliMagF& src)
{
  if (this != &src && src.fMeasuredMap) { 
    if (fMeasuredMap) delete fMeasuredMap;
    fMeasuredMap = new AliMagWrapCheb(*src.fMeasuredMap);
    SetName(src.GetName());
    fSolenoid    = src.fSolenoid;
    fBeamType    = src.fBeamType;
    fBeamEnergy  = src.fBeamEnergy;
    fInteg       = src.fInteg;
    fPrecInteg   = src.fPrecInteg;
    fFactorSol   = src.fFactorSol;
    fFactorDip   = src.fFactorDip;
    fMax         = src.fMax;
    fDipoleOFF   = src.fDipoleOFF;
    fParNames    = src.fParNames;
  }
  return *this;
}

//_______________________________________________________________________
void AliMagF::InitMachineField(BeamType_t btype, Double_t benergy)
{
  if (btype==kNoBeamField || benergy<1.) {
    fQuadGradient = fDipoleField = fCCorrField = fACorr1Field = fACorr2Field = 0.;
    return;
  }
  //
  double rigScale = benergy/7000.;   // scale according to ratio of E/Enominal
  // for ions assume PbPb (with energy provided per nucleon) and account for A/Z
  if (btype == kBeamTypeAA) rigScale *= 208./82.;
  //
  fQuadGradient = 22.0002*rigScale;
  fDipoleField  = 37.8781*rigScale;
  //
  // SIDE C
  fCCorrField   = -9.6980;
  // SIDE A
  fACorr1Field  = -13.2247;
  fACorr2Field  =  11.7905;
  //
}

//_______________________________________________________________________
void AliMagF::MachineField(const Double_t *x, Double_t *b) const
{
  // ---- This is the ZDC part
  // Compansators for Alice Muon Arm Dipole
  const Double_t kBComp1CZ = 1075., kBComp1hDZ = 260./2., kBComp1SqR = 4.0*4.0; 
  const Double_t kBComp2CZ = 2049., kBComp2hDZ = 153./2., kBComp2SqR = 4.5*4.5; 
  //  
  const Double_t kTripQ1CZ = 2615., kTripQ1hDZ = 637./2., kTripQ1SqR = 3.5*3.5;
  const Double_t kTripQ2CZ = 3408., kTripQ2hDZ = 550./2., kTripQ2SqR = 3.5*3.5;
  const Double_t kTripQ3CZ = 4130., kTripQ3hDZ = 550./2., kTripQ3SqR = 3.5*3.5;
  const Double_t kTripQ4CZ = 5015., kTripQ4hDZ = 637./2., kTripQ4SqR = 3.5*3.5;
  //
  const Double_t kDip1CZ = 6310.8,  kDip1hDZ = 945./2., kDip1SqRC = 4.5*4.5, kDip1SqRA = 3.375*3.375;
  const Double_t kDip2CZ = 12640.3, kDip2hDZ = 945./2., kDip2SqRC = 4.5*4.5, kDip2SqRA = 3.75*3.75;
  const Double_t kDip2DXC = 9.7, kDip2DXA = 9.4;
  //
  double rad2 = x[0] * x[0] + x[1] * x[1];
  //
  b[0] = b[1] = b[2] = 0;
  //
  // SIDE C **************************************************
  if(x[2]<0.){  
    if(TMath::Abs(x[2]+kBComp2CZ)<kBComp2hDZ && rad2 < kBComp2SqR){
      b[0] = fCCorrField*fFactorDip;
    } 
    else if(TMath::Abs(x[2]+kTripQ1CZ)<kTripQ1hDZ && rad2 < kTripQ1SqR){
      b[0] = fQuadGradient*x[1];
      b[1] = fQuadGradient*x[0];
    }
    else if(TMath::Abs(x[2]+kTripQ2CZ)<kTripQ2hDZ && rad2 < kTripQ2SqR){
      b[0] = -fQuadGradient*x[1];
      b[1] = -fQuadGradient*x[0];
    }
    else if(TMath::Abs(x[2]+kTripQ3CZ)<kTripQ3hDZ && rad2 < kTripQ3SqR){
      b[0] = -fQuadGradient*x[1];
      b[1] = -fQuadGradient*x[0];
    }
    else if(TMath::Abs(x[2]+kTripQ4CZ)<kTripQ4hDZ && rad2 < kTripQ4SqR){
      b[0] = fQuadGradient*x[1];
      b[1] = fQuadGradient*x[0];
    }
    else if(TMath::Abs(x[2]+kDip1CZ)<kDip1hDZ && rad2 < kDip1SqRC){
      b[1] = fDipoleField;
    }
    else if(TMath::Abs(x[2]+kDip2CZ)<kDip2hDZ && rad2 < kDip2SqRC) {
      double dxabs = TMath::Abs(x[0])-kDip2DXC;
      if ( (dxabs*dxabs + x[1]*x[1])<kDip2SqRC) {
	b[1] = -fDipoleField;
      }
    }
  }
  //
  // SIDE A **************************************************
  else{        
    if(TMath::Abs(x[2]-kBComp1CZ)<kBComp1hDZ && rad2 < kBComp1SqR) {
      // Compensator magnet at z = 1075 m 
      b[0] = fACorr1Field*fFactorDip;
    }
    //
    if(TMath::Abs(x[2]-kBComp2CZ)<kBComp2hDZ && rad2 < kBComp2SqR){
      b[0] = fACorr2Field*fFactorDip;
    }
    else if(TMath::Abs(x[2]-kTripQ1CZ)<kTripQ1hDZ && rad2 < kTripQ1SqR){
      b[0] = -fQuadGradient*x[1];
      b[1] = -fQuadGradient*x[0];
    }
    else if(TMath::Abs(x[2]-kTripQ2CZ)<kTripQ2hDZ && rad2 < kTripQ2SqR){
      b[0] =  fQuadGradient*x[1];
      b[1] =  fQuadGradient*x[0];
    }
    else if(TMath::Abs(x[2]-kTripQ3CZ)<kTripQ3hDZ && rad2 < kTripQ3SqR){
      b[0] =  fQuadGradient*x[1];
      b[1] =  fQuadGradient*x[0];
    }
    else if(TMath::Abs(x[2]-kTripQ4CZ)<kTripQ4hDZ && rad2 < kTripQ4SqR){
      b[0] = -fQuadGradient*x[1];
      b[1] = -fQuadGradient*x[0];
    }
    else if(TMath::Abs(x[2]-kDip1CZ)<kDip1hDZ && rad2 < kDip1SqRA){
      b[1] = -fDipoleField;
    }
    else if(TMath::Abs(x[2]-kDip2CZ)<kDip2hDZ && rad2 < kDip2SqRA) {
      double dxabs = TMath::Abs(x[0])-kDip2DXA;
      if ( (dxabs*dxabs + x[1]*x[1])<kDip2SqRA) {
	b[1] = fDipoleField;
      }
    }
  }
  //
}

//_______________________________________________________________________
void AliMagF::GetTPCInt(const Double_t *xyz, Double_t *b) const
{
  // Method to calculate the integral of magnetic integral from xyz to nearest cathode plane
  b[0]=b[1]=b[2]=0.0;
  if (fMeasuredMap) {
    fMeasuredMap->GetTPCInt(xyz,b);
    for (int i=3;i--;) b[i] *= fFactorSol;
  }
}

//_______________________________________________________________________
void AliMagF::GetTPCIntCyl(const Double_t *rphiz, Double_t *b) const
{
  // Method to calculate the integral of magnetic integral from point to nearest cathode plane
  // in cylindrical coordiates ( -pi<phi<pi convention )
  b[0]=b[1]=b[2]=0.0;
  if (fMeasuredMap) {
    fMeasuredMap->GetTPCIntCyl(rphiz,b);
    for (int i=3;i--;) b[i] *= fFactorSol;
  }
}
