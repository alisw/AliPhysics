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

const Double_t AliMagF::fgkSol2DipZ   =  -700.;
const Double_t AliMagF::fgkBMachineZ1 =   919.;
const Double_t AliMagF::fgkBMachineZ2 = -1972.;

//_______________________________________________________________________
AliMagF::AliMagF():
  TVirtualMagField(),
  fMeasuredMap(0),
  fMapType(k5kG),
  fSolenoid(0),
  fBeamType(kNoBeamField),
  fBeamEnergy(0),
  fCompensator(kFALSE),
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
		 BeamType_t bt, Double_t be, Bool_t compensator):
  TVirtualMagField(name),
  fMeasuredMap(0),
  fMapType(maptype),
  fSolenoid(0),
  fBeamType(bt),
  fBeamEnergy(be),
  fCompensator(compensator),
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
  fCompensator(src.fCompensator),
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
  b[0]=b[1]=b[2]=0.0;
  if (xyz[2] > fgkBMachineZ1 || xyz[2] < fgkBMachineZ2) MachineField(xyz, b);
  else if (fMeasuredMap) {
    fMeasuredMap->Field(xyz,b);
    if (xyz[2]>fgkSol2DipZ || fDipoleOFF) for (int i=3;i--;) b[i] *= fFactorSol;
    else                                  for (int i=3;i--;) b[i] *= fFactorDip;
  }
  //
}

//_______________________________________________________________________
Double_t AliMagF::GetBz(const Double_t *xyz) const
{
  // Method to calculate the field at point  xyz
  //
  if (xyz[2] <= fgkBMachineZ1 && xyz[2] >= fgkBMachineZ2) return fMeasuredMap->GetBz(xyz);
  else {
    double b[3] = {0,0,0};
    MachineField(xyz, b);
    return b[2];
  }
  //
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
    fCompensator = src.fCompensator;
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
  const double kToler = 0.1;
  if (btype==kNoBeamField) {
    fQuadGradient = fDipoleField = fCCorrField = fACorr1Field = fACorr2Field = 0.;
  }
  //
  else if (btype==kBeamTypepp && TMath::Abs(1.-benergy/5000.)<kToler ){
    // p-p @ 5+5 TeV
    fQuadGradient = 15.7145;
    fDipoleField  = 27.0558;
    // SIDE C
    fCCorrField   = 9.7017;
    // SIDE A
    fACorr1Field  = -13.2143;
    fACorr2Field  = -11.9909;
  } else if (btype == kBeamTypepp && TMath::Abs(1.-benergy/450.)<kToler) {
    // p-p 0.45+0.45 TeV
    Double_t const kEnergyRatio = benergy / 7000.;
    fQuadGradient = 22.0002 * kEnergyRatio;
    fDipoleField  = 37.8781 * kEnergyRatio;
    // SIDE C
    fCCorrField   =  9.6908;
    // SIDE A
    fACorr1Field  = -13.2014;
    fACorr2Field  = -9.6908;
  } else if ( (btype == kBeamTypepp && TMath::Abs(1.-benergy/7000.)<kToler) || 
	      (fBeamType == kBeamTypeAA) ) {
    // Pb-Pb @ 2.7+2.7 TeV or p-p @ 7+7 TeV
    fQuadGradient = 22.0002;
    fDipoleField  = 37.8781;
    // SIDE C
    fCCorrField   = 9.6908;
    // SIDE A
    fACorr1Field  = -13.2014;
    fACorr2Field  = -9.6908;
  }
  //
}

//_______________________________________________________________________
void AliMagF::MachineField(const Double_t *x, Double_t *b) const
{
  // ---- This is the ZDC part
  const Double_t kCCorrBegin = fgkBMachineZ2-0.5,kCCorrEnd = kCCorrBegin - 153., kCCorrSqRadius = 4.5*4.5;
  //
  const Double_t kCTripletBegin  = -2296.5;
  const Double_t kCQ1Begin = kCTripletBegin,        kCQ1End = kCQ1Begin-637., kCQ1SqRadius = 3.5*3.5;
  const Double_t kCQ2Begin = kCTripletBegin-908.5,  kCQ2End = kCQ2Begin-550., kCQ2SqRadius = 3.5*3.5;
  const Double_t kCQ3Begin = kCTripletBegin-1558.5, kCQ3End = kCQ3Begin-550., kCQ3SqRadius = 3.5*3.5;
  const Double_t kCQ4Begin = kCTripletBegin-2400.,  kCQ4End = kCQ4Begin-637., kCQ4SqRadius = 3.5*3.5;
  //
  const Double_t kCD1Begin = -5838.3,  kCD1End = kCD1Begin-945., kCD1SqRadius = 4.5*4.5;
  const Double_t kCD2Begin = -12167.8, kCD2End = kCD2Begin-945., kCD2SqRadius = 4.5*4.5;
  const Double_t kCD2XCentre1 = -9.7;
  const Double_t kCD2XCentre2 =  9.7;
  //
  // -> SIDE A
  // NB -> kACorr1Begin = 919. to be checked
  const Double_t kACorr1Begin = fgkBMachineZ1, kACorr1End = kACorr1Begin+260., kCCorr1SqRadius = 4.*4.;
  const Double_t kACorr2Begin = -fgkBMachineZ2 + 0.5, kACorr2End = kACorr2Begin+153., kCCorr2SqRadius = 4.5*4.5;
  const Double_t kATripletBegin  = 2296.5;
  const Double_t kAQ1Begin = kATripletBegin,	kAQ1End = kAQ1Begin+637., kAQ1SqRadius = 3.5*3.5;
  const Double_t kAQ2Begin = kATripletBegin+908.5,  kAQ2End = kAQ2Begin+550., kAQ2SqRadius = 3.5*3.5;
  const Double_t kAQ3Begin = kATripletBegin+1558.5, kAQ3End = kAQ3Begin+550., kAQ3SqRadius = 3.5*3.5;
  const Double_t kAQ4Begin = kATripletBegin+2400.,  kAQ4End = kAQ4Begin+637., kAQ4SqRadius = 3.5*3.5;
  //
  const Double_t kAD1Begin = 5838.3,  kAD1End = kAD1Begin+945., kAD1SqRadius = 3.375*3.375;
  const Double_t kAD2Begin = 12167.8, kAD2End = kAD2Begin+945., kAD2SqRadius = 3.75*3.75;
  const Double_t kAD2XCentre1 = -9.4;
  const Double_t kAD2XCentre2 =  9.4;
  //
  double rad2 = x[0] * x[0] + x[1] * x[1];
  //
  // SIDE C **************************************************
  if(x[2]<0.){  
    if(x[2] < kCCorrBegin && x[2] > kCCorrEnd && rad2 < kCCorrSqRadius){
      b[0] = fCCorrField;
      b[1] = 0.;
      b[2] = 0.;
    } 
    else if(x[2] < kCQ1Begin && x[2] > kCQ1End && rad2 < kCQ1SqRadius){
      b[0] = fQuadGradient*x[1];
      b[1] = fQuadGradient*x[0];
      b[2] = 0.;
    }
    else if(x[2] < kCQ2Begin && x[2] > kCQ2End && rad2 < kCQ2SqRadius){
      b[0] = -fQuadGradient*x[1];
      b[1] = -fQuadGradient*x[0];
      b[2] = 0.;
    }
    else if(x[2] < kCQ3Begin && x[2] > kCQ3End && rad2 < kCQ3SqRadius){
      b[0] = -fQuadGradient*x[1];
      b[1] = -fQuadGradient*x[0];
      b[2] = 0.;
    }
    else if(x[2] < kCQ4Begin && x[2] > kCQ4End && rad2 < kCQ4SqRadius){
      b[0] = fQuadGradient*x[1];
      b[1] = fQuadGradient*x[0];
      b[2] = 0.;
    }
    else if(x[2] < kCD1Begin && x[2] > kCD1End && rad2 < kCD1SqRadius){
      b[1] = fDipoleField;
      b[2] = 0.;
      b[2] = 0.;
    }
    else if(x[2] < kCD2Begin && x[2] > kCD2End){
      if(((x[0]-kCD2XCentre1)*(x[0]-kCD2XCentre1)+(x[1]*x[1]))<kCD2SqRadius
	 || ((x[0]-kCD2XCentre2)*(x[0]-kCD2XCentre2)+(x[1]*x[1]))<kCD2SqRadius){
	b[1] = -fDipoleField;
	b[2] = 0.;
	b[2] = 0.;
      }
    }
  }
  //
  // SIDE A **************************************************
  else{        
    if(fCompensator && (x[2] > kACorr1Begin && x[2] < kACorr1End) && rad2 < kCCorr1SqRadius) {
      // Compensator magnet at z = 1075 m 
      b[0] = fACorr1Field;
      b[1] = 0.;
      b[2] = 0.;
      return;
    }
    //
    if(x[2] > kACorr2Begin && x[2] < kACorr2End && rad2 < kCCorr2SqRadius){
      b[0] = fACorr2Field;
      b[1] = 0.;
      b[2] = 0.;
    }          
    else if(x[2] > kAQ1Begin && x[2] < kAQ1End && rad2 < kAQ1SqRadius){
      // First quadrupole of inner triplet de-focussing in x-direction
      b[0] = -fQuadGradient*x[1];
      b[1] = -fQuadGradient*x[0];
      b[2] = 0.;
    }
    else if(x[2] > kAQ2Begin && x[2] < kAQ2End && rad2 < kAQ2SqRadius){
      b[0] = fQuadGradient*x[1];
      b[1] = fQuadGradient*x[0];
      b[2] = 0.;
    }
    else if(x[2] > kAQ3Begin && x[2] < kAQ3End && rad2 < kAQ3SqRadius){
      b[0] = fQuadGradient*x[1];
      b[1] = fQuadGradient*x[0];
      b[2] = 0.;
    }
    else if(x[2] > kAQ4Begin && x[2] < kAQ4End && rad2 < kAQ4SqRadius){
      b[0] = -fQuadGradient*x[1];
      b[1] = -fQuadGradient*x[0];
      b[2] = 0.;
    }
    else if(x[2] > kAD1Begin && x[2] < kAD1End && rad2 < kAD1SqRadius){
      b[0] = 0.;
      b[1] = -fDipoleField;
      b[2] = 0.;
    }
    else if(x[2] > kAD2Begin && x[2] < kAD2End){
      if(((x[0]-kAD2XCentre1)*(x[0]-kAD2XCentre1)+(x[1]*x[1])) < kAD2SqRadius
	 || ((x[0]-kAD2XCentre2)*(x[0]-kAD2XCentre2)+(x[1]*x[1])) < kAD2SqRadius){
	b[1] = fDipoleField;
      }
    }
  }
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
