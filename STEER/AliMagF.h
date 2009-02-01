#ifndef ALIMAGF_H
#define ALIMAGF_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//
// Interface between the TVirtualMagField and AliMagWrapCheb: wrapper to
// the set of magnetic field data + Tosca parameterization by 
// Chebyshev polynomials
// 
// Author: ruben.shahoyan@cern.ch
//

//#include <TGeoGlobalMagField.h>
#include <TVirtualMagField.h>
class AliMagWrapCheb;

class AliMagF : public TVirtualMagField
{
 public:
  enum BMap_t      {k2kG, k5kG, k5kGUniform};
  enum BeamType_t  {kBeamTypeAA, kBeamTypepp, kNoBeamField};
  //
  AliMagF();
  AliMagF(const char *name, const char* title, Int_t integ, 
	  Double_t factorSol=1., Double_t factorDip=1., 
	  Double_t fmax=15, BMap_t maptype = k5kG,
	  const char* path="$(ALICE_ROOT)/data/maps/mfchebKGI_sym.root",
	  BeamType_t btype=kBeamTypepp, Double_t benergy=7000., Bool_t compensator=kFALSE);
  AliMagF(const AliMagF& src);             
  AliMagF& operator=(const AliMagF& src);
  virtual ~AliMagF();
  //
  virtual void Field(const Double_t *x, Double_t *b);
  void       GetTPCInt(const Double_t *xyz, Double_t *b)        const;
  void       GetTPCIntCyl(const Double_t *rphiz, Double_t *b)   const;
  Double_t   GetBz(const Double_t *xyz)                         const;
  //
  AliMagWrapCheb* GetMeasuredMap()                              const {return fMeasuredMap;}
  //
  // former AliMagF methods or their aliases
  void         SetFactorSol(Float_t fc=1.)                            {fFactorSol = fc;}
  void         SetFactorDip(Float_t fc=1.)                            {fFactorDip = fc;}
  Double_t     GetFactorSol()                                   const {return fFactorSol;}
  Double_t     GetFactorDip()                                   const {return fFactorSol;}
  Double_t     Factor()                                         const {return GetFactorSol();}
  Bool_t       IsUniform()                                      const {return fMapType == k5kGUniform;}
  //
  void         MachineField(const Double_t  *x, Double_t *b)    const;
  BMap_t       GetMapType()                                     const {return fMapType;}
  Bool_t       GetCompensator()                                 const {return fCompensator;}
  BeamType_t   GetBeamType()                                    const {return fBeamType;}
  Double_t     GetBeamEnergy()                                  const {return fBeamEnergy;}
  Double_t     Max()                                            const {return fMax;}
  Int_t        Integ()                                          const {return fInteg;}
  Int_t        PrecInteg()                                      const {return fPrecInteg;}  
  Double_t     SolenoidField()                                  const {return -fFactorSol*fSolenoid;}
  //
  Char_t*      GetDataFileName()                                const {return (Char_t*)fParNames.GetName();}
  Char_t*      GetParamName()                                   const {return (Char_t*)fParNames.GetTitle();}
  void         SetDataFileName(const Char_t* nm)                      {fParNames.SetName(nm);}
  void         SetParamName(const Char_t* nm)                         {fParNames.SetTitle(nm);}
  //
  Bool_t       LoadParameterization();
  //
 protected:
  // not supposed to be changed during the run, set only at the initialization via constructor
  void         InitMachineField(BeamType_t btype, Double_t benergy);
  void         SetBeamType(BeamType_t type)                           {fBeamType = type;}
  void         SetBeamEnergy(Float_t energy)                          {fBeamEnergy = energy;}
  void         SetCompensatorMagnet(Bool_t flag)                      {fCompensator = flag;}
  //
 protected:
  AliMagWrapCheb*  fMeasuredMap;     //! Measured part of the field map
  BMap_t           fMapType;         // field map type
  Double_t         fSolenoid;        // Solenoid field setting
  BeamType_t       fBeamType;        // Beam type: A-A (fBeamType=0) or p-p (fBeamType=1)
  Double_t         fBeamEnergy;      // Beam energy in GeV
  Bool_t           fCompensator;     // Flag for compensator magnetic field (kTrue -> ON)
  // 
  Int_t            fInteg;           // Default integration method as indicated in Geant
  Int_t            fPrecInteg;       // Alternative integration method, e.g. for higher precision
  Double_t         fFactorSol;       // Multiplicative factor for solenoid
  Double_t         fFactorDip;       // Multiplicative factor for dipole
  Double_t         fMax;             // Max Field as indicated in Geant
  Bool_t           fDipoleOFF;       // Dipole ON/OFF flag
  //
  Double_t         fQuadGradient;    // Gradient field for inner triplet quadrupoles
  Double_t         fDipoleField;     // Field value for D1 and D2 dipoles
  Double_t         fCCorrField;      // Side C 2nd compensator field
  Double_t         fACorr1Field;     // Side A 1st compensator field 
  Double_t         fACorr2Field;     // Side A 2nd compensator field
  //
  TNamed           fParNames;        // file and parameterization loadad
  //
  static const Double_t  fgkSol2DipZ;    // conventional Z of transition from L3 to Dipole field 
  static const Double_t  fgkBMachineZ1;  // Min Z of the LHC mag field range (to be checked?)
  static const Double_t  fgkBMachineZ2;  // Max Z of the LHC mag field range
  //   
  ClassDef(AliMagF, 1)           // Class for all Alice MagField wrapper for measured data + Tosca parameterization
};


#endif
