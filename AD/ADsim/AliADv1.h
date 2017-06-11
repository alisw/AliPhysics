// -*- C++ -*-
#ifndef ALIADV1_H
#define ALIADV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//////////////////////////////////////////////////
//  Manager and hits classes for set : AD       //
//////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//                  AD (ALICE Diffractive)  Detector                     //
//                                                                       //
//  This class contains the base procedures for the AD  detector         //
//  New geometry of 2014                                                 //
//  All comments should be sent to :                                     //
//                                                                       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "AliAD.h"
#include <TMath.h>
// #include <TGeoCompositeShape.h> 
// #include <TGeoMatrix.h>

class TGeoVolumeAssembly;
class TGeoCompositeShape;
class TGeoRotation;
class TGeoVolume;
class TGeoMedium;
//
class CurvedBundle; 
class PmtFiberInfo;
//
enum {
  kCSideFiberNearWLSTop, kCSideFiberNearWLSBot, kCSideFiberPmtTop, kCSideFiberPmtBot, 
  kASideFiberShort, kASideFiberLong }; 
//

class AliADv1 : public AliAD {
public:
  AliADv1();
  AliADv1(const char *name, const char *title);
  virtual ~AliADv1();
  virtual    void  AddAlignableVolumes() const;

  virtual TString  Version() { return TString("v1"); }
  virtual   Int_t  IsVersion() const { return 1; }
  virtual    void  AddHit(Int_t track, Int_t *vol, Float_t *hits);
  // virtual    void  AddDigits(Int_t* track, Int_t module, Float_t time);
  virtual    void  MakeBranch(Option_t *option);
  virtual    void  CreateGeometry();
  virtual    void  Init();
  virtual    void  StepManager();
  virtual    void  DisableTunnelStruct() { fADCstruct = kFALSE; }
  virtual    void  DisableStructuresSideC() { fADCstruct = kFALSE; }
  virtual    void  DisableStructuresSideA() { fADAstruct = kFALSE; }
  virtual    void  KeepHistory() { fKeepHistory = kTRUE; }

  enum ADCPosition_t { kADCInTunnel, kADCInCavern, kADCInBoth};

protected:

  ///////////////////////////////////////////////////////////////////////////////////
  // 
  // Geometry update
  //
  TGeoRotation * Rx180 ;  //! 
  TGeoRotation * Rz180 ;  //! 
  TGeoRotation * Ry180 ;  //! 
  TGeoRotation * Rx90  ;  //! 
  TGeoRotation * Rx90m ;  //! 
  TGeoRotation * Ry90m ;  //! 
  TGeoRotation * Ry90  ;  //! 
  TGeoRotation * Rz90  ;  //! 
   
  // Position of Fiber bundles A-Side:
  Double_t fX1FiberShort[8]; //!
  Double_t fX1FiberLong [8]; //!
  Double_t fX2Fiber     [8]; //!
  Double_t fY1Fiber     [8]; //!
  Double_t fY2Fiber     [8]; //!
  // 
  static const char * s_where[];          //!
  static const Int_t fNBundles = 8;       //!
  CurvedBundle * fTopBundles[fNBundles];  //!
  CurvedBundle * fBotBundles[fNBundles];  //!
  PmtFiberInfo * fPmtFiberTop[4];  //!
  PmtFiberInfo * fPmtFiberBot[4];  //!
  //
  // Private comunication by Arturo Tauro (2014, Apr 23)
  // According to last survey measurement done this morning, 
  // the C-Side wall is at Z = - 18959mm.
  //
  static const Double_t kZwall        ;  // Aluminium plate z position 
  static const Double_t kZendAbs      ;  // End of CC block absorber
  static const Double_t kZbegVMAOI    ;  // Begining of Warm Module
  static const Double_t kZbegValve    ;  // Begining of Valve
  static const Double_t kZbegFrontBar ;  // Begining of Front Bar
  static const Double_t kZbegCoil     ;  // Begining of compensator coil
  Double_t             RadToDeg(Double_t rad) { return rad * (180./TMath::Pi()); }
  Double_t             DegToRad(Double_t deg) { return deg * (TMath::Pi()/180.); }
  TGeoVolume         * Make_UProfile(const char * volname, const Double_t L, 
                           const TGeoMedium * medium, const Double_t W, const Double_t H, 
                           const Double_t dw, const Double_t dh);
  TGeoVolume         * Make_UProfileH(const char * volname, TGeoMedium * medium);
  TGeoVolume         * MakeVolIBeam(const char * volname, const TGeoMedium * mat, 
                           const Double_t x, const Double_t y, 
                           const Double_t dx, const Double_t dy, const Double_t dz);
  TGeoVolume         * CreateWarmModuleSupport();
  TGeoVolume         * CreateOldADA();
  TGeoVolume         * CreatePump();
  TGeoVolume         * CreateSupportZEM();
  TGeoVolumeAssembly * CreateBLM();
  TGeoVolumeAssembly * CreateADAShielding();
  TGeoVolumeAssembly * CreateVacuumChamberSupport();
  // void                 CreateCurvedBundles(TGeoVolumeAssembly * ad=0);
  TGeoVolumeAssembly * CreatePmtBoxC() ;
  TGeoCompositeShape * MakeFiberBundle (const char * shName, const Double_t xWLS, const Double_t xPMT, const Double_t yWLS, const Double_t yPMT, const Double_t Z1);
  // Bool_t               FindInArray(Int_t & , const Int_t , const Int_t nel, const Int_t * arr);
  // Float_t              GetDistanceFromFibersToWLSADC(Int_t where, Int_t arrayindex, Double_t * x);
  // Float_t              GetDistanceFromFibersToWLSADA(Int_t ADmodule, Bool_t IsShort, Float_t fHitY);
  //
public:
  ///////////////////////////////////////////////////////////////////////////////////
  // functions for ADA and ADC
  void ReadADCFromEnv(void);
  TGeoCompositeShape * MakeShapeADCpadH(const Double_t W, const Double_t H, const Double_t dz);
  virtual    void  CreateAD();

private:
  // Position of ADC: In the Tunnel, In the Cavern, or in Both
  Bool_t      fADCstruct;
  Bool_t      fADAstruct;
  ADCPosition_t fADCPosition;

  //! AD Optical parameters :
  Double_t    fADCLightYield;       //! Lightyield in BC404
  Double_t    fADALightYield;       //! Lightyield in BC404
  Double_t    fADCPhotoCathodeEfficiency;
  Double_t    fADAPhotoCathodeEfficiency;

  Bool_t      fKeepHistory;

  AliADv1(const AliAD&);
  AliADv1& operator = (const AliADv1&);

  ClassDef(AliADv1, 1)  // Class for the AD detector
};

#endif
