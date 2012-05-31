#ifndef ALI_TPC_CALIB_GLOBAL_MISALIGNMENT_H
#define ALI_TPC_CALIB_GLOBAL_MISALIGNMENT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// AliTPCCalibGlobalMisalignment class                                        //
// The class calculates the space point distortions due to simple         // 
// misalignments like shifts in caresian coordinates or a rotation        //
// of the TPC read out planes (A and C side)                              //
//                                                                        //
// date: 06/05/2010                                                       //
// Authors: Stefan Rossegger, Jim Thomas, Magnus Mager                    //
////////////////////////////////////////////////////////////////////////////

#include "AliTPCCorrection.h"
#include "TVectorD.h"
class TGeoMatrix;
class TObjArray;
class TTreeSRedirector;


class AliTPCCalibGlobalMisalignment : public AliTPCCorrection {
public:
  AliTPCCalibGlobalMisalignment();
  virtual ~AliTPCCalibGlobalMisalignment();
  
  // initialization and update functions
  //  virtual void Init();
  //  virtual void Update(const TTimeStamp &timeStamp);
  void AddAlign(const  AliTPCCalibGlobalMisalignment & add);
  // setters and getters for misalignments
  void SetXShift(Float_t xShift) {fXShift=xShift;}
  void SetYShift(Float_t yShift) {fYShift=yShift;}
  void SetZShift(Float_t zShift) {fZShift=zShift;}
  void SetRotPhiA(Float_t rotPhiA) {fRotPhiA=rotPhiA;}
  void SetRotPhiC(Float_t rotPhiC) {fRotPhiC=rotPhiC;}
  void SetdRPhiOffsetA(Float_t dRPhiOffsetA) {fdRPhiOffsetA=dRPhiOffsetA;}
  void SetdRPhiOffsetC(Float_t dRPhiOffsetC) {fdRPhiOffsetC=dRPhiOffsetC;}
  
  Float_t GetXShift() const {return fXShift;}
  Float_t GetYShift() const {return fYShift;}
  Float_t GetZShift() const {return fZShift;}
  Float_t GetRotPhiA() const {return fRotPhiA;}
  Float_t GetRotPhiC() const {return fRotPhiC;}
  Float_t GetdRPhiOffsetA() const {return fdRPhiOffsetA;}
  Float_t GetdRPhiOffsetC() const {return fdRPhiOffsetC;}
  virtual void Print(Option_t* option="") const;
  void SetQuadranAlign(const TVectorD *quadrantQ0, const TVectorD *quadrantRQ0, const TVectorD *quadrantQ1,const TVectorD *quadrantRQ1,  const TVectorD *quadrantQ2,  const TVectorD *quadrantRQ2);
  // 
  // Alignment manipulation using TGeoMatrix
  
  void SetAlignGlobal(const TGeoMatrix * matrixGlobal);
  void SetAlignGlobalDelta(const TGeoMatrix * matrixGlobalDelta);
  void SetAlignSectors(const TObjArray *arraySector);
  TGeoMatrix* GetAlignGlobal() const  {return fMatrixGlobal;}
  TGeoMatrix* GetAlignGlobalDelta() const  {return fMatrixGlobalDelta;}
  TObjArray * GetAlignSectors() const {return fArraySector;}
  //
  static AliTPCCalibGlobalMisalignment*  CreateOCDBAlign();
  static AliTPCCalibGlobalMisalignment*  CreateMeanAlign(const AliTPCCalibGlobalMisalignment *alignIn);
  static void DumpAlignment( AliTPCCalibGlobalMisalignment* align, TTreeSRedirector *pcstream, const char *name);
  //
protected:
  virtual void GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]);

private:
  Float_t fXShift;               // Shift in global X [cm]
  Float_t fYShift;               // Shift in global Y [cm]
  Float_t fZShift;               // Shift in global Z [cm]

  Float_t fRotPhiA;      // simple rotation of A side read-out plane around the Z axis [rad]
  Float_t fRotPhiC;      // simple rotation of C side read-out plane around the Z axis [rad]
  Float_t fdRPhiOffsetA;  // add a constant offset of dRPhi (or local Y) in [cm]: purely for calibration purposes!
  Float_t fdRPhiOffsetC;  // add a constant offset of dRPhi (or local Y) in [cm]: purely for calibration purposes!
  //
  // Quadrant alignment
  //
  TVectorD *fQuadrantQ0;   //OROC medium pads -delta ly+ - ly - shift (cm)
  TVectorD *fQuadrantRQ0;  //OROC medium pads -delta ly+ - ly - rotation (rad) 
  TVectorD *fQuadrantQ1;   //OROC long   pads -delta ly+ - ly - shift
  TVectorD *fQuadrantQ2;   //OROC long   pads -shift
  TVectorD *fQuadrantRQ1;  //OROC long   pads -delta ly+ - ly - rotation
  TVectorD *fQuadrantRQ2;  //OROC long   pads -rotation
  // 
  //
  // Global alignment - use native ROOT representation
  //
  TGeoMatrix * fMatrixGlobal; // global Alignment common
  TGeoMatrix * fMatrixGlobalDelta; // global Alignment common A side-C side
  TObjArray   * fArraySector; //  local Alignmnet Sector
  //
  AliTPCCalibGlobalMisalignment& operator=(const AliTPCCalibGlobalMisalignment&);
  AliTPCCalibGlobalMisalignment(const AliTPCCalibGlobalMisalignment&);
  ClassDef(AliTPCCalibGlobalMisalignment,3);
};

#endif
