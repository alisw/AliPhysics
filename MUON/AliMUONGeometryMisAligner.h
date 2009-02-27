/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup geometry
/// \class AliMUONGeometryMisAligner
/// \brief Class for misalignment of geometry transformations
//
//  Authors: Bruce Becker, Javier Castillo


#ifndef ALI_MUON_GEOMETRY_MIS_ALIGNER_H
#define ALI_MUON_GEOMETRY_MIS_ALIGNER_H

#include <TObject.h>

class AliMUONGeometryTransformer;

class TGeoCombiTrans;
class TClonesArray;

class AliMUONGeometryMisAligner:public TObject
{
 public:
  AliMUONGeometryMisAligner(Double_t cartXMisAligM, Double_t cartXMisAligW, Double_t cartYMisAligM, Double_t cartYMisAligW, Double_t angMisAligM, Double_t angMisAligW);
  AliMUONGeometryMisAligner(Double_t cartMisAligM, Double_t cartMisAligW, Double_t angMisAligM, Double_t angMisAligW);
  AliMUONGeometryMisAligner(Double_t cartMisAligW, Double_t angMisAligW);
  AliMUONGeometryMisAligner();
  virtual ~AliMUONGeometryMisAligner();
  
  //_________________________________________________________________
  // methods
  
  // return a misaligned geometry obtained from the existing one.
  AliMUONGeometryTransformer* MisAlign(const AliMUONGeometryTransformer* transformer, 
                                       Bool_t verbose = kFALSE);
  
  /// Set cartesian displacement parameters different along x, y
  void SetCartMisAlig(Double_t xmean, Double_t xwidth, Double_t ymean, Double_t ywidth, Double_t zmean = 0., Double_t zwidth = 0.)
    {fDetElemMisAlig[0][0] = xmean; fDetElemMisAlig[0][1] = xwidth; fDetElemMisAlig[1][0] = ymean; fDetElemMisAlig[1][1] = ywidth; fDetElemMisAlig[2][0] = zmean; fDetElemMisAlig[2][1] = zwidth; }

  /// Set cartesian displacement parameters, the same along x, y
  void SetCartMisAlig(Double_t mean, Double_t width)
    {fDetElemMisAlig[0][0] = mean; fDetElemMisAlig[0][1] = width; fDetElemMisAlig[1][0] = mean; fDetElemMisAlig[1][1] = width;}
  
  /// Set angular displacement
  void SetAngMisAlig(Double_t zmean, Double_t zwidth, Double_t xmean = 0., Double_t xwidth = 0., Double_t ymean = 0., Double_t ywidth = 0.)
    {fDetElemMisAlig[3][0] = xmean; fDetElemMisAlig[3][1] = xwidth; fDetElemMisAlig[4][0] = ymean; fDetElemMisAlig[4][1] = ywidth; fDetElemMisAlig[5][0] = zmean; fDetElemMisAlig[5][1] = zwidth;}
  
  /// Set cartesian displacement (Kept for backward compatibility)
  void SetMaxCartMisAlig(Double_t width) 
    {fDetElemMisAlig[0][0] = 0.0; fDetElemMisAlig[0][1] = width; fDetElemMisAlig[1][0] = 0.0; fDetElemMisAlig[1][1] = width;}
  
  /// Set angular displacement (Kept for backward compatibility)
  void SetMaxAngMisAlig(Double_t width) 
    {fDetElemMisAlig[5][0] = 0.0; fDetElemMisAlig[5][1] = width;}

  void SetXYAngMisAligFactor(Double_t factor);

  void SetZCartMisAligFactor(Double_t factor);

  /// Set option for gaussian distribution 
  void SetUseGaus(Bool_t usegaus)
    {fUseGaus=usegaus; fUseUni=!usegaus;}

  /// Set option for uniform distribution 
  void SetUseUni(Bool_t useuni)
    {fUseGaus=!useuni; fUseUni=useuni;}

  /// Set module (half chambers) cartesian displacement parameters
  void SetModuleCartMisAlig(Double_t xmean, Double_t xwidth, Double_t ymean, Double_t ywidth, Double_t zmean, Double_t zwidth) 
    {fModuleMisAlig[0][0] = xmean; fModuleMisAlig[0][1] = xwidth; fModuleMisAlig[1][0] = ymean; fModuleMisAlig[1][1] = ywidth; fModuleMisAlig[2][0] = zmean; fModuleMisAlig[2][1] = zwidth;}

  /// Set module (half chambers) cartesian displacement parameters
  void SetModuleAngMisAlig(Double_t xmean, Double_t xwidth, Double_t ymean, Double_t ywidth, Double_t zmean, Double_t zwidth) 
    {fModuleMisAlig[3][0] = xmean; fModuleMisAlig[3][1] = xwidth; fModuleMisAlig[4][0] = ymean; fModuleMisAlig[4][1] = ywidth; fModuleMisAlig[5][0] = zmean; fModuleMisAlig[5][1] = zwidth;}

  /// Set alignment resolution to misalign objects to be stored in CDB
  void SetAlignmentResolution(const TClonesArray* misAlignArray, Int_t chId=-1, Double_t chResX=-1., Double_t chResY=-1., Double_t deResX=-1., Double_t deResY=-1.);
  
 protected:
  /// Not implemented
  AliMUONGeometryMisAligner(const AliMUONGeometryMisAligner & right);
  /// Not implemented
  AliMUONGeometryMisAligner & operator =(const AliMUONGeometryMisAligner &right);
  
  
 private:
  // return a misaligned transformation
  TGeoCombiTrans MisAlignDetElem(const TGeoCombiTrans& transform) const;
  TGeoCombiTrans MisAlignModule(const TGeoCombiTrans& transform) const;
  void GetUniMisAlign(Double_t cartMisAlig[3], Double_t angMisAlig[3], const Double_t lParMisAlig[6][2]) const;
  void GetGausMisAlign(Double_t cartMisAlig[3], Double_t angMisAlig[3], const Double_t lParMisAlig[6][2]) const;

  Bool_t fUseUni;            ///< use uniform distribution for misaligmnets
  Bool_t fUseGaus;           ///< use gaussian distribution for misaligmnets
  Double_t fDetElemMisAlig[6][2]; ///< Mean and width of the displacements of the detection elements along x,y,z (translations) and about x,y,z (rotations)
  Double_t fModuleMisAlig[6][2];  ///< Mean and width of the displacements of the modules along x,y,z (translations) and about x,y,z (rotations)  

  Double_t fXYAngMisAligFactor;  ///< factor (<1) to apply to angular misalignment range since range of motion is restricted out of the xy plane
  Double_t fZCartMisAligFactor; ///< factor (<1) to apply to cartetian misalignment range since range of motion is restricted in z direction

  
  ClassDef(AliMUONGeometryMisAligner,4)	// Geometry parametrisation
};

#endif //ALI_MUON_GEOMETRY_MIS_ALIGNER_H




