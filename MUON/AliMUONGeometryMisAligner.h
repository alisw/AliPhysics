/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup geometry
/// \class AliMUONGeometryMisAligner
/// \brief Class for misalignment of geometry transformations
///
/// Authors: Bruce Becker, Javier Castillo


#ifndef ALI_MUON_GEOMETRY_MIS_ALIGNER_H
#define ALI_MUON_GEOMETRY_MIS_ALIGNER_H

#include <TObject.h>

class AliMUONGeometryTransformer;

class TRandom;
class TGeoCombiTrans;

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
  
  void SetCartMisAlig(Double_t xmean, Double_t xwidth, Double_t ymean, Double_t ywidth)
    {fCartXMisAligM = xmean; fCartXMisAligW = xwidth; fCartYMisAligM = ymean; fCartYMisAligW = ywidth;}

  void SetCartMisAlig(Double_t mean, Double_t width)
    {fCartXMisAligM = mean; fCartXMisAligW = width; fCartYMisAligM = mean; fCartYMisAligW = width;}
  
  void SetAngMisAlig(Double_t mean, Double_t width)
    {fAngMisAligM = mean; fAngMisAligW = width;}
  
  void SetMaxCartMisAlig(Double_t width) // Kept for backward compatibility
    {fCartXMisAligM = 0.0; fCartXMisAligW = width; fCartYMisAligM = 0.0; fCartYMisAligW = width;}
  
  void SetMaxAngMisAlig(Double_t width) // Kept for backward compatibility
    {fAngMisAligM = 0.0; fAngMisAligW = width;}

  void SetXYAngMisAligFactor(Double_t factor);

  void SetZCartMisAligFactor(Double_t factor);

  void SetUseGaus(Bool_t usegaus)
    {fUseGaus=usegaus; fUseUni=!usegaus;}

  void SetUseUni(Bool_t useuni)
    {fUseGaus=!useuni; fUseUni=useuni;}
  
 protected:
  AliMUONGeometryMisAligner(const AliMUONGeometryMisAligner & right);
  AliMUONGeometryMisAligner & operator =(const AliMUONGeometryMisAligner &right);
  
  
 private:
  // return a misaligned transformation
  TGeoCombiTrans MisAlign(const TGeoCombiTrans& transform) const;
  void GetUniMisAlign(Double_t cartMisAlig[3], Double_t angMisAlig[3]) const;
  void GetGausMisAlign(Double_t cartMisAlig[3], Double_t angMisAlig[3]) const;

  Bool_t fUseUni;            // use uniform distribution for misaligmnets
  Bool_t fUseGaus;            // use gaussian distribution for misaligmnets
  Double_t fCartXMisAligM;   // cartesian displacement mean along x,  (translations)
  Double_t fCartXMisAligW;   // cartesian displacement width along x,  (translations)
  Double_t fCartYMisAligM;   // cartesian displacement mean along y,  (translations)
  Double_t fCartYMisAligW;   // cartesian displacement width along y,  (translations)
  Double_t fAngMisAligM;    // Angular displacement mean (rotations)
  Double_t fAngMisAligW;    // Angular displacement range (rotations)
  Double_t fXYAngMisAligFactor; // factor (<1) to apply to angular misalignment range since range of motion is restricted out of the xy plane
  Double_t fZCartMisAligFactor; // factor (<1) to apply to cartetian misalignment range since range of motion is restricted in z direction
  TRandom *fDisplacementGenerator;  // random number generator for the displacements
  
  ClassDef(AliMUONGeometryMisAligner,3)	// Geometry parametrisation
};

#endif //ALI_MUON_GEOMETRY_MIS_ALIGNER_H




