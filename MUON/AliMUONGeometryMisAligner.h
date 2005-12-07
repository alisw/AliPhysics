/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup geometry
/// \class AliMUONGeometryMisAligner
/// \brief Class for misalignment of geometry transformations
///
/// Author: Bruce Becker


#ifndef ALI_MUON_GEOMETRY_MIS_ALIGNER_H
#define ALI_MUON_GEOMETRY_MIS_ALIGNER_H

#include <TObject.h>
#include <TGeoMatrix.h>
class TRandom;

class AliMUONGeometryTransformer;

class AliMUONGeometryMisAligner:public TObject
{
 public:
  AliMUONGeometryMisAligner(Double_t cartMisAlig, Double_t angMisAlig);
  AliMUONGeometryMisAligner();
  virtual ~AliMUONGeometryMisAligner();
  
  //_________________________________________________________________
  // methods
  
  // return a misaligned geometry obtained from the existing one.
  AliMUONGeometryTransformer* MisAlign(const AliMUONGeometryTransformer* transformer, Bool_t verbose = kFALSE);
  
  // return a misaligned transformation
  TGeoCombiTrans MisAlign(const TGeoCombiTrans& transform) const;

  void SetMaxCartMisAlig(Double_t offset)
    {fMaxCartMisAlig = offset ;}
  
  void SetMaxAngMisAlig(Double_t offset)
    {fMaxAngMisAlig = offset;}
  
  void SetXYAngMisAligFactor(Double_t factor);
  
  
 protected:
  AliMUONGeometryMisAligner(const AliMUONGeometryMisAligner & right);
  AliMUONGeometryMisAligner & operator =(const AliMUONGeometryMisAligner &right);
  
 private:
  Double_t fMaxCartMisAlig;   // cartesian displacement range, set by SetMaxCartMisAlig (translations)
  Double_t fMaxAngMisAlig;    // Angular displacement range (rotations)
  Double_t fXYAngMisAligFactor; // factor (<1) to apply to angular misalignment range since range of motion is restricted out of the xy plane
  TRandom *fDisplacementGenerator;  // random number generator for the displacements
  
  ClassDef(AliMUONGeometryMisAligner,2)	// Geometry parametrisation
};

#endif //ALI_MUON_GEOMETRY_MIS_ALIGNER_H




