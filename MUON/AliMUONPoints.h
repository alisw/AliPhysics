#ifndef ALIMUONPOINTS_H
#define ALIMUONPOINTS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004
//
/// \ingroup base
/// \class AliMUONPoints
/// \brief Class to draw detector clusters (is PolyMarker3D)

class TMarker3DBox;
#include <TMatrixFfwd.h>

class AliMUONDigit;
class AliMUONHit;
#include "AliPoints.h"

class AliMUONPoints : public AliPoints 
{
public:
  AliMUONPoints();
  AliMUONPoints(Int_t npoints);
  virtual ~AliMUONPoints();
                        
			/// Return hit index
  Int_t                 GetHitIndex() const {return fHitIndex;}
  Int_t                 GetTrackIndex() const; // *MENU*
			/// Return digit index
  Int_t                 GetDigitIndex() const {return fDigitIndex;}
			/// Return associated 3D-marker
  TMarker3DBox         *GetMarker(Int_t i) const {return fMarker[i];}
  AliMUONHit           *GetHit() const;
  AliMUONDigit         *GetDigit() const;
  virtual void          InspectHit(); // *MENU*
  virtual void          DumpHit() const; // *MENU*
  virtual void          InspectDigit(); // *MENU*
  virtual void          DumpDigit() const; // *MENU*
                        /// Set hit index
  virtual void          SetHitIndex(Int_t hitindex) {fHitIndex = hitindex;}
                        /// Set track index
  virtual void          SetTrackIndex(Int_t trackindex) {fTrackIndex = trackindex;}
                        /// Set digit index
  virtual void          SetDigitIndex(Int_t digitindex) {fDigitIndex = digitindex;}
                        /// Set associated 3D-marker
  virtual void          Set3DMarker(Int_t i,TMarker3DBox *marker) {fMarker[i] = marker;}
                        /// Set matrix
  virtual void          SetMatrix(TMatrixF *matrix) {fMatrix = matrix;}
  
protected:
  AliMUONPoints(const AliMUONPoints& points);  
  AliMUONPoints& operator = (const AliMUONPoints& rhs);

   Int_t            fHitIndex;         ///< Link to hit number 
   Int_t            fTrackIndex;       ///< Link to track number 
   Int_t            fDigitIndex;       ///< Link to digit 
   TMarker3DBox    *fMarker[3];        ///< pointer to  associated 3D-marker
   TMatrixF        *fMatrix;           ///< test
   
  ClassDef(AliMUONPoints,1) //Class to draw detector clusters (is PolyMarker3D) for MUON
};
#endif


