#ifndef ALIRECPOINT_H
#define ALIRECPOINT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Base class for Reconstructed Point        //
//  Version 0.1                               //
//  Author Yves Schutz     SUBATECH           //
//                                            //  
//                                            //
////////////////////////////////////////////////

// --- ROOT system ---

#include "TObject.h"
#include "TVector3.h"
#include "TMatrix.h"

// --- Standard library ---

#include "assert.h"
#include "iostream.h"

// --- AliRoot header files ---

#include "AliDigitNew.h"
#include "AliGeometry.h"
//#include "AliTrack.h"


class AliGeometry;

class AliRecPoint : public TObject {

public:

  AliRecPoint() ;                   // ctor            
  virtual ~AliRecPoint() ;          // dtor
 
  virtual void  AddDigit(AliDigitNew & digit) ;  // add a digit to the digits list  
  //  virtual void  AddTrack(AliTrack & track) ;  // add a track to the tracks list  
  virtual void  GetCovarianceMatrix(TMatrix & mat) ;
  virtual AliGeometry * GetGeom() const { return fGeom; } 
  virtual void  GetGlobalPosition(TVector3 & gpos, TMatrix & gmat) ; // return the global position in ALICE
  virtual int * GetDigitsList(void) const { return fDigitsList ; }
  //  virtual int * GetTracksList(void) const { return fTracksList ; }
  virtual Float_t GetEnergy() {return fAmp; } 
  virtual void  GetLocalPosition(TVector3 & pos) ;
  virtual Int_t GetDigitsMultiplicity(void) const { return fMulDigit ; }
  virtual Int_t GetMaximumDigitMultiplicity() const { return  fMaxDigit; } 
  virtual Int_t GetMaximumTrackMultiplicity() const { return  fMaxTrack; } 
  virtual Int_t GetTracksMultiplicity(void) const { return fMulTrack ; }
  virtual void  Print(Option_t * opt = "void") = 0 ; 


protected:

  Float_t       fAmp ;        // summed amplitude of digits 
  int *         fDigitsList ; // list of digits from which the point was reconstructed 
  AliGeometry * fGeom ;       // pointer to the geometry class 
  TVector3      fLocPos ;     // local position in the sub-detector coordinate
  TMatrix *     fLocPosM ;    // covariance matrix ;  
  Int_t         fMaxDigit ;   //! max initial size of digits array (not saved)
  Int_t         fMulDigit ;   // total multiplicity of digits
  Int_t         fMaxTrack ;   //! max initial size of tracks array (not saved)
  Int_t         fMulTrack ;   // total multiplicity of tracks
  int *         fTracksList ; // list of tracks to which the point was assigned 
 
public:
  

  ClassDef(AliRecPoint,1)
 
};

#endif // ALIRECPOINT_H
