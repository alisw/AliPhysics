#ifndef ALIPMDRECPOINT_H
#define ALIPMDRECPOINT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  PMD Reconstructed Point                   //
//  Version 0.1                               //
//                                            //  
//                                            //
////////////////////////////////////////////////

#include "AliRecPoint.h"

class AliPMDRecPoint : public AliRecPoint {

public:
  virtual void  AddDigit(AliDigitNew & digit) ;  // add a digit to the digit's indexes list  
  //  virtual void  AddTrack(AliTrack & track) ;  // add a track to the tracks list  
  void  Copy(AliPMDRecPoint &recp) const;
  virtual void  GetCovarianceMatrix(TMatrix & mat) ;
  virtual AliGeometry * GetGeom() const { return fGeom; } 
  virtual void  GetGlobalPosition(TVector3 & gpos, TMatrix & gmat) const ; // return global position in ALICE
  virtual int * GetDigitsList(void) const { return fDigitsList ; }
  //  virtual int * GetTracksList(void) const { return fTracksList ; }
  virtual Float_t GetEnergy() const {return fAmp; } 
  virtual void  GetLocalPosition(TVector3 & pos) const ;
  virtual Int_t GetDigitsMultiplicity(void) const { return fMulDigit ; }
  Int_t         GetIndexInList() const { return fIndexInList ; } 
  virtual Int_t GetMaximumDigitMultiplicity() const { return  fMaxDigit; } 
  virtual Int_t GetMaximumTrackMultiplicity() const { return  fMaxTrack; } 
  virtual Int_t GetTracksMultiplicity(void) const { return fMulTrack ; }
  virtual void  Print(Option_t * opt = "void") {;}
  
  AliPMDRecPoint & operator= (const AliPMDRecPoint &recp);
  void          SetIndexInList(Int_t val) { fIndexInList = val ; } 
//
  ClassDef(AliPMDRecPoint,1) // Base class for reconstructed space points
 
};

#endif // ALIPMDRECPOINT_H
