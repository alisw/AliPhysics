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

// --- Standard library ---


// --- ROOT system ---

#include <TObject.h>
#include <TMatrixFfwd.h>
#include <TVector3.h>

// --- AliRoot header files ---

class AliDigitNew;
class AliGeometry;

class AliRecPoint : public TObject {

public:

  AliRecPoint() ;                   // ctor            
  AliRecPoint(const char * opt) ;                   // ctor            
  AliRecPoint(const AliRecPoint &recp);  // copy ctor
  virtual ~AliRecPoint() ;          // dtor
 
  virtual void  AddDigit(AliDigitNew & digit) ;  // add a digit to the digit's indexes list  
  //  virtual void  AddTrack(AliTrack & track) ;  // add a track to the tracks list  
  virtual void  GetCovarianceMatrix(TMatrixF & mat) const;
  virtual AliGeometry * GetGeom() const { return fGeom; } 
  virtual void  GetGlobalPosition(TVector3 & gpos, TMatrixF & gmat) const ; // return global position in ALICE
  virtual int * GetDigitsList(void) const { return fDigitsList ; }
  //  virtual int * GetTracksList(void) const { return fTracksList ; }
  virtual Float_t GetEnergy() const {return fAmp; } 
  virtual void  GetLocalPosition(TVector3 & pos) const ;
  virtual Int_t GetDigitsMultiplicity(void) const { return fMulDigit ; }
  Int_t         GetIndexInList() const { return fIndexInList ; } 
  virtual Int_t GetMaximumDigitMultiplicity() const { return  fMaxDigit; } 
  virtual Int_t GetMaximumTrackMultiplicity() const { return  fMaxTrack; } 
  virtual Int_t GetTracksMultiplicity(void) const { return fMulTrack ; }
  AliRecPoint & operator= (const AliRecPoint &recp)
    {recp.Copy(*this); return (*this);}

  void          SetIndexInList(Int_t val) { fIndexInList = val ; } 


protected:
  void  Copy(TObject &recp) const;

  Float_t       fAmp ;        // summed amplitude of digits 
  AliGeometry * fGeom ;       //! pointer to the geometry class 
  Int_t         fIndexInList ;// the index of this RecPoint in the list stored in TreeR (to be set by analysis)
  TVector3      fLocPos ;     // local position in the sub-detector coordinate
  TMatrixF *     fLocPosM ;    // covariance matrix ;  
  Int_t         fMaxDigit ;   //! max initial size of digits array (not saved)
  Int_t         fMulDigit ;   // total multiplicity of digits
  Int_t         fMaxTrack ;   //! max initial size of tracks array (not saved)
  Int_t         fMulTrack ;   // total multiplicity of tracks
  Int_t *       fDigitsList ; //[fMulDigit] list of digit's indexes from which the point was reconstructed 
  Int_t *       fTracksList ; //[fMulTrack] list of tracks to which the point was assigned 

  ClassDef(AliRecPoint,1) // Base class for reconstructed space points
 
};

#endif // ALIRECPOINT_H
