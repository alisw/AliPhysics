#ifndef ALIPHOSRECPOINT_H
#define ALIPHOSRECPOINT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//_________________________________________________________________________
//  Base Class for PHOS Reconstructed Points  
//  A recpoint being equivalent to a cluster in encal terminology                 
//*-- Author: Gines Martinez (SUBATECH)

// --- ROOT system ---
#include <TVector3.h>

// --- AliRoot header files ---
#include "AliCluster.h"

class TClonesArray ;
class AliPHOSDigit ;
class AliDigitNew;
class TMAtrixF; 

class AliPHOSRecPoint : public AliCluster {

 public:
  
  typedef TObjArray RecPointsList ; 

  AliPHOSRecPoint() ;                   // ctor         
  AliPHOSRecPoint(const char * opt) ;   // ctor 

  AliPHOSRecPoint(const AliPHOSRecPoint &rp);
  AliPHOSRecPoint& operator= (const AliPHOSRecPoint &rp);

  
  virtual ~AliPHOSRecPoint();

//  virtual  void   AddDigit(AliDigitNew &){
//   Fatal("AddDigit", "use AddDigit(AliPHOSDigit & digit, Float_t Energy)") ; 
//  }
  virtual  void   AddDigit(AliPHOSDigit & digit, Float_t Energy, Float_t time=0) = 0 ; 
  virtual Int_t   Compare(const TObject * obj) const = 0 ;   
  virtual Int_t   DistancetoPrimitive(Int_t px, Int_t py);
  virtual void    Draw(Option_t * option="") ;
  virtual void    ExecuteEvent(Int_t event, Int_t px, Int_t py) ;
  virtual void    EvalAll(TClonesArray * digits) ;  
  void            EvalLocal2TrackingCSTransform();
  virtual void    EvalPHOSMod(AliPHOSDigit * digit) ;  
  virtual int *   GetDigitsList(void) const { return fDigitsList ; }  
  virtual Float_t GetEnergy() const {return fAmp; }
  virtual void    GetLocalPosition(TVector3 & pos) const ;   
  virtual void    GetGlobalPosition(TVector3 & gpos, TMatrixF & gmat) const ; // return global position in ALICE
  virtual Int_t   GetPHOSMod(void) const {return fPHOSMod ; }
  virtual Int_t * GetPrimaries(Int_t & number) const {number = fMulTrack ; 
                                                      return fTracksList ; }
  virtual Int_t   GetDigitsMultiplicity(void) const { return fMulDigit ; }
  Int_t           GetIndexInList() const { return fIndexInList ; }
  virtual Bool_t  IsEmc(void)const { return kTRUE ;  } 
  virtual Bool_t  IsSortable() const { 
    // tells that this is a sortable object
    return kTRUE ; 
  }  
  void            SetIndexInList(Int_t val) { fIndexInList = val ; }
  virtual void    Paint(Option_t * option="");
  virtual void    Print(Option_t *) const {
    // Print prototype
  } 

protected:
  
  Int_t     fPHOSMod ;    // PHOS Module number in which the RecPoint is found
  Int_t     fMulTrack ;   // total multiplicity of tracks to which the point was assigned
  Int_t     fMaxDigit ;   //! max initial size of digits array (not saved)
  Int_t     fMulDigit ;   // total multiplicity of digits
  Int_t     fMaxTrack ;   //! max initial size of tracks array (not saved)
  Int_t*    fDigitsList ; //[fMulDigit] list of digit's indexes from which the point was reconstructed 
  Int_t*    fTracksList ; //[fMulTrack] list of tracks to which the point was assigned 
  Float_t   fAmp ;        // summed amplitude of digits 
  Int_t     fIndexInList ;// the index of this RecPoint in the list stored in TreeR (to be set by analysis)  
  TVector3  fLocPos ;     // local position in the sub-detector coordinate


  ClassDef(AliPHOSRecPoint,3) // RecPoint for PHOS (Base Class)
 
};

#endif // AliPHOSRECPOINT_H
