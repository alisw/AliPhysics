#ifndef ALIEMCALRECPOINT_H
#define ALIEMCALRECPOINT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//_________________________________________________________________________
//  Base Class for EMCAL Reconstructed Points  
//  A recpoint being equivalent to a cluster in encal terminology                 
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

#include "AliRecPoint.h"
#include "AliEMCALDigit.h"

class AliEMCALRecPoint : public AliRecPoint {

 public:
  
  typedef TObjArray RecPointsList ; 

  AliEMCALRecPoint() ;                   // ctor         
  AliEMCALRecPoint(const char * opt) ;   // ctor 
  AliEMCALRecPoint(const AliEMCALRecPoint & rp):AliRecPoint(rp) { Fatal("cpy ctor", "not implemented") ; } 
  
  virtual ~AliEMCALRecPoint(){
    // dtor
  }
  virtual  void   AddDigit(AliDigitNew &){ Fatal("AddDigit", "use AddDigit(AliEMCALDigit & digit, Float_t Energy )") ; }
  virtual  void   AddDigit(AliEMCALDigit & digit, Float_t Energy) = 0 ; 
  virtual Int_t   Compare(const TObject * obj) const = 0 ;   
  virtual Int_t   DistancetoPrimitive(Int_t px, Int_t py);
  virtual void    Draw(Option_t * option="") ;
  virtual void    ExecuteEvent(Int_t event, Int_t, Int_t) ;
  virtual void    EvalAll(Float_t /*logWeight*/,TClonesArray * digits) ;  
  virtual void    EvalEMCALArm(AliEMCALDigit * digit) ;  
  virtual void    EvalPrimaries(TClonesArray * digits) ;  
  virtual Int_t   GetEMCALArm(void) const {return fEMCALArm ; }
  virtual void    GetGlobalPosition(TVector3 & /*gpos*/, TMatrix & /*gmat*/) const {;} // return global position in ALICE
  virtual void    GetGlobalPosition(TVector3 & gpos) const ; // return global position (r, theta, phi) in ALICE
  virtual void    GetLocalPosition(TVector3 & lpos) const ; // return loca position (x, y, z) in EMCAL
  //  virtual Int_t   GetEMCALMod(void) const {return fEMCALMod ; }
  virtual Int_t * GetPrimaries(Int_t & number) const {number = fMulTrack ; 
                                                      return fTracksList ; }
  virtual Bool_t  IsEmc(void)const { return kTRUE ;  }
  const Bool_t IsInECA(void) const { return fECASection ; } 
  const Bool_t IsInHCA(void) const { return fHCASection ; } 
  const Bool_t IsInPRE(void) const { return fPRESection ; } 
  virtual Bool_t  IsSortable() const { 
    // tells that this is a sortable object
    return kTRUE ; 
  }  
  virtual void    Paint(Option_t * option="");
  virtual void    Print(Option_t * /*opt = "void"*/) const {
    // Print prototype
  } 
  
  void SetECA() { fECASection = kTRUE ; } 
  void SetHCA() { fHCASection = kTRUE ; } 
  void SetPRE()  { fPRESection  = kTRUE ; } 
  AliEMCALRecPoint & operator = (const AliEMCALRecPoint & )  {
    Fatal("operator =", "not implemented") ;
    return *this ; 
  }

protected:
  
  Int_t fEMCALArm ; // EMCAM Arm number
  Float_t fTheta ; // theta angle in Alice
  Float_t fPhi ;   // phi angle in Alice
  Bool_t  fECASection ; // tells if the recpoint is in ECAL section 
  Bool_t  fHCASection ; // tells if the recpoint is in HCAL section 
  Bool_t  fPRESection ;  // tells if the recpoint is in PRE section 

  ClassDef(AliEMCALRecPoint,3) // RecPoint for EMCAL (Base Class)
 
};

#endif // AliEMCALRECPOINT_H
