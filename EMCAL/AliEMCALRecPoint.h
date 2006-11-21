#ifndef ALIEMCALRECPOINT_H
#define ALIEMCALRECPOINT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//_________________________________________________________________________
//  Base Class for EMCAL Reconstructed Points  
//  A recpoint being equivalent to a cluster in encal terminology                 
//*-- Author: Yves Schutz (SUBATECH)
//*-- Author: Dmitri Peressounko (RRC KI & SUBATECH)
//*-- Author: Heather Gray (LBL): merged AliEMCALRecPoint and AliEMCALTowerRecPoint 02/04

// --- ROOT system ---
class TVector3 ;  

// --- Standard library ---

// --- AliRoot header files ---

#include "AliRecPoint.h"
class AliEMCALDigit;
class AliEMCALGeometry;

class AliEMCALRecPoint : public AliRecPoint {

 public:
  
  typedef TObjArray RecPointsList ; 

  AliEMCALRecPoint() ;                   // ctor         
  AliEMCALRecPoint(const char * opt) ;   // ctor 
  AliEMCALRecPoint(const AliEMCALRecPoint & rp);
  
  virtual ~AliEMCALRecPoint();
  virtual void    AddDigit(AliDigitNew &){ Fatal("AddDigit", "use AddDigit(AliEMCALDigit & digit, Float_t Energy )") ; }
  virtual void    AddDigit(AliEMCALDigit & digit, const Float_t Energy); 
  virtual Int_t   Compare(const TObject * obj) const;   
  virtual Int_t   DistancetoPrimitive(Int_t px, Int_t py);
  virtual void    Draw(Option_t * option="") ;
  virtual void    ExecuteEvent(Int_t event, Int_t, Int_t) ;

  virtual void    SetClusterType(Int_t ver) { fClusterType = ver; }
  virtual Int_t   GetClusterType() const { return fClusterType; }

  virtual void    EvalAll(Float_t logWeight, TClonesArray * digits);
  virtual void    EvalLocalPosition(Float_t logWeight, TClonesArray * digits) ;
  //  void            EvalLocalPositionSimple(TClonesArray *digits); // ??
  virtual void    EvalPrimaries(TClonesArray * digits) ;
  virtual void    EvalParents(TClonesArray * digits) ;

  using AliRecPoint::GetGlobalPosition;
  virtual void    GetGlobalPosition(TVector3 & gpos) const; // return global position (x, y, z) in ALICE
  virtual void    GetLocalPosition(TVector3 & lpos) const;  // return local position  (x, y, z) in EMCAL SM
  virtual Int_t * GetPrimaries(Int_t & number) const {number = fMulTrack ; 
                                                      return fTracksList ; }
  virtual Int_t * GetParents(Int_t & number) const {number = fMulParent ; 
                                                      return fParentsList ; }
  Float_t         GetCoreEnergy()const {return fCoreEnergy ;}
  virtual Float_t GetDispersion()const {return fDispersion ;}
  virtual void    GetElipsAxis(Float_t * lambda)const {lambda[0] = fLambda[0]; lambda[1] = fLambda[1];};
  
  Float_t *   GetEnergiesList() const {return fEnergyList ;}       // gets the list of energies making this recpoint
  Double_t    GetPointEnergy() const;                              // gets point energy (sum of energy list)
  Float_t *   GetTimeList() const {return fTimeList ;}       // gets the list of digit times in this recpoint
  Float_t     GetMaximalEnergy(void) const ;                       // get the highest energy in the cluster
  Int_t       GetMaximumMultiplicity() const {return fMaxDigit ;}  // gets the maximum number of digits allowed
  Int_t       GetMultiplicity(void) const { return fMulDigit ; }   // gets the number of digits making this recpoint
  Int_t       GetMultiplicityAtLevel(Float_t level) const ;  // computes multiplicity of digits with 
  Int_t *     GetAbsId() const {return fAbsIdList;}
  Int_t       GetAbsId(int i) const {if(i>=0 && i<fMulDigit)return fAbsIdList[i]; else return -1;}
  // energy above relative level
  virtual Int_t GetNumberOfLocalMax(AliEMCALDigit **  maxAt, Float_t * maxAtEnergy,
                                    Float_t locMaxCut,TClonesArray * digits ) const ; 
                                                                   // searches for the local maxima 
  
  Int_t       GetPrimaryIndex() const  ;
  Float_t     GetTime(void) const{return  fTime ; }
 
  virtual Bool_t  IsEmc(void)const { return kTRUE ;  }
  virtual Bool_t  IsSortable() const { 
    // tells that this is a sortable object
    return kTRUE ; 
  }  
  virtual void    Paint(Option_t * option="");
  virtual void    Print(Option_t * option="") const ; 
  
  AliEMCALRecPoint & operator = (const AliEMCALRecPoint & )  {
    Fatal("operator =", "not implemented") ;
    return *this ; 
  }

protected:
          void  EvalCoreEnergy(Float_t logWeight,TClonesArray * digits) ;             
	  virtual void  EvalDispersion(Float_t logWeight,TClonesArray * digits) ;   // computes the dispersion of the shower
	  virtual void  EvalElipsAxis(Float_t logWeight, TClonesArray * digits );   // computes the axis of shower ellipsoide
          void  EvalTime( TClonesArray * digits );
	  virtual Bool_t AreNeighbours(AliEMCALDigit * digit1, AliEMCALDigit * digit2 ) const;
	  Float_t ThetaToEta(Float_t arg) const;  //Converts Theta (Radians) to Eta(Radians)
	  Float_t EtaToTheta(Float_t arg) const;  //Converts Eta (Radians) to Theta(Radians)

private:
          AliEMCALGeometry* fGeomPtr;  //! Pointer to geometry for utilities

          Int_t   fClusterType;    // type of cluster stored:
				   // pseudocluster or v1
	  Float_t fCoreEnergy ;       // energy in a shower core 
	  Float_t fLambda[2] ;        // shower ellipse axes
	  Float_t fDispersion ;       // shower dispersio
	  Float_t *fEnergyList ;      //[fMulDigit] energy of digits
	  Float_t *fTimeList ;        //[fMulDigit] time of digits
          Int_t   *fAbsIdList;        //[fMulDigit] absId  of digits
	  Float_t fTime ;             // Time of the digit with maximal energy deposition
	  Float_t fCoreRadius;        // The radius in which the core energy is evaluated
          Int_t fMulParent;           // Multiplicity of the parents
          Int_t fMaxParent;           // Maximum number of parents allowed
          Int_t * fParentsList;       // [fMulParent] list of the parents of the digits
          Int_t   fSuperModuleNumber; // number identifying supermodule containing recpoint

  ClassDef(AliEMCALRecPoint,7) // RecPoint for EMCAL (Base Class)
 
};

#endif // AliEMCALRECPOINT_H
