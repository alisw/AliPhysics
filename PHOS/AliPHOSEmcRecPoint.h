#ifndef ALIPHOSEMCRECPOINT_H
#define ALIPHOSEMCRECPOINT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.36  2007/04/05 10:18:58  policheh
 * Introduced distance to nearest bad crystal.
 *
 * Revision 1.35  2007/03/06 06:47:28  kharlov
 * DP:Possibility to use actual vertex position added
 *
 * Revision 1.34  2005/05/28 14:19:04  schutz
 * Compilation warnings fixed by T.P.
 *
 */

//_________________________________________________________________________
//  RecPoint implementation for PHOS-EMC 
//  An EmcRecPoint is a cluster of digits   
//           
//-- Author: Dmitri Peressounko (RRC KI & SUBATECH)

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

//#include "AliPHOSDigit.h"
#include "AliPHOSRecPoint.h"

class AliPHOSEmcRecPoint : public AliPHOSRecPoint  {

public:

  AliPHOSEmcRecPoint() ;
  AliPHOSEmcRecPoint(const char * opt) ;
  AliPHOSEmcRecPoint(const AliPHOSEmcRecPoint & rp) ; 
 
  virtual ~AliPHOSEmcRecPoint() ;  

  virtual void Clear(const Option_t* /*option*/ ="") { delete[] fEnergyList; fEnergyList=0; AliPHOSRecPoint::Clear(); }

  //This virtual function has signature different from AliPHOSRecPoint::AddDigit
  //it hides, not overrides. using - declaration should fix the problem, at least for
  //g++
//  using AliPHOSRecPoint::AddDigit;
  virtual void  AddDigit(AliPHOSDigit & digit, Float_t Energy, Float_t time=0.) ;          // add a digit to the digits list  
  Int_t       Compare(const TObject * obj) const;                         // method for sorting  

  virtual void  EvalAll(TClonesArray * digits) ; //Those tasks which can be done without vertex
  virtual void  EvalAll(Float_t logWeight, TVector3 &vtx, TClonesArray * digits) ;
          void  EvalCoreEnergy(Float_t logWeight, Float_t coreRadius, TClonesArray * digits) ;             

  //in base class this functions is non-const
  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py) /*const*/; 

  Float_t         GetCoreEnergy()const {return fCoreEnergy ;}
  virtual Float_t GetDispersion()const {return fDispersion ;}
  virtual void    GetElipsAxis(Float_t * lambda)const { lambda[0] = fLambda[0] ;
                                                        lambda[1] = fLambda[1] ; }
  Float_t *   GetEnergiesList() const {return fEnergyList ;}       // gets the list of energies making this recpoint
  Float_t     GetMaximalEnergy(void) const ;                       // get the highest energy in the cluster
  Int_t       GetMaximumMultiplicity() const {return fMaxDigit ;}  // gets the maximum number of digits allowed
  Int_t       GetMultiplicity(void) const { return fMulDigit ; }   // gets the number of digits making this recpoint
  Int_t       GetMultiplicityAtLevel(Float_t level) const ;  // computes multiplicity of digits with 
                                                                   // energy above relative level
  Short_t     GetNExMax(void) const {return fNExMax ;}             // Number of maxima found in cluster in unfolding:
                                                                   // 0: was no unfolging
                                                                   //-1: unfolding failed
  void        SetNExMax(Int_t nmax=1){fNExMax = static_cast<Short_t>(nmax) ;}
  virtual Int_t GetNumberOfLocalMax(AliPHOSDigit **  maxAt, Float_t * maxAtEnergy,
                                    Float_t locMaxCut,TClonesArray * digits ) const ; 
                                                                   // searches for the local maxima 
  //returns number of local maxima in parent cluster or -2 if unfolding failed
  Float_t     GetTime(void) const{return  fTime ; } 
  Bool_t      IsEmc(void) const { return kTRUE ; }                 // true if the recpoint is in EMC
  Bool_t      IsSortable() const {return kTRUE ; }                 // says that emcrecpoints are sortable objects 
  void        Print(Option_t *)const ; 
  void        Purify(Float_t threshold) ;                          //Removes digits below threshold

  Float_t     GetM2x()   const {return fM2x;  } // Get second X-moment
  Float_t     GetM2z()   const {return fM2z;  } // Get second Z-moment
  Float_t     GetM3x()   const {return fM3x;  } // Get third  X-moment
  Float_t     GetM4z()   const {return fM4z;  } // Get forth  Z-moment
  Float_t     GetPhixe() const {return fPhixe;} // Get angle between center gravity and eigen vector

  Float_t     GetDistanceToBadCrystal() const {return fDistToBadCrystal;}
  void        SetDistanceToBadCrystal(Float_t dist) {fDistToBadCrystal=dist;}

  AliPHOSEmcRecPoint & operator = (const AliPHOSEmcRecPoint & /*rvalue*/)  {
    Fatal("operator =", "not implemented");
    return *this ;
  }

 protected:
  virtual void  EvalLocalPosition(Float_t logWeight, TVector3 &vtx, TClonesArray * digits, TVector3 &vInc) ;// computes the position in the PHOS module 
  virtual void  EvalDispersion(Float_t logWeight, TClonesArray * digits, TVector3 &vInc) ;   // computes the dispersion of the shower
  virtual void  EvalElipsAxis(Float_t logWeight, TClonesArray * digits, TVector3 &vInc );   // computes the axis of shower ellipsoide
          void  EvalMoments(Float_t logWeight, TClonesArray * digits, TVector3 &vInc );     // computes shower moments
  virtual void  EvalPrimaries(TClonesArray * digits) ;
          void  EvalTime( TClonesArray * digits );
  virtual Bool_t AreNeighbours(AliPHOSDigit * digit1, AliPHOSDigit * digit2 ) const ;

  Float_t fCoreEnergy ;       // energy in a shower core 
  Float_t fLambda[2] ;        // shower ellipse axes
  Float_t fDispersion ;       // shower dispersion
  Float_t *fEnergyList ;      //[fMulDigit] energy of digits
  Float_t fTime ;             // Time of the digit with maximal energy deposition
  Short_t fNExMax ;           // number of (Ex-)maxima before unfolding

  Float_t fM2x;               // Second moment along X axis
  Float_t fM2z;               // Second moment along Z axis
  Float_t fM3x;               // Third  moment along X axis
  Float_t fM4z;               // Forth  moment along Z axis
  Float_t fPhixe;             // Angle between center-gravity vector and eigen vector
  Float_t fDistToBadCrystal;  // Distance to nearest bad crystal

  Int_t fDebug;               //! debug level (0 - no output)
  
  ClassDef(AliPHOSEmcRecPoint,3)  // EMC RecPoint (cluster)

};

#endif // AliPHOSEMCRECPOINT_H
