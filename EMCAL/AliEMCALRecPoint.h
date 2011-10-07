#ifndef ALIEMCALRECPOINT_H
#define ALIEMCALRECPOINT_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//_________________________________________________________________________
//  Base Class for EMCAL Reconstructed Points  
//  A recpoint being equivalent to a cluster in EMCAL terminology
//  
//  
//*-- Author: Yves Schutz (SUBATECH)
//*-- Author: Dmitri Peressounko (RRC KI & SUBATECH)
//*-- Author: Heather Gray (LBL): merged AliEMCALRecPoint and AliEMCALTowerRecPoint 02/04

// --- ROOT system ---
#include <TVector3.h>
class TGeoManager;
class TGeoPhysicalNode;
class TPad;
class TPaveText;
class TGraph;
class Riostream;
// --- Standard library ---

// --- AliRoot header files ---

#include "AliCluster.h"
class AliEMCALDigit;
class AliDigitNew;
class AliEMCALGeometry;
class AliEMCALHit;
class AliCaloCalibPedestal;

class AliEMCALRecPoint : public AliCluster {

 public:
  
  typedef TObjArray RecPointsList ; 

  AliEMCALRecPoint() ;                   // ctor         
  AliEMCALRecPoint(const char * opt) ;   // ctor 
  AliEMCALRecPoint(const AliEMCALRecPoint & rp);

  AliEMCALRecPoint& operator= (const AliEMCALRecPoint &rp);

  virtual ~AliEMCALRecPoint();

  virtual void    AddDigit(AliEMCALDigit & digit, const Float_t energy, const Bool_t shared); 
  virtual Int_t   Compare(const TObject * obj) const;  
  virtual void    Draw(Option_t * option="") ;

  virtual void    SetClusterType(Int_t ver) { fClusterType = ver ; }
  virtual Int_t   GetClusterType()    const { return fClusterType; }

  virtual void    EvalAll           (Float_t logWeight, TClonesArray * digits, const Bool_t justClusters);
  virtual void    EvalLocalPosition (Float_t logWeight, TClonesArray * digits);
  virtual void    EvalGlobalPosition(Float_t logWeight, TClonesArray * digits);

  virtual void    EvalPrimaries(TClonesArray * digits) ;
  virtual void    EvalParents  (TClonesArray * digits) ;

  void            EvalLocal2TrackingCSTransform();
  void            EvalLocalPositionFit(Double_t deff, Double_t w0, Double_t phiSlope,TClonesArray * digits);
  Bool_t          EvalLocalPosition2(TClonesArray *digits, TArrayD &ed);
  Bool_t          EvalLocalPositionFromDigits(const Double_t esum, const Double_t deff, const Double_t w0, 
                                              TClonesArray *digits, TArrayD &ed, TVector3 &locPos);
  Bool_t          EvalLocalPositionFromDigits(TClonesArray *digits, TArrayD &ed, TVector3 &locPos);
  static  void    GetDeffW0(const Double_t esum, Double_t &deff,  Double_t &w0);

  virtual void    GetGlobalPosition(TVector3 & gpos) const; // return global position (x, y, z) in ALICE
  virtual void    GetLocalPosition (TVector3 & lpos) const; // return local position  (x, y, z) in EMCAL SM
  
  virtual Int_t * GetPrimaries(Int_t & number)       const { number = fMulTrack  ; 
                                                             return fTracksList  ; }
  virtual Int_t * GetParents  (Int_t & number)       const { number = fMulParent ; 
                                                             return fParentsList ; }

  virtual Int_t   GetDigitsMultiplicity(void)  const { return fMulDigit    ; }
  Int_t           GetIndexInList()             const { return fIndexInList ; }
  virtual int *   GetDigitsList(void)          const { return fDigitsList  ; }
  virtual Float_t GetEnergy()                  const { return fAmp         ; }
  Float_t         GetCoreEnergy()              const { return fCoreEnergy  ; }
  virtual Float_t GetDispersion()              const { return fDispersion  ; }
  virtual void    GetElipsAxis(Float_t * lambda) const {lambda[0] = fLambda[0]; lambda[1] = fLambda[1];};
  Float_t *       GetEnergiesList()            const { return fEnergyList  ; } // gets the list of energies making this recpoint
  Double_t        GetPointEnergy()             const;                          // gets point energy (sum of energy list)
  Float_t         GetMaximalEnergy(void)       const ;                         // get the highest energy in the cluster
  Int_t           GetMaximalEnergyIndex(void)  const ;                         // get the index of highest energy digit
  Int_t           GetMaximumMultiplicity()     const { return fMaxDigit    ; } // gets the maximum number of digits allowed
  Int_t           GetMultiplicity(void)        const { return fMulDigit    ; } // gets the number of digits making this recpoint
  Int_t           GetMultiplicityAtLevel(Float_t level) const ;                // computes multiplicity of digits with 
  Int_t *         GetAbsId()                   const { return fAbsIdList   ; }
  Int_t           GetAbsId(Int_t i)            const { if(i>=0 && i<fMulDigit)
                                                        return fAbsIdList[i]; 
                                                        else return -1     ; }
  Int_t           GetAbsIdMaxDigit()           const { return GetAbsId(fDigitIndMax) ; }
  Int_t           GetIndMaxDigit()             const { return fDigitIndMax ; }
  void            SetIndMaxDigit(const Int_t ind)    { fDigitIndMax = ind  ; }
  void            SetIndexInList(Int_t val)          { fIndexInList = val  ; }

  virtual Int_t   GetSuperModuleNumber(void)   const { return fSuperModuleNumber;}

  // energy above relative level
  virtual Int_t   GetNumberOfLocalMax(AliEMCALDigit **  maxAt, Float_t * maxAtEnergy,
                                      Float_t locMaxCut,TClonesArray * digits ) const ; 
                                                                   // searches for the local maxima 
  // Number of local maxima found in cluster in unfolding:
  // 0: no unfolding
  //-1: unfolding failed
  Short_t         GetNExMax(void)              const { return fNExMax       ; }  // Number of maxima found in cluster in unfolding
  void            SetNExMax(Int_t nmax=1)            { fNExMax = static_cast<Short_t>(nmax) ;}
	
  Int_t           GetPrimaryIndex()            const  ;
	
  Float_t         GetTime(void)                const { return  fTime        ; }
	
  Bool_t          SharedCluster(void)          const { return  fSharedCluster ; }
  void            SetSharedCluster(Bool_t s)         { fSharedCluster = s     ; }
	
  virtual Bool_t  IsEmc(void)                  const { return kTRUE         ; }
  virtual Bool_t  IsSortable()                 const { return kTRUE         ; }  
  virtual void    Paint(Option_t * option="");
  virtual void    Print(Option_t * option="") const ; 
  
  Double_t        TmaxInCm(const Double_t e=0.0, const Int_t key=0);

  Float_t         GetDistanceToBadTower() const {return fDistToBadTower;}
  void            EvalDistanceToBadChannels(AliCaloCalibPedestal* caloped);

protected:
	  void           EvalCoreEnergy(Float_t logWeight, TClonesArray * digits) ;             
	  virtual void   EvalDispersion(Float_t logWeight, TClonesArray * digits) ;  // computes the dispersion of the shower
	  virtual void   EvalElipsAxis (Float_t logWeight, TClonesArray * digits );  // computes the axis of shower ellipsoide
	  void           EvalTime( TClonesArray * digits );
	  virtual Bool_t AreNeighbours(AliEMCALDigit * digit1, AliEMCALDigit * digit2 ) const;
	  Float_t        ThetaToEta(Float_t arg) const;  //Converts Theta (Radians) to Eta(Radians)
	  Float_t        EtaToTheta(Float_t arg) const;  //Converts Eta (Radians) to Theta(Radians)

private:

	  AliEMCALGeometry* fGeomPtr;  //! Pointer to geometry for utilities

	  Float_t  fAmp ;              // summed amplitude of digits   
	  Int_t    fIndexInList ;      // the index of this RecPoint in the
                                 // list stored in TreeR (to be set by analysis)
	  TVector3 fGlobPos ;          // global position
	  TVector3 fLocPos ;           // local  position in the sub-detector coordinate
	  Int_t    fMaxDigit ;         //! max initial size of digits array (not saved)
	  Int_t    fMulDigit ;         // total multiplicity of digits       
	  Int_t    fMaxTrack ;         //! max initial size of tracks array (not saved)
	  Int_t    fMulTrack ;         // total multiplicity of tracks
	  Int_t   *fDigitsList ;       //[fMulDigit] list of digit's indexes from which the point was reconstructed
	  Int_t   *fTracksList ;       //[fMulTrack] list of tracks to which the point was assigned

	  Int_t    fClusterType;       // type of cluster stored: v1
	  Float_t  fCoreEnergy ;       // energy in a shower core 
	  Float_t  fLambda[2] ;        // shower ellipse axes
	  Float_t  fDispersion ;       // shower dispersion
	  Float_t *fEnergyList ;       //[fMulDigit] energy of digits
	  Int_t   *fAbsIdList;         //[fMulDigit] absId  of digits
	  Float_t  fTime ;             // Time of the digit with maximal energy deposition
	  Short_t  fNExMax ;           // number of (Ex-)maxima before unfolding
	  Float_t  fCoreRadius;        // The radius in which the core energy is evaluated
	  Float_t *fDETracksList ;     //[fMulTrack] list of tracks to which the point was assigned
	  Int_t    fMulParent;         // Multiplicity of the parents
	  Int_t    fMaxParent;         // Maximum number of parents allowed
	  Int_t   *fParentsList;       // [fMulParent] list of the parents of the digits
	  Float_t *fDEParentsList;     // [fMulParent] list of the parents of the digits
	  Int_t    fSuperModuleNumber; // number identifying supermodule containing recpoint, reference is cell with maximum energy.
	  Int_t    fDigitIndMax;       // Index of digit with max energy in array fAbsIdList
	  Float_t  fDistToBadTower;    // Distance to nearest bad tower
	  Bool_t   fSharedCluster;     // States if cluster is shared by 2 SuperModules in same phi rack (0,1), (2,3) ... (10,11).
	
  ClassDef(AliEMCALRecPoint,13) // RecPoint for EMCAL (Base Class)
 
};

#endif // AliEMCALRECPOINT_H
