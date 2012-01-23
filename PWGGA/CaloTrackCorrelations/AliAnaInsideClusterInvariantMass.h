#ifndef ALIANAINSIDECLUSTERINVARIANTMASS_H
#define ALIANAINSIDECLUSTERINVARIANTMASS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
//
// Split clusters with some criteria and calculate invariant mass
// to identify them as pi0 or conversion
//
//
//-- Author: Gustavo Conesa (LPSC-Grenoble)  
//_________________________________________________________________________


// --- ROOT system ---
class TList ;
class TObjString;
class TLorentzVector;

// --- ANALYSIS system ---
#include "AliAnaCaloTrackCorrBaseClass.h"

class AliAnaInsideClusterInvariantMass : public AliAnaCaloTrackCorrBaseClass {

 public: 
  AliAnaInsideClusterInvariantMass() ; // default ctor
  virtual ~AliAnaInsideClusterInvariantMass() { ; } //virtual dtor

  Bool_t       AreNeighbours(const Int_t absId1, const Int_t absId2) const ;  
  
  TObjString * GetAnalysisCuts();
  
  TList      * GetCreateOutputObjects();
  
  TLorentzVector GetCellMomentum(const Int_t absId, Float_t energy, AliVCaloCells* cells) ;
  
  Int_t        GetNumberOfLocalMaxima(AliVCluster* cluster, AliVCaloCells* cells,
                                      Int_t *absIdList,     Float_t *maxEList)  ;
  
  void         Init();
  
  void         InitParameters();
     
  void         MakeAnalysisFillHistograms() ; 

  void         RecalibrateCellAmplitude(Float_t  & amp,  const Int_t absId);
      
  void         SetCalorimeter(TString & det)     { fCalorimeter   = det  ; }
  
  void         SetLocalMaximaCutE(Float_t cut)   { fLocMaxCutE     = cut ; }
  
  void         SetLocalMaximaCutEDiff(Float_t c) { fLocMaxCutEDiff = c   ; }
  
  void         SetM02Cut(Float_t cut)            { fM02Cut         = cut ; }

  void         SetMinNCells(Int_t cut)           { fMinNCells      = cut ; }

  void         SetPi0MassRange(Float_t min, Float_t max) { fMassPi0Min = min ; fMassPi0Max = max ; }
  void         SetEtaMassRange(Float_t min, Float_t max) { fMassEtaMin = min ; fMassEtaMax = max ; }
  void         SetConMassRange(Float_t min, Float_t max) { fMassConMin = min ; fMassConMax = max ; }
  
  void         SplitEnergy(const Int_t absId1, const Int_t absId2, 
                           AliVCluster *cluster, 
                           AliVCaloCells* cells,
                           Float_t & e1, Float_t & e2,
                           const Int_t nMax//, Int_t *absIdList, Float_t *maxEList,
);
  
  void         Print(const Option_t * opt) const;

  //For histograms
  enum mcTypes { kmcPhoton = 1, kmcConversion = 2, kmcPi0    = 3,  
                 kmcEta    = 4, kmcElectron   = 5, kmcHadron = 6 };

 private:
  
  TString      fCalorimeter ;       // Calorimeter where the gamma is searched
  Float_t      fLocMaxCutE;         // Local maxima cut must have more than this energy
  Float_t      fLocMaxCutEDiff;     // Local maxima cut, when aggregating cells, next can be a bit higher
  Float_t      fM02Cut    ;         // Study clusters with l0 larger than cut
  Int_t        fMinNCells ;         // Study clusters with ncells larger than cut
  Float_t      fMassEtaMin;         // Min Eta mass
  Float_t      fMassEtaMax;         // Max Eta mass  
  Float_t      fMassPi0Min;         // Min Pi0 mass
  Float_t      fMassPi0Max;         // Min Pi0 mass
  Float_t      fMassConMin;         // Min Conversions mass
  Float_t      fMassConMax;         // Min Conversions mass
  
  //Histograms
  
  TH2F       * fhMassNLocMax1[7] ;  //! Mass of 2 highest energy cells when 1 local max, 1-6 for different MC particle types 
  TH2F       * fhMassNLocMax2[7] ;  //! Mass of 2 cells local maxima vs E,  1-6 for different MC particle types
  TH2F       * fhMassNLocMaxN[7] ;  //! Mass of >2 cells local maxima vs E, 1-6 for different MC particle types

  TH2F       * fhNLocMax[7] ;       //! Number of maxima in cluster vs E, 1-6 for different MC particle types
  TH2F       * fhNLocMaxEMax[7] ;   //! Number of maxima in cluster vs E of each maxima, 1-6 for different MC particle types
  TH2F       * fhNLocMaxEFrac[7] ;  //! Number of maxima in cluster vs fraction of cluster E of each maxima, 1-6 for different MC particle types
  TH2F       * fhNLocMaxM02Cut[7];  //! Number of maxima in cluster vs E, 1-6 for different MC particle types, after SS cut

  TH2F       * fhM02NLocMax1[7] ;   //! M02 vs E for N max in cluster = 1, 1-6 for different MC particle types
  TH2F       * fhM02NLocMax2[7] ;   //! M02 vs E for N max in cluster = 2, 1-6 for different MC particle types
  TH2F       * fhM02NLocMaxN[7] ;   //! M02 vs E for N max in cluster > 2, 1-6 for different MC particle types

  TH2F       * fhNCellNLocMax1[7] ; //! n cells in cluster vs E for N max in cluster = 1, 1-6 for different MC particle types
  TH2F       * fhNCellNLocMax2[7] ; //! n cells in cluster vs E for N max in cluster = 2, 1-6 for different MC particle types
  TH2F       * fhNCellNLocMaxN[7] ; //! n cells in cluster vs E for N max in cluster > 2, 1-6 for different MC particle types
  
  TH2F       * fhM02Pi0LocMax1[7] ; //! M02 for Mass around pi0, N Local Maxima = 1
  TH2F       * fhM02EtaLocMax1[7] ; //! M02 for Mass around eta, N Local Maxima = 1
  TH2F       * fhM02ConLocMax1[7] ; //! M02 for Mass around close to 0, N Local Maxima = 1

  TH2F       * fhM02Pi0LocMax2[7] ; //! M02 for Mass around pi0, N Local Maxima = 2
  TH2F       * fhM02EtaLocMax2[7] ; //! M02 for Mass around eta, N Local Maxima = 2
  TH2F       * fhM02ConLocMax2[7] ; //! M02 for Mass around close to 0, N Local Maxima = 2
  
  TH2F       * fhM02Pi0LocMaxN[7] ; //! M02 for Mass around pi0, N Local Maxima > 2
  TH2F       * fhM02EtaLocMaxN[7] ; //! M02 for Mass around eta, N Local Maxima > 2
  TH2F       * fhM02ConLocMaxN[7] ; //! M02 for Mass around close to 0, N Local Maxima > 2
  
  AliAnaInsideClusterInvariantMass(              const AliAnaInsideClusterInvariantMass & g) ; // cpy ctor
  AliAnaInsideClusterInvariantMass & operator = (const AliAnaInsideClusterInvariantMass & g) ; // cpy assignment
  
  ClassDef(AliAnaInsideClusterInvariantMass,4)
  
} ;

#endif //ALIANAINSIDECLUSTERINVARIANTMASS_H



