#ifndef ALIANAINSIDECLUSTERINVARIANTMASS_H
#define ALIANAINSIDECLUSTERINVARIANTMASS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: AliAnaInsideClusterInvariantMass.h 27413 2008-07-18 13:28:12Z gconesab $ */

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
#include "AliAnaPartCorrBaseClass.h"

class AliAnaInsideClusterInvariantMass : public AliAnaPartCorrBaseClass {

 public: 
  AliAnaInsideClusterInvariantMass() ; // default ctor
  virtual ~AliAnaInsideClusterInvariantMass() { ; } //virtual dtor
 private:
  AliAnaInsideClusterInvariantMass(const AliAnaInsideClusterInvariantMass & g) ; // cpy ctor
  AliAnaInsideClusterInvariantMass & operator = (const AliAnaInsideClusterInvariantMass & g) ;//cpy assignment

 public:  
	
  Bool_t       AreNeighbours(const Int_t absId1, const Int_t absId2) const ;  
  
  TObjString * GetAnalysisCuts();
  
  TList      * GetCreateOutputObjects();
  
  TLorentzVector GetCellMomentum(const Int_t absId, AliVCaloCells* cells);
  
  Int_t        GetNumberOfLocalMaxima(AliVCluster* cluster, AliVCaloCells* cells,
                                      Int_t *absIdList,     Float_t *maxEList)  ;
  
  void         Init();
  
  void         InitParameters();
     
  void         MakeAnalysisFillHistograms() ; 

  void         RecalibrateCellAmplitude(Float_t  & amp,  const Int_t absId);
      
  void         SetCalorimeter(TString & det) { fCalorimeter = det  ; }
  
  void         SetM02Cut(Float_t cut)        { fM02Cut      = cut  ; }

  void         SetMinNCells(Float_t cut)     { fMinNCells   = cut  ; }

  void         SplitEnergy(const Int_t absId1, const Int_t absId2, AliVCaloCells* cells,
                           Float_t & e1, Float_t & e2 );
  
  void         Print(const Option_t * opt) const;

  //For histograms
  enum mcTypes { mcPhoton = 1, mcConversion = 2, mcPi0    = 3,  
                 mcEta    = 4, mcElectron   = 5, mcHadron = 6 };

 private:
  
  TString      fCalorimeter ;       // Calorimeter where the gamma is searched
  Float_t      fM02Cut  ;           // Study clusters with l0 larger than cut
  Float_t      fMinNCells ;         // Study clusters with ncells larger than cut
 
  //Histograms
  
  TH2F       * fhMassNLocMax1[7] ;  //! Mass of 2 highest energy cells when 1 local max, 1-6 for different MC particle types 
  TH2F       * fhMassNLocMax2[7] ;  //! Mass of 2 cells local maxima vs E,  1-6 for different MC particle types
  TH2F       * fhMassNLocMaxN[7] ;  //! Mass of >2 cells local maxima vs E, 1-6 for different MC particle types

  TH2F       * fhNLocMax[7] ;       //! Number of maxima in cluster vs E, 1-6 for different MC particle types
  TH2F       * fhNLocMaxM02Cut[7];  //! Number of maxima in cluster vs E, 1-6 for different MC particle types, after SS cut

  TH2F       * fhM02NLocMax1[7] ;   //! M02 vs E for N max in cluster = 1, 1-6 for different MC particle types
  TH2F       * fhM02NLocMax2[7] ;   //! M02 vs E for N max in cluster = 2, 1-6 for different MC particle types
  TH2F       * fhM02NLocMaxN[7] ;   //! M02 vs E for N max in cluster > 2, 1-6 for different MC particle types

  TH2F       * fhNCellNLocMax1[7] ; //! n cells in cluster vs E for N max in cluster = 1, 1-6 for different MC particle types
  TH2F       * fhNCellNLocMax2[7] ; //! n cells in cluster vs E for N max in cluster = 2, 1-6 for different MC particle types
  TH2F       * fhNCellNLocMaxN[7] ; //! n cells in cluster vs E for N max in cluster > 2, 1-6 for different MC particle types
  
  TH2F       * fhM02Pi0[7] ;       //! M02 for Mass around pi0
  TH2F       * fhM02Eta[7] ;       //! M02 for Mass around eta
  TH2F       * fhM02Con[7] ;       //! M02 for Mass around close to 0

  TH2F       * fhInvMassAllCells[7] ; //! Inv mass of all cells

  
  ClassDef(AliAnaInsideClusterInvariantMass,1)
  
} ;

#endif //ALIANAINSIDECLUSTERINVARIANTMASS_H



