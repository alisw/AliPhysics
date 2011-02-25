#ifndef ALIANALYSISTASKEMCALPI0CALIBSELECTION_H
#define ALIANALYSISTASKEMCALPI0CALIBSELECTION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//---------------------------------------------------------------------------// 
// Fill histograms with two-cluster invariant mass                           //
// using calibration coefficients of the previous iteration.                 //
//---------------------------------------------------------------------------//

// Root includes
class TH1F;
#include "TH2I.h"
#include "TObjArray.h"

// AliRoot includes
#include "AliAnalysisTaskSE.h"
class AliEMCALGeometry;
//class AliEMCALCalibData ;
#include "AliEMCALGeoParams.h"
class AliEMCALRecoUtils;

class AliAnalysisTaskEMCALPi0CalibSelection : public AliAnalysisTaskSE
{
public:

  AliAnalysisTaskEMCALPi0CalibSelection(const char* name);
  virtual ~AliAnalysisTaskEMCALPi0CalibSelection();

private:
  
  AliAnalysisTaskEMCALPi0CalibSelection(const AliAnalysisTaskEMCALPi0CalibSelection&); 
  AliAnalysisTaskEMCALPi0CalibSelection& operator=(const AliAnalysisTaskEMCALPi0CalibSelection&); 
  
public:
  
  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t * opt);
  virtual void LocalInit() ;
  
  void SetAsymmetryCut(Float_t asy)      {fAsyCut      = asy ;}
  void SetClusterMinEnergy(Float_t emin) {fEmin        = emin;}
  void SetClusterMaxEnergy(Float_t emax) {fEmax        = emax;}
  void SetClusterMinNCells(Int_t n)      {fMinNCells   = n   ;}
  void SetNCellsGroup(Int_t n)           {fGroupNCells = n   ;}
  void SetLogWeight(Float_t w)           {fLogWeight   = w   ;}
  
  //void SetCalibCorrections(AliEMCALCalibData* const cdata);
	
  void SwitchOnClusterCorrection()    {fCorrectClusters = kTRUE  ; }
  void SwitchOffClusterCorrection()   {fCorrectClusters = kFALSE ; }
  
  void SwitchOnSameSM()    {fSameSM = kTRUE  ; }
  void SwitchOffSameSM()   {fSameSM = kFALSE ; }
  
  Int_t  GetEMCALClusters(AliVEvent* event, TRefArray *clusters) const;
  Bool_t IsEMCALCluster(AliVCluster *clus) const;
  void SwitchOnOldAODs()   {fOldAOD = kTRUE  ; }
  void SwitchOffOldAODs()  {fOldAOD = kFALSE ; }  
  
  void SetGeometryName(TString name)                  { fEMCALGeoName = name   ; }
  TString GeometryName() const                        { return fEMCALGeoName   ; }
  void SwitchOnLoadOwnGeometryMatrices()              { fLoadMatrices = kTRUE  ; }
  void SwitchOffLoadOwnGeometryMatrices()             { fLoadMatrices = kFALSE ; }
  void SetGeometryMatrixInSM(TGeoHMatrix* m, Int_t i) { fMatrix[i]    = m      ; }

  void SetEMCALRecoUtils(AliEMCALRecoUtils * ru) {fRecoUtils = ru;}
  AliEMCALRecoUtils* GetEMCALRecoUtils() const   {return fRecoUtils;}
  
  void SetInvariantMassHistoBinRange(Int_t nBins, Float_t minbin, Float_t maxbin){
	fNbins = nBins; fMinBin = minbin; fMaxBin = maxbin; }
	  
  void GetMaxEnergyCellPosAndClusterPos(AliVCaloCells* cells, AliVCluster* clu, Int_t& iSM, Int_t& ieta, Int_t& iphi);

  void UseFilteredEventAsInput() {fFilteredInput = kTRUE;}
  void UseNormalEventAsInput()   {fFilteredInput = kFALSE;}

  
  void PrintInfo();
  
private:

  AliEMCALGeometry * fEMCALGeo;  //! EMCAL geometry
  //AliEMCALCalibData* fCalibData; // corrections to CC from the previous iteration
	
  Float_t fEmin;           // min. cluster energy
  Float_t fEmax;           // max. cluster energy
  Float_t fAsyCut;         // Asymmetry cut
  Int_t   fMinNCells;      // min. ncells in cluster
  Int_t   fGroupNCells;    // group n cells
  Float_t fLogWeight;      // log weight used in cluster recalibration
  Bool_t  fSameSM;         // Combine clusters in channels on same SM
  Bool_t  fOldAOD;         // Reading Old AODs, created before release 4.20
  Bool_t  fFilteredInput;  // Read input produced with filter.
  Bool_t  fCorrectClusters;// Correct clusters energy, position etc.
  TString fEMCALGeoName;   // Name of geometry to use.

  AliEMCALRecoUtils * fRecoUtils;  // Access to reconstruction utilities
  
  //Output histograms	
  Int_t   fNbins;  // N       mass bins of invariant mass histograms
  Float_t fMinBin; // Minimum mass bins of invariant mass histograms
  Float_t fMaxBin; // Maximum mass bins of invariant mass histograms

  TList*  fOutputContainer; //!histogram container
  TH1F*   fHmpi0[AliEMCALGeoParams::fgkEMCALModules][AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows];//! two-cluster inv. mass assigned to each cell.

  TH2F*   fHmgg;             //! two-cluster inv.mass vs pt of pair
  TH2F*   fHmggDifferentSM;  //! two-cluster inv.mass vs pt of pair, each cluster in different SM
  TH2F*   fHmggSM[AliEMCALGeoParams::fgkEMCALModules];       //! two-cluster inv.mass per SM
  TH2F*   fHmggPairSameSectorSM[AliEMCALGeoParams::fgkEMCALModules/2];   //! two-cluster inv.mass per Pair
  TH2F*   fHmggPairSameSideSM  [AliEMCALGeoParams::fgkEMCALModules-2];   //! two-cluster inv.mass per Pair

  TH2F*   fHOpeningAngle;             //! two-cluster opening angle vs pt of pair, with mass close to pi0
  TH2F*   fHOpeningAngleDifferentSM;  //! two-cluster opening angle vs pt of pair, each cluster in different SM, with mass close to pi0
  TH2F*   fHOpeningAngleSM[AliEMCALGeoParams::fgkEMCALModules];       //! two-cluster opening angle vs pt per SM,with mass close to pi0
  TH2F*   fHOpeningAnglePairSM[AliEMCALGeoParams::fgkEMCALModules];   //! two-cluster opening angle vs pt per Pair,with mass close to pi0

  TH2F*   fHIncidentAngle;             //! cluster incident angle vs pt of pair, with mass close to pi0
  TH2F*   fHIncidentAngleDifferentSM;  //! cluster incident angle vs pt of pair, each cluster in different SM, with mass close to pi0
  TH2F*   fHIncidentAngleSM[AliEMCALGeoParams::fgkEMCALModules];       //! cluster incident angle vs pt per SM,with mass close to pi0
  TH2F*   fHIncidentAnglePairSM[AliEMCALGeoParams::fgkEMCALModules];   //! cluster incident angle vs pt per Pair,with mass close to pi0
  
  TH2F*   fHAsymmetry;             //! two-cluster asymmetry vs pt of pair, with mass close to pi0
  TH2F*   fHAsymmetryDifferentSM;  //! two-cluster asymmetry vs pt of pair, each cluster in different SM, with mass close to pi0
  TH2F*   fHAsymmetrySM[AliEMCALGeoParams::fgkEMCALModules];       //! two-cluster asymmetry vs pt per SM,with mass close to pi0
  TH2F*   fHAsymmetryPairSM[AliEMCALGeoParams::fgkEMCALModules];   //! two-cluster asymmetry vs pt per Pair,with mass close to pi0
  
  TH2F*   fhTowerDecayPhotonHit[AliEMCALGeoParams::fgkEMCALModules] ;       //! Cells ordered in column/row for different module, number of times a decay photon hits
  TH2F*   fhTowerDecayPhotonEnergy[AliEMCALGeoParams::fgkEMCALModules] ;    //! Cells ordered in column/row for different module, accumulated energy in the tower by decay photons.
  TH2F*   fhTowerDecayPhotonAsymmetry[AliEMCALGeoParams::fgkEMCALModules] ; //! Cells ordered in column/row for different module, accumulated asymmetry in the tower by decay photons.

  TH1I*         fhNEvents;     //! Number of events counter histogram
  TList *       fCuts ;        //! List with analysis cuts
  Bool_t        fLoadMatrices; // Matrices set from configuration, not get from geometry.root or from ESDs/AODs
  TGeoHMatrix * fMatrix[AliEMCALGeoParams::fgkEMCALModules];    // Geometry matrices with alignments
  
  ClassDef(AliAnalysisTaskEMCALPi0CalibSelection,12);

};

#endif //ALIANALYSISTASKEMCALPI0CALIBSELECTION_H
