#ifndef ALICALORIMETERUTILS_H
#define ALICALORIMETERUTILS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliCalorimeterUtils
/// \ingroup CaloTrackCorrelationsBase
/// \brief Class with utils specific to calorimeter clusters/cells.
///
/// Class containing utility methods for calorimeters. It performs calibration,
/// cluster splitting (EMCal), geometry initialization, OADB initialization, 
/// plus other goodies.
/// For EMCal, it relys heavyly in AliEMCALRecoUtils.
///
/// More information can be found in this [twiki](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations).
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
//_________________________________________________________________________

// --- ROOT system ---
#include <TObject.h> 
#include <TString.h>
#include <TObjArray.h>
class TArrayF;  
#include <TH2I.h>
#include <TGeoMatrix.h>
class AliMCEvent;

//--- ANALYSIS system ---
class AliVEvent;
class AliVTrack;
class AliAODPWG4Particle;
class AliAODCaloCluster;
class AliVCaloCells;
class AliPHOSGeoUtils;
class AliEMCALGeometry;
class AliAODMCParticle;
class TParticle;

#include "AliEMCALRecoUtils.h"

class AliCalorimeterUtils : public TObject {

 public:   
  AliCalorimeterUtils() ; // ctor
  virtual ~AliCalorimeterUtils() ;//virtual dtor
  
  virtual void  InitParameters();
  virtual void  Print(const Option_t * opt)          const ;

  virtual Int_t GetDebug()                           const { return fDebug                ; }
  virtual void  SetDebug(Int_t d)                          { fDebug = d                   ; }
	
  //virtual void Init();
	
  // Cluster contents
  
  Bool_t        AreNeighbours(Int_t calo, Int_t absId1, Int_t absId2) const ;

  Bool_t        IsClusterSharedByTwoSuperModules(const AliEMCALGeometry * geom,
                                                 AliVCluster* cluster);
  
  Bool_t        GetFECCorrelatedCellAbsId(Int_t absId, Int_t absIdCorr[4]) const ;
  
  Bool_t        IsAbsIDsFromTCard(Int_t absId1, Int_t absId2, 
                                  Int_t & rowDiff, Int_t & colDiff) const ;
  
  Int_t         GetNumberOfLocalMaxima(AliVCluster* cluster, AliVCaloCells* cells)  ;
  
  Int_t         GetNumberOfLocalMaxima(AliVCluster* cluster, AliVCaloCells* cells,
                                       Int_t *absIdList,     Float_t *maxEList)  ;
  
  Float_t       GetLocalMaximaCutE()                 const { return fLocMaxCutE           ; }
  void          SetLocalMaximaCutE(Float_t cut)            { fLocMaxCutE     = cut        ; }
  
  Float_t       GetLocalMaximaCutEDiff()             const { return fLocMaxCutEDiff       ; }
  void          SetLocalMaximaCutEDiff(Float_t c)          { fLocMaxCutEDiff = c          ; }
  
  Int_t         GetMaxEnergyCell(AliVCaloCells* cells, AliVCluster* clu, Float_t & fraction) const ;
  
  void          SplitEnergy(Int_t absId1, Int_t absId2, AliVCluster *cluster, AliVCaloCells* cells,
                            //Float_t & e1, Float_t & e2,
                            AliAODCaloCluster *cluster1, AliAODCaloCluster *cluster2,
                            Int_t nMax, Int_t eventNumber = 0);//, Int_t *absIdList, Float_t *maxEList,
  
  void          SwitchOnClusterPlot()                      { fPlotCluster = kTRUE         ; }
  void          SwitchOffClusterPlot()                     { fPlotCluster = kFALSE        ; }

  Float_t       GetMCECellClusFracCorrection(Float_t eCell, Float_t eCluster) const ;
  void          SetMCECellClusFracCorrectionParamters(Int_t i, Float_t param) { if(i<4) fMCECellClusFracCorrParam[i] = param; }
  
  Bool_t        IsMCECellClusFracCorrectionOn()   const    { return fMCECellClusFracCorrOn   ; }
  void          SwitchOnMCECellClusFracCorrection()        { fMCECellClusFracCorrOn = kTRUE  ; }
  void          SwitchOffMCECellClusFracCorrection()       { fMCECellClusFracCorrOn = kFALSE ; }

  //------------------------------
  // Calorimeters Geometry Methods
  //------------------------------
  
  AliEMCALGeometry * GetEMCALGeometry()              const { return fEMCALGeo             ; }
  TString       EMCALGeometryName()                  const { return fEMCALGeoName         ; }  
  void          SetEMCALGeometryName(TString name)         { fEMCALGeoName = name         ; }
  void          InitEMCALGeometry() ;
  Bool_t        IsEMCALGeoMatrixSet()                const { return fEMCALGeoMatrixSet    ; }
	
  AliPHOSGeoUtils * GetPHOSGeometry()                const { return fPHOSGeo              ; }	
  TString       PHOSGeometryName()                   const { return fPHOSGeoName          ; }  
  void          SetPHOSGeometryName(TString name)          { fPHOSGeoName = name          ; }
  void          InitPHOSGeometry() ;
  Bool_t        IsPHOSGeoMatrixSet()                 const { return fPHOSGeoMatrixSet     ; }

  void          AccessGeometry(AliVEvent* inputEvent) ;
	
  void          SetImportGeometryFromFile(Bool_t import,
                                          TString path = ""){
                                                             fImportGeometryFromFile = import    ;
                                                             fImportGeometryFilePath = path      ; } // EMCAL
  
  Bool_t        IsMCParticleInCalorimeterAcceptance(Int_t calo, TParticle* particle);
  Bool_t        IsMCParticleInCalorimeterAcceptance(Int_t calo, AliAODMCParticle* particle);
  Bool_t        IsMCParticleInCalorimeterAcceptance(Int_t calo, Float_t eta, Float_t theta, Float_t phi, Int_t & absID);
  
  void          SwitchOnLoadOwnEMCALGeometryMatrices()     { fLoadEMCALMatrices = kTRUE   ; }
  void          SwitchOffLoadOwnEMCALGeometryMatrices()    { fLoadEMCALMatrices = kFALSE  ; }
  void          SetEMCALGeometryMatrixInSM(TGeoHMatrix* m, Int_t i) { fEMCALMatrix[i] = m ; }
  
  void          SwitchOnLoadOwnPHOSGeometryMatrices()      { fLoadPHOSMatrices = kTRUE    ; }
  void          SwitchOffLoadOwnPHOSGeometryMatrices()     { fLoadPHOSMatrices = kFALSE   ; }
  void          SetPHOSGeometryMatrixInSM(TGeoHMatrix* m, Int_t i) { fPHOSMatrix[i] = m   ; }
  
  void         GetEMCALSubregion(AliVCluster* clus, AliVCaloCells* cells, 
                                 Int_t & regEta, Int_t & regPhi) const ;
  
  //------------------------------
  // Bad channels
  //------------------------------

  Bool_t        IsBadChannelsRemovalSwitchedOn()     const { return fRemoveBadChannels                                ; }
  void          SwitchOnBadChannelsRemoval ()              { fRemoveBadChannels = kTRUE   ; 
                                                             fEMCALRecoUtils->SwitchOnBadChannelsRemoval(); 
                                                             if(!fPHOSBadChannelMap) InitPHOSBadChannelStatusMap()    ; }
  void          SwitchOffBadChannelsRemoval()              { fRemoveBadChannels = kFALSE ; 
                                                             fEMCALRecoUtils->SwitchOffBadChannelsRemoval()           ; }
  
  Bool_t        IsDistanceToBadChannelRecalculated() const { return  IsDistanceToBadChannelRecalculated()             ; }
  void          SwitchOnDistToBadChannelRecalculation ()   { fEMCALRecoUtils->SwitchOnDistToBadChannelRecalculation() ; }
  void          SwitchOffDistToBadChannelRecalculation()   { fEMCALRecoUtils->SwitchOffDistToBadChannelRecalculation(); }
  
  void          InitPHOSBadChannelStatusMap () ;

  Int_t         GetEMCALChannelStatus(Int_t iSM , Int_t iCol, Int_t iRow) const { 
                  return fEMCALRecoUtils->GetEMCALChannelStatus(iSM,iCol,iRow); }//Channel is ok by default

  Int_t         GetPHOSChannelStatus (Int_t imod, Int_t iCol, Int_t iRow) const { 
                  if(fPHOSBadChannelMap) return (Int_t) ((TH2I*)fPHOSBadChannelMap->At(imod))->GetBinContent(iCol,iRow); 
                  else return 0 ; }//Channel is ok by default
  
  void          SetEMCALChannelStatus(Int_t iSM , Int_t iCol, Int_t iRow, Double_t c = 1) { 
                  fEMCALRecoUtils->SetEMCALChannelStatus(iSM,iCol,iRow,c) ; }
  
  void          SetPHOSChannelStatus (Int_t imod, Int_t iCol, Int_t iRow, Double_t c = 1) {
                  if(!fPHOSBadChannelMap) InitPHOSBadChannelStatusMap() ; 
                  ((TH2I*)fPHOSBadChannelMap->At(imod))->SetBinContent(iCol,iRow,c) ; }
    
  void          SetEMCALChannelStatusMap(Int_t iSM , TH2I* h) { fEMCALRecoUtils->SetEMCALChannelStatusMap(iSM,h)      ; }
  void          SetPHOSChannelStatusMap(Int_t imod , TH2I* h) { fPHOSBadChannelMap ->AddAt(h,imod)                    ; }
  
  TH2I *        GetEMCALChannelStatusMap(Int_t iSM)  const { return fEMCALRecoUtils->GetEMCALChannelStatusMap(iSM)    ; }
  TH2I *        GetPHOSChannelStatusMap(Int_t imod)  const { return (TH2I*)fPHOSBadChannelMap->At(imod)               ; }

  void          SetEMCALChannelStatusMap(TObjArray *map)   { fEMCALRecoUtils->SetEMCALChannelStatusMap(map)           ; }
  void          SetPHOSChannelStatusMap (TObjArray *map)   { fPHOSBadChannelMap  = map                                ; }
	
  Bool_t        ClusterContainsBadChannel(Int_t calo,UShort_t* cellList, Int_t nCells);
  Bool_t        ClusterContainsBadChannel(TString /*calo*/,UShort_t* /*cellList*/, Int_t /*nCells*/) {return kFALSE;} // Stupid thing to do but just to avoid compilation break in AliTrackComparisonESD while I find the authors
	
  // Mask clusters in front of frame, EMCAL only
  Int_t         GetNMaskCellColumns()                const { return fNMaskCellColumns;}
  void          SetNMaskCellColumns(Int_t n) {
                  if(n > fNMaskCellColumns) { delete [] fMaskCellColumns ; fMaskCellColumns = new Int_t[n] ; }
                  fNMaskCellColumns = n                                                                               ; }
  void          SetMaskCellColumn(Int_t ipos, Int_t icol) { 
                  if(ipos < fNMaskCellColumns) fMaskCellColumns[ipos] = icol;
                  else printf("Not set, position larger than allocated set size first")                               ; }
  Bool_t        MaskFrameCluster(Int_t iSM, Int_t ieta) const ;
  
  //------------------------------
  // Calorimeter indexes information
  //------------------------------

  Int_t         GetModuleNumber(AliAODPWG4Particle * particle, AliVEvent* inputEvent) const;
  Int_t         GetModuleNumber(AliVCluster * cluster) const;
  Int_t         GetModuleNumberCellIndexes(Int_t absId, Int_t calo, Int_t & icol, Int_t & irow, Int_t &iRCU) const ;
  Int_t         GetModuleNumberCellIndexesAbsCaloMap(Int_t absId, Int_t calo, Int_t & icol, Int_t & irow, Int_t &iRCU, 
                                                     Int_t & icolAbs, Int_t & irowAbs) const ;
	
  //------------------------------
  // Modules fiducial region
  //------------------------------

  Bool_t        CheckCellFiducialRegion(AliVCluster* cluster, AliVCaloCells* cells) const ;
  Bool_t        CheckCellFiducialRegion(AliVCluster* cluster, AliVCaloCells* cells, AliVEvent* /**/, Int_t /**/)
                { return CheckCellFiducialRegion(cluster, cells) ; } // Stupid thing to do but just to avoid compilation break in AliTrackComparisonESD while I find the authors
  void          SetNumberOfCellsFromPHOSBorder(Int_t n)    { fNCellsFromPHOSBorder = n                                ; }
  Int_t         GetNumberOfCellsFromPHOSBorder()     const { return fNCellsFromPHOSBorder                             ; }
  void          SetNumberOfCellsFromEMCALBorder(Int_t n)   { fEMCALRecoUtils->SetNumberOfCellsFromEMCALBorder(n)      ; }
  Int_t         GetNumberOfCellsFromEMCALBorder()    const { return fEMCALRecoUtils->GetNumberOfCellsFromEMCALBorder(); }
  void          SwitchOnNoFiducialBorderInEMCALEta0()      { fEMCALRecoUtils->SwitchOnNoFiducialBorderInEMCALEta0()   ; }
  void          SwitchOffNoFiducialBorderInEMCALEta0()     { fEMCALRecoUtils->SwitchOffNoFiducialBorderInEMCALEta0()  ; }
  Bool_t        IsEMCALNoBorderAtEta0()              const { return fEMCALRecoUtils->IsEMCALNoBorderAtEta0()          ; }
  
  //------------------------------
  // Recalibration
  //------------------------------

  Bool_t        IsRecalibrationOn()                  const { return fRecalibration                                    ; }
  void          SwitchOnRecalibration()                    { fRecalibration = kTRUE ; 
                  InitPHOSRecalibrationFactors(); fEMCALRecoUtils->SwitchOnRecalibration()                            ; }
  void          SwitchOffRecalibration()                   { fRecalibration = kFALSE;
                  fEMCALRecoUtils->SwitchOffRecalibration()                                                           ; }
	
  void          InitPHOSRecalibrationFactors () ;
	
  Float_t       GetEMCALChannelRecalibrationFactor(Int_t iSM , Int_t iCol, Int_t iRow) const { 
                  return fEMCALRecoUtils->GetEMCALChannelRecalibrationFactor(iSM , iCol, iRow)                        ; }
  
  Float_t       GetPHOSChannelRecalibrationFactor (Int_t imod, Int_t iCol, Int_t iRow) const { 
                  if(fPHOSRecalibrationFactors)
                    return (Float_t) ((TH2F*)fPHOSRecalibrationFactors->At(imod))->GetBinContent(iCol,iRow); 
                  else return 1                                                                                       ; }
  
  void          SetEMCALChannelRecalibrationFactor(Int_t iSM , Int_t iCol, Int_t iRow, Double_t c = 1) { 
                  fEMCALRecoUtils->SetEMCALChannelRecalibrationFactor(iSM,iCol,iRow,c)                                ; }
	
  void          SetPHOSChannelRecalibrationFactor (Int_t imod, Int_t iCol, Int_t iRow, Double_t c = 1) {
                  if(!fPHOSRecalibrationFactors)  InitPHOSRecalibrationFactors();
                  ((TH2F*)fPHOSRecalibrationFactors->At(imod))->SetBinContent(iCol,iRow,c)                            ; }
    
  void          SetEMCALChannelRecalibrationFactors(Int_t iSM , TH2F* h) { fEMCALRecoUtils->SetEMCALChannelRecalibrationFactors(iSM,h)      ; }
  void          SetPHOSChannelRecalibrationFactors(Int_t imod , TH2F* h) { fPHOSRecalibrationFactors ->AddAt(h,imod)                        ; }
	
  TH2F *        GetEMCALChannelRecalibrationFactors(Int_t iSM)     const { return fEMCALRecoUtils->GetEMCALChannelRecalibrationFactors(iSM) ; }
  TH2F *        GetPHOSChannelRecalibrationFactors(Int_t imod)     const { return (TH2F*)fPHOSRecalibrationFactors->At(imod)                ; }
	
  void          SetEMCALChannelRecalibrationFactors(TObjArray *map)      { fEMCALRecoUtils->SetEMCALChannelRecalibrationFactors(map)        ; }
  void          SetPHOSChannelRecalibrationFactors (TObjArray *map)      { fPHOSRecalibrationFactors  = map;}

  void          RecalibrateCellTime     (Double_t & time, Int_t calo, Int_t absId, Int_t bunchCrossNumber) const ;
  void          RecalibrateCellTimeL1Phase(Double_t & time, Int_t calo, Int_t iSM, Int_t bunchCrossNumber) const;
  void          RecalibrateCellAmplitude(Float_t  & amp,  Int_t calo, Int_t absId) const ;
  Float_t       RecalibrateClusterEnergy(AliVCluster* cluster, AliVCaloCells * cells);
  Float_t       RecalibrateClusterEnergyWeightCell(AliVCluster* cluster, AliVCaloCells * cells, Float_t energyOrg);

  //------------------------------
  // Run dependent energy calibrations (EMCAL)
  //------------------------------

  void          SwitchOffRunDepCorrection()                              {  fRunDependentCorrection = kFALSE  ; }
  void          SwitchOnRunDepCorrection()                               {  fRunDependentCorrection = kTRUE   ; }
  
  //------------------------------
  // Time Recalibration (EMCAL)
  //------------------------------

  Bool_t       IsTimeRecalibrationOn()                             const { return fEMCALRecoUtils->IsTimeRecalibrationOn() ; }
  void         SwitchOffTimeRecalibration()                              { fEMCALRecoUtils->SwitchOffTimeRecalibration()   ; }
  void         SwitchOnTimeRecalibration()                               { fEMCALRecoUtils->SwitchOnTimeRecalibration()    ; }
  
  Float_t      GetEMCALChannelTimeRecalibrationFactor(Int_t bc, Int_t absID) const
  { return fEMCALRecoUtils->GetEMCALChannelTimeRecalibrationFactor(bc, absID) ; } 
	
  void         SetEMCALChannelTimeRecalibrationFactor(Int_t bc, Int_t absID, Double_t c = 0)
  { fEMCALRecoUtils->SetEMCALChannelTimeRecalibrationFactor(bc, absID, c) ; }  
  
  TH1F *       GetEMCALChannelTimeRecalibrationFactors(Int_t bc) const     { return fEMCALRecoUtils-> GetEMCALChannelTimeRecalibrationFactors(bc) ; }
  void         SetEMCALChannelTimeRecalibrationFactors(TObjArray *map)     { fEMCALRecoUtils->SetEMCALChannelTimeRecalibrationFactors(map)        ; }
  void         SetEMCALChannelTimeRecalibrationFactors(Int_t bc , TH1F* h) { fEMCALRecoUtils->SetEMCALChannelTimeRecalibrationFactors(bc , h)     ; }

  //------------------------------
  // Time Recalibration - L1 phase (EMCAL)
  //------------------------------
  Bool_t   IsL1PhaseInTimeRecalibrationOn()          const { return fEMCALRecoUtils->IsL1PhaseInTimeRecalibrationOn() ; }
  void     SwitchOffL1PhaseInTimeRecalibration()           { fEMCALRecoUtils->SwitchOffL1PhaseInTimeRecalibration()   ; }
  void     SwitchOnL1PhaseInTimeRecalibration()            { fEMCALRecoUtils->SwitchOnL1PhaseInTimeRecalibration()    ; }

  Int_t    GetEMCALL1PhaseInTimeRecalibrationForSM(Int_t iSM) const        { return fEMCALRecoUtils->GetEMCALL1PhaseInTimeRecalibrationForSM(iSM)   ; }
  void     SetEMCALL1PhaseInTimeRecalibrationForSM(Int_t iSM, Int_t c = 0) { return fEMCALRecoUtils->SetEMCALL1PhaseInTimeRecalibrationForSM(iSM,c) ; }
  
  TH1C *   GetEMCALL1PhaseInTimeRecalibrationForAllSM()const          { return fEMCALRecoUtils->GetEMCALL1PhaseInTimeRecalibrationForAllSM() ; }
  void     SetEMCALL1PhaseInTimeRecalibrationForAllSM(TObjArray *map) { fEMCALRecoUtils->SetEMCALL1PhaseInTimeRecalibrationForAllSM(map)     ; }
  void     SetEMCALL1PhaseInTimeRecalibrationForAllSM(TH1C* h)        { fEMCALRecoUtils->SetEMCALL1PhaseInTimeRecalibrationForAllSM(h)       ; }

  //------------------------------
  // EMCAL specific utils for the moment
  //------------------------------

  void          SetEMCALRecoUtils(AliEMCALRecoUtils * ru)  { fEMCALRecoUtils = ru          ; }
  AliEMCALRecoUtils* GetEMCALRecoUtils()             const { return fEMCALRecoUtils        ; }
  
  Bool_t        IsCorrectionOfClusterEnergyOn()      const { return fCorrectELinearity     ; }
  void          SwitchOnCorrectClusterLinearity()          { fCorrectELinearity = kTRUE    ; } 
  void          SwitchOffCorrectClusterLinearity()         { fCorrectELinearity = kFALSE   ; } 
  void          CorrectClusterEnergy(AliVCluster *cl);
  
  Bool_t        IsRecalculationOfClusterPositionOn() const { return fRecalculatePosition   ; }
  void          SwitchOnRecalculateClusterPosition()       { fRecalculatePosition = kTRUE  ; } 
  void          SwitchOffRecalculateClusterPosition()      { fRecalculatePosition = kFALSE ; } 
  void          RecalculateClusterPosition(AliVCaloCells* cells, AliVCluster* clu);
  void          RecalculateClusterShowerShapeParameters(AliVCaloCells* cells, AliVCluster* clu){
                  fEMCALRecoUtils->RecalculateClusterShowerShapeParameters((AliEMCALGeometry*)fEMCALGeo, cells, clu)  ; }
  
  void          RecalculateClusterDistanceToBadChannel(AliVCaloCells* cells, AliVCluster* clu){
                  fEMCALRecoUtils->RecalculateClusterDistanceToBadChannel((AliEMCALGeometry*)fEMCALGeo, cells, clu)   ; }
  
  void          RecalculateClusterPID(AliVCluster* clu)    { fEMCALRecoUtils->RecalculateClusterPID(clu)              ; }
  
  //------------------------------
  // *** Track Matching ***
  //------------------------------

  AliVTrack *   GetMatchedTrack(AliVCluster * cluster, AliVEvent * event, Int_t index = -1) const ;
  
  // Recalculation

  void          RecalculateClusterTrackMatching(AliVEvent * event, 
                                                TObjArray* clusterArray = 0x0,
                                                AliMCEvent* mc = 0x0) ;
    
  void          GetMatchedResiduals(Int_t index, Float_t &dR, Float_t &dZ) {
                  if (fRecalculateMatching) fEMCALRecoUtils->GetMatchedResiduals(index,dR,dZ)   ; }
  
  //------------------------------
  // This could be used for PHOS ...
  //------------------------------

  void          SwitchOnRecalculateClusterTrackMatching()       { fRecalculateMatching = kTRUE  ; } 
  void          SwitchOffRecalculateClusterTrackMatching()      { fRecalculateMatching = kFALSE ; } 
  Bool_t        IsRecalculationOfClusterTrackMatchingOn() const { return fRecalculateMatching   ; }
  
  Float_t       GetCutZ()                                 const { return fCutZ                  ; } // PHOS only
  void          SetCutZ(Float_t z)                              { fCutZ = z                     ; } // PHOS only

  
  Float_t       GetCutR()                                 const { return fCutR                  ; } // PHOS and EMCAL
  void          SetCutR(Float_t r)                              { fCutR = r                     ;   // PHOS and EMCA
                                                                  fEMCALRecoUtils->SetCutR(r)   ; }
  
  Float_t       GetCutEta()                               const { return fCutEta                ; } // EMCAL only
  void          SetCutEta(Float_t e)                            { fCutEta = e                   ;   // EMCAL only
                                                                  fEMCALRecoUtils->SetCutEta(e) ; }

  Float_t       GetCutPhi()                               const { return fCutPhi                ; } // EMCAL only
  void          SetCutPhi(Float_t p)                            { fCutPhi = p                   ;   // EMCAL only
                                                                  fEMCALRecoUtils->SetCutPhi(p) ; }
  
  //------------------------------
  // OADB options settings
  //------------------------------

  void          AccessOADB(AliVEvent * event) ;
  
  TString       GetPass() ;
  
  void          SwitchOnEMCALOADB()                             { fOADBForEMCAL = kTRUE         ; }
  void          SwitchOffEMCALOADB()                            { fOADBForEMCAL = kFALSE        ; }

  void          SwitchOnPHOSOADB()                              { fOADBForPHOS  = kTRUE         ; }
  void          SwitchOffPHOSOADB()                             { fOADBForPHOS  = kFALSE        ; }

  void          SetEMCALOADBFilePath(TString path)              { fOADBFilePathEMCAL  = path    ; }
  void          SetPHOSOADBFilePath (TString path)              { fOADBFilePathPHOS   = path    ; }

  //------------------------------
  // Other settings
  //------------------------------

  void          SetNumberOfSuperModulesUsed(Int_t nSM)          { fNSuperModulesUsed  = nSM     ; }
  Int_t         GetNumberOfSuperModulesUsed()             const { return fNSuperModulesUsed     ; }

  void          SetRunNumber(Int_t run)                         { fRunNumber  = run             ; }
  Int_t         GetRunNumber()                            const { return fRunNumber             ; }
  
 private:

  Int_t              fDebug;                    ///<  Debugging level.

  TString            fEMCALGeoName;             ///<  Name of geometry to use for EMCAL.
  
  TString            fPHOSGeoName;              ///<  Name of geometry to use for PHOS.	

  AliEMCALGeometry * fEMCALGeo ;                //!<! EMCAL geometry pointer.

  AliPHOSGeoUtils  * fPHOSGeo  ;                //!<! PHOS  geometry pointer. 

  Bool_t             fEMCALGeoMatrixSet;        ///<  Check if the transformation matrix is set for EMCAL.

  Bool_t             fPHOSGeoMatrixSet ;        ///<  Check if the transformation matrix is set for PHOS.

  Bool_t             fLoadEMCALMatrices;        ///<  Matrices set from configuration, not get from geometry.root or from ESDs/AODs.

  TGeoHMatrix *      fEMCALMatrix[22];          ///<  Geometry matrices with alignments.

  Bool_t             fLoadPHOSMatrices;         ///<  Matrices set from configuration, not get from geometry.root or from ESDs/AODs.

  TGeoHMatrix *      fPHOSMatrix[5];            ///<  Geometry matrices with alignments.

  Bool_t             fRemoveBadChannels;        ///<  Check the channel status provided and remove clusters with bad channels.

  TObjArray        * fPHOSBadChannelMap;        ///<  Array of histograms with map of bad channels, PHOS.

  Int_t              fNCellsFromPHOSBorder;     ///<  Number of cells from PHOS  border the cell with maximum amplitude has to be.

  Int_t              fNMaskCellColumns;         ///<  Number of masked columns.
  
  /// List of masked cells collumn index.
  Int_t            * fMaskCellColumns;          //[fNMaskCellColumns] 

  Bool_t             fRecalibration;            ///<  Switch on or off the recalibration.

  Bool_t             fRunDependentCorrection;   ///<  Switch on or off the recalibration dependent on T.

  TObjArray        * fPHOSRecalibrationFactors; ///<  Array of histograms with map of recalibration factors, PHOS.

  AliEMCALRecoUtils* fEMCALRecoUtils;           ///<  EMCAL utils for cluster rereconstruction.

  Bool_t             fRecalculatePosition;      ///<  Recalculate cluster position.
 
  Bool_t             fCorrectELinearity  ;      ///<  Correct cluster energy linearity.
  
  Bool_t             fRecalculateMatching;      ///<  Recalculate cluster position.
  
  Float_t            fCutR;                     ///<  dR cut on matching (PHOS).
  
  Float_t            fCutZ;                     ///<  dZ cut on matching (EMCAL/PHOS).
  
  Float_t            fCutEta;                   ///<  dEta cut on matching (EMCAL).
  
  Float_t            fCutPhi;                   ///<  dPhi cut on matching (EMCAL).
  
  Float_t            fLocMaxCutE;               ///<  Local maxima cut must have more than this energy.
  
  Float_t            fLocMaxCutEDiff;           ///<  Local maxima cut, when aggregating cells, next can be a bit higher.
  
  Bool_t             fPlotCluster;              ///<  Plot cluster in splitting method.
  
  Bool_t             fOADBSet ;                 ///<  AODB parameters already set.
  
  Bool_t             fOADBForEMCAL ;            ///<  Get calibration from OADB for EMCAL.
  
  Bool_t             fOADBForPHOS ;             ///<  Get calibration from OADB for PHOS.
  
  TString            fOADBFilePathEMCAL ;       ///<  Default path $ALICE_PHYSICS/OADB/EMCAL, if needed change.
  
  TString            fOADBFilePathPHOS ;        ///<  Default path $ALICE_PHYSICS/OADB/PHOS, if needed change.
  
  Bool_t             fImportGeometryFromFile;   ///<  Import geometry settings in geometry.root file.
  
  TString            fImportGeometryFilePath;   ///<  Path fo geometry.root file.

  Int_t              fNSuperModulesUsed;        ///<  Number of supermodules to be used in analysis, can be different than the real geo, to be used at initialization of histograms.
  
  Int_t              fRunNumber;                ///<  Run number of the data, take it from data itself unless set by user.

  Bool_t             fMCECellClusFracCorrOn;    ///<  Correct or not the weight of cells in cluster.
  
  Float_t            fMCECellClusFracCorrParam[4]; ///<  Parameters for the function correcting the weight of the cells in the cluster.
  
  /// Copy constructor not implemented.
  AliCalorimeterUtils(              const AliCalorimeterUtils & cu) ;
  
  /// Assignment operator not implemented.
  AliCalorimeterUtils & operator = (const AliCalorimeterUtils & cu) ; 
  
  /// \cond CLASSIMP
  ClassDef(AliCalorimeterUtils,19) ;
  /// \endcond

} ;


#endif //ALICALORIMETERUTILS_H



