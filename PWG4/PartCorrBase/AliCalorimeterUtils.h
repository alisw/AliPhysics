#ifndef ALICALORIMETERUTILS_H
#define ALICALORIMETERUTILS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  $ */

//_________________________________________________________________________
// Class utility for Calorimeter specific selection methods                ///
//
//
//
//-- Author: Gustavo Conesa (LPSC-Grenoble) 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TObject.h" 
#include "TString.h"
#include "TObjArray.h"
class TArrayF;  
#include "TH2I.h"

//--- ANALYSIS system ---
class AliVEvent;
class AliAODPWG4Particle;
class AliVCluster;
class AliVCaloCells;
#include "AliPHOSGeoUtils.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"

class AliCalorimeterUtils : public TObject {

 public:   
  AliCalorimeterUtils() ; // ctor
  virtual ~AliCalorimeterUtils() ;//virtual dtor
 private:
  AliCalorimeterUtils(const AliCalorimeterUtils & g) ; // cpy ctor
  AliCalorimeterUtils & operator = (const AliCalorimeterUtils & g) ;//cpy assignment

 public:
  
  virtual void  InitParameters();
  virtual void  Print(const Option_t * opt)          const ;

  virtual Int_t GetDebug()                           const { return fDebug                ; }
  virtual void  SetDebug(Int_t d)                          { fDebug = d                   ; }
	
  //virtual void Init();
	
  Int_t         GetMaxEnergyCell(AliVCaloCells* cells, AliVCluster* clu, Float_t & fraction) const ;
  
  //Calorimeters Geometry Methods
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

  void          SetGeometryTransformationMatrices(AliVEvent* inputEvent) ;
	
  void          SwitchOnLoadOwnEMCALGeometryMatrices()     { fLoadEMCALMatrices = kTRUE   ; }
  void          SwitchOffLoadOwnEMCALGeometryMatrices()    { fLoadEMCALMatrices = kFALSE  ; }
  void          SetEMCALGeometryMatrixInSM(TGeoHMatrix* m, Int_t i) { fEMCALMatrix[i] = m ; }
  
  void          SwitchOnLoadOwnPHOSGeometryMatrices()      { fLoadPHOSMatrices = kTRUE    ; }
  void          SwitchOffLoadOwnPHOSGeometryMatrices()     { fLoadPHOSMatrices = kFALSE   ; }
  void          SetPHOSGeometryMatrixInSM(TGeoHMatrix* m, Int_t i) { fPHOSMatrix[i] = m   ; }
  
  // Bad channels
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
                  if(fPHOSBadChannelMap)return (Int_t) ((TH2I*)fPHOSBadChannelMap->At(imod))->GetBinContent(iCol,iRow); 
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
	
  Bool_t        ClusterContainsBadChannel(TString calorimeter,UShort_t* cellList, Int_t nCells);
	
  // Mask clusters in front of frame, EMCAL only
  Int_t         GetNMaskCellColumns()                const { return fNMaskCellColumns;}
  void          SetNMaskCellColumns(Int_t n) {
                  if(n > fNMaskCellColumns) { delete [] fMaskCellColumns ; fMaskCellColumns = new Int_t[n] ; }
                  fNMaskCellColumns = n                                                                               ; }
  void          SetMaskCellColumn(Int_t ipos, Int_t icol) { 
                  if(ipos < fNMaskCellColumns) fMaskCellColumns[ipos] = icol;
                  else printf("Not set, position larger than allocated set size first")                               ; }
  Bool_t        MaskFrameCluster(const Int_t iSM, const Int_t ieta) const ;
  
  
  //Calorimeter indexes information
  Int_t         GetModuleNumber(AliAODPWG4Particle * particle, AliVEvent* inputEvent) const;
  Int_t         GetModuleNumber(AliVCluster * cluster) const;
  Int_t         GetModuleNumberCellIndexes(const Int_t absId, const TString calo, Int_t & icol, Int_t & irow, Int_t &iRCU) const ;
	
  //Modules fiducial region
  Bool_t        CheckCellFiducialRegion(AliVCluster* cluster, AliVCaloCells* cells, AliVEvent * event, Int_t iev=0) const ;
  void          SetNumberOfCellsFromPHOSBorder(Int_t n)    { fNCellsFromPHOSBorder = n                                ; }
  Int_t         GetNumberOfCellsFromPHOSBorder()     const { return fNCellsFromPHOSBorder                             ; }
  void          SetNumberOfCellsFromEMCALBorder(Int_t n)   { fEMCALRecoUtils->SetNumberOfCellsFromEMCALBorder(n)      ; }
  Int_t         GetNumberOfCellsFromEMCALBorder()    const { return fEMCALRecoUtils->GetNumberOfCellsFromEMCALBorder(); }
  void          SwitchOnNoFiducialBorderInEMCALEta0()      { fEMCALRecoUtils->SwitchOnNoFiducialBorderInEMCALEta0()   ; }
  void          SwitchOffNoFiducialBorderInEMCALEta0()     { fEMCALRecoUtils->SwitchOffNoFiducialBorderInEMCALEta0()  ; }
  Bool_t        IsEMCALNoBorderAtEta0()              const { return fEMCALRecoUtils->IsEMCALNoBorderAtEta0()          ; }
  
  // Recalibration
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

  Float_t       RecalibrateClusterEnergy(AliVCluster* cluster, AliVCaloCells * cells);

  //EMCAL specific utils for the moment
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

  //Track matching recalculation
  void          RecalculateClusterTrackMatching(AliVEvent * event)         {
                  if (fRecalculateMatching) fEMCALRecoUtils->FindMatches(event,0,fEMCALGeo)     ; }
  void          GetMatchedResiduals(Int_t index, Float_t &dR, Float_t &dZ) {
                  if (fRecalculateMatching) fEMCALRecoUtils->GetMatchedResiduals(index,dR,dZ)   ; }
  
  //This could be used for PHOS ...
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
  
 private:

  Int_t              fDebug;                 //  Debugging level
  TString            fEMCALGeoName;          //  Name of geometry to use for EMCAL.
  TString            fPHOSGeoName;           //  Name of geometry to use for PHOS.	
  AliEMCALGeometry * fEMCALGeo ;             //! EMCAL geometry pointer
  AliPHOSGeoUtils  * fPHOSGeo  ;             //! PHOS  geometry pointer  
  Bool_t             fEMCALGeoMatrixSet;     //  Check if the transformation matrix is set for EMCAL
  Bool_t             fPHOSGeoMatrixSet ;     //  Check if the transformation matrix is set for PHOS
  Bool_t             fLoadEMCALMatrices;     //  Matrices set from configuration, not get from geometry.root or from ESDs/AODs
  TGeoHMatrix *      fEMCALMatrix[10];       //  Geometry matrices with alignments
  Bool_t             fLoadPHOSMatrices;      //  Matrices set from configuration, not get from geometry.root or from ESDs/AODs
  TGeoHMatrix *      fPHOSMatrix[5];         //  Geometry matrices with alignments
  Bool_t             fRemoveBadChannels;     //  Check the channel status provided and remove clusters with bad channels
  TObjArray        * fPHOSBadChannelMap;     //  Array of histograms with map of bad channels, PHOS
  Int_t              fNCellsFromPHOSBorder;  //  Number of cells from PHOS  border the cell with maximum amplitude has to be.
  Int_t              fNMaskCellColumns;      //  Number of masked columns
  Int_t*             fMaskCellColumns;       //[fNMaskCellColumns] list of masked cell collumn
  Bool_t             fRecalibration;         //  Switch on or off the recalibration
  TObjArray        * fPHOSRecalibrationFactors;  // Array of histograms with map of recalibration factors, PHOS
  AliEMCALRecoUtils* fEMCALRecoUtils;        //  EMCAL utils for cluster rereconstruction
  Bool_t             fRecalculatePosition;   //  Recalculate cluster position
  Bool_t             fCorrectELinearity  ;   //  Correct cluster energy linearity
  Bool_t             fRecalculateMatching;   //  Recalculate cluster position
  Float_t            fCutR;                  //  dR cut on matching (PHOS)
  Float_t            fCutZ;                  //  dZ cut on matching (EMCAL/PHOS)
  Float_t            fCutEta;                //  dEta cut on matching (EMCAL)
  Float_t            fCutPhi;                //  dPhi cut on matching (EMCAL)

  
  ClassDef(AliCalorimeterUtils,7)
} ;


#endif //ALICALORIMETERUTILS_H



