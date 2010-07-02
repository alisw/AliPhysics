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
class AliAODCaloCluster;
class AliAODCaloCells;
class AliESDCaloCluster;
class AliESDCaloCells;
#include "AliPHOSGeoUtils.h"
#include "AliEMCALGeoUtils.h"

class AliCalorimeterUtils : public TObject {

 public:   
  AliCalorimeterUtils() ; // ctor
  virtual ~AliCalorimeterUtils() ;//virtual dtor
 private:
  AliCalorimeterUtils(const AliCalorimeterUtils & g) ; // cpy ctor
  AliCalorimeterUtils & operator = (const AliCalorimeterUtils & g) ;//cpy assignment

 public:
  
  virtual void InitParameters();
  virtual void Print(const Option_t * opt) const;

  virtual Int_t GetDebug()         const { return fDebug ; }
  virtual void  SetDebug(Int_t d)        { fDebug = d ; }
	
  //virtual void Init();
	
  //Calorimeters Geometry Methods
  AliEMCALGeoUtils * GetEMCALGeometry() const { return fEMCALGeo;}
  TString  EMCALGeometryName()          const { return fEMCALGeoName ; }  
  void  SetEMCALGeometryName(TString name)    { fEMCALGeoName = name ; }
  void InitEMCALGeometry() ; 
  Bool_t IsEMCALGeoMatrixSet()          const { return fEMCALGeoMatrixSet; }
	
  AliPHOSGeoUtils * GetPHOSGeometry() const { return fPHOSGeo;}	
  TString  PHOSGeometryName()         const { return fPHOSGeoName ; }  
  void  SetPHOSGeometryName(TString name)   { fPHOSGeoName = name ; }
  void InitPHOSGeometry() ; 
  Bool_t IsPHOSGeoMatrixSet()         const {return fPHOSGeoMatrixSet ; }

  void SetGeometryTransformationMatrices(AliVEvent* inputEvent) ;
	
  // Bad channels
  Bool_t IsBadChannelsRemovalSwitchedOn()  const { return fRemoveBadChannels ; }
  void SwitchOnBadChannelsRemoval ()  {fRemoveBadChannels = kTRUE  ; InitEMCALBadChannelStatusMap(); InitPHOSBadChannelStatusMap();}
  void SwitchOffBadChannelsRemoval()  {fRemoveBadChannels = kFALSE ; }
	
  void InitEMCALBadChannelStatusMap() ;
  void InitPHOSBadChannelStatusMap () ;

  Int_t GetEMCALChannelStatus(Int_t iSM , Int_t iCol, Int_t iRow) const { 
    if(fEMCALBadChannelMap) return (Int_t) ((TH2I*)fEMCALBadChannelMap->At(iSM))->GetBinContent(iCol,iRow); 
    else return 0;}//Channel is ok by default

  Int_t GetPHOSChannelStatus (Int_t imod, Int_t iCol, Int_t iRow) const { 
    if(fPHOSBadChannelMap)return (Int_t) ((TH2I*)fPHOSBadChannelMap->At(imod))->GetBinContent(iCol,iRow); 
    else return 0;}//Channel is ok by default
  
  void SetEMCALChannelStatus(Int_t iSM , Int_t iCol, Int_t iRow, Double_t c = 1) { 
    if(!fEMCALBadChannelMap)InitEMCALBadChannelStatusMap() ;
    ((TH2I*)fEMCALBadChannelMap->At(iSM))->SetBinContent(iCol,iRow,c);}
  
  void SetPHOSChannelStatus (Int_t imod, Int_t iCol, Int_t iRow, Double_t c = 1) {
	if(!fPHOSBadChannelMap) InitPHOSBadChannelStatusMap() ; 
	((TH2I*)fPHOSBadChannelMap->At(imod))->SetBinContent(iCol,iRow,c);}
    
  TH2I * GetEMCALChannelStatusMap(Int_t iSM) const {return (TH2I*)fEMCALBadChannelMap->At(iSM);}
  TH2I * GetPHOSChannelStatusMap(Int_t imod) const {return (TH2I*)fPHOSBadChannelMap->At(imod);}

  void SetEMCALChannelStatusMap(TObjArray *map) {fEMCALBadChannelMap = map;}
  void SetPHOSChannelStatusMap (TObjArray *map) {fPHOSBadChannelMap  = map;}
	
  Bool_t ClusterContainsBadChannel(TString calorimeter,UShort_t* cellList, Int_t nCells);
	
  //Calorimeter indexes information
  Int_t GetModuleNumber(AliAODPWG4Particle * particle, AliVEvent* inputEvent) const;
  Int_t GetModuleNumber(AliAODCaloCluster * cluster) const;
  Int_t GetModuleNumber(AliESDCaloCluster * cluster) const;

  Int_t GetModuleNumberCellIndexes(const Int_t absId, const TString calo, Int_t & icol, Int_t & irow, Int_t &iRCU) const ;
	
  //Modules fiducial region
  Bool_t CheckCellFiducialRegion(AliAODCaloCluster* cluster, AliAODCaloCells* cells) const ;
  Bool_t CheckCellFiducialRegion(AliESDCaloCluster* cluster, AliESDCaloCells* cells) const ;
	
  void   SetNumberOfCellsFromEMCALBorder(Int_t n) {fNCellsFromEMCALBorder = n; }
  Int_t  GetNumberOfCellsFromEMCALBorder() const  {return fNCellsFromEMCALBorder; }
  void   SetNumberOfCellsFromPHOSBorder(Int_t n)  {fNCellsFromPHOSBorder = n; }
  Int_t  GetNumberOfCellsFromPHOSBorder() const   {return fNCellsFromPHOSBorder; }
	
  void   SwitchOnNoFiducialBorderInEMCALEta0()  {fNoEMCALBorderAtEta0 = kTRUE; }
  void   SwitchOffNoFiducialBorderInEMCALEta0() {fNoEMCALBorderAtEta0 = kFALSE; }
	
  // Recalibration
  Bool_t IsRecalibrationOn()  const { return fRecalibration ; }
  void SwitchOnRecalibration()    {fRecalibration = kTRUE ; InitEMCALRecalibrationFactors(); InitPHOSRecalibrationFactors();}
  void SwitchOffRecalibration()   {fRecalibration = kFALSE ; }
	
  void InitEMCALRecalibrationFactors() ;
  void InitPHOSRecalibrationFactors () ;
	
  Float_t GetEMCALChannelRecalibrationFactor(Int_t iSM , Int_t iCol, Int_t iRow) const { 
    if(fEMCALRecalibrationFactors) return (Float_t) ((TH2F*)fEMCALRecalibrationFactors->At(iSM))->GetBinContent(iCol,iRow); 
    else return 1;}
  
  Float_t GetPHOSChannelRecalibrationFactor (Int_t imod, Int_t iCol, Int_t iRow) const { 
    if(fPHOSRecalibrationFactors)return (Float_t) ((TH2F*)fPHOSRecalibrationFactors->At(imod))->GetBinContent(iCol,iRow); 
    else return 1;}
  
  void SetEMCALChannelRecalibrationFactor(Int_t iSM , Int_t iCol, Int_t iRow, Double_t c = 1) { 
    if(!fEMCALRecalibrationFactors) InitEMCALRecalibrationFactors();
    ((TH2F*)fEMCALRecalibrationFactors->At(iSM))->SetBinContent(iCol,iRow,c);}
	
  void SetPHOSChannelRecalibrationFactor (Int_t imod, Int_t iCol, Int_t iRow, Double_t c = 1) {
    if(!fPHOSRecalibrationFactors)  InitPHOSRecalibrationFactors();
    ((TH2F*)fPHOSRecalibrationFactors->At(imod))->SetBinContent(iCol,iRow,c);}
    
  void SetEMCALChannelRecalibrationFactors(Int_t iSM , TH2F* h) {fEMCALRecalibrationFactors->AddAt(h,iSM);}
  void SetPHOSChannelRecalibrationFactors(Int_t imod , TH2F* h) {fPHOSRecalibrationFactors ->AddAt(h,imod);}
	
  TH2F * GetEMCALChannelRecalibrationFactors(Int_t iSM) const {return (TH2F*)fEMCALRecalibrationFactors->At(iSM);}
  TH2F * GetPHOSChannelRecalibrationFactors(Int_t imod) const {return (TH2F*)fPHOSRecalibrationFactors->At(imod);}
	
  void SetEMCALChannelRecalibrationFactors(TObjArray *map) {fEMCALRecalibrationFactors = map;}
  void SetPHOSChannelRecalibrationFactors (TObjArray *map) {fPHOSRecalibrationFactors  = map;}

  Float_t RecalibrateClusterEnergy(AliESDCaloCluster* cluster, AliESDCaloCells * cells);
  Float_t RecalibrateClusterEnergy(AliAODCaloCluster* cluster, AliAODCaloCells * cells);

 private:

  Int_t              fDebug;                 //  Debugging level
  TString            fEMCALGeoName;          //  Name of geometry to use for EMCAL.
  TString            fPHOSGeoName;           //  Name of geometry to use for PHOS.	
  AliEMCALGeoUtils * fEMCALGeo ;             //! EMCAL geometry pointer
  AliPHOSGeoUtils  * fPHOSGeo  ;             //! PHOS  geometry pointer  
  Bool_t             fEMCALGeoMatrixSet;     //  Check if the transformation matrix is set for EMCAL
  Bool_t             fPHOSGeoMatrixSet ;     //  Check if the transformation matrix is set for PHOS
  Bool_t             fRemoveBadChannels;     //  Check the channel status provided and remove clusters with bad channels
  TObjArray        * fEMCALBadChannelMap;    //! Array of histograms with map of bad channels, EMCAL
  TObjArray        * fPHOSBadChannelMap;     //! Array of histograms with map of bad channels, PHOS
  Int_t              fNCellsFromEMCALBorder; //  Number of cells from EMCAL border the cell with maximum amplitude has to be.
  Int_t              fNCellsFromPHOSBorder;  //  Number of cells from PHOS  border the cell with maximum amplitude has to be.
  Bool_t             fNoEMCALBorderAtEta0;   //  Do fiducial cut in EMCAL region eta = 0?
  Bool_t             fRecalibration;         //  Switch on or off the recalibration
  TObjArray        * fEMCALRecalibrationFactors; //! Array of histograms with map of recalibration factors, EMCAL
  TObjArray        * fPHOSRecalibrationFactors;  //! Array of histograms with map of recalibration factors, PHOS

  ClassDef(AliCalorimeterUtils,2)
} ;


#endif //ALICALORIMETERUTILS_H



