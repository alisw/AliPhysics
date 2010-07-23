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
class AliAODCaloCluster;
class AliAODCaloCells;
//class AliEMCALCalibData ;
#include "AliEMCALGeoParams.h"

class AliAnalysisTaskEMCALPi0CalibSelection : public AliAnalysisTaskSE
{
public:

  AliAnalysisTaskEMCALPi0CalibSelection(const char* name);
  AliAnalysisTaskEMCALPi0CalibSelection(const AliAnalysisTaskEMCALPi0CalibSelection&); 
  AliAnalysisTaskEMCALPi0CalibSelection& operator=(const AliAnalysisTaskEMCALPi0CalibSelection&); 
  virtual ~AliAnalysisTaskEMCALPi0CalibSelection();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t * opt);
  virtual void LocalInit() ;

  void SetClusterMinEnergy(Float_t emin) {fEmin=emin;}
  void SetClusterMaxEnergy(Float_t emax) {fEmax=emax;}
  void SetClusterMinNCells(Int_t n)      {fMinNCells=n;}
  void SetNCellsGroup(Int_t n)           {fGroupNCells=n;}

  void SetLogWeight(Float_t weight) {fLogWeight=weight;}
  //void SetCalibCorrections(AliEMCALCalibData* const cdata);
  void CreateAODFromESD();
  void CreateAODFromAOD();	

  void CopyAOD(Bool_t copy)   { fCopyAOD = copy ; }
  Bool_t IsAODCopied() const { return fCopyAOD ; }
	
  void SetGeometryName(TString name)   { fEMCALGeoName = name ; }
  TString GeometryName() const { return fEMCALGeoName ; }
 
  // Bad channels, copy from PWG4/PartCorrBase/AliCalorimeterUtils
  Bool_t IsBadChannelsRemovalSwitchedOn()  const { return fRemoveBadChannels ; }
  void SwitchOnBadChannelsRemoval ()  {fRemoveBadChannels = kTRUE  ; InitEMCALBadChannelStatusMap();}
  void SwitchOffBadChannelsRemoval()  {fRemoveBadChannels = kFALSE ; }
	
  void InitEMCALBadChannelStatusMap() ;
	
  Int_t GetEMCALChannelStatus(Int_t iSM , Int_t iCol, Int_t iRow) const { 
	if(fEMCALBadChannelMap) return (Int_t) ((TH2I*)fEMCALBadChannelMap->At(iSM))->GetBinContent(iCol,iRow); 
	else return 0;}//Channel is ok by default
	
  void SetEMCALChannelStatus(Int_t iSM , Int_t iCol, Int_t iRow, Double_t c = 1) { 
	if(!fEMCALBadChannelMap)InitEMCALBadChannelStatusMap() ;
	((TH2I*)fEMCALBadChannelMap->At(iSM))->SetBinContent(iCol,iRow,c);}
	
  TH2I * GetEMCALChannelStatusMap(Int_t iSM) const {return (TH2I*)fEMCALBadChannelMap->At(iSM);}
	
  void SetEMCALChannelStatusMap(TObjArray *map) {fEMCALBadChannelMap = map;}
	
  Bool_t ClusterContainsBadChannel(UShort_t* cellList, Int_t nCells);
	
  // Recalibration
  Bool_t IsRecalibrationOn()  const { return fRecalibration ; }
  void SwitchOnRecalibration()    {fRecalibration = kTRUE ; InitEMCALRecalibrationFactors();}
  void SwitchOffRecalibration()   {fRecalibration = kFALSE ; }
	
  void InitEMCALRecalibrationFactors() ;
	
  Float_t GetEMCALChannelRecalibrationFactor(Int_t iSM , Int_t iCol, Int_t iRow) const { 
	if(fEMCALRecalibrationFactors) return (Float_t) ((TH2F*)fEMCALRecalibrationFactors->At(iSM))->GetBinContent(iCol,iRow); 
	else return 1;}
	
  void SetEMCALChannelRecalibrationFactor(Int_t iSM , Int_t iCol, Int_t iRow, Double_t c = 1) { 
	if(!fEMCALRecalibrationFactors) InitEMCALRecalibrationFactors();
	((TH2F*)fEMCALRecalibrationFactors->At(iSM))->SetBinContent(iCol,iRow,c);}
	
  void SetEMCALChannelRecalibrationFactors(Int_t iSM , TH2F* h) {fEMCALRecalibrationFactors->AddAt(h,iSM);}
	
  TH2F * GetEMCALChannelRecalibrationFactors(Int_t iSM) const {return (TH2F*)fEMCALRecalibrationFactors->At(iSM);}
	
  void SetEMCALChannelRecalibrationFactors(TObjArray *map) {fEMCALRecalibrationFactors = map;}
  Float_t RecalibrateClusterEnergy(AliAODCaloCluster* cluster, AliAODCaloCells * cells);
	
  void SetInvariantMassHistoBinRange(Int_t nBins, Float_t minbin, Float_t maxbin){
	fNbins = nBins; fMinBin = minbin; fMaxBin = maxbin; }
	
private:

  void MaxEnergyCellPos(AliAODCaloCells* const cells, AliAODCaloCluster* const clu, Int_t& iSM, Int_t& ieta, Int_t& iphi);

private:

  AliEMCALGeometry * fEMCALGeo;  //! EMCAL geometry
  //AliEMCALCalibData* fCalibData; // corrections to CC from the previous iteration
	
  Float_t fEmin;          // min. cluster energy
  Float_t fEmax;          // max. cluster energy
  Int_t   fMinNCells;     // min. ncells in cluster
  Int_t fGroupNCells;     // group n cells
  Float_t fLogWeight;     // log weight used in cluster recalibration
  Bool_t  fCopyAOD;       // Copy calo information only to AOD?
  TString fEMCALGeoName;  // Name of geometry to use.
	
  Bool_t     fRemoveBadChannels;         // Check the channel status provided and remove clusters with bad channels
  TObjArray *fEMCALBadChannelMap;        // Array of histograms with map of bad channels, EMCAL
  Bool_t     fRecalibration;             // Switch on or off the recalibration
  TObjArray *fEMCALRecalibrationFactors; // Array of histograms with map of recalibration factors, EMCAL                 
 
  //Output histograms	
  Int_t   fNbins;  // N       mass bins of invariant mass histograms
  Float_t fMinBin; // Minimum mass bins of invariant mass histograms
  Float_t fMaxBin; // Maximum mass bins of invariant mass histograms

  TList*  fOutputContainer; //!histogram container
  TH1F*   fHmpi0[AliEMCALGeoParams::fgkEMCALModules][AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows];//! two-cluster inv. mass assigned to each cell.
  TH1F*   fHmgg;            //! two-cluster inv.mass
  TH1I*   fhNEvents;        //! Number of events counter histogram
  TList * fCuts ;           //! List with analysis cuts

  ClassDef(AliAnalysisTaskEMCALPi0CalibSelection,3);

};

#endif //ALIANALYSISTASKEMCALPI0CALIBSELECTION_H
