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

// AliRoot includes
#include "AliAnalysisTaskSE.h"
class AliEMCALGeometry;
class AliAODCaloCluster;
class AliAODCaloCells;
class AliEMCALCalibData ;
#include "AliEMCALGeoParams.h"

class AliAnalysisTaskEMCALPi0CalibSelection : public AliAnalysisTaskSE
{
public:

  AliAnalysisTaskEMCALPi0CalibSelection();
  AliAnalysisTaskEMCALPi0CalibSelection(const char* name);
  AliAnalysisTaskEMCALPi0CalibSelection(const AliAnalysisTaskEMCALPi0CalibSelection&); 
  AliAnalysisTaskEMCALPi0CalibSelection& operator=(const AliAnalysisTaskEMCALPi0CalibSelection&); 
  virtual ~AliAnalysisTaskEMCALPi0CalibSelection();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t * opt);
  
  void SetClusterMinEnergy(Float_t emin) {fEmin=emin;}
  void SetLogWeight(Float_t weight) {fLogWeight=weight;}
  void SetCalibCorrections(AliEMCALCalibData* const cdata);
  void CreateAODFromESD();
  void CreateAODFromAOD();	

  void CopyAOD(Bool_t copy)   { fCopyAOD = copy ; }
  Bool_t IsAODCopied() const { return fCopyAOD ; }
	
  void SetGeometryName(TString name)   { fEMCALGeoName = name ; }
  TString GeometryName() const { return fEMCALGeoName ; }

  void SetOldData(Bool_t bData) {fOldData = bData;} 
  Bool_t IsOldData() const {return fOldData;}	
 	
private:

  void MaxEnergyCellPos(AliAODCaloCells* const cells, AliAODCaloCluster* const clu, Int_t& maxId);

private:

  AliEMCALGeometry * fEMCALGeo;   // EMCAL geometry
  AliEMCALCalibData* fCalibData; // corrections to CC from the previous iteration
	
  Float_t fEmin;          // min. cluster energy
  Float_t fLogWeight;     // log weight used in cluster recalibration
  Bool_t  fCopyAOD;       // Copy calo information only to AOD?
  TString fEMCALGeoName;  // Name of geometry to use.
  Bool_t  fOldData ;      // calibrate digit amplitude, since not done in old ESDs
                         
  //Output histograms	
  TList*  fOutputContainer; //histogram container
  TH1F*   fHmpi0[AliEMCALGeoParams::fgkEMCALModules][AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows];// two-cluster inv. mass assigned to each cell.
  TH1F*   fHmgg;          // two-cluster inv.mass

  ClassDef(AliAnalysisTaskEMCALPi0CalibSelection,1);

};

#endif //ALIANALYSISTASKEMCALPI0CALIBSELECTION_H
