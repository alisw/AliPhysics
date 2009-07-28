#ifndef ALIANALYSISTASKPHOSPI0CALIBSELECTION_H
#define ALIANALYSISTASKPHOSPI0CALIBSELECTION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//---------------------------------------------------------------------------// 
// Fill histograms with two-cluster invariant mass                           //
// using calibration coefficients of the previous iteration.                 //
//---------------------------------------------------------------------------//

// Root includes
#include "TH1.h"

// AliRoot includes
#include "AliAnalysisTaskSE.h"
#include "AliPHOSGeometry.h"
#include "AliAODCaloCluster.h"
#include "AliAODCaloCells.h"
#include "AliPHOSCalibData.h"

class AliAnalysisTaskPHOSPi0CalibSelection : public AliAnalysisTaskSE
{
public:

  AliAnalysisTaskPHOSPi0CalibSelection();
  AliAnalysisTaskPHOSPi0CalibSelection(const char* name);
  AliAnalysisTaskPHOSPi0CalibSelection(const AliAnalysisTaskPHOSPi0CalibSelection&); 
  AliAnalysisTaskPHOSPi0CalibSelection& operator=(const AliAnalysisTaskPHOSPi0CalibSelection&); 
  virtual ~AliAnalysisTaskPHOSPi0CalibSelection();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t * opt);
  
  void SetClusterMinEnergy(Float_t emin) {fEmin=emin;}
  void SetLogWeight(Float_t weight) {fLogWeight=weight;}
  void SetCalibCorrections(AliPHOSCalibData* cdata);
  void CreateAODFromESD();
  void CreateAODFromAOD();	

  void CopyAOD(Bool_t copy)   { fCopyAOD = copy ; }
  Bool_t IsAODCopied() const { return fCopyAOD ; }

private:

  void MaxEnergyCellPos(AliAODCaloCells *cells, AliAODCaloCluster* clu, Int_t& maxId);

private:

  TList* fOutputContainer;
  TH1F*  fHmpi0[5][64][56];// two-cluster inv. mass assigned to each cell.

  AliPHOSGeometry * fPhosGeo;   // PHOS geometry
  AliPHOSCalibData* fCalibData; // corrections to CC from the previous iteration
	
  TH1F*   fHmgg; //two-cluster inv.mass
  Float_t fEmin; // min. cluster energy
  Float_t fLogWeight; // log weight used in cluster recalibration
  Bool_t  fCopyAOD; // Copy calo information only to AOD?
	
  ClassDef(AliAnalysisTaskPHOSPi0CalibSelection,3);

};

#endif //ALIANALYSISTASKPHOSPI0CALIBSELECTION_H
