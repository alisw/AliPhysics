#ifndef ALIANALYSISTASKPI0CALIBSELECTION_H
#define ALIANALYSISTASKPI0CALIBSELECTION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//---------------------------------------------------------------------------// 
// Fill histograms with two-cluster invariant mass                           //
// using calibration coefficients of the previous iteration.                 //
//---------------------------------------------------------------------------//


#include "AliAnalysisTaskSE.h"

class TH1F;
class AliPHOSRecoParam;
class AliPHOSGeometry; 
class AliESDCaloCluster;
class AliESDCaloCells;

class AliAnalysisTaskPi0CalibSelection : public AliAnalysisTaskSE
{
public:

  AliAnalysisTaskPi0CalibSelection();
  AliAnalysisTaskPi0CalibSelection(const char* name);
  virtual ~AliAnalysisTaskPi0CalibSelection();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t * opt);
  
  void SetClusterMinEnergy(Float_t emin) {fEmin=emin;}
  
private:

  AliAnalysisTaskPi0CalibSelection(const AliAnalysisTaskPi0CalibSelection&); 
  AliAnalysisTaskPi0CalibSelection& operator=(const AliAnalysisTaskPi0CalibSelection&); 

  void MaxEnergyCellPos(AliESDCaloCells *cells, AliESDCaloCluster* clu, Int_t& maxId);

private:

  TList* fOutputContainer;
  TH1F*  fHmpi0[5][64][56];// two-cluster inv. mass assigned to each cell.

  AliPHOSRecoParam* fRecoParam; // RecoParameters.
  AliPHOSGeometry * fPhosGeo;   // PHOS geometry

  TH1F* fHmgg; //two-cluster inv.mass
  Float_t fEmin; // min. cluster energy

  ClassDef(AliAnalysisTaskPi0CalibSelection,1);

};

#endif
