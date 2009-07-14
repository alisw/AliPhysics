#ifndef ALIANALYSISTASKPHOSPI0CALIBSELECTION_H
#define ALIANALYSISTASKPHOSPI0CALIBSELECTION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//---------------------------------------------------------------------------// 
// Fill histograms with two-cluster invariant mass                           //
// using calibration coefficients of the previous iteration.                 //
//---------------------------------------------------------------------------//


#include "AliAnalysisTaskSE.h"
#include "AliPHOSRecoParam.h"
#include "AliPHOSGeometry.h"
#include "AliAODCaloCluster.h"
#include "AliAODCaloCells.h"
#include "TH1.h"

class AliAnalysisTaskPHOSPi0CalibSelection : public AliAnalysisTaskSE
{
public:

  AliAnalysisTaskPHOSPi0CalibSelection();
  AliAnalysisTaskPHOSPi0CalibSelection(const char* name);
  virtual ~AliAnalysisTaskPHOSPi0CalibSelection();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t * opt);
  
  void SetClusterMinEnergy(Float_t emin) {fEmin=emin;}
  
private:

  AliAnalysisTaskPHOSPi0CalibSelection(const AliAnalysisTaskPHOSPi0CalibSelection&); 
  AliAnalysisTaskPHOSPi0CalibSelection& operator=(const AliAnalysisTaskPHOSPi0CalibSelection&); 

  void MaxEnergyCellPos(AliAODCaloCells *cells, AliAODCaloCluster* clu, Int_t& maxId);

private:

  TList* fOutputContainer;
  TH1F*  fHmpi0[5][64][56];// two-cluster inv. mass assigned to each cell.

  AliPHOSRecoParam* fRecoParam; // RecoParameters.
  AliPHOSGeometry * fPhosGeo;   // PHOS geometry

  TH1F* fHmgg; //two-cluster inv.mass
  Float_t fEmin; // min. cluster energy

  ClassDef(AliAnalysisTaskPHOSPi0CalibSelection,1);

};

#endif //ALIANALYSISTASKPHOSPI0CALIBSELECTION_H
