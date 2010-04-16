#ifndef ALIANALYSISTASKMUONTRACKINGEFF_H
#define ALIANALYSISTASKMUONTRACKINGEFF_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup base
/// \class AliAnalysisTaskMuonTrackingEff
/// \brief tracking chamber efficiency from data
//Author: Nicolas LE BRIS - SUBATECH Nantes


#include "AliAnalysisTask.h"


class AliCheckMuonDetEltResponse;
class AliMUONGeometryTransformer;
class AliESDEvent;
class TClonesArray;
class TH2F;

class AliAnalysisTaskMuonTrackingEff : public AliAnalysisTask
{
 public:
  AliAnalysisTaskMuonTrackingEff();
  AliAnalysisTaskMuonTrackingEff(const AliAnalysisTaskMuonTrackingEff& rhs);
  AliAnalysisTaskMuonTrackingEff& operator=(const AliAnalysisTaskMuonTrackingEff&);
  AliAnalysisTaskMuonTrackingEff(const char* name,
				 Bool_t isCosmic = kFALSE);
  virtual ~AliAnalysisTaskMuonTrackingEff();

  // Implementation of interface methods
  virtual void ConnectInputData(Option_t *option = "");
  virtual void CreateOutputObjects();
  virtual void Exec(Option_t *option);
  virtual void Terminate(Option_t *option);

  static const Int_t fgkTotNbrOfDetectionElt;    ///< The total number of detection element in the tracking system.
  static const Int_t fgkTotNbrOfChamber;

  void ComputeErrors();                        ///< Compute the error on the efficiency (see .cxx for the formula)
  
  void SetCosmic(Bool_t isCosmic) {fIsCosmicData = isCosmic;};
  Bool_t IsCosmic() const {return fIsCosmicData;};

 private:
  AliESDEvent * fESD;                             //!<ESD object
  AliMUONGeometryTransformer *fTransformer;       //!<Transformer object

  TClonesArray* fDetEltEffHistList; //!<Detetcion efficiencies histograms list. 
  TClonesArray* fDetEltTDHistList;  //!<List of histograms of the tracks detected in the detection elements. 
  TClonesArray* fDetEltTTHistList;  //!<List of histograms of the tracks which have passed through the detection elements. 
  TClonesArray* fChamberEffHistList;
  TClonesArray* fChamberTDHistList;
  TClonesArray* fChamberTTHistList;

  AliCheckMuonDetEltResponse* fChamberEff;

  Bool_t fIsCosmicData;

  ClassDef(AliAnalysisTaskMuonTrackingEff, 0)
};

#endif
