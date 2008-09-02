#ifndef ALIANALYSISTASKMUONTRACKINGEFF_H
#define ALIANALYSISTASKMUONTRACKINGEFF_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup base
/// \class AliAnalysisTaskMuonTrackingEff
/// \brief tracking chamber efficiency from data
//Author: Nicolas LE BRIS - SUBATECH Nantes


#include "AliAnalysisTask.h"
#include "AliMUONGeometryTransformer.h"
class AliESDEvent;
class TClonesArray;
class TH2F;

class AliAnalysisTaskMuonTrackingEff : public AliAnalysisTask
{
 public:
  AliAnalysisTaskMuonTrackingEff();
  AliAnalysisTaskMuonTrackingEff(const char* name,
		     const AliMUONGeometryTransformer* transformer);
  virtual ~AliAnalysisTaskMuonTrackingEff();

  // Implementation of interface methods
  virtual void ConnectInputData(Option_t *option = "");
  virtual void CreateOutputObjects();
  virtual void Exec(Option_t *option);
  virtual void Terminate(Option_t *option);

  static const Int_t fTotNbrOfDetectionElt;    ///< The total number of detection element in the tracking system.

 private:
  const AliMUONGeometryTransformer* fTransformer;
  AliESDEvent * fESD;               //!<ESD object

  TClonesArray* fDetEltEffHistList; //!<Detetcion efficiencies histograms list. 
  TClonesArray* fDetEltTDHistList;  //!<List of histograms of the tracks detected in the detection elements. 
  TClonesArray* fDetEltTTHistList;  //!<List of histograms of the tracks which have passed through the detection elements. 

  ClassDef(AliAnalysisTaskMuonTrackingEff, 1)
};

#endif
