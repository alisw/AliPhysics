#ifndef ALIANALYSISTASKMUONTRACKINGEFF_H
#define ALIANALYSISTASKMUONTRACKINGEFF_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

/// \ingroup base
/// \class AliAnalysisTaskMuonTrackingEff
/// \brief tracking chamber efficiency from ESD data
//Author: Nicolas LE BRIS - SUBATECH Nantes

#include "AliAnalysisTask.h"

class AliCheckMuonDetEltResponse;
class AliMUONGeometryTransformer;
class AliESDEvent;
class TList;
class TString;

class AliAnalysisTaskMuonTrackingEff : public AliAnalysisTask
{
 public:
  AliAnalysisTaskMuonTrackingEff();
  AliAnalysisTaskMuonTrackingEff(const AliAnalysisTaskMuonTrackingEff& rhs);
  AliAnalysisTaskMuonTrackingEff& operator=(const AliAnalysisTaskMuonTrackingEff&);
  AliAnalysisTaskMuonTrackingEff(TString name, TString path = "alien://folder=/alice/data/2010/OCDB");
  virtual ~AliAnalysisTaskMuonTrackingEff();

  // Implementation of interface methods
  virtual void ConnectInputData(Option_t *option = "");
  virtual void CreateOutputObjects();
  virtual void Exec(Option_t *option);
  virtual void Terminate(Option_t *option);

 private:
  Bool_t fIsInit;                   //!< Determine if the object has been initialized
  Bool_t fIsLoaded;                    //!< Determine if the OCDB and =geometry have been loaded
  TString fOCDBpath;                //!< OCDB path
  Int_t fUsableTracks;              //!< Number of usable tracks

  AliESDEvent * fESD;                             //< ESD object
  AliMUONGeometryTransformer *fTransformer;       //< Transformer object

  static const Int_t fgkTotNbrOfDetectionElt;    //< The total number of detection element in the tracking system.
  static const Int_t fgkTotNbrOfChamber;         //< The total number of tracking chambers in the muon spectrometer

  TList* fDetEltTDHistList;    //< List of histograms of the tracks detected in the detection elements. 
  TList* fDetEltTTHistList;    //< List of histograms of the tracks which have passed through the detection elements. 
  TList* fChamberTDHistList;   //< List of histograms of the tracks detected in the chambers.
  TList* fChamberTTHistList;   //< List of histograms of the tracks which have passed through the chambers.

  AliCheckMuonDetEltResponse* fChamberEff;  //< Class computing the efficiency

  void Init();                                 // Initialize the object  
  void LoadOCDBandGeometry();                  // Load the OCDB and the geometry

  ClassDef(AliAnalysisTaskMuonTrackingEff, 1)
};

#endif
