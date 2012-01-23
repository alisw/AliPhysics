#ifndef ALIANALYSISTASKMUONTRACKINGEFF_H
#define ALIANALYSISTASKMUONTRACKINGEFF_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

/// \ingroup base
/// \class AliAnalysisTaskMuonTrackingEff
/// \brief tracking chamber efficiency from ESD data
//Author: Nicolas LE BRIS - SUBATECH Nantes

#include "AliAnalysisTaskSE.h"

class AliMUONGeometryTransformer;
class AliESDMuonTrack;
class AliMUONTrackParam;
class TList;
class TString;
class TObjArray;

class AliAnalysisTaskMuonTrackingEff : public AliAnalysisTaskSE
{
 public:
  AliAnalysisTaskMuonTrackingEff();
  AliAnalysisTaskMuonTrackingEff(TString name);
  virtual ~AliAnalysisTaskMuonTrackingEff();

  /// Set location of the default OCDB storage (if not set use "raw://")
  void SetDefaultStorage(const char* ocdbPath) { fOCDBpath = ocdbPath; }
  
  /// set the flag to use only tracks matched with trigger or not
  void MatchTrigger(Bool_t flag = kTRUE) { fMatchTrig = flag; }
  
  /// set the flag to use only tracks passing the acceptance cuts (Rabs, eta)
  void ApplyAccCut(Bool_t flag = kTRUE) { fApplyAccCut = flag; }
  
  /// set the sigma cut value on p*DCA
  void PDCACut(Double_t cut) { fPDCACut = cut; }
  
  /// set the cut value on normalized chi2
  void Chi2Cut(Double_t cut) { fChi2Cut = cut; }
  
  /// set the cut value on minimum pt
  void PtCut(Double_t cut) { fPtCut = cut; }
  
  /// set the cut value on minimum pt
  void UseMCLabel(Bool_t flag = kTRUE) { fUseMCLabel = flag; }
  
  // Implementation of interface methods
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   NotifyRun();
  virtual void   Terminate(Option_t *);
  
  
 private:
  
  /// Not implemented
  AliAnalysisTaskMuonTrackingEff(const AliAnalysisTaskMuonTrackingEff& rhs);
  /// Not implemented
  AliAnalysisTaskMuonTrackingEff& operator = (const AliAnalysisTaskMuonTrackingEff& rhs);
  
  void TrackParamLoop(const TObjArray* trackParams);
  
  void FindAndFillMissedDetElt (const AliMUONTrackParam* trackParam,
				const Bool_t* trackFilter,
				Int_t firstMissCh, Int_t lastChamber);
  
  void CoordinatesOfMissingCluster(Double_t x1, Double_t y1, Double_t z1,
				   Double_t x2, Double_t y2, Double_t z2,
				   Double_t& x, Double_t& y) const;
  
  Bool_t CoordinatesInDetElt(Int_t DeId, Double_t x, Double_t y) const;
  
  void FillTDHistos (Int_t chamber, Int_t detElt, Double_t posXL, Double_t posYL);
  
  void FillTTHistos (Int_t chamber, Int_t detElt, Double_t posXL, Double_t posYL);

  void FillSDHistos (Int_t chamber, Int_t detElt, Double_t posXL, Double_t posYL);
  
  Int_t FromDetElt2iDet (Int_t chamber, Int_t detElt) const;
  Int_t FromDetElt2LocalId (Int_t chamber, Int_t detElt) const;
  Int_t FromLocalId2DetElt(Int_t chamber, Int_t iDet) const;
  
  
private:
  
  static const Int_t fgkNbrOfDetectionElt[10]; ///< The total number of detection element in each chamber.
  static const Int_t fgkOffset;                ///< fFirstDetectionElt[iChamber] = fOffset * (iChamber+1).
  
  Bool_t   fOCDBLoaded;           //!< Determine if the OCDB and =geometry have been loaded
  TString  fOCDBpath;             ///< OCDB path
  Bool_t   fMatchTrig;            ///< use only tracks matched with trigger
  Bool_t   fApplyAccCut;          ///< use only tracks passing the acceptance cuts (Rabs, eta)
  Double_t fPDCACut;              ///< sigma cut on p*DCA
  Double_t fChi2Cut;              ///< cut on normalized chi2
  Double_t fPtCut;                ///< cut on minimum pt
  Bool_t   fUseMCLabel;           ///< select tracks using MC label
  Float_t  fCurrentCentrality;    //!< centrality of the current event
  AliESDMuonTrack* fCurrentTrack; //!< pointer to the currently analyzed track

  AliMUONGeometryTransformer *fTransformer; //!< Transformer object

  TList* fDetEltTDHistList;  //!< List of histograms of the tracks detected in the detection elements. 
  TList* fDetEltTTHistList;  //!< List of histograms of the tracks which have passed through the detection elements. 
  TList* fDetEltSDHistList;  //!< List of histograms of the tracks only detected by one chamber of the station at the detection element level.
  TList* fChamberTDHistList; //!< List of histograms of the tracks detected in the chambers.
  TList* fChamberTTHistList; //!< List of histograms of the tracks which have passed through the chambers.
  TList* fChamberSDHistList; //!< List of histograms of the tracks only detected by one chamber of the station.
  TList* fExtraHistList;     //!< List of extra histograms.

  
  ClassDef(AliAnalysisTaskMuonTrackingEff, 4)
};

#endif

