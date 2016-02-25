/**
 * \file AliReducedHighPtEventCreator.h
 * \brief Declaration of class AliReducedHighPtEventCreator and AliReducedTrackSelectionContainer
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Apr 15, 2015
 */
#ifndef ALIREDUCEDHIGHPTEVENTCREATOR_H
#define ALIREDUCEDHIGHPTEVENTCREATOR_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskEmcal.h"
#include "AliVEvent.h"
#include <TString.h>

class TArrayI;
class TTree;
class TObjArray;

class AliEmcalTrackSelection;
class AliVCluster;
class AliVEvent;
class AliVParticle;

namespace HighPtTracks {

class AliReducedHighPtEvent;
class AliReducedPatchContainer;

/**
 * \class AliReducedTrackSelectionContainer
 * \brief Helper class mapping a virtual track selection to a cut index
 */
class AliReducedTrackSelectionContainer : public TObject{
public:
  AliReducedTrackSelectionContainer();
  AliReducedTrackSelectionContainer(Int_t index, AliEmcalTrackSelection * sel);
  virtual ~AliReducedTrackSelectionContainer();

  /**
   * Set the index of the track selection
   * \param inde Index of the track selection
   */
  void SetIndex(Int_t index) { fIndex = index; }
  /**
   * Set the underlying track selection
   * \param trackSelection Underlying track selection
   */
  void SetTrackSelection(AliEmcalTrackSelection *trackSelection) { fTrackSelection = trackSelection;}

  /**
   * Get the index of the track selection
   * \return Index of the track selection
   */
  Int_t GetIndex() const { return fIndex; }
  /**
   * Get the underlying track selection
   * \return underlying track selection
   */
  AliEmcalTrackSelection *GetTrackSelection() const { return fTrackSelection; }

protected:
  Int_t                                                           fIndex;               ///< Index of the cut
  AliEmcalTrackSelection                                          *fTrackSelection;     ///< Virtual track selection

private:
  AliReducedTrackSelectionContainer(AliReducedTrackSelectionContainer &);
  AliReducedTrackSelectionContainer &operator=(AliReducedTrackSelectionContainer &);

  /// \cond CLASSIMP
  ClassDef(AliReducedTrackSelectionContainer, 1);
  /// \endcond
};

/**
 * \class AliReducedHighPtEventCreator
 * \brief Producer task of reduced high-\f$ p_{t} \f$ events
 */
class AliReducedHighPtEventCreator: public AliAnalysisTaskEmcal {
public:
  AliReducedHighPtEventCreator();
  AliReducedHighPtEventCreator(const char *name);
  virtual ~AliReducedHighPtEventCreator();

  virtual void UserCreateOutputObjects();
  virtual Bool_t Run();

  void AddVirtualTrackSelection(AliEmcalTrackSelection * sel, Int_t index);
  /**
   * Define if trigger thresholds are swapped
   * \param doswap if true then thresholds are swapped
   */
  void SetSwapTriggerThresholds(Bool_t doswap) { fSwapTriggerThresholds = doswap; }
  /**
   * Set the energy range for clusters
   * \param minE Minimum
   * \param maxE Maximum
   */
  void SetClusterEnergyCut(Double_t minE, Double_t maxE) {
    fMinClusterE = minE;
    fMaxClusterE = maxE;
  }
  /**
   * Set the \f$ p_{t} \f$ range for particles and tracks
   * \param minpt Minimum
   * \param maxpt Maximum
   */
  void SetPtRange(Double_t minpt, Double_t maxpt){
    fMinPt = minpt;
    fMaxPt = maxpt;
  }
  /**
   * Set the \f$ \eta \f$ range for particles and tracks
   * \param mineta Minimum
   * \param maxeta Maximum
   */
  void SetEtaRange(Double_t mineta, Double_t maxeta){
    fMinEta = mineta;
    fMaxEta = maxeta;
  }

  /**
   * Set the method handling the centrality
   * \param centmethod Method handling centrality
   */
  void SetCentralityMethod(const char *centmethod){
    fCentralityMethod = centmethod;
  }
  /**
   * Set trigger bit for min. bias event selection
   * \param minbiasbit Bit number of the min. bias bit in AliVEvent
   */
  void SetMinBiasSelection(AliVEvent::EOfflineTriggerTypes minbiasbit) { fMinBiasSelection = minbiasbit; }

  /**
   * Switch for applying selection on a given centrality range in the tree production
   * \param doApply
   */
  void SetApplyCentralitySelection(Bool_t doApply = kTRUE) { fApplyCentralitySelection = doApply; }

  /**
   * Set the range of the centrality selection
   * \param mincent Minimum centrality percentile
   * \param maxcent Maximum centrality percentile
   */
  void SetCentraltityRange(Float_t mincent, Float_t maxcent) {
    fSelectCentralityRange[0] = mincent;
    fSelectCentralityRange[1] = maxcent;
  }
  /**
   * Set the EMCAL trigger selection
   * \param triggersetup
   */
  void SetTriggerSetup(const char *triggersetup){
    fTriggerSetup = triggersetup;
  }
  /**
   * Set event selection filter bits
   * @param selection Event selection filter bits
   */
  void SetBasicEventSelection(UInt_t selection) { fEventSelectionBits = selection; }

  /**
   * Keep only a given fraction of events. Events will be selected randomly
   * after full event cuts.
   * @param frac Fraction of events to keep.
   */
  void SetFractionOfEventsToKeep(Double_t frac){
    if(frac < 0) fKeepFractionEvents = 0.;
    else if(frac > 1) fKeepFractionEvents = 1.;
    else fKeepFractionEvents = frac;
  }

protected:
  Bool_t SelectEvent(AliVEvent *event) const;
  Bool_t SelectCluster(const AliVCluster *clust) const;
  Int_t  SelectTrack(AliVTrack *track, TArrayI &cutindices) const;
  void ConvertTriggerPatches(TClonesArray *patches, AliReducedPatchContainer *cont);
  void FixTrackInputEvent(AliVTrack *trk);
  Int_t GetTPCCrossedRows(const AliVTrack *trk) const;
  void GetCellEnergies(AliVCluster *emccluster, TArrayD &energies) const;
  TTree                   *fOutputTree;                   ///< Output tree
  AliReducedHighPtEvent   *fOutputEvent;                  ///< Output event
  TObjArray               *fTrackSelections;              ///< List of track selections

  UInt_t                    fEventSelectionBits;          ///< Basic event selection
  Bool_t                    fSwapTriggerThresholds;       ///< Switch for swapping of the thresholds
  Double_t                  fMinClusterE;                 ///< Min. cluster E
  Double_t                  fMaxClusterE;                 ///< Max. cluster E
  Double_t                  fMinPt;                       ///< Min. track \f$ p_{t} \f$
  Double_t                  fMaxPt;                       ///< Max. track \f$ p_{t} \f$
  Double_t                  fMinEta;                      ///< Min. track \f$ \eta \f$
  Double_t                  fMaxEta;                      ///< Max. track \f$ \eta \f$
  Float_t                   fKeepFractionEvents;          ///< Keep fraction of events (used for downscaling)
  Bool_t                    fApplyCentralitySelection;    ///< Switch for doing centrality selection
  TString                   fCentralityMethod;            ///< Method to obtain the event centrality
  TString                   fTriggerSetup;                ///< Trigger setup
  Float_t                   fSelectCentralityRange[2];    // Range for centraltity selection
  AliVEvent::EOfflineTriggerTypes fMinBiasSelection;      ///< Trigger bit for min. bias event tagging

private:
  AliReducedHighPtEventCreator(const AliReducedHighPtEventCreator &);
  AliReducedHighPtEventCreator &operator=(const AliReducedHighPtEventCreator &);

  /// \cond CLASSIMP
  ClassDef(AliReducedHighPtEventCreator, 1)
  /// \CLASSIMP
};

} /* namespace HighPtTracks */

#endif /* ALIREDUCEDHIGHPTEVENTCREATOR_H */
