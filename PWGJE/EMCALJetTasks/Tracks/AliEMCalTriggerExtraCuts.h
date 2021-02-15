/**
 * \file AliEMCalTriggerExtraCuts.h
 * \brief Declaration of class AliEMCalTriggerExtraCuts
 *
 * In this header file the class AliEMCalTriggerExtraCuts, which provides additional cuts for
 * the track selection of high-\f$ p_{t} \f$ tracks, is declared
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Apr 2, 2015
 */
#ifndef ALIEMCALTRIGGEREXTRACUTS_H
#define ALIEMCALTRIGGEREXTRACUTS_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TBits.h>
#include "AliVCuts.h"

class AliVTrack;

namespace PWGJE {
  
namespace EMCALJetTasks {

/**
 * \class AliEMCalTriggerExtraCuts
 * \brief Extra track selection cuts for the high-\f$ p_{t} \f$ track analysis
 *
 * This class implements cuts necessary for the analysis of high-\f$ p_{t} \f$ tracks
 * in triggered events which are not implemented up to now in a general track selection
 * class. Among these we have
 * - Minimum crossed rows in the TPC
 * - Minimum track length in the TPC
 */
class AliEMCalTriggerExtraCuts: public AliVCuts {
public:
  AliEMCalTriggerExtraCuts();
  virtual ~AliEMCalTriggerExtraCuts() {}

  /**
   * Setter for minimum crossed rows in the TPC
   * \param crossedRows Minimum crossed rows requested
   */
  void SetMinTPCCrossedRows(Int_t crossedRows){
    fMinCrossedRowsTPC = crossedRows;
    fRequestBitmap.SetBitNumber(kTPCCrossedRows);
  }

  /**
   * Setter for minimum track length requested
   * \param tracklength Minimum track length requested
   */
  void SetMinTPCTrackLengthCut(){
    fRequestBitmap.SetBitNumber(kTPCTrackLength);
  }

  virtual Bool_t IsSelected(TObject *o);

protected:
  /**
   * \enum CutType_t
   * \brief Bit definition for different track selection bits
   */
  enum CutType_t{
    kTPCCrossedRows         = 0,      ///< Bit for TPC crossed rows cut
    kTPCTrackLength         = 1       ///< Bit for TPC track length cut
  };

  Float_t GetTPCCrossedRows(const AliVTrack *const trk) const;
  Double_t CalculateTPCTrackLength(AliVTrack *trk) const;
  Int_t                 fMinCrossedRowsTPC;             ///< Min. number of crossed rows in the TPC
  TBits                 fRequestBitmap;                 ///< Bitmap for cuts enabled

  ClassDef(AliEMCalTriggerExtraCuts, 1);
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* PWGJE_EMCALJETTASKS_TRACKS_ALIEMCALTRIGGEREXTRACUTS_H_ */
