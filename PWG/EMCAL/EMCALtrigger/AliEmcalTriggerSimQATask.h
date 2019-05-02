/// \class AliEmcalTriggerSimQATask
///
/// Class to do QA of the EMCal Trigger Simulation
///
/// \author Michael Oliver <michael.oliver@cern.ch>, Yale University
/// \date Dec 7, 2018

#ifndef ALIEMCALTRIGGERSIMQATASK_H
#define ALIEMCALTRIGGERSIMQATASK_H

/**************************************************************************
 * Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

class TClonesArray;
class THistManager;
class AliEMCALGeometry;
//class TH1;
//class TH2;
//class TH3;

class TString;

#include <TObjArray.h>
#include <AliLog.h>

#include "THistManager.h"
#include "AliAnalysisTaskEmcal.h"
#include "AliEMCALGeometry.h"

/**
 * \class AliEmcalTriggerSimQATask
 * \brief EMCAL trigger simulation QA Task
 * Produces QA histograms for the simulated EMCAL trigger in MC.
 */
class AliEmcalTriggerSimQATask : public AliAnalysisTaskEmcal {
 public:

  enum EventEMCALTriggerType_t { // Which bit is set for the event type
    kMB = -1, // Minimum Bias
    kEL0 = 0,
    kDL0 = 1,
    kEG1 = 2,
    kDG1 = 3,
    kEG2 = 4,
    kDG2 = 5,
    kEJ1 = 6,
    kDJ1 = 7,
    kEJ2 = 8,
    kDJ2 = 9
  };

  AliEmcalTriggerSimQATask();
  AliEmcalTriggerSimQATask(const char *name);
  virtual ~AliEmcalTriggerSimQATask();

  static AliEmcalTriggerSimQATask * AddTaskEmcalTriggerSimQA();

  void SetMinAmplitude(Int_t m)             { fMinAmplitude = m; }
  void SetMinClusterEnergy(Float_t m)       { fMinClusterEnergy = m; }


 protected:
  void                                      UserCreateOutputObjects();
  void                                      ExecOnce();
  Bool_t                                    Run();
  Bool_t                                    FillHistograms();
  void                                      FillEventQA();

//  static const Int_t                        kNTriggerTypes = 6; // MB,L0,EG1,EG2,EJ1,EJ2
  static const Int_t                        kNTriggerTypes = 11; // MB,EL0,DL0,EG1,DG1,EG2,DG2,EJ1,DJ1,EJ2,DJ2
//  const TString                             fTriggerNames[kNTriggerTypes] = {"MB","L0","EG1","EG2","EJ1","EJ2"};
  const TString                             fTriggerNames[kNTriggerTypes] = {"MB","EL0","DL0","EG1","DG1","EG2","DG2","EJ1","DJ1","EJ2","DJ2"};
  EventEMCALTriggerType_t                   fTriggerTypes[kNTriggerTypes]; ///< Trigger type array
  TString                                   fTriggerPatchesName;         ///< name of input trigger array
  TClonesArray                             *fTriggerPatches;             //!<! trigger array in

  Int_t                                     fMinAmplitude;               ///< Minimum trigger patch amplitude
  Float_t                                   fMinClusterEnergy;           ///< Minimum cluster energy
  Float_t                                   fPtBinWidth;                 ///< Histogram pt bin width
  Float_t                                   fMaxPt;                      ///< Histogram pt limit

  Int_t                                     fEventTriggerBits;           //!<! Variable storing trigger bits for entire event, set by DoPatchLoop()

  // Histograms
  THistManager                              fHistManager;                ///< Histogram Manager

  void                                      DoPatchLoop();               // Loop over patches, determine trigger condition
  void                                      FillPatchHistograms(AliEMCALTriggerPatchInfo * patch,Int_t i);               // Fill histograms with patch information
  Int_t                                     GeneratePatchTriggerBits(AliEMCALTriggerPatchInfo * patch);               // Generate a trigger bits int for a single patch
  void                                      DoClusterLoop();             // Loop over clusters, fill histograms
  void                                      MatchClusterToPatches(AliVCluster * cluster); // Match cluster to trigger patches, fill histos

 private:
  AliEmcalTriggerSimQATask(const AliEmcalTriggerSimQATask&);            // not implemented
  AliEmcalTriggerSimQATask &operator=(const AliEmcalTriggerSimQATask&); // not implemented


  /// \cond CLASSIMP
  ClassDef(AliEmcalTriggerSimQATask,2);
  /// \endcond
};

#endif



