/// \class AliEMCALTriggerOfflineQAPP
/// \brief Class to generate EMCal trigger QA plots in pp collisions
///
/// This class generates QA plots for EMCal trigger in pp collisions.
/// This light version is suitable for the general QA train.
///
/// \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
/// \date May 25th, 2017

#ifndef ALIEMCALTRIGGEROFFLINELIGHTQAPP_H
#define ALIEMCALTRIGGEROFFLINELIGHTQAPP_H

/**************************************************************************
* Copyright(c) 1998-2017, ALICE Experiment at CERN, All rights reserved. *
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

// C++
#include <set>
#include <iostream>

// Root
#include <TNamed.h>
#include <TArrayF.h>

// AliRoot
#include <AliEMCALTriggerQA.h>

// AliPhysics
#include "THistManager.h"

class AliEMCALTriggerPatchInfo;
class TObjArray;
class THashList;
class AliEMCALTriggerFastOR;
class AliEMCALGeometry;
class AliVCaloCells;

class AliEMCALTriggerOfflineLightQAPP : public AliEMCALTriggerQA {
public:

  AliEMCALTriggerOfflineLightQAPP();
  AliEMCALTriggerOfflineLightQAPP(const char* name);
  AliEMCALTriggerOfflineLightQAPP(const AliEMCALTriggerOfflineLightQAPP& triggerQA);
  virtual ~AliEMCALTriggerOfflineLightQAPP();

  void   SetL0TimeRange(Int_t min, Int_t max) { fL0MinTime = min; fL0MaxTime = max; }
  void   SetFastORThresholds(Int_t l0, Int_t l1) { fMinL0FastORAmp = l0; fMinL1FastORAmp = l1; }

  // TRU bad channel
  void   AddFastORBadChannel(Short_t absId) { fBadChannels.insert(absId); }
  void   ReadFastORBadChannelFromStream(std::istream& stream);
  void   ReadFastORBadChannelFromFile(const char* fname);

  // Overloaded methods of AliEMCALTriggerQA
  void   Init();
  void   ProcessPatch(const AliEMCALTriggerPatchInfo* patch);
  void   ProcessFastor(const AliEMCALTriggerFastOR* fastor, AliVCaloCells* cells = 0);
  void   ProcessCell(const AliEMCALCellInfo& cell) {;}
  void   EventCompleted();
  TCollection* GetListOfHistograms()  { return fHistManager.GetListOfHistograms(); }
  void   EventTimeStamp(UInt_t timeStamp);


protected:
  std::set<Short_t>          fBadChannels;                 ///< Container of bad channels

  Int_t                      fL0MinTime;                   ///< Minimum L0 time
  Int_t                      fL0MaxTime;                   ///< Maximum L0 time
  Int_t                      fMinL0FastORAmp;              ///< Minimum L0 amplitude of the FastORs
  Int_t                      fMinL1FastORAmp;              ///< Minimum L1 amplitude of the FastORs
  THistManager               fHistManager;                 ///< Histogram manager

  const AliEMCALTriggerPatchInfo  *fMaxPatchEMCal[fgkNTriggerTypes][fgkNPatchTypes];         //!<! EMCal max patch (will be reset each event)
  const AliEMCALTriggerPatchInfo  *fMaxPatchDCal[fgkNTriggerTypes][fgkNPatchTypes];          //!<! DCal max patch (will be reset each event)
  Int_t                      fNTotTRU;                     //!<! Total number of TRUs
  Int_t                      fMaxFORabsId;                 //!<! Maximum FastOR abs id

private:
  AliEMCALTriggerOfflineLightQAPP &operator=(const AliEMCALTriggerOfflineLightQAPP &);

  /// \cond CLASSIMP
  ClassDef(AliEMCALTriggerOfflineLightQAPP, 1);
  /// \endcond
};

#endif
