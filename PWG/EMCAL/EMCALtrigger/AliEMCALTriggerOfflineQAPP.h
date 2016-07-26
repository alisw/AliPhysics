/// \class AliEMCALTriggerOfflineQAPP
/// \brief Class to generate EMCal trigger QA plots in pp collisions
///
/// This class generates QA plots for EMCal trigger in pp collisions
///
/// \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
/// \date Apr. 4, 2016

#ifndef ALIEMCALTRIGGEROFFLINEQAPP_H
#define ALIEMCALTRIGGEROFFLINEQAPP_H

/**************************************************************************
* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

class AliEMCALTriggerOfflineQAPP : public AliEMCALTriggerQA {
public:

  AliEMCALTriggerOfflineQAPP();
  AliEMCALTriggerOfflineQAPP(const char* name);
  AliEMCALTriggerOfflineQAPP(const AliEMCALTriggerOfflineQAPP& triggerQA);
  virtual ~AliEMCALTriggerOfflineQAPP();

  void   SetL0TimeRange(Int_t min, Int_t max) { fL0MinTime = min; fL0MaxTime = max; }
  void   SetFastORandCellThresholds(Int_t l0, Int_t l1, Double_t cell) { fMinL0FastORAmp = l0; fMinL1FastORAmp = l1; fMinCellAmp = cell; }
  void   EnableDCal(Bool_t e = kTRUE) { fDCalPlots = e; }

  // TRU bad channel
  void   AddFastORBadChannel(Short_t absId) { fBadChannels.insert(absId); }
  void   ReadFastORBadChannelFromStream(std::istream& stream);
  void   ReadFastORBadChannelFromFile(const char* fname);
  void   AddOfflineBadChannel(Short_t absId)  { fOfflineBadChannels.insert(absId)   ; }
  void   ReadOfflineBadChannelFromFile(const char* fname);
  void   ReadOfflineBadChannelFromStream(std::istream& stream);
  void   ReadFastORPedestalFromStream(std::istream& stream);
  void   ReadFastORPedestalFromFile(const char* fname);
  void   SetFastORPedestal(Short_t absId, Float_t ped);
  void   ResetFastORPedestal() { fFastORPedestal.Reset(); }

  // Overloaded methods of AliEMCALTriggerQA
  void   Init();
  void   ProcessPatch(const AliEMCALTriggerPatchInfo* patch);
  void   ProcessFastor(const AliEMCALTriggerFastOR* fastor, AliVCaloCells* cells = 0);
  void   ProcessCell(const AliEMCALCellInfo& cell);
  void   EventCompleted();
  TCollection* GetListOfHistograms()  { return fHistManager.GetListOfHistograms(); }
  void   EventTimeStamp(UInt_t timeStamp);


protected:

  std::set<Short_t>          fOfflineBadChannels;          ///< Abs ID of offline bad channels
  std::set<Short_t>          fBadChannels;                 ///< Container of bad channels
  TArrayF                    fFastORPedestal;              ///< FastOR pedestal

  Bool_t                     fDCalPlots;                   ///< Whether to add DCal QA plots
  Int_t                      fL0MinTime;                   ///< Minimum L0 time
  Int_t                      fL0MaxTime;                   ///< Maximum L0 time
  Double_t                   fMinCellAmp;                  ///< Minimum offline amplitude of the cells
  Int_t                      fMinL0FastORAmp;              ///< Minimum L0 amplitude of the FastORs
  Int_t                      fMinL1FastORAmp;              ///< Minimum L1 amplitude of the FastORs
  THistManager               fHistManager;                 ///< Histogram manager

  const AliEMCALTriggerPatchInfo  *fMaxPatchEMCal[6][3];         //!<! EMCal max patch (will be reset each event)
  const AliEMCALTriggerPatchInfo  *fMaxPatchDCal[6][3];          //!<! DCal max patch (will be reset each event)
  Double_t                   fSumOfflineEMCal;             //!<! EMCal sum of all offline energy deposition (will be reset each event)
  Int_t                      fSumL0EMCal;                  //!<! EMCal sum of all online energy deposition (will be reset each event)
  Int_t                      fSumL1EMCal;                  //!<! EMCal sum of all online energy deposition (will be reset each event)
  Double_t                   fSumOfflineDCal;              //!<! DCal sum of all offline energy deposition (will be reset each event)
  Int_t                      fSumL0DCal;                   //!<! DCal sum of all online energy deposition (will be reset each event)
  Int_t                      fSumL1DCal;                   //!<! DCal sum of all online energy deposition (will be reset each event)
  Int_t                      fNCellEMCal;                  //!<! EMCal number of offline cells (will be reset each event)
  Int_t                      fNL0EMCal;                    //!<! EMCal number of L0 FastORs (will be reset each event)
  Int_t                      fNL1EMCal;                    //!<! EMCal number of L1 FastORs (will be reset each event)
  Int_t                      fNCellDCal;                   //!<! DCal number of offline cells (will be reset each event)
  Int_t                      fNL0DCal;                     //!<! DCal number of L0 FastORs (will be reset each event)
  Int_t                      fNL1DCal;                     //!<! DCal number of L1 FastORs (will be reset each event)
  Int_t                      fNTotTRU;                     //!<! Total number of TRUs
  Int_t                      fMaxFORabsId;                 //!<! Maximum FastOR abs id

private:
  AliEMCALTriggerOfflineQAPP &operator=(const AliEMCALTriggerOfflineQAPP &);

  /// \cond CLASSIMP
  ClassDef(AliEMCALTriggerOfflineQAPP, 1);
  /// \endcond
};

#endif
