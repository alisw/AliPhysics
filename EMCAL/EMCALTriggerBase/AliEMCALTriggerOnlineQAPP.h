/// \class AliEMCALTriggerOnlineQAPP
/// \brief Class to generate EMCal trigger QA plots in pp collisions
///
/// This class generates QA plots for EMCal trigger in pp collisions
///
/// \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
/// \date Apr. 4, 2016

#ifndef ALIEMCALTRIGGERONLINEQAPP_H
#define ALIEMCALTRIGGERONLINEQAPP_H

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
#include <iostream>

// Root
#include <TList.h>

// AliRoot
#include "AliEMCALTriggerQA.h"

class AliEMCALTriggerPatchInfo;
class TObjArray;
class THashList;
class AliEMCALTriggerFastOR;
class AliEMCALGeometry;
class AliVCaloCells;
class TH1;
class TH2;

class AliEMCALTriggerOnlineQAPP : public AliEMCALTriggerQA {
public:

  AliEMCALTriggerOnlineQAPP();
  AliEMCALTriggerOnlineQAPP(const char* name);
  AliEMCALTriggerOnlineQAPP(const AliEMCALTriggerOnlineQAPP& triggerQA);
  virtual ~AliEMCALTriggerOnlineQAPP();

  void   SetL0TimeRange(Int_t min, Int_t max) { fL0MinTime = min; fL0MaxTime = max; }
  void   SetFastORandCellThresholds(Int_t l0, Int_t l1, Double_t cell) { fMinL0FastORAmp = l0; fMinL1FastORAmp = l1; fMinCellAmp = cell; }

  // Overloaded methods of AliEMCALTriggerQA
  void   Init();
  void   ProcessPatch(const AliEMCALTriggerPatchInfo* patch);
  void   ProcessFastor(const AliEMCALTriggerFastOR* fastor, AliVCaloCells* cells = 0);
  void   EventCompleted();
  TCollection* GetListOfHistograms()  { return &fHistograms; }

protected:
  static const Int_t fgkSM = 20;
  static const Int_t fgkNPatchTypes = 3;
  static const Int_t fgkNTriggerTypes = 6;
  static const Int_t fgkNDet = 2;

  Int_t                      fL0MinTime;                   ///< Minimum L0 time
  Int_t                      fL0MaxTime;                   ///< Maximum L0 time
  Double_t                   fMinCellAmp;                  ///< Minimum offline amplitude of the cells
  Int_t                      fMinL0FastORAmp;              ///< Minimum L0 amplitude of the FastORs
  Int_t                      fMinL1FastORAmp;              ///< Minimum L1 amplitude of the FastORs

  AliEMCALTriggerPatchInfo  *fMaxPatchEMCal[fgkNTriggerTypes][fgkNPatchTypes];  //!<! EMCal max patch (will be reset each event)
  AliEMCALTriggerPatchInfo  *fMaxPatchDCal[fgkNTriggerTypes][fgkNPatchTypes];   //!<! DCal max patch (will be reset each event)
  UShort_t                   fNPatches[fgkNDet][fgkNTriggerTypes][fgkNPatchTypes];//!<! Number of patches in the current event
  TList                      fHistograms;                                       //!<! List of histograms

  // General histograms
  TH1                       *fHistEvents;                      //!<! Total number of events

  // TRU histograms
  TH1                       *fHistFastORL0;                    //!<! Counter of FastOR L0 signal above zero
  TH1                       *fHistFastORL0LargeAmp;            //!<! Counter of FastOR L0 signal above some large amplitude
  TH2                       *fHistFastORL0Amp;                 //!<! Amplitude spectra of each FastOR L0 channel
  TH2                       *fHistFastORL0Time;                //!<! Time spectra of each FastOR L0 channel
  TH2                       *fHistFastORL0BySM[fgkSM];         //!<! Counter of FastOR L0 signal above zero (by SM)
  TH2                       *fHistFastORL0LargeAmpBySM[fgkSM]; //!<! Counter of FastOR L0 signal above some large amplitude (by SM)
  TH2                       *fHistFastORL0AmpBySM[fgkSM];      //!<! Integrated amplitude of each FastOR L0 channel (by SM)
  TH2                       *fHistFEEvsTRUBySM[fgkSM];         //!<! Correlation FEE vs TRU (by SM)

  // STU histograms
  TH1                       *fHistFastORL1;                    //!<! Counter of FastOR L1 signal above zero
  TH1                       *fHistFastORL1LargeAmp;            //!<! Counter of FastOR L1 signal above some large amplitude
  TH2                       *fHistFastORL1Amp;                 //!<! Amplitude spectra of each FastOR L1 channel
  TH2                       *fHistFastORL1BySM[fgkSM];         //!<! Counter of FastOR L1 signal above zero (by SM)
  TH2                       *fHistFastORL1LargeAmpBySM[fgkSM]; //!<! Counter of FastOR L1 signal above some large amplitude (by SM)
  TH2                       *fHistFastORL1AmpBySM[fgkSM];      //!<! Integrated amplitude of each FastOR L1 channel (by SM)
  TH2                       *fHistFEEvsSTUBySM[fgkSM];         //!<! Correlation FEE vs STU (by SM)

  // Trigger patch histograms
  TH1                       *fHistNPatches[fgkNDet][fgkNTriggerTypes][fgkNPatchTypes];    //!<! Spectra of patch amplitudes
  TH1                       *fHistPatchAmp[fgkNDet][fgkNTriggerTypes][fgkNPatchTypes];    //!<! Number of patches
  TH1                       *fHistMaxPatchAmp[fgkNDet][fgkNTriggerTypes][fgkNPatchTypes]; //!<! Spectra of maximum patch amplitudes
  TH2                       *fHistMaxEdgePos[fgkNTriggerTypes][fgkNPatchTypes];           //!<! Position of the maximum patch
  TH2                       *fHistAmpEdgePos[fgkNTriggerTypes][fgkNPatchTypes];           //!<! Position of the maximum patch weighted by the amplitude

private:
  AliEMCALTriggerOnlineQAPP &operator=(const AliEMCALTriggerOnlineQAPP &);

  /// \cond CLASSIMP
  ClassDef(AliEMCALTriggerOnlineQAPP, 1);
  /// \endcond
};

#endif
