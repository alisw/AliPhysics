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
/**
 * @file AliEMCALTriggerQA.cxx
 * @date Apr. 4, 2016
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 */

#include "AliEMCALTriggerPatchInfo.h"
#include "AliEMCALTriggerConstants.h"
#include "AliEMCALGeometry.h"

#include "AliEMCALTriggerQA.h"

using namespace EMCALTrigger;

/// \cond CLASSIMP
ClassImp(AliEMCALTriggerQA)
/// \endcond

const Int_t AliEMCALTriggerQA::fgkMaxPatchAmp[6] = {2000, 2000, 2000, 6000, 6000, 5000};
const TString AliEMCALTriggerQA::fgkPatchTypes[3] = {"Online", "Recalc", "Offline"};

/**
 * Default constructor for ROOT I/O
 */
AliEMCALTriggerQA::AliEMCALTriggerQA():
  TNamed(),
  fFastorL0Th(400),
  fFastorL1Th(400),
  fADCperBin(16),
  fDebugLevel(0),
  fTimeStampBinWidth(0),
  fGeom(0),
  fEventTimeStamp(0),
  fEventTimeStampBin(0)
{
  memset(fEnabledPatchTypes, 0, sizeof(Bool_t)*3);
  memset(fEnabledTriggerTypes, 0, sizeof(Bool_t)*6);
}

/**
 * Named constructor
 * \param name Name of the object
 */
AliEMCALTriggerQA::AliEMCALTriggerQA(const char* name):
  TNamed(name,name),
  fFastorL0Th(400),
  fFastorL1Th(400),
  fADCperBin(16),
  fDebugLevel(0),
  fTimeStampBinWidth(0),
  fGeom(0),
  fEventTimeStamp(0),
  fEventTimeStampBin(0)
{
  memset(fEnabledPatchTypes, 0, sizeof(Bool_t)*3);
  memset(fEnabledTriggerTypes, 0, sizeof(Bool_t)*6);
}

/**
* Copy Constructor
* \param ref Constant reference to copy from
*/
AliEMCALTriggerQA::AliEMCALTriggerQA(const AliEMCALTriggerQA& ref) :
  TNamed(ref),
  fFastorL0Th(ref.fFastorL0Th),
  fFastorL1Th(ref.fFastorL1Th),
  fADCperBin(ref.fADCperBin),
  fDebugLevel(ref.fDebugLevel),
  fTimeStampBinWidth(0),
  fGeom(0),
  fEventTimeStamp(0),
  fEventTimeStampBin(0)
{
  memcpy(fEnabledPatchTypes, ref.fEnabledPatchTypes, sizeof(Bool_t)*3);
  memcpy(fEnabledTriggerTypes, ref.fEnabledTriggerTypes, sizeof(Bool_t)*6);
}

/**
 * Destructor
 */
AliEMCALTriggerQA::~AliEMCALTriggerQA()
{
}

/**
 * Set the patch types to be plotted
 * \param type Patch type of which the status is being changed
 * \param e    Either enable or disable
 */
void AliEMCALTriggerQA::EnablePatchType(PatchTypes_t type, Bool_t e)
{
  fEnabledPatchTypes[type] = e;
}

/** Set the trigger types to be plotted
 * \param type Trigger type of which the status is being changed
 * \param e    Either enable or disable
 */
void AliEMCALTriggerQA::EnableTriggerType(EMCalTriggerType_t type, Bool_t e)
{
  fEnabledTriggerTypes[type] = e;
}

/**
 * Return the amplitude of the patch (online, recalc, offline)
 * \param patch Pointer to a AliEMCALTriggerPatchInfo object
 * \param itype Amplitude type (online, recalc, offline)
 * \return amplitude
 */
Int_t AliEMCALTriggerQA::GetAmplitude(const AliEMCALTriggerPatchInfo* patch, Int_t itype)
{
  if (!patch) return 0;
  if (itype == 0 || itype == 1) {
    return patch->GetADCAmp();
  }
  else if (itype == 2) {
    return patch->GetADCOfflineAmp();
  }
  else {
    return 0;
  }
}

/**
 * This function should be called every event to set the new time stamp.
 * It sets the time stamp in the internal field, computes the time stamp bin
 * based on fTimeStampBinWidth.
 * \param timeStamp Time stamp of the event
 */
void AliEMCALTriggerQA::EventTimeStamp(UInt_t timeStamp)
{
  fEventTimeStamp = timeStamp;

  if (fTimeStampBinWidth == 0) return;

  UInt_t timeStampBins = fEventTimeStamp / fTimeStampBinWidth;
  fEventTimeStampBin = timeStampBins*fTimeStampBinWidth;
}

/// Actions to be executed only once for the first event
void AliEMCALTriggerQA::ExecOnce()
{
  if (!fGeom) {
    fGeom = AliEMCALGeometry::GetInstance();
    if (!fGeom) {
      AliError("Could not get geometry!");
    }
  }
}
