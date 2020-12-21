/**************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                   *
 * All rights reserved.                                                               *
 *                                                                                    *
 * Redistribution and use in source and binary forms, with or without                 *
 * modification, are permitted provided that the following conditions are met:        *
 *     * Redistributions of source code must retain the above copyright               *
 *       notice, this list of conditions and the following disclaimer.                *
 *     * Redistributions in binary form must reproduce the above copyright            *
 *       notice, this list of conditions and the following disclaimer in the          *
 *       documentation and/or other materials provided with the distribution.         *
 *     * Neither the name of the <organization> nor the                               *
 *       names of its contributors may be used to endorse or promote products         *
 *       derived from this software without specific prior written permission.        *
 *                                                                                    *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND    *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED      *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE             *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY                *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES         *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;       *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND        *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT         *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 **************************************************************************************/
#ifndef ALIEMCALTRIGGERREJECTIONMAKER_H
#define ALIEMCALTRIGGERREJECTIONMAKER_H

#include "AliAnalysisTaskEmcal.h"
#include <TString.h>

class THistManager;

/**
 * @namespace PWG
 * @brief Namespace for PWG framework classes
 */
namespace PWG {

/**
 * @namespace EMCAL
 * @ingroup EMCALFW
 * @brief Namespace for EMCAL framework classes and task
 */
namespace EMCAL {

/**
 * @class AliEmcalTriggerRejectionMaker
 * @ingroup EMCALFWTASKS
 * @brief Simple task monitoring the energy spectrum of the maximum patch in the event.
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * @since Aug 30th, 2016
 *
 * With the help of the maximum trigger patch energy spectrum trigger rejection factors
 * can be calculated as function of the trigger threshold.
 */
class AliEmcalTriggerRejectionMaker: public AliAnalysisTaskEmcal {
public:

  /**
   * Dummy (I/O) constructor
   */
  AliEmcalTriggerRejectionMaker();

  /**
   * Main constructor, intializing the task
   * @param[in] name Name of the task
   */
  AliEmcalTriggerRejectionMaker(const char *name);

  /**
   * Destructor, cleaning up
   */
  virtual ~AliEmcalTriggerRejectionMaker();

  /**
   * Specify trigger to be selected. By default INT7 is used
   * as trigger.
   * @param[in] trigger Trigger (after physics selection) to be selected
   */
  void SetTriggerSelection(ULong_t trigger) { fSelectTrigger = trigger; };

  /**
   * Apply additional event selection using the trigger string
   * (for EMCAL level1 triggers)
   * @param[in] pattern Trigger string pattern to be requested
   */
  void SetTriggerPattern(const TString &pattern) { fTriggerPattern = pattern; }

protected:

  /**
   * Creating output histograms. Also initializing the
   * analysis utils in case they are not provided from
   * outside.
   */
  virtual void               UserCreateOutputObjects();

  /**
   * Perform event selection. In this task the standard
   * ALICE event selection for pPb events is used
   * @return True if the event is selected, false otherwise
   */
  virtual Bool_t             IsEventSelected();

  /**
   * Select the highest EGA and EJE patch per event and
   * fill the histogram for the corresponding spectrum with
   * the energy measurement.
   * @return Always true
   */
  virtual Bool_t             Run();

private:
  THistManager                          *fHistos;           //!<! Local Histogram handler. Histograms are stored in the AliEmcalList later
  ULong_t                               fSelectTrigger;     ///< Trigger bit selection
  TString                               fTriggerPattern;    ///< Trigger pattern string (i.e. EG1)

  AliEmcalTriggerRejectionMaker(const AliEmcalTriggerRejectionMaker &);
  AliEmcalTriggerRejectionMaker &operator=(const AliEmcalTriggerRejectionMaker &);

  /// \cond CLASSIMP
  ClassDef(AliEmcalTriggerRejectionMaker, 1);
  /// \endcond
};

} /* namespace EMCAL */

} /* namespace PWG */

#endif /* ALIEMCALTRIGGERREJECTIONMAKER_H */
