// clang-format off
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
// clang-format on
#ifndef ALIHLTEMCALDIGITSMONITORCOMPONENT
#define ALIHLTEMCALDIGITSMONITORCOMPONENT

#include <AliHLTCaloProcessor.h>

class AliHLTEMCALDigitsMonitor;
class AliHLTEMCALGeometry;

/**
 * @class   AliHLTEMCALDigitsMonitorComponent
 * @brief   HLT Component steering monitoring for EMCAL on digits level
 * @author  Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since   June 1, 2017
 *
 * This class steers the monitring for EMCAL on digits level on the HLT. As
 * HLT component is only responsible for input and output. For the monitoring
 * logics see @ref AliHLTEMCALDigitsMonitor.
 */
class AliHLTEMCALDigitsMonitorComponent : public AliHLTCaloProcessor
{
 public:
  /** Constructor */
  AliHLTEMCALDigitsMonitorComponent();

  /** Destructor */
  virtual ~AliHLTEMCALDigitsMonitorComponent();

  /** interface function, see @ref AliHLTComponent for description */
  const char* GetComponentID();

  /** interface function, see @ref AliHLTComponent for description */
  void GetInputDataTypes(std::vector<AliHLTComponentDataType>& list);

  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponentDataType GetOutputDataType();

  /** interface function, see @ref AliHLTComponent for description */
  void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);

  /** interface function, see @ref AliHLTComponent for description */
  int DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
              AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* /*outputPtr*/, AliHLTUInt32_t& /*size*/,
              std::vector<AliHLTComponentBlockData>& /*outputBlocks*/);

  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponent* Spawn();

 protected:
  /** interface function, see @ref AliHLTComponent for description */
  int DoInit(int argc, const char** argv);
  int DoDeinit() { return 0; };

  /** interface function, see @ref AliHLTComponent for description */
  virtual int Deinit(); ////////// PTH WARNING

  /** Push all histograms in the list */
  void PushHistograms(TCollection* list);

  /** Check whether a given block is of the correct input data type */
  bool CheckInputDataType(const AliHLTComponentDataType& datatype);

  /** Initialize Geometry */
  void InitialiseGeometry();

 private:
  Int_t fLocalEventCount;                   ///< Loacl event counter
  Bool_t fHistoResetOnPush;                 ///< Reset histogams after pushing
  AliHLTEMCALDigitsMonitor* fDigitsMonitor; //!<! transient
  AliHLTEMCALGeometry* fGeometry;           //!<! EMCAL geometry

  AliHLTEMCALDigitsMonitorComponent(const AliHLTEMCALDigitsMonitorComponent&);
  AliHLTEMCALDigitsMonitorComponent& operator=(const AliHLTEMCALDigitsMonitorComponent&);
};
#endif
