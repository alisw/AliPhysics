//**************************************************************************\
//* This file is property of and copyright by the ALICE Project            *\
//* ALICE Experiment at CERN, All rights reserved.                         *\
//*                                                                        *\
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *\
//*                  for The ALICE HLT Project.                            *\
//*                                                                        *\
//* Permission to use, copy, modify and distribute this software and its   *\
//* documentation strictly for non-commercial purposes is hereby granted   *\
//* without fee, provided that the above copyright notice appears in all   *\
//* copies and that both the copyright notice and this permission notice   *\
//* appear in the supporting documentation. The authors make no claims     *\
//* about the suitability of this software for any purpose. It is          *\
//* provided "as is" without express or implied warranty.                  *\
//**************************************************************************

/// \file GPUTRDTrackerComponent.h
/// \brief A TRD tracker processing component for the GPU

/// \author Ole Schmidt

#ifndef GPUTRDTRACKERCOMPONENT_H
#define GPUTRDTRACKERCOMPONENT_H

#ifndef GPUCA_ALIROOT_LIB
#define GPUCA_ALIROOT_LIB
#endif

#include "AliHLTProcessor.h"
#include "AliHLTComponentBenchmark.h"
#include "AliHLTDataTypes.h"

class TH1F;
class TList;

#include "GPUTRDDef.h"
namespace GPUCA_NAMESPACE
{
namespace gpu
{
class GPUTRDGeometry;
} // namespace gpu
} // namespace GPUCA_NAMESPACE

class GPUTRDTrackerComponent : public AliHLTProcessor
{
 public:
  /*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

  /** constructor */
  GPUTRDTrackerComponent();

  /** dummy copy constructor, defined according to effective C++ style */
  GPUTRDTrackerComponent(const GPUTRDTrackerComponent&);

  /** dummy assignment op, but defined according to effective C++ style */
  GPUTRDTrackerComponent& operator=(const GPUTRDTrackerComponent&);

  /** destructor */
  virtual ~GPUTRDTrackerComponent();

  /*
 * ---------------------------------------------------------------------------------
 * Public functions to implement AliHLTComponent's interface.
 * These functions are required for the registration process
 * ---------------------------------------------------------------------------------
 */

  /** interface function, see @ref AliHLTComponent for description */
  const char* GetComponentID();

  /** interface function, see @ref AliHLTComponent for description */
  void GetInputDataTypes(vector<AliHLTComponentDataType>& list);

  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponentDataType GetOutputDataType();

  /** @see component interface @ref AliHLTComponent::GetOutputDataType */
  int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);

  /** interface function, see @ref AliHLTComponent for description */
  void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);

  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponent* Spawn();

  int ReadConfigurationString(const char* arguments);

 protected:
  /*
 * ---------------------------------------------------------------------------------
 * Protected functions to implement AliHLTComponent's interface.
 * These functions provide initialization as well as the actual processing
 * capabilities of the component.
 * ---------------------------------------------------------------------------------
 */

  // AliHLTComponent interface functions

  /** interface function, see @ref AliHLTComponent for description */
  int DoInit(int argc, const char** argv);

  /** interface function, see @ref AliHLTComponent for description */
  int DoDeinit();

  /** interface function, see @ref AliHLTComponent for description */
  int DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks);

  /** interface function, see @ref AliHLTComponent for description */
  int Reconfigure(const char* cdbEntry, const char* chainId);

  ///////////////////////////////////////////////////////////////////////////////////

 private:
  /*
 * ---------------------------------------------------------------------------------
 * Private functions to implement AliHLTComponent's interface.
 * These functions provide initialization as well as the actual processing
 * capabilities of the component.
 * ---------------------------------------------------------------------------------
 */

  /*
 * ---------------------------------------------------------------------------------
 *                              Helper
 * ---------------------------------------------------------------------------------
 */

  /*
 * ---------------------------------------------------------------------------------
 *                             Members - private
 * ---------------------------------------------------------------------------------
 */
  GPUCA_NAMESPACE::gpu::GPUTRDTracker* fTracker; // the tracker itself
  GPUCA_NAMESPACE::gpu::GPUTRDGeometry* fGeo;    // TRD geometry needed by the tracker

  TList* fTrackList;
  bool fDebugTrackOutput;              // output GPUTRDTracks instead AliHLTExternalTrackParam
  bool fVerboseDebugOutput;            // more verbose information is printed
  bool fRequireITStrack;               // only TPC tracks with ITS match are used as seeds for tracking
  AliHLTComponentBenchmark fBenchmark; // benchmark

  ClassDef(GPUTRDTrackerComponent, 0);
};
#endif // GPUTRDTRACKERCOMPONENT_H
