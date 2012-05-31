// -*- Mode: C++ -*-
// $Id$
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file  AliHLTGlobalTrackResidualsComponent.h 
/// @author Timur Pocheptsov
/// @date  
/// @brief  A histogramming component for plotting the Y and Z track residuals
///         

/**
 * @class AliHLTGlobalTrackResidualsComponent
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b GlobalTrackResiduals  <br>
 * Library: \b libAliHLTGlobal.so         <br>
 * Input Data Types:                      <br>
 * Output Data Types:                     <br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *      
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * Configuration by component arguments.
 *
 * <h2>Default CDB entries:</h2>
  *
 * <h2>Performance:</h2>
 * The component does not process any event data.
 *
 * <h2>Memory consumption:</h2>
 * The component does not process any event data.
 *
 * @ingroup alihlt_global_components
 */

#ifndef ALIHLTGLOBALTRACKRESIDUALSCOMPONENT_H
#define ALIHLTGLOBALTRACKRESIDUALSCOMPONENT_H

#include <utility>
#include <vector>

#include <TH1F.h>

#include "AliHLTProcessor.h"

class AliHLTTPCSpacePointData;
class AliHLTGlobalBarrelTrack;

class AliHLTGlobalTrackResidualsComponent : public AliHLTProcessor {
private:
  enum EDefaults {
    kNBins = 100
  };

public:
  /** default constructor */
  AliHLTGlobalTrackResidualsComponent();

  //Overriders. These functions are required for the registration process
  /** Component's id - "name" - "GlobalTrackResiduals" */
  const char* GetComponentID();
  /** Types of input data blocks */
  void GetInputDataTypes(AliHLTComponentDataTypeList& list);
  /** The type of output data */
  AliHLTComponentDataType GetOutputDataType();
  /** Types of output data */
  int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);
  /** Approximate size of output */
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  /** "Virtual constructor" to create a component */
  AliHLTComponent* Spawn();

protected:

  //Overriders. Do component's work.
  /** Reset histograms */
  int DoInit(int argc, const char** argv);
  /** Do nothing at the moment */
  int DoDeinit();
  /** Process input - clusters and global barrel tracks */
  int DoEvent(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
  /** Supress warning from compiler about the hidden name */
  using AliHLTProcessor::DoEvent;

private:

  /** Read input data, find residuals, fill histgrams */
  void ProcessBlocks();
  /** Read input data, extract clusters */
  void ReadClusterBlocks();
  /** Extract hits' Xs for the track, sort them. */
  void SortHitsX(const AliHLTGlobalBarrelTrack& gt);
  /** Find residuals and fill histograms */
  void FillResiduals(const AliHLTGlobalBarrelTrack& gt);
  /** Clean fClustersArray and fNSpacePoints */
  void CleanClusters();
  /** Reset histograms - ranges, bins */
  void ResetHistograms();

  TH1F fResY;
  TH1F fResZ;

  std::vector<std::pair<Float_t, UInt_t> > fSortedX; //! Hits, sorted along X
  const AliHLTTPCSpacePointData*           fClustersArray[36][6]; //! Clusters for event.
  UInt_t                                   fNSpacePoints[36][6];  //! Number of points in a cluster for event.

  /** Non-copyable class */
  AliHLTGlobalTrackResidualsComponent(const AliHLTGlobalTrackResidualsComponent& rhs);
  /** Non-copyable class */
  AliHLTGlobalTrackResidualsComponent& operator = (const AliHLTGlobalTrackResidualsComponent& rhs);

  ClassDef(AliHLTGlobalTrackResidualsComponent, 0);//Component to calculate residuals.
};

#endif
