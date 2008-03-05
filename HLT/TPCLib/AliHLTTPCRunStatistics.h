//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCRUNSTATISTICS_H
#define ALIHLTTPCRUNSTATISTICS_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCRunStatistics.cxx
    @author Jochen Thaeder
    @date   
    @brief  TPC class for event statistics, derived from @see AliHLTRunStatistics
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include "TObject.h"
#include "TString.h"

#include "AliHLTRunStatistics.h"

/**
 * @class  AliHLTTPCRunStatistics
 * @brief  TPC class for event statistics, derived from @see AliHLTRunStatistics
 *
 * The run statistic classes hold information / histograms about certain 
 * characteristica of the processed events. They are devided into 3 parts. A base class 
 * @see AliHLTRunStatistics for general Information, detector specific
 * classes like @see AliHLTTPCRunStatistics for the TPC and a summary class
 * @see AliHLTRunStatisticsSummary which can hold several detector classes.
 *
 * This is the detector class for the TPC. 
 *
 * See base class @see AliHLTRunStatistics for further information.
 *
 * @ingroup alihlt_run_statistics alihlt_trigger
 */

class AliHLTTPCRunStatistics : public AliHLTRunStatistics {
  
public:
  
  /** constructor */
  AliHLTTPCRunStatistics();
  /** destructor */
  virtual ~AliHLTTPCRunStatistics();

  // -- Tracks --
   
  /** Set Total number of tracks found
   *  @param i number of tracks
   */
  void SetNTotalTracks( Int_t i  )    { fNTotalTracks = (ULong64_t) i; }

  /** Add to Total number of tracks found
   *  @param i number of tracks
   */
  void AddNTotalTracks( Int_t i  )    { fNTotalTracks += (ULong64_t) i; }

  /** Get Total number of tracks found
   *  @return number of tracks
   */
  ULong64_t GetNTotalTracks()    { return fNTotalTracks; }

  // --

  /** Set Total number of tracks found
   *  @param i number of tracks
   */
  void SetNTracksAboveClusterThreshold( Int_t i  )    { fNTracksAboveClusterThreshold = (ULong64_t) i; }

  /** Add to Tracks with ( number of cluster >= fClusterThreshold ) 
   *  @param i number of tracks
   */
  void AddNTracksAboveClusterThreshold( Int_t i  )    { fNTracksAboveClusterThreshold += (ULong64_t) i; }

  /** Get Tracks with ( number of cluster >= fClusterThreshold ) 
   *  @return number of tracks
   */
  ULong64_t GetNTracksAboveClusterThreshold()    { return fNTracksAboveClusterThreshold; }

  // -- Cluster --

  /** Get Threshold for number of clusters per track 
   *  @return threshold
   */
  Int_t GetClusterThreshold()                          { return fClusterThreshold; }

  /** Set Threshold for number of clusters per track 
   *  @param i  threshold
   */
  void SetClusterThreshold( Int_t i)                   { fClusterThreshold = i; }

  // --
  
  /** Set Total number of cluster found
   *  @param i number of cluster
   */
  void SetNTotalCluster( Int_t i  )    { fNTotalCluster = (ULong64_t) i; }

  /** Add to Total number of cluster found
   *  @param i number of cluster
   */
  void AddNTotalCluster( Int_t i  )    { fNTotalCluster += (ULong64_t) i; }

  /** Get Total number of cluster found
   *  @return number of cluster
   */
  ULong64_t GetNTotalCluster()    { return fNTotalCluster; }

  // --
  
  /** Set Number of cluster which were used in tracks
   *  @param i number of cluster
   */
  void SetNUsedCluster( Int_t i  )    { fNUsedCluster = (ULong64_t) i; }

  /** Add to Number of cluster which were used in tracks
   *  @param i number AvgClusterPerTrackof cluster
   */
  void AddNUsedCluster( Int_t i  )    {fNUsedCluster  += (ULong64_t) i; }

  /** Get Number of cluster which were used in tracks
   *  @return number of cluster
   */
  ULong64_t GetNUsedCluster()    { return fNUsedCluster; }

  // --
  
  /** Set Average of ( number of clusters per track ), floored
   *  @param i number of cluster
   */
  void SetAvgClusterPerTrack( Int_t i  )    { fAvgClusterPerTrack = (ULong64_t) i; }

  /** Add to Average of ( number of clusters per track ), floored
   *  @param i number of cluster
   */
  void AddAvgClusterPerTrack( Int_t i  )    { fAvgClusterPerTrack += (ULong64_t) i; }

  /** Get Average of ( number of clusters per track ), floored
   *  @return number of cluster
   */
  ULong64_t GetAvgClusterPerTrack()    { return fAvgClusterPerTrack; }

private:
 
  /** copy constructor prohibited */
  AliHLTTPCRunStatistics (const AliHLTTPCRunStatistics&);

  /** assignment operator prohibited */
  AliHLTTPCRunStatistics& operator= (const AliHLTTPCRunStatistics&);

  // -- Tracks --

  /** Total number of tracks found  */
  ULong64_t fNTotalTracks;                                         // see above

  /** Tracks with ( number of cluster >= fClusterThreshold ) 
   * ( for long tracks ) 
   */
  ULong64_t fNTracksAboveClusterThreshold;                         // see above

  // -- Cluster --

  /** Threshold for number of clusters per track 
   * ( for long tracks ) 
   */
  Int_t fClusterThreshold;                                        // see above

  /** Total number of cluster found  */
  ULong64_t fNTotalCluster;                                        // see above

  /** Number of cluster which were used in tracks */
  ULong64_t fNUsedCluster;                                         // see above

  /** Average of ( number of clusters per track ), floored */
  ULong64_t fAvgClusterPerTrack;                                   // see above

  ClassDef(AliHLTTPCRunStatistics, 0);

};
#endif

