//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCEVENTSTATISTICS_H
#define ALIHLTTPCEVENTSTATISTICS_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCEventStatistics.h
    @author Jochen Thaeder
    @date   
    @brief  TPC class for event statistics, derived from @see AliHLTEventStatistics
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include "TObject.h"
#include "TString.h"

#include "AliHLTEventStatistics.h"

/**
 * @class  AliHLTTPCEventStatistics
 * @brief  TPC class for event statistics, derived from @see AliHLTEventStatistics
 *
 * The event statistic classes hold information about certain characteristica 
 * of the processed events. They are devided into 3 parts. A base class 
 * @see AliHLTEventStatistics for general Information, detector specific
 * classes like @see AliHLTTPCEventStatistics for the TPC and a summary class
 * @see AliHLTEventStatisticsSummary which can hold several detector classes.
 *
 * This is the detector class for the TPC.
 *
 * See base class @see AliHLTEventStatistics for further information.
 *
 * 
 * @ingroup alihlt_run_statistics alihlt_trigger
 */

class AliHLTTPCEventStatistics : public AliHLTEventStatistics {
  
public:
  
  /** constructor */
  AliHLTTPCEventStatistics();
  /** destructor */
  virtual ~AliHLTTPCEventStatistics();

  // -- Tracks --

  /** Get Total number of tracks 
   *  @return number of tracks
   */
  Int_t GetNTotalTracks()                              { return fNTotalTracks; }

  /** Set Total number of tracks 
   *  @param i  number of tracks
   */
  void SetNTotalTracks( Int_t i)                       { fNTotalTracks = i; }

  // --

  /** Get Total number of tracks with ( number of cluster >= fClusterThreshold )
   *  @return number of tracks
   */
  Int_t GetNTracksAboveClusterThreshold()              { return fNTracksAboveClusterThreshold; }

  /** Set Total number of tracks with ( number of cluster >= fClusterThreshold ) 
   *  @param i  number of tracks
   */
  void SetNTotalTracksAboveClusterThreshold( Int_t i)  { fNTracksAboveClusterThreshold = i;}

  /** Add tracks to total number of tracks with ( number of cluster >= fClusterThreshold ) 
   *  @param i  number of tracks
   */
  void AddNTracksAboveClusterThreshold( Int_t i = 1 )  { fNTracksAboveClusterThreshold += i; }

  // --

  /** Max of ( tracks per sector ) 
   *  @return number of tracks
   */
  Int_t GetNMaxTracksPerSector()                       { return fNMaxTracksPerSector; }

  /** Max of ( tracks per sector ) 
   *  @param i  number of tracks
   */
  void SetNMaxTracksPerSector( Int_t i)                { fNMaxTracksPerSector = i; }

  // --

  /** Min of ( tracks per sector ) 
   *  @return number of tracks
   */
  Int_t GetNMinTracksPerSector()                       { return fNMinTracksPerSector; }

  /** Max of ( tracks per sector ) 
   *  @param i  number of tracks
   */
  void SetNMinTracksPerSector( Int_t i)                { fNMinTracksPerSector = i; }

  // --

  /** Average of ( tracks per sector ) 
   *  @return number of tracks
   */
  Int_t GetNAvgTracksPerSector()                       { return fNAvgTracksPerSector; }

  /** Average of ( tracks per sector ) 
   *  @param i  number of tracks
   */
  void SetNAvgTracksPerSector( Int_t i)                { fNAvgTracksPerSector = i; }

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

  /** Get Total number of cluster found
   *  @return number of cluster
   */
  Int_t GetNTotalCluster()                             { return fNTotalCluster; }

  /** Set Total number of cluster found
   *  @param i  number of cluster
   */
  void SetNTotalCluster( Int_t i)                      { fNTotalCluster = i; }

  /** Add cluster to total number of cluster found
   *  @param i  number of cluster
   */
  void AddNTotalCluster( Int_t i )                     { fNTotalCluster += i; }

  //--

  /** Get Number of cluster which were used in tracks 
   *  @return number of cluster
   */
  Int_t GetNUsedCluster()                              { return fNUsedCluster; }

  /** Set Number of cluster which were used in tracks 
   *  @param i  number of cluster
   */
  void SetNUsedCluster( Int_t i)                       { fNUsedCluster = i; }

  /** Add cluster to total number of cluster which were used in tracks 
   *  @param i  number of cluster
   */
  void AddNUsedCluster( Int_t i )                     { fNUsedCluster += i; }

  // --

  /** Get Average of ( number of clusters per track ), floored 
   *  @return number of cluster
   */
  Int_t GetAvgClusterPerTrack()                        { return fAvgClusterPerTrack; }

  /** Set Average of ( number of clusters per track ), floored 
   *  @param i  number of cluster
   */
  void SetAvgClusterPerTrack( Int_t i)                 { fAvgClusterPerTrack = i; }

private:
 
  /** copy constructor prohibited */
  AliHLTTPCEventStatistics (const AliHLTTPCEventStatistics&);

  /** assignment operator prohibited */
  AliHLTTPCEventStatistics& operator= (const AliHLTTPCEventStatistics&);

  // -- Tracks --
 
  /** Total number of tracks found  */
  Int_t fNTotalTracks;                                         // see above

  /** Tracks with ( number of cluster >= fClusterThreshold ) 
   * ( for long tracks ) 
   */
  Int_t fNTracksAboveClusterThreshold;                         // see above

  /** Max of ( tracks per sector ) */
  Int_t fNMaxTracksPerSector;                                  // see above

  /** Min of ( tracks per sector ) */
  Int_t fNMinTracksPerSector;                                  // see above

  /** Average  of ( tracks per sector ) */
  Int_t fNAvgTracksPerSector;                                  // see above

  // -- Cluster --

  /** Threshold for number of clusters per track 
   * ( for long tracks ) 
   */
  Int_t fClusterThreshold;                                     // see above

  /** Total number of cluster found  */
  Int_t fNTotalCluster;                                        // see above

  /** Number of cluster which were used in tracks */
  Int_t fNUsedCluster;                                         // see above

  /** Average of ( number of clusters per track ), floored */
  Int_t fAvgClusterPerTrack;                                   // see above

  ClassDef(AliHLTTPCEventStatistics, 0);

};
#endif

