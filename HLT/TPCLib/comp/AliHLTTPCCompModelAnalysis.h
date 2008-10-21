// XEmacs -*-C++-*-
// $Id$

#ifndef ALIHLTTPCCOMPMODELANALYSIS_H
#define ALIHLTTPCCOMPMODELANALYSIS_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCCompModelAnalysis.h
    @author J. Wagner
    @date   17-11-2007
    @brief   A HLT processing component for the Vestbo-compression-model */

#include "AliHLTTPCTrackArray.h"
#include "AliHLTTPCTrackletDataFormat.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCTrack.h"
#include "AliHLTTPCModelTrack.h"
#include "TFile.h"
#include "TH1.h"
#include "AliHLTLogging.h"

/**
 * @class AliHLTTPCCompModelAnalysis
 * @brief A HLT processing component for the Vestbo-compression-model 
 *
 * An implementiation of a class that is used to analyse the 
 * loss due to the conversion / compression in the Vestbo-model
 * resolution change due to Vestbo-model on tracking performance
 * @ingroup alihlt_tpc
 */
class AliHLTTPCCompModelAnalysis: public AliHLTLogging
{
public:


  /** type needed to build list with discarded tracks or tracks to compare */
  struct AliHLTTPCTrackList
  {
    AliHLTTPCTrackList() : 
      track(),
      pythiatrack(),
      wronglydiscarded(kFALSE),
      matchingindicator(0),
      next(NULL),
      matchingtrack(NULL)
    {}

    AliHLTTPCTrackList (const AliHLTTPCTrackList&);
    AliHLTTPCTrackList& operator= (const AliHLTTPCTrackList&);

    AliHLTTPCTrack track;       // store information of found discarded track   
    AliHLTTPCTrack pythiatrack; // store pythia information about this found discarded track
    Bool_t wronglydiscarded;    // flag to mark if track and pythia track information match together
    Int_t matchingindicator;    // only for trackanalysis the higher the number, the more probable it is that tracks match
    AliHLTTPCTrackList* next;   // pointer to next struct
    AliHLTTPCTrackList* matchingtrack; // pointer to matching track (only used in trackanalysis)
  };
  typedef struct AliHLTTPCTrackList AliHLTTPCTrackList;

  AliHLTTPCCompModelAnalysis(Bool_t modelanalysis, Bool_t trackanalysis, TString dumpfilename, TString graphfilename);
  virtual ~AliHLTTPCCompModelAnalysis();

  /** initialise track or cluster arrays depending on the analysis to be made
   * @return 0 upon success
   */  
  Int_t Init();

  /** display results for respective analysis type (track comparison or model loss analysis) 
   * @return 0 upon success
   */
  Int_t DisplayResults();

  /** function to retrieve private member variable fModelAnalysis used in ModelConverter 
   * @return 0 if fModelAnalysis is switched off
   * @return 1 if fModelAnalysis is switched on
   */
  Bool_t GetfModelAnalysis() {return fModelAnalysis;};

  /** function to retrieve private member variable fTrackAnalysis used in ModelConverter
   * @return 0 if fTrackAnalysis is switched off
   * @return 1 if fTrackAnalysis is switched on
   */
  Bool_t GetfTrackAnalysis() {return fTrackAnalysis;};

  /** fill track arrays with track data from original and secondary tracking 
   * @param tracklets           pointer to track array to be filled
   * @param fillingfirsttracks  boolean to decide which track array is to be filled (1 for first, 0 for second)
   * @return 0 upon success
   */
  Int_t SetTracks(AliHLTTPCTrackletData* tracklets, Bool_t fillingfirsttracks);

  /** fill cluster arrays with cluster data from original and secondary clusters 
   * @param clusters              pointer to cluster data array to be filled
   * @param slice                 slice number where clusters come from
   * @param patch                 patch number where clusters come from
   * @param fillingfirstclusters  boolean to decide which track array is to be filled (1 for first, 0 for second)
   * @return 0 upon success
   */
  Int_t SetClusters(AliHLTTPCClusterData* clusters, UInt_t slice, UInt_t patch, Bool_t fillingfirstclusters);

 /** store discarded tracks in model analysis to be displayed in @ref DisplayModelResults(),
  * uses @ref GetTrashTrackPythiaInfo
  * @param lowpttrack  pointer to discarded track which is to be remembered for analysis
  * @return 0 upon success
  */
  Int_t MarkTrashTrack(AliHLTTPCTrack* lowpttrack);

  /** store discarded cluster not assigned to any track in model analysis to be displayed in @ref DisplayModelResults(),
   * uses @ref GetClusterPythiaInfo
   * @param discardedcluster  cluster data of discarded cluster
   * @param slice             slice where discarded cluster occurrs
   * @param patch             patch where discarded cluster occurrs
   * @return 0 upon success
   */
  Int_t MarkTrashCluster(AliHLTTPCClusterData *discardedcluster, UInt_t slice, UInt_t patch);
 
private:
  /** copy constructor prohibited */
  AliHLTTPCCompModelAnalysis (const AliHLTTPCCompModelAnalysis&); 

  /** assignment operator prohibited */
  AliHLTTPCCompModelAnalysis& operator= (const AliHLTTPCCompModelAnalysis&);

  /** private function to display results from model loss analysis
   * @return 0 upon success
   */
  Int_t DisplayModelResults();

  /** private function to display results from track comparison 
   * @return 0 upon success
   */
  Int_t DisplayTrackResults();

 /** compare tracks and store differences to be displayed in @ref DisplayTrackResults() 
  * @return 0 upon success
  */
  Int_t CompareTracks();

  /** compare clusters for all tracks and create graphs out of the differences 
   *  function used in @ref DisplayTrackResults()
   * @return 0 upon success
   */
  Int_t CompareClusters(Bool_t relativedifferences = 1);

  /** get Pythia information about tracks in track comparison
   * @param comparabletrack track to look for pythia information
   * @return pythiatrack    track information from pythia lookup
   */
  AliHLTTPCTrack GetComparableTrackPythiaInfo(AliHLTTPCTrack comparabletrack);

  /** compare discarded track parameters with parameters from Pythia event
   * @param discardedtrack  pointer to a discarded track (usually with low pt)
   * @return 0 upon correct decision (track with low pt accoridng to Pythia, i.e. track = delta-electron or similar noise)
   * @return 1 upon wrong decision (track wrongly discarded, should by taken into account according to Pythia information)
   */
  Bool_t GetTrashTrackPythiaInfo(AliHLTTPCTrack* discardedtrack);

  /** compare information of a cluster not assigned to any track with its Pythia information
   * @param discardedcluster  pointer to discarded cluster
   * @return 0 upon correct decision (cluster not assigned to any track is true in Pythia, i.e. cluster = noise cluster)
   * @return 1 upon wrong decision (cluster wrongly discarded, i.e. it belongs to a valuable track according to Pythia)
   */
  Bool_t GetClusterPythiaInfo(AliHLTTPCClusterData* discardedcluster);

  /** compare two tracks in order to find if the match
   * @param firsttracklistelement   track from orignal tracking, stored in AliHLTTPCTrackList
   * @param secondtracklistelement  track from secondary tracking, stored in AliHLTTPCTrackList
   * @return matchingindicator      if all parameters are equal: matchingindicator = 10
   */
  Int_t CompareTrackInfo(AliHLTTPCTrackList* firsttracklistelement, AliHLTTPCTrackList* secondtracklistelement);

  /** compare two pythia tracks in order to find if the match
   * @param firsttracklistelement   track from orignal tracking, stored in AliHLTTPCTrackList
   * @param secondtracklistelement  track from secondary tracking, stored in AliHLTTPCTrackList
   * @return matchingindicator      if all pythia parameters are equal: matchingindicator = 10 (should be so!)
   */
  Int_t ComparePythiaTrackInfo(AliHLTTPCTrackList* firsttracklistelement, AliHLTTPCTrackList* secondtracklistelement);

  /** if -graphfile filename.root is given as input parameter, histrograms are created
   * @param relativedifferences boolean to decide whether to plot histograms 
   *                            with relative differences in track paramters (1) or not (0), 1 by default
   * @return 0 upon success
   */
  Int_t CreateGraphs(Bool_t relativedifferences = 1);

  /** flag to decide wheter to do track or model loss analysis */
  Bool_t fModelAnalysis; // switch on model analysis
  Bool_t fTrackAnalysis; // switch on track analysis

  /** name of humanly readable dump file for analysis results */
  TString fDumpFileName;
  /** name of root file name to store histograms of track comparison */
  TString fGraphFileName;
 
  /** members for track analysis:  */
  /** array with original tracks from first tracking */
  AliHLTTPCTrackArray fFirstTrackArray;  // array to store tracks of first tracking
  /** array with tracks from secondary tracking with Vestbo-decompressed clusters */
  AliHLTTPCTrackArray fSecondTrackArray; // array to store tracks of second tracking (after compression/expansion of model)
  /** pointer to first element of first track list containing tracks and their pythia information  */
  AliHLTTPCTrackList* fFirstTrackList;
  /** pointer to first element of second track list containing tracks and their pythia information */
  AliHLTTPCTrackList* fSecondTrackList;
  /** array of original clusters for deviation analysis to secondary clusters */
  AliHLTTPCClusterData* fOriginalClusters[36][6];
  /** array of secondary clusters for deviation analysis to origin clusters */
  AliHLTTPCClusterData* fSecondaryClusters[36][6];
  /** number of tracks with pt < 0.1GeV in first array */
  Int_t fFirstTrashTracks;
  /** number of tracks with pt < 0.1GeV in second array */
  Int_t fSecondTrashTracks;
  /** total number of compared tracks */
  Int_t fTotalComparedTracks;
  /** number of matched tracks under 0.1 GeV in first array */
  Int_t fMatchedFirstTrashTracks;
  /** number of matched tracks under 0.1 GeV in second array */
  Int_t fMatchedSecondTrashTracks;
  /** number of original tracks that do not have secondary track matches */
  Int_t fFirstUnmatchedTracks;
  /** number of secondary tracks that do not have original track matches */
  Int_t fSecondUnmatchedTracks;
  /** tolerance limits of relative deviation of compared track quantities, set in Init() */
  Float_t fToleranceDeviation;  // (firstparameter - secondparameter)/firstparameter <= fToleranceDeviation

  /** members for model analysis: */
  /** pointer to list of discarded tracks */
  AliHLTTPCTrackList* fTrackListPointer; // pointer to list of discarded tracks
  /** array of discarded clusters */
  AliHLTTPCClusterData* fDiscardedClusters[36][6]; // array of discarded valuable clusters
  
  /** total number of discarded clusters not assigned to any track per event */
  Int_t fTotalDiscardedClusters;    // number of discarded clusters
  /** total number of wrongly discarded clusters according to Pythia information per event */
  Int_t fValuableDiscardedClusters; // number of discarded clusters which should be assigned to a track 
  /** total number of discarded tracks per event */
  Int_t fTrashTracks;               // number of discarded tracks (pt lower 0.1 GeV)

  ClassDef(AliHLTTPCCompModelAnalysis, 1)
  
};
#endif
