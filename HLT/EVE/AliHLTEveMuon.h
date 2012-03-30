/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliHLTEveCalo.h
/// @author Svein Lindal
/// @brief  Muon Instance of Eve display processor


#ifndef ALIHLTEVEMUON_H
#define ALIHLTEVEMUON_H

#include "AliHLTEveBase.h"
class AliHLTHOMERBlockDesc;
class TEveStraightLineSet;
class TEvePointSet;
class TEveTrackList;
class AliMUONTrack;
class AliHLTMUONTrackStruct;

class AliHLTEveMuon : public AliHLTEveBase {

public:
  
  /** Constructor  **/
  AliHLTEveMuon();

  /** Destructor **/
 ~AliHLTEveMuon();

  /** Inherited form AliHLTEveBase */
  void ProcessBlock(AliHLTHOMERBlockDesc * block);

  /** inherited from AliHLTEveBase */
  void UpdateElements();
  
  /** inherited from AliHLTEveBase */
  void ResetElements();

private:
  
  /** copy constructor prohibited */
  AliHLTEveMuon(const AliHLTEveMuon&);
  /** assignment operator prohibited */
  AliHLTEveMuon& operator = (const AliHLTEveMuon &);

  /** Inherited from AliHLTEveBase */
  void ProcessHistogram(AliHLTHOMERBlockDesc * block );

  /** Process block containing clusters */
  void ProcessClusters(AliHLTHOMERBlockDesc * block, TEvePointSet * clusters);
  /** Process block containing Manso tracks */
  void ProcessTracks(AliHLTHOMERBlockDesc * block, TEveStraightLineSet * tracks);
  /** Process block containing Full Tracks **/
  Int_t ProcessFullTracks(AliHLTHOMERBlockDesc * block, TEveTrackList * tracks);
  
  /** Convert muon Full Tracks block to Muon tracks **/
  int MakeMUONTrack(AliMUONTrack *muonTrack, const AliHLTMUONTrackStruct *muonHLTTrack);

  /** create the cluster pointset**/
  TEvePointSet * CreateClusters();
  /** create the Manso tracks lineset **/
  TEveStraightLineSet * CreateTrackSet();
  /** create the tracks for Full Tracker **/
  TEveTrackList * CreateFullTrackList();

  TEveTrackList *  fFullTrackList; //The track elements for Full Tracker
  TEveStraightLineSet * fTracks; //The track elements for Manso Tracker
  TEvePointSet * fClusters; //The cluster elements

  ClassDef(AliHLTEveMuon, 0);
};

#endif
