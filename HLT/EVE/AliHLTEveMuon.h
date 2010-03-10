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
  /** Process block containing tracks */
  void ProcessTracks(AliHLTHOMERBlockDesc * block, TEveStraightLineSet * tracks);
  
  /** create the cluster pointset*/
  TEvePointSet * CreateClusters();
  /** create the tracks lineset */
  TEveStraightLineSet * CreateTrackSet();

  TEveStraightLineSet * fTracks; //The track elements
  TEvePointSet * fClusters; //The cluster elements

  ClassDef(AliHLTEveMuon, 0);
};

#endif
