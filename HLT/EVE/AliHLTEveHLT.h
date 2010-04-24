/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliHLTEveCalo.h
/// @author Svein Lindal
/// @brief  HLT Instance of Eve display processor
// Author: Svein Lindal <slindal@fys.uio.no>

#ifndef ALIHLTEVEHLT_H
#define ALIHLTEVEHLT_H

#include "AliHLTEveBase.h"
class TEvePointSet;
class AliESDtrack;
class AliEveTrack;
class TEveTrackList;
class TEveTrackPropagator;
class TString;
class AliExternalTrackParam;
class TH1F;

class AliHLTEveHLT : public AliHLTEveBase {

public:
  
  /** Constructor  **/
  AliHLTEveHLT();

  /** Destructor **/
 ~AliHLTEveHLT();
  
  /** Inherited from AliHLTEveBase */
  void ProcessBlock(AliHLTHOMERBlockDesc * block);

  /** inherited from AliHLTEveBase */
  void UpdateElements();
  
  /** inherited from AliHLTEveBase */
  void ResetElements();

private:

  /** copy constructor prohibited */
  AliHLTEveHLT(const AliHLTEveHLT&);
  /** assignment operator prohibited */
  AliHLTEveHLT& operator = (const AliHLTEveHLT );

  /*Create the pointset for the display */
  void CreateTrackList();

  // Make a standard track representation and put it into given container.
  // Choose which parameters to use a track's starting point.
  // If gkFixFailedITSExtr is TRUE (FALSE by default) and
  // if ITS refit failed, take track parameters at inner TPC radius.
  AliEveTrack * MakeEsdTrack(AliESDtrack *at, TEveTrackList* cont);

  // Process the ESD block and call the functions necessary to fill the tracklist
  void ProcessEsdBlock( AliHLTHOMERBlockDesc * block, TEveTrackList * cont );

  //Set up the track propagator
  void SetUpTrackPropagator(TEveTrackPropagator* trkProp, Float_t magF, Float_t maxR);

  //Create a title for the track
  TString CreateTrackTitle(AliESDtrack* t);

  ///Create the pointset to display primary vertex
  void CreateVertexPointSet();

  //Add track param to AliEveTrack
  void AddTrackParamToTrack(AliEveTrack* track, const AliExternalTrackParam* tp);

  //Process histogram block
  void ProcessHistograms(AliHLTHOMERBlockDesc * block, TCanvas * canvas);

  //Create tpc qa histograms
  void CreateHistograms();

  //Draw tpc qa histograms
  void DrawHistograms();


  Bool_t fTrueField;        //Use true field?
  Bool_t fUseIpOnFailedITS; // Use IP as origin if ITS refit fails?
  Bool_t fUseRkStepper;    // Use Runge Kutta for something something?

  TEveTrackList * fTrackList;  //Eve tracklist 
  TEvePointSet * fPointSetVertex;      //Display primary vertex

  TCanvas * fTrCanvas;  //Canvas for track qa histos

  TH1F * fHistPt;      //Pt histo 
  TH1F * fHistP;       //p histo
  TH1F * fHistEta;     //Eta histo
  TH1F * fHistTheta;   //Theta histo
  TH1F * fHistPhi;     //Phi histo
  TH1F * fHistnClusters;//nClusters histo
  TH1F * fHistMult;    //Mult histo

  ClassDef(AliHLTEveHLT, 0);
};

#endif
