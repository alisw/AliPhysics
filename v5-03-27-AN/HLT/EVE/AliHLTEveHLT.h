//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTEVEHLT_H
#define ALIHLTEVEHLT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTEveHLT.h
/// @author Svein Lindal <slindal@fys.uio.no>
/// @brief  HLT Instance of Eve display processor

#include "AliHLTEveBase.h"
class TEvePointSet;
class AliESDtrack;
class AliEveTrack;
class TEveTrackList;
class TEveTrackPropagator;
class TString;
class AliExternalTrackParam;
class TH1F;
class TH1;
class AliESDEvent;

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

  void ProcessEsdEvent( AliESDEvent * esd);

  void DestroyOldTrackList();
  static void * DestroyGarbage(void * arg);

private:

  /** copy constructor prohibited */
  AliHLTEveHLT(const AliHLTEveHLT&);
  /** assignment operator prohibited */
  AliHLTEveHLT& operator = (const AliHLTEveHLT& );

  /*Create the pointset for the display */
  void CreateTrackList();

  // Make a standard track representation and put it into given container.
  // Choose which parameters to use a track's starting point.
  // If gkFixFailedITSExtr is TRUE (FALSE by default) and
  // if ITS refit failed, take track parameters at inner TPC radius.
  AliEveTrack * MakeEsdTrack(AliESDtrack *at, TEveTrackList* cont);

  // Process the ESD block and call the functions necessary to fill the tracklist
  void ProcessEsdBlock( AliHLTHOMERBlockDesc * block );

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

  //Process trigger block
  void ProcessGlobalTrigger(AliHLTHOMERBlockDesc * block);

  //Create tpc qa histograms
  void CreateHistograms();

  //Create eve tracks and put them in track list
  void FillTrackList(AliESDtrack * esdTrack);

  //Create track lists
  void CreateTrackLists();

  //Draw tpc qa histograms
  void DrawHistograms();
  void FillHistograms(AliESDtrack * esdTrack);


  //Add histograms to a canvas in pad number cdCount
  void AddHistogramToCanvas(TH1 * histogram, TCanvas * canvas, Int_t &cdCount);

  //Get color from pt
  Color_t GetColor(Float_t pt);
  //Get color bin track belongs to
  Int_t GetColorBin(Float_t pt);


  Bool_t fTrueField;        //Use true field?
  Bool_t fUseIpOnFailedITS; // Use IP as origin if ITS refit fails?
  Bool_t fUseRkStepper;    // Use Runge Kutta for something something?

  TEveTrackList * fTrackList;  //Eve tracklist 
  TEveElementList * fTrackLists; //Holder for tracklists
  TEveTrackList * fOldTrackList;  //Eve tracklist 
  TEvePointSet * fPointSetVertex;      //Display primary vertex

  TCanvas * fTrCanvas;  //Canvas for track qa histos
  TCanvas * fVertexCanvas; //Other canvas
  
  TH1F * fHistEta;     //Eta histo
  TH1F * fHistPhi;     //Phi histo
  TH1F * fHistnClusters;//nClusters histo
  TH1F * fHistMult;    //Mult histo
  TH1F * fHistDCAr;    //DCA r histo
  Int_t fTrCount;
  Int_t fVCount;
  
  Int_t fNTrackBins;
  
  ClassDef(AliHLTEveHLT, 0);
};

#endif
