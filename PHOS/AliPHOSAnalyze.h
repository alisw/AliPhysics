#ifndef ALIPHOSANALYZE_H
#define ALIPHOSANALYZE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Algorythm class to analyze PHOS events    //
//  Yves Schutz            SUBATECH           //
//                                            //  
//                                            //
////////////////////////////////////////////////

// --- ROOT system ---

#include "TFile.h"

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSv0.h"
#include "AliPHOSGeometry.h"

class AliPHOSAnalyze : public TObject {

public:

  AliPHOSAnalyze() ;              // ctor
  AliPHOSAnalyze(Text_t * name) ; // ctor
  virtual ~AliPHOSAnalyze() ;     // dtor

  void AnalyzeOneEvent(Int_t evt = -999) ;  // analyzes a single event ;
  Bool_t Init(Int_t evt) ;                  // does various initialisations
  void DisplayKineEvent(Int_t evt = -999) ; // displays the Kine events in ALICE coordinate 
  void DisplayRecParticles() ;              // displays RecParticles in ALICE coordinate  
  void DisplayRecPoints() ;                 // displays RecPoints in module coordinate  
  void DisplayTrackSegments() ;             // displays TrackSegments in module coordinate  
  Bool_t OpenRootFile(Text_t * name) ;      // opens the root file
 
private:
  
  AliPHOSClusterizer * fClu ;       // a clusterizer 
  Int_t fEvt ;                      // the evt number being processed 
  AliPHOSGeometry * fGeom;          // the PHOS Geometry object
  AliPHOSv0 * fPHOS ;               // the PHOS object from the root file 
  AliPHOSParticleGuesser * fPag ;   // a particle guesser
  AliPHOSReconstructioner * fRec ;  // a reconstructioner  
  TFile * fRootFile ;               // the root file that contains the data
  AliPHOSTrackSegmentMaker * fTrs ; // a tracksegmentmaker ;
public:

ClassDef(AliPHOSAnalyze,1)  // PHOS event analyzis , version 1

};

#endif // AliPHOSANALYZE_H
