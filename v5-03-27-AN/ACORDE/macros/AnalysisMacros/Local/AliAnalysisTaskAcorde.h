/////////////////////////////////////////////////////////////////////////////
//
// AliAnalysisTaskAcorde class
//
// Description:
//
//	Reads the information of ACORDE-ESD
//	Also it is implemented a simple matching between tracks
//	to look for the extrapolation of them to ACORDE modules
//
//	Create some fHistograms and one tree with the information
//	of the matching tracks
//
//  Author: Mario Rodríguez Cahuantzi
//		<mario.rocah@gmail.com>
//		<mrodrigu@mail.cern.ch>
//
//  Created: June 30th. 2010 @ FCFM - BUAP, Puebla, MX
//  Last update: created
//
/////////////////////////////////////////////////////////////////////////////

#ifndef ALIANALYSISTASKACORDE_H
#define ALIANALYSISTASKACORDE_H

class TH2F;
class TH1F;
class AliESDEvent;
class TArrayI;
class TList;
#include "AliAnalysisTask.h"

class AliAnalysisTaskAcorde : public AliAnalysisTask {
 public:
  AliAnalysisTaskAcorde(const char *name = "AliAnalysisTaskAcorde");
  virtual ~AliAnalysisTaskAcorde();		//! Destructor fo task
  
  virtual void   ConnectInputData(Option_t *); 	//! Connects input data to class analysis
  virtual void   CreateOutputObjects(); 	//! Creates output object (fCosmicTree)
  virtual void   Exec(Option_t *option); 	//! Execution class
  virtual void   Terminate(Option_t *); 	//! Terminate class 

 private:


  AliESDEvent *fESD;    		//! ESD object
  TArrayI      *fPair;			//! Pair track connected (up & down)
  TTree *fCosmicTree;			//! TTree with some information of matched tracks
  Int_t fnTracks;			//! # of recontructed tracks
  Int_t fNMatchTracks;			//! # of matched tracks


  // Cut definitions

  const Int_t fkMinTPCclusters; 		//! cut in clusters
  const Float_t fkMinTrackDist;		//! cut in distance
  const Float_t fkMinCutDir;		//! minimum cut

  Float_t fXAco;				//! x-coordinate of extrapolated track to ACORDE
  Float_t fYAco;				//! y-coordinate of extrapolated track to ACORDE
  Float_t fZAco;				//! z-coordinate of extrapolated track to ACORDE
  Int_t fTrigger;			//! fTrigger label
  TString fActiveTriggerDetector;	//! detector string	

  Int_t fNSingleTracks;			//! no. of single track
  Int_t fNMatchedTracks;			//! no. of matched track
  

  TList *fHisto;				//! list of fHistograms
  TH1F *fAcordeHitsAll;			//! hits of acorde
  TH1F *fAcordeMultiAll;			//! multi. of acorde modules
  TH1F *fAcordeHitsTPC;			//! hits of acorde (if track)
  TH1F *fAcordeMultiTPC;			//! multi. of acorde modules (if track)
  TH1F *fNTracksMatched;			//! matched tracks
  TH1F *fNTracks;			//! no. of tracks
  TH2F *fTracksToAcorde;		//! tracks extrapolated to acorde.

  AliAnalysisTaskAcorde(const AliAnalysisTaskAcorde&); // not implemented
  AliAnalysisTaskAcorde& operator=(const AliAnalysisTaskAcorde&); // not implemented
  
  ClassDef(AliAnalysisTaskAcorde, 1); // example of analysis
};

#endif
