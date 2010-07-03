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
//	Create some histograms and one tree with the information
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

#ifndef AliAnalysisTaskAcorde_cxx
#define AliAnalysisTaskAcorde_cxx

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
  virtual void   CreateOutputObjects(); 	//! Creates output object (cosmicTree)
  virtual void   Exec(Option_t *option); 	//! Execution class
  virtual void   Terminate(Option_t *); 	//! Terminate class 

 private:


  AliESDEvent *fESD;    		//! ESD object
  TArrayI      *fPair;			//! Pair track connected (up & down)
  TTree *cosmicTree;			//! TTree with some information of matched tracks
  Int_t nTracks;			//! # of recontructed tracks
  Int_t nMatchTracks;			//! # of matched tracks


  // Cut definitions

  const Int_t minTPCclusters; 		//! cut in clusters
  const Float_t minTrackDist;		//! cut in distance
  const Float_t minCutDir;		//! minimum cut

  Float_t xAco;				//! x-coordinate of extrapolated track to ACORDE
  Float_t yAco;				//! y-coordinate of extrapolated track to ACORDE
  Float_t zAco;				//! z-coordinate of extrapolated track to ACORDE
  Int_t trigger;			//! trigger label
  TString ActiveTriggerDetector;	//! detector string	

  Int_t nSingleTracks;			//! no. of single track
  Int_t nMatchedTracks;			//! no. of matched track
  

  TList *histo;				//! list of histograms
  TH1F *acordeHitsAll;			//! hits of acorde
  TH1F *acordeMultiAll;			//! multi. of acorde modules
  TH1F *acordeHitsTPC;			//! hits of acorde (if track)
  TH1F *acordeMultiTPC;			//! multi. of acorde modules (if track)
  TH1F *nTracksMatched;			//! matched tracks
  TH1F *ntracks;			//! no. of tracks
  TH2F *fTracksToAcorde;		//! tracks extrapolated to acorde.

  AliAnalysisTaskAcorde(const AliAnalysisTaskAcorde&); // not implemented
  AliAnalysisTaskAcorde& operator=(const AliAnalysisTaskAcorde&); // not implemented
  
  ClassDef(AliAnalysisTaskAcorde, 1); // example of analysis
};

#endif
