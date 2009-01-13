// $Id$
#ifndef ALIHLTGLOBALTRACKMERGER_H
#define ALIHLTGLOBALTRACKMERGER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTGlobalTrackMerger.h
    @author Jacek Otwinowski (Jacek.Otwinowski@gsi.de)
    @date   
    @brief  The HLT TPC merger base class
*/


class AliHLTTPCTrack;
class AliHLTTPCTrackSegmentData;
class AliHLTTPCVertex;
class AliHLTTPCTrackArray;
class AliTRDtrackV1;
class AliESDEvent;
class AliESDVertex;
class AliExternalTrackParam;

class TClonesArray;
class TTreeStream;
class TTreeSRedirector;

#include "AliHLTLogging.h"
#include "AliESDtrack.h"

/** 
 * @class AliHLTGlobalTrackMerger
 * Global track merger for the barrel section.
 *
 * @ingroup alihlt_global
 * @author Jacek.Otwinowski@gsi.de
 */
class AliHLTGlobalTrackMerger : public AliHLTLogging {
public:
  AliHLTGlobalTrackMerger();
  /** destructor */
  virtual ~AliHLTGlobalTrackMerger();

  // load tracks
  Bool_t LoadTracks(TClonesArray *aTRDTracks,AliESDEvent *esdEvent=0);
  Bool_t LoadTracks(AliHLTTPCTrackArray *aTPCTracks,AliESDEvent *esdEvent=0);

  // set matching parameters
  void SetParameter(Double_t maxy=1., Double_t maxz=1., Double_t maxsnp=0.05, Double_t maxtgl=0.1, Double_t signed1Pt=0.001);

  // match tracks
  Bool_t MatchTracks(AliExternalTrackParam *extTPC=0, AliESDtrack *trackTRD=0);

  // merge tracks
  Bool_t Merge(AliESDEvent *esdEvent=0);
  Bool_t MergeTracks(AliESDtrack *trackTPC=0, AliESDtrack *trackTRD=0, AliESDEvent *esdEvent=0);

  // create AliESDtrack objects
  void FillTPCESD(AliHLTTPCTrack* tpcTrack=0, ULong_t flags=AliESDtrack::kTPCin, AliESDEvent* esdEvent=0); 
  void FillTRDESD(AliTRDtrackV1*  trdTrack=0, ULong_t flags=AliESDtrack::kTRDin, AliESDEvent* esdEvent=0);

  // propagate tracks to DCA to primary vertex
  void PropagateTracksToDCA(AliESDEvent *esdEvent=0);

  // Smooth track parameters 
  // (origin Sergey Gorbunov)
  /*
  Bool_t SmoothTracks( Double_t T1[], Double_t C1[], Double_t Chi21, Int_t NDF1,
	       	       Double_t T2[], Double_t C2[], Double_t Chi22, Int_t NDF2,
		       Double_t T [], Double_t C [], Double_t &Chi2, Int_t &NDF,
		       Int_t N );
  */
 
  Bool_t SmoothTracks( const Double_t T1[], const Double_t C1[], Double_t Chi21, Int_t NDF1,
	       	       const Double_t T2[], const Double_t C2[], Double_t Chi22, Int_t NDF2,
		       Double_t T [], Double_t C [], Double_t &Chi2, Int_t &NDF,
		       Int_t N );
  Int_t IndexS(Int_t i, Int_t j) {
     return ( j<=i ) ? i*(i+1)/2+j :j*(j+1)/2+i;
  }
  void MultSSQ( const Double_t *A, const Double_t *B, Double_t *C, Int_t N );
  Bool_t InvertS( Double_t A[], Int_t N );  

private:

  // TRD-TPC matching parameters
  Double_t fMaxY;    //! max Y track (cm)
  Double_t fMaxZ;    //! max Z track (cm)
  Double_t fMaxSnp;  //! max local sine of the track momentum azimuthal angle
  Double_t fMaxTgl;  //! max tangent of the track momentum dip angle
  Double_t fMaxSigned1Pt; //! max 1/pt (1/(GeV/c)) 

  // primary vertex
  AliESDVertex *fVertex; //! event vertex (needed to propagate all tracks to DCA) 

  TTreeSRedirector *fDebugStreamer; //!
  
  /** copy constructor prohibited */
  AliHLTGlobalTrackMerger(const AliHLTGlobalTrackMerger&);
  /** assignment operator prohibited */
  AliHLTGlobalTrackMerger& operator=(const AliHLTGlobalTrackMerger&);

  ClassDef(AliHLTGlobalTrackMerger,1) //Merging base class
};

#endif
