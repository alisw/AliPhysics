#ifndef ALIMUONTRACK_H
#define ALIMUONTRACK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/
// Revision of includes 07/05/2004

/// \ingroup rec
/// \class AliMUONTrack
/// \brief Reconstructed track in ALICE dimuon spectrometer
///
////////////////////////////////////////////////////
/// Reconstructed track in ALICE dimuon spectrometer
////////////////////////////////////////////////////

#include "AliMUONTrackParam.h" // object belongs to the class
#include <TClonesArray.h>

class AliMUONHitForRec;
class AliMUONObjectPair;

class AliMUONTrack : public TObject 
{
 public:
  AliMUONTrack(); // Default constructor
  AliMUONTrack(AliMUONObjectPair *segment); // Constructor from a segment
  virtual ~AliMUONTrack(); // Destructor
  AliMUONTrack (const AliMUONTrack& track); // copy constructor
  AliMUONTrack& operator=(const AliMUONTrack& track); // assignment operator


	/// return pointeur to track parameters at vertex
  AliMUONTrackParam*         GetTrackParamAtVertex() {return &fTrackParamAtVertex;}
	/// set track parameters at vertex
  void                       SetTrackParamAtVertex(const AliMUONTrackParam* trackParam) {fTrackParamAtVertex = *trackParam;}

	/// return array of track parameters at hit
  TClonesArray*              GetTrackParamAtHit() const {return fTrackParamAtHit;}
	/// reset array of track parameters at hit
  void                       ResetTrackParamAtHit() { fTrackParamAtHit->Delete(); }
  void                       AddTrackParamAtHit(const AliMUONTrackParam *trackParam, AliMUONHitForRec *hitForRec); 
  void                       RemoveTrackParamAtHit(AliMUONTrackParam *trackParam);
  void                       UpdateTrackParamAtHit();
  void                       UpdateCovTrackParamAtHit();
  
	/// return array of hitForRec at hit
  TClonesArray*              GetHitForRecAtHit() const {return fHitForRecAtHit;}
	/// reset array of hitForRec at hit
  void                       ResetHitForRecAtHit() { fHitForRecAtHit->Delete(); }
  void                       AddHitForRecAtHit(const AliMUONHitForRec *hitForRec); 

	/// return the number of hits attached to the track
  Int_t                      GetNTrackHits() const {return fNTrackHits;}
	/// set the number of hits attached to the track
  void                       SetNTrackHits(Int_t nTrackHits) {fNTrackHits = nTrackHits;}

	/// return kTrue if the vertex must be used to constrain the fit, kFalse if not
  Bool_t                     GetFitWithVertex() const {return fFitWithVertex;}
	/// set the flag telling whether the vertex must be used to constrain the fit or not
  void                       SetFitWithVertex(Bool_t fitWithVertex) { fFitWithVertex = fitWithVertex; }
	/// return the vertex used during the tracking procedure
  AliMUONHitForRec*          GetVertex() const {return fVertex;}
  void                       SetVertex(const AliMUONHitForRec* vertex);

	/// return kTrue if the multiple scattering must be accounted for in the fit, kFalse if not
  Bool_t                     GetFitWithMCS() const {return fFitWithMCS;}
	/// set the flag telling whether the multiple scattering must be accounted for in the fit or not
  void                       SetFitWithMCS(Bool_t fitWithMCS) {fFitWithMCS = fitWithMCS;}
  
  Bool_t    ComputeHitWeights(TMatrixD* mcsCovariances = 0);
  Bool_t    ComputeLocalChi2(Bool_t accountForMCS);
  Double_t  ComputeGlobalChi2(Bool_t accountForMCS);

	/// return the minimum value of the function minimized by the fit
  Double_t                   GetFitFMin() const {return fGlobalChi2;}
	/// set the minimum value of the function minimized by the fit
  void                       SetFitFMin(Double_t chi2) { fGlobalChi2 = chi2;}
  
  	/// return kTRUE if the track has been improved
  Bool_t                     IsImproved() const {return fImproved;}
	/// set the flag telling whether the track has been improved or not
  void                       SetImproved(Bool_t improved) { fImproved = improved;}
  
        /// return 1,2,3 if track matches with trigger track, 0 if not
  Int_t                      GetMatchTrigger(void) const {return fMatchTrigger;}
  /// returns the local trigger number corresponding to the trigger track 
  Int_t                      GetLoTrgNum(void) const {return floTrgNum;}
	/// set the flag telling whether track matches with trigger track or not
   void			     SetMatchTrigger(Int_t matchTrigger) {fMatchTrigger = matchTrigger;}
   /// set the local trigger number corresponding to the trigger track
   void			     SetLoTrgNum(Int_t loTrgNum) {floTrgNum = loTrgNum;}
	/// return the chi2 of trigger/track matching 
  Double_t                   GetChi2MatchTrigger(void) const {return fChi2MatchTrigger;}
	/// set the chi2 of trigger/track matching 
  void                       SetChi2MatchTrigger(Double_t chi2MatchTrigger) {fChi2MatchTrigger = chi2MatchTrigger;}

  Int_t                      HitsInCommon(AliMUONTrack* track) const;

  Double_t                   GetNormalizedChi2() const;

  Bool_t*                    CompatibleTrack(AliMUONTrack* track, Double_t sigma2Cut) const; // return array of compatible chamber
  
	/// return track number in TrackRefs
  Int_t                      GetTrackID() const {return fTrackID;}
	/// set track number in TrackRefs
  void                       SetTrackID(Int_t trackID) {fTrackID = trackID;}

        /// set word telling which trigger chambers where hit by track
  UShort_t                   GetHitsPatternInTrigCh() const {return fHitsPatternInTrigCh;}
        /// set word telling which trigger chambers where hit by track
  void                       SetHitsPatternInTrigCh(UShort_t hitsPatternInTrigCh) {fHitsPatternInTrigCh = hitsPatternInTrigCh;}

  /// set local trigger information for the matched trigger track
  void SetLocalTrigger(Int_t loCirc, Int_t loStripX, Int_t loStripY, Int_t loDev, Int_t loLpt, Int_t loHpt);
  /// return local trigger information for the matched trigger track
  Int_t GetLocalTrigger(void) const { return fLocalTrigger;              }
  /// number of triggering circuit
  Int_t LoCircuit(void) const { return fLocalTrigger & 0xFF;       }
  /// x-strip local trigger 
  Int_t LoStripX(void) const  { return fLocalTrigger >>  8 & 0x1F; }
  /// y-strip local trigger 
  Int_t LoStripY(void) const  { return fLocalTrigger >> 13 & 0x0F; }
  /// deviation local trigger 
  Int_t LoDev(void)    const  { return fLocalTrigger >> 17 & 0x1F; }
  /// low pt decision local trigger 
  Int_t LoLpt(void)    const  { return fLocalTrigger >> 22 & 0x03; }
  /// high pt decision local trigger 
  Int_t LoHpt(void)    const  { return fLocalTrigger >> 24 & 0x03; }

  void                       RecursiveDump(void) const; // Recursive dump (with track hits)

  virtual void               Print(Option_t* opt="") const;

  virtual void Clear(Option_t* opt="");


 private:
 
  AliMUONTrackParam fTrackParamAtVertex; //!< Track parameters at vertex
  TClonesArray *fTrackParamAtHit; ///< Track parameters at hit
  TClonesArray *fHitForRecAtHit; ///< Cluster parameters at hit
  Int_t fNTrackHits; ///< Number of hits attached to the track
  
  Bool_t fFitWithVertex; //!< kTRUE if using the vertex to constrain the fit, kFALSE if not
  AliMUONHitForRec *fVertex; //!< Vertex used during the tracking procedure if required
  
  Bool_t fFitWithMCS; //!< kTRUE if accounting for multiple scattering in the fit, kFALSE if not
  
  TMatrixD* fHitWeightsNonBending; //!< weights matrix, in non bending direction, of hits attached to the track
  				   //!< (accounting for multiple scattering and hits resolution)
  TMatrixD* fHitWeightsBending;    //!< weights matrix, in bending direction, of hits attached to the track
  				   //!< (accounting for multiple scattering and hits resolution)
  
  Double_t fGlobalChi2; ///< Global chi2 of the track
  
  Bool_t fImproved; //!< kTRUE if the track has been improved
  
  Int_t fMatchTrigger;  ///<  0 track does not match trigger
                        ///<  1 track match but does not pass pt cut
                        ///<  2 track match Low pt cut
                        ///<  3 track match High pt cut
  Int_t floTrgNum; ///< the number of the corresponding loTrg, -1 if no matching
  Double_t fChi2MatchTrigger; ///< chi2 of trigger/track matching 
  
  Int_t fTrackID; ///< track ID = track number in TrackRefs
  UShort_t fHitsPatternInTrigCh; ///< Word containing info on the hits left in trigger chambers

  Int_t fLocalTrigger;    ///< packed local trigger information
  
  
  // methods
  Bool_t ComputeHitWeights(TMatrixD& hitWeightsNB, TMatrixD& hitWeightsB, TMatrixD* mcsCovariances = 0, AliMUONHitForRec* discardedHit = 0) const;
  void   ComputeMCSCovariances(TMatrixD& mcsCovariances) const;
  
  
  ClassDef(AliMUONTrack, 7) // Reconstructed track in ALICE dimuon spectrometer
};
	
#endif
