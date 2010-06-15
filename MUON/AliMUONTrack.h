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

#include <TClonesArray.h>
#include <TMatrixD.h>

class AliMUONVCluster;
class AliMUONObjectPair;
class AliMUONTrackParam;

class AliMUONTrack : public TObject 
{
 public:
  AliMUONTrack(); // Default constructor
  AliMUONTrack(AliMUONObjectPair *segment, Double_t bendingVertexDispersion); // Constructor from a segment
  virtual ~AliMUONTrack(); // Destructor
  AliMUONTrack (const AliMUONTrack& track); // copy constructor
  AliMUONTrack& operator=(const AliMUONTrack& track); // assignment operator

  void Reset();
  
  TClonesArray* GetTrackParamAtCluster() const;
  void          AddTrackParamAtCluster(const AliMUONTrackParam &trackParam, AliMUONVCluster &cluster, Bool_t copy = kFALSE); 
  void          RemoveTrackParamAtCluster(AliMUONTrackParam *trackParam);
  Bool_t        UpdateTrackParamAtCluster();
  Bool_t        UpdateCovTrackParamAtCluster();
  
  Bool_t IsValid(UInt_t requestedStationMask, Bool_t request2ChInSameSt45 = kFALSE);
  
  void TagRemovableClusters(UInt_t requestedStationMask);
  
  /// return the number of clusters attached to the track
  Int_t GetNClusters() const {return fTrackParamAtCluster ? fTrackParamAtCluster->GetEntriesFast() : 0;}

  /// return kTrue if the vertex must be used to constrain the fit, kFalse if not
  Bool_t FitWithVertex() const {return fFitWithVertex;}
  /// set the flag telling whether the vertex must be used to constrain the fit or not
  void   FitWithVertex(Bool_t fitWithVertex) { fFitWithVertex = fitWithVertex; }
  /// return the vertex resolution square used during the tracking procedure
  void   GetVertexErrXY2(Double_t &nonBendingErr2, Double_t &bendingErr2) const
         { nonBendingErr2 = fVertexErrXY2[0]; bendingErr2 = fVertexErrXY2[1]; }
  /// set the vertex resolution square used during the tracking procedure
  void   SetVertexErrXY2(Double_t nonBendingErr2, Double_t bendingErr2)
	 { fVertexErrXY2[0] = nonBendingErr2; fVertexErrXY2[1] = bendingErr2; }

  /// return kTrue if the multiple scattering must be accounted for in the fit, kFalse if not
  Bool_t FitWithMCS() const {return fFitWithMCS;}
  /// set the flag telling whether the multiple scattering must be accounted for in the fit or not
  void   FitWithMCS(Bool_t fitWithMCS) {fFitWithMCS = fitWithMCS;}
  
  Bool_t   ComputeClusterWeights(TMatrixD* mcsCovariances = 0);
  Bool_t   ComputeLocalChi2(Bool_t accountForMCS);
  Double_t ComputeGlobalChi2(Bool_t accountForMCS);

  /// return the minimum value of the function minimized by the fit
  Double_t GetGlobalChi2() const {return fGlobalChi2;}
  /// set the minimum value of the function minimized by the fit
  void     SetGlobalChi2(Double_t chi2) { fGlobalChi2 = chi2;}
  
  /// return kTRUE if the track has been improved
  Bool_t IsImproved() const {return fImproved;}
  /// set the flag telling whether the track has been improved or not
  void   SetImproved(Bool_t improved) { fImproved = improved;}
  
  /// return 1,2,3 if track matches with trigger track, 0 if not
  Int_t    GetMatchTrigger(void) const {return fMatchTrigger;}
  /// returns the local trigger number corresponding to the trigger track (obsolete)
  Int_t    GetLoTrgNum(void) const {return LoCircuit();}
  /// set the flag telling whether track matches with trigger track or not
  void	   SetMatchTrigger(Int_t matchTrigger) {fMatchTrigger = matchTrigger;}
  /// return the chi2 of trigger/track matching 
  Double_t GetChi2MatchTrigger(void) const {return fChi2MatchTrigger;}
  /// set the chi2 of trigger/track matching 
  void     SetChi2MatchTrigger(Double_t chi2MatchTrigger) {fChi2MatchTrigger = chi2MatchTrigger;}

  Int_t ClustersInCommon(AliMUONTrack* track, Int_t stMin = 0, Int_t stMax = 4) const;

  Int_t    GetNDF() const;
  Double_t GetNormalizedChi2() const;

  Int_t  FindCompatibleClusters(const AliMUONTrack &track, Double_t sigma2Cut, Bool_t compatibleCluster[10]) const;
  Bool_t Match(AliMUONTrack &track, Double_t sigma2Cut, Int_t &nMatchClusters) const;
  
  /// return pointer to track parameters at vertex (can be 0x0)
  AliMUONTrackParam* GetTrackParamAtVertex() const {return fTrackParamAtVertex;}
  void               SetTrackParamAtVertex(const AliMUONTrackParam* trackParam);

  /// set word telling which trigger chambers where hit by track
  UShort_t GetHitsPatternInTrigCh() const {return fHitsPatternInTrigCh;}
  /// set word telling which trigger chambers where hit by track
  void     SetHitsPatternInTrigCh(UShort_t hitsPatternInTrigCh) {fHitsPatternInTrigCh = hitsPatternInTrigCh;}

  /// set local trigger information for the matched trigger track
  void  SetLocalTrigger(Int_t loCirc, Int_t loStripX, Int_t loStripY, Int_t loDev, Int_t loLpt, Int_t loHpt, UChar_t respWithoutChamber=0);
  /// return local trigger information for the matched trigger track
  Int_t GetLocalTrigger(void) const { return fLocalTrigger;        }
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
  /// Word stating if trigger would be fired without one chamber
  Int_t GetTriggerWithoutChamber(void) const { return fLocalTrigger >> 26 & 0xF; }
  /// Check if trigger would be fired without chamber (ich [0,3])
  Bool_t TriggerFiredWithoutChamber(Int_t ich) const { return GetTriggerWithoutChamber() >> (3 - ich) & 0x1; }

  void  FindMCLabel();
  /// set the corresponding MC track number
  void  SetMCLabel(Int_t label) {fTrackID = label;}
  /// return the corresponding MC track number
  Int_t GetMCLabel() const {return fTrackID;}
  
  void RecursiveDump(void) const; // Recursive dump (with associated clusters)

  virtual void Print(Option_t* opt="") const;

  virtual void Clear(Option_t* opt="");
  
  /// return the maximum chi2 above which the track can be considered as abnormal (due to extrapolation failure, ...)
  static Double_t MaxChi2() {return fgkMaxChi2;}
  
  /// set the flag which is kTRUE if that track shares cluster(s) with another
  void   Connected(Bool_t flag = kTRUE) {fConnected = flag;}
  /// return the flag which is kTRUE if that track shares cluster(s) with another
  Bool_t IsConnected() const {return fConnected;}


 private:
  
  static const Double_t fgkMaxChi2; ///< maximum chi2 above which the track can be considered as abnormal
  
  mutable TClonesArray* fTrackParamAtCluster; ///< Track parameters at cluster
  
  Bool_t   fFitWithVertex;   //!< kTRUE if using the vertex to constrain the fit, kFALSE if not
  Double_t fVertexErrXY2[2]; //!< Vertex resolution square used during the tracking procedure if required
  
  Bool_t fFitWithMCS; //!< kTRUE if accounting for multiple scattering in the fit, kFALSE if not
  
  TMatrixD* fClusterWeightsNonBending; //!< weights matrix, in non bending direction, of clusters attached to the track
				       //!< (accounting for multiple scattering and cluster resolution)
  TMatrixD* fClusterWeightsBending;    //!< weights matrix, in bending direction, of clusters attached to the track
				       //!< (accounting for multiple scattering and cluster resolution)
  
  Double_t fGlobalChi2; ///< Global chi2 of the track
  
  Bool_t fImproved; //!< kTRUE if the track has been improved
  
  Int_t fMatchTrigger;  ///<  0 track does not match trigger
                        ///<  1 track match but does not pass pt cut
                        ///<  2 track match Low pt cut
                        ///<  3 track match High pt cut
  Double_t fChi2MatchTrigger; ///< chi2 of trigger/track matching 
  
  Int_t fTrackID; ///< Point to the corresponding MC track
  
  AliMUONTrackParam* fTrackParamAtVertex; //!< Track parameters at vertex
  
  UShort_t fHitsPatternInTrigCh; ///< Word containing info on the hits left in trigger chambers

  Int_t fLocalTrigger;    ///< packed local trigger information
  
  Bool_t fConnected; ///< kTRUE if that track shares cluster(s) with another
  
  // methods
  Bool_t ComputeClusterWeights(TMatrixD& clusterWeightsNB, TMatrixD& clusterWeightsB,
			       TMatrixD* mcsCovariances = 0, const AliMUONVCluster* discardedCluster = 0) const;
  void   ComputeMCSCovariances(TMatrixD& mcsCovariances) const;
  
  
  ClassDef(AliMUONTrack, 9) // Reconstructed track in ALICE dimuon spectrometer
};
	
#endif
