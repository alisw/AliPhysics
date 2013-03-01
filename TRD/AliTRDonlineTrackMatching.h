#ifndef ALITRDONLINETRACKMATCHING_H
#define ALITRDONLINETRACKMATCHING_H

//
// Track matching between TRD online tracks and ESD tracks.
// Author: Felix Rettig <rettig@compeng.uni-frankfurt.de>

#include "TObject.h"
class TH1;
class AliESDEvent;
class AliExternalTrackParam;
class AliTRDgeometry;
class AliESDtrack;
class AliESDTrdTrack;

//#define TRD_TM_DEBUG

class AliTRDonlineTrackMatching {
 public:
  AliTRDonlineTrackMatching();
  ~AliTRDonlineTrackMatching();

  inline static Short_t TrdLsiSec(Short_t lsi) { return (lsi < 0) ? -1 : (lsi/5); }           // convert linear stack index 0-89 to TRD sector 0-17
  inline static Short_t TrdLsiSi(Short_t lsi) { return (lsi < 0) ? -1 : (lsi%5); }            // convert stack index 0-89 to local stack index 0-4
  inline static Short_t TrdDetLsi(Short_t det) { return det/6; }                              // convert TRD detector/chamber 0-539 index to linear stack index 0-89
  inline static Short_t TrdDetSi(Short_t det) { return (det%30)/6; }                          // convert detector (=chamber) number 0-539 to local stack index 0-4
  inline static Short_t TrdDetLyr(Short_t det) { return det%6; }                              // convert detector (=chamber) number 0-539 to local layer 0-5
  inline static Short_t TrdSecSiLsi(Short_t sec, Short_t si) { return 5*sec + si; }          // convert sector (0-17) and local stack index (0-4) to linear stack index 0-89

  static Short_t EstimateSector(const Double_t globalCoords[3]);
  static Short_t EstimateLayer(Double_t radius);
  static Short_t EstimateLocalStack(const Double_t globalCoords[3]);
  static Short_t EstimateStack(const Double_t globalCoords[3]);

  static Bool_t StackToTrack(const AliExternalTrackParam* track, Short_t &stack, UShort_t &layersWithTracklet, Double_t magFieldinKiloGauss);
  static Bool_t StackToTrack(const AliESDtrack* track, Short_t &stack, UShort_t &layersWithTracklet, Double_t magFieldinKiloGauss);

  static Bool_t TrackPlaneIntersect(AliExternalTrackParam *trk, Double_t pnt[3], Double_t norm[3], Double_t mag);
  Int_t EstimateTrackDistance(AliESDtrack *esd_track, AliESDTrdTrack* gtu_track, Double_t mag, Double_t *ydist, Double_t *zdist);
  static Double_t RateTrackMatch(Double_t distY, Double_t distZ, Double_t rpt, Double_t gpt);

  static void SetEsdTrackCutMinTPCrows(Double_t setting) { fEsdTrackCutMinTPCrows = setting; };
  static void SetEsdTrackCutMinRatioRowsFindableClusters(Double_t setting) { fEsdTrackCutMinRatioRowsFindableClusters = setting; };
  static void SetEsdTrackCutMaxChi2TPCclusters(Float_t setting) { fEsdTrackCutMaxChi2TPCclusters = setting; };
  static void SetEsdTrackCutMaxChi2ITSclusters(Float_t setting) { fEsdTrackCutMaxChi2ITSclusters = setting; };
  static void SetEsdTrackCutMaxDCAtoVertexXY(Float_t setting) { fEsdTrackCutMaxDCAtoVertexXY = setting; };
  static void SetEsdTrackCutMaxDCAtoVertexZ( Float_t setting) { fEsdTrackCutMaxDCAtoVertexZ = setting; };
  static void SetEsdTrackCutITSlayerMask(UShort_t setting) { fEsdTrackCutsITSlayerMask = setting; };
  static void SetEsdTrackCutsChi2TPCconstrainedVsGlobal(Float_t setting) { fEsdTrackVCutsChi2TPCconstrainedVsGlobal = setting; };
  static void SetEsdTrackCutMinimal(Bool_t setting) { fEsdTrackCutMinimal = setting; };
  static void SetEsdTrackCutRequireITSrefit(Bool_t setting) { fEsdTrackCutRequireITSrefit = setting; };
  static void SetEsdTrackCutRequireTPCrefit(Bool_t setting) { fEsdTrackCutRequireTPCrefit = setting; };
  static void SetEsdTrackCutPrim(Bool_t setting) { fEsdTrackCutPrim = setting; };

  static void SetEsdTrackDefaultCuts(const char* cutIdent = "minimal");

  void SetMinMatchRating(Float_t setting) { fMinMatchRating = setting; };
  void SetHistMatchRating(TH1* hist) { fHistMatchRating = hist; };

  static Bool_t AcceptTrack(const AliESDtrack* esdTrack, const AliESDEvent* esdEvent);
  Bool_t ProcessEvent(AliESDEvent *esdEvent);

 protected:

  inline static Double_t PtDiffRel(Double_t refPt, Double_t gtuPt);

  static const unsigned int fgkMaxEsdTracksPerStack = 20000; // max no. of tracks per stack
  static const unsigned int fgkTrdStacks = 90;		     // no. of TRD stacks (0-89)
  static const Float_t     fgkSaveInnerRadius;		     // inner radius for extrapolation
  static const Float_t     fgkSaveOuterRadius;		     // outer radius for extrapolation

  static Float_t fEsdTrackCutMinTPCrows;                   // track cut
  static Float_t fEsdTrackCutMinRatioRowsFindableClusters; // track cut
  static Float_t fEsdTrackCutMaxChi2TPCclusters;	   // track cut
  static Float_t fEsdTrackCutMaxChi2ITSclusters;	   // track cut
  static Float_t fEsdTrackCutMaxDCAtoVertexXY;		   // track cut
  static Float_t fEsdTrackCutMaxDCAtoVertexZ;		   // track cut
  static UShort_t fEsdTrackCutsITSlayerMask;		   // track cut
  static Float_t fEsdTrackVCutsChi2TPCconstrainedVsGlobal; // track cut
  static Float_t fEsdTrackCutPtDCAOfs;			   // track cut
  static Float_t fEsdTrackCutPtDCACoeff;		   // track cut
  static Bool_t fEsdTrackCutMinimal;                       // activate/deactive minimal track cuts
  static Bool_t fEsdTrackCutRequireITSrefit;		   // track cut
  static Bool_t fEsdTrackCutRequireTPCrefit;		   // track cut
  static Bool_t fEsdTrackCutPrim;			   // select primaries or keep secondaries for matching

  AliTRDgeometry* fTRDgeo;           //! TRD geometry instance
  Float_t fMinMatchRating;	     //! min rating to accept pair of matched tracks
  TH1* fHistMatchRating;	     //! optional histogram for match rating for all TRD tracks

  AliTRDonlineTrackMatching(const AliTRDonlineTrackMatching &c);
  AliTRDonlineTrackMatching& operator=(const AliTRDonlineTrackMatching& rhs);
};

#endif
