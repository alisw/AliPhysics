#ifndef ALITRDTRACKERV1_H
#define ALITRDTRACKERV1_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */ 
/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  The TRD tracker                                                       //
//                                                                        //
//  Authors:                                                              //
//    Marian Ivanov <M.Ivanov@gsi.de>                                     //
//    Alex Bercuci <A.Bercuci@gsi.de>                                     //
//    Jouri Belikov <J.Belikov@cern.ch>                                   //
//    Markus Fasel <M.Fasel@gsi.de>                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALITRACKER_H
#include "AliTracker.h"
#endif

#ifndef ALITRDTRACKINGSECTOR_H
#include "AliTRDtrackingSector.h"
#endif

/**************************************************************************
* Class Status see source file                                           *
**************************************************************************/

class TFile;
class TTreeSRedirector;
class TClonesArray;
class TLinearFitter;

class AliRieman;
class AliESDEvent;
class AliCluster;
class AliTrackPoint;

class AliTRDcluster;
class AliTRDseedV1;
class AliTRDtrackingChamber;
class AliTRDchamberTimeBin;
class AliTRDtrackerFitter;
class AliTRDtrackV1;
class AliTRDReconstructor;
class AliTRDtrackerV1 : public AliTracker
{
public:
  enum{
    kOwner = BIT(14)
  };
  enum{
    kMaxLayersPerSector   = 1000
    , kMaxTimeBinIndex    = 216
    , kTrackingSectors    = 18
    , kNTimeBins          = 35
    , kNPlanes            = 6
    , kNSeedPlanes        = 4
    , kMaxTracksStack     = 100
    , kNConfigs           = 15
  };
  AliTRDtrackerV1(AliTRDReconstructor *rec = 0x0);
  virtual ~AliTRDtrackerV1();
  
  //temporary
  AliTRDtrackingSector* GetTrackingSector(Int_t sec) {return &fTrSec[sec];}
  
  Int_t           Clusters2Tracks(AliESDEvent *esd);
  AliCluster*     GetCluster(Int_t index) const;
  AliTRDseedV1*   GetTracklet(Int_t index) const;
  AliKalmanTrack* GetTrack(Int_t index) const;
  TClonesArray*   GetListOfClusters() const  { return fClusters;}
  TClonesArray*   GetListOfTracklets() const { return fTracklets;}
  TClonesArray*   GetListOfTracks() const    { return fTracks;}
  static Int_t    GetNTimeBins()             { return fgNTimeBins;}
  static void     GetExtrapolationConfig(Int_t iconfig, Int_t planes[2]);
  static void     GetSeedingConfig(Int_t iconfig, Int_t planes[4]);
  static TLinearFitter*  GetTiltedRiemanFitter();
  static TLinearFitter*  GetTiltedRiemanFitterConstraint();
  static AliRieman*      GetRiemanFitter();
  static void     FitRieman(AliTRDcluster **clusters, Double_t chi2[2]);
  static Float_t  FitRieman(AliTRDseedV1 *tracklets, Double_t *chi2, Int_t *planes = 0x0);
  static Float_t  FitTiltedRiemanConstraint(AliTRDseedV1 *tracklets, Double_t zVertex);
  static Float_t  FitTiltedRieman(AliTRDseedV1 *tracklets, Bool_t sigError);
  
  static Double_t FitRiemanTilt(const AliTRDtrackV1 *trk, AliTRDseedV1 *tracklets = 0x0, Bool_t err=0, Int_t np = 0, AliTrackPoint *points = 0x0);
  static Double_t FitLine(const AliTRDtrackV1 *trk, AliTRDseedV1 *tracklets = 0x0, Bool_t err=0, Int_t np = 0, AliTrackPoint *points = 0x0);
  static Double_t FitKalman(AliTRDtrackV1 *trk, AliTRDseedV1 *tracklets = 0x0, Bool_t up=0, Int_t np = 0, AliTrackPoint *points = 0x0);

  Bool_t          IsClustersOwner() const    { return TestBit(kOwner);}
  void            SetClustersOwner(Bool_t own=kTRUE) {SetBit(kOwner, own); if(!own) fClusters = 0x0;}

  Int_t           FollowBackProlongation(AliTRDtrackV1 &t);
  Int_t           FollowProlongation(AliTRDtrackV1 &t);
  Int_t           LoadClusters(TTree *cTree);
  Int_t           LoadClusters(TClonesArray *clusters);
  Int_t           PropagateBack(AliESDEvent *event);
  static Int_t    PropagateToX(AliTRDtrackV1 &t, Double_t xToGo, Double_t maxStep);
  Int_t           ReadClusters(TClonesArray* &array, TTree *in) const;
  Int_t           RefitInward(AliESDEvent *event);
  static void     SetNTimeBins(Int_t nTimeBins){fgNTimeBins = nTimeBins; }
  void            SetReconstructor(const AliTRDReconstructor *rec){ fReconstructor = rec; }
  void            UnloadClusters();
//   void            UseClusters(const AliKalmanTrack *t, Int_t from = 0) const;

  class AliTRDLeastSquare{
  public:
    AliTRDLeastSquare();
    ~AliTRDLeastSquare(){};
    
    void		AddPoint(Double_t *x, Double_t y, Double_t sigmaY);
    void		RemovePoint(Double_t *x, Double_t y, Double_t sigmaY);
    void		Eval();
    
    Double_t	GetFunctionParameter(Int_t ParNumber) const {return fParams[ParNumber];}
    Double_t	GetFunctionValue(Double_t *xpos) const;
    void		GetCovarianceMatrix(Double_t *storage) const;
  private:
          AliTRDLeastSquare(const AliTRDLeastSquare &);
    AliTRDLeastSquare& operator=(const AliTRDLeastSquare &);
    Double_t 	fParams[2];						// Fitparameter	
    Double_t 	fCovarianceMatrix[3];			// Covariance Matrix
    Double_t 	fSums[6];						// Sums
  };

protected:
  static Bool_t  AdjustSector(AliTRDtrackV1 *track); 
  Double_t       BuildSeedingConfigs(AliTRDtrackingChamber **stack, Int_t *configs);
  Int_t          BuildTrackingContainers();
  static Float_t CalculateChi2Z(AliTRDseedV1 *tracklets, Double_t offset, Double_t slope, Double_t xref);
  Int_t          Clusters2TracksSM(Int_t sector, AliESDEvent *esd);
  Int_t          Clusters2TracksStack(AliTRDtrackingChamber **stack, TClonesArray *esdTrackList);
  AliTRDseedV1*  GetTracklet(AliTRDtrackV1 *trk, Int_t plane, Int_t &idx);
  Bool_t         GetTrackPoint(Int_t index, AliTrackPoint &p) const;	
  Float_t        GetR4Layer(Int_t ly) const { return fR[ly];}
  Int_t          MakeSeeds(AliTRDtrackingChamber **stack, AliTRDseedV1 *sseed, Int_t *ipar);
  AliTRDtrackV1* MakeTrack(AliTRDseedV1 *seeds, Double_t *params);
  AliTRDtrackV1* SetTrack(AliTRDtrackV1 *track);
  AliTRDseedV1*  SetTracklet(AliTRDseedV1 *tracklet);

private:
  AliTRDtrackerV1(const AliTRDtrackerV1 &tracker);
  AliTRDtrackerV1 &operator=(const AliTRDtrackerV1 &tracker);
  Double_t       	CookLikelihood(AliTRDseedV1 *cseed, Int_t planes[4]);
  Double_t       	CalculateTrackLikelihood(AliTRDseedV1 *tracklets, Double_t *chi2);
  Int_t          	ImproveSeedQuality(AliTRDtrackingChamber **stack, AliTRDseedV1 *tracklet);
  static Float_t	CalculateReferenceX(AliTRDseedV1 *tracklets);
  void        ResetSeedTB();
  Float_t     GetChi2Y(AliTRDseedV1 *tracklets) const;
  Float_t     GetChi2Z(AliTRDseedV1 *tracklets) const;

private:
  const AliTRDReconstructor *fReconstructor;            // reconstructor manager
  AliTRDgeometry      *fGeom;                           // Pointer to TRD geometry
  AliTRDtrackingSector fTrSec[kTrackingSectors];        // Array of tracking sectors;    
  TClonesArray        *fClusters;                       // List of clusters
  TClonesArray        *fTracklets;                      // List of tracklets
  TClonesArray        *fTracks;                         // List of tracks
  Float_t              fR[kNPlanes];                    //! rough radial position of each TRD layer

  // should go to the recoParam
  static const Double_t    fgkMaxChi2;                  // Max increment in track chi2 
  static const Float_t     fgkMinClustersInTrack;       // Min number of clusters in track
  static const Float_t     fgkLabelFraction;            // Min fraction of same label
  static const Double_t    fgkMaxSnp;                   // Maximal snp for tracking
  static const Double_t    fgkMaxStep;                  // Maximal step for tracking  
  
  // stand alone tracking
  static Double_t      fgTopologicQA[kNConfigs];        //  Topologic quality
  Double_t             fTrackQuality[kMaxTracksStack];  //  Track quality 
  Int_t                fSeedLayer[kMaxTracksStack];     //  Seed layer
  AliTRDchamberTimeBin *fSeedTB[kNSeedPlanes]; // seeding time bin planes
  Int_t                fSieveSeeding;                   //! Seeding iterator
  
  static const Double_t fgkX0[kNPlanes];                // default values for the position of anode wire
  static Int_t         fgNTimeBins;                     // Timebins per plane in track prolongation 
  static TLinearFitter *fgTiltedRieman;                 //  Fitter for the tilted Rieman fit without vertex constriant
  static TLinearFitter *fgTiltedRiemanConstrained;      //  Fitter for the tilted Rieman fit with vertex constraint	
  static AliRieman     *fgRieman;                       //  Fitter for the untilted Rieman fit
  
  ClassDef(AliTRDtrackerV1, 3)                          //  TRD tracker - tracklet based tracking

};
#endif
