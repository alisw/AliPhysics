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

//#ifndef ALITRACKER_H
#include "AliTracker.h"
//#endif

//#ifndef ALITRDTRACKINGSECTOR_H
#include "AliTRDtrackingSector.h"
//#endif

//#ifndef ROOT_TMatrixDfwd
#include <TMatrixDfwd.h>
//#endif

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
class AliTRDrecoParam;
class AliTRDtrackerV1 : public AliTracker
{
public:
  enum{
    kOwner            = BIT(14) // owner of clusters
   ,kRemoveContainers = BIT(15) // delete containers after usage
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
  AliTRDtrackerV1(AliTRDReconstructor *rec = NULL);
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
  static Float_t  FitRieman(AliTRDseedV1 *tracklets, Double_t *chi2, Int_t *const planes = NULL);
  static Float_t  FitTiltedRiemanConstraint(AliTRDseedV1 *tracklets, Double_t zVertex);
  static Float_t  FitTiltedRieman(AliTRDseedV1 *tracklets, Bool_t sigError);
  static Double_t FitTiltedRiemanV1(AliTRDseedV1 *tracklets);
  
  static Double_t FitRiemanTilt(const AliTRDtrackV1 *trk, AliTRDseedV1 *tracklets = NULL, Bool_t err=0, Int_t np = 0, AliTrackPoint *points = NULL);
  static Double_t FitLine(const AliTRDtrackV1 *trk, AliTRDseedV1 *tracklets = NULL, Bool_t err=0, Int_t np = 0, AliTrackPoint *points = NULL);
  static Double_t FitKalman(AliTRDtrackV1 *trk, AliTRDseedV1 * const tracklets = NULL, Bool_t up=0, Int_t np = 0, AliTrackPoint *points = NULL);

  Bool_t          IsClustersOwner() const    { return TestBit(kOwner);}
  Bool_t          HasRemoveContainers() const    { return TestBit(kRemoveContainers);}
  void            SetClustersOwner(Bool_t own=kTRUE) {SetBit(kOwner, own); if(!own) fClusters = NULL;}
  void            SetRemoveContainers(Bool_t rm=kTRUE) {SetBit(kRemoveContainers, rm);}

  Int_t           FollowBackProlongation(AliTRDtrackV1 &t);
  Int_t           FollowProlongation(AliTRDtrackV1 &t);
  Int_t           LoadClusters(TTree *cTree);
  Int_t           LoadClusters(TClonesArray *const clusters);
  Int_t           PropagateBack(AliESDEvent *event);
  static Int_t    PropagateToX(AliTRDtrackV1 &t, Double_t xToGo, Double_t maxStep);
  Bool_t          ReadClusters(TTree *in);
  Int_t           RefitInward(AliESDEvent *event);
  static void     SetNTimeBins(Int_t nTimeBins){fgNTimeBins = nTimeBins; }
  void            SetReconstructor(const AliTRDReconstructor *rec) {fkReconstructor = rec;}
  void            UnloadClusters();

  class AliTRDLeastSquare{
  public:
    AliTRDLeastSquare();
    ~AliTRDLeastSquare(){};
    
    void          AddPoint(const Double_t * const x, Double_t y, Double_t sigmaY);
    void          RemovePoint(const Double_t * const x, Double_t y, Double_t sigmaY);
    Bool_t        Eval();
    void          Reset();
    
    Double_t      GetFunctionParameter(Int_t ParNumber) const {return fParams[ParNumber];}
    Double_t      GetFunctionValue(const Double_t * const xpos) const;
    void          GetCovarianceMatrix(Double_t *storage) const;
  private:
    AliTRDLeastSquare(const AliTRDLeastSquare &);
    AliTRDLeastSquare& operator=(const AliTRDLeastSquare &);
    Double_t      fParams[2];           // Fitparameter	
    Double_t      fCovarianceMatrix[3]; // Covariance Matrix
    Double_t      fSums[6];             // Sums
  };

  class AliTRDtrackFitterRieman{
    public:
      AliTRDtrackFitterRieman();
      ~AliTRDtrackFitterRieman();

      Double_t Eval();
      void Reset();

      Double_t GetYat(Double_t x) const;
      Double_t GetDyDxAt(Double_t x) const;
      Double_t GetZat(Double_t x) const;
      Double_t GetDzDx() const { return fParameters[4]; };
      Double_t GetCurvature() const;
      void GetCovAt(Double_t x, Double_t *cov) const;

      void SetRiemanFitter(TLinearFitter *const fitter) { fTrackFitter = fitter; }
      void SetTracklet(Int_t il, AliTRDseedV1 * const tracklet);
      void SetSysClusterError(Double_t err) { fSysClusterError = err; };
    private:
      AliTRDtrackFitterRieman(const AliTRDtrackFitterRieman &);
      AliTRDtrackFitterRieman &operator=(const AliTRDtrackFitterRieman &);
      void UpdateFitters(AliTRDseedV1 * const tracklet);
      Bool_t CheckAcceptable(Double_t offset, Double_t slope);
      Double_t CalculateReferenceX();

      TLinearFitter *fTrackFitter;        // Fitter for linearized track model
      AliTRDLeastSquare *fZfitter;        // Linear fitter in z-Direction
      AliTRDseedV1 *fTracklets[kNPlanes]; // Tracklet container
      TMatrixD *fCovarPolY;               // Polynomial Covariance Matrix Estimation (y-direction)
      TMatrixD *fCovarPolZ;               // Polynomial Covariance Matrix Estimation (z-direction)
      Double_t fXref;                     // Reference x position for fit in z-Direction
      Double_t fSysClusterError;          // Systematic cluster Error
      Double_t fParameters[5];            // Track Model Parameter
      Double_t fSumPolY[5];               // Sums for polynomial Covariance Matrix Estimation (y-direction)
      Double_t fSumPolZ[3];               // Sums for polynomial Covariance Matrix Estimation (z-direction)
  };

protected:
  static Bool_t  AdjustSector(AliTRDtrackV1 *const track); 
  Double_t       BuildSeedingConfigs(AliTRDtrackingChamber **stack, Int_t *configs);
  Int_t          BuildTrackingContainers();
  static Float_t CalculateChi2Z(const AliTRDseedV1 *tracklets, Double_t offset, Double_t slope, Double_t xref);
  Int_t          Clusters2TracksSM(Int_t sector, AliESDEvent *esd);
  Int_t          Clusters2TracksStack(AliTRDtrackingChamber **stack, TClonesArray * const esdTrackList);
  AliTRDseedV1*  GetTracklet(AliTRDtrackV1 *const trk, Int_t plane, Int_t &idx);
  Bool_t         GetTrackPoint(Int_t index, AliTrackPoint &p) const;	
  Float_t        GetR4Layer(Int_t ly) const { return fR[ly];}
  Int_t          MakeSeeds(AliTRDtrackingChamber **stack, AliTRDseedV1 * const sseed, const Int_t * const ipar);
  AliTRDtrackV1* MakeTrack(AliTRDseedV1 * const tracklet);
  AliTRDtrackV1* SetTrack(const AliTRDtrackV1 * const track);
  AliTRDseedV1*  SetTracklet(const AliTRDseedV1 * const tracklet);
  void           UnsetTrackletsTrack(const AliTRDtrackV1 * const track);

private:
  AliTRDtrackerV1(const AliTRDtrackerV1 &tracker);
  AliTRDtrackerV1 &operator=(const AliTRDtrackerV1 &tracker);
  Double_t       CookLikelihood(AliTRDseedV1 *cseed, Int_t planes[4]);
  Double_t       CalculateTrackLikelihood(Double_t *chi2);
  Bool_t         ImproveSeedQuality(AliTRDtrackingChamber **stack, AliTRDseedV1 *tracklet, Double_t &chi2);
  static Float_t	CalculateReferenceX(const AliTRDseedV1 *const tracklets);
  void        ResetSeedTB();
  Float_t     GetChi2Y(const AliTRDseedV1 *const tracklets) const;
  Float_t     GetChi2Z(const AliTRDseedV1 *const tracklets) const;
  Float_t     GetChi2Phi(const AliTRDseedV1 *const tracklets) const;

  const AliTRDReconstructor *fkReconstructor;           // reconstructor manager
  const AliTRDrecoParam     *fkRecoParam;               // reco param for the current event
  AliTRDgeometry      *fGeom;                           // Pointer to TRD geometry
  AliTRDtrackingSector fTrSec[kTrackingSectors];        // Array of tracking sectors;    
  TClonesArray        *fClusters;                       // List of clusters
  TClonesArray        *fTracklets;                      // List of tracklets
  TClonesArray        *fTracks;                         // List of tracks
  TClonesArray        *fTracksESD;                      // List of ESD tracks in current SM
  Float_t              fR[kNPlanes];                    //! rough radial position of each TRD layer

  // stand alone tracking
  static Double_t      fgTopologicQA[kNConfigs];        //  Topologic quality
  Double_t             fTrackQuality[kMaxTracksStack];  //  Track quality 
  Int_t                fSeedLayer[kMaxTracksStack];     //  Seed layer
  AliTRDchamberTimeBin *fSeedTB[kNSeedPlanes]; // seeding time bin planes
  Int_t                fSieveSeeding;                   //! Seeding iterator
  Int_t                fEventInFile;                    //! event in file being tracked (debug purposes)
  
  static const Double_t fgkX0[kNPlanes];                // default values for the position of anode wire
  static Int_t         fgNTimeBins;                     // Timebins per plane in track prolongation 
  static TLinearFitter *fgTiltedRieman;                 //  Fitter for the tilted Rieman fit without vertex constriant
  static TLinearFitter *fgTiltedRiemanConstrained;      //  Fitter for the tilted Rieman fit with vertex constraint	
  static AliRieman     *fgRieman;                       //  Fitter for the untilted Rieman fit
  
  ClassDef(AliTRDtrackerV1, 7)                          //  TRD tracker - tracklet based tracking

};
#endif
