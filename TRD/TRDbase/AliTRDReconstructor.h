#ifndef ALITRDRECONSTRUCTOR_H
#define ALITRDRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class for TRD reconstruction                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliReconstructor.h"
#include "AliRecoParam.h"
#include "AliDetectorRecoParam.h"
#include "AliTRDpidUtil.h"
#include "AliTRDrecoParam.h"
#include "AliTRDdigitsParam.h"
#include "AliESDTrdTrigger.h"

class TClonesArray;
class TTreeSRedirector;
class AliRawReader;
class AliTRDclusterizer;
class AliTRDonlineTrackMatching;
class AliTRDReconstructor: public AliReconstructor 
{
public:
  enum ETRDReconstructorSteer {
    kDigitsConversion= BIT(0)
    ,kWriteClusters  = BIT(1)
    ,kWriteTracklets = BIT(2)
    ,kSeeding        = BIT(3)
    ,kHLT            = BIT(4)
    ,kProcTracklets  = BIT(5) // process online tracklets
    ,kDebug          = BIT(6)
    ,kClRadialCorr   = BIT(7) // toggle radial correction in clusters
    ,kOwner          = BIT(8)
    ,kNsteer         = 8      // number of tasks
  };
  AliTRDReconstructor();
  virtual ~AliTRDReconstructor();
	
	virtual void        Init();

  virtual void        ConvertDigits(AliRawReader *rawReader, TTree *digitsTree) const;
  virtual AliTracker* CreateTracker() const;
  TTreeSRedirector*   GetDebugStream(AliTRDrecoParam::ETRDReconstructionTask task) const { return task < AliTRDrecoParam::kTRDreconstructionTasks ? fDebugStream[task] : NULL; }

  virtual void        FillESD(AliRawReader *, TTree *clusterTree, AliESDEvent *esd) const { FillESD((TTree * )NULL, clusterTree, esd);                    }
  virtual void        FillESD(TTree *digitsTree, TTree *clusterTree, AliESDEvent *esd) const;
  static TClonesArray* GetClusters();
  static TClonesArray* GetTracklets(const char *trkltype = "");
  static TClonesArray* GetTracks();
  static Int_t        GetNTimeBins()             { return fgNTimeBins;}
  Int_t               GetNdEdxSlices() const     { return (Int_t)AliTRDpidUtil::GetNdEdxSlices(GetPIDMethod());}
  AliTRDpidUtil::ETRDPIDMethod       GetPIDMethod() const       { return GetRecoParam()->IsPIDNeuralNetwork() ? AliTRDpidUtil::kNN : AliTRDpidUtil::kLQ;}
  static const AliTRDrecoParam* GetRecoParam()   { return dynamic_cast<const AliTRDrecoParam*>(AliReconstructor::GetRecoParam(2)); }
  static Float_t      GetMinClustersInTrack()    { return fgkMinClustersInTrack;}
  static Float_t      GetLabelFraction()         { return fgkLabelFraction;}
  static Double_t     GetMaxChi2()               { return fgkMaxChi2;}
  static Double_t     GetMaxSnp()                { return fgkMaxSnp;}
  static Double_t     GetMaxStep()               { return fgkMaxStep;}
  static Double_t     GetEpsilon()               { return fgkEpsilon;}

  virtual Bool_t      HasDigitConversion() const { return fSteerParam&kDigitsConversion;  };
  Bool_t              IsCosmic() const           { return GetRecoParam()->GetEventSpecie() & AliRecoParam::kCosmic;}
  Bool_t              IsWritingClusters() const  { return fSteerParam&kWriteClusters;}
  Bool_t              IsWritingTracklets() const { return fSteerParam&kWriteTracklets;}
  Bool_t              IsHLT() const              { return fSteerParam&kHLT;}
  Bool_t              IsSeeding() const          { return fSteerParam&kSeeding;}
  Bool_t              IsProcessingTracklets() const { return fSteerParam&kProcTracklets;}
  Bool_t              IsDebugStreaming() const   { return (fSteerParam&kDebug || AliTRDReconstructor::GetStreamLevel()>0);}
  Bool_t              UseClusterRadialCorrection() const { return fSteerParam&kClRadialCorr;}

  static void         Options(UInt_t steer=0);
  virtual void        Reconstruct(AliRawReader *rawReader, TTree *clusterTree) const;
  virtual void        Reconstruct(TTree *digitsTree, TTree *clusterTree) const;

 static void         SetClusters(TClonesArray *clusters)  { fgClusters = clusters;}
 static void         SetTracklets(TClonesArray *tracklets) { fgTracklets = tracklets;}
 static void         SetTracks(TClonesArray *tracks) { fgTracks = tracks;}
  void	              SetOption(Option_t *opt);
  static Int_t GetStreamLevel()               { return fgStreamLevel;}
  static void  SetStreamLevel(Int_t level) { fgStreamLevel = level;}
  //
  static void SetExtraMaxClPerLayer(Int_t n)      {if (n>0) fgExtraMaxClPerLayer = n;}
  static void SetExtraBoundaryTolerance(double v) {fgExtraBoundaryTolerance =v;}
  static void SetExtraRoadY(double v)             {fgExtraRoadY = v;}
  static void SetExtraRoadZ(double v)             {fgExtraRoadZ = v;}
  static void SetExtraChi2Out(double v)           {fgExtraChi2Out = v;}
  //
  static Int_t    GetExtraMaxClPerLayer()     {return fgExtraMaxClPerLayer;}
  static Double_t GetExtraBoundaryTolerance() {return fgExtraBoundaryTolerance;}
  static Double_t GetExtraRoadY()             {return fgExtraRoadY;}
  static Double_t GetExtraRoadZ()             {return fgExtraRoadZ;}
  static Double_t GetExtraChi2Out()           {return fgExtraChi2Out;}
  //
private:
  AliTRDReconstructor(const AliTRDReconstructor &r); //Not implemented
  AliTRDReconstructor& operator = (const AliTRDReconstructor&); //Not implemented
  void                ResetContainers() const;
  static Int_t               fgStreamLevel; // flag for streaming      - for TRD reconstruction

  static Char_t const *fgSteerNames[kNsteer];//! steering names
  static Char_t const *fgSteerFlags[kNsteer];//! steering flags
  static Char_t const   *fgTaskNames[AliTRDrecoParam::kTRDreconstructionTasks]; //! tasks names
  static Char_t const   *fgTaskFlags[AliTRDrecoParam::kTRDreconstructionTasks]; //! tasks flags
  UInt_t            fSteerParam;          // steering bits
  // configuration vars for tracking
  static const Double_t    fgkMaxChi2;                  // Max increment in track chi2
  static const Float_t     fgkMinClustersInTrack;       // Min number of clusters in track
  static const Float_t     fgkLabelFraction;            // Min fraction of same label
  static const Double_t    fgkMaxSnp;                   // Maximal snp for tracking
  static const Double_t    fgkMaxStep;                  // Maximal step for tracking
  static const Double_t    fgkEpsilon;                  // Precision of radial coordinate

  TTreeSRedirector *fDebugStream[AliTRDrecoParam::kTRDreconstructionTasks];// Debug Streamer container;
 
  static TClonesArray *fgClusters;    //  list of clusters for local reconstructor
  static TClonesArray *fgTracklets;   //  list of online tracklets for local reconstructor
  static TClonesArray *fgTracks;      //  list of GTU tracks for local reconstructor
  static Int_t         fgNTimeBins;   //  number of time bins as given by the clusterizer
  AliTRDclusterizer   *fClusterizer;  //! instance of TRD clusterizer
  static AliTRDonlineTrackMatching fgOnlineTrackMatcher; // track matcher between on-line and off-line track
  static AliESDTrdTrigger fgTriggerFlags; //  L1 trigger flags
  //
  static Double_t fgExtraBoundaryTolerance; // additional tolerance for boundary check
  static Double_t fgExtraRoadY;             // additional road in Y
  static Double_t fgExtraRoadZ;             // additional road in Z
  static Double_t fgExtraChi2Out;           // additional chi2 on backpropagation
  static Int_t    fgExtraMaxClPerLayer;     // additional cl. per layer allowed
  //
  ClassDef(AliTRDReconstructor, 6)    //  Class for the TRD reconstruction

};

#endif

