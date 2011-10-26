#ifndef ALITRDCALIBTASK_H
#define ALITRDCALIBTASK_H

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration task for offline calibration                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


class TList;
class TObject;
class TObjArray;
class TH2F;
class TH1F;
class TH1I;
class TH2S;
class TProfile2D;
class TH2I;
class TTree;
class AliESDEvent;
class AliESDfriend;
class AliESDtrack;
class AliESDfriendTrack;
class AliTRDtrackV1;
class AliTRDseedV1;
class AliTRDCalibraFillHisto;
class AliTRDcluster;
class AliESDtrackCuts;
class AliTRDCalDet;
class AliTRDCalibChamberStatus;

#include "TObjString.h"
#include "AliAnalysisTaskSE.h" 
#include "TMath.h"

class AliTRDCalibTask : public AliAnalysisTaskSE {
 public:
  AliTRDCalibTask(const char *name = "AliTRDCalibTask");
  virtual ~AliTRDCalibTask();
  
  //  virtual void   ConnectInputData(Option_t *);
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *);
  virtual Bool_t Load(const Char_t *filename);
  virtual Bool_t Load(TList *lister);
  virtual Long64_t  Merge(TCollection *li);
  void           AddTask(const AliTRDCalibTask * calibTask);
  Bool_t         FindP1TrackPHtrackletV1Test(const AliTRDseedV1 *tracklet, Int_t nbclusters);
  TList          *GetList() const {return fListHist;};

  void SetOnInstance(Bool_t onInstance)                             {fOnInstance=onInstance;};
  void SetHisto2d(Bool_t histo2d)                                   {fHisto2d=histo2d;};
  void SetVector2d(Bool_t vector2d)                                 {fVector2d=vector2d;};
  void SetVdriftLinear(Bool_t vdriftLinear)                         {fVdriftLinear = vdriftLinear;};
  void SetNbTimeBins(Int_t nbTimeBins)                              {fNbTimeBins=nbTimeBins;};  
  
  void SetNz(Short_t nz, Int_t i)                                      {fNz[i]=nz;};
  void SetNrphi(Short_t nrphi, Int_t i)                                {fNrphi[i]=nrphi;};
  
  void AddSelectedTriggerClass(const char*name)                        {fSelectedTrigger->Add(new TObjString(name));};
  void SetReject(Bool_t rejected)                                      {fRejected = rejected;};
  
  void SetESDtrackCuts(AliESDtrackCuts * const esdtrackCuts)           {fEsdTrackCuts = esdtrackCuts;};
  void SetRequirePrimaryVertex(Bool_t requirePrimaryVertex)            {fRequirePrimaryVertex = requirePrimaryVertex;};
  void SetUseTPCVertex()                                               {fVtxTPC=kTRUE ; fVtxSPD=kFALSE;} 
  void SetUseSPDVertex()                                               {fVtxTPC=kFALSE; fVtxSPD=kTRUE ;} 
  void SetMinNbOfContributors(Int_t minNbOfContributors)               {fMinNbContributors = minNbOfContributors;};  
  void SetRangePrimaryVertexZ(Double_t rangePrimaryVertexZ)            {fRangePrimaryVertexZ = TMath::Abs(rangePrimaryVertexZ);}; 
  void SetMinNbTracks(Int_t minNbTracks)                               {fMinNbTracks = minNbTracks;};
  void SetMaxNbTracks(Int_t maxNbTracks)                               {fMaxNbTracks = maxNbTracks;};
 
  void SetVersionGainUsed(Int_t versionGainUsed)                       { fVersionGainUsed = versionGainUsed;   }
  void SetSubVersionGainUsed(Int_t subVersionGainUsed)                 { fSubVersionGainUsed = subVersionGainUsed;   }
  void SetVersionGainLocalUsed(Int_t versionGainLocalUsed)             { fVersionGainLocalUsed = versionGainLocalUsed;   }
  void SetSubVersionGainLocalUsed(Int_t subVersionGainLocalUsed)       { fSubVersionGainLocalUsed = subVersionGainLocalUsed;   }
  void SetVersionVdriftUsed(Int_t versionVdriftUsed)                   { fVersionVdriftUsed = versionVdriftUsed;   }
  void SetSubVersionVdriftUsed(Int_t subVersionVdriftUsed)             { fSubVersionVdriftUsed = subVersionVdriftUsed;   }
  
  void SetLow(Int_t low)                                               {fLow=low;};
  void SetHigh(Int_t high)                                             {fHigh=high;};
  void SetFillZero(Bool_t fillZero)                                    {fFillZero =  fillZero;};
  void SetNormalizeNbOfCluster(Bool_t normalizeNbOfCluster =  kTRUE)   {fNormalizeNbOfCluster = normalizeNbOfCluster;};
  void SetMaxCluster(Float_t maxCluster)                               {fMaxCluster =  maxCluster; }; 
  void SetNbMaxCluster(Short_t nbMaxCluster)                           {fNbMaxCluster =  nbMaxCluster; }; 
  void SetOfflineTracks()                                              {fOfflineTracks=kTRUE; fStandaloneTracks=kFALSE; };
  void SetStandaloneTracks()                                           {fStandaloneTracks=kTRUE; fOfflineTracks=kFALSE; };
  
  void SetCalDetGain(AliTRDCalDet * const calDetGain)                  {fCalDetGain = calDetGain;};  

  void SetMaxEvent(Int_t nbevents)                                     { fMaxEvent = nbevents; };
  void SetDebug(Int_t debug)                                           { fDebug = debug; };

 private:
  Bool_t SetVersionSubversion();
  Bool_t ParticleGood(int i) const;

  AliESDEvent  *fESD;                            //! ESD object
  AliESDfriend *fESDfriend;                      //! ESD friend
  const AliESDtrack *fkEsdTrack;                  //! ESD track
  AliESDfriendTrack *fFriendTrack;               //! ESD friend track
  TObject *fCalibObject;                         //! calibration objects attached to the ESD friend
  AliTRDtrackV1 *fTrdTrack;                      //! trdtrack
  AliTRDcluster *fCl;                            //! cluster
  
  TList       *fListHist;                        //! list of histograms

  AliTRDCalibraFillHisto *fTRDCalibraFillHisto;  //! calibration analyse object
  AliTRDCalibChamberStatus *fTRDChamberStatus;   //! calibration chamber status

  TH1I        *fNEvents;                         //! counter  
  TH1I        *fNEventsInput;                    //! counter  
  
  TH1F        *fNbTRDTrack;                      //! nb ESD tracks with TRD clusters
  TH1F        *fNbTRDTrackOffline;               //! nb ESD tracks with TRD clusters
  TH1F        *fNbTRDTrackStandalone;            //! nb ESD tracks with TRD clusters
  TH2F        *fNbTPCTRDtrack;                   //! nb TPC and TRD tracks when problems
  TH2F        *fNbGoodTracks;                    //! nb of good tracks
   
  TH1F        *fNbTimeBin;                       //! nb Time Bin
  TH1F        *fNbTimeBinOffline;                //! nb Time Bin offline
  TH1F        *fNbTimeBinStandalone;             //! nb Time Bin standalone
  TH1F        *fNbClusters;                      //! nb Clusters
  TH1F        *fNbClustersOffline;               //! nb Clusters offline
  TH1F        *fNbClustersStandalone;            //! nb Clusters standalone
  TH1F        *fNbTracklets;                     //! nb Tracklets
  TH1F        *fNbTrackletsOffline;              //! nb Tracklets offline
  TH1F        *fNbTrackletsStandalone;           //! nb Tracklets standalone
  
  TH2F        *fAbsoluteGain;                    //! Absolute Gain without AliESDfriend
  TH2I        *fCH2dSum;                         //! CH2d charge all
  TProfile2D  *fPH2dSum;                         //! PH2d PH all
  TH2I        *fCH2dSM;                          //! CH2d per SM
  TProfile2D  *fPH2dSM;                          //! PH2d per SM
  TH2I        *fCH2dTest;                        //! CH2d for test
  TProfile2D  *fPH2dTest;                        //! PH2d for test
  TH2S *fLinearVdriftTest;                       //! VdriftLinear for test

  Bool_t      fOnInstance;                       // On Instance
  Bool_t      fHisto2d;                          // histo
  Bool_t      fVector2d;                         // vector
  Bool_t      fVdriftLinear;                     // vdrift Linear

  Int_t       fNbTimeBins;                       // number of timebins 

  Short_t     fNz[3];                            // Nz mode 
  Short_t     fNrphi[3];                         // Nrphi mode
   
  TObjArray   *fSelectedTrigger;                 // Trigger class names accepted/rejected
  Bool_t      fRejected;                         // Reject the selected trigger class
  
  AliESDtrackCuts *fEsdTrackCuts;                // Quality cut on the AliESDtrack
  Bool_t      fRequirePrimaryVertex;             // Primary Vertex
  Bool_t      fVtxTPC;                           // Flag for use of TPC vertex
  Bool_t      fVtxSPD;                           // Flag for use of SPD vertex
  Int_t       fMinNbContributors;                // Min number of contributors
  Double_t    fRangePrimaryVertexZ;              // Were the primary vertex is
  Int_t       fMinNbTracks;                      // Min Nb Tracks
  Int_t       fMaxNbTracks;                      // Max Nb Tracks
  
  Int_t       fLow;                              // lower limit of nb of TRD clusters per tracklet
  Int_t       fHigh;                             // higher limit of nb of TRD clusters per tracklet
  Bool_t      fFillZero;                         // fill zero
  Bool_t      fNormalizeNbOfCluster;             // normalize with number of clusters (per default not) 
  Float_t     fRelativeScale;                    // relative scale for gas gain
  Float_t     fMaxCluster;                       // Maxcluster 
  Short_t     fNbMaxCluster;                     // Number of tb at the end
  Bool_t      fOfflineTracks;                    // Only Offline refitted tracks
  Bool_t      fStandaloneTracks;                 // Take only standalone tracks

  Int_t       fFirstRunGain;                     // FirstRunGainUsed 
  Int_t       fVersionGainUsed;                  // VersionGainUsed 
  Int_t       fSubVersionGainUsed;               // SubVersionGainUsed
  Int_t       fFirstRunGainLocal;                // FirstRunGainLocalUsed 
  Int_t       fVersionGainLocalUsed;             // VersionGainLocalUsed 
  Int_t       fSubVersionGainLocalUsed;          // SubVersionGainLocalUsed
  Int_t       fFirstRunVdrift;                   // FirstRunVdriftUsed 
  Int_t       fVersionVdriftUsed;                // VersionVdriftUsed 
  Int_t       fSubVersionVdriftUsed;             // SubVersionVdriftUsed
  Int_t       fFirstRunExB;                      // FirstRunExBUsed 
  Int_t       fVersionExBUsed;                   // VersionExBUsed 
  Int_t       fSubVersionExBUsed;                // SubVersionExBUsed

  AliTRDCalDet *fCalDetGain;                     // Calib object gain

  Int_t       fMaxEvent;                         // max events
  Int_t       fCounter;                          // max events
  Int_t       fDebug;                            // fDebug

  AliTRDCalibTask(const AliTRDCalibTask&); 
  AliTRDCalibTask& operator=(const AliTRDCalibTask&); 

  ClassDef(AliTRDCalibTask, 1); 
};

#endif


