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
  void SetExbAlt(Bool_t exbalt)                                     {fExbAlt = exbalt;};
  void SetNbTimeBins(Int_t nbTimeBins)                              {fNbTimeBins=nbTimeBins;};  
  void SetNumberBinCharge(Short_t nbBinCharge)                      {fNumberBinCharge=nbBinCharge;};  
  void SetRangeCharge(Float_t rangeCharge)                          {fRangeCharge=rangeCharge;};
  void SetVdBindx(Short_t vdBindx)                                  {fVdBindx=vdBindx;};  
  void SetVdBindy(Short_t vdBindy)                                  {fVdBindy=vdBindy;};    
  void SetVdRangedx(Double_t vdRangex)                               {fVdRangex=vdRangex;};  
  void SetVdRangedy(Double_t vdRangey)                               {fVdRangey=vdRangey;};    
  void SetDebugLevelTRDCalibraFillHisto(Short_t debugLevelTRDCalibraFillHisto) {fDebugLevelTRDCalibraFillHisto = debugLevelTRDCalibraFillHisto;};

  
  void SetNz(Short_t nz, Int_t i)                                      {fNz[i]=nz;};
  void SetNrphi(Short_t nrphi, Int_t i)                                {fNrphi[i]=nrphi;};
  
  void SetSelectTrigger(Bool_t selectTrigger)                          {fSelectTrigger = selectTrigger;};
  void AddSelectedTriggerClass(const char*name)                        {fSelectedTrigger->Add(new TObjString(name));};
  void SetReject(Bool_t rejected)                                      {fRejected = rejected;};
  
  void SetESDtrackCuts(AliESDtrackCuts * const esdtrackCuts)           {fEsdTrackCuts = esdtrackCuts;};
  void SetRequirePrimaryVertex(Bool_t requirePrimaryVertex)            {fRequirePrimaryVertex = requirePrimaryVertex;};
  void SetUseTPCVertex()                                               {fVtxTPC=kTRUE ; fVtxSPD=kFALSE;} 
  void SetUseSPDVertex()                                               {fVtxTPC=kFALSE; fVtxSPD=kTRUE ;} 
  void SetMinNbOfContributors(Int_t minNbOfContributors)               {fMinNbContributors = minNbOfContributors;};  
  void SetRangePrimaryVertexZ(Double_t rangePrimaryVertexZ)            {fRangePrimaryVertexZ = TMath::Abs(rangePrimaryVertexZ);}; 
  void SetRejectPileUpWithSPD(Bool_t rejectPileUpWithSPD)              {fRejectPileUpWithSPD = rejectPileUpWithSPD;};
  void SetRejectPileUpWithTOF(Bool_t rejectPileUpWithTOF)              {fRejectPileUpWithTOF = rejectPileUpWithTOF;};
  void SetRejectPileUpWithTOFOrITS(Bool_t rejectPileUpWithTOFOrITS)    {fRejectPileUpWithTOFOrITS = rejectPileUpWithTOFOrITS;};
  void SetMinNbTracks(Int_t minNbTracks)                               {fMinNbTracks = minNbTracks;};
  void SetMaxNbTracks(Int_t maxNbTracks)                               {fMaxNbTracks = maxNbTracks;};
  void SetCutWithVdriftCalib(Bool_t cutWithVdriftCalib)                {fCutWithVdriftCalib = cutWithVdriftCalib;};
  void SetMinNbTRDtracklets(Int_t minNbTRDtracklets)                   {fMinNbTRDtracklets = minNbTRDtracklets;};
  void SetMinTRDMometum(Double_t minTRDMomentum)                       {fMinTRDMomentum = minTRDMomentum;};
  void SetScaleGainWithTPCSignal(Bool_t scaleGainWithTPCSignal)        {fScaleGainWithTPCSignal = scaleGainWithTPCSignal;}; 

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

  Bool_t IsPHQon() const {return fPHQon;}
  void SetPHQon(const Bool_t kphq){ fPHQon = kphq; }

 private:
  Bool_t SetVersionSubversion();
  Bool_t ParticleGood(int i) const;

  AliESDEvent  *fESD;                            //! ESD object
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
  TH1I        *fNEventsTrigger;                  //! counter trigger   
  
  TH1F        *fNbTRDTrack;                      //! nb ESD tracks with TRD clusters
  TH1F        *fNbTRDTrackOffline;               //! nb ESD tracks with TRD clusters
  TH1F        *fNbTRDTrackStandalone;            //! nb ESD tracks with TRD clusters
  TH2F        *fNbTPCTRDtrack;                   //! nb TPC and TRD tracks when problems
  TH2F        *fNbGoodTracks;                    //! nb of good tracks versus centrality
  TH1F        *fNbGoodTracks1D;                  //! nb of good tracks
   
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
  TH2F        *fTOFbc;                           //! Check TOF branch crossing
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
  Bool_t      fExbAlt;                           // alternative exb calculation

  Short_t     fDebugLevelTRDCalibraFillHisto;    // Debug Level Fill Histo
  Int_t       fNbTimeBins;                       // number of timebins 
  Short_t     fNumberBinCharge;                  // Number of bins for the gain factor
  Float_t     fRangeCharge;                      // Range Charge
  Short_t     fVdBindx;                          // Nb of bin in vd histos x
  Short_t     fVdBindy;                          // Nb of bin in vd histos y
  Double_t    fVdRangex;                         // Range vd histos x
  Double_t    fVdRangey;                         // Range vd histos y

  Short_t     fNz[3];                            // Nz mode 
  Short_t     fNrphi[3];                         // Nrphi mode
   
  Bool_t      fSelectTrigger;                    // Select trigger
  TObjArray   *fSelectedTrigger;                 // Trigger class names accepted/rejected
  Bool_t      fRejected;                         // Reject the selected trigger class
  
  AliESDtrackCuts *fEsdTrackCuts;                // Quality cut on the AliESDtrack
  Bool_t      fRequirePrimaryVertex;             // Primary Vertex
  Bool_t      fVtxTPC;                           // Flag for use of TPC vertex
  Bool_t      fVtxSPD;                           // Flag for use of SPD vertex
  Int_t       fMinNbContributors;                // Min number of contributors
  Double_t    fRangePrimaryVertexZ;              // Were the primary vertex is
  Bool_t      fRejectPileUpWithSPD;              // Reject pile-up events with SPD
  Bool_t      fRejectPileUpWithTOF;              // Reject pile-up tracks with TOF
  Bool_t      fRejectPileUpWithTOFOrITS;         // Reject pile-up tracks with TOF or ITS
  Int_t       fMinNbTracks;                      // Min Nb Tracks
  Int_t       fMaxNbTracks;                      // Max Nb Tracks
  Bool_t      fCutWithVdriftCalib;               // CutWithVdriftCalib for the gain and PH
  Int_t       fMinNbTRDtracklets;                // Min number of TRD tracklets
  Float_t     fMinTRDMomentum;                   // Min TRD momentum  
  Bool_t      fScaleGainWithTPCSignal;           // Scale the TPC gain with the TPC signal

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

  Bool_t fPHQon;                                  //switch of phq

  AliTRDCalibTask(const AliTRDCalibTask&); 
  AliTRDCalibTask& operator=(const AliTRDCalibTask&); 

  ClassDef(AliTRDCalibTask, 8); 
};

#endif


