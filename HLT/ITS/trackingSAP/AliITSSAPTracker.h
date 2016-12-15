#ifndef ALIITSSAPTRACKER_H
#define ALIITSSAPTRACKER_H

#include <TH1.h>
#include <TH2.h>
#include <TBits.h>
#include <TStopwatch.h>
#include <TObjArray.h>
#include <algorithm>
#include <vector>
#include "AliExternalTrackParam.h"
#include "AliESDVertex.h"

//------- compilation options, comment out all for best performance ------
//#define _TIMING_                            // print timing info, use only in offline mode
//#define _CONTROLH_                          // fill control histos, use only in offline mode
//#define _DEBUG_                             // print debug info, use only in offline mode
//------------------------------------------------------------------------

class AliITSRecPoint;
class AliESDVertex;
class AliITSSAPLayer;
class AliHLTITSTrackPoint;

class AliITSSAPTracker : public TObject
{
 public:
  enum {kALrSPD1,kALrSPD2, kALrSDD1,kALrSDD2, kALrSSD1,kALrSSD2,kNLrActive};
  enum {kLrBeamPime, kLrSPD1,kLrSPD2, kLrShield1, kLrSDD1,kLrSDD2, kLrShield2, kLrSSD1,kLrSSD2,
	kMaxLrITS,kNLrPassive=kMaxLrITS-kNLrActive};
  enum {kInvalidBit=BIT(14),kVtUsedBit=BIT(15)};
  //
  struct SPDtracklet {
  SPDtracklet() : id1(0),id2(0),label(0),dphi(0),dtht(0),chi2(0) {}
    int id1;
    int id2;
    int label;
    float dphi;
    float dtht;
    float chi2;
  };
  typedef struct SPDtracklet SPDtracklet_t;
  //
  struct ITStrack {
  ITStrack() : paramOut(),paramInw(),chi2(0),ncl(0),nmiss(0),label(0),trackletID(0) {}
    AliExternalTrackParam paramOut;
    AliExternalTrackParam paramInw;
    float chi2;
    short ncl;
    short nmiss;
    int clID[6];
    int label;
    int trackletID;
  };
  typedef struct ITStrack ITStrack_t;
  //
  AliITSSAPTracker();
  virtual ~AliITSSAPTracker();
  //
  void ProcessEvent();
  void Init();
  void Clear(Option_t *opt="");
  void ClearTracklets();
  void ClearTracks()                               {fTracks.clear();}
  //
  void SetSPDVertex(const AliESDVertex* v)         {fSPDVertex = v;}
  const AliESDVertex* GetSPDVertex()  const        {return fSPDVertex;}
  void AddCluster(AliITSRecPoint* cl);
  void SetBz(float v)                              {fBz = v;}
  void SetMaxRSPDVtx(float v=1.5)                  {fMaxRSPDVtx = v>0.5 ? v:0.5;}
  //
  // methods for trackleting ---------------->>>
  Bool_t FindTracklets();
  Int_t  AssociateClusterOfL2(int icl2);
  Bool_t IsBlacklisted(int id1,int id2);
  void   Blacklist(int id1,int id2);
  //
  void SetPhiShift(float v=0.0045)                  {fPhiShift = v;}
  void SetSigThetaTracklet(float v=0.025)           {fSigThetaTracklet = v;}
  void SetSigPhiTracklet(float v=0.08)              {fSigPhiTracklet = v;}
  void SetChi2CutTracklet(float v=1.5)              {fChi2CutTracklet = v;}
  void SetMaxTrackletsToRunTracking(int n=9999)     {fMaxTrackletsToRunTracking = n;}
  Int_t GetMaxTrackletsToRunTracking() const        {return fMaxTrackletsToRunTracking;}
  //
  Double_t GetClSystYErr2(Int_t lr)    const        {return fgkClSystYErr2[lr];}
  Double_t GetClSystZErr2(Int_t lr)    const        {return fgkClSystZErr2[lr];}
  //
  int  GetNTracklets()                 const        {return (int)fTracklets.size();}
  int  GetNTracks()                    const        {return fNTracks;}
  void PrintTracklets()                const;
  void PrintTracklet(Int_t itr)        const;
  const AliITSSAPTracker::SPDtracklet_t& GetTracklet(int i) const {return fTracklets[i];}
  // methods for trackleting ----------------<<<
  //
  // methods for track reconstruction ------->>>
  Float_t GetMinPt()                   const        {return fMinPt;}
  void    SetMinPt(Float_t v=0.3)                   {fMinPt = v<0.2 ? 0.2 : v;}
  void    CalcAuxTracking();
  Bool_t  CreateTrack(AliITSSAPTracker::ITStrack_t& track, AliITSSAPTracker::SPDtracklet_t& trlet);
  void    Tracklets2Tracks();
  AliITSSAPLayer* GetLayer(int i)         const        {return (AliITSSAPLayer*)fLayers[i];}
  Int_t   GetActiveLayerID(int i)      const        {return fgkLr2Active[i];}
  Float_t GetChi2TotCut(int ncl)       const;
  Bool_t  CrossPassiveLayer(AliExternalTrackParam& track, Int_t lrID);
  Bool_t  FollowToLayer(AliITSSAPTracker::ITStrack_t& track, Int_t lrID);
  Double_t GetXatLabRLin(AliExternalTrackParam& track, double r);
  void    CookLabel(AliITSSAPTracker::ITStrack_t& track);
  void    CookLabel(AliITSSAPTracker::SPDtracklet_t& tracklet);
  void    PrintTrack(const AliITSSAPTracker::ITStrack_t& track) const;
  Bool_t  GetSkipLayer(int lr)                    const  {return fSkipLayer[lr];}
  void    SetSkipLayer(int lr, Bool_t v=kTRUE)           {fSkipLayer[lr] = v;}
  Bool_t  IsAcceptableTrack(const AliITSSAPTracker::ITStrack_t& track) const;
  void    PrintTracks()                const;
  Int_t   GetTrackletMCTruth(AliITSSAPTracker::SPDtracklet_t& trlet) const;
  void    RefitInward();
  Bool_t  RefitInward(int itr);
  void    SetMaxMissedLayers(int n=0)  { fMaxMissedLayers = n;}
  Int_t   GetMaxMissedLayers()    const  { return fMaxMissedLayers;}
  AliITSSAPTracker::ITStrack_t& GetTrack(int i) const {return (ITStrack_t &)fTracks[i];}
  // methods for track reconstruction -------<<<
  //
  // methods for vertex reconstruction ------>>>
  Bool_t  FitTrackVertex();
  AliESDVertex& GetTrackVertex()      const       {return (AliESDVertex&)fTrackVertex;}
  void      SetMaxVtxIter(Int_t n=5)              {fMaxVtxIter = n;}
  Int_t     GetMaxVtxIter()           const       {return fMaxVtxIter;}
  void      SetStopScaleChange(float v=0.8)       {fStopScaleChange = v;}
  Float_t   GetStopScaleChange()           const  {return fStopScaleChange;}
  // methods for vertex reconstruction ------<<<
  //

  // track point calculation

  Int_t GetTrackPoint( Int_t iLayer, Int_t clusterIndex, AliHLTITSTrackPoint& p ) const ;

 protected:
  //
  AliITSSAPLayer* fLayers[kNLrActive];
  Bool_t    fSkipLayer[kNLrActive];                     //! skip layer
  Int_t     fNClusters[kNLrActive];                     //! number of clusters per event
  //
  // for SPD trackleting ----------------- >>>
  std::vector<bool> fSPD2Discard;                       //! status of SPD2 clusters during trackleting
  std::vector<SPDtracklet_t> fTracklets;                //! found tracklets
  std::vector<int> fSPD1Tracklet;                       //! id+1 of traclet using this SPD1 cluster
  TBits*   fBlacklist;                            //! blacklisted combinations
  Float_t  fPhiShift;                             //! Phi shift reference value (at 0.5 T)
  Float_t  fSigThetaTracklet;                     //! sigTheta for tracklets
  Float_t  fSigPhiTracklet;                       //! sigPhi for tracklets
  Float_t  fChi2CutTracklet;                      //! cut on tracklet total chi2
  Float_t  fPhiShiftSc;                           //! dPhi offset to account for bending
  Float_t  fDThetaTrackletSc;                     //! max dTheta for tracklets with scaling from chi2 cut
  Float_t  fDPhiTrackletSc;                       //! max dPhi for tracklets with scaling from chi2 cut
  Float_t  fBz;                                   //! Bz field in ITS
  Float_t  fMaxRSPDVtx;                           //! max allowed R of SPD vertex
  //
  // auxilary precomputed stuff
  Float_t  fDPhiTol;                              //! tolerance on phi, accounting for bending
  Float_t  fDThSig2Inv;                           //! inverse of sigma dTheta
  Float_t  fDPhSig2Inv;                           //! inverse of sigma dPhi
  // for SPD trackleting ----------------- <<<
  //
  // for track reconstruction ------------ >>>
  Float_t  fMinPt;                                //! user pt cutoff
  Float_t  fCurvMax;                              //! max curvature to reconstruct
  Float_t  fZSPD2CutMin;                          //! min Z of tracklet SPD2 to consider tracking
  Float_t  fZSPD2CutMax;                          //! max Z of tracklet SPD2 to consider tracking
  Float_t  fMaxChi2Tr2Cl;                         //! cut on cluster-to-track chi2
  Float_t  fAddErr2YspdVtx;                       //! additional error to Y of the SPD vertex in track fit
  Float_t  fAddErr2ZspdVtx;                       //! additional error to Z of the SPD vertex in track fit
  Float_t  fChi2TotCut[kNLrActive];               //! cut on total chi2 depending on track length
  //
  Float_t  fNSigma2[kNLrActive];                  //! N^2 sigma margin for cluster search
  Float_t  fYToler2[kNLrActive];                  //! Y additional margin^2 for cluster search
  Float_t  fZToler2[kNLrActive];                  //! Z additional margin^2 for cluster search
  Float_t  fMaxDRPhi;                             //! maximal search window in rphi
  Float_t  fMaxDZ;                                //! maximal search window in z

  Float_t  fMSDist[kNLrActive];                   //! shift due to the MS for 1 GeV particle
  Float_t  fMSPhi[kNLrActive];                    //! dphi due to the MS for 1 GeV particle
  Float_t  fTolPhiCrude[kNLrActive];              //! tolerance in dphi for particle of unknown momentum
  Float_t  fTolZCrude[kNLrActive];                //! tolerance in Z for particle of unknown momentum
  Float_t  fMissChi2Penalty;                      //! penalize missed clusters
  Int_t    fMaxMissedLayers;                      //! allow to miss at most this number of layers
  Int_t    fNTracks;                              //! n found tracks
  Int_t    fMaxTrackletsToRunTracking;            //! skip tracking if too many SPD tracklets
  //
  Int_t    fMaxVtxIter;                           //! max number of vertexing iterations
  Float_t  fStopScaleChange;                      //! stop vertexing if sigma scaling change is above this
  //
  std::vector<ITStrack_t> fTracks;                //! found tracks container
  AliESDVertex fTrackVertex;                      //! fitted track vertex
  Bool_t    fFitVertex;                           //! fit vertex with tracks
  // for track reconstruction ------------ <<<
  //		      
		      
  const AliESDVertex* fSPDVertex;                  //! external vertex

  static const Float_t fgkRhoLITS[kMaxLrITS];      // <rho*L> for each material layer
  static const Float_t fgkX2X0ITS[kMaxLrITS];      // <x/x0> for each material layer
  static const Float_t fgkRLayITS[kMaxLrITS];     // radii of material layers
  static const Float_t fgkRSpanITS[kMaxLrITS];    // half R span of the material layer
  static const Float_t fgkZSpanITS[kMaxLrITS];    // half Z span of the material layer
  static const Int_t   fgkPassivLrITS[kNLrPassive];  // list of passive layer enums
  static const Int_t   fgkActiveLrITS[kNLrActive]; // list of active layer enums
  static const Double_t fgkClSystYErr2[kNLrActive]; // syst error^2 for Y direction
  static const Double_t fgkClSystZErr2[kNLrActive]; // syst error^2 for Y direction

  static const Int_t   fgkLr2Active[kMaxLrITS]; // conversion from LrID to ActiveLr ID
  static const Int_t   fgkLrDefBins[kNLrActive][2]; // default binning for cluster navigator
  static const Int_t   fgkDummyLabel;               // dummy MC label
  static const Float_t fgkDefMass;                  // default mass for tracking
  //
#ifdef _TIMING_
 public:
  enum {kSWTotal,kSWTracklets,kSWTracks,kSWVertex,kNSW};
  void PrintTiming();
  const TStopwatch& GetStopwatch(int i)     const {return fSW[i];}
  const char*       GetStopwatchName(int i) const {return fgkSWNames[i];}
 protected:
  static const char* fgkSWNames[kNSW];
  TStopwatch fSW[kNSW];
#endif
  //
#ifdef _CONTROLH_
 protected:
  TObjArray fArrHisto;
  TH2F *fHTrackletMC,*fHTrackletAll,*fHTrackletFake,*fHTrackMC,*fHTrackAll,*fHTrackFake,*fHVtxDiffXY
    ,*fHVtxDiffXMlt,*fHVtxDiffYMlt,*fHVtxDiffZMlt
    ,*fHVtxPullXMlt,*fHVtxPullYMlt,*fHVtxPullZMlt
    ,*fHVtxMCSPDDiffXY
    ,*fHVtxMCSPDDiffXMlt,*fHVtxMCSPDDiffYMlt,*fHVtxMCSPDDiffZMlt
    ,*fHVtxMCSPDPullXY
    ,*fHVtxMCSPDPullXMlt,*fHVtxMCSPDPullYMlt,*fHVtxMCSPDPullZMlt
    ,*fHChi2NDFvsPT,*fHChi2vsNC;
  TH1F *fHVtxMltRef,*fHVtxOKMlt,*fHVtxDiffZ,*fHVtxMCSPDDiffZ;
  //
  void FillRecoStat();
  void FillTrackingControlHistos(int lrID,int lbl,const AliExternalTrackParam* bestTr,
				 const double cpar[2],const double ccov[3],const AliITSRecPoint* bestCl);
  void BookHistos();
  Double_t* DefLogAx(double xMn,double xMx, int nbin);
 public:
  void SaveHistos(const char* outFName="XXXITSTrackerControlH.root");
  TObjArray* GetHistos() const {return (TObjArray*)&fArrHisto;}
  enum {kHResidY,kHPullY,kHResidZ,kHPullZ,kHChi2Cl};
#endif

 private:
  AliITSSAPTracker(const AliITSSAPTracker&);
  AliITSSAPTracker& operator=(const AliITSSAPTracker&);
  ClassDef(AliITSSAPTracker,0) // ITS SA primaries tracker/vertexer
};


//______________________________________________
inline Bool_t AliITSSAPTracker::IsBlacklisted(int id1,int id2)
{
  // check if this combination is blacklisted
  return fBlacklist->TestBitNumber(UInt_t(id1*fNClusters[0])+id2);
}

//______________________________________________
inline void AliITSSAPTracker::Blacklist(int id1,int id2)
{
  // blacklist this combination
  return fBlacklist->SetBitNumber(UInt_t(id1*fNClusters[0])+id2);
}

//______________________________________________
inline Float_t AliITSSAPTracker::GetChi2TotCut(int ncl) const
{
  // return chi2 cut for given number of clusters. Min ncl=3
  return fChi2TotCut[ncl-2];
}

#endif
