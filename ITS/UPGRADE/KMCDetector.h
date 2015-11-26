#ifndef DETECTORK_H
#define DETECTORK_H

#include <TNamed.h>
#include <TList.h>
#include <TClonesArray.h>
#include <TArrayD.h>
#include "AliExternalTrackParam.h"

//-------------------------------------------------------------------------
// Current support and development: Ruben Shahoyan (Ruben.Shahoyan@cern.ch) 
//-------------------------------------------------------------------------

class KMCLayer;
class KMCCluster;
class TH2F;
class TH1F;
class TGraph;
class TArrayI;

//////////////////////////////////////////////////////////////////////////////////////////////
//--------------------------------------------------------------------------------------------
class KMCProbe : public AliExternalTrackParam {
 public:
  enum {kBitKilled=BIT(14)};
  enum {kNDOF=5};
  enum {kY2=0,kZ2=2,kSnp2=5,kTgl2=9,kPtI2=14};
  enum {kY,kZ,kSnp,kTgl,kPtI};
  //
 KMCProbe() : fMass(0.14),fChi2(0),fHits(0),fFakes(0),fNHits(0),fNHitsITS(0),fNHitsITSFake(0),fInnLrCheck(fgNITSLayers) {for (int i=kMaxITSLr;i--;) fClID[i] =-2;
}
  KMCProbe(KMCProbe& src);
  KMCProbe& operator=(const KMCProbe& src);
  virtual ~KMCProbe() {}
  virtual   void   Print(Option_t* option = "") const;
  virtual   void   Reset() {fMass=0.14; fChi2=0; fHits=fFakes=0;  for (int i=kMaxITSLr;i--;) fClID[i]=-2; AliExternalTrackParam::Reset();}
  virtual   Bool_t IsSortable()                     const {return kTRUE;}
  virtual   Int_t  Compare(const TObject* obj)      const;
  void      ResetCovMat();
  //
  void      Kill(Bool_t v=kTRUE)                          {SetBit(kBitKilled,v);}
  Bool_t    IsKilled()                              const {return TestBit(kBitKilled);}
  //
  Bool_t    CorrectForMeanMaterial(const KMCLayer* lr, Bool_t inward=kTRUE);
  Bool_t    GetXatLabR(Double_t r,Double_t &x, Double_t bz, Int_t dir=0) const;
  Bool_t    PropagateToR(double r, double b, int dir=0);
  
  //
  void      SetMass(double m=0.14)                      {fMass = m;}
  Double_t  GetMass()                             const {return fMass;}
  Double_t  GetChi2()                             const {return fChi2;}
  void      SetChi2(double chi2)                        {fChi2 = chi2;}
  void      AddChi2(double chi2)                        {fChi2 += chi2;}
  void      SetInnerLrChecked(Int_t n)                  {if (n<fgNITSLayers) fInnLrCheck = n;}
  UShort_t  GetNITSHits()                         const {return fNHitsITS;}
  UShort_t  GetNFakeITSHits()                     const {return fNHitsITSFake;}
  UShort_t  GetNHits()                            const {return fNHits;}
  UShort_t  GetInnerLayerChecked()                const {return fInnLrCheck;}
  UInt_t&   GetHitsPatt()                               {return fHits;}
  UInt_t&   GetFakesPatt()                              {return fFakes;}
  Double_t  GetNormChi2(Bool_t penalize=kFALSE)   const;
  void      AddHit(Int_t lr, double chi2, Int_t clID=-1);
  void      ResetHit(Int_t lr);
  Bool_t    IsHit(Int_t lr)                       const {return (lr<fgNITSLayers) ? IsWBit(fHits,lr)  : kFALSE;}
  Bool_t    IsHitFake(Int_t lr)                   const {return (lr<fgNITSLayers) ? IsWBit(fFakes,lr) : kFALSE;}
  Bool_t    PropagateToCluster(KMCCluster* cl, Double_t b);
  // protected: 
  static void   SetWBit(UInt_t &patt,UInt_t bit)               {patt |= 0x1<<bit;}
  static void   ResetWBit(UInt_t &patt,UInt_t bit)             {patt &= ~(0x1<<bit);}
  static Bool_t IsWBit(const UInt_t &patt,const UInt_t bit)    {return patt&(0x1<<bit);}
  static void   SetNITSLayers(Int_t n)                         {fgNITSLayers = n;}
  static int    GetNITSLayers()                                {return fgNITSLayers;}
  //
  static Double_t GetMissingHitPenalty()                        {return fgMissingHitPenalty;}
  static void     SetMissingHitPenalty(double p=2.)             {fgMissingHitPenalty = p;}
  //
 public:
  enum {kMaxITSLr=12};
  Float_t fMass;   // mass
  Float_t fChi2;   // total chi2
  UInt_t  fHits;   // pattern on hits (max 32!)
  UInt_t  fFakes;  // pattern of fakes among hits
  UShort_t fNHits;    // total hits
  UShort_t fNHitsITS; // total ITS hits
  UShort_t fNHitsITSFake; // number of fake ITS hits
  UShort_t fInnLrCheck;   // lowest active layer where update was checked
  Int_t    fClID[kMaxITSLr];
  //
  static Int_t    fgNITSLayers;
  static Double_t fgMissingHitPenalty;
  ClassDef(KMCProbe,2);  
};

//_______________________________________
inline Double_t KMCProbe::GetNormChi2(Bool_t penalize) const
{
  // normalized chi2, penilized for missing hits
  //  if (fNHitsITS<3) return 0;
  double chi2 = fChi2;
  if (penalize) {
    int nMiss = fgNITSLayers - fNHitsITS - fInnLrCheck;
    chi2 = fChi2 + fgMissingHitPenalty*nMiss;
  }
  return chi2/( (fNHitsITS<<1)-kNDOF);
}

//_______________________________________
inline void KMCProbe::AddHit(Int_t lr, double chi2, Int_t clID) {
  // note: lr is active layer ID
  if (lr<0) return;
  fNHits++;
  if (lr<fgNITSLayers) {
    fChi2 += chi2;
    SetWBit(fHits,lr); 
    fNHitsITS++;
    if (clID>-1) {
      SetWBit(fFakes,lr);
      fNHitsITSFake++;
    }
    fClID[lr] = clID;
    //else ResetWBit(fFakes,lr);
  }
}

//_______________________________________
inline void KMCProbe::ResetHit(Int_t lr) {
  // note: lr is active layer ID
  if (lr>=fgNITSLayers) return; 
  if (IsWBit(fHits,lr))  {fNHitsITS--;     ResetWBit(fHits,lr);}
  if (IsWBit(fFakes,lr)) {fNHitsITSFake--; ResetWBit(fFakes,lr);}
}

//////////////////////////////////////////////////////////////////////////////////////////////
//--------------------------------------------------------------------------------------------
class KMCCluster : public TObject {
 public:
  //
  enum {kBitKilled=BIT(14)};
 KMCCluster(Float_t y=0, Float_t z=0, Float_t x=0, Float_t phi=0) : fY(y),fZ(z),fX(x),fPhi(phi) {}
  KMCCluster(KMCCluster &src);
  KMCCluster& operator=(const KMCCluster& src);
  virtual ~KMCCluster() {}
  void Reset()               {Clear();}
  //
  Double_t GetY()    const   {return fY;}
  Double_t GetX()    const   {return fX;}
  Double_t GetZ()    const   {return fZ;}
  Double_t GetPhi()  const   {return fPhi;}
  //
  void    Kill(Bool_t v=kTRUE)          {SetBit(kBitKilled,v);}
  Bool_t  IsKilled()              const {return TestBit(kBitKilled);}
  Float_t fY; 
  Float_t fZ; 
  Float_t fX;
  Float_t fPhi;
  void Set(Float_t y, Float_t z, Float_t x, Float_t phi) {fY=y; fZ=z; fX=x; fPhi=phi; ResetBit(kBitKilled);}
  virtual void Print(Option_t * = 0) const;
  ClassDef(KMCCluster,1);
};

//////////////////////////////////////////////////////////////////////////////////////////////
//--------------------------------------------------------------------------------------------
class KMCLayer : public TNamed {
public:
  enum {kTypeNA=-1,kITS,kTPC,kBitVertex=BIT(15)};
  KMCLayer(char *name);
  Float_t GetRadius()   const {return fR;}
  Float_t GetRadL()     const {return fx2X0;}
  Float_t GetXTimesRho() const {return fXRho;}
  Float_t GetPhiRes()   const {return fPhiRes;}
  Float_t GetZRes()     const {return fZRes;}
  Float_t GetLayerEff() const {return fEff;}
  Int_t   GetActiveID() const {return fActiveID;}
  virtual void  Print(Option_t* option = "") const;
  //
  Bool_t IsDead()       const {return fIsDead;}
  Bool_t IsITS()        const {return fType==kITS;}
  Bool_t IsTPC()        const {return fType==kTPC;}
  Bool_t IsVertex()     const {return TestBit(kBitVertex);}
  //
  Int_t    AddBgCluster(double y,double z,double x,double phi) {int n=GetNBgClusters(); new (fClBg[n]) KMCCluster(y,z,x,phi); return ++n;}
  KMCCluster* GetBgCluster(Int_t i) {return (KMCCluster*)fClBg[i];}
  KMCCluster* GetMCCluster()        {return (KMCCluster*)&fClMC;}
  //
  KMCProbe*       GetAnProbe()    const {return (KMCProbe*)&fTrCorr;}
  Int_t           GetNMCTracks()  const {return fTrMC.GetEntries();}
  KMCProbe*       AddMCTrack(KMCProbe* src=0);
  KMCProbe*       GetMCTrack(Int_t it) const {return (KMCProbe*)fTrMC[it];}
  KMCProbe*       GetWinnerMCTrack()  {if (!fTrMC.IsSorted()) fTrMC.Sort(); return fTrMC.GetEntries() ? (KMCProbe*)fTrMC[0]:0;}
  TClonesArray*         GetMCTracks()        const {return (TClonesArray*)&fTrMC;}
  void Reset();
  void ResetMC() { fClBg.Clear(); fTrMC.Clear();}
  void ResetBgClusters() { fClBg.Clear(); }
  void ResetMCTracks()   { fTrMC.Clear(); }
  Int_t GetNBgClusters() const { return fClBg.GetEntries();}
  static Double_t GetDefEff()   {return fgDefEff;}
  static void     SetDefEff(double eff=1) {fgDefEff = eff>1. ? 1.: (eff<0? 0:eff);}
  //
  //
  Float_t fR; 
  Float_t fx2X0;
  Float_t fXRho;    // x*density
  Float_t fPhiRes; 
  Float_t fZRes;   
  Float_t fEff;
  Bool_t  fIsDead;
  Int_t   fType;   // its, tpc etc
  Int_t   fActiveID;   // active layer id
  Float_t fSig2EstD;
  Float_t fSig2EstZ;
  //
  KMCCluster   fClCorr;     // ideal cluster
  KMCCluster   fClMC;       // MC cluster (from MS scattered track)
  TClonesArray fClBg; // bg clusters for MC
  //
  KMCProbe     fTrCorr;   // ideal track
  TClonesArray fTrMC;      // MC tracks
  //
  static Double_t fgDefEff;
  ClassDef(KMCLayer,1);
};

//////////////////////////////////////////////////////////////////////////////////////////////
//--------------------------------------------------------------------------------------------
class KMCDetector : public TNamed {
 public:
  enum {kUtilHisto=BIT(14)};
  KMCDetector();
  KMCDetector(char *name,char *title);
  virtual ~KMCDetector();

  void AddLayer(char *name, Float_t radius, Float_t radL, Float_t xrho=0, Float_t phiRes=999999, Float_t zRes=999999, Float_t eff=-1);
  Int_t GetLayerID(Int_t actID) const;
  void KillLayer(char *name);
  void SetRadius(char *name, Float_t radius);
  void SetRadiationLength(char *name, Float_t radL);
  void SetResolution(char *name, Float_t phiRes=999999, Float_t zRes=999999);
  void SetLayerEfficiency(char *name, Float_t eff=1.0);
  void RemoveLayer(char *name);

  Float_t GetRadius(char *name);
  Float_t GetRadiationLength(char *name);
  Float_t GetResolution(char *name, Int_t axis=0);
  Float_t GetLayerEfficiency(char *name);

  void PrintLayout(); 
  void PlotLayout(Int_t plotDead = kTRUE);
  
  void MakeAliceAllNew(Bool_t flagTPC =0,  Bool_t flagMon=1,int setVer=0);
  void MakeAliceCurrent(Bool_t flagTPC =0, Int_t AlignResiduals = 0);
  void AddTPC(Float_t phiResMean, Float_t zResMean, Int_t skip=1);
  void RemoveTPC();

  void SetBField(Float_t bfield) {fBFieldG = bfield*10; }
  Float_t GetBField() const {return fBFieldG/10; }
  void SetLhcUPCscale(Float_t lhcUPCscale) {fLhcUPCscale = lhcUPCscale; }
  Float_t GetLhcUPCscale() const { return fLhcUPCscale; }
  void SetParticleMass(Float_t particleMass) {fParticleMass = particleMass; }
  Float_t GetParticleMass() const { return fParticleMass; }
  void SetIntegrationTime(Float_t integrationTime) {fIntegrationTime = integrationTime; }
  Float_t GetIntegrationTime() const { return fIntegrationTime; }
  void SetMaxRadiusOfSlowDetectors(Float_t maxRadiusSlowDet) {fMaxRadiusSlowDet =  maxRadiusSlowDet; }
  Float_t GetMaxRadiusOfSlowDetectors() const { return fMaxRadiusSlowDet; }
  void SetAvgRapidity(Float_t avgRapidity) {fAvgRapidity = avgRapidity; }
  Float_t GetAvgRapidity() const { return fAvgRapidity; }
  void SetConfidenceLevel(Float_t confLevel) {fConfLevel = confLevel; }
  Float_t GetConfidenceLevel() const { return fConfLevel; }
  void SetAtLeastCorr(Int_t atLeastCorr ) {fAtLeastCorr = atLeastCorr; }
  Float_t GetAtLeastCorr() const { return fAtLeastCorr; }
  void SetAtLeastFake(Int_t atLeastFake ) {fAtLeastFake = atLeastFake; }
  Float_t GetAtLeastFake() const { return fAtLeastFake; }

  void SetdNdEtaCent(Int_t dNdEtaCent ) {fdNdEtaCent = dNdEtaCent; }
  Float_t GetdNdEtaCent() const { return fdNdEtaCent; }

  Int_t GetNLayers()          const {return fLayers.GetEntries(); }
  Int_t GetNActiveLayers()    const {return fNActiveLayers; }
  Int_t GetNActiveITSLayers() const {return fNActiveITSLayers; }

  Bool_t SolveSingleTrackViaKalman(Double_t mass, Double_t pt, Double_t eta);
  Bool_t SolveSingleTrackViaKalmanMC(int offset=6);
  Bool_t SolveSingleTrack(Double_t mass, Double_t pt, Double_t eta, TObjArray* sumArr=0, int nMC=10000,int offset=6);
  KMCProbe* KalmanSmooth(int actLr, int actMin,int actMax) const;
  KMCProbe* KalmanSmoothFull(int actLr, int actMin,int actMax) const; //TBD
  void   EliminateUnrelated();
  //
  KMCProbe* PrepareKalmanTrack(double pt, double lambda, double mass, int charge, double phi=0,double x=0,double y=0,double z=0);
  Bool_t TransportKalmanTrackWithMS(KMCProbe *probTr, int maxLr) const;
  void   ApplyMS(KMCProbe* trc,  double x2x0) const;
  Bool_t PropagateToLayer(KMCProbe* trc, KMCLayer* lr, int dir) const;
  Bool_t UpdateTrack(KMCProbe* trc, KMCLayer* lr, KMCCluster* cl, Bool_t goToCluster=kTRUE) const;

  // Helper functions
  Double_t ThetaMCS                 ( Double_t mass, Double_t RadLength, Double_t momentum ) const;
  Double_t ProbGoodHit              ( Double_t radius, Double_t searchRadiusRPhi, Double_t searchRadiusZ )   ; 
  Double_t ProbGoodChiSqHit         ( Double_t radius, Double_t searchRadiusRPhi, Double_t searchRadiusZ )   ; 
  Double_t ProbGoodChiSqPlusConfHit ( Double_t radius, Double_t leff, Double_t searchRadiusRPhi, Double_t searchRadiusZ )   ; 
  Double_t ProbNullChiSqPlusConfHit ( Double_t radius, Double_t leff, Double_t searchRadiusRPhi, Double_t searchRadiusZ )   ; 
 
  // Howard W. hit distribution and convolution integral
  Double_t Dist              ( Double_t Z, Double_t radius ) ;  
  Double_t HitDensity        ( Double_t radius )   ;
  Double_t UpcHitDensity     ( Double_t radius )   ;
  Double_t IntegratedHitDensity  ( Double_t multiplicity, Double_t radius )   ;
  Double_t OneEventHitDensity    ( Double_t multiplicity, Double_t radius ) const   ;
  Double_t D0IntegratedEfficiency( Double_t pt, Double_t corrEfficiency[][20] ) const ;
  
  TGraph* GetGraphMomentumResolution(Int_t color, Int_t linewidth=1);
  TGraph* GetGraphPointingResolution(Int_t axis,Int_t color, Int_t linewidth=1);
  TGraph* GetGraphPointingResolutionTeleEqu(Int_t axis,Int_t color, Int_t linewidth=1);

  TGraph* GetGraphImpactParam(Int_t mode, Int_t axis, Int_t color, Int_t linewidth=1);

  TGraph* GetGraphRecoEfficiency(Int_t particle, Int_t color, Int_t linewidth=1); 
  TGraph* GetGraphRecoFakes(Int_t particle,Int_t color, Int_t linewidth);

  TGraph* GetGraph(Int_t number, Int_t color, Int_t linewidth=1);

  void MakeStandardPlots(Bool_t add =0, Int_t color=1, Int_t linewidth=1,Bool_t onlyPionEff=0);

  // method to extend AliExternalTrackParam functionality
  Bool_t IsZero(double val, double tol=1e-9) const {return TMath::Abs(val)<tol;}
  TList *GetLayers()                   const {return (TList*)&fLayers;}
  KMCLayer* GetLayer(Int_t i)          const {return (KMCLayer*) fLayers.At(i);}
  KMCLayer* GetActiveLayer(Int_t actID)    const {int pid=GetLayerID(actID); return pid<0 ? 0:GetLayer(pid);}
  KMCLayer* GetLayer(const char* name) const {return (KMCLayer*) fLayers.FindObject(name);}
  KMCProbe* GetProbeTrack()       const {return (KMCProbe*)&fProbe;}
  void   ClassifyLayers();
  void   ResetMCTracks(int maxLr=-1);                      
  //
  Bool_t   GetUseBackground()               const {return fUseBackground;}
  void     SetUseBackground(Bool_t v=kTRUE)       {fUseBackground = v;}
  void     CheckTrackProlongations(KMCProbe *probe, KMCLayer* lr, KMCLayer* lrP);
  void     ResetSearchLimits() {fBgYMin=fBgZMin=1e6; fBgYMax=fBgZMax=-1e6; fNBgLimits=0;}
  void     UpdateSearchLimits(KMCProbe* probe, KMCLayer* lr);
  Int_t    GenBgClusters(KMCLayer* lr);
  Bool_t   NeedToKill(KMCProbe* probe) const;
  Double_t PropagateBack(KMCProbe* trc);
  //
  // definition of reconstructable track
  void     RequireMaxChi2Cl(double cut=25.)           {fMaxChi2Cl = cut>0 ? cut:9; fMaxChi2ClSQ = TMath::Sqrt(fMaxChi2Cl);}
  void     RequireMinITSHits(Int_t n=4)               {fMinITSHits = n;}
  void     RequireMaxNormChi2NDF(double cut=5.)       {fMaxNormChi2NDF = cut>0 ? cut:9;}
  void     RequirePattern(UInt_t *patt, int groups);
  //
  Double_t GetMaxChi2Cl()                      const {return fMaxChi2Cl;}
  Double_t GetMaxNormChi2NDFusterKMC()              const {return fMaxNormChi2NDF;}
  Int_t    GetMinITSHits()                     const {return fMinITSHits;}
  //
  Double_t GetUpdCalls()                       const {return fUpdCalls;}
  void     InitMCWatchHistos();
  TH2F*    GetHMCLrResidRPhi()                 const {return fHMCLrResidRPhi;}
  TH2F*    GetHMCLrResidZ()                    const {return fHMCLrResidZ;}
  TH2F*    GetHMCLrChi2()                      const {return fHMCLrChi2;}
  //
  void     PrintITS(Option_t* opt="") const {for (int i=0;i<=fLastActiveITSLayer;i++) if (!GetLayer(i)->IsDead()) GetLayer(i)->Print(opt);}
  static void SetVtxConstraint(double d=-1, double z=-1) {fgVtxConstraint[0]=d; fgVtxConstraint[1]=z;}
  //
  void CalcHardSearchLimits(double dzv);
  void SetMaxSeedToPropagate(Int_t n=300) {fMaxSeedToPropagate = n;}

 protected:
 
  Int_t fNLayers;        // total number of layers in the model
  Int_t fNActiveLayers;  // number of active layers in the model
  Int_t fNActiveITSLayers;  // number of active ITS layers in the model
  Int_t fLastActiveLayer;       // id of last active layer
  Int_t fLastActiveITSLayer;    // id of last active ITS layer
  Int_t fLastActiveLayerTracked;    // id of last active layer really used for tracking of given pt
  TList fLayers;                // List of layer pointers
  Float_t fBFieldG;             // Magnetic Field in Gauss (set in Tesla)
  Float_t fLhcUPCscale;         // UltraPeripheralElectrons: scale from RHIC to LHC
  Float_t fIntegrationTime;     // electronics integration time
  Float_t fConfLevel;           // Confidence Level for the tracking
  Float_t fAvgRapidity;         // rapidity of the track (= mean)
  Float_t fParticleMass;        // Particle used for tracking. Standard: mass of pion
  Double_t fMaxRadiusSlowDet;   // Maximum radius for slow detectors.  Fast detectors 
                                // and only fast detectors reside outside this radius.
  Int_t fAtLeastCorr;     // min. number of correct hits for the track to be "good"
  Int_t fAtLeastFake;     // min. number of fake hits for the track to be "fake"

  Int_t fdNdEtaCent;       // Multiplicity
  //
  // reconstruction settings
  Double_t fMaxChi2Cl;   // max cluster-track chi2 
  Double_t fMaxNormChi2NDF;// max chi2/NDF to accept
  Int_t    fMinITSHits;  // min ITS hits in track to accept
  //
  Double_t fMaxChi2ClSQ; // precalculated sqrt(chi2);
  //
  Int_t    fMaxSeedToPropagate; // take 1st fMaxSeedToPropagate seeds from each layer
  // background settings
  Bool_t   fUseBackground; // do we want to simulate background?
  Double_t fBgYMin;      // min Y for current layer bg generation
  Double_t fBgYMax;      // max Y for current layer bg generation
  Double_t fBgZMin;      // min Z for current layer bg generation
  Double_t fBgZMax;      // max Z for current layer bg generation
  TArrayD  fBgYMinTr;    // min Y for each seed at current layer
  TArrayD  fBgYMaxTr;    // max Y for each seed at current layer
  TArrayD  fBgZMinTr;    // min Z for each seed at current layer
  TArrayD  fBgZMaxTr;    // max Z for each seed at current layer
  Int_t    fNBgLimits;   // depth of limits array
  //
  enum {kMaxNDetectors = 200};
 
  Double_t fTransMomenta[40];                          // array of transverse momenta
  Double_t fMomentumRes[40];                           // array of momentum resolution
  Double_t fResolutionRPhi[40];                        // array of rphi resolution
  Double_t fResolutionZ[40];                           // array of z resolution
  Double_t fDetPointRes[kMaxNDetectors][40];           // array of rphi resolution per layer
  Double_t fDetPointZRes[kMaxNDetectors][40];          // array of z resolution per layer
  Double_t fEfficiency[3][40];                         // efficiency for different particles
  Double_t fFake[3][40];                               // fake prob for different particles
  //
  KMCProbe fProbe;
  Double_t fDensFactorEta;                             // density scaling for non-0 eta (precalculated
  //
  Double_t fUpdCalls;                  // number of kalman updates
  TH2F*  fHMCLrResidRPhi;              // Residials on lr, rphi
  TH2F*  fHMCLrResidZ;                 // Residials on lr, rphi
  TH2F*  fHMCLrChi2;                   // track to clusyer chi2 on each lr, rphi
  TArrayI fPattITS;                     // bit pattern of each group of active layers where hit is required
  //
  static Double_t fgVtxConstraint[2];  // if both positive, the vertex is used as constraint (accounted in chi2 but not in update)
  ClassDef(KMCDetector,1);
};

//____________________________________________________________________________
inline Bool_t KMCProbe::PropagateToCluster(KMCCluster* cl, double b)
{
  // propagate track to cluster frame
  if (!Rotate(cl->GetPhi()) || !PropagateTo(cl->GetX(),b)) {
    AliDebug(2,Form("Failed to propager track to cluster at phi=%.3f X=%.3f",cl->GetPhi(),cl->GetX()));
    if (AliLog::GetGlobalDebugLevel()>1) Print();
    return kFALSE;
  }
  return kTRUE;
}

//////////////////////////////////////////////////////////////////////////////////////////////
//--------------------------------------------------------------------------------------------
//
// Class to collect summary of certrain track performance
// The track can be selected according to min and max ITS clusters to accept (total, fake and or correct)
// and according to presence of (any and/or correct)  or absence of fakes clusters in certain groups
// of ITS layers.

class KMCTrackSummary : public TNamed
{
 public:
  KMCTrackSummary(const char* name=0,const char* title=0,int nlr=0);
  virtual ~KMCTrackSummary();
  //
  void Add(const KMCTrackSummary* src, double scl=1);

  virtual void  Print(Option_t* option = "") const;
  KMCTrackSummary* MakeCopy(const char* pref) const;
  void   SetNamePrefix(const char* pref="");

  Bool_t CheckTrack(KMCProbe* trc);
  void   AddTrack(KMCProbe* trc);
  //
  TH1*   GetHMCChi2()               const {return fHMCChi2;}
  TH1*   GetHMCSigDCARPhi()         const {return fHMCSigDCARPhi;}
  TH1*   GetHMCSigDCAZ()            const {return fHMCSigDCAZ;}
  TH1*   GetHMCSigPt()              const {return fHMCSigPt;}
  TH2*   GetHMCNCl()                const {return fHMCNCl;}
  TH1*   GetHMCFakePatt()           const {return fHMCFakePatt;}
  //
  Double_t GetCountAcc()            const {return fCountAcc;}
  Double_t GetCountTot()            const {return fCountTot;}
  Double_t GetEff()                 const {return fCountTot>0 ? fCountAcc/fCountTot : 0.;}
  Double_t GetEffErr()              const;
  KMCProbe& GetAnProbe()      const {return (KMCProbe&)fAnProbe;}
  KMCProbe& GetRefProbe()     const {return (KMCProbe&)fRefProbe;}
  void SetAnProbe(KMCProbe* pr)     {if (pr) fAnProbe = *pr;}
  void SetRefProbe(KMCProbe* pr)    {if (pr) fRefProbe = *pr;}
  //
  void SetMinMaxClITS(Int_t mn=0,Int_t mx=999)     {fMinMaxClITS[0]=mn; fMinMaxClITS[1]=mx;}
  void SetMinMaxClITSFake(Int_t mn=0,Int_t mx=999) {fMinMaxClITSFake[0]=mn; fMinMaxClITSFake[1]=mx;}
  void SetMinMaxClITSCorr(Int_t mn=0,Int_t mx=999) {fMinMaxClITSCorr[0]=mn; fMinMaxClITSCorr[1]=mx;}
  //
  Double_t GetUpdCalls()           const {return fCountTot>0 ? fUpdCalls/fCountTot : 0;}
  void     AddUpdCalls(double v)         {fUpdCalls += v;}
  //
  void AddPatternITS(Int_t patt)           {AddPattern(patt,fPattITS);}
  void AddPatternITSFakeExcl(Int_t patt)   {AddPattern(patt,fPattITSFakeExcl);}
  void AddPatternITSCorr(Int_t patt)       {AddPattern(patt,fPattITSCorr);}
  static UInt_t Bits(Bool_t l0=0, Bool_t l1=0, Bool_t l2=0, Bool_t l3=0, Bool_t l4=0, Bool_t l5=0, Bool_t l6=0, Bool_t l7=0,Bool_t  l8=0, Bool_t l9=0,
		     Bool_t l10=0,Bool_t l11=0,Bool_t l12=0,Bool_t l13=0,Bool_t l14=0,Bool_t l15=0,Bool_t l16=0,Bool_t l17=0,Bool_t l18=0,Bool_t l19=0,
		     Bool_t l20=0,Bool_t l21=0,Bool_t l22=0,Bool_t l23=0,Bool_t l24=0,Bool_t l25=0,Bool_t l26=0,Bool_t l27=0,Bool_t l28=0,Bool_t l29=0,
		     Bool_t l30=0,Bool_t l31=0);
  //
  // protected:
  void   AddPattern(Int_t patt,TArrayI& arr) {int n=arr.GetSize(); arr.Set(n+1); arr[n]=patt;}
  Bool_t CheckPattern(UInt_t patt, Int_t cond) {return ((UInt_t)cond)&patt;}
  Int_t  OutOfRange(Int_t n, Int_t range[2])  {if (n<range[0]) return -1; if (n>range[1]) return 1; return 0;}
  static void   PrependTNamed(TNamed* obj, const char* nm=0, const char* tit=0);
  //
 private:
 KMCTrackSummary(const KMCTrackSummary& src) :TNamed(src) {}
 protected:
  //
  Int_t   fNITSLayers;                  // number of its layers
  // track selection conditions
  Int_t   fMinMaxClITS[2];              // min and max N ITS clusters (total)
  Int_t   fMinMaxClITSFake[2];          // min and max N ITS fake clusters  
  Int_t   fMinMaxClITSCorr[2];          // min and max N ITS correct clusters
  //
  TArrayI fPattITS;                     // bit pattern of each group of layers where hit is required
  TArrayI fPattITSFakeExcl;             // bit pattern of each group of layers where fake hit is not tolerated
  TArrayI fPattITSCorr;                 // bit pattern of each group of layers where correc hit is required
  //
  TH1F*  fHMCChi2;                     // track chi2
  TH1F*  fHMCSigDCARPhi;               // DCA distribution in RPhi from MC
  TH1F*  fHMCSigDCAZ;                  // DCA distribution in Z for from MC
  TH1F*  fHMCSigPt;                    // Pt difference distribution from MC
  TH2F*  fHMCNCl;                      // number of clusters
  TH1F*  fHMCFakePatt;                 // Fakes pattern
  //
  Double_t fCountAcc;                  // track counter of accepted tracks
  Double_t fCountTot;                  // track counter of all tracks
  Double_t fUpdCalls;                  // total update calls
  //
  KMCProbe fRefProbe;                  // reference probe
  KMCProbe fAnProbe;                   // analital probe
  static Int_t fgSumCounter;           // global counter of bins

  ClassDef(KMCTrackSummary,1)
};

inline Double_t KMCTrackSummary::GetEffErr() const {
  //
  if (fCountTot<1) return 0;
  double r = fCountAcc/fCountTot;
  double err = r*(1.-r)/fCountTot;
  return err>0 ? TMath::Sqrt(err) : 0;
}

inline KMCTrackSummary* KMCTrackSummary::MakeCopy(const char* pref) const {
  // create a copy, prepending all names with prefix
  KMCTrackSummary* sm = (KMCTrackSummary*) Clone(this->GetName()); 
  sm->SetNamePrefix(pref); 
  return sm;
} 

#endif
