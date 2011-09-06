#ifndef ALITRACKLETTASKUNI_H
#define ALITRACKLETTASKUNI_H

///////////////////////////////////////////////////////////////////////////
// Class AliTrackletTask                                                 //
// Analysis task to study performance of tracklet reconstruction         //
// algorithm and combinatorial background                                //
// Author:  M. Nicassio (INFN Bari)                                      //
// Contact: Maria.Nicassio@ba.infn.it, Domenico.Elia@ba.infn.it          //
///////////////////////////////////////////////////////////////////////////

class TH1F; 
class TH2F;
class AliESDEvent;
class TList;
class TNtuple;

class AliMCParticle;
class AliITSMultRecBg;

#include "../ITS/AliITSsegmentationSPD.h"
#include "AliAnalysisTaskSE.h"
#include "AliTriggerAnalysis.h" 

class AliTrackletTaskUni : public AliAnalysisTaskSE {
 public:
  enum {kData,kBgInj,kBgRot,kBgMix,kMC};
  //
  enum {  // define here id's of the standard histos in corresponding TObjArray* fHistosTr...
    kHEtaZvDist,      // 3 d sparse histo with dist  (uncut) vs zv vs eta
    kHEtaZvDPhiS,     // 3 d sparse histo with dphiS (uncut) vs zv vs eta
    kHEtaZvCut,       // zv vs eta with strict cut on tracklets applied (dist<1 or |dPhi|<narrowWindow)
    kHDPhiDTheta,     // measured dTheta vs dPhi
    kHDPhiSDThetaX,   // dTheta (1/sin^2 scaled if needed) vs dPhi (bending subtracted)
    kHEtaDPhiS,       // dPhi (bending subtracted) vs eta
    kHEtaDThetaX,     // dTheta (1/sin^2 scaled if needed) vs eta
    kHEtaDist,        // Weighted distance vs eta
    kHZvDPhiS,        // dPhi (bending subtracted) vs Zv
    kHZvDThetaX,      // dTheta (1/sin^2 scaled if needed) vs Zv
    kHZvDist          // Weighted distance vs Zv
  };
  enum { // define here id's of any custom histos to be added to fHistosCustom
    kHStat,            // job info (meaning of bins defined in the enum below)
    kHZVEtaPrimMC,     // Zv vs eta for all primary tracks (true MC multiplicity)
    //
    kHZVtxNoSel,       // Z vertex distribution before event selection
    kHNTrackletsNoSel, // N tracklets before event selection
    kHNClSPD1NoSel,    // N clusters on SPD1 before event selection
    kHNClSPD2NoSel,    // N clusters on SPD2 before event selection
    kHV0NoSel,         // V0 mult before selection
    kHV0NClSPD2NoSel,  // V0 - nspd2 correlation
    //
    kHZVtx,            // Z vertex distribution
    kHNTracklets,      // N tracklets
    kHNClSPD1,         // N clusters on SPD1
    kHNClSPD2,         // N clusters on SPD2
    kHV0,              // V0 mult after selection
    //
    kHZVtxMixDiff,     // difference in Z vtx of mixed events
    kHNTrMixDiff,      // difference in N tracklets of mixed events
    //
    kHPrimPDG,         // PDG code of prim tracklet
    kHSecPDG,          // PDG code of sec tracklet
    kHPrimParPDG,      // PDG code of prim tracklet parent
    kHSecParPDG        // PDG code of sec tracklet parent
  }; // custom histos

  // bins for saved parameters
  enum {kDummyBin,
	kEvTot0,      // events read
	kEvTot,       // events read after vertex quality selection
	kEvTotPlp,    // events with pile-up
	kEvProcData,  // events with data mult.object (ESD or reco)
	kEvProcInj,   // events Injected
	kEvProcRot,   // events Rotated
	kEvProcMix,   // events Mixed
	//
	kDPhi,        // dphi window
	kDTht,        // dtheta window
	kNStd,        // N.standard deviations to keep
	kPhiShift,    // bending shift
	kThtS2,       // is dtheta scaled by 1/sin^2
	kThtCW,       // on top of w.dist cut cut also on 1 sigma dThetaX
	kPhiOvl,      // overlap params
	kZEtaOvl,     // overlap params
	kNoOvl,       // flag that overlap are suppressed
	//
	kPhiRot,      // rotation phi
	kInjScl,      // injection scaling
	kEtaMin,      // eta cut
	kEtaMax,      // eta cut
	kZVMin,       // min ZVertex to process
	kZVMax,       // max ZVertex to process
	kTrcMin,      // min mult to process
	kTrcMax,      // max mult to process
	//
	kMCV0Scale,   // scaling value for V0 in MC	
	//
	kOneUnit=49,  // just 1 to track mergings
	kNWorkers=50, // n workers
	kNStatBins
  };

  //
  AliTrackletTaskUni(const char *name = "AliTrackletTaskUni");
  virtual ~AliTrackletTaskUni(); 
  
  virtual void  UserCreateOutputObjects();
  virtual void  UserExec(Option_t *option);
  virtual void  Terminate(Option_t *);
  void          RegisterStat();

  void       SetUseMC(Bool_t mc = kFALSE)              {fUseMC = mc;}
  void       SetCheckReconstructables(Bool_t c=kFALSE) {fCheckReconstructables = c;}
  TObjArray* BookHistosSet(const char* pref, UInt_t selHistos=0xffffffff);
  TObjArray* BookCustomHistos();
  void       AddHisto(TObjArray* histos, TObject* h, Int_t at=-1);
  void       FillHistosSet(TObjArray* histos, double phi,double theta,double dphi,double dtheta,double dist);
  // RS
  void       SetNStdDev(Float_t f=1.)           {fNStdDev = f<1e-5 ? 1e-5:f;}
  void       SetScaleDThetaBySin2T(Bool_t v=kFALSE) {fScaleDTBySin2T = v;}
  void       SetCutOnDThetaX(Bool_t v=kFALSE)   {fCutOnDThetaX = v;}
  void       SetPhiWindow(float w=0.08)         {fDPhiWindow   = w<1e-5 ? 1e-5:w;}
  void       SetThetaWindow(float w=0.025)      {if (w<0) fCutOnDThetaX=kTRUE; fDThetaWindow = TMath::Abs(w)<1e-5 ? 1e-5:TMath::Abs(w);}
  void       SetPhiShift(float w=0.0045)        {fDPhiShift = w;}
  void       SetPhiOverlapCut(float w=0.005)    {fPhiOverlapCut = w;}
  void       SetZetaOverlapCut(float w=0.05)    {fZetaOverlap = w;}
  void       SetPhiRot(float w=0)               {fPhiRot = w;}
  void       SetInjScale(Float_t s=1.)          {fInjScale = s>0? s:1.;}
  void       SetRemoveOverlaps(Bool_t w=kFALSE) {fRemoveOverlaps = w;}
  void       SetScaleMCV0(Float_t s=1.0)        {fMCV0Scale = s;}  
  //
  void       SetEtaCut(Float_t etaCut)          {fEtaMax = TMath::Abs(etaCut); fEtaMin= -fEtaMax;}
  void       SetEtaMin(Float_t etaMin)          {fEtaMin = etaMin;}
  void       SetEtaMax(Float_t etaMax)          {fEtaMax = etaMax;}
  void       SetZVertexMin(Float_t z)           {fZVertexMin = z;}
  void       SetZVertexMax(Float_t z)           {fZVertexMax = z;}
  void       SetMultCutMin(Int_t n=0)           {fMultCutMin = n;}
  void       SetMultCutMax(Int_t n=99999)       {fMultCutMax = n;}
  //
  Bool_t     GetDoNormalReco()             const {return fDoNormalReco;}
  Bool_t     GetDoInjection()              const {return fDoInjection;}
  Bool_t     GetDoRotation()               const {return fDoRotation;}
  Bool_t     GetDoMixing()                 const {return fDoMixing;}
  //
  void       SetDoNormalReco(Bool_t v=kTRUE)    {fDoNormalReco = v;}
  void       SetDoInjection(Bool_t v=kTRUE)     {fDoInjection = v;}
  void       SetDoRotation(Bool_t v=kTRUE)      {fDoRotation = v;}
  void       SetDoMixing(Bool_t v=kTRUE)        {fDoMixing = v;}
  //
  void       SetDontMerge(Bool_t v=kTRUE)       {fDontMerge = v;}
  /*
  void       SetTrigger(AliTriggerAnalysis::Trigger trigger)  { fTrigger = trigger; }
  void       SetMCCentralityBin(MCCentralityBin mccentrbin)   { fMCCentralityBin = mccentrbin;}
  void       SetCentralityLowLim(Float_t centrlowlim)         { fCentrLowLim = centrlowlim;}
  void       SetCentralityUpLim(Float_t centruplim)           { fCentrUpLim = centruplim;}
  void       SetCentralityEst(TString centrest)               { fCentrEst = centrest;}
  */
  //
 protected:
  void       InitMultReco();
  Bool_t     HaveCommonParent(const float* clLabs0,const float* clLabs1);
  void       FillHistos(Int_t type, const AliMultiplicity* mlt);
  void       FillMCPrimaries(TH2F* hetaz);
  void       FillSpecies(Int_t primsec, Int_t id, Double_t dist);
  Int_t      GetPdgBin(Int_t pdgCode);
  void       CheckReconstructables();
  //
 protected:
  TList*       fOutput;                   // output list send on output slot 1 
  //
  Bool_t       fDoNormalReco;              // do normal reco
  Bool_t       fDoInjection;               // do injection
  Bool_t       fDoRotation;                // do rotation
  Bool_t       fDoMixing;                  // do mixing
  //
  Bool_t       fUseMC; 
  Bool_t       fCheckReconstructables;
  //
  TObjArray*   fHistosTrData;              //! all tracklets in data
  TObjArray*   fHistosTrInj;               //! injected
  TObjArray*   fHistosTrRot;               //! rotated
  TObjArray*   fHistosTrMix;               //! mixed
  //
  TObjArray*   fHistosTrPrim;              //! primary
  TObjArray*   fHistosTrSec;               //! secondary
  TObjArray*   fHistosTrComb;              //! combinatorials
  TObjArray*   fHistosTrCombU;             //! combinatorials uncorrelated
  //
  TObjArray*   fHistosTrRcblPrim;          //! Primary Reconstructable
  TObjArray*   fHistosTrRcblSec;           //! Secondary Reconstructable
  TObjArray*   fHistosCustom;              //! custom histos
  //
  // Settings for the reconstruction
  // tracklet reco settings
  Float_t      fEtaMin;                    // histos filled only for this eta range
  Float_t      fEtaMax;                    // histos filled only for this eta range
  Float_t      fZVertexMin;                // min Z vtx to process
  Float_t      fZVertexMax;                // max Z vtx to process
  Int_t        fMultCutMin;                // min mult in ESD to process?
  Int_t        fMultCutMax;                // max mult in ESD to process?
  Float_t      fMCV0Scale;                 // scaling factor for V0 in MC
  //
  Bool_t       fScaleDTBySin2T;            // request dTheta scaling by 1/sin^2(theta)
  Bool_t       fCutOnDThetaX;              // if true, apart from NStdDev cut apply also the cut on dThetaX
  Float_t      fNStdDev;                   // cut on weighted distance
  Float_t      fDPhiWindow;                // max dPhi
  Float_t      fDThetaWindow;              // max dTheta
  Float_t      fDPhiShift;                 // mean bend
  Float_t      fPhiOverlapCut;             // overlaps cut in phi
  Float_t      fZetaOverlap;               // overlaps cut in Z
  Float_t      fPhiRot;                    // rotate L1 wrt L2
  Float_t      fInjScale;                  // scaling factor for injection
  Bool_t       fRemoveOverlaps;            // request overlaps removal
  //
  AliITSMultRecBg *fMultReco;              //! mult.reco object
  TTree*       fRPTree;                    //! tree of recpoints
  TTree*       fRPTreeMix;                 //! tree of recpoints for mixing
  AliStack*    fStack;                     //! MC stack
  AliMCEvent*  fMCEvent;                   //! MC Event
  Float_t      fESDVtx[3];                 //  ESD vertex
  //
  /*
  AliTriggerAnalysis::Trigger fTrigger;    // requested trigger
  MCCentralityBin fMCCentralityBin;        // to select MC centrality bin in which corrections are calculated
  Float_t      fCentrLowLim;               // to select centrality bin on data
  Float_t      fCentrUpLim;                // to select centrality bin on data
  TString      fCentrEst;                  // to select centrality estimator
  */
  Bool_t fDontMerge;                       // no merging requested
  static const char*  fgkPDGNames[];                //!pdg names
  static const Int_t  fgkPDGCodes[];                //!pdg codes
  //
 private:    
  AliTrackletTaskUni(const AliTrackletTaskUni&); // not implemented
  AliTrackletTaskUni& operator=(const AliTrackletTaskUni&); // not implemented 
  
  ClassDef(AliTrackletTaskUni, 1);  
};


#endif
