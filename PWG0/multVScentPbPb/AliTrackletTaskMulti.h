#ifndef ALITRACKLETTASKMULTI_H
#define ALITRACKLETTASKMULTI_H

///////////////////////////////////////////////////////////////////////////
// Class AliTrackletTaskMulti                                            //
// Analysis task to produce data and MC histos needed for tracklets      //
// dNdEta extraction in multiple bins in one go                          //
// Author:  ruben.shahoyan@cern.ch                                       //
///////////////////////////////////////////////////////////////////////////

class TH1F; 
class TH2F;
class TH3F;
class AliESDEvent;
class TList;
class TNtuple;

class AliMCParticle;
class AliITSMultRecBg;
class AliESDTrackCuts;

#include "../ITS/AliITSsegmentationSPD.h"
#include "AliAnalysisTaskSE.h"
#include "AliTriggerAnalysis.h" 
#include <TMath.h>

class AliTrackletTaskMulti : public AliAnalysisTaskSE {
 public:
  enum {kData,kBgInj,kBgRot,kBgMix,kMC};
  enum {kCentV0M,kCentFMD,kCentTRK,kCentTKL,kCentCL0,kCentCL1,kCentV0MvsFMD,kCentTKLvsV0,kCentZEMvsZDC,kNCentTypes}; // what is used to define centrality
  //
  enum {  // define here id's of the standard histos in corresponding TObjArray* fHistosTr...
    kHEtaZvCut,       // histo zv vs eta for tracklets passing final selection (dist<1 or |dPhi|<narrowWindow ...)
    kHDPhiDTheta,     // measured dTheta vs dPhi
    kHDPhiSDThetaX,   // dTheta (1/sin^2 scaled if needed) vs dPhi (bending subtracted)
    kHWDist,          // Weighted distance 
    kNStandardH       // number of standard histos per centrality bin
  };
  enum { // define here id's of any custom histos to be added to fHistosCustom
    kHStat,            // job info (meaning of bins defined in the enum below)
    //
    kHStatCent,        // events per centrality bin with real values on the axis
    kHStatCentBin,     // events per centrality bin
    //
    kHNPrimMeanMC,     // <n> primaries per mult bin
    kHNPrim2PartMC,    // <n> prim per part.pair per mult bin
    kHNPrim2BCollMC,   // <n> prim per bin.coll per mult bin
    kHNPrim2PartNpMC,  // <n> prim per n part vs npart
    kHNPrim2BCollNpMC, // <n> prim per n part vs npart
    kHNPartMC,         // n.part.pairs according to MC
    kHNPartMeanMC,     // <n> part pairs per mult bin
    kHNBCollMC,        // n.bin.colls according to MC
    kHNBCollMeanMC,
    //
    kHZVtxNoSel,       // Z vertex distribution before event selection
    kHV0NoSel,         // V0 before selection
    kHNClSPD2NoSel,    // NSPD2 before selection
    kHZDCZEMNoSel,     // ZDC ZEM before selection
    //
    kHZVtx,            // Z vertex distribution
    kHV0,              // V0 before selection
    kHNClSPD2,         // NSPD2 before selection
    kHZDCZEM,          // ZDC ZEM before selection
    

    kHZVtxMixDiff,     // difference in Z vtx of mixed events
    kHNTrMixDiff,      // difference in N tracklets of mixed events
    //
    kHPrimPDG,         // PDG code of prim tracklet
    kHSecPDG,          // PDG code of sec tracklet
    kHPrimParPDG,      // PDG code of prim tracklet parent
    kHSecParPDG,       // PDG code of sec tracklet parent
    //
    kHClUsedInfoL0,    // used clusters of lr0
    kHClUsedInfoL1,    // used clusters of lr1
    kHClAllInfoL0,     // all clusters of lr0
    kHClAllInfoL1,     // all clusters of lr1
    //
    // This MUST be last one: this is just beginning of many histos (one per bin)
    kHZVEtaPrimMC      // Zv vs eta for all primary tracks (true MC multiplicity)
  }; // custom histos

  // bins for saved parameters
  enum {kDummyBin,
	kEvTot0,      // events read
	kEvTot,       // events read after vertex quality selection
	kOneUnit,     // just 1 to track primate merges
	kNWorkers,    // n workers
	//
	kCentVar,     // cetrality var. used
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
	//
	kDPiSCut,     // cut on dphi used to extract signal (when WDist is used in analysis, put it equal to kDPhi
	kNStdCut,     // cut on weighted distance (~1) used to extract signal 
	//
	kMCV0Scale,   // scaling value for V0 in MC
	//
	// here we put entries for each mult.bin
	kBinEntries = 50,
	kEvProcData,  // events with data mult.object (ESD or reco)
	kEvProcInj,   // events Injected, total
	kEvProcRot,   // events Rotated
	kEvProcMix,   // events Mixed
	kEntriesPerBin
  };

  //
  AliTrackletTaskMulti(const char *name = "AliTrackletTaskMulti");
  virtual ~AliTrackletTaskMulti(); 
  
  virtual void  UserCreateOutputObjects();
  virtual void  UserExec(Option_t *option);
  virtual void  Terminate(Option_t *);

  void       SetUseCentralityVar(Int_t v=kCentV0M)     {fUseCentralityVar = v;}
  void       SetUseMC(Bool_t mc = kFALSE)              {fUseMC = mc;}
  void       SetCheckReconstructables(Bool_t c=kFALSE) {fCheckReconstructables = c;}
  TObjArray* BookHistosSet(const char* pref, UInt_t selHistos=0xffffffff);
  TObjArray* BookCustomHistos();
  void       AddHisto(TObjArray* histos, TObject* h, Int_t at=-1);
  void       FillHistosSet(TObjArray* histos, double eta, /*double phi,double theta,*/double dphi,double dtheta,double dthetaX,double dist);
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
  //
  void       SetDPhiSCut(Float_t c=0.06)        {fDPhiSCut = c;}
  void       SetNStdCut(Float_t c=1.0)          {fNStdCut = c;}
  void       SetScaleMCV0(Float_t s=1.0)        {fMCV0Scale = s;}  
  //
  void       SetEtaCut(Float_t etaCut)          {fEtaMax = TMath::Abs(etaCut); fEtaMin= -fEtaMax;}
  void       SetEtaMin(Float_t etaMin)          {fEtaMin = etaMin;}
  void       SetEtaMax(Float_t etaMax)          {fEtaMax = etaMax;}
  void       SetZVertexMin(Float_t z)           {fZVertexMin = z;}
  void       SetZVertexMax(Float_t z)           {fZVertexMax = z;}
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
  //
 protected:
  void       InitMultReco();
  Bool_t     HaveCommonParent(const float* clLabs0,const float* clLabs1);
  void       FillHistos(Int_t type, const AliMultiplicity* mlt);
  void       FillMCPrimaries();
  void       FillSpecies(Int_t primsec, Int_t id);
  void       FillClusterInfo();
  Int_t      GetPdgBin(Int_t pdgCode);
  void       CheckReconstructables();
  Int_t      GetCentralityBin(Float_t percentile) const;
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
  Float_t      fDPhiSCut;                  // cut on signal dphiS
  Float_t      fNStdCut;                   // cut on signal weighted distance
  Float_t      fMCV0Scale;                 // scaling factor for V0 in MC
  //
  AliITSMultRecBg *fMultReco;              //! mult.reco object
  TTree*       fRPTree;                    //! tree of recpoints
  TTree*       fRPTreeMix;                 //! tree of recpoints for mixing
  AliStack*    fStack;                     //! MC stack
  AliMCEvent*  fMCEvent;                   //! MC Event
  Float_t      fESDVtx[3];                 //  ESD vertex
  //
  Float_t fNPart;                          // number of participant pairs from MC
  Float_t fNBColl;                         // number of bin. collision from MC
  Int_t  fCurrCentBin;                     // current centrality bin
  Int_t  fNCentBins;                       // N of mult bins
  Int_t  fUseCentralityVar;                // what is used to determine the centrality
  //
  static const Float_t fgkCentPerc[];               //! centrality in percentiles
  //
  static const char*  fgCentSelName[];              //!centrality types
  static const char*  fgkPDGNames[];                //!pdg names
  static const Int_t  fgkPDGCodes[];                //!pdg codes
  //
 private:    
  AliTrackletTaskMulti(const AliTrackletTaskMulti&); // not implemented
  AliTrackletTaskMulti& operator=(const AliTrackletTaskMulti&); // not implemented 
  
  ClassDef(AliTrackletTaskMulti, 1);  
};


#endif
