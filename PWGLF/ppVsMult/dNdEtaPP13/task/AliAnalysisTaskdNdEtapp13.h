#ifndef AliAnalysisTaskdNdEtapp13_H
#define AliAnalysisTaskdNdEtapp13_H

///////////////////////////////////////////////////////////////////////////
// Class AliAnalysisTaskdNdEtapp13                                            //
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
class AliStack;
class AliMultiplicity;

class AliPPVsMultUtils;

#include "AliITSsegmentationSPD.h"
#include "AliAnalysisTaskSE.h"
#include "AliTriggerAnalysis.h"
#include "AliVEvent.h"
#include <TMath.h>

class AliAnalysisTaskdNdEtapp13 : public AliAnalysisTaskSE {
public:
  enum {kData,kBgInj,kBgRot,kMC};
  enum {kCentV0M,kCentV0A,kCentV0C,kCentFMD,kCentTRK,kCentTKL,kCentCL0,kCentCL1,kCentV0MvsFMD,kCentZNA,kCentTKLvsV0,kCentZEMvsZDC,kCentV0A123,kCentV0A0,kCentV0S,kCentMB,kCentRef,kCentV0av,kNCentTypes}; // what is used to define centrality
  //
  enum {  // define here id's of the standard histos in corresponding TObjArray* fHistosTr...
  kHEtaZvCut,       // histo zv vs eta for tracklets passing final selection (dist<1 or |dPhi|<narrowWindow ...)
  kHDPhiDTheta,     // measured dTheta vs dPhi
  kHDPhiSDThetaX,   // dTheta (1/sin^2 scaled if needed) vs dPhi (bending subtracted)
  kHWDist,          // Weighted distance
  //    kHEtaZvSPD1,      // histo zv vs eta for SPD1 single clusters
  kHWDvEta,         // WDist vs eta
  kNStandardH       // number of standard histos per centrality bin
};
enum { // define here id's of any custom histos to be added to fHistosCustom
kHStat,            // job info (meaning of bins defined in the enum below)
//
kHStatCent,        // events per centrality bin with real values on the axis
kHStatCentBin,     // events per centrality bin
kHCentDistNoSel,   // events per centrality percentile before selection
kHCentDistTrig,   // events per centrality percentile after trigger
kHCentDist,        // events per centrality percentile
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
kHNPrimEta2All,    // NPrim |eta|<2 all gen events
kHNPrimEta2Vt,     // NPrim |eta|<2 events with vtx
kHNPrimEta2Sel,    // NPrim |eta|<2 events after selection
kHNPrimEta2SelVt,  // NPrim |eta|<2 events passing selection with vertex
kHNCorrMCEta2,     // Ntracklets |eta|<2 vs NPrim |eta|<2
//
//
kHZVtxNoSel,       // Z vertex distribution before event selection
kHV0NoSel,         // V0 before selection
kHV0Trig,         // V0 after trigger
kHNClSPD2NoSel,    // NSPD2 before selection
kHZDCZEMNoSel,     // ZDC ZEM before selection
//
kHTotalNchNoPhSel,
kHProcessMCNoPhSel,
//
kHZVtxMCNoPhSel,   // Z vertex distribution MC before physics selection
kHZVtxMCNoVtSel,   // Z vertex distribution MC before rec.vertex selection
kHZVtxMC,          // Z vertex distribution MC after all selections
kHZVtx,            // Z vertex distribution
kHV0,              // V0 before selection
kHNClSPD2,         // NSPD2 before selection
kHZDCZEM,          // ZDC ZEM before selection
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
kHMltEstTrITSTPC,  // mult distribution with kTrackletsITSTPC
kHMltEstTrITSSA,   // mult distribution with kTrackletsITSSA
kHMltEstTr,        // mult distribution with tracklets
kHMltMC,           // mult distribution for MC primaries
//
kHEffMatrix,
//
kHEtaPhi,          // traklets eta/phi
// These MUST be last one: this is just beginning of many histos (one per bin)
kHZVEtaPrimMC,                        // Zv gen vs eta for all primary tracks (true MC multiplicity) in all events
kHZVrEtaPrimMC=kHZVEtaPrimMC+50,      // Zv rec vs eta for all primary tracks (true MC multiplicity) in sel. events
kHZVResMC=kHZVrEtaPrimMC+50,      // zv resolution

kHCorrMatrix=kHZVResMC+50, // correlation matrix

kHCorrMatrixSel = kHCorrMatrix+50, // correlation matrix with events sel

kHCorrMatrixSel2 = kHCorrMatrixSel+50 // correlation matrix 2 with tracklets (triggered + Selected)


}; // custom histos

// bins for saved parameters
enum {kDummyBin,
  kEvTot0,      // events read
  kEvTot,       // events read after vertex quality selection
  kEvAfterPhysSel  , // events after pileup rejection
  kEvAfterPileUp  , // events after pileup rejection
  kEvAfterClsVsTrk, // events after the cluster vs tracklet cut
  kEvAfterAsymCut , // events after V0 asymmetry cut
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
    kNEvSDMC,     // number of pure SD event from MC
    // here we put entries for each mult.bin
    kBinEntries = 50,
    kEvInMltBin = 0,
    kEvProcData,  // events with data mult.object (ESD or reco) passing all selections
    kEvProcInj,   // events Injected, total
    kEvProcRot,   // events Rotated
    kEvProcMix,   // events Mixed
    kEvCentBin,   // events in centrality bin (w/o any selection)
    kEvPassPS,    // events passed ev.sel + trigger
    kEvPassVtx,   // and having vertex (not necessarily in needed range)
    kEntriesPerBin = 10
  };

  //
  AliAnalysisTaskdNdEtapp13(const char *name = "AliAnalysisTaskdNdEtapp13");
  virtual ~AliAnalysisTaskdNdEtapp13();

  virtual void  UserCreateOutputObjects();
  virtual void  UserExec(Option_t *option);
  virtual void  Terminate(Option_t *);
  void       RegisterStat();
  //
  void       SetTriggerSelection(UInt_t sel=AliVEvent::kINT7) {fTrigSel = sel;}
  void       SetCentPercentiles(const Float_t *arr, Int_t nbins);
  void       SetCentPercentiles(const Double_t *arr, Int_t nbins);
  //
  void       CheckCentralityVar(const char* var);
  void       SetUseCentralityVar(Int_t v)              {fUseCentralityVar = fgCentSelName[v];}
  void       SetUseCentralityVar(const char* v="V0M")  {CheckCentralityVar(v); fUseCentralityVar = v;}
  void       SetUseMC(Bool_t mc = kFALSE)              {fUseMC = mc;}
  void       SetCheckReconstructables(Bool_t c=kFALSE) {fCheckReconstructables = c;}
  TObjArray* BookHistosSet(const char* pref, UInt_t selHistos=0xffffffff);
  TObjArray* BookCustomHistos();
  void       AddHisto(TObjArray* histos, TObject* h, Int_t at=-1);
  Bool_t     FillHistosSet(TObjArray* histos, double eta, /*double phi,double theta,*/double dphi,double dtheta,double dthetaX,double dist);
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
  void       SetEtaBin(Float_t deta)            {fEtaBin = deta>0.01 ? deta : 0.01;}
  Float_t    GetEtaBin()                  const {return fEtaBin;}

  void       SetZVertexMin(Float_t z)           {fZVertexMin = z;}
  void       SetZVertexMax(Float_t z)           {fZVertexMax = z;}
  void       SetZVBin(Float_t dz=1)             {fZVBin = dz>0.1 ? dz : 0.1;}
  Float_t    GetZVBin()                    const  {return fZVBin;}
  //
  Bool_t     GetDoNormalReco()             const {return fDoNormalReco;}
  Bool_t     GetDoInjection()              const {return fDoInjection;}
  Bool_t     GetDoRotation()               const {return fDoRotation;}
  //
  void       SetDoNormalReco(Bool_t v=kTRUE)    {fDoNormalReco = v;}
  void       SetDoInjection(Bool_t v=kTRUE)     {fDoInjection = v;}
  void       SetDoRotation(Bool_t v=kTRUE)      {fDoRotation = v;}
  //
  Bool_t GetUseSpecialOutput()            const {return fUseSpecialOutput;}
  void   SetUseSpecialOutput(Bool_t v=kTRUE)    {fUseSpecialOutput=v;}
  void       SetUseBCMod(Bool_t bc = kFALSE)              {fUseBCMod = bc;}
  void       SetBCMod(Int_t mod=2)        {fBCMod4 = mod;}
  void       SetCutOnPhi(Bool_t p=kFALSE)   {fCutOnPhi = p;}
  void SetCalibfilePath(TString CalibfilePath) {fCalibfilePath=CalibfilePath;}
  void SetCalibHisto(TH1D *calibhisto=0x0);



  //
protected:
  void       InitMultReco();
  Bool_t     HaveCommonParent(const float* clLabs0,const float* clLabs1);
  void       FillHistos(Int_t type, const AliMultiplicity* mlt);
  void       FillMCPrimaries();
  void       FillSpecies(Int_t primsec, Int_t id);
  void       FillClusterInfo();
  void       FillClusterInfoFromMult(const AliMultiplicity* mlt, double zVertex);
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
  //
  Bool_t       fUseMC;
  Bool_t       fCheckReconstructables;
  //
  TObjArray*   fHistosTrData;              //! all tracklets in data
  TObjArray*   fHistosTrInj;               //! injected
  TObjArray*   fHistosTrRot;               //! rotated
  //
  TObjArray*   fHistosTrPrim;              //! primary
  TObjArray*   fHistosTrSec;               //! secondary
  TObjArray*   fHistosTrComb;              //! combinatorials
  TObjArray*   fHistosTrCombU;             //! combinatorials uncorrelated
  //
  TObjArray*   fHistosTrRcblPrim;          //! Primary Reconstructable
  TObjArray*   fHistosTrRcblSec;           //! Secondary Reconstructable
  TObjArray*   fHistosCustom;              //! custom histos
  TH1D*        fHcalib3;

  //
  // Settings for the reconstruction
  // tracklet reco settings
  Float_t      fEtaMin;                    // histos filled only for this eta range
  Float_t      fEtaMax;                    // histos filled only for this eta range
  Float_t      fEtaBin;                    // eta bin size
  Float_t      fZVertexMin;                // min Z vtx to process
  Float_t      fZVertexMax;                // max Z vtx to process
  Float_t      fZVBin;                     // bin size in Zv
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
  UInt_t       fTrigSel; // requested trigger (to be changed to AliBits)

  AliITSMultRecBg *fMultReco;              //! mult.reco object
  TTree*       fRPTree;                    //! tree of recpoints
  AliStack*    fStack;                     //! MC stack
  AliMCEvent*  fMCEvent;                   //! MC Event
  Float_t      fESDVtx[3];                 //! ESD vertex
  Float_t      fVtxMC[3];                  //! MC gen vertex
  //
  Int_t   fNPrimMCeta2;                    //! N of primaries |eta|<2
  Int_t   fNTreta2;                        //! N of signal tracklets |eta|<2
  Float_t fNPart;                          // number of participant pairs from MC
  Float_t fNBColl;                         // number of bin. collision from MC
  Int_t  fCurrCentBin;                     // current centrality bin
  Int_t  fNCentBins;                       // N of mult bins
  TArrayF fCentPerc;                       // centrality in percentiles
  TString fUseCentralityVar;                // what is used to determine the centrality
  Bool_t fIsSelected;                      //! did current event pass phys.sel.?
  Bool_t fVtxOK;                           //! rec.vertex is good
  Bool_t fUseSpecialOutput;                // flag to open special output
  Bool_t fUseBCMod;                         // flag to use bunch crossing mod 4 events
  Int_t      fBCMod4;                        // Select BC Mod4
  Bool_t fCutOnPhi;                          // Set cut on affected phi regions
  TString fCalibfilePath;                     // Set V0M MC calibration file name path



  //
  static const char*  fgCentSelName[];              //!centrality types
  static const char*  fgkPDGNames[];                //!pdg names
  static const Int_t  fgkPDGCodes[];                //!pdg codes
  //
  Double_t fWeight;
  AliPPVsMultUtils *fPPVsMultUtils; //!
private:
  AliAnalysisTaskdNdEtapp13(const AliAnalysisTaskdNdEtapp13&); // not implemented
  AliAnalysisTaskdNdEtapp13& operator=(const AliAnalysisTaskdNdEtapp13&); // not implemented

  ClassDef(AliAnalysisTaskdNdEtapp13, 1);
};


#endif
