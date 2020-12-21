#ifndef ALISINGLETRACKEFFICUTS_H
#define ALISINGLETRACKEFFICUTS_H

#include <TString.h>
#include "TObject.h"

#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliAnalysisCuts.h"
#include "AliPID.h"
#include "AliAODTrack.h"
#include "AliPIDCombined.h"


class AliSingleTrackEffCuts : public AliAnalysisCuts 
{
 public:

  enum {
      kNoBayesianPID=0,         // no Bayesian PID applied (default)
      kMaximumBayesianProb=1,   // maximum Bayesian PID probability
      kThresholdBayesianProb=2  // Bayesian PID with probability Threshold
  };

  AliSingleTrackEffCuts();
  AliSingleTrackEffCuts(const char* name, const char* title);

  AliSingleTrackEffCuts(const AliSingleTrackEffCuts& source);
  AliSingleTrackEffCuts& operator=(const AliSingleTrackEffCuts& source);

  virtual ~AliSingleTrackEffCuts();

  Bool_t IsSelected(TList* list) { 
    AliWarning(Form(" Function not implemented, list having %d entries",list->GetEntries())); 
    return kFALSE; 
  } // not implemented

  //
  // Event and track selections
  //
  // Event selection at MC level (z-vtx)
  Bool_t IsMCEventSelected(TObject *obj);
  // Event selection at reconstructed level (trigger, zvtx)
  Bool_t IsRecoEventSelected(TObject *obj);
  //
  // Check generated particles (primary, charged, pdgcode)
  Bool_t IsMCParticleGenerated(TObject *obj);
  // Check generated particles (eta, pt)
  Bool_t IsMCParticleInKineAcceptance(TObject *obj);
  // Check if particle has left enough hits in the detectors (only at ESD level)
  Bool_t IsMCParticleInReconstructable(TObject *obj);
  // Check if reconstructed particle is in the acceptance (eta, pt)
  Bool_t IsRecoParticleKineAcceptance(TObject *obj);
  // Check if reconstructed particle passes the PID selections
  //  (only at AOD level for now) if TPC or TOF accept
  Bool_t IsRecoParticlePID(TObject *obj);

  //
  // Setters
  //
  // Set maximum distance of particle origin from interaction point
  void SetMaxRadiusOfParticleOrigin(Float_t rmax){fMaxProdRadius=rmax;}
  // Set eta range for acceptance cut (both MC and reco level)
  void SetEtaRange(Float_t etamin, Float_t etamax){ fEtaMin=etamin; fEtaMax=etamax; }
  // Set rapidity range for acceptance cut (both MC and reco level)
  void SetYRange(Float_t ymin, Float_t ymax){ fYMin=ymin; fYMax=ymax; }
  // Set pt range for acceptance cut (both MC and reco level)
  void SetPtRange(Float_t ptmin, Float_t ptmax){ fPtMin=ptmin; fPtMax=ptmax; }
  // Set PDG code to be checked (by default =0, no check)
  void SetPdgCode(Int_t pdgCode){ fPdgCode = pdgCode; fIsPdgCode=kTRUE; }
  Int_t GetPdgCode() { return fPdgCode; }
  // Set if check for charged particle
  void SetIsCharged(Bool_t charge){ fIsCharged=charge; }

  // Set/get for AOD/ESD analysis
  void SetIsAOD(Bool_t flag){ fisAOD = flag; }
  Bool_t IsAOD(){ return fisAOD; }

  // Set for vertex type (0: not cut; 1: SPDZ; 2: SPD3D; 3: Tracks)
  void SetMinVtxType(Int_t type=3) { fMinVtxType=type; }
  void SetUseEventsWithOnlySPDVertex(Bool_t flag=kTRUE){ 
    if(flag) fMinVtxType=1;
    else fMinVtxType=3;
  }
  // Set minimum vertex contributors
  void SetMinVtxContr(Int_t contr=1) { fMinVtxContr=contr; }
  // Set maximum z-vtx cut
  void SetMaxVtxZ(Float_t z=1e6) { fMaxVtxZ=z; }
  // Appky or not cut on difference SPD-TPC vtx (0: no cut, 1: |zvtx-SPD - zvtx-TPC|<0.5cm)
  void SetCutOnZVertexSPD(Int_t cut) { fCutOnZVertexSPD=cut; }

  // Select event trigger mask
  void SetTriggerMask(ULong64_t mask=0) { fTriggerMask=mask; }
  UInt_t GetTriggerMask(){ return fTriggerMask; }

  // Set minimum number of ITS, TPC, TOF or MUON clusters
  void SetNumberOfClusters(Int_t nITS, Int_t nTPC, Int_t nTOF, Int_t nMUON){
    fnClusITS = nITS; fnClusTPC = nTPC; fnClusTOF = nTOF; fnClusMUON = nMUON;
  }

  // PID setters (flag, particle specie, detector, limits)
  void SetUsePid(Bool_t flag=kTRUE) { fusePid=flag; }
  void SetParticleSpecie(AliPID::EParticleType type=AliPID::kPion) { fParticlePid=type; }
  void SetUseTPCPid(Bool_t flag=kTRUE) { fuseTPCPid=flag; }
  void SetUseTOFPid(Bool_t flag=kTRUE) { fuseTOFPid=flag; }
  // set given number of sigma cut per P bin
  void SetTPCSigmaPtBins(Int_t nPtBins, Float_t *pBinLimits, Float_t *sigmaBin);
  void SetTOFSigmaPtBins(Int_t nPtBins, Float_t *pBinLimits, Float_t *sigmaBin);
  // maximum momentum to use PID
  void SetMaximumPTPC(Float_t p) { fPmaxTPC=p; }
  void SetMaximumPTOF(Float_t p) { fPmaxTOF=p; }
  // Bayesian PID settings (0=no, 1=maximum)
  void SetUseCombinPID(Int_t flag=0){ fuseCombinPid=flag; }
  void SetPIDThreshold(Float_t value=0.3){ fThreshold=value; }
  //
  // Getters
  //
  Bool_t GetUsePid() const{ return fusePid; }
  Int_t  GetParticleSpecie() const { return fParticlePid; }
  Float_t *GetPTPCBinLimits() const { return fPTPCBinLimits; }
  Int_t  GetNPTPCBins() const {return fnPTPCBins; }
  Float_t *GetPTOFBinLimits() const { return fPTOFBinLimits; }
  Int_t  GetNPTOFBins() const { return fnPTOFBins; }
  Int_t GetUseCombinPID(){ return fuseCombinPid; }
  Float_t GetPIDThreshold(){ return fThreshold; }

 protected:

  Bool_t IsVertexSelected(AliVEvent *event);
  Bool_t CheckTPCPIDStatus(AliAODTrack *track) const;
  Bool_t CheckTOFPIDStatus(AliAODTrack *track) const;

  Bool_t fisAOD;  // flag wether it is AOD:1 or ESD:0 analysis

  Bool_t fIsPdgCode; // flag to check pdg code
  Int_t fPdgCode;    // particle pdg code
  Float_t fMaxProdRadius;  // maximum radius (cm) for primary selection

  Float_t fEtaMin;   // minimum eta cut
  Float_t fEtaMax;   // maximum eta cut
  Float_t fYMin;     // minimum Y cut
  Float_t fYMax;     // maximum Y cut
  Float_t fPtMin;    // minimum Pt cut
  Float_t fPtMax;    // maximum Pt cut
  Bool_t  fIsCharged; // check if particle is charged (MC level)

  UInt_t  fTriggerMask;   // event trigger mask
  Int_t   fMinVtxType;    // 0: not cut; 1: SPDZ; 2: SPD3D; 3: Tracks
  Int_t   fMinVtxContr;   // minimum vertex contributors
  Float_t fMaxVtxZ;       // maximum |z| of primary vertex
  Int_t fCutOnZVertexSPD; // 0: no cut, 1: |zvtx-SPD - zvtx-TPC|<0.5cm

  Int_t fnClusITS;   // minimum number of ITS clusters
  Int_t fnClusTPC;   // minimum number of TPC clusters
  Int_t fnClusTOF;   // minimum number of TOF clusters
  Int_t fnClusMUON;  // minimum number of MUON clusters

  Bool_t    fusePid;          // flag to use or not Pid
  Int_t     fParticlePid;     // integer to define the particle specie to check
  //
  Bool_t    fuseTPCPid;       // flag to use TPC Pid
  Int_t     fnPTPCBins;       // "number of limits", that is fnPBins+1
  Int_t     fnPTPCBinLimits;  // "number of limits", that is fnPBins+1
  Float_t*  fPTPCBinLimits;   //[fnPTPCBinLimits]  p bins
  Float_t*  fnSigmaTPC;       //[fnPTPCBins]
  Float_t   fPmaxTPC;         // maximum TPC P to use Pid
  //
  Bool_t    fuseTOFPid;       // flag to use TOF Pid
  Int_t     fnPTOFBins;       // "number of limits", that is fnPBins+1
  Int_t     fnPTOFBinLimits;  // "number of limits", that is fnPBins+1
  Float_t*  fPTOFBinLimits;   //[fnPTOFBinLimits]  p bins
  Float_t*  fnSigmaTOF;       //[fnPTOFBins]
  Float_t   fPmaxTOF;         // maximum TOF P to use Pid
  Int_t 	fuseCombinPid;	  // Bayesiand PID (0=no, 1=maximum, 2=Threshold)
  Float_t	fThreshold;		  // Threshold for PID Combined

  ClassDef(AliSingleTrackEffCuts,3)  // base class for cuts on single tracks
 };

#endif
