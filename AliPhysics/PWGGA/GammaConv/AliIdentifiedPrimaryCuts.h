#ifndef ALIIDENTIFIEDPRIMARYCUTS_H
#define ALIIDENTIFIEDPRIMARYCUTS_H

// Class handling all kinds of selection cuts for primary

// Authors: Svein Lindal, Daniel Lohner												*
// Authors: Ana Marin for Identified charged

#include "AliAODpidUtil.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliMCEvent.h"
#include "AliAnalysisCuts.h"
#include "AliESDtrackCuts.h"
#include "TH1F.h"

class AliESDEvent;
class AliAODEvent;
class AliConversionPhotonBase;
class AliKFVertex;
class AliKFParticle;
class TH1F;
class TH2F;
class AliPIDResponse;
class AliAnalysisCuts;
class iostream;
class TList;
class AliAnalysisManager;


using namespace std;

class AliIdentifiedPrimaryCuts : public AliAnalysisCuts {
		
 public: 
  
  
  enum cutIds {
    kParticleTypeCut,
    kEtaCut,
    kClsITSCut,
    kClsTPCCut,
    kDCACut,
    kPtCut,
    kPidedxSigmaITSCut,
    kPidedxSigmaTPCCut,
    kPiTOFSigmaPID,
    kMassCut,
    kNCuts
  };
  
  
  enum identifiedCuts {
    kIdentifiedIn=0,
    kNoTracks,
    kTrackCuts,
    kdEdxCuts,
    kIdentifiedOut
  };
  
  
  Bool_t SetCutIds(TString cutString); 
  Int_t fCuts[kNCuts];
  Bool_t SetCut(cutIds cutID, Int_t cut);
  Bool_t UpdateCutString();
  static const char * fgkCutNames[kNCuts];
  
  
  Bool_t InitializeCutsFromCutString(const TString analysisCutSelection); 
  
  
  AliIdentifiedPrimaryCuts(const char *name="IdentifiedCuts", const char * title="Identified Cuts");
  virtual ~AliIdentifiedPrimaryCuts();                            //virtual destructor
  
  virtual Bool_t IsSelected(TObject* /*obj*/){return kTRUE;}
  virtual Bool_t IsSelected(TList* /*list*/) {return kTRUE;}
  
  TString GetCutNumber();
  
  // Cut Selection
  Bool_t IdentifiedIsSelectedMC(Int_t labelParticle,AliMCEvent *mcEvent);
  Bool_t IdentifiedIsSelectedAODMC(Int_t labelParticle, TClonesArray *AODMCTrackArray);
  Bool_t TrackIsSelected(AliESDtrack* lTrack);
  Bool_t TrackIsSelectedAOD(AliAODTrack* lTrack);
  Bool_t IdentifiedIsSelected(AliESDtrack* lTrack);
  Bool_t IdentifiedIsSelectedAOD(AliAODTrack* lTrack);
  static AliIdentifiedPrimaryCuts * GetStandardCuts2010PbPb();
  static AliIdentifiedPrimaryCuts * GetStandardCuts2010pp();
  Bool_t InitPIDResponse();
  
  void SetPIDResponse(AliPIDResponse * pidResponse) {fPIDResponse = pidResponse;}
  AliPIDResponse * GetPIDResponse() { return fPIDResponse;}
  
  void PrintCuts();
  void PrintCutsWithValues();
  
  void    SetLightOutput( Bool_t flag ){fDoLightOutput = flag; return;}
  
  void InitCutHistograms(TString name="",Bool_t preCut = kTRUE,TString cutName="");
  void SetFillCutHistograms(TString name="",Bool_t preCut = kTRUE,TString cutName=""){if(!fHistograms){InitCutHistograms(name,preCut,cutName);};}
  TList *GetCutHistograms(){return fHistograms;}
  
  static AliVTrack * GetTrack(AliVEvent * event, Int_t label);
  
  ///Cut functions
  Bool_t dEdxCuts(AliVTrack * track);
  
  Bool_t SetTPCdEdxCutIdentifiedLine(Int_t pidedxSigmaCut);
  Bool_t SetITSdEdxCutIdentifiedLine(Int_t ededxSigmaCut);
  Bool_t SetITSClusterCut(Int_t clsITSCut);
  Bool_t SetTPCClusterCut(Int_t clsTPCCut);
  Bool_t SetParticleTypeCut(Int_t particleTypeCut);
  Bool_t SetEtaCut(Int_t etaCut);
  Bool_t SetPtCut(Int_t ptCut);
  Bool_t SetDCACut(Int_t dcaCut);
  void SetEtaShift(Double_t etaShift){fEtaShift = etaShift;}
  Bool_t SetTOFIdentifiedPIDCut(Int_t TOFelectronPID);
  Bool_t SetMassCut(Int_t massCut);
  Double_t GetMassCut(){return fMassCut;}
  
  // Request Flags
  Int_t GetParticleTypeCut(){ return  fParticleTypeCut;}
  Double_t GetEtaCut(){ return  fEtaCut;}
  Double_t GetNFindableClustersTPC(AliVTrack* lTrack);
  Bool_t   DoWeights(){return fDoWeights;}
  Bool_t   DoMassCut(){return fDoMassCut;}
  
 protected:
  
  TList           *fHistograms;
  Bool_t          fDoLightOutput;             ///< switch for running light output, kFALSE -> normal mode, kTRUE -> light mode
  AliPIDResponse  *fPIDResponse;
  AliESDtrackCuts *fEsdTrackCuts;
  
  Int_t    fParticleTypeCut;
  Double_t fEtaCut; //eta cutç
  Double_t fEtaShift;
  Bool_t   fDoEtaCut;
  Double_t fPtCut;
  Double_t fMinClsTPC; // minimum clusters in the TPC
  Double_t fChi2PerClsTPC; // maximum Chi2 per cluster in the TPC
  Bool_t   fRequireTPCRefit; // require a refit in the TPC
  Double_t fMinClsTPCToF; // minimum clusters to findable clusters
  Bool_t   fDodEdxSigmaITSCut; // flag to use the dEdxCut ITS based on sigmas
  Bool_t   fDodEdxSigmaTPCCut; // flag to use the dEdxCut TPC based on sigmas
  Bool_t   fDoTOFsigmaCut; // flag to use TOF pid cut RRnewTOF
  Double_t fPIDnSigmaAboveIdentifiedLineITS; // sigma cut ITS
  Double_t fPIDnSigmaBelowIdentifiedLineITS; // sigma cut ITS
  Double_t fPIDnSigmaAboveIdentifiedLineTPC; // sigma cut TPC
  Double_t fPIDnSigmaBelowIdentifiedLineTPC; // sigma cut TPC
  Double_t fPIDnSigmaAboveIdentifiedLineTOF; // sigma cut TOF
  Double_t fPIDnSigmaBelowIdentifiedLineTOF; // sigma cut TOF 
  Bool_t   fUseCorrectedTPCClsInfo; // flag to use corrected tpc cl info
  Bool_t   fUseTOFpid; // flag to use tof pid
  Bool_t   fRequireTOF; //flg to analyze only tracks with TOF signal
  Bool_t   fDoMassCut;
  Double_t fMassCut;	
  Bool_t   fDoWeights;
  Double_t fMaxDCAToVertexZ;
  
  
  
  // Histograms
  TObjString *fCutString; // cut number used for analysis
  TString fCutStringRead;
  TH1F *fHistCutIndex; // bookkeeping for cuts
  TH1F *fHistdEdxCuts;  // bookkeeping for dEdx cuts
  TH2F *fHistITSdEdxbefore; // ITS dEdx before cuts
  TH2F *fHistITSdEdxafter;
  TH2F *fHistTPCdEdxbefore; // TPC dEdx before cuts
  TH2F *fHistTPCdEdxafter; // TPC dEdx after cuts
  TH2F *fHistTPCdEdxSignalbefore; //TPC dEdx signal before
  TH2F *fHistTPCdEdxSignalafter; //TPC dEdx signal  after
  TH2F *fHistTOFbefore; // TOF after cuts
  TH2F *fHistTOFafter; // TOF after cuts
  TH2F *fHistTrackDCAxyPtbefore;
  TH2F *fHistTrackDCAxyPtafter;
  TH2F *fHistTrackDCAzPtbefore;
  TH2F *fHistTrackDCAzPtafter;
  TH2F *fHistTrackNFindClsPtTPCbefore;
  TH2F *fHistTrackNFindClsPtTPCafter;
  
  TString fStringITSClusterCut;
  
 private:
  
  AliIdentifiedPrimaryCuts(const AliIdentifiedPrimaryCuts&); // not implemented
  AliIdentifiedPrimaryCuts& operator=(const AliIdentifiedPrimaryCuts&); // not implemented
  
  
  ClassDef(AliIdentifiedPrimaryCuts,1)
    };

#endif
