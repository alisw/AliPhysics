#ifndef ALIPRIMARYPIONCUTS_H
#define ALIPRIMARYPIONCUTS_H

// Class handling all kinds of selection cuts for primary

// Authors: Svein Lindal, Daniel Lohner												*


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
class TObjString;

using namespace std;

class AliPrimaryPionCuts : public AliAnalysisCuts {
		
	public: 


	enum cutIds {
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


	enum pionCuts {
		kPionIn=0,
		kNoTracks,
		kTrackCuts,
		kdEdxCuts,
		kPionOut
	};


	Bool_t SetCutIds(TString cutString); 
	Int_t fCuts[kNCuts];
	Bool_t SetCut(cutIds cutID, Int_t cut);
	Bool_t UpdateCutString();
	static const char * fgkCutNames[kNCuts];


	Bool_t InitializeCutsFromCutString(const TString analysisCutSelection); 
	

	AliPrimaryPionCuts(const char *name="PionCuts", const char * title="Pion Cuts");
	virtual ~AliPrimaryPionCuts();                            //virtual destructor

	virtual Bool_t IsSelected(TObject* /*obj*/){return kTRUE;}
	virtual Bool_t IsSelected(TList* /*list*/) {return kTRUE;}

	TString GetCutNumber();

		// Cut Selection
    Bool_t PionIsSelectedMC(Int_t labelParticle,AliMCEvent *mcEvent);
    Bool_t PionIsSelectedAODMC(Int_t labelParticle, TClonesArray *AODMCTrackArray);
    Bool_t TrackIsSelected(AliESDtrack* lTrack);
    Bool_t TrackIsSelectedAOD(AliAODTrack* lTrack);
    Bool_t PionIsSelected(AliESDtrack* lTrack);
    Bool_t PionIsSelectedAOD(AliAODTrack* lTrack);
	static AliPrimaryPionCuts * GetStandardCuts2010PbPb();
	static AliPrimaryPionCuts * GetStandardCuts2010pp();
	void  SetHybridTrackCutsAODFiltering(Int_t runflag);
	void  SetPtDepDCACuts(Double_t pt);
	Bool_t IsDCACutAccepted(AliAODTrack* lTrack);
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

	Bool_t SetTPCdEdxCutPionLine(Int_t pidedxSigmaCut);
	Bool_t SetITSdEdxCutPionLine(Int_t ededxSigmaCut);
	Bool_t SetITSClusterCut(Int_t clsITSCut);
	Bool_t SetTPCClusterCut(Int_t clsTPCCut);
	Bool_t SetEtaCut(Int_t etaCut);
	Bool_t SetPtCut(Int_t ptCut);
	Bool_t SetDCACut(Int_t dcaCut);
	void SetEtaShift(Double_t etaShift){fEtaShift = etaShift;}
	Bool_t SetTOFPionPIDCut(Int_t TOFelectronPID);
	Bool_t SetMassCut(Int_t massCut);
	void SetPeriodName(TString periodName){fPeriodName = periodName;}
	Double_t GetMassCut(){return fMassCut;}
    Double_t GetMassCut_WithNDM(){return fMassCut_WithNDM;}
	void SetPrefilterRunFlag(Int_t runflag){fRunFlag = runflag;}
	
	// Request Flags
	Double_t GetEtaCut(){ return  fEtaCut;}
    Double_t GetNFindableClustersTPC(AliVTrack* lTrack);
	Bool_t   DoWeights(){return fDoWeights;}
	Bool_t 	 DoMassCut(){return fDoMassCut;}
    Bool_t 	 DoMassCut_WithNDM(){return fDoMassCut_WithNDM;}
	Bool_t 	 RequireVertexConstrain(){return fRequireVertexConstrain;}
	Bool_t 	 Use4VecForMass(){return fUse4VecForMass;}
	
	protected:

    TList           *fHistograms;
    Bool_t          fDoLightOutput;             ///< switch for running light output, kFALSE -> normal mode, kTRUE -> light mode
    AliPIDResponse  *fPIDResponse;
	AliESDtrackCuts *fEsdTrackCuts;
	AliESDtrackCuts *fEsdTrackCutsGC;

	Double_t fEtaCut; //eta cutç
	Double_t fEtaShift;
	Bool_t   fDoEtaCut;
	Double_t fPtCut;
	Double_t fMinClsTPC; // minimum clusters in the TPC
    Double_t fChi2PerClsTPC; // maximum Chi2 per cluster in the TPC
    Bool_t   fRequireTPCRefit; // require a refit in the TPC
	Double_t fMinClsTPCToF; // minimum clusters to findable clusters
    Double_t fMaxSharedClsTPCFrac; // maximum fraction of shared clusters to TPCnClus
	Double_t fMinClsITS; // minimum clustersin the ITS
	Bool_t   fDodEdxSigmaITSCut; // flag to use the dEdxCut ITS based on sigmas
	Bool_t   fDodEdxSigmaTPCCut; // flag to use the dEdxCut TPC based on sigmas
	Bool_t   fDoTOFsigmaCut; // flag to use TOF pid cut RRnewTOF
	Double_t fPIDnSigmaAbovePionLineITS; // sigma cut ITS
	Double_t fPIDnSigmaBelowPionLineITS; // sigma cut ITS
	Double_t fPIDnSigmaAbovePionLineTPC; // sigma cut TPC
	Double_t fPIDnSigmaBelowPionLineTPC; // sigma cut TPC
	Double_t fPIDnSigmaAbovePionLineTOF; // sigma cut TOF
	Double_t fPIDnSigmaBelowPionLineTOF; // sigma cut TOF 
	Bool_t   fUseCorrectedTPCClsInfo; // flag to use corrected tpc cl info
	Bool_t   fUseTOFpid; // flag to use tof pid
	Bool_t   fRequireTOF; //flg to analyze only tracks with TOF signal
	Bool_t   fDoMassCut;
    Bool_t   fDoMassCut_WithNDM;
	Double_t fMassCut;	
    Double_t fMassCut_WithNDM;
	Bool_t fUse4VecForMass; // use only momentum 4vector to calculate inv mass
	Bool_t fRequireVertexConstrain; // require contrain to primary vertex (only for AOD)
	Bool_t   fDoWeights;
    Double_t fMaxDCAToVertexZ;
    Double_t fMaxDCAToVertexXY;
    Bool_t fUsePtDepXYDCA;
    Bool_t fUseDCAToVertex2D;
    TString fMaxDCAToVertexXYPtDep;
	Int_t  fRunFlag; // runflag used to set track prefiltering
	


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
	TH1F *fHistTrackSelectedEta;
	TH1F *fHistTrackSelectedPhi;
	TH1F *fHistTrackSelectedPt;
	TH1F *fHistTrackSelectedPtWithoutITS;
	
	TString fStringITSClusterCut;
	TString fPeriodName;
	
	private:

	AliPrimaryPionCuts(const AliPrimaryPionCuts&); // not implemented
	AliPrimaryPionCuts& operator=(const AliPrimaryPionCuts&); // not implemented


    ClassDef(AliPrimaryPionCuts,13)
};

#endif
