#ifndef AliAnalysisTaskAntipd_cxx
#define AliAnalysisTaskAntipd_cxx

// example of analysis task

class TH1F;
class TH1I;
class TH1;
class AliVEvent;
class AliAODEvent;
// class AliESDtrackCuts;
class TProfile;
class TList;
class AliPIDResponse; // main class for PID analysis

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskAntipd : public AliAnalysisTaskSE {
public:
	AliAnalysisTaskAntipd();
	AliAnalysisTaskAntipd(const char *name);
	virtual ~AliAnalysisTaskAntipd();

	// Double_t Psi_pair(AliESDtrack *neg, AliESDtrack *pos);
    void   UserCreateOutputObjects();
    void   UserExec(Option_t *option);
    void   Terminate(Option_t *);

    // void     SetTrackCuts(AliESDtrackCuts* const cuts) { fTrackCuts = cuts; }
    void     SetFilterBit(Int_t fBit)            {fFilterBit = fBit;}
    void     SetLowPCut(Float_t lowP)            {fLowPCut = lowP;}
    void     SetHighPCut(Float_t highP)          {fHighPCut = highP;}
    void     SetEtaCut(Float_t EtaCut)           {fEtaCut = EtaCut;}
    void     SetMinNITSCl(Int_t MinNclIts)       {fMinClIts = MinNclIts;}
    void     SetMaxDCAxy(Float_t maxDCAxy)       {fMaxDCAxy = maxDCAxy;}
    void     SetMaxDCAz(Float_t maxDCAz)         {fMaxDCAz = maxDCAz;}
    // PID
    void     SetMaxTPCnSigma(Float_t maxTPCnSigma)       {fMaxTPCnSigma = maxTPCnSigma;}
    void     SetMaxTOFnSigma(Float_t maxTOFnSigma)       {fMaxTOFnSigma = maxTOFnSigma;}
    void     SetUseTOFPidCut(Bool_t useTOFPidCut)        {fUseTOFPidCut = useTOFPidCut;}
    void     SetMomForTOFanaProt(Float_t momTOFprot)     {fMomTOFProt = momTOFprot;}
    void     SetMomForTOFanaDeut(Float_t momTOFdeut)     {fMomTOFDeut = momTOFdeut;}
    void     SetITSnSigmaRange(Float_t minITSnSigma, Float_t maxITSnSigma)
    {fMinITSnSigma = minITSnSigma; fMaxITSnSigma = maxITSnSigma;}

    void     CreateHistosTrack(std::vector<TH1*> &histos);
    void     FillHistosTrack(std::vector<TH1*> &histos, AliAODTrack *track);
    Float_t  GetTOFBeta(AliAODTrack *);
    Float_t  GetMass2TOF(Float_t, AliAODTrack *);

private:
    // transient members are not streamed
    AliVEvent   *fVEvent;           //! vEvent
    AliAODEvent	*fAODEvent;         //! AOD event
    TList       *fOutputList;       //! Output list
    TList       *fOutputEvent;      //! Event list
    TList       *fOutputProtons;    //! protons list
    TList       *fOutputAProtons;   //! anti-protons list
    TList       *fOutputDeuterons;  //! deuterons list
    TList       *fOutputADeuterons; //! anti-deuterons list
		TList       *fOutputHe3;        //! He3 list
		TList       *fOutputAHe3;        //! anti-He3 list

    TProfile    *fHistTrackCuts;    //! track cuts config

    // Event histograms
    TH1I        *fEventStat;       //! Event statistics
    TH1F        *fVertexZ;         //! Vertex Z coordinate
    TH1F        *fMultPercentileV0M; // all percentile range200 measuring multiplicity
    TH1F        *fMultPercentileV0MZoomed; //zoomed in 0.1% where the multiplicity
    // track histograms
    std::vector<TH1*> fHistsProton; // vector of proton histograms
    std::vector<TH1*> fHistsAProton; // vector of anti-proton histograms

    std::vector<TH1*> fHistsDeuteron; // vector of deuteron histograms
    std::vector<TH1*> fHistsADeuteron; // vector of anti-deuteron histograms

    std::vector<TH1*> fHistsHe3;      // vector of He3 histograms
		std::vector<TH1*> fHistsAHe3;      // vector of anti-He3 histograms

    AliPIDResponse *fPIDResponse;   //! PID response object

    // persistent members are streamed (copied/stored)
    // AliESDtrackCuts *fTrackCuts; // Track cuts
    Int_t         fFilterBit;
    Float_t       fLowPCut;
    Float_t       fHighPCut;
    Float_t       fEtaCut;
    Int_t         fMinClIts;
    Float_t       fMaxDCAxy;
    Float_t       fMaxDCAz;
    // PID
    Float_t       fMaxTPCnSigma;
    Float_t       fMaxTOFnSigma;
    Bool_t        fUseTOFPidCut;
    Float_t       fMinITSnSigma;
    Float_t       fMaxITSnSigma;
    Float_t       fMomTOFProt;
    Float_t       fMomTOFDeut;

    enum {kSelectedEvents=0, kINT7selected, kDAQincomplete, kV0timing, kClusterTrackletCut, kVertexZ, kVtxNcontrib, kSPDPileUp, kVtxDisplace, kVtxRes, kNbinsEvent};

	AliAnalysisTaskAntipd(const AliAnalysisTaskAntipd&); // not implemented
	AliAnalysisTaskAntipd& operator=(const AliAnalysisTaskAntipd&); // not implemented

	ClassDef(AliAnalysisTaskAntipd, 1);
};

#endif
