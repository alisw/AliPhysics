// ----------------------------------------------------------------------------
// AliAnalysisTaskDiHadronPID.h
// ----------------------------------------------------------------------------
// Author: Misha Veldhoen (misha.veldhoen@cern.ch)
// Last Rev.: May 2nd 2012. (v 8.00)
// ----------------------------------------------------------------------------

#ifndef ALIANALYSISTASKDIHADRONPID_H
#define ALIANALYSISTASKDIHADRONPID_H

#include <iostream>
#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"
#include "TMath.h"

using namespace std;

class TH1F;
class TH2F;
class TH3F;
class TList;
class TObjArray;
class TClonesArray;
class TString;

class AliAODTrack;
class AliAODEvent;
class AliAODVertex;

class AliPIDResponse;

class AliAnalysisTaskDiHadronPID: public AliAnalysisTaskSE {

public:
	// Required functions.
	AliAnalysisTaskDiHadronPID();
	AliAnalysisTaskDiHadronPID(const char *name);
	virtual ~AliAnalysisTaskDiHadronPID();
	
	virtual void UserCreateOutputObjects();
	virtual void UserExec(Option_t *option);
	virtual void Terminate(Option_t *);
    
    // Setters
    void SetVerbose(Int_t verbose) {fVerbose=verbose;}
    void SetPrintBufferSize(Bool_t printbuffersize=kTRUE) {fPrintBufferSize=printbuffersize;}
    void SetCalculateMixedEvents(Bool_t mixedevents=kTRUE) {fCalculateMixedEvents=mixedevents;}
    void SetMC(Bool_t mc=kTRUE) {fMC = mc;}
    void SetBeamType(TString beamtype) {
		if ((beamtype!="pp")&&(beamtype!="PbPb")) {
			cout<<"SetBeamType -> Beamtype must be pp or PbPb"<<endl;
			return;
		}
		fBeamType=beamtype;
	}
    
	void SetMaxEta(Double_t maxeta) {
		if (TMath::Abs(maxeta)>0.9) {
			cout<<"SetMaxEta -> |eta| must be < 0.9"<<endl;
			return;
		}
        fMaxEta = maxeta;
	}
    
    void SetMaxRapidityInInclusiveSpectra(Double_t maxrap) {fMaxRap=maxrap;}
        
    void SetMaxPlotEta(Double_t maxploteta) {
		if (TMath::Abs(maxploteta)>1.0) {
			cout<<"SetMaxPlotEta -> |eta| must be < 1.0"<<endl;
			return;
		}
        fMaxPlotEta = maxploteta;
	}
    
    void SetMaxPt(Double_t maxpt) {
        if (maxpt<5.) {
            cout<<"SetMaxPt -> Maximum pT must be > 5.0 GeV/c."<<endl;
            return;
        }
        fMaxPt = maxpt;
    }
    
    void SetNEtaBins(Int_t netabins) {
        if (netabins<1||netabins>72) {
            cout<<"SetNEtaBins -> Number of bins must be between 1 and 72"<<endl;
            return;
        }
        fNEtaBins = netabins;
    }
    
    void SetNPhiBins(Int_t nphibins) {
        if (nphibins<1||nphibins>72) {
            cout<<"SetNPhiBins -> Number of bins must be between 1 and 72"<<endl;
            return;
        }
        fNPhiBins = nphibins;

    }
    
    void SetVertexZMixedEvents(Double_t vertexzmixedevents) {
        if (vertexzmixedevents<0.||vertexzmixedevents>10.) {
            cout<<"SetVertexZMixedEvents -> must be 0 < z < 10"<<endl;
            return;
        } 
        fVertexZMixedEvents=vertexzmixedevents;
    } 
    
    void SetZoomed(Bool_t zoomed=kTRUE) {fZoomed=zoomed;}
	void SetDoDCACut(Bool_t dodcacut=kTRUE) {fDoDCACut=dodcacut;}
    void SetDoITSCut(Bool_t doitscut=kTRUE) {fDoITSCut=doitscut;}
    void SetDemandNoMismatch(Bool_t demandnomismatch=kTRUE) {fDemandNoMismatch=demandnomismatch;}
    void SetTrigBufferMaxSize(Int_t trigbuffermaxsize) {
        if (trigbuffermaxsize<10||trigbuffermaxsize>25000) {
            cout<<"SetTrigBufferMaxSize -> Max buffer size must be between 10 and 25000."<<endl;
            return;
        }
        fTrigBufferMaxSize=trigbuffermaxsize;
    }
    
    void SetCentralityCut(Double_t centralitycutmax, Double_t centralitycutmin) {
        if (centralitycutmax<0.) {
            cout<<"SetCentralityCut -> Centrality cannot be lower than 0."<<endl;
            return;
        }
        if (centralitycutmin<centralitycutmax) {
            cout<<"SetCentralityCut -> Maximum centrality needs to be smaller than the minimum centrality. (It's confusing I know)"<<endl;
            return;
        }
        if (centralitycutmin>100.) {
            cout<<"SetCentralityCut -> Minimum centrality cannot exceed 100%."<<endl;
            return;
        }
        fCentralityCutMax=centralitycutmax;
        fCentralityCutMin=centralitycutmin;
    
    }
	
    // Getters
    Bool_t GetVerbose() {return fVerbose;}
    Bool_t GetPrintBufferSize() {return fPrintBufferSize;}
    Bool_t GetCalculateMixedEvents() {return fCalculateMixedEvents;}
    TString GetBeamType() {return fBeamType;}
    Double_t GetMaxEta() {return fMaxEta;}
    Double_t GetMaxPlotEta() {return fMaxPlotEta;}
    Double_t GetMaxPt() {return fMaxPt;}
    Int_t GetNEtaBins() {return fNEtaBins;}
    Int_t GetNPhiBins() {return fNPhiBins;}
    Double_t GetVertexZMixedEvents() {return fVertexZMixedEvents;}
    Bool_t GetZoomed() {return fZoomed;}
    Bool_t GetDoDCACut() {return fDoDCACut;}
    Bool_t GetDoITSCut() {return fDoITSCut;}
    Bool_t GetDemandNoMismatch() {return fDemandNoMismatch;}
	Int_t GetTrigBufferMaxSize() {return fTrigBufferMaxSize;}
    Double_t GetCentralityCutMax() {return fCentralityCutMax;}
    Double_t GetCentralityCutMin() {return fCentralityCutMin;}
    
private:
	// Private Functions.	
	AliAnalysisTaskDiHadronPID(const AliAnalysisTaskDiHadronPID&); // NOT IMPLEMENTED.
	AliAnalysisTaskDiHadronPID& operator=(const AliAnalysisTaskDiHadronPID&); // NOT IMPLEMENTED.

    void FillGlobalTracksArray();
	AliAODTrack* GetGlobalTrack(AliAODTrack* track);
    
	Bool_t SelectEvent(AliAODVertex *vertex);
    Int_t ClassifyTrack(AliAODTrack* track);

	Double_t PhiRange(Double_t DPhi);
    Int_t ConvertPdgCode(Int_t pdgcode);
	
private:
	// PID object.
	AliPIDResponse		*fPIDResponse;              //! PID Response Handler.
	
	// Event and Track related objects.
	AliAODEvent			*fAODEvent;                 //! The AOD Event.
	AliAODHeader		*fAODHeader;                //! The AOD Header.
	AliAODVertex		*fAODVertex;                //! The AOD Vertex.
	
	AliAODTrack			*fAODTrack;                 //! Current AOD Track.

	TObjArray			*fGlobalTracks;             //! Partner Tracks.
	TClonesArray		*fMCTracks;					//! MC tracks, indexed by their track label.
	
	// HISTOGRAMS.

	// Event QA plots.
	TH1F			    *fCentrality;               //! Centrality Histogram.
	TH1F			    *fVertexZ;                  //! Vertex Z position.

    // Track QA plots.
    TH2F                *fDCA;                      //! DCA XY vs Z before DCA cut.
    TH2F                *fDCAZoomed;                //!
    TH2F                *fDCAZoomedTwice;           //!
    TH2F                *fDCACut;                   //! DCA XY vs Z after DCA cut (if performed!).
    TH2F                *fDCAZoomedCut;             //!
    TH2F                *fDCAZoomedTwiceCut;        //!
    
    TH1F                *fITSHits;                  //! 3 bins, [no hits in first 2 layers, 1 hit, 2 hits]
    
    TH1F                *fTrackCutsCount;           //! Counts of used tracks after cuts
    TH2F                *fTrackCutsPt;              //! pT spectrum after cuts.
    TH2F                *fTrackCutsEta;             //! eta spectrum after cuts.
    TH2F                *fTrackCutsPhi;             //! phi spectrum after cuts.
    
    TH1F                *fEtaSpectrumTrig;          //! eta spectrum of triggers (pT > 5.0 tracks, trigger track cuts.)
    TH2F                *fEtaSpectrumAssoc;         //! eta spectrum of associateds as a function of pT
    TH2F                *fPhiSpectrumAssoc;         //! phi spectrum of associateds as a function of pT
    
	// PID QA plots.
	TH2F				*fTPCnSigmaProton;          //! TPC nSigma plot for Protons.
	TH2F				*fTOFnSigmaProton;          //! TOF nSigma plot for Protons.
	TH2F				*fTPCnSigmaPion;            //! TPC nSigma plot for Pions.
	TH2F				*fTOFnSigmaPion;            //! TOF nSigma plot for Pions.	
	TH2F				*fTPCnSigmaKaon;            //! TPC nSigma plot for Kaons.
	TH2F				*fTOFnSigmaKaon;            //! TOF nSigma plot for Kaons.
	
	TH3F				*fTPCSignal;				//! TPC signal (pt,eta).
	TH3F				*fTOFSignal;				//! TOF signal (pt,eta).
    TH3F                *fInclusiveTPCTOF[3][10];   //! inclusive TPC-TOF histogram as a function of pT (and eta)
    TH3F                *fInclusiveTPCTOFRap[3][10];//! inclusive TPC-TOF histogram as a function of pT and rapidity, with additional rapidity cut.
    
	// Efficiency Plots (Monte Carlo)
	TH2F				*fPtEtaDistrDataPrim[6];	//! pT distribution of physical primaries [species: pi+,pi-,K+,K-,p,pbar].
	TH2F				*fPtEtaDistrDataSec[6];		//! 

	TH2F				*fPtRapDistrDataPrimRapCut[6];//! pT distribution of physical primaries [species: pi+,pi-,K+,K-,p,pbar].
	TH2F				*fPtRapDistrDataSecRapCut[6];//! with an additional rapidity cut.
	    
	TH2F				*fPtEtaDistrMCPrim[6];		//! pT distribution of MCParticles. [species: pi+,pi-,K+,K-,p,pbar].
	TH2F				*fPtEtaDistrMCSec[6];		//! 

    TH2F				*fPtRapDistrMCPrimRapCut[6];//! pT distribution of MCParticles. [species: pi+,pi-,K+,K-,p,pbar].
	TH2F				*fPtRapDistrMCSecRapCut[6];	//! with an additional rapidity cut.
    
	TH3F				*fDiHadronMC[6];			//! DPhiDEta Plot per species [species]
	
	// Di-Hadron Correlations.
	TH3F				*fDiHadron;					//! regular di-hadron correlation, accepting all associateds.
	THnSparseF			*fDiHadronTPCTOF[3][10];	//! Di-Hadron correlations with both TPC and TOF signal.

    // Mixed Events.
    TH3F				*fMixedEvents;              //! Mixed Events, associated track cuts.
    //TH3F 				*fMixedEventsTPCTOFCut[3];	//! For every species seperately we keep track of mixed events.
	THnSparseF			*fMixedEventsTPCTOF[3][10];	//!
	
	// List of Histograms.
	TList			    *fHistoList;                //! List of Histograms.
	
	// Analysis Task Configuration Variables.
    Bool_t              fCalculateMixedEvents;      // 
	TString				fBeamType;                  // pp or PbPb
    Bool_t              fMC;                        // runs over MC.
	Double_t			fMaxEta;                    // Q: Do we need to take extra care of the binning?
    Double_t            fMaxPlotEta;                //
    Double_t            fMaxRap;                    // Max rapidity, applied to the inclusive spectra.
    Double_t            fMaxPt;                     //
    Int_t               fNEtaBins;                  // Number of bins in eta
    Int_t               fNPhiBins;                  // Number of bins in phi
    Double_t            fVertexZMixedEvents;        // Events with a vertex z difference smaller than 
                                                    // this number (standard 2cm) will be mixed.
    
    Double_t            fCentralityCutMax;          // Maximum centrality (standard 0%)
    Double_t            fCentralityCutMin;          // Minimum centrality (standard 10%)
    Bool_t              fZoomed;                    //
    Bool_t              fDoITSCut;                  // Cut the tracks with not at least one SPD hit.
	Bool_t              fDoDCACut;                  // Perform a DCA cut to get rid of secondaries.
    Bool_t              fDemandNoMismatch;          // 
    
    Int_t               fTrackCutLabelNumbers[8];   // Track Cut labels.
    
    // Level of verbal output.
    //  0 -> No output.
    //  1 -> Only error messages.
    //  2 -> Information about output creation (beginning of the job)
    //  3 -> Event information.
    //  4 -> Track information.
    Int_t               fVerbose;                   //
    Bool_t              fPrintBufferSize;           //

    // Trigger buffer.
    Double_t            fTrigBuffer[25000][4];      //!
	Int_t               fTrigBufferIndex;           //!
    Int_t               fTrigBufferSize;            //!
	Int_t				fTrigBufferMaxSize;			//!
    
    
	ClassDef(AliAnalysisTaskDiHadronPID,1);
	
};

#endif

