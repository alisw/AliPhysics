// ----------------------------------------------------------------------------
// Header file.
// ----------------------------------------------------------------------------
// Last Rev.: 17th of February 2012. (v08)
// ----------------------------------------------------------------------------

#ifndef ALIANALYSISTASKDIHADRONPID_H
#define ALIANALYSISTASKDIHADRONPID_H

#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"

class TH1F;
class TH2F;
class TH3F;
//class THnSparseF;
class TList;
class TObjArray;
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
    void SetVerbose(Bool_t verbose=kTRUE) {fVerbose=verbose;}
    void SetCalculateMixedEvents(Bool_t mixedevents=kTRUE) {fCalculateMixedEvents=mixedevents;}
    
    // Getters
    Bool_t GetVerbose() {return fVerbose;}
    Bool_t GetCalculateMixedEvents() {return fCalculateMixedEvents;}
    
private:
	// Private Functions.	
	AliAnalysisTaskDiHadronPID(const AliAnalysisTaskDiHadronPID&); // NOT IMPLEMENTED.
	AliAnalysisTaskDiHadronPID& operator=(const AliAnalysisTaskDiHadronPID&); // NOT IMPLEMENTED.
	
	Bool_t SelectTrack(AliAODTrack *track, Int_t cuts=0);
	Bool_t SelectEvent(AliAODVertex *vertex);
	
	void FillPIDPartnersArray();
	AliAODTrack* GetTrackPartner(AliAODTrack* track);

	Double_t PhiRange(Double_t DPhi);
	
private:
	// PID object.
	AliPIDResponse		*fPIDResponse;              //! PID Response Handler.
	
	// Event and Track related objects.
	AliAODEvent			*fAODEvent;                 //! The AOD Event.
	AliAODHeader		*fAODHeader;                //! The AOD Header.
	AliAODVertex		*fAODVertex;                //! The AOD Vertex.
	
	AliAODTrack			*fAODTrack;                 //! Current AOD Track.

	TObjArray			*fPIDPartners;              //! Partner Tracks.
	
	// HISTOGRAMS.

	// Event Sample Plots.
	TH1F			    *fCentrality;               //! Centrality Histogram.
	TH1F			    *fVertexZ;                  //! Vertex Z position.

    // Track cuts counts.
    TH1F                *fTrackCuts;                //! Different Track cuts.
    
	// Unidentified Spectra.
	TH1F			    *fPtSpectrum;               //! Unidentified Pt-spectrum (after standard track cuts).
    TH3F                *fAssociatedDistribution;   //! Associated distribution. 

	// QA plots PID.
	TH2F				*fTPCnSigmaProton;          //! TPC nSigma plot for Protons.
	TH2F				*fTOFnSigmaProton;          //! TOF nSigma plot for Protons.
	TH2F				*fTPCnSigmaPion;            //! TPC nSigma plot for Pions.
	TH2F				*fTOFnSigmaPion;            //! TOF nSigma plot for Pions.	
	TH2F				*fTPCnSigmaKaon;            //! TPC nSigma plot for Kaons.
	TH2F				*fTOFnSigmaKaon;            //! TOF nSigma plot for Kaons.
	
    // PID signals as function of pT and Eta
	TH3F				*fTPCSignal;				//! TPC signal (pt,eta).
	TH3F				*fTOFSignal;				//! TOF signal (pT,eta).
    
	// Unidentified Di-Hadron Correlations & Mixed Events.
	TH3F				*fDiHadron;                 //! Di-Hadron correlations.	
	TH3F				*fMixedEvents;              //! Mixed Events.
	
	// Di-Hadron Correlations with TPC and TOF signals.
    TH3F                *fDiHadronTPC[3][10];       //! Di-Hadron correlations with TPC signal.
    TH3F                *fDiHadronTOF[3][10];       //! Di-Hadron correlations with TOF signal.
	THnSparseF			*fDiHadronTPCTOF[3][10];	//! Di-Hadron correlations with both TPC and TOF signal.

	// List of Histograms.
	TList			    *fHistoList;                //! List of Histograms.
	
	// Other data members.
	Bool_t				fVerbose;                   // Generates verbal output if kTRUE.
	UInt_t				fMask;                      // PID efficiency graphs will be made for this mask.
    Bool_t              fCalculateMixedEvents;      // 
    
    // Trigger buffer.
    Double_t            fTrigBuffer[25000][4];      //!
	Int_t               fTrigBufferIndex;           //!
    Int_t               fTrigBufferSize;            //!
    
    
	ClassDef(AliAnalysisTaskDiHadronPID,1);
	
};

#endif

