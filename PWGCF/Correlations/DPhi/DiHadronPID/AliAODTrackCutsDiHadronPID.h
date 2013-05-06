#ifndef ALIAODTRACKCUTSDIHADRONPID_H
#define ALIAODTRACKCUTSDIHADRONPID_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. * 
* See cxx source for full Copyright notice */ 
/* $Id$ */

#include "TFormula.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"

class AliAODTrackCutsDiHadronPID : public TNamed 

{

public:
	enum HistoClass {kAllCharged = 0, kPositive = 1, kNegative = 2,
		kAllPion = 3, kPosPion = 4, kNegPion = 5, 
		kAllKaon = 6, kPosKaon = 7, kNegKaon = 8, 
		kAllProton = 9, kPosProton = 10, kNegProton = 11};

public:
	AliAODTrackCutsDiHadronPID();					// Default Constructor
	AliAODTrackCutsDiHadronPID(const char* name);	// Named Constructor
	virtual ~AliAODTrackCutsDiHadronPID();			// Destructor
	virtual Long64_t Merge(TCollection* list);		// Merger

private:
	AliAODTrackCutsDiHadronPID(const AliAODTrackCutsDiHadronPID&);
	AliAODTrackCutsDiHadronPID& operator=(const AliAODTrackCutsDiHadronPID&);

// -----------------------------------------------------------------------
//  Interface, methods used to get information about the track cuts, and to
//  retrieve filled histograms. 
// -----------------------------------------------------------------------

public:
	void PrintCuts() const;

	// List of QA histograms.
	TList* GetListOfDataQAHistos() const;
	TList* GetListOfPrimRecMCTrackQAHistos() const;
	TList* GetListOfPrimGenMCTrackQAHistos() const;
	TList* GetListOfSecRecMCTrackQAHistos() const;
	TList* GetListOfSecGenMCTrackQAHistos() const;

	// Return a specific member of one of the lists of histograms.
	TObject* GetHistData(const char* name) const {return fDataTrackQAHistos->FindObject(name);} 
	TObject* GetHistPrimRecMC(const char* name) const {return fPrimRecMCTrackQAHistos->FindObject(name);}
	TObject* GetHistPrimGenMC(const char* name) const {return fPrimGenMCTrackQAHistos->FindObject(name);}
	TObject* GetHistSecRecMC(const char* name) const {return fSecRecMCTrackQAHistos->FindObject(name);}
	TObject* GetHistSecGenMC(const char* name) const {return fSecGenMCTrackQAHistos->FindObject(name);}

	// Return a projection of signal or mismatch onto the TOF axis (FIXME: works only data, not yet MC)
	TH1F* GetHistDataTOFProjection(Int_t charge, Int_t species, Int_t ptbin);
	TObjArray* GetDataTOFProjection(Int_t charge, Int_t species);
 	TH1F* GetHistDataTOFMismatch(Int_t charge, Int_t species, Int_t ptbin);
 	TObjArray* GetDataTOFMismatch(Int_t charge, Int_t species);
	
	// Return a projection of signal or mismatch onto the TOF axis (FIXME: works only data, not yet MC)	
 	TH2F* GetHistDataTPCTOFProjection(Int_t charge, Int_t species, Int_t ptbin);
 	TObjArray* GetDataTPCTOFProjection(Int_t charge, Int_t species);
 	TH2F* GetHistDataTPCTOFMismatch(Int_t charge, Int_t species, Int_t ptbin);
 	TObjArray* GetDataTPCTOFMismatch(Int_t charge, Int_t species);

	// Note that there are two p_T axes, one for PID histograms, and one for other histograms.
	// Methods regarding the "other" p_T axis.
	Int_t GetNPtBins() const {return fNPtBins;}
	Double_t* GetPtAxis() {return fPtAxis;}
	Double_t GetPtMin(const Int_t bin) const;
	Double_t GetPtMax(const Int_t bin) const;
	Double_t GetPtBinWidth(Int_t bin) const {return (GetPtMax(bin) - GetPtMin(bin)); }

	// Methods regarding the PID p_T axis.
	Int_t GetNPtBinsPID(const Int_t ptclass = -1) const;
	Double_t* GetPtAxisPID() const;
	Double_t GetPtMinPID(const Int_t bin) const;
	Double_t GetPtMaxPID(const Int_t bin) const;
	Double_t GetPtBinWidthPID(Int_t bin) const {return (GetPtMaxPID(bin) - GetPtMinPID(bin));}
	Double_t GetPtClassMin(const Int_t ptclass) const;
	Double_t GetPtClassMax(const Int_t ptclass) const;

	Int_t GetNTOFbins(const Int_t ptclass, const Int_t species) const {return fTOFbins[ptclass][species];}
	Double_t GetTOFmin(const Int_t ptclass, const Int_t species) const {return fTOFLowerBound[ptclass][species];}
	Double_t GetTOFmax(const Int_t ptclass, const Int_t species) const {return fTOFUpperBound[ptclass][species];}

	Int_t GetNTPCbins(const Int_t ptclass, const Int_t species) const {return fTPCbins[ptclass][species];}
	Double_t GetTPCmin(const Int_t ptclass, const Int_t species) const {return fTPCLowerBound[ptclass][species];}
	Double_t GetTPCmax(const Int_t ptclass, const Int_t species) const {return fTPCUpperBound[ptclass][species];}

	// Getters (Cuts)
	UInt_t GetFilterMask() const {return fFilterMask;}
	Double_t GetMaxEta() const {return fMaxEta;}
	ULong_t GetDemandedFlags() const {return fDemandedFlags;}
	TFormula* GetPtDeptDCACutFormula() const {return fPtDeptDCAxyCutFormula;}
	Double_t GetDCAzCut() const {return fDCAzCut;}
	UInt_t GetMinSPDHitsForPtDeptDCACut() const {return fMinSPDHitsForPtDeptDCAcut;}

	// Getters (Settings)
	Bool_t GetIsMC() const {return fIsMC;}
	Int_t GetDebugLevel() const {return fDebug;}

// -----------------------------------------------------------------------
//  Methods used to configure the track cuts object, to be called at
//  initialization, i.e., before the object is added to an analysis task.  
// -----------------------------------------------------------------------

public:

	// Request Certain QA histograms being filled.
	Bool_t RequestQAHistos(const Int_t histoclass, const Bool_t Enable3DSpectra = kFALSE, const Bool_t EnablePIDHistos = kFALSE);
	
	// Setters (Cuts)
	void SetPtRange(Double_t minpt, Double_t maxpt) {
		fMinPt = minpt;
		fMaxPt = maxpt;
		fTestPt = kTRUE;
	}
	void SetFilterMask(UInt_t filtermask) {
		fFilterMask = filtermask;
		fTestFilterMask = kTRUE;
	}
	void SetMaxEta(Double_t maxeta) {
		fMaxEta = maxeta;
		fTestMaxEta = kTRUE;
	}
	void SetMaxRapidity(Double_t maxrapidity) {
		fMaxRapidity = maxrapidity;
		fTestMaxRapidity = kTRUE;
	}
	void SetDemandNoMismatch() {
		fTestTOFmismatch = kTRUE;
	}
	void SetDemandFlags(ULong_t demandedflags) {
		fDemandedFlags = demandedflags;
		fTestFlags = kTRUE;
	}
	void SetPtDeptDCACut(TFormula* DCAxyCutFormula, Double_t DCAzCut, UInt_t MinSPDHits = 1) {
		fPtDeptDCAxyCutFormula = DCAxyCutFormula;
		fDCAzCut = DCAzCut;
		fMinSPDHitsForPtDeptDCAcut = MinSPDHits;
		fTestPtDeptDCAcut = kTRUE;
	}

// Setters (Settings)
	void SetIsMC(Bool_t ismc = kTRUE) {fIsMC = ismc;}

	void SetDebugLevel(Int_t debuglevel) {fDebug = debuglevel;}

// -----------------------------------------------------------------------
//  Methods called by the analysis task. 
// -----------------------------------------------------------------------

public:

// These two functions signal the beginning and the end of a new event.
	void StartNewEvent();				// Some things are set to zero.
	void EventIsDone(Bool_t IsMC);		// Some final histograms are filled.
	void CreateHistos(); 				// Should be called by the UserCreateOutput() function of the analysis task.

// Is Selected, for different types of tracks.
	Bool_t IsSelectedData(AliTrackDiHadronPID* track, Double_t randomhittime = -1.e20);
	Bool_t IsSelectedGeneratedMC(AliAODMCParticle* particle);
	Bool_t IsSelectedReconstructedMC(AliTrackDiHadronPID* track);

// -----------------------------------------------------------------------
//  Internal methods.
// -----------------------------------------------------------------------

public:

// For PID histograms we have a certain number of bins in pT, spread out over five
// large histograms, i.e., one for the lowest pT, and the biggest range in TOF/TPC,
// one for the higher pT and smaller range in TOF/TPC, etc. The following methods
// are a mapping between the total pT bin (what the user uses), and the pt bin
// within one of the five histograms (what's used internally)
	Int_t GetPtClass(const Int_t ptbin) const;
	Int_t GetBinInPtClass(const Int_t ptbin) const;

private:

// Checks, return kTRUE if track passes the cut.
	Bool_t CheckPt(Double_t pt) const {
		// TODO: TO IMPLEMENTATION.		
		if (!fTestPt) return kTRUE;
		if ((pt > fMinPt) && (pt < fMaxPt)) return kTRUE;
		return kFALSE;
	}
	Bool_t CheckMaxEta(Double_t eta) const {
		if (!fTestMaxEta) return kTRUE;				// Accepted if there is no check on this parameter.
		if (TMath::Abs(eta) < fMaxEta) return kTRUE;
		return kFALSE;
	}
	Bool_t CheckRapidity(Double_t rap) const {
		if (!fTestMaxRapidity) return kTRUE;
		if (TMath::Abs(rap) < fMaxRapidity) return kTRUE;
		return kFALSE;
	}
	Bool_t CheckFilterMask(UInt_t filtermap) const {
		if (!fTestFilterMask) return kTRUE;
		if ((fFilterMask & filtermap) == fFilterMask) return kTRUE;
		return kFALSE;
	}
	Bool_t CheckFlags(ULong_t flags) const {
		if (!fTestFlags) return kTRUE;
		if ((flags & fDemandedFlags) == fDemandedFlags) return kTRUE;
		return kFALSE;
	}
	Bool_t CheckTOFmismatch(Bool_t ismismatch) const {
		if (!fTestTOFmismatch) return kTRUE; // if we're not cutting on mismatch, then it's accepted.
		if (!ismismatch) return kTRUE; 		// so if the track is not a mismatch, then it is accepted.
		return kFALSE; 						// if it is a mismatch, then it's not accepted.
	}
	Bool_t CheckPtDeptDCACut(Double_t dcaz, Double_t dcaxy, Double_t pt, UInt_t SPDhits) const {
		if (!fTestPtDeptDCAcut) return kTRUE;
		if (SPDhits < fMinSPDHitsForPtDeptDCAcut) return kTRUE; // If there are not enough SPD hits to do the cut.
		if ((dcaz < fDCAzCut) && (dcaxy < fPtDeptDCAxyCutFormula->Eval(pt))) return kTRUE;
		return kFALSE;
	}

// Filling QA histograms.
	Bool_t FillDataHistos(Int_t histoclass, AliTrackDiHadronPID* track);
	Bool_t FillTOFMismatchHistos(Int_t histoclass, AliTrackDiHadronPID* track, Double_t randomhittime);
	Bool_t FillGenMCHistos(Int_t histoclass, AliAODMCParticle* particle);
	Bool_t FillRecMCHistos(Int_t histoclass, AliTrackDiHadronPID* track);

// Initializing QA histograms.
	Bool_t InitializeDataHistos(Int_t histoclass);
	Bool_t InitializeGenMCHistos(Int_t histoclass);
	Bool_t InitializeRecMCHistos(Int_t histoclass);

	void InitializeDefaultHistoNamesAndAxes();

	TH1F* InitializePtSpectrum(const char* name, Int_t histoclass);
	TH3F* InitializePhiEtaPt(const char* name, Int_t histoclass);
	TH1F* InitializeNTracksHisto(const char* name, Int_t histoclass);
	TH1F* InitializeDCAxyHisto(const char* name, Int_t histoclass);
	TH1F* InitializeDCAzHisto(const char* name, Int_t histoclass);
	TH3F* InitializeAcceptanceHisto(const char* /*name*/, Int_t /*histoclass*/); // TO BE IMPLEMENTED.
	TH2F* InitializeDCASpectrum(const char* name, Int_t histoclass);


	TH3F* InitializePIDHisto(const char* name, Int_t histoclass, Int_t expspecies, Int_t ptclass);
	TH2F* InitializeTOFMismatchHisto(const char* name, Int_t histoclass, Int_t expspecies, Int_t ptclass);

// -----------------------------------------------------------------------
//  Data members.
// -----------------------------------------------------------------------

private:
// Track Cuts
	Double_t				fMinPt;
	Double_t				fMaxPt;
	UInt_t 					fFilterMask;					// FilterMask to-be-checked.
	Double_t 				fMaxEta;						// Max Eta of the track.
	Double_t				fMaxRapidity;					// Rapidity cut (only done for PID plots!!)
	ULong_t 				fDemandedFlags;					// Flags demanded on the track.
	UInt_t 					fMinSPDHitsForPtDeptDCAcut;		// Required number of SPD hits for performing Pt-Dept DCA cut.
	TFormula* 				fPtDeptDCAxyCutFormula;			// Formula for the Pt-Dept DCA cut.
	Double_t				fDCAzCut;						// Max z at DCA.

// Settings
	Bool_t					fIsMC;							// Is the current event MC or not.

// Requested Histograms;
	Bool_t					fHistRequested[12];				//
	Bool_t					f3DSpectraEnabeled[12];			//
	Bool_t					fPIDHistosEnabeled[12];			//

// Which Track Cuts will be tested.
	Bool_t					fTestPt;						//
	Bool_t 					fTestFilterMask;				//
	Bool_t 					fTestMaxEta;					//
	Bool_t					fTestMaxRapidity;				//
	Bool_t 					fTestFlags;						//
	Bool_t 					fTestTOFmismatch;				//
	Bool_t 					fTestPtDeptDCAcut;				//

// QA histograms for Data.
	TList*					fDataTrackQAHistos;				// 
	TH1F*					fHistDataPt[3];					//! Pt distribution of tracks passing this cut.
	TH3F*					fHistDataPhiEtaPt[3];			//! Pt, Eta, Phi distribution.
	TH1F*					fHistDataNTracks[3];			//! Number of tracks passing the cut per event (filling by EventIsDone()).
	TH1F*					fHistDataDCAxy[3];				//! DCA_{xy} distribution.
	TH1F*					fHistDataDCAz[3];				//! DCA_{z} distribution
	TH2F* 					fHistDataDCAxyOneSigma[12];		//! DCA_{xy} distribution of particles as identified by 1 sigma method.
	Int_t					fNTracks[12];					//! Number of tracks

	TH3F* 					fHistDataPID[3][3][5];			//! TPC/TOF v.s. pT, [charge][mass assumption][ptclass]
	TH2F*					fHistTOFMismatch[3][3][5];		//! TOF Mismatch histograms, [charge][mass assumption][ptclass]
	TH3F*					fHistTPCTOFMismatch[3][3][5];	//! TPC/TOF mismatch histograms (Same as TOF, but now the TPC hit of the track is included.)

// QA histograms for Primary Reconstructed MC tracks.
	TList*					fPrimRecMCTrackQAHistos;		//
	TH1F*					fHistPrimRecMCPt[12];			//! Pt distribution of reconstructed MC track passing this cut.
	TH1F*					fHistPrimRecNTracks[12];		//!
	TH2F*					fHistPrimRecMCDCA[12];			//! DCA_xy distribution of reconstructed MC track passing this cut.

// QA histograms for Primary Generated MC particles.
	TList*					fPrimGenMCTrackQAHistos;		//
	TH1F*					fHistPrimGenMCPt[12];			//! Pt distribution of generated MC particles passing this cut.

// QA histograms for Secondary Reconstructed MC tracks.
	TList*					fSecRecMCTrackQAHistos;			//
	TH1F*					fHistSecRecMCPt[12];			//! Pt distribution of reconstructed MC track passing this cut.
	TH2F*					fHistSecRecMCDCAMat[12];		//! DCA_xy distribution of material decay particles.
	TH2F*					fHistSecRecMCDCAWeak[12];		//! DCA_xy distribution of weak decay.

// QA histograms for Secondary Generated MC particles.
	TList*					fSecGenMCTrackQAHistos;			//
	TH1F*					fHistSecGenMCPt[12];			//! Pt distribution of generated MC particles passing this cut.

// Binning of all the histograms.
	Double_t				fPtAxis[57];					// Pt axis used in all histograms, except PID and Mismatch histograms.
	Int_t					fNPtBins;						// Number of bins in the pt-axis.
	Int_t					fNEtaBins;						// 
	Int_t					fNPhiBins;						//

	Double_t				fPtBoundaryPID[6];				// There are five different PID histo's. This array gives the pT range of these histograms.
	Int_t					fNPtBinsPID[5];					// This array gives the number of pT bins for each of these histograms.

	Double_t				fTOFLowerBound[5][3];			// These arrays give the lower and upper bound of the TOF axes,
	Double_t				fTOFUpperBound[5][3];			// for each species, as well as the number of bins. The numbers
	Int_t					fTOFbins[5][3];					// size of the array is [ptrange][species].

	Double_t				fTPCLowerBound[5][3];			// The same, but now for TPC.
	Double_t				fTPCUpperBound[5][3];
	Int_t					fTPCbins[5][3];

// Naming conventions of the histograms.
	TString					fHistoName[12];					// Names of the histogram classes.
	TString					fHistoLatex[12];				// Names of the histogram classes in LaTeX.
	TString					fParticleName[3];				// Names of the particles (Pion, Kaon, Proton)
	TString					fPtClassName[5];				// Names of the ptclasses (should only be for internal use)

	Int_t 					fDebug;							// Debug flag.

	ClassDef(AliAODTrackCutsDiHadronPID,5);

};

#endif
