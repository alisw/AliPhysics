#ifndef ALIAODTRACKCUTSDIHADRONPID_H
#define ALIAODTRACKCUTSDIHADRONPID_H

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

// -------------------------------------------------------------------------
//  Interface, methods used to get information about the track cuts, and to
//  retrieve filled histograms. 
// -------------------------------------------------------------------------

public:
	void PrintCuts();								// Gives an overview of the cuts.

// Return the list of QA histos
	TList* GetListOfDataQAHistos() const{
		if (fDataTrackQAHistos) {return fDataTrackQAHistos;}
		else return 0x0;
	}
	TList* GetListOfPrimRecMCTrackQAHistos() const{
		if (fPrimRecMCTrackQAHistos) {return fPrimRecMCTrackQAHistos;}
		else return 0x0;
	}
	TList* GetListOfPrimGenMCTrackQAHistos() const {
		if (fPrimGenMCTrackQAHistos) {return fPrimGenMCTrackQAHistos;}
		else return 0x0;
	}
	TList* GetListOfSecRecMCTrackQAHistos() const {
		if (fSecRecMCTrackQAHistos) return fSecRecMCTrackQAHistos;
		else return 0x0;
	}
	TList* GetListOfSecGenMCTrackQAHistos() const {
		if (fSecGenMCTrackQAHistos) return fSecGenMCTrackQAHistos;
		else return 0x0;
	}

// Note that the PID histograms have in principle a different number of pT bins.
// FIXME: This is not very nice...
	Int_t GetNPtBins() const {return fNPtBins;}
	Int_t GetNPtBinsPID(Int_t ptclass = -1) const {
		// if class = -1, then return the sum.
		if (ptclass == -1) {
			Int_t nptbinspid = 0;
			for (Int_t iPtClass = 0; iPtClass < 5; iPtClass++) {
				nptbinspid += fNPtBinsPID[iPtClass];
			}
			return nptbinspid;
		} else if (ptclass >= 0 && ptclass < 5) {
			return fNPtBinsPID[ptclass];
		} else {return -999;}
	}	

// Returns the Pt axis for PID and other histograms.
	Double_t* GetPtAxis() {return fPtAxis;}
	Float_t* GetPtAxisPID() const {
		const Int_t nptbinspid = GetNPtBinsPID();
		Float_t* ptaxis = new Float_t[nptbinspid];
		for (Int_t iPtBin = 0; iPtBin < nptbinspid; iPtBin++) {
			ptaxis[iPtBin] = GetPtMinPID(iPtBin + 1); 
		}
		ptaxis[nptbinspid] = GetPtMaxPID(nptbinspid);
		return ptaxis;
	}

// Return data histogram with a specific name. Since the histograms are not streamed, and only the
// TList containing them is, we have to retrieve them from the list.
	TObject* GetHistData(const char* name) const {return fDataTrackQAHistos->FindObject(name);} 
	TObject* GetHistPrimRecMC(const char* name) const {return fPrimRecMCTrackQAHistos->FindObject(name);}
	TObject* GetHistPrimGenMC(const char* name) const {return fPrimGenMCTrackQAHistos->FindObject(name);}
	TObject* GetHistSecRecMC(const char* name) const {return fSecRecMCTrackQAHistos->FindObject(name);}
	TObject* GetHistSecGenMC(const char* name) const {return fSecGenMCTrackQAHistos->FindObject(name);}

// Since we will often want to have TOF histograms, here are a few methods which return the 
// appropriate projections. The class does not own these projections, and the user must take care of them.
	TH1F* GetHistDataTOFProjection(Int_t charge, Int_t species, Int_t ptbin);
 	TH1F* GetHistDataTOFMismatch(Int_t charge, Int_t species, Int_t ptbin);
	Double_t GetPtMinPID(Int_t bin) const {
		Int_t ptclass = GetPtClass(bin);
		Int_t bininptclass = GetBinInPtClass(bin);
		Double_t minpt = fPtBoundaryPID[ptclass];
		Double_t maxpt = fPtBoundaryPID[ptclass+1];
		Double_t ptres = (maxpt - minpt)/((Double_t)fNPtBinsPID[ptclass]);
		return (minpt + ptres * ((Double_t)(bininptclass - 1)) );
	}
	Double_t GetPtMaxPID(Int_t bin) const {
		Int_t ptclass = GetPtClass(bin);
		Int_t bininptclass = GetBinInPtClass(bin);
		Double_t minpt = fPtBoundaryPID[ptclass];
		Double_t maxpt = fPtBoundaryPID[ptclass+1];
		Double_t ptres = (maxpt - minpt)/((Double_t)fNPtBinsPID[ptclass]);
		return (minpt + ptres * ((Double_t)(bininptclass)) );
	}
	Double_t GetPtClassMin(Int_t ptclass) const {
		if (ptclass >= 0 && ptclass < 5) {
			return fPtBoundaryPID[ptclass];
		} else {return -999;}
	}
	Double_t GetPtClassMax(Int_t ptclass) const {
		if (ptclass >= 0 && ptclass < 5) {
			return fPtBoundaryPID[ptclass+1];
		} else {return -999;}
	}

	Double_t GetPtBinWidthPID(Int_t bin) const {return (GetPtMaxPID(bin) - GetPtMinPID(bin)); }

// BE CAREFUL! Following methods do not apply to the ptaxis of the PID histograms! For that, call
// the GetPtMinPID and GetPtMaxPID methods. 
	Double_t GetPtMin(Int_t bin) const {
		if ((bin < 1) || (bin > fNPtBins + 1)) {cout<<"Bin is out of range..."<<endl; return -999.;}
		else {return fPtAxis[bin - 1];}
	}
	Double_t GetPtMax(Int_t bin) const {
		if ((bin < 1) || (bin > fNPtBins + 1)) {cout<<"Bin is out of range..."<<endl; return -999.;}
		else {return fPtAxis[bin];}
	}
	Double_t GetPtBinWidth(Int_t bin) const {return (GetPtMax(bin) - GetPtMin(bin)); }

	Int_t GetNTOFbins(Int_t ptclass, Int_t species) const {return fTOFbins[ptclass][species];}
	Double_t GetTOFmin(Int_t ptclass, Int_t species) const {return fTOFLowerBound[ptclass][species];}
	Double_t GetTOFmax(Int_t ptclass, Int_t species) const {return fTOFUpperBound[ptclass][species];}

	Int_t GetNTPCbins(Int_t ptclass, Int_t species) const {return fTPCbins[ptclass][species];}
	Double_t GetTPCmin(Int_t ptclass, Int_t species) const {return fTPCLowerBound[ptclass][species];}
	Double_t GetTPCmax(Int_t ptclass, Int_t species) const {return fTPCUpperBound[ptclass][species];}

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

// -------------------------------------------------------------------------
//  Methods used to configure the track cuts object, to be called at
//  initialization, i.e., before the object is added to an analysis task.  
// -------------------------------------------------------------------------

public:

// Request Certain QA histograms being filled.
	Bool_t RequestQAHistos(Int_t histoclass, Bool_t Enable3DSpectra = kFALSE, Bool_t EnablePIDHistos = kFALSE) {
		if ((histoclass > -1) && (histoclass < 12)) {
			fHistRequested[histoclass] = kTRUE;
			f3DSpectraEnabeled[histoclass] = Enable3DSpectra;
			fPIDHistosEnabeled[histoclass] = EnablePIDHistos;
			//cout<<"histoclass: "<<histoclass<<" requested: "<<fHistRequested[histoclass]<<endl;
			return kTRUE;
		} else { 
			return kFALSE;
		}
	}
	
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

// -------------------------------------------------------------------------
//  Methods called by the analysis task. 
// -------------------------------------------------------------------------

public:

// These two functions signal the beginning and the end of a new event.
	void StartNewEvent();				// Some things are set to zero.
	void EventIsDone(Bool_t IsMC);		// Some final histograms are filled.
	void CreateHistos(); 				// Should be called by the UserCreateOutput() function of the analysis task.

// Is Selected, for different types of tracks.
	Bool_t IsSelectedData(AliTrackDiHadronPID* track, Double_t randomhittime = -1.e20);
	Bool_t IsSelectedGeneratedMC(AliAODMCParticle* particle);
	Bool_t IsSelectedReconstructedMC(AliTrackDiHadronPID* track);

// -------------------------------------------------------------------------
//  Internal methods.
// -------------------------------------------------------------------------

public:

// For PID histograms we have a certain number of bins in pT, spread out over five
// large histograms, i.e., one for the lowest pT, and the biggest range in TOF/TPC,
// one for the higher pT and smaller range in TOF/TPC, etc. The following methods
// are a mapping between the total pT bin (what the user uses), and the pt bin
// within one of the five histograms (what's used internally)
	Int_t GetPtClass(const Int_t ptbin) const {
		
		// Returns a number [0..4]
		Int_t currentptclass = 0;
		Int_t currentptbin = fNPtBinsPID[0];
		while (currentptbin < ptbin) {
			currentptclass++;
			if (currentptclass == 5) {break;}
			currentptbin += fNPtBinsPID[currentptclass];
		}
		if (currentptclass == 5) {cout<<"GetPtClass -> ptbin out of range!"<<endl; return -1;}
		return currentptclass;
	}
	Int_t GetBinInPtClass(const Int_t ptbin) const {

		// Returns a number [1..Nbins]
		Int_t ptclass = GetPtClass(ptbin);
		if (ptclass == -1) {return -1;}

		Int_t ptbinout = ptbin;
		for (Int_t iPtClass = 0; iPtClass < ptclass; iPtClass++) {ptbinout -= fNPtBinsPID[iPtClass];}

		return ptbinout;

	}

private:

// Checks, return kTRUE if track passes the cut.
	Bool_t CheckPt(Double_t pt) const {
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
	TH1F* InitializeNTracksHisto(const char* name, Int_t histoclass);
	TH1F* InitializeDCAxyHisto(const char* name, Int_t histoclass);
	TH1F* InitializeDCAzHisto(const char* name, Int_t histoclass);
	TH3F* InitializeAcceptanceHisto(const char* /*name*/, Int_t /*histoclass*/); // TO BE IMPLEMENTED.
	TH2F* InitializeDCASpectrum(const char* name, Int_t histoclass);


	TH3F* InitializePIDHisto(const char* name, Int_t histoclass, Int_t expspecies, Int_t ptclass);
	TH2F* InitializeTOFMismatchHisto(const char* name, Int_t histoclass, Int_t expspecies, Int_t ptclass);

// -------------------------------------------------------------------------
//  Data members.
// -------------------------------------------------------------------------

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
	TH1F*					fHistDataNTracks[3];			//! Number of tracks passing the cut per event (filling by EventIsDone()).
	TH1F*					fHistDataDCAxy[3];				//! DCA_{xy} distribution.
	TH1F*					fHistDataDCAz[3];				//! DCA_{z} distribution
	TH2F* 					fHistDataDCAxyOneSigma[12];		//! DCA_{xy} distribution of particles as identified by 1 sigma method.
	Int_t					fNTracks[3];					//! Number of tracks

	TH3F* 					fHistDataPID[3][3][5];			//! TPC/TOF v.s. pT, [charge][mass assumption][ptclass]
	TH2F*					fHistTOFMismatch[3][3][5];		//! TOF Mismatch histograms, [charge][mass assumption][ptclass]

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
	Double_t				fPtAxis[47];					// Pt axis used in all histograms, except PID and Mismatch histograms.
	Int_t					fNPtBins;						// Number of bins in the pt-axis.

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

	ClassDef(AliAODTrackCutsDiHadronPID,4);

};

#endif
