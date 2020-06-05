#include "AliESDtrackCuts.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AlidNdPtTools.h"
#include "TGraph.h"
#include "THnSparse.h"
#include "TPDGCode.h"

class AlidNdPtTools;

using namespace std;

/// \cond CLASSIMP
ClassImp(AlidNdPtTools)
	/// \endcond

	//____________________________________________________________________________

	THnSparseF* AlidNdPtTools::fSparseTmp = 0;
	TGraph* AlidNdPtTools::fGsscale = 0;
	TGraph* AlidNdPtTools::fGsscale1 = 0;
	TGraph* AlidNdPtTools::fGsscale2 = 0;
	TGraph* AlidNdPtTools::fGsscaleB = 0;
	TGraph* AlidNdPtTools::fGsscaleB1 = 0;
	TGraph* AlidNdPtTools::fGsscaleB2 = 0;

	//____________________________________________________________________________

	/// Function to fill THnSparse or THn with up to 12 dimensions
	///
	/// Ugly, but works and makes filling easier
	///
	/// \param s   Pointer to histogram to be filled
	/// \param x0  x of dimension 0
	/// \param x1  x of dimension 1
	/// \param x2  x of dimention 2
	/// \param x3  x of dimention 3
	/// \param x4  x of dimention 4
	/// \param x5  x of dimention 5
	/// \param x6  x of dimention 6
	/// \param x7  x of dimention 7
	/// \param x8  x of dimention 8
	/// \param x9  x of dimention 9
	/// \param x10 x of dimention 10
	/// \param x11 x of dimention 11
	///
	/// \return return value of THnBase->Fill(...) or 0 in case of error

	Long64_t AlidNdPtTools::FillHist(THnBase* s, Double_t x0, Double_t x1,
			Double_t x2, Double_t x3, Double_t x4,
			Double_t x5, Double_t x6, Double_t x7,
			Double_t x8, Double_t x9, Double_t x10,
			Double_t x11) {
		if (s->GetNdimensions() > 12) {
			return 0;
		}
		Double_t vals[12] = {x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11};
		return s->Fill(vals);
	}

//____________________________________________________________________________

/// Function to fill THnSparse or THn with up to 12 dimensions with weight
///
/// \param w   weight to be used for Histogram filling
/// \param s   Pointer to histogram to be filled
/// \param x0  x of dimension 0
/// \param x1  x of dimension 1
/// \param x2  x of dimention 2
/// \param x3  x of dimention 3
/// \param x4  x of dimention 4
/// \param x5  x of dimention 5
/// \param x6  x of dimention 6
/// \param x7  x of dimention 7
/// \param x8  x of dimention 8
/// \param x9  x of dimention 9
/// \param x10 x of dimention 10
/// \param x11 x of dimention 11
///
/// \return return value of THnBase->Fill(...) or 0 in case of error

Long64_t AlidNdPtTools::FillHist(Double_t w, THnBase* s, Double_t x0,
		Double_t x1, Double_t x2, Double_t x3,
		Double_t x4, Double_t x5, Double_t x6,
		Double_t x7, Double_t x8, Double_t x9,
		Double_t x10, Double_t x11) {
	if (s->GetNdimensions() > 12) {
		return 0;
	}
	Double_t vals[12] = {x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11};
	return s->Fill(vals, w);
}

//____________________________________________________________________________

/// Function to fill THnSparse or THn with weight
///
/// \param weight   weight to be used for Histogram filling
/// \param s   Pointer to histogram to be filled
/// \param val vector of values
///
/// \return return value of THnBase->Fill(...) or 0 in case of error

Long64_t AlidNdPtTools::FillHistWeighted(THnBase* s, std::vector<double> const& val, double weight) {
    if(val.size() < 1) return 0;
    return s->Fill(val.data(), weight);
}

//____________________________________________________________________________

/// Add an Axis (Dimension) to the THnSparseF
///
/// function to add a user defined axes
/// with normal (linear) binning
/// number of bins, upper and lower range
/// has to be provided
///
/// \param label axis label (short version)
/// \param title axis title (long version with latex code)
/// \param nbins number of bins
/// \param xmin lower edge of first bin
/// \param xmax upper edge of last bin
/// \param option not used and currently ignored
///
/// \return the total number of dimensions after adding this axis, 0 in case of
/// error

Int_t AlidNdPtTools::AddAxis(const char* label, const char* title, Int_t nbins,
		Double_t xmin, Double_t xmax, const char* option) {
	Int_t n = 1;
	if (fSparseTmp) {
		n += fSparseTmp->GetNdimensions();
	}
	TString s;
	TArrayI bin(n);
	TArrayD min(n);
	TArrayD max(n);
	for (int i = 0; i < n - 1; i++) {
		bin[i] = fSparseTmp->GetAxis(i)->GetNbins();
		min[i] = fSparseTmp->GetAxis(i)->GetXmin();
		max[i] = fSparseTmp->GetAxis(i)->GetXmax();
		s += fSparseTmp->GetAxis(i)->GetName();
		s += ":";
	}
	bin[n - 1] = nbins;
	min[n - 1] = xmin;
	max[n - 1] = xmax;
	s += label;
	THnSparseF* h = new THnSparseF("fSparseTmp", s.Data(), n, bin.GetArray(),
			min.GetArray(), max.GetArray());
	for (int i = 0; i < n - 1; i++) {
		if (fSparseTmp->GetAxis(i)->GetXbins() &&
				fSparseTmp->GetAxis(i)->GetXbins()->GetSize()) {
			h->SetBinEdges(i, fSparseTmp->GetAxis(i)->GetXbins()->GetArray());
		}
		h->GetAxis(i)->SetTitle(fSparseTmp->GetAxis(i)->GetTitle());
		h->GetAxis(i)->SetName(fSparseTmp->GetAxis(i)->GetName());
	}
	h->GetAxis(n - 1)->SetTitle(title);
	h->GetAxis(n - 1)->SetName(label);
	if (fSparseTmp) {
		delete fSparseTmp;
	}
	fSparseTmp = h;
	return fSparseTmp->GetNdimensions();
}

//____________________________________________________________________________

/// Add an Axis (Dimension) to the THnSparseF
///
/// function to add a user defined axes
/// with normal (linear) binning
/// number of bins, upper and lower range
/// has to be provided
///
/// \param label this is used as label AND title for the new axis
/// \param nbins number of bins
/// \param xmin lower edge of first bin
/// \param xmax upper edge of last bin
/// \param option not used and currently ignored
///
/// \return the total number of dimensions after adding this axis, 0 in case of
/// error

Int_t AlidNdPtTools::AddAxis(const char* label, Int_t nbins, Double_t xmin,
		Double_t xmax, const char* option) {
	return AddAxis(label, label, nbins, xmin, xmax, option);
}

//____________________________________________________________________________

/// Add an Axis (Dimension) to the THnSparseF
///
/// function to add a user defined binning with
/// the option of variable bin size
/// number of bins and and an array defining the bin edges
/// has to be provided
///
/// \param label axis label (short version)
/// \param title axis title (long version with latex code)
/// \param nbins number of bins
/// \param xbins array of length nbins+1 containing the bin edges
/// \param option not used and currently ignored
///
/// \return the total number of dimensions after adding this axis, 0 in case of
/// error

Int_t AlidNdPtTools::AddAxis(const char* label, const char* title, Int_t nbins,
		Double_t* xbins, const char* option) {
	Int_t n = AddAxis(label, title, nbins, xbins[0], xbins[nbins], option);
	fSparseTmp->SetBinEdges(n - 1, xbins);
	return n;
}

//____________________________________________________________________________

/// Add an Axis (Dimension) to the THnSparseF
///
/// function to add a user defined binning with
/// the option of variable bin size
/// number of bins and and an array defining the bin edges
/// has to be provided
///
/// \param label this is used as label AND title for the new axis
/// \param nbins number of bins
/// \param xbins array of length nbins+1 containing the bin edges
/// \param option not used and currently ignored
///
/// \return the total number of dimensions after adding this axis, 0 in case of
/// error

Int_t AlidNdPtTools::AddAxis(const char* label, Int_t nbins, Double_t* xbins,
		const char* option) {
	return AddAxis(label, label, nbins, xbins, option);
}

//____________________________________________________________________________

/// Add an Axis (Dimension) to the THnSparseF
///
/// function to add from a series of pre-defined options
/// option supplied in not case-senstitiv
///
/// currently the following ones are available
/// "pt"                standard pT axis
/// "ptfew"             reduced pt binning
/// "ptveryfew"         much reduced pt binning
/// "ptmario"           marios pt binning
/// "cent"              standard centrality binning
/// "varsig35"
/// "mult6kfine"        multiplicity bining 0-6000 in fine bins
/// "mult6kcoarse"      multiplicity bining 0-6000 in coarse bins
/// "mult100kcoarse"    multiplicity bining 0-100000 in coarse bins
///
/// \param label axis label (short version)
/// \param title axis title (long version with latex code)
/// \param option string to steer the binning
///
/// \return the total number of dimensions after adding this axis, 0 in case of
/// error

Int_t AlidNdPtTools::AddAxis(const char* label, const char* title,
		const char* option) {
	TString o(option);
	o.ToLower();
	if (o.Contains("ptfew")) {
		const Int_t nbins = 21;
		Double_t xbins[22] = {0.0, 0.1, 0.2,  0.3,  0.4,  0.5,  0.6, 0.7,
			0.8, 0.9, 1.0,  1.1,  1.2,  1.3,  1.4, 1.5,
			2.0, 5.0, 10.0, 20.0, 50.0, 200.0};
		return AddAxis(label, title, nbins, xbins);
	}
	if (o.Contains("ptveryfew")) {
		const Int_t nbins = 8;
		Double_t xbins[9] = {0.0, 0.15, 0.5, 1.0, 2.0, 5.0, 10, 25.0, 200.0};
		return AddAxis(label, title, nbins, xbins);
	}
	if (o.Contains("ptmario")) {
		const Int_t nbins = 52;
		Double_t xbins[53] = {
			0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45, 0.5, 0.55, 0.6,
			0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95, 1.0,  1.1, 1.2,  1.3,
			1.4,  1.5,  1.6,  1.7,  1.8,  1.9,  2.0,  2.2,  2.4, 2.6,  2.8,
			3.0,  3.2,  3.4,  3.6,  3.8,  4.0,  4.5,  5.0,  5.5, 6.0,  6.5,
			7.0,  8.0,  9.0,  10.0, 20.0, 30.0, 40.0, 50.0, 60.0};
		return AddAxis(label, title, nbins, xbins);
	}
	if (o.Contains("pt")) {
		const Int_t nbins = 81;
		Double_t xbins[82] = {
			0.0,   0.05, 0.1,  0.15,  0.2,   0.25,  0.3,   0.35,  0.4,   0.45,
			0.5,   0.55, 0.6,  0.65,  0.7,   0.75,  0.8,   0.85,  0.9,   0.95,
			1.0,   1.1,  1.2,  1.3,   1.4,   1.5,   1.6,   1.7,   1.8,   1.9,
			2.0,   2.2,  2.4,  2.6,   2.8,   3.0,   3.2,   3.4,   3.6,   3.8,
			4.0,   4.5,  5.0,  5.5,   6.0,   6.5,   7.0,   8.0,   9.0,   10.0,
			11.0,  12.0, 13.0, 14.0,  15.0,  16.0,  18.0,  20.0,  22.0,  24.0,
			26.0,  28.0, 30.0, 32.0,  34.0,  36.0,  40.0,  45.0,  50.0,  60.0,
			70.0,  80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0,
			180.0, 200.0};
		return AddAxis(label, title, nbins, xbins);
	}
	if (o.Contains("cent")) {
		const Int_t nbins = 11;
		Double_t xbins[12] = {0.,  5.,  10., 20., 30., 40.,
			50., 60., 70., 80., 90., 100.};
		return AddAxis(label, title, nbins, xbins);
	}
	if (o.Contains("varsig35")) {
		const Int_t nbins = 35;
		Double_t xbins[36] = {-1, 0,  1,   2,   3,   4,   5,   6,    7,
			8,  9,  10,  11,  12,  13,  14,  15,   16,
			17, 18, 19,  20,  30,  40,  50,  60,   70,
			80, 90, 100, 200, 300, 400, 500, 1000, 2000};
		return AddAxis(label, title, nbins, xbins);
	}
	if (o.Contains("mult6kfine")) {
		// variable mult binning total 0-6000
		// 1-width 0-100            101
		// 10-width 100-1000         90
		// 100-width 1000-6000       50
		const Int_t nbins = 241;
		Double_t xbins[242];
		xbins[0] = -0.5;
		int i = 0;
		for (; i <= 100; i++) {
			xbins[i + 1] = xbins[i] + 1;
		}
		for (; i <= 100 + 90; i++) {
			xbins[i + 1] = xbins[i] + 10;
		}
		for (; i <= 100 + 90 + 50; i++) {
			xbins[i + 1] = xbins[i] + 100;
		}
		return AddAxis(label, title, nbins, xbins);
	}
	if (o.Contains("mult6kcoarse")) {
		// variable mult binning total 0-6000
		// 1-width 0-10              11
		// 10-width 10-100            9
		// 100-width 100-1000         9
		// 200-width 1000-6000       25
		const Int_t nbins = 54;
		Double_t xbins[55];
		xbins[0] = -0.5;
		int i = 0;
		for (; i <= 10; i++) {
			xbins[i + 1] = xbins[i] + 1;
		}
		for (; i <= 10 + 9; i++) {
			xbins[i + 1] = xbins[i] + 10;
		}
		for (; i <= 10 + 9 + 9; i++) {
			xbins[i + 1] = xbins[i] + 100;
		}
		for (; i <= 10 + 9 + 9 + 25; i++) {
			xbins[i + 1] = xbins[i] + 200;
		}
		return AddAxis(label, title, nbins, xbins);
	}
	if (o.Contains("mult100kcoarse")) {
		// variable mult binning total 0-100000
		// 1-width 0-10             11
		// 10-width 10-100           9
		// 100-width 100-1000        9
		// 1000-width 1000-10000     9
		// 10000-width 10000-100000  9
		const Int_t nbins = 47;
		Double_t xbins[48];
		xbins[0] = -0.5;
		int i = 0;
		for (; i <= 10; i++) {
			xbins[i + 1] = xbins[i] + 1;
		}
		for (; i <= 10 + 9; i++) {
			xbins[i + 1] = xbins[i] + 10;
		}
		for (; i <= 10 + 9 + 9; i++) {
			xbins[i + 1] = xbins[i] + 100;
		}
		for (; i <= 10 + 9 + 9 + 9; i++) {
			xbins[i + 1] = xbins[i] + 1000;
		}
		for (; i <= 10 + 9 + 9 + 9 + 9; i++) {
			xbins[i + 1] = xbins[i] + 10000;
		}
		return AddAxis(label, title, nbins, xbins);
	}

	return 0;
}

//____________________________________________________________________________

/// Add an Axis (Dimension) to the THnSparseF
///
/// function to add from a series of pre-defined options
/// option supplied in not case-senstitiv
///
/// currently the following ones are available
/// "pt"                standard pT axis
/// "ptfew"             reduced pt binning
/// "ptveryfew"         much reduced pt binning
/// "ptmario"           marios pt binning
/// "cent"              standard centrality binning
/// "varsig35"
/// "mult6kfine"        multiplicity bining 0-6000 in fine bins
/// "mult6kcoarse"      multiplicity bining 0-6000 in coarse bins
/// "mult100kcoarse"    multiplicity bining 0-100000 in coarse bins
///
/// \param label this is used as label AND title for the new axis
/// \param option string to steer the binning
///
/// \return the total number of dimensions after adding this axis, 0 in case of
/// error

Int_t AlidNdPtTools::AddAxis(const char* label, const char* option) {
	return AddAxis(label, label, option);
}

//____________________________________________________________________________

/// Add an Axis (Dimension) to the THnSparseF
///
/// function to add from a series of pre-defined options
/// option supplied in not case-senstitiv
///
/// currently the following ones are available
/// "pt"                standard pT axis
/// "ptfew"             reduced pt binning
/// "ptveryfew"         much reduced pt binning
/// "ptmario"           marios pt binning
/// "mcpt"              standard pT axis for MC
/// "mcptfew"           reduced pt binning for MC
/// "mcptveryfew"       much reduced pt binning
/// "mcptmario"         marios pt binning
/// "cent"              standard centrality binning
/// "mult6kfine"        multiplicity bining 0-6000 in fine bins
/// "mult6kcoarse"      multiplicity bining 0-6000 in coarse bins
/// "mult100kcoarse"    multiplicity bining 0-100000 in coarse bins
///
/// \param option string to steer the binning
///
/// \return the total number of dimensions after adding this axis, 0 in case of
/// error

Int_t AlidNdPtTools::AddAxis(const char* option) {
	TString o(option);
	o.ToLower();
	if (o.EqualTo("pt"))
		return AddAxis("pT", "#it{p}_{T} (GeV/#it{c})", "pt");
	if (o.EqualTo("mcpt"))
		return AddAxis("MCpT", "it{p}_{T,MC} (GeV/#it{c})", "pt");
	if (o.EqualTo("ptfew"))
		return AddAxis("pT", "it{p}_{T} (GeV/#it{c})", "ptfew");
	if (o.EqualTo("mcptfew"))
		return AddAxis("MCpT", "it{p}_{T,MC} (GeV/#it{c})", "ptfew");
	if (o.EqualTo("ptveryfew"))
		return AddAxis("pT", "it{p}_{T} (GeV/#it{c})", "ptveryfew");
	if (o.EqualTo("mcptveryfew"))
		return AddAxis("MCpT", "it{p}_{T,MC} (GeV/#it{c})", "ptveryfew");
	if (o.EqualTo("ptmario"))
		return AddAxis("pT", "it{p}_{T} (GeV/#it{c})", "ptmario");
	if (o.EqualTo("mcptmario"))
		return AddAxis("MCpT", "it{p}_{T,MC} (GeV/#it{c})", "ptmario");
	if (o.EqualTo("mult6kfine"))
		return AddAxis("mult", "N", "mult6kfine");
	if (o.EqualTo("mult6kcoarse"))
		return AddAxis("mult", "N", "mult6kcoarse");
	if (o.EqualTo("mult1000kcoarse"))
		return AddAxis("mult", "N", "mult1000kcoarse");
	if (o.EqualTo("cent"))
		return AddAxis("cent", "centrality (%)", "cent");
	return 0;
}

//____________________________________________________________________________

/// Create a THnSparseF histogram
///
/// Before this function actually creates a histogram
/// axis have to be added using the various AddAxis() functions
///
/// \param name  name of the histogram
///
/// \return newly created THnSparseF histogram or 0 in case fo error

THnSparseF* AlidNdPtTools::CreateHist(const char* name) {
	if (!fSparseTmp)
		return 0;
	THnSparseF* h = fSparseTmp;
	h->SetName(name);
	fSparseTmp = 0;
	return h;
}

//____________________________________________________________________________

/// Create a TH1D histogram for Logging purposes
///
/// \param name  name of the histogram
/// \param title title of the histogram
///
/// \return newly created TH1D histogram

TH1D* AlidNdPtTools::CreateLogHist(const char* name, const char* title) {
	TH1D* h = 0;
	if (title) {
		h = new TH1D(name, title, 200, 0, 200);
	} else {
		h = new TH1D(name, name, 200, 0, 200);
	}
	return h;
}

//____________________________________________________________________________

/// Create a TH1D histogram for Logging purposes
///
/// \param name  name and title of the histogram
///
/// \return newly created TH1D histogram

TH1D* AlidNdPtTools::CreateLogHist(const char* name) {
	return CreateLogHist(name, name);
}

//____________________________________________________________________________

/// Function to create AliESDtrackCuts with various settings
///
/// specfiy the type of track cuts to be created in the option string
/// this is not case senstitiv
///
/// currently the folowing options are implemented:
/// ""                          identical to "TPCITSgeo"
/// "default"                   identical to "TPCITSgeo"
/// "TPCITSgeo"                 all standard track cuts for TPC-ITS primaries
/// including the golden chi2 cut and the geometric cut "TPCITSgeoNoDCAr" same
/// as "TPCITSgeo" but without the DCAr cut "TPCITSforDCArStudy"
/// "TPCITSgeoNoDCAr" but without the golden chi2 cut (for DCAr fits) "TPCgeo"
/// all default TPC cuts, including the geometric length cut "TPCgeoNoDCAr" same
/// as "TPCgeo" but without the cut on DCAr (for DCAr fits) "TPCgeo+ITShit" same
/// as "TPCgeo" but in addition require a hit in the ITS (for Matching Efficieny
/// studies) "TPCgeo+ITSrefit"           same as "TPCgeo" but in addition
/// require the ITS refit, i.e. 2 hitsin the ITS (for Matching Efficieny
/// studies) "TPCgeo+SPDhit"             same as "TPCgeo" but in addition
/// require a hit in any layer of the SPD (for Matching Efficieny studies)
/// "TPCgeo+SPDhit+ITSrefit"    identical to "TPCgeo+ITSrefit+SPDhit"
/// "TPCgeo+ITSrefit+SPDhit"    same as "TPCgeo" but in addition require a hit
/// in any layer of the SPD AND a ITS refit (for Matching Efficieny studies)
/// "TPConlyMinimal"            minimal tpc only cuts
///
/// all these cuts do not include the eta-cut, to add add one  of the following
/// options: "Eta05" "Eta08" "Eta10"
///
/// For variation of geometrical length cut put one of the follow options in
/// addition:
///  "GeoCutVar1"
///  "GeoCutVar2"
///  "GeoCutVar3"
///  "GeoCutVar4"
///  "GeoCutVar5"
///  "GeoCutVar6"
///  see
///  https://alice-notes.web.cern.ch/system/files/notes/analysis/472/2018-Mar-29-analysis_note-Charged_particle_spectra_pp5TeV_v4.pdf
///  for futher information
///
/// \param option string to select the type of track cuts to be created
///
/// \return the newly created AliESDtrackCuts or 0 in case of error

AliESDtrackCuts* AlidNdPtTools::CreateESDtrackCuts(const char* option, int _cutMode, bool _SaveHistos) {
	TString o(option);
	AliESDtrackCuts* cuts = new AliESDtrackCuts(o.Data());
	o.ToLower();
    
	// if eta ranges is provided set the eta range
	// and remove the part of the string containting the eta range
	if (o.Contains("eta05")) {
		cuts->SetEtaRange(-0.5, 0.5);
		o.ReplaceAll("eta05", "");
	} else if (o.Contains("eta08")) {
		cuts->SetEtaRange(-0.8, 0.8);
		o.ReplaceAll("eta08", "");
	} else if (o.Contains("eta10")) {
		cuts->SetEtaRange(-1.0, 1.0);
		o.ReplaceAll("eta10", "");
	}

	// as default use the cuts with geometric length cut
	if ((o.EqualTo("")) || (o.EqualTo("default")) || 100!=_cutMode) {
		o = "tpcitsgeo";
	}

	if (o.EqualTo("tpcitsgeo")) { // default
		//         cuts = new AliESDtrackCuts("default TPCITS with geo L cut");
		cuts->SetRequireTPCRefit(kTRUE);
		cuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
		cuts->SetMaxChi2PerClusterTPC(4);
		cuts->SetMaxFractionSharedTPCClusters(0.4);
		cuts->SetRequireITSRefit(kTRUE);
		cuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
				AliESDtrackCuts::kAny);
		cuts->SetMaxChi2PerClusterITS(36.);
		cuts->SetDCAToVertex2D(kFALSE);
		cuts->SetRequireSigmaToVertex(kFALSE);
		cuts->SetMaxDCAToVertexZ(2.0);
		// 7*(0.0026+0.0050/pt^1.01)
		cuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
		cuts->SetAcceptKinkDaughters(kFALSE);
		// tpcc cut
		cuts->SetMaxChi2TPCConstrainedGlobal(36.);
		// Geometrical-Length Cut
		cuts->SetCutGeoNcrNcl(3, 130, 1.5, 0.85, 0.7);
		//         cuts->SetEtaRange(-0.8,0.8);

		if(_cutMode == 101) {cuts->SetMaxChi2PerClusterITS(25.);}
		if(_cutMode == 102) {cuts->SetMaxChi2PerClusterITS(49.);}

		if(_cutMode == 103) {cuts->SetMaxChi2PerClusterTPC(3.0); }
		if(_cutMode == 104) {cuts->SetMaxChi2PerClusterTPC(5.0); }

		if(_cutMode == 105) {cuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);}
		if(_cutMode == 106) {cuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);}

		if(_cutMode == 107) {cuts->SetMaxFractionSharedTPCClusters(0.2);}
		if(_cutMode == 108) {cuts->SetMaxFractionSharedTPCClusters(1.0);}

		if(_cutMode == 109) {cuts->SetMaxChi2TPCConstrainedGlobal(25.);}
		if(_cutMode == 110) {cuts->SetMaxChi2TPCConstrainedGlobal(49.);}

		if(_cutMode == 111) {cuts->SetMaxDCAToVertexXYPtDep("0.0104+0.0200/pt^1.01");}
		if(_cutMode == 112) {cuts->SetMaxDCAToVertexXYPtDep("0.0260+0.0500/pt^1.01");}

		if(_cutMode == 113) {cuts->SetMaxDCAToVertexZ(1.0);}
		if(_cutMode == 114) {cuts->SetMaxDCAToVertexZ(5.0);}

		if(_cutMode == 115) {cuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);}

		if(_cutMode == 116) {cuts->SetCutGeoNcrNcl(3,120,1.5,0.85,0.7);}
		if(_cutMode == 117) {cuts->SetCutGeoNcrNcl(3,140,1.5,0.85,0.7);}

		if(_cutMode == 118) {cuts->SetCutGeoNcrNcl(4,130,1.5,0.85,0.7);}
		if(_cutMode == 119) {cuts->SetCutGeoNcrNcl(2,130,1.5,0.85,0.7);}

	} else if (o.EqualTo("tpcitsnogeo")) {
		cuts->SetRequireTPCRefit(kTRUE);
		cuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
		cuts->SetMaxChi2PerClusterTPC(4);
		cuts->SetMaxFractionSharedTPCClusters(0.4);
		cuts->SetRequireITSRefit(kTRUE);
		cuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
				AliESDtrackCuts::kAny);
		cuts->SetMaxChi2PerClusterITS(36.);
		cuts->SetDCAToVertex2D(kFALSE);
		cuts->SetRequireSigmaToVertex(kFALSE);
		cuts->SetMaxDCAToVertexZ(2.0);
		// 7*(0.0026+0.0050/pt^1.01)
		cuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
		cuts->SetAcceptKinkDaughters(kFALSE);
		// tpcc cut
		cuts->SetMaxChi2TPCConstrainedGlobal(36.);
		// Geometrical-Length Cut
		//        cuts->SetCutGeoNcrNcl(3,130,1.5,0.85,0.7);
	} else if (o.EqualTo("tpcitsnogeonogold")) {
        cuts->SetRequireTPCRefit(kTRUE);
        cuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
        cuts->SetMaxChi2PerClusterTPC(4);
        cuts->SetMaxFractionSharedTPCClusters(0.4);
        cuts->SetRequireITSRefit(kTRUE);
        cuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                AliESDtrackCuts::kAny);
        cuts->SetMaxChi2PerClusterITS(36.);
        cuts->SetDCAToVertex2D(kFALSE);
        cuts->SetRequireSigmaToVertex(kFALSE);
        cuts->SetMaxDCAToVertexZ(2.0);
        // 7*(0.0026+0.0050/pt^1.01)
        cuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
        cuts->SetAcceptKinkDaughters(kFALSE);
        // tpcc cut
//        cuts->SetMaxChi2TPCConstrainedGlobal(36.);
        // Geometrical-Length Cut
        //        cuts->SetCutGeoNcrNcl(3,130,1.5,0.85,0.7);
    } else if (o.EqualTo("tpcitsgeonodcar")) {
		//         cuts = new AliESDtrackCuts("default TPCITS with geo L cut
		//         without DCAr");
		cuts->SetRequireTPCRefit(kTRUE);
		cuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
		cuts->SetMaxChi2PerClusterTPC(4);
		cuts->SetMaxFractionSharedTPCClusters(0.4);
		cuts->SetRequireITSRefit(kTRUE);
		cuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
				AliESDtrackCuts::kAny);
		cuts->SetMaxChi2PerClusterITS(36.);
		cuts->SetDCAToVertex2D(kFALSE);
		cuts->SetRequireSigmaToVertex(kFALSE);
		cuts->SetMaxDCAToVertexZ(2.0);
		cuts->SetAcceptKinkDaughters(kFALSE);
		// tpcc cut
		cuts->SetMaxChi2TPCConstrainedGlobal(36.);
		// Geometrical-Length Cut
		cuts->SetCutGeoNcrNcl(3, 130, 1.5, 0.85, 0.7);
		//         cuts->SetEtaRange(-0.8,0.8);

	} else if (o.EqualTo("tpcitsfordcarstudy")) {
		//         cuts = new AliESDtrackCuts("default TPCITS with geo L cut
		//         without DCAr and Chi2 TPCc vs. Global");
		cuts->SetRequireTPCRefit(kTRUE);
		cuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
		cuts->SetMaxChi2PerClusterTPC(4);
		cuts->SetMaxFractionSharedTPCClusters(0.4);
		cuts->SetRequireITSRefit(kTRUE);
		cuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
				AliESDtrackCuts::kAny);
		cuts->SetMaxChi2PerClusterITS(36.);
		cuts->SetDCAToVertex2D(kFALSE);
		cuts->SetRequireSigmaToVertex(kFALSE);
		cuts->SetMaxDCAToVertexZ(2.0);
		cuts->SetAcceptKinkDaughters(kFALSE);
		// tpcc cut
		// cuts->SetMaxChi2TPCConstrainedGlobal(36.);
		// Geometrical-Length Cut
		cuts->SetCutGeoNcrNcl(3, 130, 1.5, 0.85, 0.7);
		//         cuts->SetEtaRange(-0.8,0.8);

	} else if (o.EqualTo("tpcgeo")) {
		//         cuts = new AliESDtrackCuts("TPConly with geo L");
		cuts->SetRequireTPCRefit(kTRUE);
		cuts->SetAcceptKinkDaughters(kFALSE);
		cuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
		cuts->SetMaxChi2PerClusterTPC(4);
		cuts->SetMaxFractionSharedTPCClusters(0.4);
		cuts->SetMaxDCAToVertexZ(3);
		cuts->SetMaxDCAToVertexXY(3);
		// Geometrical-Length Cut
		cuts->SetCutGeoNcrNcl(3, 130, 1.5, 0.85, 0.7);
		//         cuts->SetEtaRange(-0.8,0.8);

	} else if (o.EqualTo("tpcgeonodcar")) {
		cuts = new AliESDtrackCuts("TPConly with geo L without DCAr");
		cuts->SetRequireTPCRefit(kTRUE);
		cuts->SetAcceptKinkDaughters(kFALSE);
		cuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
		cuts->SetMaxChi2PerClusterTPC(4);
		cuts->SetMaxFractionSharedTPCClusters(0.4);
		cuts->SetMaxDCAToVertexZ(3);
		// Geometrical-Length Cut
		cuts->SetCutGeoNcrNcl(3, 130, 1.5, 0.85, 0.7);
		//         cuts->SetEtaRange(-0.8,0.8);

	} else if (o.EqualTo("tpcgeo+itshit")) {
		//         cuts = new AliESDtrackCuts("TPConly with geo L + hit in
		//         ITS");
		cuts->SetRequireTPCRefit(kTRUE);
		cuts->SetAcceptKinkDaughters(kFALSE);
		cuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
		cuts->SetMaxChi2PerClusterTPC(4);
		cuts->SetMaxFractionSharedTPCClusters(0.4);
		cuts->SetMaxDCAToVertexZ(3);
		cuts->SetMaxDCAToVertexXY(3);
		// Geometrical-Length Cut
		cuts->SetCutGeoNcrNcl(3, 130, 1.5, 0.85, 0.7);
		// its hit
		cuts->SetMinNClustersITS(1);
		//         cuts->SetEtaRange(-0.8,0.8);

	} else if (o.EqualTo("tpcgeo+itsrefit")) {
		//         cuts = new AliESDtrackCuts("TPC with geo L + ITSrefit");
		cuts->SetRequireTPCRefit(kTRUE);
		cuts->SetAcceptKinkDaughters(kFALSE);
		cuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
		cuts->SetMaxChi2PerClusterTPC(4);
		cuts->SetMaxFractionSharedTPCClusters(0.4);
		cuts->SetMaxDCAToVertexZ(3);
		cuts->SetMaxDCAToVertexXY(3);
		// Geometrical-Length Cut
		cuts->SetCutGeoNcrNcl(3, 130, 1.5, 0.85, 0.7);
		// its refit
		cuts->SetRequireITSRefit(kTRUE);

	} else if (o.EqualTo("tpcgeo+spdhit")) {
		//         cuts = new AliESDtrackCuts("TPC with geo L + hit in SPD");
		cuts->SetRequireTPCRefit(kTRUE);
		cuts->SetAcceptKinkDaughters(kFALSE);
		cuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
		cuts->SetMaxChi2PerClusterTPC(4);
		cuts->SetMaxFractionSharedTPCClusters(0.4);
		cuts->SetMaxDCAToVertexZ(3);
		cuts->SetMaxDCAToVertexXY(3);
		// Geometrical-Length Cut
		cuts->SetCutGeoNcrNcl(3, 130, 1.5, 0.85, 0.7);
		// spd hit
		cuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
				AliESDtrackCuts::kAny);
		//         cuts->SetEtaRange(-0.8,0.8);

	} else if (o.EqualTo("tpcgeo+itsrefit+spdhit") ||
			o.EqualTo("tpcgeo+spdhit+itsrefit")) {
		//         cuts = new AliESDtrackCuts("TPC with geo L + ITSrefit + hit
		//         in SPD");
		cuts->SetRequireTPCRefit(kTRUE);
		cuts->SetAcceptKinkDaughters(kFALSE);
		cuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
		cuts->SetMaxChi2PerClusterTPC(4);
		cuts->SetMaxFractionSharedTPCClusters(0.4);
		cuts->SetMaxDCAToVertexZ(3);
		cuts->SetMaxDCAToVertexXY(3);
		// Geometrical-Length Cut
		cuts->SetCutGeoNcrNcl(3, 130, 1.5, 0.85, 0.7);
		// its refit + spd hit
		cuts->SetRequireITSRefit(kTRUE);
		cuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
				AliESDtrackCuts::kAny);
		//         cuts->SetEtaRange(-0.8,0.8);

	} else if (o.EqualTo("tpconlyminimal")) {
		//         cuts = new AliESDtrackCuts("minimal TPC only cuts");
		cuts->SetRequireTPCRefit(kTRUE);
		cuts->SetMaxChi2PerClusterTPC(4);
		cuts->SetMaxFractionSharedTPCClusters(0.4);
		cuts->SetMinNClustersTPC(50);
		//         cuts->SetEtaRange(-1.0,1.0);
	} else {
		// in case of no valid argument for cuts supplied
		// return 0 to avoid mistakes
		delete cuts;
		cuts = 0;
	}

    if(_SaveHistos){
        cuts->DefineHistograms(kRed);
    }
    
	return cuts;
}

//____________________________________________________________________________

/// Retrieve the scaling factor for MC primaries and secondaries
///
/// WARNING! this is only for LHC17pq for now!
///
/// This method looks up the proper scaling factors and returns them
/// for online use
///
/// there is the systflag parameter
/// to determine the type of correction.
/// currenlty implemented:
///  0: nominal (average correction)
/// -1: systematic varation down (minimal correction)
/// +1: systematic varation up (maximal correction)
///
/// period index
/// currenlty two periods supported:
/// periodindex 0: pp 13 TeV, LHC18b  -- this is the (new) default
/// periodindex 1: pp 5TeV, LHC17pq   -- previoulsy that was the default! be
/// aware!
///
/// TODO: add other fits, implement automatic switching
///
/// \param particle MC particle
/// \param event    ESD event
/// \param systflag Flag for syst variation
/// \param periodindex Temporary solution to select the period for the secondary
/// scaling
///
/// \return scaling factor accoring to the supplied arguments

Double_t AlidNdPtTools::MCScalingFactor(AliMCParticle* particle,
		AliMCEvent* event, Int_t systflag,
		Int_t periodindex) {
	// event multiplicity is ignored for now
	// TODO add multipclity dependence

	// protection
	if (!particle) {
		return 1.0;
	}
	if (!event) {
		return 1.0;
	}

	// protection
	// apply correction only for mid-rapidity eta<1.5 and charged particles
	if (TMath::Abs(particle->Eta()) > 1.5) {
		return 1.0;
	}
	if (particle->Charge() == 0) {
		return 1.0;
	}

	// get all particle id, prodcution type and pt
	Double_t mcpt = particle->Pt();
	ProductionType prod = kUnknown;

	if (event->IsSecondaryFromMaterial(particle->GetLabel())) {
		prod = kSecDecay;
	}
	if (event->IsSecondaryFromWeakDecay(particle->GetLabel())) {
		prod = kSecMaterial;
	}
	if (event->IsPhysicalPrimary(particle->GetLabel())) {
		prod = kPrim;
	}

	ParticleType ptype = ParticleTypeFromPDG(particle->PdgCode());

	// for now use hard coded values
	if (periodindex == 1) {
		// period = pp 5TeV, LHC17pq
		if (prod == kSecMaterial || prod == kSecDecay) {
			if (systflag == 0) {
				if (!fGsscaleB) {
					Double_t x[17] = {0.05, 0.125, 0.175, 0.225, 0.275, 0.35,
						0.45, 0.55,  0.65,  0.8,   1,     1.2,
						1.4,  1.75,  3.5,   27.5,  125};
					Double_t y[17] = {
						1.07631, 1.07631, 1.07631, 1.14635, 1.22387, 1.28186,
						1.37162, 1.40387, 1.42242, 1.36301, 1.35462, 1.38724,
						1.39539, 1.46115, 1.48689, 1.48689, 1.48689};
					fGsscaleB = new TGraph(17, x, y);
				}
				return fGsscaleB->Eval(mcpt);
			}
			if (systflag == 1) {
				if (!fGsscaleB1) {
					Double_t x[17] = {0.05, 0.125, 0.175, 0.225, 0.275, 0.35,
						0.45, 0.55,  0.65,  0.8,   1,     1.2,
						1.4,  1.75,  3.5,   27.5,  125};
					Double_t y[17] = {
						1.11695, 1.11695, 1.11695, 1.18939, 1.28676, 1.36386,
						1.45975, 1.46592, 1.46357, 1.39391, 1.38146, 1.43343,
						1.4839,  1.56826, 1.60185, 1.60185, 1.60185};
					fGsscaleB1 = new TGraph(17, x, y);
				}
				return fGsscaleB1->Eval(mcpt);
			}
			if (systflag == -1) {
				if (!fGsscaleB2) {
					Double_t x[17] = {0.05, 0.125, 0.175, 0.225, 0.275, 0.35,
						0.45, 0.55,  0.65,  0.8,   1,     1.2,
						1.4,  1.75,  3.5,   27.5,  125};
					Double_t y[17] = {
						1.04527, 1.04527, 1.04527, 1.09517, 1.16257, 1.19152,
						1.289,   1.37694, 1.39118, 1.34851, 1.2972,  1.30618,
						1.26338, 1.30091, 1.30519, 1.30519, 1.30519};
					fGsscaleB2 = new TGraph(17, x, y);
				}
				return fGsscaleB2->Eval(mcpt);
			}
		}
	} else {
		// preiod pp 13 TeV, LHC18b
		// this is the default
		if (prod == kSecMaterial || prod == kSecDecay) {
			if (systflag == 0) {
				if (!fGsscale) {
					Double_t x[17] = {0.05, 0.125, 0.175, 0.225, 0.275, 0.35,
						0.45, 0.55,  0.65,  0.8,   1,     1.2,
						1.4,  1.75,  3.5,   27.5,  125};
					Double_t y[17] = {
						1.04546, 1.04546, 1.04546, 1.10722, 1.17168, 1.20625,
						1.26868, 1.30039, 1.31913, 1.27285, 1.25653, 1.28177,
						1.29794, 1.35358, 1.38495, 1.38495, 1.38495};
					fGsscale = new TGraph(17, x, y);
				}
				return fGsscale->Eval(mcpt);
			}
			if (systflag == 1) {
				if (!fGsscale1) {
					Double_t x[17] = {0.05, 0.125, 0.175, 0.225, 0.275, 0.35,
						0.45, 0.55,  0.65,  0.8,   1,     1.2,
						1.4,  1.75,  3.5,   27.5,  125};
					Double_t y[17] = {
						1.08392, 1.08392, 1.08392, 1.15106, 1.21342, 1.28518,
						1.34486, 1.3514,  1.36498, 1.30593, 1.2892,  1.34095,
						1.39871, 1.46542, 1.51347, 1.51347, 1.51347};
					fGsscale1 = new TGraph(17, x, y);
				}
				return fGsscale1->Eval(mcpt);
			}
			if (systflag == -1) {
				if (!fGsscale2) {
					Double_t x[17] = {0.05, 0.125, 0.175, 0.225, 0.275, 0.35,
						0.45, 0.55,  0.65,  0.8,   1,     1.2,
						1.4,  1.75,  3.5,   27.5,  125};
					Double_t y[17] = {
						1.00192, 1.00192, 1.00192, 1.06146, 1.12252, 1.1249,
						1.20334, 1.2597,  1.27133, 1.23842, 1.19989, 1.18686,
						1.15164, 1.18588, 1.19261, 1.19261, 1.19261};
					fGsscale2 = new TGraph(17, x, y);
				}
				return fGsscale2->Eval(mcpt);
			}
		}
	}

	// internally use the dummy function
	return MCScalingFactor(prod, ptype, mcpt);
}

//____________________________________________________________________________

/// Retrieve the scaling factor for MC primaries and secondaries
///
/// WARNING! only dummy function - use only for testing and for
///
/// This method looks up the proper scaling factors and returns them
/// for offline and -- especially -- online use
///
/// For now this is only a dummy until proper implementation of the
/// functionality is available
///
///
/// \param prod MC production type (primary, secondary from material or
/// secondary from decay) \param part MC particle species (pion, proton, kaon,
/// sigma+, sigma-, xi-, omega-, electron, muon, other) \param pt   MC pT of the
/// particle
///
/// \return scaling factor accoring to the supplied arguments

Double_t AlidNdPtTools::MCScalingFactor(ProductionType prod, ParticleType part,
		Double_t pt) {
	// if prod or part type not set return scaling of 1
	if (prod == kUnknown || part == kUndefined)
		return 1.0;
	// dummy function for testing, scaling sigmas up by a factor two
	// secondaries from decay by factor 1.5
	// and leaves the rest unchanged
	// TODO this should call the ALiMCSpectra weights once they are ready
	// TODO and the corresponding solution for secondariy scaling
	if (prod == kSecDecay || (prod == kPrim))
		return 1.0;
	if (prod == kPrim) {
		// TODO add AliMCSpectraWeights here
		if ((part == kSigmaP) || (part == kSigmaM))
			return 1.0;
	}
	// in other case return scaling factor 1
	return 1.0;
}

//____________________________________________________________________________

/// Convert the PDG code to a ParticleType
///
/// \param pdgCode PDG code
///
/// \return AlidNdPtTools::ParticleType (pion, proton, kaon, sigma+, sigma-,
/// xi-, omega-, electron, muon, other)

AlidNdPtTools::ParticleType AlidNdPtTools::ParticleTypeFromPDG(Int_t pdgCode) {
	if (pdgCode == kElectron || pdgCode == kPositron) {
		return kEl;
	}
	if (pdgCode == kMuonMinus || pdgCode == kMuonPlus) {
		return kMu;
	}
	if (pdgCode == kPiPlus || pdgCode == kPiMinus) {
		return kPi;
	}
	if (pdgCode == kProton || pdgCode == kProtonBar) {
		return kPr;
	}
	if (pdgCode == kKPlus || pdgCode == kKMinus) {
		return kKa;
	}
	if (pdgCode == kSigmaPlus || pdgCode == kSigmaBarMinus) {
		return kSigmaP;
	}
	if (pdgCode == kSigmaMinus || pdgCode == kSigmaBarPlus) {
		return kSigmaM;
	}
	if (pdgCode == kXiMinus || pdgCode == kXiPlusBar) {
		return kXi;
	}
	if (pdgCode == kOmegaMinus || pdgCode == kOmegaPlusBar) {
		return kOmega;
	}
	return kOther;
}

//____________________________________________________________________________
