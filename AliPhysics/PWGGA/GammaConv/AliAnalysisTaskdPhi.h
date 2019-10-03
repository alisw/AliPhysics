////////////////////////////////////////////////
//--------------------------------------------- 
// Class doing conversion gamma dPhi correlations
// Gamma Conversion analysis
//---------------------------------------------
////////////////////////////////////////////////

#ifndef AliAnalysisTaskdPhi_cxx
#define AliAnalysisTaskdPhi_cxx

#include "AliAnalysisTaskSE.h"

#include <TAxis.h>
#include <TH3I.h>
#include <THnSparse.h>
#include <THn.h>

//#include <AliAnalysisFilter.h>
#include <iostream>
//#include <AliAnaConvCorrBase.h>
#include <AliLog.h>
#include <AliAnalysisCuts.h>
//class AliAnaConvCorrPion;
class AliConvEventCuts;
class AliConversionPhotonCuts;
class AliConversionMesonCuts;
class AliV0ReaderV1;
class TList;
class TH2I;
//class THnSparseF;

using namespace std;

class AliAnalysisTaskdPhi : public AliAnalysisTaskSE {
  
public:
	AliAnalysisTaskdPhi(const char *name="slindal_dPhi");
	virtual ~AliAnalysisTaskdPhi();
	
	virtual void   UserCreateOutputObjects();
	//virtual void   SetUpCorrObjects();
	//  virtual void   SetUpCorrAxes(AliAnaConvCorrBase * corr);
	virtual void   SetUpBins();
	virtual void   UserExec(Option_t *option);
	virtual void   Terminate(Option_t *);
	
	TAxis& GetAxistPt()     { return fAxistPt;      }
	TAxis& GetAxiscPt()     { return fAxiscPt;      }
	TAxis& GetAxisdEta()    { return fAxisdEta;     }
	TAxis& GetAxisTrigEta() { return fAxisTrigEta;  }
	TAxis& GetAxisAssEta()  { return fAxisAssEta;   }
	TAxis& GetAxisPhi()     { return fAxisdPhi;     }
	TAxis& GetAxisZ()       { return fAxisZ;        }
	TAxis& GetAxisCent()    { return fAxisCent;     }
	TAxis& GetAxisPiMass()  { return fAxisPiM;      }
	
	void SetV0Filter(AliConvEventCuts * filterEvent, AliConversionPhotonCuts * filter) { fV0FilterEvent=filterEvent, fV0FilterPhoton = filter; }
	void AddEventFilter(TObject * filter, Bool_t high = kTRUE) { fEventFilters[high].AddLast(filter); }
	void AddV0Filter(TObject * filter, Bool_t high = kTRUE) { fV0Filters[high].AddLast(filter); }
	void AddMesonFilter(TObject * filter, Bool_t high = kTRUE) { fMesonFilters[high].AddLast(filter); }
	void AddTrackFilter(TObject * filter, Bool_t high = kTRUE) { fTrackFilters[high].AddLast(filter); }
	void SetCorrectionMap(THnF * map)  { fCorrectionMap = map; }

	void SetEventFilter(AliConvEventCuts * filter) { fEventFilter = filter; }
	void SetMesonFilter(AliConversionMesonCuts * filter) { fMesonFilter = filter; }
	void SetPhotonFilter(AliConversionPhotonCuts * filter) { fPhotonFilter = filter; }
	void SetV0Reader(AliV0ReaderV1 * reader) { fV0Reader = reader; }
	void SaveReaderHists(Bool_t save = kTRUE) { fSaveReaderHists = save; }
	AliAnalysisCuts * GetTrackCuts() const { return fTrackFilter; }
	void SetTrackFilter( AliAnalysisCuts * cuts) { if (fTrackFilter) delete fTrackFilter; fTrackFilter = cuts; }
  
protected:
	
	TClonesArray * GetConversionGammas(Bool_t isAOD);
  
private:

	//void CorrelateWithTracks(AliAODConversionParticle * particle, TObjArray tracks[], Int_t ntrackfilters, Bool_t ** lowtrackmap, Int_t nltf, Int_t const tIDs[4], Double_t dphiValues[]);
	//void FillCounters(TObjArray * particles, TObjArray tracks[], Int_t ntrackfilters, Float_t cent, Float_t vtxz);

	Double_t GetTrackCorrection(Double_t vtxz, AliVTrack * track);

	
	///Get the distance in phi between trigger particle and correlated particle
	Float_t GetDPhi(Float_t dPhi) { 
		if ( dPhi < 3*TMath::PiOver2() && dPhi > - TMath::PiOver2() ) return dPhi;
		else return ( (dPhi>0)? dPhi - TMath::TwoPi() : dPhi + TMath::TwoPi() ); 
	}

	THnSparseF * CreateSparse(TString nameString, TString titleString, TList * axesList);
	Int_t GetBin(TAxis &axis, Double_t value);
	THnSparseF * GetMEHistogram(Int_t binz, Int_t binc, TObjArray * array);
	//  AliAnaConvCorrBase * GetCorrObject(Int_t binz, Int_t binc, TObjArray * array);
	void Process(TObjArray * gammas, TObjArray * tracks, Float_t cent, Float_t vtxz);
	void FindDeltaAODBranchName(AliVEvent * event);
	

	///Members
	TList 							*fHistograms; 					//histograms
	THnSparseF 						*fCorrSparse;
	THnSparseF 						*fTrigSparse;
	THnSparseF 						*fTrackSparse;
	THnSparseF 						*fMassSparse;

	AliV0ReaderV1 					*fV0Reader; 					// V0 reader
	Bool_t 							fSaveReaderHists; 				// save histograms from v0 reader
	AliConvEventCuts				*fV0FilterEvent; 				// additional v0 filter on top of v0 reader
	AliConversionPhotonCuts			*fV0FilterPhoton; 				// additional v0 filter on top of v0 reader
	TObjArray 						fV0Filters[2]; 					// Array of v0 filters, increasingly loose ideally.
	TObjArray 						fEventFilters[2]; 				// Array of event filters, increasingly loose ideally.
	AliConvEventCuts 				*fEventFilter;	 				// additional v0 filter for events only
	AliConversionPhotonCuts 		*fPhotonFilter; 				// additional v0 filter for photons only
	AliConversionMesonCuts 			*fMesonFilter; 					// additional meson filter behind fv0filter
	TObjArray 						fMesonFilters[2]; 				// Array of Meson filters
	AliAnalysisCuts 				*fTrackFilter; 					// Cuts for corr tracks
	TObjArray 						fTrackFilters[2]; 				// Array of track cuts

	TObjArray 						fGammas;
	TObjArray 						fTracks;
	//Int_t 								fTFBin; 						// Axis for track filters

	TH2I 							*hMEvents; 						// event histrogam
	TH2I 							*hTrackCent; 					// event histrogam
	TH3F 							*hTrigPt; 						// trigger pt
	TH2F 							*hTrackPt; 						// track pt
	TH1F 							*hTrigPhi; 						// trigger phi

	//AliAnaConvCorrBase * 				fPhotonCorr;	 				// photon
	//AliAnaConvCorrPion * 				fPionCorr; 						// pion

	TString 						fDeltaAODBranchName; 			// comment

	TAxis 							fAxistPt;						// comment
	TAxis 							fAxiscPt;						// comment
	TAxis 							fAxisdEta;						// comment
	TAxis 							fAxisTrigEta;					// comment
	TAxis 							fAxisAssEta;					// comment
	TAxis 							fAxisdPhi;						// comment
	TAxis 							fAxisCent;						// comment
	TAxis 							fAxisZ;							// comment
	TAxis 							fAxisPiM;						// comment
	TAxis 							fAxisTrackFilters; 				// comment
	TAxis 							fAxisV0Filters; 				// comment
	TAxis 							fAxisMesonFilters; 				// comment
	Bool_t 							fkTrackAxis; 					// on or off
	Bool_t 							fkV0Axis; 						// on or off
	Bool_t 							fkPionAxis; 					// on or off

	TList 							fAxesList; 						// dphi axes list
	TList 							fTrigAxesList; 					// Trigger axes list
	TList 							fTrackAxesList; 				// Trackociated particles axes list
	TList 							fMassAxesList; 					// Mass vs pt sparse

	Bool_t 							fDoPhoton; 						// do photon analysis?
	THnF * 							fCorrectionMap;

	
	AliAnalysisTaskdPhi(const AliAnalysisTaskdPhi&); // not implemented
	AliAnalysisTaskdPhi& operator=(const AliAnalysisTaskdPhi&); // not implemented
	
	ClassDef(AliAnalysisTaskdPhi, 11); 
};

inline THnSparseF * AliAnalysisTaskdPhi::GetMEHistogram(Int_t binz, Int_t binc, TObjArray * array) {
	///Get Mixed Event histogram
	if(binz < 0 || binz > fAxisZ.GetNbins()) {
		cout << "error out of z axis range: " << binz << endl; 
		return NULL;
	}  
	if(binc < 0 || binc >= fAxisCent.GetNbins()) {
		cout << "error out of centraliy axis range: " << binc << endl; 
		return NULL;
	}  
	
	TObjArray * arrayc = static_cast<TObjArray*>(array->At(binc));
	THnSparseF * histogram = static_cast<THnSparseF*>(arrayc->At(binz));
	return histogram;
}


// inline AliAnaConvCorrBase * AliAnalysisTaskdPhi::GetCorrObject(Int_t binz, Int_t binc, TObjArray * array) {
//   ///Get correlation object
//   if(binc < 0 || binz < 0) {
// 	  AliError("We have a bad bin!!!");
// 	  return NULL;
// 	}

//   TObjArray * arrayc = static_cast<TObjArray*>(array->At(binc));
//   AliAnaConvCorrBase * corrmaker = static_cast<AliAnaConvCorrBase*>(arrayc->At(binz));
//   return corrmaker;
// }

inline Int_t AliAnalysisTaskdPhi::GetBin(TAxis & axis, Double_t value) {
	//Return bin - 1 if within range, else return -1
	Int_t bin = axis.FindFixBin(value);
	bin = (bin > 0 && bin <= axis.GetNbins()) ? bin -1 : -1;
	return bin;
}

#endif

