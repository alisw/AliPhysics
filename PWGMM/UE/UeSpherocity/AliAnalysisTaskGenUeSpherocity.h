/**  
 *  
 *  \class AliAnalysisTaskGenUeSpherocity
 *   
 *  This macro is intended to provide pT spectra vs multiplicity (two estimators) and spherocity (pT weighted and modified version proposed by Lund) author: Antonio Ortiz, original version by Gyula Bencedi 
 *  \date: 	June 4, 2017
 * 
 */
#ifndef ALIANALYSISTASKGENUESPHEROCITY_H
#define ALIANALYSISTASKGENUESPHEROCITY_H

// ROOT includes

#include <TObject.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TH1F.h>
#include <THnSparse.h>
#include <TArrayD.h>
#include <TParticle.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TObjArray.h>

// AliRoot includes

#include <AliAnalysisTaskSE.h>
#include <AliMCEvent.h>
#include <AliStack.h>
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliInputEventHandler.h"

class TList;
class TH1F;
class TH1I;
class TH2D;
class TParticle;
class AliStack;
class AliVVertex;
class AliVParticle;

class AliAnalysisTaskGenUeSpherocity : public AliAnalysisTaskSE { //

	public:

		AliAnalysisTaskGenUeSpherocity();
		AliAnalysisTaskGenUeSpherocity(const char *name);

		virtual ~AliAnalysisTaskGenUeSpherocity();

		virtual void UserCreateOutputObjects();
		virtual void Init();
		virtual void UserExec(Option_t *option);
		virtual void Terminate(Option_t *);
		virtual void SetYRange(Float_t y){ fY=y; }
		virtual void SetGenerator(TString generator){fGenerator=generator;}
		virtual void SetMinPtLeading(Double_t minptl){fMinPtLeading=minptl;}

	private:

		Int_t   GetIndexLeading(Bool_t fIsPseudoRec);
		Int_t   GetMultipliciy(Bool_t fIsPseudoRec, std::vector<Int_t> &mult, std::vector<Float_t> &ptArray,  std::vector<Float_t> &etaArray, std::vector<Float_t> &phiArray);
		void    MakeRTAnalysis(Bool_t fIsPseudoRec);
		void    MakeUeSoNch08Analysis(std::vector<Int_t> &mult);
		void 	ParticleSel(Bool_t fIsPseudoRec, const std::vector<Int_t> &mult);
		Float_t GetSpherocity(
				Int_t nch_so,
				const std::vector<Float_t> &pt,
				const std::vector<Float_t> &eta,
				const std::vector<Float_t> &phi,
				const Bool_t isPtWeighted
				);
		Int_t	GetPidCode(Int_t pdgCode) const;
		Bool_t  IsGoodTrack(Int_t pid, Int_t pdgCode,Double_t pt);
		void    MakeAnaGen(Int_t nch_gen, std::vector<Int_t> &multArray, std::vector<Float_t> &ptArray,  std::vector<Float_t> &etaArray, std::vector<Float_t> &phiArray);
		void    MakeAnaRec(Int_t nch_gen, std::vector<Int_t> &multArray, std::vector<Float_t> &ptArray,  std::vector<Float_t> &etaArray, std::vector<Float_t> &phiArray);
		Int_t    GetMultBin(Bool_t fIsPseudoRec,Int_t mult_int, Int_t mult_select);
		Double_t DeltaPhi(Double_t phia, Double_t phib,
				Double_t rangeMin = -TMath::Pi()/2,
				Double_t rangeMax = 3*TMath::Pi()/2);

	protected:

		Bool_t IsMCEventSelected(TObject* obj);

		AliMCEvent*              fMcEvent;    //!<! MC event
		AliInputEventHandler*    fMcHandler;  //!<! MCEventHandler

		AliStack*   fStack;
		TString     fGenerator;
		Int_t       fIndexLeadingGen;
		Int_t       fIndexLeadingRec;
		Double_t    fMinPtLeading;
		Float_t     fSizeStep;
		//Int_t       fNso_gen;
		//Int_t       fNso_rec;
		Float_t     fspherocity_gen_ptWeighted;
		Float_t     fspherocity_gen;
		Float_t     fspherocity_rec_ptWeighted;
		Float_t     fspherocity_rec;
		Int_t       fbinPerc0_gen;
		Int_t       fbinPerc1_gen;
		Int_t       fbinPerc2_gen;

		Int_t       fbinPerc0_rec;
		Int_t       fbinPerc1_rec;
		Int_t       fbinPerc2_rec;

		Float_t	fY;     	///< rapidity cut

		void FillHisto(const char* objkey, Double_t x);
		void FillHisto(const char* objkey, Double_t x, Double_t y);

		template <class T> T* InitHisto(const char* hname = "hname", const char* htitle = "htitle", Int_t nxbins = 100, Double_t xmin = 0., Double_t xmax = 20., const char* xtitle = "xtitle", const char* ytitle = "ytitle");
		template <class T> T* InitHisto(const char* hname = "hname", const char* htitle = "htitle", Int_t nxbins = 100, Double_t xmin = 0., Double_t xmax = 20., Int_t nybins = 100, Double_t ymin = 0., Double_t ymax = 20., const char* xtitle = "xtitle", const char* ytitle = "ytitle");


		TH1I * fHistEvt;	 	//!<! 	QA of event properties
		TH1F * fHistEta;
		TH1I* fHistPart;	 	//!<! 	QA of particle properties
		TH1D * fDphiNS;
		TH1D * fDphiAS;
		TH1D * fDphiTS;
		TH1D * fMultTS;
		TH1D * fDphiNSRec;
		TH1D * fDphiASRec;
		TH1D * fDphiTSRec;
		TH1D * fMultTSRec;
		TH2D * fSoWeighedVsNchPtL;

		TH2D * fHistPtVsNchNS[11];
		TH2D * fHistPtVsNchAS[11];
		TH2D * fHistPtVsNchTS[11];
		TH2D * fHistPtVsNchNSRec[11];
		TH2D * fHistPtVsNchASRec[11];
		TH2D * fHistPtVsNchTSRec[11];

		TH1F * fHistPt[11];
		TH1F * fHistPtRec[11];
		TH2D * fHistPtVsSoPtW0[11][7];
		TH2D * fHistPtVsSoPtWRec0[11][7];
		TH2D * fHistPtVsSoPtW1[11][7];
		TH2D * fHistPtVsSoPtWRec1[11][7];
		TH2D * fHistPtVsSoPtW2[11][7];
		TH2D * fHistPtVsSoPtWRec2[11][7];

		TH2D * fHistPtVsSo0[11][7];
		TH2D * fHistPtVsSoRec0[11][7];
		TH2D * fHistPtVsSo1[11][7];
		TH2D * fHistPtVsSoRec1[11][7];
		TH2D * fHistPtVsSo2[11][7];
		TH2D * fHistPtVsSoRec2[11][7];

		TH2D * fMult[3];
		TH2D * fSoVsNch[3];
		TH2D * fSoWeighedVsNch[3];
		TH1F * fNch[3];
		TH1F * fNchSoSel[3];
		TH2D * fMultRec[3];
		TH2D * fSoVsNchRec[3];
		TH2D * fSoWeighedVsNchRec[3];
		TH1F * fNchRec[3];
		TH1F * fNchSoSelRec[3];

		TH2D * hDphiSoIS[122];
		TH2D * hPtVsSoIS[122];
		TH2D * hPtVsSoNS[122];
		TH2D * hPtVsSoAS[122];
		TH2D * hPtVsSoTS[122];

		TList* fListOfObjects;	//!<! Output list of objects


		AliAnalysisTaskGenUeSpherocity(const AliAnalysisTaskGenUeSpherocity&);            // not implemented
		AliAnalysisTaskGenUeSpherocity& operator=(const AliAnalysisTaskGenUeSpherocity&); // not implemented

		ClassDef(AliAnalysisTaskGenUeSpherocity, 1);	// Analysis task for LF spectra analysis  
};

#endif
