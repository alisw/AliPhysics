/**  
 *  
 *  \class AliAnalysisTaskGenRT
 *  
 *  \brief Analysis task to be used for on-the-fly simulations in LEGO trains
 * 
 *  Minimal analysis task meant to be used for on-the-fly simulations in LEGO trains (for generating light flavor particle pt spectra)

 *  \author:	Gyula Bencedi  <Gyula.Bencedi@cern.ch>, WIGNER RCP
 * 
 *  \date: 	June 4, 2017
 * 
 */
#ifndef ALIANALYSISTASKGENRT_H
#define ALIANALYSISTASKGENRT_H

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
class TParticle;
class AliStack;
class AliVVertex;
class AliVParticle;

class AliAnalysisTaskGenRT : public AliAnalysisTaskSE { //

	public:

		AliAnalysisTaskGenRT();
		AliAnalysisTaskGenRT(const char *name);

		virtual ~AliAnalysisTaskGenRT();

		virtual void UserCreateOutputObjects();
		virtual void Init();
		virtual void UserExec(Option_t *option);
		virtual void Terminate(Option_t *);

		void IsXiAllowed(Bool_t allow){ AllowXi = allow; }
		void SetYRange(Float_t y){ fY=y; }
		void SetMinPtLeading(Double_t minpt ){ fMinPtLeading = minpt; }
		void SetMaxPtLeading(Double_t maxpt ){ fMaxPtLeading = maxpt; }

	private:

		Int_t	GetIndexLeading();
		void	MakeRTAnalysis();
		Int_t   GetMultipliciy(std::vector<Int_t> &mult);
		Double_t DeltaPhi(Double_t phia, Double_t phib,
				Double_t rangeMin = -TMath::Pi()/2,
				Double_t rangeMax = 3*TMath::Pi()/2);

		Int_t	GetPidCode(Int_t pdgCode) const;

	protected:

		Bool_t IsMCEventSelected(TObject* obj);

		AliMCEvent*              fMcEvent;    //!<! MC event
		AliInputEventHandler*    fMcHandler;  //!<! MCEventHandler

		AliStack* 	fStack;

		Bool_t AllowXi;
		Float_t	fY;     	///< rapidity cut
		Int_t fIndexLeadingGen;
		Double_t fMinPtLeading;
		Double_t fMaxPtLeading;


		void FillHisto(const char* objkey, Double_t x);

		template <class T> T* InitHisto(const char* hname = "hname", const char* htitle = "htitle", Int_t nxbins = 100, Double_t xmin = 0., Double_t xmax = 20., const char* xtitle = "xtitle", const char* ytitle = "ytitle");

		TH1I* fHistEvt;	 	//!<! 	QA of event properties
		TH1I* fHistPart;	 	//!<! 	QA of particle properties	
		TH1F* hPtLeadingGen;
		TH2F* hMultTSvsPtLeading;
		TH1F* fDphiNS;
		TH1F* fDphiAS;
		TH1F* fDphiTS;
		TH1F* fMultTS;
		TH2F* fHistPtVsNchNS[11];
		TH2F* fHistPtVsNchAS[11];
		TH2F* fHistPtVsNchTS[11];
		TH1F* fMult[6];
		TH1F* fMult2[6];


		TList* fListOfObjects;	//!<! Output list of objects


		AliAnalysisTaskGenRT(const AliAnalysisTaskGenRT&);            // not implemented
		AliAnalysisTaskGenRT& operator=(const AliAnalysisTaskGenRT&); // not implemented

		ClassDef(AliAnalysisTaskGenRT, 1);	// Analysis task for LF spectra analysis  
};

#endif
