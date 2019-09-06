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

		void SetYRange(Float_t y){ fY=y; }

	private:

		void	EventSel(TObject* obj);
		Int_t       GetMultipliciy(std::vector<Int_t> &mult, std::vector<Float_t> &ptArray,  std::vector<Float_t> &etaArray, std::vector<Float_t> &phiArray, TObject* obj);
		void 	ParticleSel(const std::vector<Int_t> &mult, TObject *obj);
		Float_t     GetSpherocity(
				const std::vector<Float_t> &pt,
				const std::vector<Float_t> &eta,
				const std::vector<Float_t> &phi,
				const Bool_t isPtWeighted
				);
		Int_t	GetPidCode(Int_t pdgCode) const;

	protected:

		Bool_t IsMCEventSelected(TObject* obj);

		AliMCEvent*              fMcEvent;    //!<! MC event
		AliInputEventHandler*    fMcHandler;  //!<! MCEventHandler

		AliStack* 	fStack;
		Float_t     fSizeStep;
		Int_t       fNrec;
		Float_t	fY;     	///< rapidity cut

		void FillHisto(const char* objkey, Double_t x);
		void FillHisto(const char* objkey, Double_t x, Double_t y);

		template <class T> T* InitHisto(const char* hname = "hname", const char* htitle = "htitle", Int_t nxbins = 100, Double_t xmin = 0., Double_t xmax = 20., const char* xtitle = "xtitle", const char* ytitle = "ytitle");
		template <class T> T* InitHisto(const char* hname = "hname", const char* htitle = "htitle", Int_t nxbins = 100, Double_t xmin = 0., Double_t xmax = 20., Int_t nybins = 100, Double_t ymin = 0., Double_t ymax = 20., const char* xtitle = "xtitle", const char* ytitle = "ytitle");



		TH1I* fHistEvt;	 	//!<! 	QA of event properties
		TH1I* fHistPart;	 	//!<! 	QA of particle properties

		TList* fListOfObjects;	//!<! Output list of objects


		AliAnalysisTaskGenUeSpherocity(const AliAnalysisTaskGenUeSpherocity&);            // not implemented
		AliAnalysisTaskGenUeSpherocity& operator=(const AliAnalysisTaskGenUeSpherocity&); // not implemented

		ClassDef(AliAnalysisTaskGenUeSpherocity, 1);	// Analysis task for LF spectra analysis  
};

#endif
