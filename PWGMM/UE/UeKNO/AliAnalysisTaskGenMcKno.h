/**  
 *  
 *  \class AliAnalysisTaskGenMcKno
 *  
 *  Antonio Ortiz (ICN-UNAM), antonio.ortiz@nucleares.unam.mx 
 *  First version: 	April 23, 2020
 * 
 */
#ifndef ALIANALYSISTASKGENMCKNO_H
#define ALIANALYSISTASKGENMCKNO_H

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

class AliAnalysisTaskGenMcKno : public AliAnalysisTaskSE { //

	public:

		AliAnalysisTaskGenMcKno();
		AliAnalysisTaskGenMcKno(const char *name);

		virtual ~AliAnalysisTaskGenMcKno();

		virtual void UserCreateOutputObjects();
		virtual void Init();
		virtual void UserExec(Option_t *option);
		virtual void Terminate(Option_t *);
		virtual void SetPtLeading(Double_t minptl, Double_t maxptl){fMinPtLeading=minptl; fMaxPtLeading=maxptl;}

	private:

		Int_t   GetIndexLeading();
		Int_t   GetMultipliciy(std::vector<Int_t> &mult);
		void    MakeRTAnalysis(std::vector<Int_t> &mult);
		Int_t	GetPidCode(Int_t pdgCode) const;
		Double_t DeltaPhi(Double_t phia, Double_t phib,
				Double_t rangeMin = -TMath::Pi()/2,
				Double_t rangeMax = 3*TMath::Pi()/2);

	protected:

		Bool_t IsMCEventSelected(TObject* obj);

		AliMCEvent*              fMcEvent;    //!<! MC event
		AliInputEventHandler*    fMcHandler;  //!<! MCEventHandler

		AliStack*   fStack;
		Int_t       fIndexLeading;
		Double_t    fMinPtLeading;
		Double_t    fMaxPtLeading;	
		
		TH1D * fDphiNS;
		TH1D * fDphiAS;
		TH1D * fDphiTS;
		TH1D * fMultTS;
		TH1D * fptL;

		TH2D * fHistPtVsNchNS[12][5];
		TH2D * fHistPtVsNchAS[12][5];
		TH2D * fHistPtVsNchTS[12][5];
		TH1D * fMult[5];
		TH1D * fMult2[5];

		TList* fListOfObjects;	//!<! Output list of objects


		AliAnalysisTaskGenMcKno(const AliAnalysisTaskGenMcKno&);            // not implemented
		AliAnalysisTaskGenMcKno& operator=(const AliAnalysisTaskGenMcKno&); // not implemented

		ClassDef(AliAnalysisTaskGenMcKno, 1);	// Analysis task for LF spectra analysis  
};

#endif
