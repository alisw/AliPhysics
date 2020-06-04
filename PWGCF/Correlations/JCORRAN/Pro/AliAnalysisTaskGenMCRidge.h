#ifndef ALIANALYSISTASKGENMCRIDGE_H
#define ALIANALYSISTASKGENMCRIDGE_H

using namespace std;

class AliMultSelection;
class AliVMultiplicity;
class TClonesArray;
class AliJJetTask;
class AliDirList;
class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;
class AliEmcalTrackSelection;
class AliAnalysisUtils;
class AliCalorimeterUtils;
class AliMultSelection;
class TLorentzVector;

#include "TFile.h"
#include <TSystem.h>
#include <TGrid.h>
#include "THnSparse.h"
#include "AliAnalysisTaskSE.h"
#include "AliTriggerAnalysis.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliPIDCombined.h"
#include "THistManager.h"
#include "AliAODMCParticle.h"
#include <deque>
#include <iostream>
#include <fstream>
#include "AliAnalysisUtils.h"
#include <vector>
#include <TVector.h>
#include <TRandom.h>
#include <TString.h>
#include <TLorentzVector.h>
#include "AliJJetTask.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "AliJetContainer.h"
#include "AliInputEventHandler.h"

class AliAnalysisTaskGenMCRidge : public AliAnalysisTaskEmcalJet {

	pubilc:

	typedef std::vector<Bool_t> Bool_1d;
	typedef std::vector<Double_t> Double1D;
	typedef std::vector<int>      Int1D;
	typedef std::vector<Double1D> Double2D;
	typedef std::vector<TLorentzVector>   TLorentzVector1D;



        AliAnalysisTaskGenMCRidge();
        AliAnalysisTaskGenMCRidge
        (
              const char *name
              , const char *option
        );

        AliAnalysisTaskGenMCRidge
          (
           const AliAnalysisTaskGenMCRidge& ap
          );

        AliAnalysisTaskGenMCRidge& operator =
        (
              const AliAnalysisTaskGenMCRidge& ap
        );

        ~AliAnalysisTaskGenMCRidge();

        virtual void    UserCreateOutputObjects();
        virtual void    Exec(Option_t *);
        virtual void    FinishTaskOutput();
        virtual void    Terminate(Option_t *);


	void RhoSparse(AliJetContainer *ktContainer, AliJetContainer *aktContainer, Int_t numberofexcludingjets);
	Bool_t isOverlapping(AliEmcalJet* jet1, AliEmcalJet* jet2);


	private:

	TAxis                           binCent; //! 
	TAxis                           binTPt; //!
	TAxis                           binAPt; //!
	TAxis                           binPhi; //!
	TAxis                           binEta; //!
 	TAxis                           binLHPt; //!
	TAxis                           binJetpT; //!

	AliInputEventHandler*           fMcHandler;
	AliMCEvent*              	fMcEvent;
	AliStack*			fStack;

	ClassDef(AliAnalysisTaskGenMCRidge, 1);
}
#endif
