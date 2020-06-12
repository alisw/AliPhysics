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
#include <TParticle.h>
#include "AliJFJTask.h"

class AliAnalysisTaskGenMCRidge : public AliAnalysisTaskEmcalJet {
    public:
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

        double GetMultiplicity( AliStack* stack );
        bool GetProperTracks( AliStack* stack );
        void GetCorrelations();
        void SetPtHardMin( double pthardmin ){fPtHardMin = pthardmin; };
        void SetPtHardMax( double pthardmax ){fPtHardMax = pthardmax; };
        void SetJFJTaskName(TString name){ fJFJTaskName=name; } 
        void ESETagging(int iESE, double lpPT);
        TAxis AxisFix( TString name, int nbin, Double_t xmin, Double_t xmax);
        TAxis AxisVar( TString name, std::vector<Double_t> bin );
        TAxis AxisLog( TString name, int nbin, Double_t xmin, Double_t xmax
            , Double_t xmin0);
        TAxis AxisStr( TString name, std::vector<TString> bin );
        THnSparse * CreateTHnSparse(TString name, TString title
            , Int_t ndim, std::vector<TAxis> bins, Option_t * opt="");
        THnSparse * CreateTHnSparse(TString name, TString title
            , TString templ, Option_t * opt="");
        Long64_t FillTHnSparse( TString name, std::vector<Double_t> x, Double_t w=1.);
        Long64_t FillTHnSparse( THnSparse *h, std::vector<Double_t> x, Double_t w=1.);


    private:
	typedef std::vector < TParticle* > tracklist;
	typedef std::deque<tracklist> eventpool;
	typedef std::vector<vector<eventpool> > mixingpool;

	mixingpool			fEMpool; //!

	TAxis                           binCent; //! 
	TAxis                           binTPt; //!
	TAxis                           binAPt; //!
	TAxis                           binPhi; //!
	TAxis                           binEta; //!
 	TAxis                           binLHPt; //!
	TAxis                           binJetPt; //!
	TAxis				binNtrig; //!
	TAxis				binZ; //!

	AliInputEventHandler*           fMcHandler=nullptr; //!
	AliMCEvent*              	fMcEvent=nullptr; //!
	AliStack*			fStack=nullptr; //!

	AliDirList*                     fOutput=nullptr; //!
	THistManager*                   fHistos=nullptr; //!
	TString                         fOption; 

	std::vector < TParticle* >      goodTracks; //!
	std::vector < Double_t >        NTracksPerPtBin; //!

	Double_t			RHO;

	Double_t			fMult;
	Double_t			fZ;
	Double_t			fLHPt;
	Double_t			fJetPt;

	Int_t				fbookingsize = 7;
  
  double fPtHardMin;
  double fPtHardMax;
  Bool_t TagThisEvent[5];
  AliJFJTask        *fJFJTask;
  TString fJFJTaskName;



    ClassDef(AliAnalysisTaskGenMCRidge, 1);
};

#endif
