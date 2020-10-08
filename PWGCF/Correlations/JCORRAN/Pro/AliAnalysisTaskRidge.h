/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

 // Short comment describing what this class does needed!

 //===========================================================
 // Dummy comment, should be replaced by a real one
 //===========================================================

#ifndef ALIANALYSISTASKRIDGE_H
#define ALIANALYSISTASKRIDGE_H

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

class AliAnalysisTaskRidgeRunTable {
    public:
        enum {kPP,kPA,kAA,kUnknownCollType};
        AliAnalysisTaskRidgeRunTable();
        AliAnalysisTaskRidgeRunTable(Int_t runnumber);
        ~AliAnalysisTaskRidgeRunTable();

        Bool_t IsPP(){
            return fCollisionType==kPP; 
        }
        Bool_t IsPA(){
            return fCollisionType==kPA; 
        }
        Bool_t IsAA(){
            return fCollisionType==kAA; 
        }
    private:
        Int_t  fCollisionType=kPP; //! Is proton-proton collisions?
};


class AliAnalysisTaskRidge : public AliAnalysisTaskEmcalJet {
    public:
      typedef std::vector<Bool_t> Bool_1d;
      typedef std::vector<Double_t> Double1D;
      typedef std::vector<int>      Int1D;
      typedef std::vector<Double1D> Double2D;
      typedef std::vector<TLorentzVector>   TLorentzVector1D;

        enum {  kSD=0, kDD, kND, kCD, kAllProc};
        //PN = unlike sign, PP and NN are like signs

        AliAnalysisTaskRidge();
        AliAnalysisTaskRidge
        ( 
              const char *name
              , const char *option
        );

        AliAnalysisTaskRidge
          (
           const AliAnalysisTaskRidge& ap
          );

        AliAnalysisTaskRidge& operator = 
        (
              const AliAnalysisTaskRidge& ap
        );

        ~AliAnalysisTaskRidge();
        

        virtual void    UserCreateOutputObjects();
        virtual void    Exec(Option_t *);
        virtual void    FinishTaskOutput();
        virtual void    Terminate(Option_t *);

        void SetOption(char * option) {fOption = option;}
        void SetFilterBit(UInt_t filterbit) {fFilterBit = filterbit;}

	void SetEfficiencyFile(char* fname) { TGrid::Connect("alien://"); fefficiencyFile = TFile::Open(fname,"READ"); }
	void SetEfficiency3DFile(char* fname) { TGrid::Connect("alien://"); fefficiency3DFile = TFile::Open(fname,"READ"); }

        Bool_t  GoodTracksSelection(int trk);
        Bool_t  GoodTrackletSelection();
	Bool_t  GoodTracksSelectionMC();

        void FillTracks();
        void FillTracklets();
	void FillTracksMC();

        void SetMixing (Bool_t setmixing) {fsetmixing = setmixing;}
        void SetIsAA (Bool_t isaa) {IsAA = isaa;}
        void SetIsMC (Bool_t ismc) {IsMC = ismc;}
        void SetParticleType(Int_t partype) {fParticleType = partype;}
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

	bool ALICEAccp(double pt, double eta);
	bool CMSAccp(double pt, double eta);

        struct Tracklet{
          Double_t eta;
          Double_t phi;
        };
	struct Particlelet{
		Double_t eta;
		Double_t phi;
		Double_t pt;
		Int_t IsTrackRecon;
	};

	void RhoSparse(AliJetContainer *ktContainer, AliJetContainer *aktContainer, Int_t numberofexcludingjets);

    protected:
	void MeasureBgDensity(AliJetContainer* ktContainer);
	Bool_t isOverlapping(AliEmcalJet* jet1, AliEmcalJet* jet2);
	
    private:
        typedef std::vector<AliVTrack*> tracklist;
        typedef std::deque<tracklist>  eventpool;
        typedef std::vector<vector<eventpool> > mixingpool; 

        typedef std::vector<Tracklet> trackletlist;
        typedef std::deque<trackletlist>  eventpooltracklet;
        typedef std::vector<vector<eventpooltracklet> > mixingpooltracklet; 

	typedef std::vector<Particlelet> trackMClist;	
	typedef std::deque<trackMClist>  eventpoolMC;
	typedef std::vector<vector<eventpoolMC> > mixingpoolMC;

      
        TString                         fOption;
	AliDirList*				fOutput=nullptr; //!

	TFile*				fefficiencyFile= TFile::Open("EffOut.root","read"); //
	TFile*				fefficiency3DFile=nullptr; //

        AliTriggerAnalysis*             fTrigger=nullptr; //!
        AliESDtrackCuts*                fTrackCuts=nullptr; //!
        AliVEvent*                      fEvt=nullptr; //!
	AliJJetTask*			fJetTask=nullptr; //!

        UInt_t                          fFilterBit=0x300;
        Bool_t                          IsFirstEvent=kTRUE;
        AliAnalysisTaskRidgeRunTable*   fRunTable=nullptr; //!

        Float_t                         fCent;
        Double_t                        fZ;
	Double_t			fZ_gen;
 
        std::vector < UInt_t >          goodtrackindices; //!
        std::vector < UInt_t >          goodtrackindicesMCALICE; //!
	std::vector < UInt_t >          goodtrackindicesMCCMS; //!

	std::vector < Double_t > 	NTracksPerPtBin;
	std::vector < Double_t >        NTracksPerPtBinMCALICE;
	std::vector < Double_t >        NTracksPerPtBinMCCMS;

        mixingpool                      fEMpool; //!
        mixingpooltracklet              fEMpooltracklet; //!
	mixingpoolMC			fEMpoolMCALICE; //!
	mixingpoolMC                    fEMpoolMCCMS; //!  

        TAxis                           binCent; //! 
        TAxis                           binZ; //!
	TAxis				binTPt; //!
	TAxis				binAPt; //!
	TAxis				binPhi; //!
	TAxis				binEta; //!
	TAxis				binMCEta; //!
	TAxis				binLtpt; //!
	TAxis				binJetpT; //!
	TAxis				binRho; //!

	TAxis				binUnipT; //!
	TAxis				binTrig; //!
	TAxis				binV0Amp; //!
	TAxis				binTrkEff; //!

	TAxis				binPt; //!
	TAxis				binPt1; //!
	TAxis				binNtrig; //!

	TAxis				binPhiTrack; //!
	TAxis				binEtaTrack; //!

        Int_t                           centbin = -1 ;
        Int_t                           zbin = -1 ;
        Double_t                        fptcut = 0.2;
        Double_t                        fetacut = 0.9;
	Double_t			fLT_pT;
	Double_t			fLT_pT_MCALICE;
	Double_t			fLT_pT_MCCMS;
	Double_t			fJetPt;
	Double_t			fJetPtMC;

	Bool_t                          fsetmixing = kTRUE;
        Bool_t                          IsDGV0=kFALSE;
        Bool_t                          IsDGV0FMD=kFALSE;
        Bool_t                          IsAA=kFALSE;
        Bool_t                          IsMC=kFALSE;
        THistManager*                   fHistos=nullptr; //!
        TClonesArray*                   fMCArray=nullptr; //!
        Int_t                           fParticleType = 99999;
        Int_t                           fNTracks = 0;
        AliMultSelection               *sel=nullptr;//! 
        Int_t                           bookingsize = 7;
        AliVMultiplicity*               fMultiplicity=nullptr;//!
	Int_t				bookingsizeMC = 7;
	std::vector< std::vector< double > > Eff;
	std::vector< std::vector< std::vector< double > > > Eff3D;

	Double1D EffpT;

	Int_t fEff_npT_step = 40;
	Double_t fEff_pT_min = 0.2;
	Double_t fEff_pT_max = 20.0;
	Double_t fEff_pT_l = 0.1;

	Int_t fEff_neta_step = 18;
	Double_t fEff_eta_min = -0.9;
	Double_t fEff_eta_max = 0.9;
	Double_t fEff_eta_l = 0.1;

	Int_t fEff_nphi_step = 180;

        Int_t ITS_fEff_neta_step = 26;
        Double_t ITS_fEff_eta_min = -1.3;
        Double_t ITS_fEff_eta_max = 1.3;
        Double_t ITS_fEff_eta_l = 0.1;

	Double_t AbsZmax = 8.0;
	Double_t V0M_mean;

	Double_t RHO;
	Double_t RHOM;


    ClassDef(AliAnalysisTaskRidge, 1);
};

#endif

