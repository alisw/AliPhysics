#ifndef AliAnalysisTaskInclusivef0f2_H
#define AliAnalysisTaskInclusivef0f2_H

#include "THnSparse.h"
#include "AliAnalysisTaskSE.h"
#include "AliTriggerAnalysis.h"
#include "AliESDtrackCuts.h"
#include "THistManager.h"
#include "AliPIDResponse.h"
#include <deque>
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliPIDCombined.h"
#include "AliPID.h"
#include "AliStack.h"
#include "AliEventCuts.h"


class AliMultSelection;
class AliVMultiplicity;
class AliStack;

class AliAnalysisTaskInclusivef0f2RunTable {
    public:
        enum {kPP,kPA,kAA,kUnknownCollType};
        AliAnalysisTaskInclusivef0f2RunTable();
        AliAnalysisTaskInclusivef0f2RunTable(Int_t runnumber);
        ~AliAnalysisTaskInclusivef0f2RunTable();

        Bool_t IsPP(){
            return fCollisionType==kPP;
        }
        Bool_t IsPA(){
            return fCollisionType==kPA;
        }
        Bool_t IsAA(){
            return fCollisionType==kAA;
        }
	void SetColl(int i){
		fCollisionType = i;
	}
    private:
        Int_t  fCollisionType; //! Is proton-proton collisions?
};
class AliAnalysisTaskInclusivef0f2 : public AliAnalysisTaskSE{
 public:
	typedef std::vector<Double_t> Double1D;

	AliAnalysisTaskInclusivef0f2();
        AliAnalysisTaskInclusivef0f2(const char *name, const char *option);
	~AliAnalysisTaskInclusivef0f2();

	AliAnalysisTaskInclusivef0f2(const AliAnalysisTaskInclusivef0f2& ap);
	AliAnalysisTaskInclusivef0f2& operator=(const AliAnalysisTaskInclusivef0f2& ap); 
        virtual void UserCreateOutputObjects();
        virtual void UserExec(Option_t* option);
	virtual void FinishTaskOutput();
        virtual void Terminate(Option_t *);
	bool GoodTracksSelection(int iTrackCut, double TPCsig, double TOFsig, double TPCalonesig);
        void FillTracks();
	int GetPID(AliPIDResponse *pid, const AliVTrack *trk);

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

	void SetIsMC (Bool_t ismc) {IsMC = ismc;}

	AliEventCuts fEventCuts;  // Event cuts

 private:
        typedef std::vector<AliVTrack*> tracklist;
        typedef std::deque<tracklist> eventpool;
        typedef std::vector<std::vector<eventpool> > mixingpool;
	typedef std::vector<mixingpool> mixingpooltrkcut;
	
        mixingpool                      fEMpool; //!
	mixingpooltrkcut		fEMpooltrk;

	AliAnalysisTaskInclusivef0f2RunTable*   fRunTable=nullptr; //!
        AliAODEvent*            	fAOD;           //! input event
	TString                         fOption;
        TList*                  	fOutput=nullptr;    //! output list
        TH1F*                   	fHistPt;        //! dummy histogram
	AliTriggerAnalysis*             fTrigger=nullptr; //!
	AliESDtrackCuts*                fTrackCuts=nullptr; //!
	THistManager*                   fHistos=nullptr; //!
	Bool_t                          IsFirstEvent=kTRUE;
	AliVEvent*                      fEvt=nullptr; //!
	AliPIDResponse                 *fPIDResponse=nullptr; //!
	AliPIDCombined                 *fPIDCombined=nullptr;
	std::vector<UInt_t>             goodtrackindices; //!
	UInt_t                          fFilterBit;
	Int_t                           fParticleType;
        AliMultSelection               *sel=nullptr;//! 
        Int_t                           bookingsize = 20;
        AliVMultiplicity*               fMultiplicity=nullptr;//!
	TClonesArray*                   fMCArray=nullptr;
	bool				IsMC=kFALSE;
	AliStack*			fmcstack=nullptr;
        TAxis                           binCent; //! 
	TAxis				binCentForMC;
        TAxis                           binZ; //!
        TAxis                           binPt; //!
	TAxis                           binType; //!
        TAxis                           binMass; //!
	TAxis				binPtTrack;
	TAxis				binPhiTrack;
	TAxis				binEtaTrack;
	TAxis				binCharge;
	TAxis				binTrackCutBit;
	TAxis				binPID;
	TAxis				binExKaonNum;
	TAxis				binTrackPt;
	TAxis				binSigma;
	TAxis				binEta;
	TAxis				binSwitch;

        Double_t                        fCent=-1;
        Double_t                        fZ=-30;
        Double_t                        fptcut = 0.15;
        Double_t                        fetacut = 0.8;
        Int_t                           centbin = -1 ;
        Int_t                           zbin = -1 ;
	Int_t				trkbin = -1;
	Int_t				pidbin = -1;

	std::vector<Double_t>		EffpT;
	Int_t fEff_npT_step = 140;
	Double_t fEff_pT_max = 14.0;

        ClassDef(AliAnalysisTaskInclusivef0f2, 1);
};

#endif

