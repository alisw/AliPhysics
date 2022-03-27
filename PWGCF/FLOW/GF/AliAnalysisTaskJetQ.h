#ifndef AliAnalysisTaskJetQ_H
#define AliAnalysisTaskJetQ_H
#define C_PI_HALF 1.5707963
#define C_PI_TH 4.7123890
#define C_TWOPI 6.2831853
#include "TH3D.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TList.h"
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliAODTrack.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliMultSelection.h"
#include "AliAODInputHandler.h"
#include "AliMCEvent.h"
#include "AliMultSelection.h"
#include "AliVEvent.h"
#include "TList.h"
#include "TMath.h"
#include "AliEventPoolManager.h"
#include "AliBasicParticle.h"
#include "TAxis.h"
using namespace std;
class AliAnalysisTaskJetQ : public AliAnalysisTaskSE
{
    public:
        AliAnalysisTaskJetQ();
        AliAnalysisTaskJetQ(const char *name);
        virtual ~AliAnalysisTaskJetQ();
        virtual void UserCreateOutputObjects();
        virtual void UserExec(Option_t* option);
        virtual void Terminate(Option_t* option);
        void SetCentralityBins(Int_t nBins, Double_t *bins) {setupAxis(nBins,bins,fCentBins,fCentAxis); };
        void SetVtxZBins(Int_t nBins, Double_t *bins) {setupAxis(nBins,bins,fVzBins,fVzAxis); };
        void SetPtBins(Int_t nBins, Double_t *bins) {setupAxis(nBins,bins,fPtBins,fPtAxis); };
        void SetEventMixingCapacity(Int_t nTotEv, Int_t nTotTr, Int_t frReady, Int_t nMinEv) { fEvMixPars[0]=nTotEv; fEvMixPars[1]=nTotTr; fEvMixPars[2]=frReady; fEvMixPars[3]=nMinEv; }; //Fraction is given in %, so should be an integer number!
    private:
        AliAnalysisTaskJetQ(const AliAnalysisTaskJetQ&); // not implemented
        AliAnalysisTaskJetQ& operator=(const AliAnalysisTaskJetQ&); // not implemented
        Bool_t CheckTrigger(Double_t);
        Bool_t AcceptAOD(AliAODEvent*, Double_t lvtxXYZ[3]);
        Int_t FindGivenPt(const Double_t &ptMin, const Double_t &ptMax);
        Int_t FillCorrelations(Int_t &triggerIndex, const Double_t &ptAsMin, const Double_t &ptAsMax);
        Int_t FillMixedEvent(Int_t &triggerIndex, AliEventPool *inpool);
        void fill2DHist(TH1 *&inh, Double_t &xval, Double_t &yval) { ((TH2*)inh)->Fill(xval,yval); };
        void fill3DHist(TH1 *&inh, Double_t &xval, Double_t &yval, Double_t &zval) { ((TH3*)inh)->Fill(xval,yval,zval); };

        void fixPhi(Double_t &inPhi) { if(inPhi<-C_PI_HALF) inPhi+=C_TWOPI; else if(inPhi>C_PI_TH) inPhi-=C_TWOPI; };
        void setupAxis(Int_t &nBins, Double_t *&bins, vector<Double_t> &binCont, TAxis *&ax) {if(binCont.size()>0) binCont.clear(); for(Int_t i=0;i<=nBins;i++) binCont.push_back(bins[i]); if(ax) delete ax; ax = new TAxis(nBins, bins); };
        AliAODEvent *fAOD; //! do not store
        AliMCEvent *fMCEvent; //! do not store
        AliEventPoolManager *fPoolMgr; //! do not store
        TObjArray *fPoolTrackArray; //! do not store
        UInt_t fTriggerType;
        TList *fOutList;
        TAxis *fCentAxis;
        TAxis *fVzAxis;
        TAxis *fPtAxis;
        TH2D *fNormCounter; //!
        TH1 *fCorrPlot; //!
        TH1 *fMixCorrPlot; //!
        vector<Double_t> fCentBins;
        vector<Double_t> fVzBins;
        vector<Double_t> fPtBins;
        Bool_t fPtDif;
        AliEventCuts fEventCuts;
        Int_t fEvMixPars[4];
        ClassDef(AliAnalysisTaskJetQ, 1);
};
#endif
