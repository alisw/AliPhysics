#ifndef ALIDNDPTHELPER_H
#define ALIDNDPTHELPER_H

//
// static dNdPt helper functions
//
// Origin: Jan Fiete Grosse-Oetringhaus
// Modified and Extended: Jacek Otwinowski 19/11/2009
// last change: 2013-02-05 by M.Knichel
//

#include <TObject.h>
#include <THnSparse.h>
#include "AliTriggerAnalysis.h"
#include "AliPWG0Helper.h"

class TF1;
class TH1;
class TH2;
class TH3;
class TTree;
class THnSparse;
class TParticle;
class AliHeader;
class AliGenEventHeader;
class AliStack;
class AliESD;
class AliESDEvent;
class AliESDtrack;
class AliMCEvent;
class AliESDVertex;
class AliESDtrackCuts;
class AlidNdPtAcceptanceCuts;
class AlidNdPtEventCuts;
class AliPWG0Helper;
class AliOfflineTrigger;
class AliMultiplicity;

class AliGenDPMjetEventHeader;
class AliGenCocktailEventHeader;
class AliGenPythiaEventHeader;
class AliVertexerTracks;
class AliLog;
class AliHeader;

class AlidNdPtHelper : public TObject
{
  public:
    enum AnalysisMode { kInvalid = -1, kSPD = 0, kTPC, kTPCITS, kTPCSPDvtx, kTPCSPDvtxUpdate, kTPCTrackSPDvtx, kTPCTrackSPDvtxUpdate, kTPCITSHybrid, kTPCITSHybridTrackSPDvtx, kTPCITSHybridTrackSPDvtxDCArPt, kITSStandAloneTrackSPDvtx,kITSStandAloneTPCTrackSPDvtx, kMCRec };

    enum ParticleMode { kAllPart = 0, kMCPion, kMCKaon, kMCProton, kPlus, kMinus, kCosmic, kBackgroundTrack, kMCRest, kVZEROCase1, kVZEROCase2};

    enum OutputObject { kInvalidObject = -1, kCutAnalysis = 0, kAnalysis, kAnalysisPbPb, kCorrection, kSystematics, kCutAnalysisPbPb };

    enum TrackObject  { kInvalidTrackObject = -1, kAllTracks = 0, kAccTracks, kRecTracks, kMCTracks };
    enum EventObject  { kInvalidEventObject = -1, kAllEvents = 0, kTriggeredEvents, kAccEvents, kRecEvents, kMCEvents };
    enum CutSteps     { kCutSteps = 3 };

    static const AliESDVertex* GetVertex(AliESDEvent* const aEsd, const AlidNdPtEventCuts *const evtCuts,const  AlidNdPtAcceptanceCuts *const accCuts,const  AliESDtrackCuts *const trackCuts,  AnalysisMode analysisMethod, Bool_t debug = kFALSE,Bool_t bRedoTPC = kFALSE, Bool_t bUseMeanVertex = kFALSE);

    static const AliESDVertex* GetTPCVertexZ(const AliESDEvent* const aEsd, const AlidNdPtEventCuts *const evtCuts, const AlidNdPtAcceptanceCuts *const accCuts, const AliESDtrackCuts *const trackCuts, Float_t fraction=0.8, Int_t ntracksMin=2);

    static Bool_t TestRecVertex(const AliESDVertex* vertex, const AliESDVertex* vertexSPD, AnalysisMode analysisMode, Bool_t debug = kFALSE);

    static Bool_t IsPrimaryParticle(AliStack *const stack, Int_t idx, ParticleMode particleMode);
    static Bool_t IsCosmicTrack(AliESDtrack *const track1, AliESDtrack *const track2);
    static Bool_t IsGoodImpPar(const AliESDtrack *const track);
    static Int_t ConvertPdgToPid(const TParticle *const particle);

    static Bool_t SelectEvent(const AliESDEvent* const aEsd, AliESDtrackCuts* const esdTrackCuts);
    static Bool_t SelectMCEvent(AliMCEvent* const mcEvent);

    static TObjArray *GetAllChargedTracks(AliESDEvent *const esdEvent, AnalysisMode analysisMode);

    static TH1F* MakeResol(TH2F * const his, Int_t integ, Bool_t type, Bool_t drawBins, Int_t minHistEntries);
    static TH1F* CreateResHisto(TH2F* const hRes2, TH1F **phMean, Int_t integ,  Bool_t drawBinFits, Int_t minHistEntries);

    static Int_t GetTPCMBTrackMult(const AliESDEvent* const esdEvent, const AlidNdPtEventCuts *const evtCuts,const  AlidNdPtAcceptanceCuts *const accCuts, const AliESDtrackCuts *const trackCuts);
    static Int_t GetTPCMBPrimTrackMult(const AliESDEvent* const esdEvent, AliStack * const stack,const AlidNdPtEventCuts *const evtCuts, const AlidNdPtAcceptanceCuts *const accCuts, const AliESDtrackCuts *const trackCuts);

    static Int_t GetSPDMBTrackMult(const AliESDEvent* const esdEvent, Float_t deltaThetaCut =0.025, Float_t deltaPhiCut = 0.08);
    static Int_t GetSPDMBPrimTrackMult(const AliESDEvent* const esdEvent, AliStack *const  stack, Float_t deltaThetaCut =0.025, Float_t deltaPhiCut = 0.08);
    static Int_t GetMCTrueTrackMult(AliMCEvent *const mcEvent, AlidNdPtEventCuts *const evtCuts, AlidNdPtAcceptanceCuts *const accCuts);

    static AliESDtrack* GetTPCOnlyTrackSPDvtx(const AliESDEvent* const esdEvent, Int_t iTrack, Bool_t bUpdate);
    static AliESDtrack* GetTPCOnlyTrackTrackSPDvtx(const AliESDEvent* const esdEvent, Int_t iTrack, Bool_t bUpdate);
    static AliESDtrack* GetTrackSPDvtx(const AliESDEvent* const esdEvent, Int_t iTrack, Bool_t bUpdate);
    static AliESDtrack* GetTrackTrackSPDvtx(const AliESDEvent* const esdEvent, Int_t iTrack, Bool_t bUpdate);

    static void PrintConf(AnalysisMode analysisMode, AliTriggerAnalysis::Trigger trigger);
    static void PrintMCInfo(AliStack *const pStack, Int_t label);

    static TH1* ScaleByBinWidth(TH1 *const hist=0);
    static TH1* GetContCorrHisto(TH1 *const hist=0);
    static TH1* CalcRelativeDifference(const TH1 *const hist1=0, const TH1 *const hist2=0);
    static TH1* CalcRelativeDifferenceFun(const TH1 *const hist=0, TF1 *const fun=0);
    static TH1* NormalizeToEvent(const TH2 *const hist1=0, const TH1 *const hist2=0);

    //static THnSparse* GenerateCorrMatrix(THnSparse *const hist1, const THnSparse *const hist2, char *const name);
    //static TH2* GenerateCorrMatrix(TH2 *const hist1, TH2 *const hist2, char *const name);
    //static TH1* GenerateCorrMatrix(TH1 *const hist1, TH1 *const hist2, char *const name);

    //static THnSparse* GenerateContCorrMatrix(THnSparse *const hist1, THnSparse *const hist2, char *const name);
    //static TH2* GenerateContCorrMatrix(TH2 *const hist1, TH2 *const hist2, char *const name);
    //static TH1* GenerateContCorrMatrix(TH1 *const hist1, TH1 *const hist2, char *const name);

    static THnSparse* GenerateCorrMatrix(THnSparse *const hist1, const THnSparse *const hist2, const char * name);
    static TH2* GenerateCorrMatrix(TH2 *const hist1, TH2 *const hist2, const char* name);
    static TH1* GenerateCorrMatrix(TH1 *const hist1, TH1 *const hist2, const char* name);

    static THnSparse* GenerateContCorrMatrix(THnSparse *const hist1, const THnSparse *const hist2, const char* name);
    static TH2* GenerateContCorrMatrix(TH2 *const hist1, TH2 *const hist2, const char* name);
    static TH1* GenerateContCorrMatrix(TH1 *const hist1, TH1 *const hist2, const char* name);

    static Double_t GetStrangenessCorrFactor(const Double_t pt);
    static Double_t GetStrangenessCorrFactorPbPb(const Double_t pt);    
    static Double_t GetLinearInterpolationValue(const Double_t x1, const Double_t y1, const Double_t x2, const Double_t y2, const Double_t pt);

    // function to rebin THnSparse, the content of hist1 will be rebinned, hist2 serves as a protoype for the binning
    static THnSparse* RebinTHnSparse(const THnSparse* hist1, THnSparse* hist2, const Char_t* newname = "",  Option_t* option = "");
    
    // function to get processtype (kSD, kND) from DPMJET header in PA
    static AliPWG0Helper::MCProcessType GetEventProcessTypePA(AliHeader* aHeader, Bool_t adebug = kFALSE);
    static AliPWG0Helper::MCProcessType GetDPMjetEventProcessTypePA(AliGenEventHeader* aHeader, Bool_t adebug = kFALSE);

    ClassDef(AlidNdPtHelper, 2);

  private:
    AlidNdPtHelper(const AlidNdPtHelper&);
    AlidNdPtHelper& operator=(const AlidNdPtHelper&);
};

#endif

