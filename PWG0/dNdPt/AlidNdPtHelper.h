#ifndef ALIDNDPTHELPER_H
#define ALIDNDPTHELPER_H

//
// static dNdPt helper functions
//
// Origin: Jan Fiete Grosse-Oetringhaus
// Modified and Extended: Jacek Otwinowski 19/11/2009
//

#include <TObject.h>
#include <THnSparse.h>
#include "AliTriggerAnalysis.h"

class TF1;
class TH1;
class TH2;
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

    enum OutputObject { kInvalidObject = -1, kCutAnalysis = 0, kAnalysis, kCorrection, kSystematics };
    enum TrackObject  { kInvalidTrackObject = -1, kAllTracks = 0, kAccTracks, kRecTracks, kMCTracks };
    enum EventObject  { kInvalidEventObject = -1, kAllEvents = 0, kTriggeredEvents, kAccEvents, kRecEvents, kMCEvents };
    enum CutSteps     { kCutSteps = 3 };

    static const AliESDVertex* GetVertex(AliESDEvent* const aEsd, AlidNdPtEventCuts *const evtCuts, AlidNdPtAcceptanceCuts *const accCuts, AliESDtrackCuts *const trackCuts,  AnalysisMode analysisMethod, Bool_t debug = kFALSE,Bool_t bRedoTPC = kFALSE, Bool_t bUseMeanVertex = kFALSE);

    static const AliESDVertex* GetTPCVertexZ(AliESDEvent* const aEsd, AlidNdPtEventCuts *const evtCuts, AlidNdPtAcceptanceCuts *const accCuts, AliESDtrackCuts *const trackCuts, Float_t fraction=0.8, Int_t ntracksMin=2);

    static Bool_t TestRecVertex(const AliESDVertex* vertex, const AliESDVertex* vertexSPD, AnalysisMode analysisMode, Bool_t debug = kFALSE);

    static Bool_t IsPrimaryParticle(AliStack *const stack, Int_t idx, ParticleMode particleMode);
    static Bool_t IsCosmicTrack(AliESDtrack *const track1, AliESDtrack *const track2);
    static Bool_t IsGoodImpPar(AliESDtrack *const track);
    static Int_t ConvertPdgToPid(TParticle *const particle);

    static TObjArray *GetAllChargedTracks(AliESDEvent *const esdEvent, AnalysisMode analysisMode);

    static TH1F* MakeResol(TH2F * const his, Int_t integ, Bool_t type, Bool_t drawBins, Int_t minHistEntries);
    static TH1F* CreateResHisto(TH2F* const hRes2, TH1F **phMean, Int_t integ,  Bool_t drawBinFits, Int_t minHistEntries);

    static Int_t GetTPCMBTrackMult(AliESDEvent* const esdEvent, AlidNdPtEventCuts *const evtCuts, AlidNdPtAcceptanceCuts *const accCuts, AliESDtrackCuts *const trackCuts);
    static Int_t GetTPCMBPrimTrackMult(AliESDEvent* const esdEvent, AliStack * const stack,AlidNdPtEventCuts *const evtCuts, AlidNdPtAcceptanceCuts *const accCuts, AliESDtrackCuts *const trackCuts);

    static Int_t GetSPDMBTrackMult(AliESDEvent* const esdEvent, Float_t deltaThetaCut =0.025, Float_t deltaPhiCut = 0.08);
    static Int_t GetSPDMBPrimTrackMult(AliESDEvent* const esdEvent, AliStack *const  stack, Float_t deltaThetaCut =0.025, Float_t deltaPhiCut = 0.08);
    static Int_t GetMCTrueTrackMult(AliMCEvent *const mcEvent, AlidNdPtEventCuts *const evtCuts, AlidNdPtAcceptanceCuts *const accCuts);

    static AliESDtrack* GetTPCOnlyTrackSPDvtx(AliESDEvent* const esdEvent, Int_t iTrack, Bool_t bUpdate);
    static AliESDtrack* GetTPCOnlyTrackTrackSPDvtx(AliESDEvent* const esdEvent, Int_t iTrack, Bool_t bUpdate);
    static AliESDtrack* GetTrackSPDvtx(AliESDEvent* const esdEvent, Int_t iTrack, Bool_t bUpdate);
    static AliESDtrack* GetTrackTrackSPDvtx(AliESDEvent* const esdEvent, Int_t iTrack, Bool_t bUpdate);

    static void PrintConf(AnalysisMode analysisMode, AliTriggerAnalysis::Trigger trigger);
    static void PrintMCInfo(AliStack *const pStack, Int_t label);

    static TH1* ScaleByBinWidth(TH1 *const hist=0);
    static TH1* GetContCorrHisto(TH1 *const hist=0);
    static TH1* CalcRelativeDifference(TH1 *const hist1=0, TH1 *const hist2=0);
    static TH1* CalcRelativeDifferenceFun(TH1 *const hist=0, TF1 *const fun=0);
    static TH1* NormalizeToEvent(TH2 *const hist1=0, TH1 *const hist2=0);

    static THnSparse* GenerateCorrMatrix(THnSparse *const hist1, THnSparse *const hist2, char *const name);
    static TH2* GenerateCorrMatrix(TH2 *const hist1, TH2 *const hist2, char *const name);
    static TH1* GenerateCorrMatrix(TH1 *const hist1, TH1 *const hist2, char *const name);

    static THnSparse* GenerateContCorrMatrix(THnSparse *const hist1, THnSparse *const hist2, char *const name);
    static TH2* GenerateContCorrMatrix(TH2 *const hist1, TH2 *const hist2, char *const name);
    static TH1* GenerateContCorrMatrix(TH1 *const hist1, TH1 *const hist2, char *const name);

    static Double_t GetStrangenessCorrFactor(const Double_t pt);
    static Double_t GetLinearInterpolationValue(const Double_t x1, const Double_t y1, const Double_t x2, const Double_t y2, const Double_t pt);

    ClassDef(AlidNdPtHelper, 0);

  private:
    AlidNdPtHelper(const AlidNdPtHelper&);
    AlidNdPtHelper& operator=(const AlidNdPtHelper&);
};

#endif

