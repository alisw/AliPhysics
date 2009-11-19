#ifndef ALIDNDPTHELPER_H
#define ALIDNDPTHELPER_H

//
// static dNdPt helper functions
//
// Origin: Jan Fiete Grosse-Oetringhaus
// Modified and Extended: Jacek Otwinowski 19/11/2009
//

#include <TObject.h>

class TParticle;
class TH1;
class TH1F;
class TH2F;
class TH2;
class TH3;
class TF1;
class TTree;
class TArrayI;

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

#include "AliPWG0Helper.h"
#include "THnSparse.h"

class AlidNdPtHelper : public TObject
{
  public:
    enum AnalysisMode { kInvalid = -1, kSPD = 0, kTPC, kTPCITS, kTPCSPDvtx, kMCRec, kMCPion, kMCKaon, kMCProton, kPlus, kMinus };
    static const AliESDVertex* GetVertex(AliESDEvent* aEsd, AlidNdPtEventCuts *evtCuts, AlidNdPtAcceptanceCuts *accCuts, AliESDtrackCuts *trackCuts,  AnalysisMode analysisMethod, Bool_t debug = kFALSE,Bool_t bRedoTPC = kFALSE, Bool_t bUseMeanVertex = kFALSE);

    static const AliESDVertex* GetTPCVertexZ(AliESDEvent* aEsd, AlidNdPtEventCuts *evtCuts, AlidNdPtAcceptanceCuts *accCuts, AliESDtrackCuts *trackCuts, Float_t fraction=0.8, Int_t ntracksMin=2);

    static Bool_t TestRecVertex(const AliESDVertex* vertex, AnalysisMode analysisMode, Bool_t debug = kFALSE);

    static Bool_t IsPrimaryParticle(AliStack *stack, Int_t idx, AnalysisMode analysisMode);
    static void PrintConf(AnalysisMode analysisMode, AliPWG0Helper::Trigger trigger);
    static Int_t ConvertPdgToPid(TParticle *particle);

    enum OutputObject { kInvalidObject = -1, kCutAnalysis = 0, kAnalysis, kCorrection, kSystematics };
    enum TrackObject  { kInvalidTrackObject = -1, kAllTracks = 0, kAccTracks, kRecTracks, kMCTracks };
    enum EventObject  { kInvalidEventObject = -1, kAllEvents = 0, kTriggeredEvents, kAccEvents, kRecEvents, kMCEvents };
    enum CutSteps     { kCutSteps = 3 };

    static TObjArray *GetAllChargedTracks(AliESDEvent *esdEvent, AnalysisMode analysisMode);

    static TH1F* MakeResol(TH2F * his, Int_t integ, Bool_t type, Bool_t drawBins, Int_t minHistEntries);
    static TH1F* CreateResHisto(TH2F* hRes2, TH1F **phMean, Int_t integ,  Bool_t drawBinFits, Int_t minHistEntries);

    static Int_t GetTPCMBTrackMult(AliESDEvent* esdEvent, AlidNdPtEventCuts *evtCuts, AlidNdPtAcceptanceCuts *accCuts, AliESDtrackCuts *trackCuts);
    static Int_t GetTPCMBPrimTrackMult(AliESDEvent* esdEvent, AliStack * stack,AlidNdPtEventCuts *evtCuts, AlidNdPtAcceptanceCuts *accCuts, AliESDtrackCuts *trackCuts);

    static Int_t GetSPDMBTrackMult(AliESDEvent* esdEvent, Float_t deltaThetaCut =0.025, Float_t deltaPhiCut = 0.08);
    static Int_t GetSPDMBPrimTrackMult(AliESDEvent* esdEvent, AliStack * stack, Float_t deltaThetaCut =0.025, Float_t deltaPhiCut = 0.08);
    static Int_t GetMCTrueTrackMult(AliMCEvent *mcEvent, AlidNdPtEventCuts *evtCuts, AlidNdPtAcceptanceCuts *accCuts);

    static void PrintMCInfo(AliStack *pStack, Int_t label);

    static TH1* ScaleByBinWidth(TH1 *hist=0);
    static TH1* GetContCorrHisto(TH1 *hist=0);
    static TH1* CalcRelativeDifference(TH1 *hist1=0, TH1 *hist2=0);
    static TH1* CalcRelativeDifferenceFun(TH1 *hist=0, TF1 *fun=0);
    static TH1* NormalizeToEvent(TH2 *hist1=0, TH1 *hist2=0);

    static THnSparse* GenerateCorrMatrix(THnSparse *hist1, THnSparse *hist2, char *name);
    static TH2* GenerateCorrMatrix(TH2 *hist1, TH2 *hist2, char *name);
    static TH1* GenerateCorrMatrix(TH1 *hist1, TH1 *hist2, char *name);

    static THnSparse* GenerateContCorrMatrix(THnSparse *hist1, THnSparse *hist2, char *name);
    static TH2* GenerateContCorrMatrix(TH2 *hist1, TH2 *hist2, char *name);
    static TH1* GenerateContCorrMatrix(TH1 *hist1, TH1 *hist2, char *name);

    ClassDef(AlidNdPtHelper, 0);

  private:
    AlidNdPtHelper(const AlidNdPtHelper&);
    AlidNdPtHelper& operator=(const AlidNdPtHelper&);
};

#endif

