/* $Id: AlidNdPtHelper.h 28655 2008-09-10 12:57:42Z jgrosseo $ */

#ifndef ALIDNDPTHELPER_H
#define ALIDNDPTHELPER_H

//
// static helper functions
// origin PWG0 (Jan Fiete, CKB) and extended by Jacek Otwinowski (JO)
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

#include "THnSparse.h"
class AlidNdPtHelper : public TObject
{
  public:
    enum Trigger { kMB1 = 0, kMB2, kSPDFASTOR }; // definition from ALICE-INT-2005-025
    enum AnalysisMode { kInvalid = -1, kSPD = 0, kTPC, kTPCITS, kTPCSPDvtx, kMCRec, kMCPion, kMCKaon, kMCProton, kPlus, kMinus };
    // in case we want to use bitmaps...
    // kDiffractiveProcess is artifficial
    enum MCProcessType { kInvalidProcess = -1, kND = 0x1, kDD = 0x2, kSD = 0x4, kDiffractiveProcess = 0x9 }; 

    static Bool_t IsEventTriggered(const AliESD* aEsd, Trigger trigger = kMB2);
    static Bool_t IsEventTriggered(ULong64_t triggerMask, Trigger trigger = kMB2);
    static const AliESDVertex* GetVertex(AliESDEvent* aEsd, AlidNdPtEventCuts *evtCuts, AlidNdPtAcceptanceCuts *accCuts, AliESDtrackCuts *trackCuts,  AnalysisMode analysisMethod, Bool_t debug = kFALSE,Bool_t bRedoTPC = kFALSE, Bool_t bUseMeanVertex = kFALSE);
    //static const AliESDVertex* GetVertex(AliESDEvent* aEsd, AnalysisMode analysisMethod, Bool_t debug = kFALSE,Bool_t bRedoTPC = kFALSE, Bool_t bUseMeanVertex = kFALSE);

    static const AliESDVertex* GetTPCVertexZ(AliESDEvent* aEsd, Float_t sigmaXYcut=3., Float_t distXYcut=3., Float_t distZcut=30., Int_t nclCut=50, Float_t fraction=0.8, Int_t ntracksMin=2);

    static const Bool_t TestVertex(const AliESDVertex* vertex, AnalysisMode analysisMode, Bool_t debug = kFALSE);

    static Bool_t IsPrimaryCharged(TParticle* aParticle, Int_t aTotalPrimaries, Bool_t adebug = kFALSE);
    static Bool_t IsPrimaryParticle(AliStack *stack, Int_t idx, AnalysisMode analysisMode);

    static AlidNdPtHelper::MCProcessType GetEventProcessType(AliHeader* aHeader, Bool_t adebug = kFALSE);
    static AlidNdPtHelper::MCProcessType GetPythiaEventProcessType(AliGenEventHeader* aHeader, Bool_t adebug = kFALSE);
    static AlidNdPtHelper::MCProcessType GetDPMjetEventProcessType(AliGenEventHeader* aHeader, Bool_t adebug = kFALSE);
    static Int_t GetLastProcessType() { return fgLastProcessType; }

    static TParticle* FindPrimaryMother(AliStack* stack, Int_t label);
    static Int_t FindPrimaryMotherLabel(AliStack* stack, Int_t label);

    static void CreateProjections(TH3* hist, Bool_t save = kFALSE);
    static void CreateDividedProjections(TH3* hist, TH3* hist2, const char* axis = 0, Bool_t putErrors = kFALSE, Bool_t save = kFALSE);
    static const char* GetAxisTitle(TH3* hist, const char axis);

    static void NormalizeToBinWidth(TH1* hist);
    static void NormalizeToBinWidth(TH2* hist);

    static void PrintConf(AnalysisMode analysisMode, Trigger trigger);

    // added by JO
    static Int_t ConvertPdgToPid(TParticle *particle);

    enum OutputObject { kInvalidObject = -1, kCutAnalysis = 0, kAnalysis, kCorrection, kSystematics };
    enum TrackObject  { kInvalidTrackObject = -1, kAllTracks = 0, kAccTracks, kRecTracks, kMCTracks };
    enum EventObject  { kInvalidEventObject = -1, kAllEvents = 0, kTriggeredEvents, kAccEvents, kRecEvents, kMCEvents };
    enum CutSteps     { kCutSteps = 3 };

    //static TObjArray *GetAllChargedTracks(AliESDEvent *esdEvent, AnalysisMode analysisMode);
    static TObjArray *GetAllChargedTracks(AliESDEvent *esdEvent, const AliESDVertex *vtx, AnalysisMode analysisMode);
    static AliESDtrack* GetTPCOnlyTrack(AliESDEvent* esd, const AliESDVertex *vtx, Int_t iTrack);

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

  protected:
    static Int_t fgLastProcessType;    // stores the raw value of the last process type extracnted
 
    ClassDef(AlidNdPtHelper, 0);

  private:
    AlidNdPtHelper(const AlidNdPtHelper&);
    AlidNdPtHelper& operator=(const AlidNdPtHelper&);
};

#endif

