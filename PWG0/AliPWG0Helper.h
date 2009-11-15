/* $Id$ */

#ifndef ALIPWG0HELPER_H
#define ALIPWG0HELPER_H

#include <TObject.h>

// static helper functions

class AliESDEvent;
class AliESDVertex;
class TParticle;
class TH1;
class TH2;
class TH3;
class AliHeader;
class AliGenEventHeader;
class AliStack;
class TTree;
class AliOfflineTrigger;

class AliPWG0Helper : public TObject
{
  public:
    enum Trigger { kAcceptAll = 1, kMB1 = 2, kMB2, kMB3, kSPDGFO, kV0A, kV0C, kZDC, kZDCA, kZDCC, kFMDA, kFMDC, kFPANY, kStartOfFlags = 0x0100, kOfflineFlag = 0x8000 }; // MB1, MB2, MB3 definition from ALICE-INT-2005-025
    enum AnalysisMode { kInvalid = -1, kSPD = 0x1, kTPC = 0x2, kTPCITS = 0x4, kFieldOn = 0x8 };
    // in case we want to use bitmaps...
    enum MCProcessType { kInvalidProcess = -1, kND = 0x1, kDD = 0x2, kSD = 0x4 };

    static Bool_t IsEventTriggered(const AliESDEvent* aEsd, Trigger trigger);
    static Bool_t IsEventTriggered(ULong64_t triggerMask, Trigger trigger);
    static const AliESDVertex* GetVertex(AliESDEvent* aEsd, AnalysisMode analysisMethod, Bool_t debug = kFALSE, Bool_t bRedoTPC = kFALSE);
    static Bool_t TestVertex(const AliESDVertex* vertex, AnalysisMode analysisMode, Bool_t debug = kFALSE);
    
    static Bool_t IsPrimaryCharged(TParticle* aParticle, Int_t aTotalPrimaries, Bool_t adebug = kFALSE);

    static AliPWG0Helper::MCProcessType GetEventProcessType(AliHeader* aHeader, Bool_t adebug = kFALSE);
    static AliPWG0Helper::MCProcessType GetPythiaEventProcessType(AliGenEventHeader* aHeader, Bool_t adebug = kFALSE);
    static AliPWG0Helper::MCProcessType GetDPMjetEventProcessType(AliGenEventHeader* aHeader, Bool_t adebug = kFALSE);
    static Int_t GetLastProcessType() { return fgLastProcessType; }

    static TParticle* FindPrimaryMother(AliStack* stack, Int_t label);
    static Int_t FindPrimaryMotherLabel(AliStack* stack, Int_t label);

    static void CreateProjections(TH3* hist, Bool_t save = kFALSE);
    static void CreateDividedProjections(TH3* hist, TH3* hist2, const char* axis = 0, Bool_t putErrors = kFALSE, Bool_t save = kFALSE);
    static const char* GetAxisTitle(TH3* hist, const char axis);

    static void NormalizeToBinWidth(TH1* hist);
    static void NormalizeToBinWidth(TH2* hist);

    static const char* GetTriggerName(Trigger trigger);
    static void PrintConf(AnalysisMode analysisMode, Trigger trigger);
    
    static AliOfflineTrigger* GetOfflineTrigger();

  protected:
    static Int_t fgLastProcessType;    // stores the raw value of the last process type extracnted
    static AliOfflineTrigger* fgOfflineTrigger;  // class that implemenents the offline trigger logic
 
    ClassDef(AliPWG0Helper, 0)

  private:
    AliPWG0Helper(const AliPWG0Helper&);
    AliPWG0Helper& operator=(const AliPWG0Helper&);
};

#endif

