/* $Id$ */

#ifndef ALIPWG0HELPER_H
#define ALIPWG0HELPER_H

#include <TObject.h>
#include <AliTriggerAnalysis.h>

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
class AliMultiplicity;

class AliPWG0Helper : public TObject
{
  public:
    // kTPCSPD: tracks + tracklets not using any cluster associated to tracks 
    enum AnalysisMode { kInvalid = -1, kSPD = 0x1, kTPC = 0x2, kTPCITS = 0x4, kFieldOn = 0x8, kSPDOnlyL0 = 0x10, kTPCSPD = 0x20}; 
    // in case we want to use bitmaps...
    enum MCProcessType { kInvalidProcess = -1, kND = 0x1, kDD = 0x2, kSD = 0x4, kOnePart = 0x8 };
    enum DiffTreatment { kMCFlags = 0, kUA5Cuts = 1, kE710Cuts, kALICEHadronLevel };

    static const AliESDVertex* GetVertex(const AliESDEvent* aEsd, AnalysisMode analysisMethod, Bool_t debug = kFALSE);
    static Bool_t TestVertex(const AliESDVertex* vertex, AnalysisMode analysisMode, Bool_t debug = kFALSE);
    
    static Bool_t IsPrimaryCharged(TParticle* aParticle, Int_t aTotalPrimaries, Bool_t adebug = kFALSE);

    static MCProcessType GetEventProcessType(AliESDEvent* esd, AliHeader* header, AliStack* stack, DiffTreatment diffTreatment);
    static MCProcessType GetEventProcessType(AliHeader* aHeader, Bool_t adebug = kFALSE);
    static MCProcessType GetPythiaEventProcessType(AliGenEventHeader* aHeader, Bool_t adebug = kFALSE);
    static MCProcessType GetDPMjetEventProcessType(AliGenEventHeader* aHeader, Bool_t adebug = kFALSE);
    static Int_t GetLastProcessType() { return fgLastProcessType; }

    static TParticle* FindPrimaryMother(AliStack* stack, Int_t label);
    static Int_t FindPrimaryMotherLabel(AliStack* stack, Int_t label);

    static void CreateProjections(TH3* hist, Bool_t save = kFALSE);
    static void CreateDividedProjections(TH3* hist, TH3* hist2, const char* axis = 0, Bool_t putErrors = kFALSE, Bool_t save = kFALSE);
    static const char* GetAxisTitle(TH3* hist, const char axis);

    static void NormalizeToBinWidth(TH1* hist);
    static void NormalizeToBinWidth(TH2* hist);

    static void PrintConf(AnalysisMode analysisMode, AliTriggerAnalysis::Trigger trigger, DiffTreatment diffTreatment);
    
    static Double_t Rapidity(Double_t pt, Double_t pz, Double_t m);
    static Bool_t IsHadronLevelSingleDiffractive(AliStack* stack, Float_t cms, Float_t xiMin, Float_t xiMax);
    
  protected:
    static Int_t fgLastProcessType;    // stores the raw value of the last process type extracted
 
    ClassDef(AliPWG0Helper, 0)

  private:
    AliPWG0Helper(const AliPWG0Helper&);
    AliPWG0Helper& operator=(const AliPWG0Helper&);
};

#endif

