/* $Id$ */

#ifndef ALIPWG0HELPER_H
#define ALIPWG0HELPER_H

#include <TObject.h>

// static helper functions

class AliESD;
class AliESDEvent;
class AliESDVertex;
class TParticle;
class TH1;
class TH2;
class TH3;
class AliHeader;
class AliStack;
class TTree;

class AliPWG0Helper : public TObject
{
  public:
    enum Trigger { kMB1 = 0, kMB2 }; // definition from ALICE-INT-2005-025
    enum AnalysisMode { kSPD = 0, kTPC, kTPCITS };

    static Bool_t IsEventTriggered(const AliESD* aEsd, Trigger trigger = kMB2);
    static Bool_t IsEventTriggered(ULong64_t triggerMask, Trigger trigger = kMB2);
    static Bool_t IsVertexReconstructed(const AliESD* aEsd);
    static Bool_t IsVertexReconstructed(const AliESDVertex* vtxESD);
    static const AliESDVertex* GetVertex(const AliESDEvent* aEsd, AnalysisMode analysisMethod, Bool_t debug = kFALSE);

    static Bool_t IsPrimaryCharged(TParticle* aParticle, Int_t aTotalPrimaries, Bool_t adebug = kFALSE);

    static Int_t GetPythiaEventProcessType(AliHeader* aHeader, Bool_t adebug = kFALSE);
    static TParticle* FindPrimaryMother(AliStack* stack, Int_t label);
    static Int_t FindPrimaryMotherLabel(AliStack* stack, Int_t label);

    static void CreateProjections(TH3* hist, Bool_t save = kFALSE);
    static void CreateDividedProjections(TH3* hist, TH3* hist2, const char* axis = 0, Bool_t putErrors = kFALSE, Bool_t save = kFALSE);
    static const char* GetAxisTitle(TH3* hist, const char axis);

    static void SetBranchStatusRecursive(TTree* tree, char *bname, Bool_t status, Bool_t debug = kFALSE);

    static void NormalizeToBinWidth(TH1* hist);
    static void NormalizeToBinWidth(TH2* hist);

  protected:
    ClassDef(AliPWG0Helper, 0)

  private:
    AliPWG0Helper(const AliPWG0Helper&);
    AliPWG0Helper& operator=(const AliPWG0Helper&);
};

#endif

