#ifndef AliPHOSJetJetMC_cxx
#define AliPHOSJetJetMC_cxx

//Author: Daiki Sekihata (Hiroshima University)

#include "TObject.h"
#include "TString.h"
#include "AliVEvent.h"
#include "AliGenPythiaEventHeader.h"

class AliVEvent;
class AliGenPythiaEventHeader;

class AliPHOSJetJetMC : public TObject {

  public:
    AliPHOSJetJetMC();
    AliPHOSJetJetMC(Int_t pThardbin);
    virtual ~AliPHOSJetJetMC();
 
    Int_t GetPtHardBin()  {return fPtHardBin;}
    //Float_t GetPtHard()   {return fPtHard;}
    //Float_t GetXsection() {return fXsection;}
    //Int_t GetNTrials()    {return fNTrials;}

    AliGenPythiaEventHeader* GetPythiaEventHeader(AliVEvent *event);
    Bool_t ComparePtHardWithJet(AliVEvent *event);
    Bool_t ComparePtHardWithSingleParticle(AliVEvent *event);
    Bool_t ComparePtHardBin(AliVEvent *event);

    Int_t GetFirstJetIndex() {return fFirstJetIndex;}
    Int_t GetLastJetIndex()  {return fLastJetIndex;}
    Int_t GetGeneratorJetIndex() {return fGenJetID;}
    Int_t GetFirstUEIndex() {return fFirstJetIndex;}
    Int_t GetLastUEIndex()  {return fLastJetIndex;}
    Int_t GetGeneratorUEIndex() {return fGenUEID;}
    void ConfigureJetJetMC(AliVEvent *event);

    void SetJetPtFactor(Double_t factor) {fPtHardAndJetPtFactor = factor;}
    void SetSingleParticlePtFactor(Double_t factor) {fPtHardAndSinglePtFactor = factor;}

  protected:
    TList *GetGenHeaderList(AliVEvent *event);
    TString GetProductionName();    

  
  private:
    Int_t fPtHardBin;
    Float_t fPtHard;
    Float_t fXsection;
    Int_t fNTrials;
    Double_t fPtHardAndJetPtFactor;
    Double_t fPtHardAndSinglePtFactor;
    Int_t fFirstJetIndex;
    Int_t fLastJetIndex;
    Int_t fGenJetID;
    Int_t fFirstUEIndex;
    Int_t fLastUEIndex;
    Int_t fGenUEID;//generator index of underlying event. (HIJING or PYTHIA)

    AliPHOSJetJetMC(const AliPHOSJetJetMC&);
    AliPHOSJetJetMC& operator=(const AliPHOSJetJetMC&);

    ClassDef(AliPHOSJetJetMC, 8);

};

#endif

