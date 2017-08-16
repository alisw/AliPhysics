#ifndef AliPHOSJetJetMC_cxx
#define AliPHOSJetJetMC_cxx

//Author: Daiki Sekihata (Hiroshima University)

#include "TObject.h"
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
    void ConfigureJetJetMC(AliVEvent *event);

  protected:
    TList *GetGenHeaderList(AliVEvent *event);
      
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

    AliPHOSJetJetMC(const AliPHOSJetJetMC&);
    AliPHOSJetJetMC& operator=(const AliPHOSJetJetMC&);

    ClassDef(AliPHOSJetJetMC, 5);

};

#endif

