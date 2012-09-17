#ifndef AliAnalysisConversionQA_cxx
#define AliAnalysisConversionQA_cxx

#include "AliAnalysisTaskSE.h"
#include "AliConversionPhotonBase.h"
#include "TH1.h"
#include "TH2.h"
#include "TTreeStream.h"
#include "AliLog.h"
#include <vector>
#include "AliV0ReaderV1.h"
#include "AliConversionCuts.h"
#include "TList.h"
#include "AliStack.h"
#include "TClonesArray.h"


using namespace std;


class AliAnalysisTaskConversionQA : public AliAnalysisTaskSE{

public:

    AliAnalysisTaskConversionQA(const char *name);
    virtual ~AliAnalysisTaskConversionQA();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);

    void SetV0Reader(AliV0ReaderV1 *v0Reader){fV0Reader=v0Reader;}
    void SetConversionCuts(AliConversionCuts* conversionCuts,Bool_t IsHeavyIon ){
       fConversionCuts=conversionCuts;
       fIsHeavyIon = IsHeavyIon;
    }
    
private:

    void ProcessQA();

    AliV0ReaderV1 *fV0Reader;
    TClonesArray *fConversionGammas; //Reconstructed Photons;
    AliConversionCuts *fConversionCuts; // Cuts used by the V0Reader
    TTreeSRedirector *fStreamQA;
    Bool_t fIsHeavyIon;
    TList *fOutputList;

    ClassDef(AliAnalysisTaskConversionQA, 0);
};

#endif

