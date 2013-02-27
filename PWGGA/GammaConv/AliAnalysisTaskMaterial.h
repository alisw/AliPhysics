#ifndef AliAnalysisTaskMaterial_cxx
#define AliAnalysisTaskMaterial_cxx

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


class AliAnalysisTaskMaterial : public AliAnalysisTaskSE{

public:

    AliAnalysisTaskMaterial(const char *name);
    virtual ~AliAnalysisTaskMaterial();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);

    void SetV0Reader(AliV0ReaderV1 *v0Reader){fV0Reader=v0Reader;}
    void SetConversionCuts(AliConversionCuts* conversionCuts,Bool_t IsHeavyIon ){
       fConversionCuts=conversionCuts;
       fIsHeavyIon = IsHeavyIon;
    }
    
private:

    void ProcessPhotons();
	 void ProcessMCPhotons();
	 void FillMCTree(Int_t stackPos);
	 Int_t CountESDTracks14();
	 Int_t CountESDTracks0914();
	 Int_t CountESDTracks09();
	 
    AliV0ReaderV1 *fV0Reader;
    TClonesArray *fConversionGammas; //Reconstructed Photons;
    AliConversionCuts *fConversionCuts; // Cuts used by the V0Reader
    TTreeSRedirector *fStreamMaterial;
	 TTreeSRedirector *fStreamResolution;
    Bool_t fIsHeavyIon;
    TList *fOutputList;
	 AliESDEvent *fESDEvent;
	 AliMCEvent *fMCEvent;
	
    ClassDef(AliAnalysisTaskMaterial, 0);
};

#endif

