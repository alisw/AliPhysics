#ifndef AliAnalysisTaskResolution_cxx
#define AliAnalysisTaskResolution_cxx

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


class AliAnalysisTaskResolution : public AliAnalysisTaskSE{

 public:
   AliAnalysisTaskResolution();
   AliAnalysisTaskResolution(const char *name);
   virtual ~AliAnalysisTaskResolution();

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
   Int_t CountTracks0914();
   Int_t CountTracks09();

   AliV0ReaderV1 *fV0Reader;
   TClonesArray *fConversionGammas; //Reconstructed Photons;
   AliConversionCuts *fConversionCuts; // Cuts used by the V0Reader
   TTree *fTreeEvent;
   TTree *fTreeResolution;
   Float_t fPrimVtxZ;
   Int_t fNContrVtx;
   Int_t fNESDtracksEta09;
   Int_t fNESDtracksEta0914;
   Int_t fNESDtracksEta14;
   TVectorF fGammaRecCoords;
   TVectorF fGammaMCCoords;
   Float_t fChi2ndf;
   Bool_t fIsHeavyIon;
   TList *fOutputList;
   TList *fEventList;
   TList *fResolutionList;
   AliESDEvent *fESDEvent;
   AliMCEvent *fMCEvent;

   AliAnalysisTaskResolution(const AliAnalysisTaskResolution&); // not implemented
   AliAnalysisTaskResolution& operator=(const AliAnalysisTaskResolution&); // not implemented


   ClassDef(AliAnalysisTaskResolution, 0);
};

#endif
