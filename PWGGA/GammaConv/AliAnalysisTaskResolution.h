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
#include "AliConvEventCuts.h"
#include "AliConversionPhotonCuts.h"
#include "TList.h"
#include "TClonesArray.h"


using namespace std;


class AliAnalysisTaskResolution : public AliAnalysisTaskSE{

	public:
		AliAnalysisTaskResolution();
		AliAnalysisTaskResolution(const char *name);
		virtual ~AliAnalysisTaskResolution();

		virtual void   UserCreateOutputObjects();
    virtual Bool_t Notify                   ();
		virtual void   UserExec(Option_t *option);
		virtual void   Terminate(Option_t *);

    void SetV0ReaderName(TString name){fV0ReaderName=name; return;}
		void SetV0Reader(AliV0ReaderV1 *v0Reader){fV0Reader=v0Reader;}
		void SetConversionCuts(AliConversionPhotonCuts* conversionCuts,Int_t IsHeavyIon ){
			fConversionCuts=conversionCuts;
			fIsHeavyIon = IsHeavyIon;
		}
		void SetEventCuts(AliConvEventCuts* conversionCuts,Int_t IsHeavyIon ){
			fEventCuts=conversionCuts;
			fIsHeavyIon = IsHeavyIon;
		}
		void SetIsMC(Bool_t isMC){fIsMC=isMC;}
	
	private:

		void ProcessPhotons();
		Int_t CountTracks0914();
		Int_t CountTracks09();

		AliV0ReaderV1 *fV0Reader;
        TString fV0ReaderName;
		TClonesArray *fConversionGammas; //Reconstructed Photons;
		AliConvEventCuts *fEventCuts; // Cuts used by the V0Reader
		AliConversionPhotonCuts *fConversionCuts; // Cuts used by the V0Reader
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
		Int_t fIsHeavyIon;
		Bool_t fIsMC;
		TList *fOutputList;
		TList *fEventList;
		TList *fResolutionList;
		AliESDEvent *fESDEvent;
		AliMCEvent *fMCEvent;

		AliAnalysisTaskResolution(const AliAnalysisTaskResolution&); // not implemented
		AliAnalysisTaskResolution& operator=(const AliAnalysisTaskResolution&); // not implemented


    ClassDef(AliAnalysisTaskResolution, 3);
};

#endif
