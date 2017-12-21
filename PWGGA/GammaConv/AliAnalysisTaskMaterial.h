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
#include "AliConversionPhotonCuts.h"
#include "AliConvEventCuts.h"
#include "TList.h"
#include "TClonesArray.h"


using namespace std;


class AliAnalysisTaskMaterial : public AliAnalysisTaskSE{

	public:

		AliAnalysisTaskMaterial();
		AliAnalysisTaskMaterial(const char *name);
		virtual ~AliAnalysisTaskMaterial();

		virtual void   UserCreateOutputObjects();
    virtual Bool_t Notify                   ();
		virtual void   UserExec(Option_t *option);
		virtual void   Terminate(Option_t *);

		void SetIsMC(Bool_t isMC){fIsMC=isMC;}
		void SetV0Reader(AliV0ReaderV1 *v0Reader){fV0Reader=v0Reader;}
        void SetV0ReaderName(TString name){fV0ReaderName=name; return;}
		void SetConversionCuts(AliConversionPhotonCuts* conversionCuts,Int_t IsHeavyIon ){
			fConversionCuts=conversionCuts;
			fIsHeavyIon = IsHeavyIon;
		}
		void SetEventCuts(AliConvEventCuts* conversionCuts,Int_t IsHeavyIon ){
			fEventCuts=conversionCuts;
			fIsHeavyIon = IsHeavyIon;
		}
		
	private:
		
		void ProcessPhotons();
		void ProcessMCPhotons();
        void FillMCTree(Int_t eventPos);
		Int_t CountTracks0914();
		Int_t CountTracks09();

		AliV0ReaderV1 				*fV0Reader;					// 
        TString                     fV0ReaderName;
		TClonesArray 				*fConversionGammas; 		// Reconstructed Photons;
		AliConversionPhotonCuts 	*fConversionCuts; 			// Cuts used by the V0Reader
		AliConvEventCuts 			*fEventCuts; 				// Cuts used by the V0Reader
		TList 						*fOutputList; 				//
		TList 						*fEventList;				//
		TList 						*fRecGammaList;				//
		TList 						*fAllMCGammaList;			//
		TList 						*fAllMCConvGammaList;		//
		TTree* 						fTreeEvent;					//
		TTree* 						fTreeMaterialRec;			//
		TTree* 						fTreeMaterialAllGamma;		//
		TTree* 						fTreeMaterialConvGamma;		//
		Float_t 					fPrimVtxZ;					//
		Int_t 						fNContrVtx;					//
		Int_t 						fNESDtracksEta09;			//
		Int_t 						fNESDtracksEta0914;			//
		Int_t 						fNESDtracksEta14;			//
		Float_t 					fGammaMCPt;					//
		Float_t 					fGammaMCTheta;				//
		Float_t 					fGammaMCConvPt;				//
		Float_t 					fGammaMCConvTheta;			//
		TVectorF 					fMCConvCords;				//
		TVectorF 					fMCConvDaughterProp;		//
		Float_t 					fGammaPt;					//
		Float_t 					fGammaTheta;				//
		Float_t 					fGammaChi2NDF;				//
		TVectorF 					fRecCords;					//
		TVectorF 					fDaughterProp;				//
		UChar_t 					fKind;						//
		Int_t 						fIsHeavyIon;				//
		Bool_t 						fIsMC;						//
		AliESDEvent 				*fESDEvent;					//
		AliMCEvent 					*fMCEvent;					//

		AliAnalysisTaskMaterial(const AliAnalysisTaskMaterial&); // not implemented
		AliAnalysisTaskMaterial& operator=(const AliAnalysisTaskMaterial&); // not implemented


        ClassDef(AliAnalysisTaskMaterial, 3);
};

#endif
