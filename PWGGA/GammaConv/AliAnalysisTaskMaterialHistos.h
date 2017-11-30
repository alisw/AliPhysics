#ifndef AliAnalysisTaskMaterialHistos_cxx
#define AliAnalysisTaskMaterialHistos_cxx

#include "AliAnalysisTaskSE.h"
#include "AliConversionPhotonBase.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTreeStream.h"
#include "AliLog.h"
#include <vector>
#include "AliV0ReaderV1.h"
#include "AliConversionPhotonCuts.h"
#include "AliConvEventCuts.h"
#include "TList.h"
#include "TClonesArray.h"


using namespace std;


class AliAnalysisTaskMaterialHistos : public AliAnalysisTaskSE{

	public:

		AliAnalysisTaskMaterialHistos();
		AliAnalysisTaskMaterialHistos(const char *name);
		virtual ~AliAnalysisTaskMaterialHistos();

		virtual void   UserCreateOutputObjects();
		virtual void   UserExec(Option_t *option);
		virtual void   Terminate(Option_t *);
		void SetIsHeavyIon(Int_t flag)                                { fIsHeavyIon                 = flag    ;}
		void SetIsMC(Int_t isMC){fIsMC=isMC;}
		void SetV0Reader(AliV0ReaderV1 *v0Reader){fV0Reader=v0Reader;}
        void SetV0ReaderName(TString name){fV0ReaderName=name; return;}

    void SetEventCutList(Int_t nCuts, TList *CutArray)          { fnCuts                        = nCuts     ;
                                                                  fEventCutArray                = CutArray  ;}
    void SetConversionCutList(Int_t nCuts, TList *CutArray)     { fnCuts                        = nCuts     ;
                                                                  fConversionCutArray           = CutArray  ;}
 


		/* void SetDoTreesForMaterial(Bool_t flag)                              { fDoTreesForMaterial               = flag    ;} */
		/* void SetDoHistosForMaterial(Bool_t flag)                             { fDoHistosForMaterial              = flag    ;} */
	private:
		
		void ProcessPhotons();
		void ProcessMCPhotons();
        void FillMCTree(Int_t eventPos);
		Int_t CountTracks0914();
		Int_t CountTracks09();

		AliV0ReaderV1 		*fV0Reader;					// 
        TString              fV0ReaderName;
		TClonesArray 			*fConversionGammas; 		// Reconstructed Photons;
    TList*                            fGammaCandidates;                           //
    TList*            fConversionCutArray;
    TList*            fEventCutArray;         //
    TList**           fCutFolder;                                 //
    TList**           fESDList;                                   //
    TList**           fTrueList;                                  //
    TList**           fMCList;                                    //
		TList 						*fOutputList; 				//
		TList 						*fAllMCGammaList;			//
		TList 						*fAllMCConvGammaList;		//
		Float_t 					fPrimVtxZ;					//
		Int_t 						fNContrVtx;					//
		Int_t 						fNESDtracksEta09;			//
		Int_t 						fNESDtracksEta0914;			//
		Int_t 						fNESDtracksEta14;			//
		Float_t 					fGammaMCPt;					//
		Float_t 					fGammaMCTheta;				//
		Float_t 					fGammaMCConvPt;				//
		Float_t 					fGammaMCConvTheta;			//
		Float_t 					fGammaPt;					//
		Float_t 					fGammaTheta;				//
		Float_t 					fGammaChi2NDF;				//
		UChar_t 					fKind;						//
		Int_t 						fIsHeavyIon;			//
		Int_t 						fIsMC;					           	//
		AliVEvent*        fInputEvent;                                //
		AliMCEvent*       fMCEvent;				         //
		Int_t             fnCuts;                      //
		Int_t             fiCut;                       //
		TH1F**            hNEvents         ;           //!
		TH1F**            hNGoodESDTracks09;           //!
		TH1F**            hNGoodESDTracks14;           //!
		TH1F**            hNGoodESDTracks09_14;        //!
		TH2F**            hESDConversionMappingRPhi;      //!
		TH2F**            hESDConversionMappingRZ;        //!
		TH1F**            hESDConversionR;                //! 
		TH2F**            hESDConversionPtvsR;                //!   
		TH2F**            hESDConversionAsymP;            //!  
		TH1F**            hESDConversionMidPtR;           //!
		TH1F**            hESDConversionHighPtR;           //!
		TH1F**            hESDConversionEtaPR;           //!
		TH1F**            hESDConversionEtaNR;           //!
		TH3F**            hESDConversionRInBins;          //!
		TH3F**            hESDConversionPhiInBins;        //!
		TH1F**            hESDConversionEta;              //!  
		TH1F**            hESDConversionMidPtEta;         //!
		TH1F**            hESDConversionHighPtEta;         //!
		TH1F**            hESDConversionPt;              //!
		TH1F**            hESDConversionPt5cm;              //!
		TH1F**           	hESDConversionDCA;             //!
		TH1F**           	hESDConversionMidPtDCA;        //!
		TH1F**           	hESDConversionHighPtDCA;        //!
		TH1F**           	hESDConversionPsiPair;         //!
		TH1F**           	hESDConversionChi2;            //!
		TH1F**           	hESDConversionMass;            //!

 		TH2F**            hMCConversionMappingRPhi;      //!
		TH1F**            hMCConversionR;                //!  
		TH2F**            hMCConversionPtvsR;                //!  
		TH1F**            hMCConversionMidPtR;           //!  
		TH1F**            hMCConversionHighPtR;           //! 
		TH1F**            hMCConversionEtaPR;           //!   
		TH1F**            hMCConversionEtaNR;           //!   
		TH1F**            hMCConversionEta;              //!  
		TH1F**            hMCConversionMidPtEta;         //!
		TH1F**            hMCConversionHighPtEta;         //!
		TH1F**            hMCConversionPt;               //!
		TH1F**            hMCAllGammaPt;                 //!
		TH2F**            hMCTrueConversionMappingRPhi;      //!
		TH2F**            hMCTrueConversionMappingRZ;        //!
		TH1F**            hMCTrueConversionR;                //!  
		TH2F**            hMCTrueConversionPtvsR;                //!  
		TH1F**            hMCTrueConversionMidPtR;           //!  
		TH1F**            hMCTrueConversionHighPtR;           //!  
		TH1F**            hMCTrueConversionEtaPR;           //!  
		TH1F**            hMCTrueConversionEtaNR;           //!  
		TH1F**            hMCTrueConversionEta;              //!  
		TH1F**            hMCTrueConversionMidPtEta;         //!
		TH1F**            hMCTrueConversionHighPtEta;         //!
		TH1F**            hMCTrueConversionPt;               //!
		TH1F**            hMCTrueConversionPt5cm;               //!
		TH2F**            hMCTrueConversionAsymP;            //!
		TH1F**            hMCTrueConversionDCA;              //!
		TH1F**            hMCTrueConversionMidPtDCA;         //!
		TH1F**            hMCTrueConversionHighPtDCA;         //!
		TH1F**            hMCTrueConversionPsiPair;          //!
		TH1F**            hMCTrueConversionChi2;             //!
		TH1F**            hMCTrueConversionMass;             //!

		TH1F**            hMCTruePi0DalConversionR;                //!  
		TH1F**            hMCTruePi0DalConversionMidPtR;           //!  
		TH1F**            hMCTruePi0DalConversionHighPtR;           //!  
		TH1F**            hMCTruePi0DalConversionEta;              //!  
		TH1F**            hMCTruePi0DalConversionMidPtEta;         //!
		TH1F**            hMCTruePi0DalConversionHighPtEta;         //!
		TH1F**            hMCTruePi0DalConversionPt;              //!
		TH1F**            hMCTrueEtaDalConversionR;                //!  
		TH1F**            hMCTrueEtaDalConversionMidPtR;           //!  
		TH1F**            hMCTrueEtaDalConversionHighPtR;           //!  
		TH1F**            hMCTrueEtaDalConversionEta;              //!  
		TH1F**            hMCTrueEtaDalConversionMidPtEta;         //!
		TH1F**            hMCTrueEtaDalConversionHighPtEta;         //!
		TH1F**            hMCTrueEtaDalConversionPt;              //!
		TH1F**            hMCTrueCombConversionR;                //!  
		TH1F**            hMCTrueCombConversionMidPtR;           //!  
		TH1F**            hMCTrueCombConversionHighPtR;           //!  
		TH1F**            hMCTrueCombConversionEta;              //!  
		TH1F**            hMCTrueCombConversionMidPtEta;         //!
		TH1F**            hMCTrueCombConversionHighPtEta;         //!
		TH1F**            hMCTrueCombConversionPt;              //!
		//		Bool_t            fDoHistosForMaterial;     // flag for using Trees for Material Budget evaluation

		AliAnalysisTaskMaterialHistos(const AliAnalysisTaskMaterialHistos&); // not implemented
		AliAnalysisTaskMaterialHistos& operator=(const AliAnalysisTaskMaterialHistos&); // not implemented


        ClassDef(AliAnalysisTaskMaterialHistos, 10);
};

#endif
