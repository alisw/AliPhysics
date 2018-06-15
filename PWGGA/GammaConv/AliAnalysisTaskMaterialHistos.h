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
        virtual void   SetLogBinningXTH2(TH2* histoRebin);
		virtual void   Terminate(Option_t *);
		void SetIsHeavyIon(Int_t flag)                                {fIsHeavyIon                 = flag;}
		void SetIsMC(Int_t isMC)                                      {fIsMC=isMC;}
		void SetV0Reader(AliV0ReaderV1 *v0Reader)                     {fV0Reader=v0Reader;}
        void SetV0ReaderName(TString name)                            {fV0ReaderName=name; return;}

        void SetEventCutList(Int_t nCuts, TList *CutArray)            {fnCuts                      = nCuts;
                                                                       fEventCutArray              = CutArray;}
        void SetConversionCutList(Int_t nCuts, TList *CutArray)       {fnCuts                      = nCuts;
                                                                       fConversionCutArray         = CutArray;}



		/* void SetDoTreesForMaterial(Bool_t flag)                              { fDoTreesForMaterial               = flag    ;} */
		/* void SetDoHistosForMaterial(Bool_t flag)                             { fDoHistosForMaterial              = flag    ;} */
	private:

		void ProcessPhotons();
		void ProcessMCPhotons();
        void FillMCHistograms(Int_t eventPos);
		Int_t CountTracks09();
		Int_t CountTracks0914();

		AliV0ReaderV1 	  *fV0Reader;	                //
        TString           fV0ReaderName;
		TClonesArray 	  *fConversionGammas; 		    // Reconstructed Photons;
        TList*            fGammaCandidates;             //
        TList*            fConversionCutArray;
        TList*            fEventCutArray;               //
        TList**           fCutFolder;                   //
        TList**           fESDList;                     //
        TList**           fTrueList;                    //
        TList**           fMCList;                      //
		TList 			  *fOutputList; 	 			//
		TList 			  *fAllMCGammaList;	   		    //
		TList 			  *fAllMCConvGammaList;		    //
		Float_t 		  fPrimVtxZ;					//
		Int_t 			  fNContrVtx;					//
		Int_t 			  fNESDtracksEta09;		      	//
		Int_t 			  fNESDtracksEta0914;			//
		Int_t 			  fNESDtracksEta14;		      	//
		Float_t 		  fGammaMCPt;					//
		Float_t 		  fGammaMCTheta;				//
		Float_t 		  fGammaMCConvPt;				//
		Float_t 		  fGammaMCConvTheta;			//
		Float_t 		  fGammaPt;				     	//
		Float_t 		  fGammaTheta;		            //
		Float_t 		  fGammaChi2NDF;				//
		UChar_t 		  fKind;						//
		Int_t 			  fIsHeavyIon;		         	//
		Int_t 	          fIsMC;		               	//
		AliVEvent*        fInputEvent;                  //
		AliMCEvent*       fMCEvent;				        //
		Int_t             fnCuts;                       //
		Int_t             fiCut;                        //
		TH1F**            hNEvents;                     //!
		TH1F**            hNGoodESDTracksEta09;         //!
		TH1F**            hNGoodESDTracksEta14;         //!
		TH1F**            hNGoodESDTracksEta09_14;      //!
		TH2F**            hESDConversionRPhi;           //!
		TH2F**            hESDConversionRZ;             //!
		TH2F**            hESDConversionRPt;            //!
		TH2F**            hESDConversionREta;           //!
		TH1F**            hESDConversionDCA;            //!
		TH1F**            hESDConversionPsiPair;        //!
		TH1F**            hESDConversionChi2;           //!
		TH1F**            hESDConversionMass;           //!
		TH1F**            hESDConversionRRejSmall;      //!
		TH1F**            hESDConversionRRejLarge;      //!
		TH2F**            hESDConversionAsymP;          //!  
        TH2F**            hElectronRdEdx;               //!
        TH2F**            hElectronRNSigmadEdx;         //!
        TH2F**            hPositronRdEdx;               //!
        TH2F**            hPositronRNSigmadEdx;         //!

 		TH2F**            hMCConversionRPhi;            //!
		TH2F**            hMCConversionRPt;             //!
		TH2F**            hMCConversionREta;            //!
 		TH1F**            hMCConversionRRejSmall;       //!
 		TH1F**            hMCConversionRRejLarge;       //!
		TH1F**            hMCAllGammaPt;                //!


		TH2F**            hMCTrueConversionRPhi;        //!
		TH2F**            hMCTrueConversionRZ;          //!
		TH2F**            hMCTrueConversionRPt;         //!
		TH2F**            hMCTrueConversionREta;        //!
		TH1F**            hMCTrueConversionDCA;         //!
		TH1F**            hMCTrueConversionPsiPair;     //!
		TH1F**            hMCTrueConversionChi2;        //!
		TH1F**            hMCTrueConversionMass;        //!
		TH2F**            hMCTrueConversionAsymP;          //!
 		TH1F**            hMCTrueConversionRRejSmall;   //!
 		TH1F**            hMCTrueConversionRRejLarge;   //!

		TH2F**            hMCTruePi0DalConversionRPt;   //!
		TH1F**            hMCTruePi0DalConversionEta;   //!
		TH2F**            hMCTrueEtaDalConversionRPt;   //!
		TH1F**            hMCTrueEtaDalConversionEta;   //!
		TH2F**            hMCTrueCombinatorialConversionRPt;     //!
		TH1F**            hMCTrueCombinatorialConversionEta;     //!
		//Bool_t            fDoHistosForMaterial;             // flag for using Trees for Material Budget evaluation

		AliAnalysisTaskMaterialHistos(const AliAnalysisTaskMaterialHistos&); // not implemented
		AliAnalysisTaskMaterialHistos& operator=(const AliAnalysisTaskMaterialHistos&); // not implemented


        ClassDef(AliAnalysisTaskMaterialHistos, 13);
};

#endif
