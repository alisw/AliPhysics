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
  virtual Bool_t Notify();
  virtual void   UserExec(Option_t *option);
  virtual void   SetLogBinningXTH2(TH2* histoRebin);
  virtual void   Terminate(Option_t *);
  void SetIsHeavyIon(Int_t flag)                                {fIsHeavyIon                 = flag;}
  void SetIsMC(Int_t isMC)                                      {fIsMC=isMC;}
  void SetV0Reader(AliV0ReaderV1 *v0Reader)                     {fV0Reader=v0Reader;}
  void SetV0ReaderName(TString name)                            {fV0ReaderName=name; return;}
  void SetDoDeDxMaps(Int_t flag)                          { fDoDeDxMaps           = flag    ;}
  void SetDoMultWeights(Int_t flag)                          { fDoMultWeights     = flag    ;}
  
  void SetEventCutList(Int_t nCuts, TList *CutArray)            {fnCuts                      = nCuts;
    fEventCutArray              = CutArray;}
  void SetConversionCutList(Int_t nCuts, TList *CutArray)       {fnCuts                      = nCuts;
    fConversionCutArray         = CutArray;}
  void SetDoMaterialBudgetWeightingOfGammasForTrueMesons(Bool_t flag) {fDoMaterialBudgetWeightingOfGammasForTrueMesons = flag;}
  void SetDoSelectBCnumbers(Int_t flag )                        { fDoSelectBCNumber     = flag;}

  /* void SetDoTreesForMaterial(Bool_t flag)                              { fDoTreesForMaterial               = flag    ;} */
  /* void SetDoHistosForMaterial(Bool_t flag)                             { fDoHistosForMaterial              = flag    ;} */
 private:
  
  void ProcessPhotons();
  void ProcessMCPhotons();
  void FillMCHistograms(Int_t eventPos);
  Int_t CountTracks08();
  Int_t CountTracks08pt200();
  Int_t CountTracks08pt300();
  Int_t CountTracks08pt400();
  Int_t CountTracks0814();
  void DoSelectBCNumbers();
  
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
  TList**           fDeDxMapList;                 //
  TList 			  *fOutputList; 	 			//
  TList 			  *fAllMCGammaList;	//
  TList 			  *fAllMCConvGammaList;	//
  Float_t 		  fPrimVtxZ;		//
  Int_t 			  fNContrVtx;		//
  Int_t 			  fNESDtracksEta08;	//
  Int_t 			  fNESDtracksEta08pt200;	//
  Int_t 			  fNESDtracksEta08pt300;	//
  Int_t 			  fNESDtracksEta08pt400;	//
  Int_t 			  fNESDtracksEta0814;	//
  Int_t 			  fNESDtracksEta14;	//
  Float_t 		  fGammaMCPt;		//
  Float_t 		  fGammaMCTheta;	//
  Float_t 		  fGammaMCConvPt;	//
  Float_t 		  fGammaMCConvTheta;	//
  Float_t 		  fGammaPt;		//
  Float_t 		  fGammaTheta;		//
  Float_t 		  fGammaChi2NDF;	//
  UChar_t 		  fKind;		//
  Int_t 			  fIsHeavyIon;		//
  Int_t 	          fIsMC;		        //
  AliVEvent*        fInputEvent;                  //
  AliMCEvent*       fMCEvent;	                //
  Int_t             fnCuts;                       //
  Int_t             fiCut;                        //
  Int_t             fDoDeDxMaps;                  //
  Int_t             fDoMultWeights;               //
  Double_t          fWeightMultMC;                //
  TH1F**            hNEvents;                     //!
  TH1F**            hBCNumber;                    //!
  TH1F**            hBCNumberSelected;            //!
  TH1F**            hNGoodESDTracksEta08;         //!
  TH1F**            hNGoodESDTracksWeightedEta08; //!
  TH1F**            hNGoodESDTracksEta08pt200;         //!
  TH1F**            hNGoodESDTracksWeightedEta08pt200; //!
  TH1F**            hNGoodESDTracksEta08pt300;         //!
  TH1F**            hNGoodESDTracksWeightedEta08pt300; //!
  TH1F**            hNGoodESDTracksEta08pt400;         //!
  TH1F**            hNGoodESDTracksWeightedEta08pt400; //!
  TH1F**            hNGoodESDTracksEta14;         //!
  TH1F**            hNGoodESDTracksEta08_14;      //!
  TH1F**            fHistoNV0Tracks;              //!
  TH1F**            fHistoNV0TracksWeighted;      //!
  
  TH2F**            hESDConversionRPhi;           //!
  TH2F**            hESDConversionRPhiFromConv;   //!
  TH2F**            hESDConversionRZ;             //!
  TH2F**            hESDConversionRPt;            //!
  TH2F**            hESDConversionWOWeightRPt;    //!
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
  TH2F**            hMCConversionRPhiFromConv;    //!
  TH2F**            hMCConversionRPt;             //!
  TH2F**            hMCConversionWOWeightRPt;             //!
  TH2F**            hMCConversionREta;            //!
  TH1F**            hMCConversionRRejSmall;       //!
  TH1F**            hMCConversionRRejLarge;       //!
  TH1F**            hMCAllGammaPt;                //!
  TH1F**            hMCAllGammaWOWeightPt;        //!
  TH2F**            hMCAllSecondaryGammaPt;       //!
  TH3F**            hMCSecondaryConvGammaPtR;      //!
  
  
  TH2F**            hMCTrueConversionRPhi;        //!
  TH2F**            hMCTrueConversionRPhiFromConv;//!
  TH2F**            hMCTrueConversionRZ;          //!
  TH2F**            hMCTrueConversionRPt;         //!
  TH2F**            hMCTrueConversionWOWeightRPt;         //!
  TH2F**            hMCTrueConversionRPtMCRPt;    //!
  TH2F**            hMCTrueConversionWOWeightRPtMCRPt;    //!
  TH2F**            hMCTrueConversionREta;        //!
  TH1F**            hMCTrueConversionDCA;         //!
  TH1F**            hMCTrueConversionPsiPair;     //!
  TH1F**            hMCTrueConversionChi2;        //!
  TH1F**            hMCTrueConversionMass;        //!
  TH2F**            hMCTrueConversionAsymP;          //!
  TH1F**            hMCTrueConversionRRejSmall;   //!
  TH1F**            hMCTrueConversionRRejLarge;   //!
  TH2F**            hMCTruePrimConversionRPt;     //!
  TH2F**            hMCTruePrimConversionWOWeightRPt;     //!
  TH2F**            hMCTrueSecConversionRPt;      //!
  TH3F**            hMCTrueSecondaryConvGammaRPt;//!
  TH3F**            hMCTrueSecondaryConvGammaMCRPt;//!
  TH2F**            hMCTruePi0DalConversionRPt;   //!
  TH1F**            hMCTruePi0DalConversionEta;   //!
  TH2F**            hMCTrueEtaDalConversionRPt;   //!
  TH1F**            hMCTrueEtaDalConversionEta;   //!
  TH2F**            hMCTrueCombinatorialConversionRPt;     //!
  TH1F**            hMCTrueCombinatorialConversionEta;     //!
  TH3F**            hPositrondEdxMapsR0;         //!
  TH3F**            hElectrondEdxMapsR0;         //!
  TH3F**            hPositrondEdxMapsR1;         //!
  TH3F**            hElectrondEdxMapsR1;         //!
  TH3F**            hPositrondEdxMapsR2;         //!
  TH3F**            hElectrondEdxMapsR2;         //!
  TH3F**            hPositrondEdxMapsR3;         //!
  TH3F**            hElectrondEdxMapsR3;         //!
  Bool_t            fDoMaterialBudgetWeightingOfGammasForTrueMesons;
  //Bool_t            fDoHistosForMaterial;             // flag for using Trees for Material Budget evaluation
  Bool_t            fDoSelectBCNumber;
  UShort_t          fBCNumber;
  Int_t             fRunNumber;


  AliAnalysisTaskMaterialHistos(const AliAnalysisTaskMaterialHistos&); // not implemented
  AliAnalysisTaskMaterialHistos& operator=(const AliAnalysisTaskMaterialHistos&); // not implemented
  
  
  ClassDef(AliAnalysisTaskMaterialHistos, 25);
};

#endif
