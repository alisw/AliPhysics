#ifndef AliAnalysisTaskTrackQA_cxx
#define AliAnalysisTaskTrackQA_cxx

#include "AliAnalysisTaskSE.h"
#include "AliConversionPhotonBase.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTreeStream.h"
#include "AliLog.h"
#include <vector>
#include "AliV0ReaderV1.h"
#include "AliConvEventCuts.h"
#include "AliIdentifiedPrimarySelector.h"
//#include "AliPrimaryPionSelector.h"
//#include "AliPrimaryKaonSelector.h"
//#include "AliPrimaryProtonSelector.h"
//#include "AliPrimaryDeuteronSelector.h"
#include "TList.h"
#include "TClonesArray.h"


using namespace std;


class AliAnalysisTaskTrackQA : public AliAnalysisTaskSE{

 public:

  AliAnalysisTaskTrackQA();
  AliAnalysisTaskTrackQA(const char *name);
  virtual ~AliAnalysisTaskTrackQA();
  
  virtual void   UserCreateOutputObjects();
  virtual Bool_t Notify();
  virtual void   UserExec(Option_t *option);
  virtual void   SetLogBinningXTH2(TH2* histoRebin);
  virtual void   Terminate(Option_t *);
  void SetIsHeavyIon(Int_t flag)                                {fIsHeavyIon                 = flag;}
  void SetIsMC(Int_t isMC)                                      {fIsMC=isMC;}
  void SetV0Reader(AliV0ReaderV1 *v0Reader)                     {fV0Reader=v0Reader;}
  void SetV0ReaderName(TString name)                            {fV0ReaderName=name; return;}
//  void SetIdentifiedPrimarySelectorName(TString name)                        {fPionSelectorName=name; return;}
  void SetPionSelectorName(TString name)                        {fPionSelectorName=name; return;}
  void SetKaonSelectorName(TString name)                        {fKaonSelectorName=name; return;}
  void SetProtonSelectorName(TString name)                      {fProtonSelectorName=name; return;}
  void SetDeuteronSelectorName(TString name)                    {fDeuteronSelectorName=name; return;}

  //  void SetDoDeDxMaps(Int_t flag)                          { fDoDeDxMaps           = flag    ;}
  //  void SetDoMultWeights(Int_t flag)                          { fDoMultWeights     = flag    ;}

  void SetEventCutList(Int_t nCuts, TList *CutArray)            {fnCuts                      = nCuts;
    fEventCutArray              = CutArray;}
  //  void SetConversionCutList(Int_t nCuts, TList *CutArray)       {fnCuts                      = nCuts;
  //  fConversionCutArray         = CutArray;}
  void SetPionCutList(TList *CutArray){ fPionCutArray = CutArray;}
  void SetKaonCutList(TList *CutArray){ fKaonCutArray = CutArray;}
  void SetProtonCutList(TList *CutArray){ fProtonCutArray = CutArray;}
  void SetDeuteronCutList(TList *CutArray){ fDeuteronCutArray = CutArray;}



 private:
  void ProcessPionCandidates();
  void ProcessKaonCandidates();
  void ProcessProtonCandidates();
  void ProcessDeuteronCandidates();


  Int_t CountTracks08();

  
  AliV0ReaderV1*    fV0Reader;	                //
  TString           fV0ReaderName;
  AliIdentifiedPrimarySelector*     fPionSelector;                  //!<! primary charged pion selector, basic selection of pi+,pi-
  TString           fPionSelectorName; 
  AliIdentifiedPrimarySelector*     fKaonSelector;                  //!<! primary charged Kaon selector, basic selection of K+,K-
  TString           fKaonSelectorName; 
  AliIdentifiedPrimarySelector*   fProtonSelector;                //!<! primary charged proton selector, basic selection of p+,p-
  TString           fProtonSelectorName; 
  AliIdentifiedPrimarySelector* fDeuteronSelector;              //!<! primary charged deuteron selector, basic selection of d+,d-
  TString           fDeuteronSelectorName; 
  TList*            fEventCutArray;              ///<
  TList*            fPionCutArray;               ///<
  TList*            fKaonCutArray;               ///<
  TList*            fProtonCutArray;             ///<
  TList*            fDeuteronCutArray;           ///<
  vector<Int_t>     fSelectorNegPionIndex;                    //!<! array with pion indices for negative pions from fPionSelector
  vector<Int_t>     fSelectorPosPionIndex;                    //!<! array with pion indices for positive pions from fPionSelector
  vector<Int_t>     fSelectorNegKaonIndex;                    //!<! array with pion indices for negative pions from fPionSelector
  vector<Int_t>     fSelectorPosKaonIndex;                    //!<! array with pion indices for positive pions from fPionSelector
  vector<Int_t>     fSelectorNegProtonIndex;                  //!<! array with pion indices for negative pions from fPionSelector
  vector<Int_t>     fSelectorPosProtonIndex;                  //!<! array with pion indices for positive pions from fPionSelector
  vector<Int_t>     fSelectorNegDeuteronIndex;                //!<! array with pion indices for negative pions from fPionSelector
  vector<Int_t>     fSelectorPosDeuteronIndex;                //!<! array with pion indices for positive pions from fPionSelector
  TList*            fPosPionCandidates;                       //!<! good positive pion candidates
  TList*            fNegPionCandidates;                       //!<! good negative pion candidates
  TList*            fPosKaonCandidates;                       //!<! good positive pion candidates
  TList*            fNegKaonCandidates;                       //!<! good negative pion candidates
  TList*            fPosProtonCandidates;                       //!<! good positive pion candidates
  TList*            fNegProtonCandidates;                       //!<! good negative pion candidates
  TList*            fPosDeuteronCandidates;                       //!<! good positive pion candidates
  TList*            fNegDeuteronCandidates;                       //!<! good negative pion candidates
  TList**           fCutFolder;                   //
  TList**           fESDList;                     //
  TList**           fTrueList;                    //
  TList**           fMCList;                      //
  TList*            fOutputContainer;                                   //!<! output container
  Int_t 	    fNESDtracksEta08;	         //
  Int_t    	    fNContrVtx;		       //
  TH2F**            hESDPosPionDCAxy;         //!
  TH2F**            hESDNegPionDCAxy;         //!
  TH2F**            hESDPosKaonDCAxy;         //!
  TH2F**            hESDNegKaonDCAxy;         //!
  TH2F**            hESDPosProtonDCAxy;       //!
  TH2F**            hESDNegProtonDCAxy;       //!
  TH2F**            hESDPosDeuteronDCAxy;     //!
  TH2F**            hESDNegDeuteronDCAxy;     //!
  TH2F**            hESDPosPionDCAz;         //!
  TH2F**            hESDNegPionDCAz;         //!
  TH2F**            hESDPosKaonDCAz;         //!
  TH2F**            hESDNegKaonDCAz;         //!
  TH2F**            hESDPosProtonDCAz;       //!
  TH2F**            hESDNegProtonDCAz;       //!
  TH2F**            hESDPosDeuteronDCAz;     //!
  TH2F**            hESDNegDeuteronDCAz;     //!
  Int_t 	    fIsHeavyIon;		//
  Int_t 	    fIsMC;		        //
  AliVEvent*        fInputEvent;                  //
  AliMCEvent*       fMCEvent;	                //
  Int_t             fnCuts;                       //
  Int_t             fiCut;                        //
  TH1F**            hNEvents;                     //!
  TH1F**            hNGoodESDTracksEta08;         //!
 
 


  AliAnalysisTaskTrackQA(const AliAnalysisTaskTrackQA&); // not implemented
  AliAnalysisTaskTrackQA& operator=(const AliAnalysisTaskTrackQA&); // not implemented


  ClassDef(AliAnalysisTaskTrackQA, 1);
};

#endif
