#include "TH1.h"
#include "THnSparse.h"
#include "TList.h"
#include "AliESDEvent.h"
#include "AliESD.h"
#include "TLorentzVector.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliESDtrackCuts.h"
class AliESDtrackCuts;
class TChain;

class AliAnalysisTaskV0QA : public AliAnalysisTask {
 public:
  //  AliAnalysisTaskV0QA() : AliAnalysisTask(), fESD(0), fChain(0) {}
  AliAnalysisTaskV0QA(const char *name);
  virtual ~AliAnalysisTaskV0QA();
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void InspectListOfChargedParticles();
  void InspectListOfV0s();
  void FillHnSparseGamma();
  void FillHnSparseK0();
  void FillHnSparseL();
  void FillHnSparseAL();

  //  void getPID(AliESDtrack *esdTrack, Stat_t &fpid, Stat_t &fweight);


  //  AliESDtrackCuts* fEsdTrackCuts;           // Object containing the parameters of the esd track cuts
  //  void SetESDtrackCuts();




  Int_t GetTPCReference(Int_t label);




 private:



  AliESDEvent *fESD; //ESD object
  AliStack * stack;
  AliMCEventHandler* mctruth; 
 
  TChain * fChain;
  TList * fOutputContainer; // ! output data container

  THnSparse   *fSparseV0;
  THnSparse   *fSparseK0;
  THnSparse   *fSparseL;
  THnSparse   *fSparseAL;

  //////////////////////////////////

  Int_t nEv;

  Int_t nConvGamGeant;

  Int_t * gConvGamGeantIndex;
  Int_t * eNegConvGamGeantIndex;
  Int_t * ePosConvGamGeantIndex;
  Float_t * eNegConvGamGeantLength;
  Float_t * ePosConvGamGeantLength;


  Int_t * eNegConvGamSingleRecIndex;
  Int_t * ePosConvGamSingleRecIndex;
  Int_t * eNegConvGamV0RecIndex;
  Int_t * ePosConvGamV0RecIndex;
  Int_t * ConvGamV0RecIndexPos;
  Int_t * ConvGamV0RecIndexNeg;


  Int_t gDim;
  // Lambda
  Int_t nDecayLGeant;
  Int_t * lDecayLGeantIndex;
  Int_t * piNegDecayLGeantIndex;
  Int_t * pPosDecayLGeantIndex;
  Float_t * piNegDecayLGeantLength;
  Float_t * pPosDecayLGeantLength;

  Int_t * piNegDecayLSingleRecIndex;
  Int_t * pPosDecayLSingleRecIndex;
  Int_t * piNegDecayLV0RecIndex;
  Int_t * pPosDecayLV0RecIndex;
  Int_t * DecayLV0RecIndexPos;
  Int_t * DecayLV0RecIndexNeg;


  // AntiLambda
  Int_t nDecayALGeant;
  Int_t * alDecayALGeantIndex;
  Int_t * piPosDecayALGeantIndex;
  Int_t * apNegDecayALGeantIndex;
  Float_t * piPosDecayALGeantLength;
  Float_t * apNegDecayALGeantLength;

  Int_t * piPosDecayALSingleRecIndex;
  Int_t * apNegDecayALSingleRecIndex;
  Int_t * piPosDecayALV0RecIndex;
  Int_t * apNegDecayALV0RecIndex;
  Int_t * DecayALV0RecIndexPos;
  Int_t * DecayALV0RecIndexNeg;


  // K0S
  Int_t nDecayK0Geant;
  Int_t * K0DecayK0GeantIndex;
  Int_t * piNegDecayK0GeantIndex;
  Int_t * piPosDecayK0GeantIndex;
  Float_t * piNegDecayK0GeantLength;
  Float_t * piPosDecayK0GeantLength;

  Int_t * piNegDecayK0SingleRecIndex;
  Int_t * piPosDecayK0SingleRecIndex;
  Int_t * piNegDecayK0V0RecIndex;
  Int_t * piPosDecayK0V0RecIndex;
  Int_t * DecayK0V0RecIndexPos;
  Int_t * DecayK0V0RecIndexNeg;

  Int_t piPosK0Index;
  Int_t piNegK0Index;

  
  Int_t nTracksPrim;


  Int_t tpcRefit;
  Int_t itsRefit;
  Int_t trdRefit;
  Int_t trdOut;




  Double_t probabilityPos[AliPID::kSPECIES];
  Double_t probabilityNeg[AliPID::kSPECIES];

  Double_t * fValueL;
  Double_t * fValueAL;
  Double_t * fValueK0;
  Double_t * fValueV0;

  Double_t * xminV0;
  Double_t * xmaxV0;
  Int_t    * binsV0;
  Int_t    fDim;



  TObjArray* fRefTPC;
  int fLabelsTPC[100000];


  TClonesArray *clRefsN;
  TClonesArray *clRefsP;

  // MC variables

  AliAnalysisTaskV0QA (const AliAnalysisTaskV0QA & v0QA );
  AliAnalysisTaskV0QA & operator=(const AliAnalysisTaskV0QA & v0QA);



  ClassDef(AliAnalysisTaskV0QA, 1); // example of analysis
};
