#ifndef AliAnalysisTaskV0QA_h
#define AliAnalysisTaskV0QA_h
//----------------------------------
// Class to check the V0 method efficiency for 
// Author A. Marin   revision 18/10/2009
//----------------------------------

#include "THnSparse.h"
#include "TList.h"
#include "AliPID.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

class TH1;
class THnSparse;
class TList;
class AliESDEvent;
class AliESD;
class TLorentzVector;
class AliAnalysisManager;
class AliAnalysisDataContainer;
class AliESDtrackCuts;
class AliMCEventHandler;
class AliStack;
class TChain;

class AliAnalysisTaskV0QA : public AliAnalysisTaskSE {
 public:
  //  AliAnalysisTaskV0QA() : AliAnalysisTask(), fESD(0), fChain(0) {}
    AliAnalysisTaskV0QA();
    AliAnalysisTaskV0QA(const char *name);
  
  virtual ~AliAnalysisTaskV0QA();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
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



  AliESDEvent *fESD;                // ESD object
  AliStack * fStack;                // The MC Stack
  AliMCEventHandler* fMCtruth;      // The mc info
 
  TChain * fChain;                  // Input chian
  TList * fOutputContainer; // ! output data container

  THnSparse   *fSparseV0;   // THnSparse with Gamma info 
  THnSparse   *fSparseK0;   // THnSparse with K0 info
  THnSparse   *fSparseL;    // THnSparse with L info
  THnSparse   *fSparseAL;   // THnSparse with antiL info

  //////////////////////////////////

  Int_t fnEv;               // Number of event to analyse
  Int_t fgDim;              // Dimension of the THnSparse

  Int_t fnConvGamGeant;     // number of conversions in mc 

  Int_t * fgConvGamGeantIndex; //[fgDim] index of conversions in mc 
  Int_t * feNegConvGamGeantIndex; //[fgDim] index of e- from conversions in mc 
  Int_t * fePosConvGamGeantIndex; //[fgDim] index of e+ from conversions in mc 
  Float_t * feNegConvGamGeantLength; //[fgDim]  length of the e- from conv 
  Float_t * fePosConvGamGeantLength; //[fgDim]  length of the e+ from conv 


  Int_t * feNegConvGamSingleRecIndex; //[fgDim] index of e- from conversions reconstructed single
  Int_t * fePosConvGamSingleRecIndex; //[fgDim] index of e+ from conversions reconstructed single
  Int_t * feNegConvGamV0RecIndex;     //[fgDim] index of e- from conversions reconstructed in V0
  Int_t * fePosConvGamV0RecIndex;     //[fgDim] index of e- from conversions reconstructed in V0
  Int_t * fConvGamV0RecIndexPos;      //[fgDim] index of V0 from conversions reconstructed
  Int_t * fConvGamV0RecIndexNeg;      //[fgDim] index of V0 from conversions reconstructed



   // Lambda  
  Int_t fnDecayLGeant;                // number of Lambda in mc 
  Int_t * flDecayLGeantIndex;         //[fgDim] index of Lambda in MC
  Int_t * fpiNegDecayLGeantIndex;     //[fgDim] index of pi- from L in MC
  Int_t * fpPosDecayLGeantIndex;      //[fgDim] index of proton from L in MC 
  Float_t * fpiNegDecayLGeantLength;  //[fgDim] length of the pi- from MC
  Float_t * fpPosDecayLGeantLength;   //[fgDim] length of the proton from MC

  Int_t * fpiNegDecayLSingleRecIndex; //[fgDim] index of pi- from L reconstr. single
  Int_t * fpPosDecayLSingleRecIndex;  //[fgDim] index of proton from L reconstr. single
  Int_t * fpiNegDecayLV0RecIndex;     //[fgDim] index of pi- from L reconstr. V0
  Int_t * fpPosDecayLV0RecIndex;      //[fgDim] index of proton from L reconstr. V0
  Int_t * fDecayLV0RecIndexPos;       //[fgDim] index of pi- from L reconstr. V0
  Int_t * fDecayLV0RecIndexNeg;       //[fgDim] index of proton from L reconstr.


  // AntiLambda
  Int_t fnDecayALGeant;               // number of Lambdabar in mc 
  Int_t * falDecayALGeantIndex;       //[fgDim] index of Lambdabar in MC
  Int_t * fpiPosDecayALGeantIndex;    //[fgDim] index of pi+ from AL in MC
  Int_t * fapNegDecayALGeantIndex;    //[fgDim] index of antiproton from AL in MC
  Float_t * fpiPosDecayALGeantLength; //[fgDim] Length of pi+ in MC 
  Float_t * fapNegDecayALGeantLength; //[fgDim] Length of antiproton in MC 

  Int_t * fpiPosDecayALSingleRecIndex; //[fgDim] index of pi+ from AL reconstr. single
  Int_t * fapNegDecayALSingleRecIndex; //[fgDim] index of antiproton from AL reconstr. single
  Int_t * fpiPosDecayALV0RecIndex;     //[fgDim] index of pi+ from AL reconstr. V0 
  Int_t * fapNegDecayALV0RecIndex;     //[fgDim] index of antiproton from AL reconstr. V0  
  Int_t * fDecayALV0RecIndexPos;       //[fgDim] index of pi+ V0
  Int_t * fDecayALV0RecIndexNeg;       //[fgDim] index of antiproton V0


  // K0S
  Int_t fnDecayK0Geant;                // number of K0s in mc  
  Int_t * fK0DecayK0GeantIndex;        //[fgDim] index of K0S in MC
  Int_t * fpiNegDecayK0GeantIndex;     //[fgDim] index of pi- from K0s in MC
  Int_t * fpiPosDecayK0GeantIndex;     //[fgDim] index of pi+ from K0s in MC
  Float_t * fpiNegDecayK0GeantLength;  //[fgDim] length of the pi- from K0s in MC
  Float_t * fpiPosDecayK0GeantLength;  //[fgDim] length of the pi+ from K0s in MC

  Int_t * fpiNegDecayK0SingleRecIndex;  //[fgDim] index of pi- from K0S reconstr. single
  Int_t * fpiPosDecayK0SingleRecIndex;  //[fgDim] index of pi+ from K0S reconstr. single 
  Int_t * fpiNegDecayK0V0RecIndex;      //[fgDim] index of pi- from K0S reconstr. V0    
  Int_t * fpiPosDecayK0V0RecIndex;      //[fgDim] index of pi+ from K0S reconstr. V0
  Int_t * fDecayK0V0RecIndexPos;        //[fgDim] index of pi+ V0 
  Int_t * fDecayK0V0RecIndexNeg;        //[fgDim] index of pi- V0  

  Int_t fpiPosK0Index;    //      
  Int_t fpiNegK0Index;    //

  
  Int_t fnTracksPrim;    // number of primary tracks contributing to vertex


  Int_t ftpcRefit;    // tpcRefit condition
  Int_t fitsRefit;    // itsRefit condition
  Int_t ftrdRefit;    // trdRefit condition
  Int_t ftrdOut;      // trdOut condition




  Double_t fprobabilityPos[AliPID::kSPECIES];
  Double_t fprobabilityNeg[AliPID::kSPECIES];

  Int_t    fDim;        // number of dimensions THnSparse
  Double_t * fValueL;   //[fDim] values to THnSparse for L
  Double_t * fValueAL;  //[fDim] values to THnSparse for AL
  Double_t * fValueK0;  //[fDim] values to THnSparse for K0
  Double_t * fValueV0;  //[fDim] values to THnSparse for Gamma

  Double_t * fxminV0;   //[fDim] min value to THnSparse
  Double_t * fxmaxV0;   //[fDim] max value to THnSparse
  Int_t    * fbinsV0;   //[fDim] number of bins to THnSparse
  Int_t	 fCentralityC;  // centrality



  TObjArray* fRefTPC;   // references for the tpc
  int fLabelsTPC[100000]; // labels for the tpc


  TClonesArray *fclRefsN;  // negative references for the tpc
  TClonesArray *fclRefsP;  // positive references for the tpc

  // MC variables

  AliAnalysisTaskV0QA (const AliAnalysisTaskV0QA & v0QA );
  AliAnalysisTaskV0QA & operator=(const AliAnalysisTaskV0QA & v0QA);



  ClassDef(AliAnalysisTaskV0QA, 1); // example of analysis
};


#endif
