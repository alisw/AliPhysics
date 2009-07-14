#define AliAnalysisTaskV0QA_cxx
#include "Riostream.h"
#include "TROOT.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"

#include "TCanvas.h"
#include "TSystem.h"
#include "TLorentzVector.h"

#include "AliAnalysisTask.h"
#include "AliTrackReference.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliVertexerTracks.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "AliESD.h"
#include "AliLog.h"

#include "AliAnalysisTaskV0QA.h"

ClassImp(AliAnalysisTaskV0QA)

//________________________________________________________________________
AliAnalysisTaskV0QA::AliAnalysisTaskV0QA(const char *name) :AliAnalysisTask(name,""), 
fESD(0), 
stack(0),
mctruth(0),
fChain(0),
fOutputContainer(0),
fSparseV0(0),
fSparseK0(0),
fSparseL(0),
fSparseAL(0),
nEv(0),
nConvGamGeant(-1),
gConvGamGeantIndex(0),
eNegConvGamGeantIndex(0),
ePosConvGamGeantIndex(0),
eNegConvGamGeantLength(0),
ePosConvGamGeantLength(0),
eNegConvGamSingleRecIndex(0),
ePosConvGamSingleRecIndex(0),
eNegConvGamV0RecIndex(0),
ePosConvGamV0RecIndex(0),
ConvGamV0RecIndexPos(0),
ConvGamV0RecIndexNeg(0),
gDim(50),
nDecayLGeant(-1),
lDecayLGeantIndex(0),
piNegDecayLGeantIndex(0),
pPosDecayLGeantIndex(0),
piNegDecayLGeantLength(0),
pPosDecayLGeantLength(0),
piNegDecayLSingleRecIndex(0),
pPosDecayLSingleRecIndex(0),
piNegDecayLV0RecIndex(0),
pPosDecayLV0RecIndex(0),
DecayLV0RecIndexPos(0),
DecayLV0RecIndexNeg(0),
nDecayALGeant(-1),
alDecayALGeantIndex(0),
piPosDecayALGeantIndex(0),
apNegDecayALGeantIndex(0),
piPosDecayALGeantLength(0),
apNegDecayALGeantLength(0),
piPosDecayALSingleRecIndex(0),
apNegDecayALSingleRecIndex(0),
piPosDecayALV0RecIndex(0),
apNegDecayALV0RecIndex(0),
DecayALV0RecIndexPos(0),
DecayALV0RecIndexNeg(0),
nDecayK0Geant(-1),
K0DecayK0GeantIndex(0),
piNegDecayK0GeantIndex(0),
piPosDecayK0GeantIndex(0),
piNegDecayK0GeantLength(0),
piPosDecayK0GeantLength(0),
piNegDecayK0SingleRecIndex(0),
piPosDecayK0SingleRecIndex(0),
piNegDecayK0V0RecIndex(0),
piPosDecayK0V0RecIndex(0),
DecayK0V0RecIndexPos(0),
DecayK0V0RecIndexNeg(0),
piPosK0Index(-1),
piNegK0Index(-1),
nTracksPrim(-1),
tpcRefit(0),
itsRefit(0),
trdRefit(0),
trdOut(0),
fValueL(0),
fValueAL(0),
fValueK0(0),
fValueV0(0),
xminV0(0),
xmaxV0(0),
binsV0(0),
fDim(37),
fRefTPC(0),
clRefsN(0),
clRefsP(0)

 {
  // Constructor.
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  //  DefineOutput(0,TObjArray::Class());
  DefineOutput(0,TList::Class());

  // Reconstructed arrays

  nEv=0;
  fDim=37; 
  fValueK0 = new Double_t[fDim];
  fValueL = new Double_t[fDim];
  fValueAL = new Double_t[fDim];
  fValueV0 = new Double_t[fDim];
  xminV0 = new Double_t[fDim];
  xmaxV0 = new Double_t[fDim];
  binsV0 = new Int_t[fDim];


  gDim=50; 
  gConvGamGeantIndex = new Int_t[gDim];
  eNegConvGamGeantIndex = new Int_t[gDim];
  ePosConvGamGeantIndex = new Int_t[gDim];
  eNegConvGamGeantLength = new Float_t[gDim];
  ePosConvGamGeantLength = new Float_t[gDim];

  eNegConvGamSingleRecIndex = new Int_t[gDim];
  ePosConvGamSingleRecIndex = new Int_t[gDim];

  eNegConvGamV0RecIndex = new Int_t[gDim];
  ePosConvGamV0RecIndex = new Int_t[gDim];

  ConvGamV0RecIndexPos = new Int_t[gDim];
  ConvGamV0RecIndexNeg = new Int_t[gDim];

  // Lambda to proton pi-
  lDecayLGeantIndex = new Int_t[gDim];
  piNegDecayLGeantIndex = new Int_t[gDim];
  pPosDecayLGeantIndex = new Int_t[gDim];
  piNegDecayLGeantLength = new Float_t[gDim];
  pPosDecayLGeantLength = new Float_t[gDim];

  piNegDecayLSingleRecIndex = new Int_t[gDim];
  pPosDecayLSingleRecIndex = new Int_t[gDim];

  piNegDecayLV0RecIndex = new Int_t[gDim];
  pPosDecayLV0RecIndex = new Int_t[gDim];

  DecayLV0RecIndexPos = new Int_t[gDim];
  DecayLV0RecIndexNeg = new Int_t[gDim];

  //K0S to pi+ pi-
  K0DecayK0GeantIndex = new Int_t[gDim];
  piNegDecayK0GeantIndex = new Int_t[gDim];
  piPosDecayK0GeantIndex = new Int_t[gDim];
  piNegDecayK0GeantLength = new Float_t[gDim];
  piPosDecayK0GeantLength = new Float_t[gDim];

  piNegDecayK0SingleRecIndex = new Int_t[gDim];
  piPosDecayK0SingleRecIndex = new Int_t[gDim];

  piNegDecayK0V0RecIndex = new Int_t[gDim];
  piPosDecayK0V0RecIndex = new Int_t[gDim];

  DecayK0V0RecIndexPos = new Int_t[gDim];
  DecayK0V0RecIndexNeg = new Int_t[gDim];

  //Antilambda to antiproton piplus
  alDecayALGeantIndex = new Int_t[gDim];
  piPosDecayALGeantIndex = new Int_t[gDim];
  apNegDecayALGeantIndex = new Int_t[gDim];
  piPosDecayALGeantLength = new Float_t[gDim];
  apNegDecayALGeantLength = new Float_t[gDim];

  piPosDecayALSingleRecIndex = new Int_t[gDim];
  apNegDecayALSingleRecIndex = new Int_t[gDim];

  piPosDecayALV0RecIndex = new Int_t[gDim];
  apNegDecayALV0RecIndex = new Int_t[gDim];

  DecayALV0RecIndexPos = new Int_t[gDim];
  DecayALV0RecIndexNeg = new Int_t[gDim];


  //////////////////////////////////////////////
  clRefsP = new TClonesArray("AliTrackReference");
  clRefsN = new TClonesArray("AliTrackReference");

  //  SetESDtrackCuts();


  AliLog::SetGlobalLogLevel(AliLog::kError);

}

//________________________________________________________________________
void AliAnalysisTaskV0QA::ConnectInputData(Option_t *) {
  printf("   ConnectInputData %s\n", GetName());


  
  AliESDInputHandler* esdH = (AliESDInputHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  fESD = esdH->GetEvent();
  
  fChain = (TChain*)GetInputData(0);

}

//_____________________________________________________
AliAnalysisTaskV0QA::~AliAnalysisTaskV0QA()
{
  // Remove all pointers


  delete [] clRefsP;
  delete [] clRefsN;


  delete [] fValueK0;
  delete [] fValueL;
  delete [] fValueAL;
  delete [] fValueV0;
  delete [] binsV0;
  delete [] xminV0;
  delete [] xmaxV0;

  delete [] gConvGamGeantIndex;
  delete [] eNegConvGamGeantIndex;
  delete [] ePosConvGamGeantIndex;

  delete [] eNegConvGamSingleRecIndex;
  delete [] ePosConvGamSingleRecIndex;

  delete [] eNegConvGamV0RecIndex;
  delete [] ePosConvGamV0RecIndex;
  delete [] ConvGamV0RecIndexPos;
  delete [] ConvGamV0RecIndexNeg;

  delete [] lDecayLGeantIndex;
  delete [] piNegDecayLGeantIndex;
  delete [] pPosDecayLGeantIndex;

  delete [] piNegDecayLGeantLength;
  delete [] pPosDecayLGeantLength;
  delete [] piNegDecayLSingleRecIndex;
  delete [] pPosDecayLSingleRecIndex;

  delete [] piNegDecayLV0RecIndex;
  delete [] pPosDecayLV0RecIndex;
  delete [] DecayLV0RecIndexPos;
  delete [] DecayLV0RecIndexNeg;

  delete [] alDecayALGeantIndex;
  delete [] piPosDecayALGeantIndex;
  delete [] apNegDecayALGeantIndex;

  delete [] piPosDecayALGeantLength;
  delete [] apNegDecayALGeantLength;
  delete [] piPosDecayALSingleRecIndex;
  delete [] apNegDecayALSingleRecIndex;

  delete [] piPosDecayALV0RecIndex;
  delete [] apNegDecayALV0RecIndex;
  delete [] DecayALV0RecIndexPos;
  delete [] DecayALV0RecIndexNeg;


  delete [] piNegDecayK0GeantIndex;
  delete [] piPosDecayK0GeantIndex;

  delete [] piNegDecayK0GeantLength;
  delete [] piPosDecayK0GeantLength;
  delete [] piNegDecayK0SingleRecIndex;
  delete [] piPosDecayK0SingleRecIndex;

  delete [] piNegDecayK0V0RecIndex;
  delete [] piPosDecayK0V0RecIndex;

  delete [] DecayK0V0RecIndexPos;
  delete [] DecayK0V0RecIndexNeg;

}


//________________________________________________________________________
void AliAnalysisTaskV0QA::CreateOutputObjects() {

  for(Int_t d=0;d<fDim;d++){
    binsV0[d]=70;
  }
  xminV0[0]=   0;     // 1/sqrt(pt) Gamma geant
  xmaxV0[0]=   8;


  xminV0[1]=-2.5;   // eta Gamma Geant
  xmaxV0[1]= 1.5;


  xminV0[2]=-2*TMath::Pi();   // phi Gamma geant
  xmaxV0[2]= TMath::Pi();


  xminV0[3]=   0;     // r geant
  xmaxV0[3]= 200;


  xminV0[4]=-250;   // z geant
  xmaxV0[4]= 250;


  xminV0[5]=   0;     // 1/sqrt(pt) Geant Pos
  xmaxV0[5]=   8;


  xminV0[6]=-2.5;   // eta geant Pos
  xmaxV0[6]= 1.5;


  xminV0[7]=-2*TMath::Pi();   // phi Geant Pos
  xmaxV0[7]= TMath::Pi();

  
  xminV0[8]=0;   // Track Length TPC Geant Pos
  xmaxV0[8]= 200;


  xminV0[9]=   0;     // 1/sqrt(pt) Geant Neg
  xmaxV0[9]=   8;


  xminV0[10]=-2.5;   // eta Geant Neg
  xmaxV0[10]= 1.5;


  xminV0[11]=-2*TMath::Pi();   // phi Geant Neg
  xmaxV0[11]= TMath::Pi();


  xminV0[12]=0;   // Track Length TPC Geant Neg
  xmaxV0[12]= 200;



  //-----------Rec single variables

  xminV0[13]=   -0.5;     // (pt-ptGeant)/ptGeant rec Pos
  xmaxV0[13]=   0.5;


  xminV0[14]=-2.5;   // eta  rec Pos
  xmaxV0[14]= 1.5;


  xminV0[15]=-2*TMath::Pi();   // phi rec Pos
  xmaxV0[15]= TMath::Pi();


  xminV0[16]=   0;     // Impact parameter rec Pos
  xmaxV0[16]=   100;


  xminV0[17]=   0;     // nsigmas Impact parameter rec Pos
  xmaxV0[17]=   100;



  xminV0[18]=   -1;     // Ncls ITS rec Pos
  xmaxV0[18]=   6;


  xminV0[19]=   -1;     // Ncls TPC rec Pos
  xmaxV0[19]=   180;


  xminV0[20]=   -2;     // Status Single  TPC rec Pos
  xmaxV0[20]=   2;


  xminV0[21]=   -0.5;     // (pt-ptGeant)/ptGeant rec Neg
  xmaxV0[21]=   0.5;


  xminV0[22]=-2.5;   // eta  rec Neg
  xmaxV0[22]= 1.5;


  xminV0[23]=-2*TMath::Pi();   // phi rec Neg
  xmaxV0[23]= TMath::Pi();


  xminV0[24]=   0;     // Impact parameter rec Neg
  xmaxV0[24]=   100;


  xminV0[25]=   0;     // Sigmas Impact parameter rec Neg
  xmaxV0[25]=   100;



  xminV0[26]=   -1;     // Ncls ITS rec Neg
  xmaxV0[26]=   6;


  xminV0[27]=   -1;     // Ncls TPC rec Neg
  xmaxV0[27]=   180;


  xminV0[28]=   -2;     // Status Single  TPC rec Neg
  xmaxV0[28]=   2;

  // ------------------Rec V0 variables



  xminV0[29]=   -0.5;     // (pt-ptGeant)/ptGeant rec V0 Pos
  xmaxV0[29]=   0.5;


  xminV0[30]=-2.5;   // eta  rec V0 Pos
  xmaxV0[30]= 1.5;


  xminV0[31]=-2*TMath::Pi();   // phi rec V0 Pos
  xmaxV0[31]= TMath::Pi();


  xminV0[32]=   -2;     // Status V0 TPC rec Pos
  xmaxV0[32]=   2;




  xminV0[33]=   -0.5;     // 1/sqrt(pt) rec V0 Neg
  xmaxV0[33]=   0.5;


  xminV0[34]=-2.5;   // eta  rec V0 Neg
  xmaxV0[34]= 1.5;


  xminV0[35]=-2*TMath::Pi();   // phi rec V0 Neg
  xmaxV0[35]= TMath::Pi();



  xminV0[36]=   -2;     // Status V0 TPC rec Neg
  xmaxV0[36]=   2;



  TString  axisName[37]={"ptGammaGeant", 
			 "etaGammaGeant",
			 "phiGammaGeant",
			 "rGeant",
			 "zGeant",
			 "ptEPlusGeant",  
			 "etaEPlusGeant",
			 "phiEPlusGeant",
			 "TPCTrackLengthEPlusGeant",
			 "ptEMinusGeant", 
			 "etaEMinusGeant",
			 "phiEMinusGeant",
			 "TPCTrackLengthEMinusGeant",
			 "ptResEPlusRecSingle", 
			 "etaEPlusRecSingle",
			 "phiEPlusRecSingle",
			 "bXYZEPlusRecSingle",
			 "sigbXYZEPlusRecSingle",
			 "NclsITSEPlusRecSingle",
			 "NclsTPCEPlusRecSingle",
			 "statusRecSinglePos",
			 "ptResEMinusRecSingle", 
			 "etaEMinusRecSingle",
			 "phiEMinusRecSingle",
			 "bXYZEMinusRecSingle",
			 "sigbXYZEMinusRecSingle",
			 "NclsITSEMinusRecSingle",
			 "NclsTPCEMinusRecSingle",
			 "statusRecSingleNeg",
			 "ptResEPlusRecV0", 
			 "etaEPlusRecV0",
			 "phiEPlusRecV0",
			 "statusV0SinglePos",
			 "ptResEMinusRecV0", 
			 "etaEMinusRecV0",
			 "phiEMinusRecV0",
                         "statusV0SingleNeg"};


  fSparseV0= new THnSparseF("sparseV0","sparseV0",fDim,binsV0,xminV0,xmaxV0);

  for (Int_t iaxis=0; iaxis<fDim; iaxis++){
   fSparseV0->GetAxis(iaxis)->SetName(axisName[iaxis]);
   fSparseV0->GetAxis(iaxis)->SetTitle(axisName[iaxis]);
  }

  TString  axisNameK0[37]={"ptK0Geant", 
			 "etaK0Geant",
			 "phiK0Geant",
			 "rGeant",
			 "zGeant",
			 "ptPiPlusGeant",  
			 "etaPiPlusGeant",
			 "phiPiPlusGeant",
			 "TPCTrackLengthPiPlusGeant",
			 "ptPiMinusGeant", 
			 "etaPiMinusGeant",
			 "phiPiMinusGeant",
			 "TPCTrackLengthPiMinusGeant",
			 "ptResPiPlusRecSingle", 
			 "etaPiPlusRecSingle",
			 "phiPiPlusRecSingle",
			 "bXYZPiPlusRecSingle",
			 "sigbXYZPiPlusRecSingle",
			 "NclsITSPiPlusRecSingle",
			 "NclsTPCPiPlusRecSingle",
			 "statusRecSinglePos",
			 "ptResPiMinusRecSingle", 
			 "etaPiMinusRecSingle",
			 "phiPiMinusRecSingle",
			 "bXYZPiMinusRecSingle",
			 "sigbXYZPiMinusRecSingle",
			 "NclsITSPiMinusRecSingle",
			 "NclsTPCPiMinusRecSingle",
			 "statusRecSingleNeg",
			 "ptResPiPlusRecV0", 
			 "etaPiPlusRecV0",
			 "phiPiPlusRecV0",
			 "statusRecV0Pos",
			 "ptResPiMinusRecV0", 
			 "etaPiMinusRecV0",
			 "phiPiMinusRecV0",
                         "statusRecV0Neg"};



  fSparseK0= new THnSparseF("sparseK0","sparseK0",fDim,binsV0,xminV0,xmaxV0);
  for (Int_t iaxis=0; iaxis<fDim; iaxis++){
   fSparseK0->GetAxis(iaxis)->SetName(axisNameK0[iaxis]);
   fSparseK0->GetAxis(iaxis)->SetTitle(axisNameK0[iaxis]);
  }

  TString  axisNameL[37]={"ptLGeant", 
			 "etaLGeant",
			 "phiLGeant",
			 "rGeant",
			 "zGeant",
			 "ptPPlusGeant",  
			 "etaPPlusGeant",
			 "phiPPlusGeant",
			 "TPCTrackLengthPPlusGeant",
			 "ptPiMinusGeant", 
			 "etaPiMinusGeant",
			 "phiPiMinusGeant",
			 "TPCTrackLengthPiMinusGeant",
			 "ptResPPlusRecSingle", 
			 "etaPPlusRecSingle",
			 "phiPPlusRecSingle",
			 "bXYZPPlusRecSingle",
			 "sigbXYZPPlusRecSingle",
			 "NclsITSPPlusRecSingle",
			 "NclsTPCPPlusRecSingle",
			 "statusRecSinglePos",
			 "ptResPiMinusRecSingle", 
			 "etaPiMinusRecSingle",
			 "phiPiMinusRecSingle",
			 "bXYZPiMinusRecSingle",
			 "sigbXYZPiMinusRecSingle",
			 "NclsITSPiMinusRecSingle",
			 "NclsTPCPiMinusRecSingle",
                         "statusRecSingleNeg",
			 "ptResPPlusRecV0", 
			 "etaPPlusRecV0",
			 "phiPPlusRecV0",
			 "statusRecV0Pos",
			 "ptResPiMinusRecV0", 
			 "etaPiMinusRecV0",
			 "phiPiMinusRecV0",
                         "statusRecV0Neg"};


  fSparseL= new THnSparseF("sparseL","sparseL",fDim,binsV0,xminV0,xmaxV0);
  for (Int_t iaxis=0; iaxis<fDim; iaxis++){
   fSparseL->GetAxis(iaxis)->SetName(axisNameL[iaxis]);
   fSparseL->GetAxis(iaxis)->SetTitle(axisNameL[iaxis]);
  }

  TString  axisNameAL[37]={"ptALGeant", 
			 "etaALGeant",
			 "phiALGeant",
			 "rGeant",
			 "zGeant",
			 "ptPiPluusGeant",  
			 "etaPiPlusGeant",
			 "phiPiPlusGeant",
			 "TPCTrackLengthPiPlusGeant",
			 "ptAPMinusGeant", 
			 "etaAPMinusGeant",
			 "phiAPMinusGeant",
			 "TPCTrackLengthAPMinusGeant",
			 "ptResPiPlusRecSingle", 
			 "etaPiPlusRecSingle",
			 "phiPiPlusRecSingle",
			 "bXYZPiPlusRecSingle",
			 "sigbXYZPiPlusRecSingle",
			 "NclsITSPiPlusRecSingle",
			 "NclsTPCPiPlusRecSingle",
			 "statusRecSinglePos",
			 "ptResAPMinusRecSingle", 
			 "etaAPMinusRecSingle",
			 "phiAPMinusRecSingle",
			 "bXYZAPMinusRecSingle",
			 "sigbXYZAPMinusRecSingle",
			 "NclsITSAPMinusRecSingle",
			 "NclsTPCAPMinusRecSingle",
                         "statusRecSingleNeg",
			 "ptResPiPlusRecV0", 
			 "etaPiPlusRecV0",
			 "phiPiPlusRecV0",
			 "statusRecV0Pos",
			 "ptResAPMinusRecV0", 
			 "etaAPMinusRecV0",
			 "phiAPMinusRecV0",
                         "statusRecV0Neg"};


  fSparseAL= new THnSparseF("sparseAL","sparseAL",fDim,binsV0,xminV0,xmaxV0);
  for (Int_t iaxis=0; iaxis<fDim; iaxis++){
   fSparseAL->GetAxis(iaxis)->SetName(axisNameAL[iaxis]);
   fSparseAL->GetAxis(iaxis)->SetTitle(axisNameAL[iaxis]);
  }

  // create output container
 
  fOutputContainer = new TList() ;
  fOutputContainer->SetName(GetName()) ;


  fOutputContainer->Add(fSparseV0);
  fOutputContainer->Add(fSparseK0);
  fOutputContainer->Add(fSparseL);
  fOutputContainer->Add(fSparseAL);
  
}

//________________________________________________________________________
void AliAnalysisTaskV0QA::Exec(Option_t *) {

  if (!fESD) {
    cout<< "not a tree"<< endl;
    return;
  }

  nEv++;


  //Get MC data 
  mctruth = (AliMCEventHandler*)((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());

  //  Double_t vertex[3];
  Double_t MaxVertex=150.;
  Double_t MaxEta=1.2;
  Double_t LineCutZRSlope=0.662486;
  Double_t LineCutZValue=7.;
  Int_t elecGIndex=-1;
  Int_t posiGIndex=-1;
  Int_t pPosLIndex=-1;
  Int_t piNegLIndex=-1;

  Int_t apNegALIndex=-1;
  Int_t piPosALIndex=-1;

  nConvGamGeant=-1;
  nDecayK0Geant=-1;
  nDecayLGeant=-1;
  nDecayALGeant=-1;

  for(Int_t i=0; i<gDim;i++){
    gConvGamGeantIndex[i] = -1;
    eNegConvGamGeantIndex[i] = -1;
    ePosConvGamGeantIndex[i] = -1;
    
    eNegConvGamSingleRecIndex[i] = -1;
    ePosConvGamSingleRecIndex[i] = -1;

    eNegConvGamV0RecIndex[i] = -1;
    ePosConvGamV0RecIndex[i] = -1;
    ConvGamV0RecIndexPos[i] = -1;
    ConvGamV0RecIndexNeg[i] = -1;


    K0DecayK0GeantIndex[i] = -1;
    piNegDecayK0GeantIndex[i] = -1;
    piPosDecayK0GeantIndex[i] = -1;
    
    piNegDecayK0SingleRecIndex[i] = -1;
    piPosDecayK0SingleRecIndex[i] = -1;

    piNegDecayK0V0RecIndex[i] = -1;
    piPosDecayK0V0RecIndex[i] = -1;
    DecayK0V0RecIndexPos[i] = -1;
    DecayK0V0RecIndexNeg[i] = -1;


    lDecayLGeantIndex[i] = -1;
    piNegDecayLGeantIndex[i] = -1;
    pPosDecayLGeantIndex[i] = -1;
    
    piNegDecayLSingleRecIndex[i] = -1;
    pPosDecayLSingleRecIndex[i] = -1;

    piNegDecayLV0RecIndex[i] = -1;
    pPosDecayLV0RecIndex[i] = -1;
    DecayLV0RecIndexPos[i] = -1;
    DecayLV0RecIndexNeg[i] = -1;

    // Antilambda
    alDecayALGeantIndex[i] = -1;
    piPosDecayALGeantIndex[i] = -1;
    apNegDecayALGeantIndex[i] = -1;
    
    piPosDecayALSingleRecIndex[i] = -1;
    apNegDecayALSingleRecIndex[i] = -1;

    piPosDecayALV0RecIndex[i] = -1;
    apNegDecayALV0RecIndex[i] = -1;
    DecayALV0RecIndexPos[i] = -1;
    DecayALV0RecIndexNeg[i] = -1;


  }

  Int_t doMC=1;

  AliKFVertex primVtx(*(fESD->GetPrimaryVertex()));
  nTracksPrim=primVtx.GetNContributors();


  if(mctruth && nTracksPrim>0){

   stack = mctruth->MCEvent()->Stack();


   if ( doMC){

    for (Int_t iTracks = 0; iTracks < mctruth->MCEvent()->GetNumberOfTracks(); iTracks++) {
      

     TParticle* particle = stack->Particle(iTracks);



     if (!particle) {
       Printf("ERROR: Could not receive particle %d (mc loop)", iTracks);
       continue;
     }
     
     if(particle->Pt()<0.050) continue;
     if(TMath::Abs(particle->Eta())> 1.2) continue;

     
     if (particle->GetPdgCode()== 22){
       
       
       if(particle->GetMother(0) >-1 && stack->Particle(particle->GetMother(0))->GetPdgCode() == 22){
	 continue; // no photon as mothers!
       }
       
       if(particle->GetMother(0) >= stack->GetNprimary()){
	 continue; // the gamma has a mother, and it is not a primary particle
       }
       
       TParticle* ePos = NULL;
       TParticle* eNeg = NULL;
       elecGIndex=-1;
       posiGIndex=-1;
  
       if(particle->GetNDaughters() >= 2){
	 for(Int_t daughterIndex=particle->GetFirstDaughter();daughterIndex<=particle->GetLastDaughter();daughterIndex++){
	   TParticle *tmpDaughter = stack->Particle(daughterIndex);
	   if(tmpDaughter->GetUniqueID() == 5){
	     if(tmpDaughter->GetPdgCode() == 11){
	       eNeg = tmpDaughter;
	       elecGIndex=daughterIndex;
	     }
	     else if(tmpDaughter->GetPdgCode() == -11){
	       ePos = tmpDaughter;
	       posiGIndex=daughterIndex;
	     }
	   }
	 }
       }
       
       
       if(ePos == NULL || eNeg == NULL){ // means we do not have two daughters from pair production
	 continue;
       }
       
       if(TMath::Abs(ePos->Eta())> MaxEta || TMath::Abs(eNeg->Eta())> MaxEta){
	 continue;
       }	
       
       if(ePos->R()> MaxVertex ){
	 continue; // cuts on distance from collision point
       }
      
      
       if( (TMath::Abs(ePos->Vz()) * LineCutZRSlope - LineCutZValue)  > ePos->R() ){
	 continue;               // line cut to exclude regions where we do not reconstruct
       }		

       
       // Looking at the existance of TPC references

       TParticle* ePosTPC;
       mctruth->MCEvent()->GetParticleAndTR(posiGIndex,ePosTPC,clRefsP);

       AliMCParticle *mcParticlePos = mctruth->MCEvent()->GetTrack(posiGIndex);
       if(!mcParticlePos) continue;

       Int_t counter; 
       Float_t tpcTrackLengthePos = mcParticlePos->GetTPCTrackLength(fESD->GetMagneticField(),0.05,counter,3.0); 

 

       int nPointsP =  clRefsP->GetEntries();

       if (fRefTPC) delete fRefTPC;
       fRefTPC = new TObjArray();

       for(int iPoint=0; iPoint<nPointsP; iPoint++) {
	 AliTrackReference *ref = (AliTrackReference*)clRefsP->At(iPoint);
	 if (ref->DetectorId() == AliTrackReference::kTPC) fRefTPC->Add(new AliTrackReference(*ref));	 
       }

       fRefTPC->Sort();

       for(int i=0; i<fRefTPC->GetEntries(); i++) {
	 AliTrackReference *ref = (AliTrackReference*)(*fRefTPC)[i]; 
	 fLabelsTPC[i] = ref->GetTrack();
       }
       int labelPosRefs=0;
       for(int iPoint=GetTPCReference(iTracks);iPoint<fRefTPC->GetEntries();iPoint++){
	 AliTrackReference* aRef = (AliTrackReference*)(*fRefTPC)[iPoint];
	 if (aRef->GetTrack() != posiGIndex ) break;			
	 ++labelPosRefs;
       }
     



       TParticle* eNegTPC;
       mctruth->MCEvent()->GetParticleAndTR(elecGIndex,eNegTPC,clRefsN);

       AliMCParticle *mcParticleNeg = mctruth->MCEvent()->GetTrack(elecGIndex);
       if(!mcParticleNeg) continue;

       Int_t counterN; 
       Float_t tpcTrackLengtheNeg = mcParticleNeg->GetTPCTrackLength(fESD->GetMagneticField(),0.05,counterN,3.0); 
       int nPointsN =  clRefsN->GetEntries();

       if (fRefTPC) delete fRefTPC;
       fRefTPC = new TObjArray();

       for(int iPoint=0; iPoint<nPointsN; iPoint++) {
	 AliTrackReference *ref = (AliTrackReference*)clRefsN->At(iPoint);
	 if (ref->DetectorId() == AliTrackReference::kTPC) fRefTPC->Add(new AliTrackReference(*ref));	 
       }

       fRefTPC->Sort();

       for(int i=0; i<fRefTPC->GetEntries(); i++) {
	 AliTrackReference *ref = (AliTrackReference*)(*fRefTPC)[i]; 
	 fLabelsTPC[i] = ref->GetTrack();
       }
       int labelNegRefs=0;
       for(int iPoint=GetTPCReference(iTracks);iPoint<fRefTPC->GetEntries();iPoint++){
	 AliTrackReference* aRef = (AliTrackReference*)(*fRefTPC)[iPoint];
	 if (aRef->GetTrack() != elecGIndex ) break;			
	 ++labelNegRefs;
       }


       
       if ( labelNegRefs==0 || labelPosRefs==0) continue; // if e+/e- do not have a TPC ref continue;
       ////////////////////////////////////////////////////////////////////
       
       
       nConvGamGeant++;
       gConvGamGeantIndex[nConvGamGeant]=iTracks;
       eNegConvGamGeantIndex[nConvGamGeant] = elecGIndex;
       ePosConvGamGeantIndex[nConvGamGeant] = posiGIndex;
       
       eNegConvGamGeantLength[nConvGamGeant] = tpcTrackLengtheNeg;
       ePosConvGamGeantLength[nConvGamGeant] = tpcTrackLengthePos;

     }

     
     TParticle* piPos = NULL;
     TParticle* piNeg = NULL;
     piPosK0Index=-1;
     piNegK0Index=-1;

     if (particle->GetPdgCode()== 310){          // k0short
       if(particle->GetNDaughters() == 2){
	 for(Int_t daughterIndex=particle->GetFirstDaughter();daughterIndex<=particle->GetLastDaughter();daughterIndex++){
	   TParticle *tmpDaughter = stack->Particle(daughterIndex);
	   if(tmpDaughter->GetPdgCode() == 211){
	     piPos= tmpDaughter;
	     piPosK0Index=daughterIndex;
	   }
	   else if(tmpDaughter->GetPdgCode() == -211){
	     piNeg = tmpDaughter;
	     piNegK0Index=daughterIndex;
	   }
	 }
       }

       if(piPos == NULL || piNeg == NULL){ // means we do not have two daughters from K0short decay
	 continue;
       }
       
       if(TMath::Abs(piPos->Eta())> MaxEta || TMath::Abs(piNeg->Eta())> MaxEta){
	 continue;
       }	
       
       if(piPos->R()> MaxVertex ){
	 continue; // cuts on distance from collision point
       }
      
      
       if( (TMath::Abs(piPos->Vz()) * LineCutZRSlope - LineCutZValue)  > piPos->R() ){
	 continue;               // line cut to exclude regions where we do not reconstruct
       }		

      // Looking at the existance of TPC references

       TParticle* ePosTPC;
       mctruth->MCEvent()->GetParticleAndTR(piPosK0Index,ePosTPC,clRefsP);

       AliMCParticle *mcParticlePos = mctruth->MCEvent()->GetTrack(piPosK0Index);
       if(!mcParticlePos) continue;

       Int_t counter; 
       Float_t tpcTrackLengthePos = mcParticlePos->GetTPCTrackLength(fESD->GetMagneticField(),0.05,counter,3.0); 
 

       int nPointsP =  clRefsP->GetEntries();
       if (fRefTPC) delete fRefTPC;
       fRefTPC = new TObjArray();

       for(int iPoint=0; iPoint<nPointsP; iPoint++) {
	 AliTrackReference *ref = (AliTrackReference*)clRefsP->At(iPoint);
	 if (ref->DetectorId() == AliTrackReference::kTPC) fRefTPC->Add(new AliTrackReference(*ref));	 
       }

       fRefTPC->Sort();

       for(int i=0; i<fRefTPC->GetEntries(); i++) {
	 AliTrackReference *ref = (AliTrackReference*)(*fRefTPC)[i]; 
	 fLabelsTPC[i] = ref->GetTrack();
       }
       int labelPosRefs=0;
       for(int iPoint=GetTPCReference(iTracks);iPoint<fRefTPC->GetEntries();iPoint++){
	 AliTrackReference* aRef = (AliTrackReference*)(*fRefTPC)[iPoint];
	 if (aRef->GetTrack() != piPosK0Index ) break;			
	 ++labelPosRefs;
       }
     



       TParticle* eNegTPC;
       mctruth->MCEvent()->GetParticleAndTR(piNegK0Index,eNegTPC,clRefsN);

       AliMCParticle *mcParticleNeg = mctruth->MCEvent()->GetTrack(piNegK0Index);
       if(!mcParticleNeg) continue;

       Int_t counterN; 
       Float_t tpcTrackLengtheNeg = mcParticleNeg->GetTPCTrackLength(fESD->GetMagneticField(),0.05,counterN,3.0); 

       int nPointsN =  clRefsN->GetEntries();
       if (fRefTPC) delete fRefTPC;
       fRefTPC = new TObjArray();

       for(int iPoint=0; iPoint<nPointsN; iPoint++) {
	 AliTrackReference *ref = (AliTrackReference*)clRefsN->At(iPoint);
	 if (ref->DetectorId() == AliTrackReference::kTPC) fRefTPC->Add(new AliTrackReference(*ref));	 
       }

       fRefTPC->Sort();

       for(int i=0; i<fRefTPC->GetEntries(); i++) {
	 AliTrackReference *ref = (AliTrackReference*)(*fRefTPC)[i]; 
	 fLabelsTPC[i] = ref->GetTrack();
       }
       int labelNegRefs=0;
       for(int iPoint=GetTPCReference(iTracks);iPoint<fRefTPC->GetEntries();iPoint++){
	 AliTrackReference* aRef = (AliTrackReference*)(*fRefTPC)[iPoint];
	 if (aRef->GetTrack() != piNegK0Index ) break;			
	 ++labelNegRefs;
       }
       
       if ( labelNegRefs==0 || labelPosRefs==0) continue; // if pi+/pi- do not have a TPC ref continue;
       ////////////////////////////////////////////////////////////////////
       
       nDecayK0Geant++;

       K0DecayK0GeantIndex[nDecayK0Geant]=iTracks;
       piNegDecayK0GeantIndex[nDecayK0Geant]=piNegK0Index;
       piPosDecayK0GeantIndex[nDecayK0Geant]=piPosK0Index;
       piNegDecayK0GeantLength[nDecayK0Geant]=tpcTrackLengtheNeg;
       piPosDecayK0GeantLength[nDecayK0Geant]=tpcTrackLengthePos;
      
     }    


     TParticle* pPos = NULL;
     TParticle* piNegL = NULL;
     pPosLIndex=-1;
     piNegLIndex=-1;

 
     if (particle->GetPdgCode()== 3122){        //lambda
      
       if(particle->GetNDaughters() == 2){
	 for(Int_t daughterIndex=particle->GetFirstDaughter();daughterIndex<=particle->GetLastDaughter();daughterIndex++){
	   TParticle *tmpDaughter = stack->Particle(daughterIndex);
	   if(tmpDaughter->GetPdgCode() == 2212){
	     pPos= tmpDaughter;
	     pPosLIndex=daughterIndex;
	   }
	   else if(tmpDaughter->GetPdgCode() == -211){
	     piNegL = tmpDaughter;
	     piNegLIndex=daughterIndex;
	   }
	 }
       }

       if(pPos == NULL || piNegL == NULL){ // means we do not have two daughters from lambda decay
	 continue;
       }

       if(TMath::Abs(pPos->Eta())> MaxEta || TMath::Abs(piNegL->Eta())> MaxEta){
	 continue;
       }	
       
       if(pPos->R()> MaxVertex ){
	 continue; // cuts on distance from collision point
       }
      
      
       if( (TMath::Abs(pPos->Vz()) * LineCutZRSlope - LineCutZValue)  > pPos->R() ){
	 continue;               // line cut to exclude regions where we do not reconstruct
       }


     // Looking at the existance of TPC references

       TParticle* ePosTPC;
       mctruth->MCEvent()->GetParticleAndTR(pPosLIndex,ePosTPC,clRefsP);

       AliMCParticle *mcParticlePos = mctruth->MCEvent()->GetTrack(pPosLIndex);
       if(!mcParticlePos) continue;

       Int_t counter; 
       Float_t tpcTrackLengthePos = mcParticlePos->GetTPCTrackLength(fESD->GetMagneticField(),0.05,counter,3.0); 
 

       int nPointsP =  clRefsP->GetEntries();
       if (fRefTPC) delete fRefTPC;
       fRefTPC = new TObjArray();

       for(int iPoint=0; iPoint<nPointsP; iPoint++) {
	 AliTrackReference *ref = (AliTrackReference*)clRefsP->At(iPoint);
	 if (ref->DetectorId() == AliTrackReference::kTPC) fRefTPC->Add(new AliTrackReference(*ref));	 
       }

       fRefTPC->Sort();

       for(int i=0; i<fRefTPC->GetEntries(); i++) {
	 AliTrackReference *ref = (AliTrackReference*)(*fRefTPC)[i]; 
	 fLabelsTPC[i] = ref->GetTrack();
       }
       int labelPosRefs=0;
       for(int iPoint=GetTPCReference(iTracks);iPoint<fRefTPC->GetEntries();iPoint++){
	 AliTrackReference* aRef = (AliTrackReference*)(*fRefTPC)[iPoint];
	 if (aRef->GetTrack() != pPosLIndex ) break;			
	 ++labelPosRefs;
       }
     



       TParticle* eNegTPC;
       mctruth->MCEvent()->GetParticleAndTR(piNegLIndex,eNegTPC,clRefsN);

       AliMCParticle *mcParticleNeg = mctruth->MCEvent()->GetTrack(piNegLIndex);
       if(!mcParticleNeg) continue;

       Int_t counterN; 
       Float_t tpcTrackLengtheNeg = mcParticleNeg->GetTPCTrackLength(fESD->GetMagneticField(),0.05,counterN,3.0); 

       int nPointsN =  clRefsN->GetEntries();
       if (fRefTPC) delete fRefTPC;
       fRefTPC = new TObjArray();

       for(int iPoint=0; iPoint<nPointsN; iPoint++) {
	 AliTrackReference *ref = (AliTrackReference*)clRefsN->At(iPoint);
	 if (ref->DetectorId() == AliTrackReference::kTPC) fRefTPC->Add(new AliTrackReference(*ref));	 
       }

       fRefTPC->Sort();

       for(int i=0; i<fRefTPC->GetEntries(); i++) {
	 AliTrackReference *ref = (AliTrackReference*)(*fRefTPC)[i]; 
	 fLabelsTPC[i] = ref->GetTrack();
       }
       int labelNegRefs=0;
       for(int iPoint=GetTPCReference(iTracks);iPoint<fRefTPC->GetEntries();iPoint++){
	 AliTrackReference* aRef = (AliTrackReference*)(*fRefTPC)[iPoint];
	 if (aRef->GetTrack() != piNegLIndex ) break;			
	 ++labelNegRefs;
       }
       
       if ( labelNegRefs==0 || labelPosRefs==0) continue; // if proton/pi- do not have a TPC ref continue;
       ////////////////////////////////////////////////////////////////////
     
       nDecayLGeant++;

       lDecayLGeantIndex[nDecayLGeant]=iTracks;
     
       piNegDecayLGeantIndex[nDecayLGeant]=piNegLIndex;
       pPosDecayLGeantIndex[nDecayLGeant]=pPosLIndex;

       piNegDecayLGeantLength[nDecayLGeant]=tpcTrackLengtheNeg;
       pPosDecayLGeantLength[nDecayLGeant]=tpcTrackLengthePos;

     
     }

     /////////////////////////////////////////////////////

     // AntiLambda    
     TParticle* apNeg = NULL;
     TParticle* piPosAL = NULL;
     apNegALIndex=-1;
     piPosALIndex=-1;

 
     if (particle->GetPdgCode()== -3122){        //antilambda

       if(particle->GetNDaughters() == 2){
	 for(Int_t daughterIndex=particle->GetFirstDaughter();daughterIndex<=particle->GetLastDaughter();daughterIndex++){
	   TParticle *tmpDaughter = stack->Particle(daughterIndex);
	   if(tmpDaughter->GetPdgCode() == -2212){
	     apNeg= tmpDaughter;
	     apNegALIndex=daughterIndex;
	   }
	   else if(tmpDaughter->GetPdgCode() == 211){
	     piPosAL = tmpDaughter;
	     piPosALIndex=daughterIndex;
	   }
	 }
       }

       if(apNeg == NULL || piPosAL == NULL){ // means we do not have two daughters from antilambda decay
	 continue;
       }

       if(TMath::Abs(apNeg->Eta())> MaxEta || TMath::Abs(piPosAL->Eta())> MaxEta){
	 continue;
       }	
       
       if(apNeg->R()> MaxVertex ){
	 continue; // cuts on distance from collision point
       }
      
      
       if( (TMath::Abs(apNeg->Vz()) * LineCutZRSlope - LineCutZValue)  > apNeg->R() ){
	 continue;               // line cut to exclude regions where we do not reconstruct
       }


     // Looking at the existance of TPC references

       TParticle* ePosTPC;
       mctruth->MCEvent()->GetParticleAndTR(piPosALIndex,ePosTPC,clRefsP);

       AliMCParticle *mcParticlePos = mctruth->MCEvent()->GetTrack(piPosALIndex);
       if(!mcParticlePos) continue;

       Int_t counter; 
       Float_t tpcTrackLengthePos = mcParticlePos->GetTPCTrackLength(fESD->GetMagneticField(),0.05,counter,3.0); 

       int nPointsP =  clRefsP->GetEntries();
       if (fRefTPC) delete fRefTPC;
       fRefTPC = new TObjArray();

       for(int iPoint=0; iPoint<nPointsP; iPoint++) {
	 AliTrackReference *ref = (AliTrackReference*)clRefsP->At(iPoint);
	 if (ref->DetectorId() == AliTrackReference::kTPC) fRefTPC->Add(new AliTrackReference(*ref));	 
       }

       fRefTPC->Sort();

       for(int i=0; i<fRefTPC->GetEntries(); i++) {
	 AliTrackReference *ref = (AliTrackReference*)(*fRefTPC)[i]; 
	 fLabelsTPC[i] = ref->GetTrack();
       }
       int labelPosRefs=0;
       for(int iPoint=GetTPCReference(iTracks);iPoint<fRefTPC->GetEntries();iPoint++){
	 AliTrackReference* aRef = (AliTrackReference*)(*fRefTPC)[iPoint];
	 if (aRef->GetTrack() != piPosALIndex ) break;			
	 ++labelPosRefs;
       }
     

       TParticle* eNegTPC;
       mctruth->MCEvent()->GetParticleAndTR(apNegALIndex,eNegTPC,clRefsN);

       AliMCParticle *mcParticleNeg = mctruth->MCEvent()->GetTrack(apNegALIndex);
       if(!mcParticleNeg) continue;

       Int_t counterN; 
       Float_t tpcTrackLengtheNeg = mcParticleNeg->GetTPCTrackLength(fESD->GetMagneticField(),0.05,counterN,3.0); 

       int nPointsN =  clRefsN->GetEntries();
       if (fRefTPC) delete fRefTPC;
       fRefTPC = new TObjArray();

       for(int iPoint=0; iPoint<nPointsN; iPoint++) {
	 AliTrackReference *ref = (AliTrackReference*)clRefsN->At(iPoint);
	 if (ref->DetectorId() == AliTrackReference::kTPC) fRefTPC->Add(new AliTrackReference(*ref));	 
       }

       fRefTPC->Sort();

       for(int i=0; i<fRefTPC->GetEntries(); i++) {
	 AliTrackReference *ref = (AliTrackReference*)(*fRefTPC)[i]; 
	 fLabelsTPC[i] = ref->GetTrack();
       }
       int labelNegRefs=0;
       for(int iPoint=GetTPCReference(iTracks);iPoint<fRefTPC->GetEntries();iPoint++){
	 AliTrackReference* aRef = (AliTrackReference*)(*fRefTPC)[iPoint];
	 if (aRef->GetTrack() != apNegALIndex ) break;			
	 ++labelNegRefs;
       }


       
       if ( labelNegRefs==0 || labelPosRefs==0) continue; // if proton/pi- do not have a TPC ref continue;
       ////////////////////////////////////////////////////////////////////
     
       nDecayALGeant++;
       alDecayALGeantIndex[nDecayALGeant]=iTracks;
     
       piPosDecayALGeantIndex[nDecayALGeant]=piPosALIndex;
       apNegDecayALGeantIndex[nDecayALGeant]=apNegALIndex;

       piPosDecayALGeantLength[nDecayALGeant]=tpcTrackLengthePos;
       apNegDecayALGeantLength[nDecayALGeant]=tpcTrackLengtheNeg;

     
     }  // AntiLambda    

    } //track loop 
   }    

  }


  AliKFParticle::SetField(fESD->GetMagneticField());

  const AliESDVertex *pvertex = fESD->GetPrimaryVertex();
  Double_t xyzVtx[3];
  pvertex->GetXYZ(xyzVtx);

  if(nTracksPrim>0) {

    InspectListOfChargedParticles();
    InspectListOfV0s();


    if(nConvGamGeant>-1){
      FillHnSparseGamma();
    }

    if(nDecayK0Geant>-1){
      FillHnSparseK0();
    }

    if(nDecayLGeant>-1){
      FillHnSparseL();
    }

    if(nDecayALGeant>-1){
      FillHnSparseAL();
    }

  }

 
  PostData(0,fOutputContainer );
  

}

void AliAnalysisTaskV0QA::Terminate(Option_t *) {
  // Draw some histogram at the end.

}


Int_t AliAnalysisTaskV0QA::GetTPCReference(Int_t label) {

  int start = TMath::BinarySearch(fRefTPC->GetEntries(), fLabelsTPC, label);

  while (start >= 0) {
    AliTrackReference *ref = (AliTrackReference*)(*fRefTPC)[start];
    if (ref->GetTrack() != label) return start+1;
    start--;
  }

  return 0;
}





void AliAnalysisTaskV0QA::InspectListOfChargedParticles(){


  for(Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++){

    AliESDtrack* curTrack = fESD->GetTrack(iTracks);

    if(!curTrack){
      continue;
    }


//     if( !(curTrack->GetStatus() & AliESDtrack::kTPCrefit)){
//       continue;
//     }


    Int_t labelMC = TMath::Abs(curTrack->GetLabel());

    if ( labelMC > stack->GetNtrack() ) continue;


    TParticle* curParticle = stack->Particle(labelMC);
    if(curParticle->GetMother(0)==-1){
      continue;
    }

     
    if(TMath::Abs(curParticle->GetPdgCode()) == 11){ // e+/e-
      
      if( stack->Particle(curParticle->GetMother(0))->GetPdgCode()==22 ){  // e+/e- from gamma
	if( curParticle->GetUniqueID()!=5 ){ // e+/e- from gamma conversion
	  continue;
	}

	for(Int_t iGamConv=0;iGamConv<nConvGamGeant+1;iGamConv++ ){
	  if(curTrack->GetSign()>0){
	    if (labelMC== ePosConvGamGeantIndex[iGamConv]){
	      ePosConvGamSingleRecIndex[iGamConv]=iTracks;
	    }
	  }else{
	    if (labelMC== eNegConvGamGeantIndex[iGamConv]){
	      eNegConvGamSingleRecIndex[iGamConv]=iTracks;
	    }
	  }
	} // loop over geant converted gammas

      }
    } // condition to select reconstructed electrons
  


    if(TMath::Abs(curParticle->GetPdgCode()) == 211 || TMath::Abs(curParticle->GetPdgCode())==2212 ){ // pi+/pi-
      
      if( stack->Particle(curParticle->GetMother(0))->GetPdgCode()==310 || 
	  stack->Particle(curParticle->GetMother(0))->GetPdgCode()==3122 || 
	  stack->Particle(curParticle->GetMother(0))->GetPdgCode()==-3122 ){  // pi+/proton/pi- from K0/Lambda

	for(Int_t iK0Dec=0;iK0Dec<nDecayK0Geant+1;iK0Dec++ ){
	  if(curTrack->GetSign()>0){
	    if (labelMC== piPosDecayK0GeantIndex[iK0Dec]){
	      piPosDecayK0SingleRecIndex[iK0Dec]=iTracks;
	    }
	  }else{
	    if (labelMC== piNegDecayK0GeantIndex[iK0Dec]){
	      piNegDecayK0SingleRecIndex[iK0Dec]=iTracks;
	    }
	  }
	} // loop over geant decay K0

	for(Int_t iLDec=0;iLDec<nDecayLGeant+1;iLDec++ ){
	  if(curTrack->GetSign()>0){
	    if (labelMC== pPosDecayLGeantIndex[iLDec]){
	      pPosDecayLSingleRecIndex[iLDec]=iTracks;
	    }
	  }else{
	    if (labelMC== piNegDecayLGeantIndex[iLDec]){
	      piNegDecayLSingleRecIndex[iLDec]=iTracks;
	    }
	  }
	} // loop over geant decay Lambda

	for(Int_t iALDec=0;iALDec<nDecayALGeant+1;iALDec++ ){
	  if(curTrack->GetSign()<0){
	    if (labelMC== apNegDecayALGeantIndex[iALDec]){
	      apNegDecayALSingleRecIndex[iALDec]=iTracks;
	    }
	  }else{
	    if (labelMC== piPosDecayALGeantIndex[iALDec]){
	      piPosDecayALSingleRecIndex[iALDec]=iTracks;
	    }
	  }
	} // loop over geant decay antiLambda
      }
    } // condition to select reconstructed electrons
   } // all reconstructed track


}

void AliAnalysisTaskV0QA::InspectListOfV0s(){

  AliESDtrack* trackPos= NULL;
  AliESDtrack* trackNeg= NULL;
  Int_t grandMotherPos=-1;
  Int_t grandMotherNeg=-1;
  Int_t motherPos=-1;
  Int_t motherNeg=-1;
  Int_t pIndex=-1;
  Int_t nIndex=-1;

  for(Int_t iV0MI = 0; iV0MI < fESD->GetNumberOfV0s(); iV0MI++) {

    AliESDv0 * fV0MIs = fESD->GetV0(iV0MI);
    

    if ( !fV0MIs->GetOnFlyStatus() ){
	continue;
    }
    if(nTracksPrim<=0) {
      continue;
    }


    AliESDtrack* trackPosTest = fESD->GetTrack(fV0MIs->GetPindex());
    AliESDtrack* trackNegTest = fESD->GetTrack(fV0MIs->GetNindex());


    if ( trackPosTest->GetSign() == trackNegTest->GetSign()){
     continue;
    }

    //To avoid ghosts

//     if( !(trackPosTest->GetStatus() & AliESDtrack::kTPCrefit)){
//       continue;
//     }

//     if( !(trackNegTest->GetStatus() & AliESDtrack::kTPCrefit)){
//       continue;
//     }

    if( trackPosTest->GetSign() ==1){
      trackPos =fESD->GetTrack(fV0MIs->GetPindex());
      trackNeg =fESD->GetTrack(fV0MIs->GetNindex());
      pIndex=fV0MIs->GetPindex();
      nIndex=fV0MIs->GetNindex();
    }

    if( trackPosTest->GetSign() ==-1){
      trackPos =fESD->GetTrack(fV0MIs->GetNindex());
      trackNeg =fESD->GetTrack(fV0MIs->GetPindex());
      pIndex=fV0MIs->GetNindex();
      nIndex=fV0MIs->GetPindex();

    }

    Int_t labelNeg=TMath::Abs(trackNeg->GetLabel());
    if(labelNeg > stack->GetNtrack() ) continue;
    TParticle * particleNeg= stack->Particle(labelNeg);

    Int_t labelPos=TMath::Abs(trackPos->GetLabel());
    if(labelPos > stack->GetNtrack() ) continue;
    TParticle * particlePos= stack->Particle(labelPos);


    if(particlePos->GetMother(0)>-1){
      grandMotherPos=stack->Particle(particlePos->GetMother(0))->GetMother(0);
      motherPos=particlePos->GetMother(0);
    }
    
    if(particleNeg->GetMother(0)>-1){
      grandMotherNeg=stack->Particle(particleNeg->GetMother(0))->GetMother(0);
      motherNeg=particleNeg->GetMother(0);
    }

    if(motherPos == motherNeg &&  motherPos!=-1 ){
      if( particlePos->GetPdgCode() ==-11  &&   particleNeg->GetPdgCode()==11 ){
	for(Int_t iGamConv=0;iGamConv<nConvGamGeant+1;iGamConv++ ){
	  if (labelPos== ePosConvGamGeantIndex[iGamConv]){
	    ePosConvGamV0RecIndex[iGamConv]=pIndex;
	    ConvGamV0RecIndexPos[iGamConv]=iV0MI;
	  }
	  if (labelNeg== eNegConvGamGeantIndex[iGamConv]){
	    eNegConvGamV0RecIndex[iGamConv]=nIndex;
	    ConvGamV0RecIndexNeg[iGamConv]=iV0MI;
	  }

	} // loop over geant converted gammas
      }

      if( particlePos->GetPdgCode()==211  &&   particleNeg->GetPdgCode()==-211 ){
	for(Int_t iK0Dec=0;iK0Dec<nDecayK0Geant+1;iK0Dec++ ){
	  if (labelPos== piPosDecayK0GeantIndex[iK0Dec]){
	    piPosDecayK0V0RecIndex[iK0Dec]=pIndex;
	    DecayK0V0RecIndexPos[iK0Dec]=iV0MI;
	  }
	  if (labelNeg== piNegDecayK0GeantIndex[iK0Dec]){
	    piNegDecayK0V0RecIndex[iK0Dec]=nIndex;
	    DecayK0V0RecIndexNeg[iK0Dec]=iV0MI;
	  }

	} // loop over geant K0
      }

      if( particlePos->GetPdgCode()==2212  && particleNeg->GetPdgCode()==-211 ){
	for(Int_t iLDec=0;iLDec<nDecayLGeant+1;iLDec++ ){
	  if (labelPos== pPosDecayLGeantIndex[iLDec]){
	    pPosDecayLV0RecIndex[iLDec]=pIndex;
	    DecayLV0RecIndexPos[iLDec]=iV0MI;
	  }
	  if (labelNeg== piNegDecayLGeantIndex[iLDec]){
	    piNegDecayLV0RecIndex[iLDec]=nIndex;
	    DecayLV0RecIndexNeg[iLDec]=iV0MI;
	  }

	} // loop over geant Lambda
      }

      if( particleNeg->GetPdgCode()==-2212  && particlePos->GetPdgCode()==211 ){
	for(Int_t iALDec=0;iALDec<nDecayALGeant+1;iALDec++ ){
	  if (labelNeg== apNegDecayALGeantIndex[iALDec]){
	    apNegDecayALV0RecIndex[iALDec]=nIndex;
	    DecayALV0RecIndexNeg[iALDec]=iV0MI;
	  }
	  if (labelPos== piPosDecayALGeantIndex[iALDec]){
	    piPosDecayALV0RecIndex[iALDec]=pIndex;
	    DecayALV0RecIndexPos[iALDec]=iV0MI;
	  }

	} // loop over geant antiLambda
      }


    }
    
  }
  for(Int_t iGamConv=0;iGamConv<nConvGamGeant+1;iGamConv++ ){
    if ( ConvGamV0RecIndexNeg[iGamConv]!=  ConvGamV0RecIndexPos[iGamConv]){
      ePosConvGamV0RecIndex[iGamConv]=-1;
      eNegConvGamV0RecIndex[iGamConv]=-1; 
      ConvGamV0RecIndexNeg[iGamConv]=-1;
      ConvGamV0RecIndexPos[iGamConv]=-1;

    }
  }

  for(Int_t iLDec=0;iLDec<nDecayLGeant+1;iLDec++ ){
    if(DecayLV0RecIndexPos[iLDec] !=  DecayLV0RecIndexNeg[iLDec]){
      piNegDecayLV0RecIndex[iLDec]=-1;
      pPosDecayLV0RecIndex[iLDec]=-1;
      DecayLV0RecIndexNeg[iLDec]=-1;
      DecayLV0RecIndexPos[iLDec]=-1;
    }
  }

  for(Int_t iALDec=0;iALDec<nDecayALGeant+1;iALDec++ ){
    if(DecayALV0RecIndexPos[iALDec] !=  DecayALV0RecIndexNeg[iALDec]){
      piPosDecayALV0RecIndex[iALDec]=-1;
      apNegDecayALV0RecIndex[iALDec]=-1;
      DecayALV0RecIndexNeg[iALDec]=-1;
      DecayALV0RecIndexPos[iALDec]=-1;
    }
  }

  for(Int_t iK0Dec=0;iK0Dec<nDecayK0Geant+1;iK0Dec++ ){
    if(DecayK0V0RecIndexPos[iK0Dec] !=  DecayK0V0RecIndexNeg[iK0Dec]){
      piNegDecayK0V0RecIndex[iK0Dec]=-1;
      piPosDecayK0V0RecIndex[iK0Dec]=-1;
      DecayK0V0RecIndexNeg[iK0Dec]=-1;
      DecayK0V0RecIndexPos[iK0Dec]=-1;
    }
  }
  

}
void AliAnalysisTaskV0QA::FillHnSparseGamma()
{
  Double_t massE=0.00051099892;
  Double_t ppSgl[3];
  Double_t pmSgl[3];
  Float_t bPosSgl[2];
  Float_t bNegSgl[2];
  Float_t bPosCov[3];
  Float_t bNegCov[3];
  
  Double_t ppV0[3];
  Double_t pmV0[3];
  Double_t xrG[3];
  
  TLorentzVector posSglTrack;
  TLorentzVector negSglTrack;
  Double_t posPt,posEta,posPhi;
  Double_t negPt,negEta,negPhi;

  TLorentzVector posV0Track;
  TLorentzVector negV0Track;
  Double_t posV0Pt,posV0Eta,posV0Phi;
  Double_t negV0Pt,negV0Eta,negV0Phi;
  
  Float_t nClsITSPos=-1;
  Float_t nClsITSNeg=-1;

  Float_t nClsTPCPos=-1;
  Float_t nClsTPCNeg=-1;

  Int_t statusSingPos=-1;
  Int_t statusSingNeg=-1;

  Int_t statusV0Pos=-1;
  Int_t statusV0Neg=-1;


  for(Int_t i=0;i<nConvGamGeant+1;i++){
    TParticle* gamPart = stack->Particle(gConvGamGeantIndex[i]);
    TParticle* ePosPart = stack->Particle(ePosConvGamGeantIndex[i]);
    TParticle* eNegPart = stack->Particle(eNegConvGamGeantIndex[i]);
    if (ePosConvGamSingleRecIndex[i]!=-1){
      AliESDtrack * ePosSglTrack = fESD->GetTrack(ePosConvGamSingleRecIndex[i]);
      ePosSglTrack->GetPxPyPz(ppSgl); 
      posSglTrack.SetXYZM(ppSgl[0],ppSgl[1],ppSgl[2],massE);
      posPt  = posSglTrack.Pt();
      posEta = posSglTrack.Eta();
      posPhi = posSglTrack.Phi();
      ePosSglTrack->GetImpactParameters(bPosSgl,bPosCov);
      nClsITSPos=ePosSglTrack->GetNcls(0);
      nClsTPCPos=ePosSglTrack->GetNcls(1);
      statusSingPos=1;
    }else{
      posPt  = 1000000;
      posEta = -2.;
      posPhi = -2*TMath::Pi();
      bPosSgl[0]=-100.;
      bPosSgl[1]=-100.;
      bPosCov[0]=-100;
      bPosCov[2]=-100;
      nClsITSPos=-1;
      nClsTPCPos=-1;
      statusSingPos=-1;
    }
    
    if (eNegConvGamSingleRecIndex[i]!=-1){
      AliESDtrack * eNegSglTrack = fESD->GetTrack(eNegConvGamSingleRecIndex[i]);
      eNegSglTrack->GetPxPyPz(pmSgl); 
      negSglTrack.SetXYZM(pmSgl[0],pmSgl[1],pmSgl[2],massE);
      negPt  = negSglTrack.Pt();
      negEta = negSglTrack.Eta();
      negPhi = negSglTrack.Phi();
      eNegSglTrack->GetImpactParameters(bNegSgl,bNegCov);
      nClsITSNeg=eNegSglTrack->GetNcls(0);
      nClsTPCNeg=eNegSglTrack->GetNcls(1);
      statusSingNeg=1;
    }else{
      negPt  = 1000000;
      negEta = -2.;
      negPhi = -2*TMath::Pi();
      bNegSgl[0]=-100.;
      bNegSgl[1]=-100.;
      bNegCov[0]=-100;
      bNegCov[2]=-100;
      nClsITSNeg=-1;
      nClsTPCNeg=-1;
      statusSingNeg=-1;
    }

    posV0Pt  = 1000000;
    posV0Eta = -2.;
    posV0Phi = -2*TMath::Pi();
    negV0Pt  = 1000000;
    negV0Eta = -2.;
    negV0Phi = -2*TMath::Pi();
    
    if(ConvGamV0RecIndexPos[i]!=-1){
      AliESDv0 * fV0MIs = fESD->GetV0(ConvGamV0RecIndexPos[i]);
      AliESDtrack* trackPosTest = fESD->GetTrack(fV0MIs->GetPindex());
      AliESDtrack* trackNegTest = fESD->GetTrack(fV0MIs->GetNindex());

      if (ePosConvGamV0RecIndex[i]!=-1 ){
	//	AliESDtrack * ePosV0Track = fESD->GetTrack(ePosConvGamV0RecIndex[i]);
	if ( trackPosTest->GetSign()==1 ) {
	  fV0MIs->GetPPxPyPz(ppV0[0],ppV0[1],ppV0[2]);
	}else{
	  fV0MIs->GetNPxPyPz(ppV0[0],ppV0[1],ppV0[2]);
	}
	posV0Track.SetXYZM(ppV0[0],ppV0[1],ppV0[2],massE);

	posV0Pt  = posV0Track.Pt();
	posV0Eta = posV0Track.Eta();
	posV0Phi = posV0Track.Phi();
	statusV0Pos=1;
      }else{
	posV0Pt  = 1000000;
	posV0Eta = -2.;
	posV0Phi = -2*TMath::Pi();
	statusV0Pos=-1;
      }
      
      if (eNegConvGamV0RecIndex[i]!=-1 ){
	//	AliESDtrack * eNegV0Track = fESD->GetTrack(eNegConvGamV0RecIndex[i]);
	if ( trackNegTest->GetSign()==-1 ) {
	  fV0MIs->GetNPxPyPz(pmV0[0],pmV0[1],pmV0[2]);
	}else{
	  fV0MIs->GetPPxPyPz(pmV0[0],pmV0[1],pmV0[2]);
	}
	negV0Track.SetXYZM(pmV0[0],pmV0[1],pmV0[2],massE);
	
	negV0Pt  = negV0Track.Pt();
	negV0Eta = negV0Track.Eta();
	negV0Phi = negV0Track.Phi();
	statusV0Neg=1;
      }else{
	negV0Pt  = 1000000;
	negV0Eta = -2.;
	negV0Phi = -2*TMath::Pi();
	statusV0Neg=-1;
      }
    }
    
    xrG[0] = ePosPart->Vx();
    xrG[1] = ePosPart->Vy();
    xrG[2] = ePosPart->Vz();
    
    //--------- Geant variables ----------------------
    fValueV0[0] = 1./TMath::Sqrt(gamPart->Pt());
    fValueV0[1] = gamPart->Eta();

    Double_t tmpGPhi=gamPart->Phi();
    if( gamPart->Phi()>TMath::Pi()){
      tmpGPhi=gamPart->Phi()-2*TMath::Pi();
    }
    fValueV0[2] = tmpGPhi;

    fValueV0[3] = TMath::Sqrt(xrG[0]*xrG[0]+xrG[1]*xrG[1]);
    fValueV0[4] = xrG[2];
 

    fValueV0[5] = 1./TMath::Sqrt(ePosPart->Pt());
    fValueV0[6] = ePosPart->Eta();

    Double_t tmpPPhi=ePosPart->Phi();
    if( ePosPart->Phi()>TMath::Pi()){
      tmpPPhi = ePosPart->Phi()-2*TMath::Pi();
    }
     fValueV0[7] = tmpPPhi;
     fValueV0[8] = ePosConvGamGeantLength[i];

    fValueV0[9] = 1./TMath::Sqrt(eNegPart->Pt());
    fValueV0[10] = eNegPart->Eta();

    Double_t tmpNPhi=eNegPart->Phi();
    if( eNegPart->Phi()>TMath::Pi()){
      tmpNPhi = eNegPart->Phi()-2*TMath::Pi();
    }
    fValueV0[11] = tmpNPhi;
    fValueV0[12] = eNegConvGamGeantLength[i];    

    //---- Single track variables----------------------

    fValueV0[13] = (posPt-ePosPart->Pt())/ePosPart->Pt();
    fValueV0[14] = posEta;
    fValueV0[15] = posPhi;
    fValueV0[16] = TMath::Sqrt( bPosSgl[0]* bPosSgl[0] +  bPosSgl[1]* bPosSgl[1] );
    fValueV0[17] = TMath::Sqrt( bPosSgl[0]* bPosSgl[0] +  bPosSgl[1]* bPosSgl[1] )/TMath::Sqrt(bPosCov[0]*bPosCov[0]+bPosCov[2]*bPosCov[2]);
    fValueV0[18] = nClsITSPos;
    fValueV0[19] = nClsTPCPos;
    fValueV0[20] = statusSingPos;


    fValueV0[21] = (negPt-eNegPart->Pt())/eNegPart->Pt();
    fValueV0[22] = negEta;
    fValueV0[23] = negPhi;
    fValueV0[24] = TMath::Sqrt( bNegSgl[0]* bNegSgl[0] +  bNegSgl[1]* bNegSgl[1] );
    fValueV0[25] = TMath::Sqrt( bNegSgl[0]* bNegSgl[0] +  bNegSgl[1]* bNegSgl[1] )/TMath::Sqrt(bNegCov[0]*bNegCov[0]+bNegCov[2]*bNegCov[2]);

    fValueV0[26] = nClsITSNeg;
    fValueV0[27] = nClsTPCNeg;
    fValueV0[28] = statusSingNeg;

    
    //---- V0 track variables----------------------

    fValueV0[29] = (posV0Pt-ePosPart->Pt())/ePosPart->Pt();
    fValueV0[30] = posV0Eta;
    fValueV0[31] = posV0Phi;
    fValueV0[32] = statusV0Pos;

    fValueV0[33] = (negV0Pt-eNegPart->Pt())/eNegPart->Pt();
    fValueV0[34] = negV0Eta;
    fValueV0[35] = negV0Phi;
    fValueV0[36] = statusV0Neg;

    fSparseV0->Fill(fValueV0);
  }


}

void AliAnalysisTaskV0QA::FillHnSparseK0()
{

  Double_t massPi=0.13957018;
  Double_t ppSgl[3];
  Double_t pmSgl[3];
  Float_t bPosSgl[2];
  Float_t bNegSgl[2];
  Float_t bPosCov[3];
  Float_t bNegCov[3];
  
  Double_t ppV0[3];
  Double_t pmV0[3];
  Double_t xrG[3];
  
  TLorentzVector posSglTrack;
  TLorentzVector negSglTrack;
  Double_t posPt,posEta,posPhi;
  Double_t negPt,negEta,negPhi;

  TLorentzVector posV0Track;
  TLorentzVector negV0Track;
  Double_t posV0Pt,posV0Eta,posV0Phi;
  Double_t negV0Pt,negV0Eta,negV0Phi;
  
  Float_t nClsITSPos=-1;
  Float_t nClsITSNeg=-1;

  Float_t nClsTPCPos=-1;
  Float_t nClsTPCNeg=-1;

  Int_t statusSingPos=-1;
  Int_t statusSingNeg=-1;

  Int_t statusV0Pos=-1;
  Int_t statusV0Neg=-1;
  
  for(Int_t i=0;i<nDecayK0Geant+1;i++){
    TParticle* K0Part = stack->Particle(K0DecayK0GeantIndex[i]);
    TParticle* ePosPart = stack->Particle(piPosDecayK0GeantIndex[i]);
    TParticle* eNegPart = stack->Particle(piNegDecayK0GeantIndex[i]);
    if (piPosDecayK0SingleRecIndex[i]!=-1){
      AliESDtrack * ePosSglTrack = fESD->GetTrack(piPosDecayK0SingleRecIndex[i]);
      ePosSglTrack->GetPxPyPz(ppSgl); 
      posSglTrack.SetXYZM(ppSgl[0],ppSgl[1],ppSgl[2],massPi);
      posPt  = posSglTrack.Pt();
      posEta = posSglTrack.Eta();
      posPhi = posSglTrack.Phi();
      ePosSglTrack->GetImpactParameters(bPosSgl,bPosCov);
      nClsITSPos=ePosSglTrack->GetNcls(0);
      nClsTPCPos=ePosSglTrack->GetNcls(1);
      statusSingPos=1;
    }else{
      posPt  = 1000000;
      posEta = -2.;
      posPhi = -2*TMath::Pi();
      bPosSgl[0]=-100.;
      bPosSgl[1]=-100.;
      bPosCov[0]=-100;
      bPosCov[2]=-100;

      nClsITSPos=-1;
      nClsTPCPos=-1;
      statusSingPos=-1;
    }

    if (piNegDecayK0SingleRecIndex[i]!=-1){
      AliESDtrack * eNegSglTrack = fESD->GetTrack(piNegDecayK0SingleRecIndex[i]);
      eNegSglTrack->GetPxPyPz(pmSgl); 
      negSglTrack.SetXYZM(pmSgl[0],pmSgl[1],pmSgl[2],massPi);
      negPt  = negSglTrack.Pt();
      negEta = negSglTrack.Eta();
      negPhi = negSglTrack.Phi();
      eNegSglTrack->GetImpactParameters(bNegSgl,bNegCov);
      nClsITSNeg=eNegSglTrack->GetNcls(0);
      nClsTPCNeg=eNegSglTrack->GetNcls(1);
      statusSingNeg=1;
    }else{
      negPt  = 1000000;
      negEta = -2.;
      negPhi = -2*TMath::Pi();
      bNegSgl[0]=-100.;
      bNegSgl[1]=-100.;
      bNegCov[0]=-100;
      bNegCov[2]=-100;
      nClsITSNeg=-1;
      nClsTPCNeg=-1;
      statusSingNeg=-1;
    }

    posV0Pt  = 1000000;
    posV0Eta = -2.;
    posV0Phi = -2*TMath::Pi();
    negV0Pt  = 1000000;
    negV0Eta = -2.;
    negV0Phi = -2*TMath::Pi();

    if(DecayK0V0RecIndexPos[i]!=-1){
      AliESDv0 * fV0MIs = fESD->GetV0(DecayK0V0RecIndexPos[i]);
      AliESDtrack* trackPosTest = fESD->GetTrack(fV0MIs->GetPindex());
      AliESDtrack* trackNegTest = fESD->GetTrack(fV0MIs->GetNindex());

      if (piPosDecayK0V0RecIndex[i]!=-1 ){
	//	AliESDtrack * ePosV0Track = fESD->GetTrack(piPosDecayK0V0RecIndex[i]);
	if ( trackPosTest->GetSign()==1 ) {
	  fV0MIs->GetPPxPyPz(ppV0[0],ppV0[1],ppV0[2]);
	}else{
	  fV0MIs->GetNPxPyPz(ppV0[0],ppV0[1],ppV0[2]);
	}
	posV0Track.SetXYZM(ppV0[0],ppV0[1],ppV0[2],massPi);

	posV0Pt  = posV0Track.Pt();
	posV0Eta = posV0Track.Eta();
	posV0Phi = posV0Track.Phi();
	statusV0Pos=1;
      }else{
	posV0Pt  = 1000000;
	posV0Eta = -2.;
	posV0Phi = -2*TMath::Pi();
	statusV0Pos=-1;
      }
      
      if (piNegDecayK0V0RecIndex[i]!=-1 ){
	//	AliESDtrack * eNegV0Track = fESD->GetTrack(piNegDecayK0V0RecIndex[i]);
	if ( trackNegTest->GetSign()==-1 ) {
	  fV0MIs->GetNPxPyPz(pmV0[0],pmV0[1],pmV0[2]);
	}else{
	  fV0MIs->GetPPxPyPz(pmV0[0],pmV0[1],pmV0[2]);
	}
	negV0Track.SetXYZM(pmV0[0],pmV0[1],pmV0[2],massPi);
	
	negV0Pt  = negV0Track.Pt();
	negV0Eta = negV0Track.Eta();
	negV0Phi = negV0Track.Phi();
	statusV0Neg=1;
      }else{
	negV0Pt  = 1000000;
	negV0Eta = -2.;
	negV0Phi = -2*TMath::Pi();
	statusV0Neg=-1;
      }
    }

    xrG[0] = ePosPart->Vx();
    xrG[1] = ePosPart->Vy();
    xrG[2] = ePosPart->Vz();
    
    
    //--------- Geant variables ----------------------
    fValueK0[0] = 1./TMath::Sqrt(K0Part->Pt());
    fValueK0[1] = K0Part->Eta();

    Double_t tmpGPhi=K0Part->Phi();
    if( K0Part->Phi()>TMath::Pi()){
      tmpGPhi=K0Part->Phi()-2*TMath::Pi();
    }
    fValueK0[2] = tmpGPhi;

    fValueK0[3] = TMath::Sqrt(xrG[0]*xrG[0]+xrG[1]*xrG[1]);
    fValueK0[4] = xrG[2];
 

    fValueK0[5] = 1./TMath::Sqrt(ePosPart->Pt());
    fValueK0[6] = ePosPart->Eta();

    Double_t tmpPPhi=ePosPart->Phi();
    if( ePosPart->Phi()>TMath::Pi()){
      tmpPPhi = ePosPart->Phi()-2*TMath::Pi();
    }
     fValueK0[7] = tmpPPhi;
     fValueK0[8] = piPosDecayK0GeantLength[i];

    fValueK0[9] = 1./TMath::Sqrt(eNegPart->Pt());
    fValueK0[10] = eNegPart->Eta();

    Double_t tmpNPhi=eNegPart->Phi();
    if( eNegPart->Phi()>TMath::Pi()){
      tmpNPhi = eNegPart->Phi()-2*TMath::Pi();
    }
    fValueK0[11] = tmpNPhi;
    fValueK0[12] = piNegDecayK0GeantLength[i];    
    //---- Single track variables----------------------

    fValueK0[13] = (posPt-ePosPart->Pt())/ePosPart->Pt() ;
    fValueK0[14] = posEta;
    fValueK0[15] = posPhi;
    fValueK0[16] = TMath::Sqrt( bPosSgl[0]* bPosSgl[0] +  bPosSgl[1]* bPosSgl[1] );
    fValueK0[17] = TMath::Sqrt( bPosSgl[0]* bPosSgl[0] +  bPosSgl[1]* bPosSgl[1] )/TMath::Sqrt(bPosCov[0]*bPosCov[0]+bPosCov[2]*bPosCov[2]);
   
    fValueK0[18] = nClsITSPos;
    fValueK0[19] = nClsTPCPos;
    fValueK0[20] = statusSingPos;

    fValueK0[21] = (negPt-eNegPart->Pt())/eNegPart->Pt();
    fValueK0[22] = negEta;
    fValueK0[23] = negPhi;
    fValueK0[24] = TMath::Sqrt( bNegSgl[0]* bNegSgl[0]+  bNegSgl[1]* bNegSgl[1] );
    fValueK0[25] = TMath::Sqrt( bNegSgl[0]* bNegSgl[0]+  bNegSgl[1]* bNegSgl[1] )/TMath::Sqrt(bNegCov[0]*bNegCov[0]+bNegCov[2]*bNegCov[2]);
    fValueK0[26] = nClsITSNeg;
    fValueK0[27] = nClsTPCNeg;
    fValueK0[28] = statusSingNeg;

    
    //---- V0 track variables----------------------

    fValueK0[29] = (posV0Pt-ePosPart->Pt())/ePosPart->Pt();
    fValueK0[30] = posV0Eta;
    fValueK0[31] = posV0Phi;
    fValueK0[32] = statusV0Pos;

    fValueK0[33] = (negV0Pt-eNegPart->Pt())/eNegPart->Pt();
    fValueK0[34] = negV0Eta;
    fValueK0[35] = negV0Phi;
    fValueK0[36] = statusV0Neg;

    fSparseK0->Fill(fValueK0);
  }
  

}
void AliAnalysisTaskV0QA::FillHnSparseL()
{

  Double_t massPi=0.13957018;
  Double_t massP=0.93827203;


  Double_t ppSgl[3];
  Double_t pmSgl[3];
  Float_t bPosSgl[2];
  Float_t bNegSgl[2];
  Float_t bPosCov[3];
  Float_t bNegCov[3];
  
  Double_t ppV0[3];
  Double_t pmV0[3];
  Double_t xrG[3];
  
  TLorentzVector posSglTrack;
  TLorentzVector negSglTrack;
  Double_t posPt,posEta,posPhi;
  Double_t negPt,negEta,negPhi;

  TLorentzVector posV0Track;
  TLorentzVector negV0Track;
  Double_t posV0Pt,posV0Eta,posV0Phi;
  Double_t negV0Pt,negV0Eta,negV0Phi;
  
  Float_t nClsITSPos=-1;
  Float_t nClsITSNeg=-1;

  Float_t nClsTPCPos=-1;
  Float_t nClsTPCNeg=-1;

  Int_t statusSingPos=-1;
  Int_t statusSingNeg=-1;

  Int_t statusV0Pos=-1;
  Int_t statusV0Neg=-1;

  for(Int_t i=0;i<nDecayLGeant+1;i++){
    TParticle* lPart = stack->Particle(lDecayLGeantIndex[i]);
    TParticle* ePosPart = stack->Particle(pPosDecayLGeantIndex[i]);
    TParticle* eNegPart = stack->Particle(piNegDecayLGeantIndex[i]);
    if (pPosDecayLSingleRecIndex[i]!=-1){
      AliESDtrack * ePosSglTrack = fESD->GetTrack(pPosDecayLSingleRecIndex[i]);
      ePosSglTrack->GetPxPyPz(ppSgl); 
      posSglTrack.SetXYZM(ppSgl[0],ppSgl[1],ppSgl[2],massP);
      posPt  = posSglTrack.Pt();
      posEta = posSglTrack.Eta();
      posPhi = posSglTrack.Phi();
      ePosSglTrack->GetImpactParameters(bPosSgl,bPosCov);
      nClsITSPos=ePosSglTrack->GetNcls(0);
      nClsTPCPos=ePosSglTrack->GetNcls(1);
      statusSingPos=1;
    }else{
      posPt  = 1000000;
      posEta = -2.;
      posPhi = -2*TMath::Pi();
      bPosSgl[0]=-100.;
      bPosSgl[1]=-100.;
      bPosCov[0]=-100;
      bPosCov[2]=-100;
      nClsITSPos=-1;
      nClsTPCPos=-1;
      statusSingPos=-1;
    }
    
    if (piNegDecayLSingleRecIndex[i]!=-1){
      AliESDtrack * eNegSglTrack = fESD->GetTrack(piNegDecayLSingleRecIndex[i]);
      eNegSglTrack->GetPxPyPz(pmSgl); 
      negSglTrack.SetXYZM(pmSgl[0],pmSgl[1],pmSgl[2],massPi);
      negPt  = negSglTrack.Pt();
      negEta = negSglTrack.Eta();
      negPhi = negSglTrack.Phi();
      eNegSglTrack->GetImpactParameters(bNegSgl,bNegCov);
      nClsITSNeg=eNegSglTrack->GetNcls(0);
      nClsTPCNeg=eNegSglTrack->GetNcls(1);
      statusSingNeg=1;
    }else{
      negPt  = 1000000;
      negEta = -2.;
      negPhi = -2*TMath::Pi();
      bNegSgl[0]=-100.;
      bNegSgl[1]=-100.;
      bNegCov[0]=-100;
      bNegCov[2]=-100;
      nClsITSNeg=-1;
      nClsTPCNeg=-1;
      statusSingNeg=-1;
    }

    posV0Pt  = 1000000;
    posV0Eta = -2.;
    posV0Phi = -2*TMath::Pi();
    negV0Pt  = 1000000;
    negV0Eta = -2.;
    negV0Phi = -2*TMath::Pi();
    
    if(DecayLV0RecIndexPos[i]!=-1){
      AliESDv0 * fV0MIs = fESD->GetV0(DecayLV0RecIndexPos[i]);
      AliESDtrack* trackPosTest = fESD->GetTrack(fV0MIs->GetPindex());
      AliESDtrack* trackNegTest = fESD->GetTrack(fV0MIs->GetNindex());

      if (pPosDecayLV0RecIndex[i]!=-1 ){
	//	AliESDtrack * ePosV0Track = fESD->GetTrack(pPosDecayLV0RecIndex[i]);
	if ( trackPosTest->GetSign()==1 ) {
	  fV0MIs->GetPPxPyPz(ppV0[0],ppV0[1],ppV0[2]);
	}else{
	  fV0MIs->GetNPxPyPz(ppV0[0],ppV0[1],ppV0[2]);
	}
	posV0Track.SetXYZM(ppV0[0],ppV0[1],ppV0[2],massP);

	posV0Pt  = posV0Track.Pt();
	posV0Eta = posV0Track.Eta();
	posV0Phi = posV0Track.Phi();
	statusV0Pos=1;
      }else{
	posV0Pt  = 1000000;
	posV0Eta = -2.;
	posV0Phi = -2*TMath::Pi();
	statusV0Pos=-1;
      }
      
      if (piNegDecayLV0RecIndex[i]!=-1 ){
	//	AliESDtrack * eNegV0Track = fESD->GetTrack(piNegDecayLV0RecIndex[i]);
	if ( trackNegTest->GetSign()==-1 ) {
	  fV0MIs->GetNPxPyPz(pmV0[0],pmV0[1],pmV0[2]);
	}else{
	  fV0MIs->GetPPxPyPz(pmV0[0],pmV0[1],pmV0[2]);
	}
	negV0Track.SetXYZM(pmV0[0],pmV0[1],pmV0[2],massPi);
	
	negV0Pt  = negV0Track.Pt();
	negV0Eta = negV0Track.Eta();
	negV0Phi = negV0Track.Phi();
	statusV0Neg=1;
      }else{
	negV0Pt  = 1000000;
	negV0Eta = -2.;
	negV0Phi = -2*TMath::Pi();
	statusV0Neg=-1;
      }
    }
    
    xrG[0] = ePosPart->Vx();
    xrG[1] = ePosPart->Vy();
    xrG[2] = ePosPart->Vz();
    
    //--------- Geant variables ----------------------
    fValueL[0] = 1./TMath::Sqrt(lPart->Pt());
    fValueL[1] = lPart->Eta();

    Double_t tmpGPhi=lPart->Phi();
    if( lPart->Phi()>TMath::Pi()){
      tmpGPhi=lPart->Phi()-2*TMath::Pi();
    }
    fValueL[2] = tmpGPhi;

    fValueL[3] = TMath::Sqrt(xrG[0]*xrG[0]+xrG[1]*xrG[1]);
    fValueL[4] = xrG[2];
 

    fValueL[5] = 1./TMath::Sqrt(ePosPart->Pt());
    fValueL[6] = ePosPart->Eta();

    Double_t tmpPPhi=ePosPart->Phi();
    if( ePosPart->Phi()>TMath::Pi()){
      tmpPPhi = ePosPart->Phi()-2*TMath::Pi();
    }
     fValueL[7] = tmpPPhi;
     fValueL[8] = pPosDecayLGeantLength[i];

    fValueL[9] = 1./TMath::Sqrt(eNegPart->Pt());
    fValueL[10] = eNegPart->Eta();

    Double_t tmpNPhi=eNegPart->Phi();
    if( eNegPart->Phi()>TMath::Pi()){
      tmpNPhi = eNegPart->Phi()-2*TMath::Pi();
    }
    fValueL[11] = tmpNPhi;
    fValueL[12] = piNegDecayLGeantLength[i];    
    //---- Single track variables----------------------

    fValueL[13] = (posPt-ePosPart->Pt())/ePosPart->Pt();
    fValueL[14] = posEta;
    fValueL[15] = posPhi;
    fValueL[16] = TMath::Sqrt( bPosSgl[0]* bPosSgl[0] +  bPosSgl[1]* bPosSgl[1]);
    fValueL[17] = TMath::Sqrt( bPosSgl[0]* bPosSgl[0] +  bPosSgl[1]* bPosSgl[1] )/TMath::Sqrt(bPosCov[0]*bPosCov[0]+bPosCov[2]*bPosCov[2]);   
    fValueL[18] = nClsITSPos;
    fValueL[19] = nClsTPCPos;
    fValueL[20] = statusSingPos;

    fValueL[21] = (negPt-eNegPart->Pt())/eNegPart->Pt() ;
    fValueL[22] = negEta;
    fValueL[23] = negPhi;
    fValueL[24] = TMath::Sqrt( bNegSgl[0]* bNegSgl[0] +  bNegSgl[1]* bNegSgl[1] );
    fValueL[25] = TMath::Sqrt( bNegSgl[0]* bNegSgl[0] +  bNegSgl[1]* bNegSgl[1] )/TMath::Sqrt(bNegCov[0]*bNegCov[0]+bNegCov[2]*bNegCov[2]);
    fValueL[26] = nClsITSNeg;
    fValueL[27] = nClsTPCNeg;
    fValueL[28] = statusSingNeg;


    
    //---- V0 track variables----------------------

    fValueL[29] = (posV0Pt-ePosPart->Pt())/ePosPart->Pt();
    fValueL[30] = posV0Eta;
    fValueL[31] = posV0Phi;
    fValueL[32] = statusV0Pos;


    fValueL[33] = (negV0Pt-eNegPart->Pt())/eNegPart->Pt();
    fValueL[34] = negV0Eta;
    fValueL[35] = negV0Phi;
    fValueL[36] = statusV0Neg;

    fSparseL->Fill(fValueL);
  }


}

void AliAnalysisTaskV0QA::FillHnSparseAL()
{

  Double_t massPi=0.13957018;
  Double_t massP=0.93827203;


  Double_t ppSgl[3];
  Double_t pmSgl[3];
  Float_t bPosSgl[2];
  Float_t bNegSgl[2];
  Float_t bPosCov[3];
  Float_t bNegCov[3];
  
  Double_t ppV0[3];
  Double_t pmV0[3];
  Double_t xrG[3];
  
  TLorentzVector posSglTrack;
  TLorentzVector negSglTrack;
  Double_t posPt,posEta,posPhi;
  Double_t negPt,negEta,negPhi;

  TLorentzVector posV0Track;
  TLorentzVector negV0Track;
  Double_t posV0Pt,posV0Eta,posV0Phi;
  Double_t negV0Pt,negV0Eta,negV0Phi;
  
  Float_t nClsITSPos=-1;
  Float_t nClsITSNeg=-1;

  Float_t nClsTPCPos=-1;
  Float_t nClsTPCNeg=-1;

  Int_t statusSingPos=-1;
  Int_t statusSingNeg=-1;

  Int_t statusV0Pos=-1;
  Int_t statusV0Neg=-1;


  for(Int_t i=0;i<nDecayALGeant+1;i++){
    TParticle* alPart = stack->Particle(alDecayALGeantIndex[i]);
    TParticle* eNegPart = stack->Particle(apNegDecayALGeantIndex[i]);
    TParticle* ePosPart = stack->Particle(piPosDecayALGeantIndex[i]);
    if (piPosDecayALSingleRecIndex[i]!=-1){
      AliESDtrack * ePosSglTrack = fESD->GetTrack(piPosDecayALSingleRecIndex[i]);
      ePosSglTrack->GetPxPyPz(ppSgl); 
      posSglTrack.SetXYZM(ppSgl[0],ppSgl[1],ppSgl[2],massPi);
      posPt  = posSglTrack.Pt();
      posEta = posSglTrack.Eta();
      posPhi = posSglTrack.Phi();
      ePosSglTrack->GetImpactParameters(bPosSgl,bPosCov);
      nClsITSPos=ePosSglTrack->GetNcls(0);
      nClsTPCPos=ePosSglTrack->GetNcls(1);
      statusSingPos=1;
    }else{
      posPt  = 1000000;
      posEta = -2.;
      posPhi = -2*TMath::Pi();
      bPosSgl[0]=-100.;
      bPosSgl[1]=-100.;
      bPosCov[0]=-100;
      bPosCov[2]=-100;
      nClsITSPos=-1;
      nClsTPCPos=-1;
      statusSingPos=-1;
    }
    
    if (apNegDecayALSingleRecIndex[i]!=-1){
      AliESDtrack * eNegSglTrack = fESD->GetTrack(apNegDecayALSingleRecIndex[i]);
      eNegSglTrack->GetPxPyPz(pmSgl); 
      negSglTrack.SetXYZM(pmSgl[0],pmSgl[1],pmSgl[2],massP);
      negPt  = negSglTrack.Pt();
      negEta = negSglTrack.Eta();
      negPhi = negSglTrack.Phi();
      eNegSglTrack->GetImpactParameters(bNegSgl,bNegCov);
      nClsITSNeg=eNegSglTrack->GetNcls(0);
      nClsTPCNeg=eNegSglTrack->GetNcls(1);
      statusSingNeg=1;
    }else{
      negPt  = 1000000;
      negEta = -2.;
      negPhi = -2*TMath::Pi();
      bNegSgl[0]=-100.;
      bNegSgl[1]=-100.;
      bNegCov[0]=-100;
      bNegCov[2]=-100;
      nClsITSNeg=-1;
      nClsTPCNeg=-1;
      statusSingNeg=-1;
    }

    posV0Pt  = 1000000;
    posV0Eta = -2.;
    posV0Phi = -2*TMath::Pi();
    negV0Pt  = 1000000;
    negV0Eta = -2.;
    negV0Phi = -2*TMath::Pi();
    
    if(DecayALV0RecIndexPos[i]!=-1){
      AliESDv0 * fV0MIs = fESD->GetV0(DecayALV0RecIndexPos[i]);
      AliESDtrack* trackPosTest = fESD->GetTrack(fV0MIs->GetPindex());
      AliESDtrack* trackNegTest = fESD->GetTrack(fV0MIs->GetNindex());

      if (piPosDecayALV0RecIndex[i]!=-1 ){
	//	AliESDtrack * ePosV0Track = fESD->GetTrack(piPosDecayALV0RecIndex[i]);
	if ( trackPosTest->GetSign()==1 ) {
	  fV0MIs->GetPPxPyPz(ppV0[0],ppV0[1],ppV0[2]);
	}else{
	  fV0MIs->GetNPxPyPz(ppV0[0],ppV0[1],ppV0[2]);
	}
	posV0Track.SetXYZM(ppV0[0],ppV0[1],ppV0[2],massPi);

	posV0Pt  = posV0Track.Pt();
	posV0Eta = posV0Track.Eta();
	posV0Phi = posV0Track.Phi();
	statusV0Pos=1;
      }else{
	posV0Pt  = 1000000;
	posV0Eta = -2.;
	posV0Phi = -2*TMath::Pi();
	statusV0Pos=-1;
      }
      
      if (apNegDecayALV0RecIndex[i]!=-1 ){
	//	AliESDtrack * eNegV0Track = fESD->GetTrack(apNegDecayALV0RecIndex[i]);
	if ( trackNegTest->GetSign()==-1 ) {
	  fV0MIs->GetNPxPyPz(pmV0[0],pmV0[1],pmV0[2]);
	}else{
	  fV0MIs->GetPPxPyPz(pmV0[0],pmV0[1],pmV0[2]);
	}
	negV0Track.SetXYZM(pmV0[0],pmV0[1],pmV0[2],massP);
	
	negV0Pt  = negV0Track.Pt();
	negV0Eta = negV0Track.Eta();
	negV0Phi = negV0Track.Phi();
	statusV0Neg=1;
      }else{
	negV0Pt  = 1000000;
	negV0Eta = -2.;
	negV0Phi = -2*TMath::Pi();
	statusV0Neg=-1;
      }
    }
    
    xrG[0] = ePosPart->Vx();
    xrG[1] = ePosPart->Vy();
    xrG[2] = ePosPart->Vz();
    
    //--------- Geant variables ----------------------
    fValueAL[0] = 1./TMath::Sqrt(alPart->Pt());
    fValueAL[1] = alPart->Eta();

    Double_t tmpGPhi=alPart->Phi();
    if( alPart->Phi()>TMath::Pi()){
      tmpGPhi=alPart->Phi()-2*TMath::Pi();
    }
    fValueAL[2] = tmpGPhi;

    fValueAL[3] = TMath::Sqrt(xrG[0]*xrG[0]+xrG[1]*xrG[1]);
    fValueAL[4] = xrG[2];
 

    fValueAL[5] = 1./TMath::Sqrt(ePosPart->Pt());
    fValueAL[6] = ePosPart->Eta();

    Double_t tmpPPhi=ePosPart->Phi();
    if( ePosPart->Phi()>TMath::Pi()){
      tmpPPhi = ePosPart->Phi()-2*TMath::Pi();
    }
     fValueAL[7] = tmpPPhi;
     fValueAL[8] = piPosDecayALGeantLength[i];

    fValueAL[9] = 1./TMath::Sqrt(eNegPart->Pt());
    fValueAL[10] = eNegPart->Eta();

    Double_t tmpNPhi=eNegPart->Phi();
    if( eNegPart->Phi()>TMath::Pi()){
      tmpNPhi = eNegPart->Phi()-2*TMath::Pi();
    }
    fValueAL[11] = tmpNPhi;
    fValueAL[12] = apNegDecayALGeantLength[i];    
    //---- Single track variables----------------------

    fValueAL[13] = (posPt-ePosPart->Pt())/ePosPart->Pt();
    fValueAL[14] = posEta;
    fValueAL[15] = posPhi;
    fValueAL[16] = TMath::Sqrt( bPosSgl[0]* bPosSgl[0] +  bPosSgl[1]* bPosSgl[1]);
    fValueAL[17] = TMath::Sqrt( bPosSgl[0]* bPosSgl[0] +  bPosSgl[1]* bPosSgl[1] )/TMath::Sqrt(bPosCov[0]*bPosCov[0]+bPosCov[2]*bPosCov[2]);   
    fValueAL[18] = nClsITSPos;
    fValueAL[19] = nClsTPCPos;
    fValueAL[20] = statusSingPos;

    fValueAL[21] = (negPt-eNegPart->Pt())/eNegPart->Pt() ;
    fValueAL[22] = negEta;
    fValueAL[23] = negPhi;
    fValueAL[24] = TMath::Sqrt( bNegSgl[0]* bNegSgl[0] +  bNegSgl[1]* bNegSgl[1] );
    fValueAL[25] = TMath::Sqrt( bNegSgl[0]* bNegSgl[0] +  bNegSgl[1]* bNegSgl[1] )/TMath::Sqrt(bNegCov[0]*bNegCov[0]+bNegCov[2]*bNegCov[2]);
    fValueAL[26] = nClsITSNeg;
    fValueAL[27] = nClsTPCNeg;
    fValueAL[28] = statusSingNeg;


    
    //---- V0 track variables----------------------

    fValueAL[29] = (posV0Pt-ePosPart->Pt())/ePosPart->Pt();
    fValueAL[30] = posV0Eta;
    fValueAL[31] = posV0Phi;
    fValueAL[32] = statusV0Pos;


    fValueAL[33] = (negV0Pt-eNegPart->Pt())/eNegPart->Pt();
    fValueAL[34] = negV0Eta;
    fValueAL[35] = negV0Phi;
    fValueAL[36] = statusV0Neg;

    fSparseAL->Fill(fValueAL);
  }
}


// void AliAnalysisTaskV0QA::SetESDtrackCuts()
// {

//   fEsdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");

//   fEsdTrackCuts->SetRequireTPCRefit(kTRUE);
//   fEsdTrackCuts->SetAcceptKinkDaughters(kFALSE);



// }
