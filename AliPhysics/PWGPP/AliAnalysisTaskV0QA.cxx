//------------------------------------------------
// Implementation of AliAnalysisTaskV0QA class.
// Calculates the "on the fly" V0 method efficiency 
// for Gamma, K0s, lambda, antilambda
// Needs MC information
// Author: A. Marin  Revision 18/10/09
//-------------------------------------------------
#define AliAnalysisTaskV0QA_cxx

#include "TChain.h"
#include "TH1.h"
#include "TH2.h"


//#include "TLorentzVector.h"

#include "AliAnalysisTaskSE.h"
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
#include "AliCentrality.h"
#include "AliAnalysisTaskV0QA.h"

ClassImp(AliAnalysisTaskV0QA)
AliAnalysisTaskV0QA::AliAnalysisTaskV0QA() :AliAnalysisTaskSE(), 
fESD(0), 
fStack(0),
fMCtruth(0),
fChain(0),
fOutputContainer(0),
fSparseV0(0),
fSparseK0(0),
fSparseL(0),
fSparseAL(0),
fnEv(0),
fgDim(500),
fnConvGamGeant(-1),
fgConvGamGeantIndex(0),
feNegConvGamGeantIndex(0),
fePosConvGamGeantIndex(0),
feNegConvGamGeantLength(0),
fePosConvGamGeantLength(0),
feNegConvGamSingleRecIndex(0),
fePosConvGamSingleRecIndex(0),
feNegConvGamV0RecIndex(0),
fePosConvGamV0RecIndex(0),
fConvGamV0RecIndexPos(0),
fConvGamV0RecIndexNeg(0),
fnDecayLGeant(-1),
flDecayLGeantIndex(0),
fpiNegDecayLGeantIndex(0),
fpPosDecayLGeantIndex(0),
fpiNegDecayLGeantLength(0),
fpPosDecayLGeantLength(0),
fpiNegDecayLSingleRecIndex(0),
fpPosDecayLSingleRecIndex(0),
fpiNegDecayLV0RecIndex(0),
fpPosDecayLV0RecIndex(0),
fDecayLV0RecIndexPos(0),
fDecayLV0RecIndexNeg(0),
fnDecayALGeant(-1),
falDecayALGeantIndex(0),
fpiPosDecayALGeantIndex(0),
fapNegDecayALGeantIndex(0),
fpiPosDecayALGeantLength(0),
fapNegDecayALGeantLength(0),
fpiPosDecayALSingleRecIndex(0),
fapNegDecayALSingleRecIndex(0),
fpiPosDecayALV0RecIndex(0),
fapNegDecayALV0RecIndex(0),
fDecayALV0RecIndexPos(0),
fDecayALV0RecIndexNeg(0),
fnDecayK0Geant(-1),
fK0DecayK0GeantIndex(0),
fpiNegDecayK0GeantIndex(0),
fpiPosDecayK0GeantIndex(0),
fpiNegDecayK0GeantLength(0),
fpiPosDecayK0GeantLength(0),
fpiNegDecayK0SingleRecIndex(0),
fpiPosDecayK0SingleRecIndex(0),
fpiNegDecayK0V0RecIndex(0),
fpiPosDecayK0V0RecIndex(0),
fDecayK0V0RecIndexPos(0),
fDecayK0V0RecIndexNeg(0),
fpiPosK0Index(-1),
fpiNegK0Index(-1),
fnTracksPrim(-1),
ftpcRefit(0),
fitsRefit(0),
ftrdRefit(0),
ftrdOut(0),
fDim(38),
fValueL(0),
fValueAL(0),
fValueK0(0),
fValueV0(0),
fxminV0(0),
fxmaxV0(0),
fbinsV0(0),
fCentralityC(-1),
fRefTPC(0),
fclRefsN(0),
fclRefsP(0)

 {
  // Default Constructor.
  for(Int_t i=0;i<100000;i++) fLabelsTPC[i] = 0;
 }

//________________________________________________________________________
AliAnalysisTaskV0QA::AliAnalysisTaskV0QA(const char *name) :AliAnalysisTaskSE(name), 
fESD(0), 
fStack(0),
fMCtruth(0),
fChain(0),
fOutputContainer(0),
fSparseV0(0),
fSparseK0(0),
fSparseL(0),
fSparseAL(0),
fnEv(0),
fgDim(500),
fnConvGamGeant(-1),
fgConvGamGeantIndex(0),
feNegConvGamGeantIndex(0),
fePosConvGamGeantIndex(0),
feNegConvGamGeantLength(0),
fePosConvGamGeantLength(0),
feNegConvGamSingleRecIndex(0),
fePosConvGamSingleRecIndex(0),
feNegConvGamV0RecIndex(0),
fePosConvGamV0RecIndex(0),
fConvGamV0RecIndexPos(0),
fConvGamV0RecIndexNeg(0),
fnDecayLGeant(-1),
flDecayLGeantIndex(0),
fpiNegDecayLGeantIndex(0),
fpPosDecayLGeantIndex(0),
fpiNegDecayLGeantLength(0),
fpPosDecayLGeantLength(0),
fpiNegDecayLSingleRecIndex(0),
fpPosDecayLSingleRecIndex(0),
fpiNegDecayLV0RecIndex(0),
fpPosDecayLV0RecIndex(0),
fDecayLV0RecIndexPos(0),
fDecayLV0RecIndexNeg(0),
fnDecayALGeant(-1),
falDecayALGeantIndex(0),
fpiPosDecayALGeantIndex(0),
fapNegDecayALGeantIndex(0),
fpiPosDecayALGeantLength(0),
fapNegDecayALGeantLength(0),
fpiPosDecayALSingleRecIndex(0),
fapNegDecayALSingleRecIndex(0),
fpiPosDecayALV0RecIndex(0),
fapNegDecayALV0RecIndex(0),
fDecayALV0RecIndexPos(0),
fDecayALV0RecIndexNeg(0),
fnDecayK0Geant(-1),
fK0DecayK0GeantIndex(0),
fpiNegDecayK0GeantIndex(0),
fpiPosDecayK0GeantIndex(0),
fpiNegDecayK0GeantLength(0),
fpiPosDecayK0GeantLength(0),
fpiNegDecayK0SingleRecIndex(0),
fpiPosDecayK0SingleRecIndex(0),
fpiNegDecayK0V0RecIndex(0),
fpiPosDecayK0V0RecIndex(0),
fDecayK0V0RecIndexPos(0),
fDecayK0V0RecIndexNeg(0),
fpiPosK0Index(-1),
fpiNegK0Index(-1),
fnTracksPrim(-1),
ftpcRefit(0),
fitsRefit(0),
ftrdRefit(0),
ftrdOut(0),
fDim(38),
fValueL(0),
fValueAL(0),
fValueK0(0),
fValueV0(0),
fxminV0(0),
fxmaxV0(0),
fbinsV0(0),
fCentralityC(-1),
fRefTPC(0),
fclRefsN(0),
fclRefsP(0)

 {

   fnEv=0;
   fDim=38; 

   fValueK0 = new Double_t[fDim];
   fValueL = new Double_t[fDim];
   fValueAL = new Double_t[fDim];
   fValueV0 = new Double_t[fDim];
   fxminV0 = new Double_t[fDim];
   fxmaxV0 = new Double_t[fDim];
   fbinsV0 = new Int_t[fDim];

  for(Int_t i=0;i<100000;i++) fLabelsTPC[i] = 0;

  fgDim=500; 
  fgConvGamGeantIndex = new Int_t[fgDim];
  feNegConvGamGeantIndex = new Int_t[fgDim];
  fePosConvGamGeantIndex = new Int_t[fgDim];
  feNegConvGamGeantLength = new Float_t[fgDim];
  fePosConvGamGeantLength = new Float_t[fgDim];

  feNegConvGamSingleRecIndex = new Int_t[fgDim];
  fePosConvGamSingleRecIndex = new Int_t[fgDim];

  feNegConvGamV0RecIndex = new Int_t[fgDim];
  fePosConvGamV0RecIndex = new Int_t[fgDim];

  fConvGamV0RecIndexPos = new Int_t[fgDim];
  fConvGamV0RecIndexNeg = new Int_t[fgDim];

  // Lambda to proton pi-
  flDecayLGeantIndex = new Int_t[fgDim];
  fpiNegDecayLGeantIndex = new Int_t[fgDim];
  fpPosDecayLGeantIndex = new Int_t[fgDim];
  fpiNegDecayLGeantLength = new Float_t[fgDim];
  fpPosDecayLGeantLength = new Float_t[fgDim];

  fpiNegDecayLSingleRecIndex = new Int_t[fgDim];
  fpPosDecayLSingleRecIndex = new Int_t[fgDim];

  fpiNegDecayLV0RecIndex = new Int_t[fgDim];
  fpPosDecayLV0RecIndex = new Int_t[fgDim];

  fDecayLV0RecIndexPos = new Int_t[fgDim];
  fDecayLV0RecIndexNeg = new Int_t[fgDim];

  //K0S to pi+ pi-
  fK0DecayK0GeantIndex = new Int_t[fgDim];
  fpiNegDecayK0GeantIndex = new Int_t[fgDim];
  fpiPosDecayK0GeantIndex = new Int_t[fgDim];
  fpiNegDecayK0GeantLength = new Float_t[fgDim];
  fpiPosDecayK0GeantLength = new Float_t[fgDim];

  fpiNegDecayK0SingleRecIndex = new Int_t[fgDim];
  fpiPosDecayK0SingleRecIndex = new Int_t[fgDim];

  fpiNegDecayK0V0RecIndex = new Int_t[fgDim];
  fpiPosDecayK0V0RecIndex = new Int_t[fgDim];

  fDecayK0V0RecIndexPos = new Int_t[fgDim];
  fDecayK0V0RecIndexNeg = new Int_t[fgDim];

  //Antilambda to antiproton piplus
  falDecayALGeantIndex = new Int_t[fgDim];
  fpiPosDecayALGeantIndex = new Int_t[fgDim];
  fapNegDecayALGeantIndex = new Int_t[fgDim];
  fpiPosDecayALGeantLength = new Float_t[fgDim];
  fapNegDecayALGeantLength = new Float_t[fgDim];

  fpiPosDecayALSingleRecIndex = new Int_t[fgDim];
  fapNegDecayALSingleRecIndex = new Int_t[fgDim];

  fpiPosDecayALV0RecIndex = new Int_t[fgDim];
  fapNegDecayALV0RecIndex = new Int_t[fgDim];

  fDecayALV0RecIndexPos = new Int_t[fgDim];
  fDecayALV0RecIndexNeg = new Int_t[fgDim];


 
  fclRefsP = new TClonesArray("AliTrackReference");
  fclRefsN = new TClonesArray("AliTrackReference");

  //  SetESDtrackCuts();


  AliLog::SetGlobalLogLevel(AliLog::kError);
//
  DefineOutput(1, TList::Class());
}

//_____________________________________________________
AliAnalysisTaskV0QA::~AliAnalysisTaskV0QA()
{
  // Remove all pointers


  delete [] fclRefsP;
  delete [] fclRefsN;


  delete [] fValueK0;
  delete [] fValueL;
  delete [] fValueAL;
  delete [] fValueV0;
  delete [] fbinsV0;
  delete [] fxminV0;
  delete [] fxmaxV0;

  delete [] fgConvGamGeantIndex;
  delete [] feNegConvGamGeantIndex;
  delete [] fePosConvGamGeantIndex;

  delete [] feNegConvGamSingleRecIndex;
  delete [] fePosConvGamSingleRecIndex;

  delete [] feNegConvGamV0RecIndex;
  delete [] fePosConvGamV0RecIndex;
  delete [] fConvGamV0RecIndexPos;
  delete [] fConvGamV0RecIndexNeg;

  delete [] flDecayLGeantIndex;
  delete [] fpiNegDecayLGeantIndex;
  delete [] fpPosDecayLGeantIndex;

  delete [] fpiNegDecayLGeantLength;
  delete [] fpPosDecayLGeantLength;
  delete [] fpiNegDecayLSingleRecIndex;
  delete [] fpPosDecayLSingleRecIndex;

  delete [] fpiNegDecayLV0RecIndex;
  delete [] fpPosDecayLV0RecIndex;
  delete [] fDecayLV0RecIndexPos;
  delete [] fDecayLV0RecIndexNeg;

  delete [] falDecayALGeantIndex;
  delete [] fpiPosDecayALGeantIndex;
  delete [] fapNegDecayALGeantIndex;

  delete [] fpiPosDecayALGeantLength;
  delete [] fapNegDecayALGeantLength;
  delete [] fpiPosDecayALSingleRecIndex;
  delete [] fapNegDecayALSingleRecIndex;

  delete [] fpiPosDecayALV0RecIndex;
  delete [] fapNegDecayALV0RecIndex;
  delete [] fDecayALV0RecIndexPos;
  delete [] fDecayALV0RecIndexNeg;


  delete [] fpiNegDecayK0GeantIndex;
  delete [] fpiPosDecayK0GeantIndex;

  delete [] fpiNegDecayK0GeantLength;
  delete [] fpiPosDecayK0GeantLength;
  delete [] fpiNegDecayK0SingleRecIndex;
  delete [] fpiPosDecayK0SingleRecIndex;

  delete [] fpiNegDecayK0V0RecIndex;
  delete [] fpiPosDecayK0V0RecIndex;

  delete [] fDecayK0V0RecIndexPos;
  delete [] fDecayK0V0RecIndexNeg;

}


//________________________________________________________________________
void AliAnalysisTaskV0QA::UserCreateOutputObjects() {
  // Create Ouptut objects

  for(Int_t d=0;d<fDim;d++){
    fbinsV0[d]=70;
  }
  fxminV0[0]=   0;     // 1/sqrt(pt) Gamma geant
  fxmaxV0[0]=   8;


  fxminV0[1]=-2.5;   // eta Gamma Geant
  fxmaxV0[1]= 1.5;


  fxminV0[2]=-2*TMath::Pi();   // phi Gamma geant
  fxmaxV0[2]= TMath::Pi();


  fxminV0[3]=   0;     // r geant
  fxmaxV0[3]= 200;


  fxminV0[4]=-250;   // z geant
  fxmaxV0[4]= 250;


  fxminV0[5]=   0;     // 1/sqrt(pt) Geant Pos
  fxmaxV0[5]=   8;


  fxminV0[6]=-2.5;   // eta geant Pos
  fxmaxV0[6]= 1.5;


  fxminV0[7]=-2*TMath::Pi();   // phi Geant Pos
  fxmaxV0[7]= TMath::Pi();

  
  fxminV0[8]=0;   // Track Length TPC Geant Pos
  fxmaxV0[8]= 200;


  fxminV0[9]=   0;     // 1/sqrt(pt) Geant Neg
  fxmaxV0[9]=   8;


  fxminV0[10]=-2.5;   // eta Geant Neg
  fxmaxV0[10]= 1.5;


  fxminV0[11]=-2*TMath::Pi();   // phi Geant Neg
  fxmaxV0[11]= TMath::Pi();


  fxminV0[12]=0;   // Track Length TPC Geant Neg
  fxmaxV0[12]= 200;



  //-----------Rec single variables

  fxminV0[13]=   -0.5;     // (pt-ptGeant)/ptGeant rec Pos
  fxmaxV0[13]=   0.5;


  fxminV0[14]=-2.5;   // eta  rec Pos
  fxmaxV0[14]= 1.5;


  fxminV0[15]=-2*TMath::Pi();   // phi rec Pos
  fxmaxV0[15]= TMath::Pi();


  fxminV0[16]=   0;     // Impact parameter rec Pos
  fxmaxV0[16]=   100;


  fxminV0[17]=   0;     // nsigmas Impact parameter rec Pos
  fxmaxV0[17]=   100;



  fxminV0[18]=   -1;     // Ncls ITS rec Pos
  fxmaxV0[18]=   6;


  fxminV0[19]=   -1;     // Ncls TPC rec Pos
  fxmaxV0[19]=   180;


  fxminV0[20]=   -2;     // Status Single  TPC rec Pos
  fxmaxV0[20]=   2;


  fxminV0[21]=   -0.5;     // (pt-ptGeant)/ptGeant rec Neg
  fxmaxV0[21]=   0.5;


  fxminV0[22]=-2.5;   // eta  rec Neg
  fxmaxV0[22]= 1.5;


  fxminV0[23]=-2*TMath::Pi();   // phi rec Neg
  fxmaxV0[23]= TMath::Pi();


  fxminV0[24]=   0;     // Impact parameter rec Neg
  fxmaxV0[24]=   100;


  fxminV0[25]=   0;     // Sigmas Impact parameter rec Neg
  fxmaxV0[25]=   100;



  fxminV0[26]=   -1;     // Ncls ITS rec Neg
  fxmaxV0[26]=   6;


  fxminV0[27]=   -1;     // Ncls TPC rec Neg
  fxmaxV0[27]=   180;


  fxminV0[28]=   -2;     // Status Single  TPC rec Neg
  fxmaxV0[28]=   2;

  // ------------------Rec V0 variables



  fxminV0[29]=   -0.5;     // (pt-ptGeant)/ptGeant rec V0 Pos
  fxmaxV0[29]=   0.5;


  fxminV0[30]=-2.5;   // eta  rec V0 Pos
  fxmaxV0[30]= 1.5;


  fxminV0[31]=-2*TMath::Pi();   // phi rec V0 Pos
  fxmaxV0[31]= TMath::Pi();


  fxminV0[32]=   -2;     // Status V0 TPC rec Pos
  fxmaxV0[32]=   2;




  fxminV0[33]=   -0.5;     // 1/sqrt(pt) rec V0 Neg
  fxmaxV0[33]=   0.5;


  fxminV0[34]=-2.5;   // eta  rec V0 Neg
  fxmaxV0[34]= 1.5;


  fxminV0[35]=-2*TMath::Pi();   // phi rec V0 Neg
  fxmaxV0[35]= TMath::Pi();



  fxminV0[36]=   -2;     // Status V0 TPC rec Neg
  fxmaxV0[36]=   2;

  fxminV0[37]=   -1;     // Centrality
  fxmaxV0[37]=   10;   



  TString  axisName[38]={"ptGammaGeant", 
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
                         "statusV0SingleNeg",
	                 "centrality"};


  fSparseV0= new THnSparseF("sparseV0","sparseV0",fDim,fbinsV0,fxminV0,fxmaxV0);

  for (Int_t iaxis=0; iaxis<fDim; iaxis++){
   fSparseV0->GetAxis(iaxis)->SetName(axisName[iaxis]);
   fSparseV0->GetAxis(iaxis)->SetTitle(axisName[iaxis]);
  }

  TString  axisNameK0[38]={"ptK0Geant", 
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
                         "statusRecV0Neg",	                 
                         "centrality"};



  fSparseK0= new THnSparseF("sparseK0","sparseK0",fDim,fbinsV0,fxminV0,fxmaxV0);
  for (Int_t iaxis=0; iaxis<fDim; iaxis++){
   fSparseK0->GetAxis(iaxis)->SetName(axisNameK0[iaxis]);
   fSparseK0->GetAxis(iaxis)->SetTitle(axisNameK0[iaxis]);
  }

  TString  axisNameL[38]={"ptLGeant", 
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
                         "statusRecV0Neg",
	                 "centrality"};


  fSparseL= new THnSparseF("sparseL","sparseL",fDim,fbinsV0,fxminV0,fxmaxV0);
  for (Int_t iaxis=0; iaxis<fDim; iaxis++){
   fSparseL->GetAxis(iaxis)->SetName(axisNameL[iaxis]);
   fSparseL->GetAxis(iaxis)->SetTitle(axisNameL[iaxis]);
  }

  TString  axisNameAL[38]={"ptALGeant", 
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
                         "statusRecV0Neg",
	                 "centrality"};



  fSparseAL= new THnSparseF("sparseAL","sparseAL",fDim,fbinsV0,fxminV0,fxmaxV0);
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
void AliAnalysisTaskV0QA::UserExec(Option_t *) {
  // Execution of the Task

    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    
  if (!fESD) {
    //cout<< "not a tree"<< endl;
    return;
  }

  fnEv++;


  //Get MC data 
  fMCtruth = (AliMCEventHandler*)((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());

  //  Double_t vertex[3];
  Double_t maxVertex=150.;
  Double_t maxEta=0.9;
  Double_t lineCutZRSlope=tan(2*atan(exp(-maxEta)));
  Double_t lineCutZValue=7.;
  Int_t elecGIndex=-1;
  Int_t posiGIndex=-1;
  Int_t pPosLIndex=-1;
  Int_t piNegLIndex=-1;

  Int_t apNegALIndex=-1;
  Int_t piPosALIndex=-1;

  fnConvGamGeant=-1;
  fnDecayK0Geant=-1;
  fnDecayLGeant=-1;
  fnDecayALGeant=-1;

  for(Int_t i=0; i<fgDim;i++){
    fgConvGamGeantIndex[i] = -1;
    feNegConvGamGeantIndex[i] = -1;
    fePosConvGamGeantIndex[i] = -1;
    
    feNegConvGamSingleRecIndex[i] = -1;
    fePosConvGamSingleRecIndex[i] = -1;

    feNegConvGamV0RecIndex[i] = -1;
    fePosConvGamV0RecIndex[i] = -1;
    fConvGamV0RecIndexPos[i] = -1;
    fConvGamV0RecIndexNeg[i] = -1;


    fK0DecayK0GeantIndex[i] = -1;
    fpiNegDecayK0GeantIndex[i] = -1;
    fpiPosDecayK0GeantIndex[i] = -1;
    
    fpiNegDecayK0SingleRecIndex[i] = -1;
    fpiPosDecayK0SingleRecIndex[i] = -1;

    fpiNegDecayK0V0RecIndex[i] = -1;
    fpiPosDecayK0V0RecIndex[i] = -1;
    fDecayK0V0RecIndexPos[i] = -1;
    fDecayK0V0RecIndexNeg[i] = -1;


    flDecayLGeantIndex[i] = -1;
    fpiNegDecayLGeantIndex[i] = -1;
    fpPosDecayLGeantIndex[i] = -1;
    
    fpiNegDecayLSingleRecIndex[i] = -1;
    fpPosDecayLSingleRecIndex[i] = -1;

    fpiNegDecayLV0RecIndex[i] = -1;
    fpPosDecayLV0RecIndex[i] = -1;
    fDecayLV0RecIndexPos[i] = -1;
    fDecayLV0RecIndexNeg[i] = -1;

    // Antilambda
    falDecayALGeantIndex[i] = -1;
    fpiPosDecayALGeantIndex[i] = -1;
    fapNegDecayALGeantIndex[i] = -1;
    
    fpiPosDecayALSingleRecIndex[i] = -1;
    fapNegDecayALSingleRecIndex[i] = -1;

    fpiPosDecayALV0RecIndex[i] = -1;
    fapNegDecayALV0RecIndex[i] = -1;
    fDecayALV0RecIndexPos[i] = -1;
    fDecayALV0RecIndexNeg[i] = -1;


  }

  Int_t doMC=1;

  AliKFVertex primVtx(*(fESD->GetPrimaryVertex()));
  fnTracksPrim=primVtx.GetNContributors();


  if(fMCtruth && fnTracksPrim>0){

   fStack = fMCtruth->MCEvent()->Stack();


   if ( doMC){

    for (Int_t iTracks = 0; iTracks < fMCtruth->MCEvent()->GetNumberOfTracks(); iTracks++) {
      

     TParticle* particle = fStack->Particle(iTracks);



     if (!particle) {
       Printf("ERROR: Could not receive particle %d (mc loop)", iTracks);
       continue;
     }
     
     if(particle->Pt()<0.050) continue;
     if(TMath::Abs(particle->Eta())> maxEta) continue;

     
     if (particle->GetPdgCode()== 22){
       
       
       if(particle->GetMother(0) >-1 && fStack->Particle(particle->GetMother(0))->GetPdgCode() == 22){
	 continue; // no photon as mothers!
       }
       
       if(particle->GetMother(0) >= fStack->GetNprimary()){
	 continue; // the gamma has a mother, and it is not a primary particle
       }
       
       TParticle* ePos = NULL;
       TParticle* eNeg = NULL;
       elecGIndex=-1;
       posiGIndex=-1;
  
       if(particle->GetNDaughters() >= 2){
	 for(Int_t daughterIndex=particle->GetFirstDaughter();daughterIndex<=particle->GetLastDaughter();daughterIndex++){
	   TParticle *tmpDaughter = fStack->Particle(daughterIndex);
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
       
       if(TMath::Abs(ePos->Eta())> maxEta || TMath::Abs(eNeg->Eta())> maxEta){
	 continue;
       }	
       
       if(ePos->R()> maxVertex ){
	 continue; // cuts on distance from collision point
       }
      
      
       if( (TMath::Abs(ePos->Vz()) * lineCutZRSlope - lineCutZValue)  > ePos->R() ){
	 continue;               // line cut to exclude regions where we do not reconstruct
       }		

       
       // Looking at the existance of TPC references

       TParticle* ePosTPC;
       fMCtruth->MCEvent()->GetParticleAndTR(posiGIndex,ePosTPC,fclRefsP);

       AliMCParticle *mcParticlePos = (AliMCParticle*) (fMCtruth->MCEvent()->GetTrack(posiGIndex));
       if(!mcParticlePos) continue;

       Int_t counter; 
       Float_t tpcTrackLengthePos = mcParticlePos->GetTPCTrackLength(fESD->GetMagneticField(),0.05,counter,3.0); 

 

       int nPointsP =  fclRefsP->GetEntries();

       if (fRefTPC) delete fRefTPC;fRefTPC=NULL;
       fRefTPC = new TObjArray();

       for(int iPoint=0; iPoint<nPointsP; iPoint++) {
	 AliTrackReference *ref = (AliTrackReference*)fclRefsP->At(iPoint);
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
       fMCtruth->MCEvent()->GetParticleAndTR(elecGIndex,eNegTPC,fclRefsN);

       AliMCParticle *mcParticleNeg = (AliMCParticle*) (fMCtruth->MCEvent()->GetTrack(elecGIndex));
       if(!mcParticleNeg) continue;

       Int_t counterN; 
       Float_t tpcTrackLengtheNeg = mcParticleNeg->GetTPCTrackLength(fESD->GetMagneticField(),0.05,counterN,3.0); 
       int nPointsN =  fclRefsN->GetEntries();

       if (fRefTPC) delete fRefTPC; fRefTPC=NULL;
       fRefTPC = new TObjArray();

       for(int iPoint=0; iPoint<nPointsN; iPoint++) {
	 AliTrackReference *ref = (AliTrackReference*)fclRefsN->At(iPoint);
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
       
       
       fnConvGamGeant++;
       fgConvGamGeantIndex[fnConvGamGeant]=iTracks;
       feNegConvGamGeantIndex[fnConvGamGeant] = elecGIndex;
       fePosConvGamGeantIndex[fnConvGamGeant] = posiGIndex;
       
       feNegConvGamGeantLength[fnConvGamGeant] = tpcTrackLengtheNeg;
       fePosConvGamGeantLength[fnConvGamGeant] = tpcTrackLengthePos;

     }

     
     TParticle* piPos = NULL;
     TParticle* piNeg = NULL;
     fpiPosK0Index=-1;
     fpiNegK0Index=-1;

     if (particle->GetPdgCode()== 310){          // k0short
       if(particle->GetNDaughters() == 2){
	 for(Int_t daughterIndex=particle->GetFirstDaughter();daughterIndex<=particle->GetLastDaughter();daughterIndex++){
	   TParticle *tmpDaughter = fStack->Particle(daughterIndex);
	   if(tmpDaughter->GetPdgCode() == 211){
	     piPos= tmpDaughter;
	     fpiPosK0Index=daughterIndex;
	   }
	   else if(tmpDaughter->GetPdgCode() == -211){
	     piNeg = tmpDaughter;
	     fpiNegK0Index=daughterIndex;
	   }
	 }
       }

       if(piPos == NULL || piNeg == NULL){ // means we do not have two daughters from K0short decay
	 continue;
       }
       
       if(TMath::Abs(piPos->Eta())> maxEta || TMath::Abs(piNeg->Eta())> maxEta){
	 continue;
       }	
       
       if(piPos->R()> maxVertex ){
	 continue; // cuts on distance from collision point
       }
      
      
       if( (TMath::Abs(piPos->Vz()) * lineCutZRSlope - lineCutZValue)  > piPos->R() ){
	 continue;               // line cut to exclude regions where we do not reconstruct
       }		

      // Looking at the existance of TPC references

       TParticle* ePosTPC;
       fMCtruth->MCEvent()->GetParticleAndTR(fpiPosK0Index,ePosTPC,fclRefsP);

       AliMCParticle *mcParticlePos = (AliMCParticle*) (fMCtruth->MCEvent()->GetTrack(fpiPosK0Index));
       if(!mcParticlePos) continue;

       Int_t counter; 
       Float_t tpcTrackLengthePos = mcParticlePos->GetTPCTrackLength(fESD->GetMagneticField(),0.05,counter,3.0); 
 

       int nPointsP =  fclRefsP->GetEntries();
       if (fRefTPC) delete fRefTPC; fRefTPC=NULL;
       fRefTPC = new TObjArray();

       for(int iPoint=0; iPoint<nPointsP; iPoint++) {
	 AliTrackReference *ref = (AliTrackReference*)fclRefsP->At(iPoint);
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
	 if (aRef->GetTrack() != fpiPosK0Index ) break;			
	 ++labelPosRefs;
       }
     



       TParticle* eNegTPC;
       fMCtruth->MCEvent()->GetParticleAndTR(fpiNegK0Index,eNegTPC,fclRefsN);

       AliMCParticle *mcParticleNeg = (AliMCParticle*) (fMCtruth->MCEvent()->GetTrack(fpiNegK0Index));
       if(!mcParticleNeg) continue;

       Int_t counterN; 
       Float_t tpcTrackLengtheNeg = mcParticleNeg->GetTPCTrackLength(fESD->GetMagneticField(),0.05,counterN,3.0); 

       int nPointsN =  fclRefsN->GetEntries();
       if (fRefTPC) delete fRefTPC; fRefTPC=NULL;
       fRefTPC = new TObjArray();

       for(int iPoint=0; iPoint<nPointsN; iPoint++) {
	 AliTrackReference *ref = (AliTrackReference*)fclRefsN->At(iPoint);
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
	 if (aRef->GetTrack() != fpiNegK0Index ) break;			
	 ++labelNegRefs;
       }
       
       if ( labelNegRefs==0 || labelPosRefs==0) continue; // if pi+/pi- do not have a TPC ref continue;
       ////////////////////////////////////////////////////////////////////
       
       fnDecayK0Geant++;

       fK0DecayK0GeantIndex[fnDecayK0Geant]=iTracks;
       fpiNegDecayK0GeantIndex[fnDecayK0Geant]=fpiNegK0Index;
       fpiPosDecayK0GeantIndex[fnDecayK0Geant]=fpiPosK0Index;
       fpiNegDecayK0GeantLength[fnDecayK0Geant]=tpcTrackLengtheNeg;
       fpiPosDecayK0GeantLength[fnDecayK0Geant]=tpcTrackLengthePos;
      
     }    


     TParticle* pPos = NULL;
     TParticle* piNegL = NULL;
     pPosLIndex=-1;
     piNegLIndex=-1;

 
     if (particle->GetPdgCode()== 3122){        //lambda
      
       if(particle->GetNDaughters() == 2){
	 for(Int_t daughterIndex=particle->GetFirstDaughter();daughterIndex<=particle->GetLastDaughter();daughterIndex++){
	   TParticle *tmpDaughter = fStack->Particle(daughterIndex);
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

       if(TMath::Abs(pPos->Eta())> maxEta || TMath::Abs(piNegL->Eta())> maxEta){
	 continue;
       }	
       
       if(pPos->R()> maxVertex ){
	 continue; // cuts on distance from collision point
       }
      
      
       if( (TMath::Abs(pPos->Vz()) * lineCutZRSlope - lineCutZValue)  > pPos->R() ){
	 continue;               // line cut to exclude regions where we do not reconstruct
       }


     // Looking at the existance of TPC references

       TParticle* ePosTPC;
       fMCtruth->MCEvent()->GetParticleAndTR(pPosLIndex,ePosTPC,fclRefsP);

       AliMCParticle *mcParticlePos = (AliMCParticle*) (fMCtruth->MCEvent()->GetTrack(pPosLIndex));
       if(!mcParticlePos) continue;

       Int_t counter; 
       Float_t tpcTrackLengthePos = mcParticlePos->GetTPCTrackLength(fESD->GetMagneticField(),0.05,counter,3.0); 
 

       int nPointsP =  fclRefsP->GetEntries();
       if (fRefTPC) delete fRefTPC; fRefTPC=NULL;
       fRefTPC = new TObjArray();

       for(int iPoint=0; iPoint<nPointsP; iPoint++) {
	 AliTrackReference *ref = (AliTrackReference*)fclRefsP->At(iPoint);
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
       fMCtruth->MCEvent()->GetParticleAndTR(piNegLIndex,eNegTPC,fclRefsN);

       AliMCParticle *mcParticleNeg = (AliMCParticle*) (fMCtruth->MCEvent()->GetTrack(piNegLIndex));
       if(!mcParticleNeg) continue;

       Int_t counterN; 
       Float_t tpcTrackLengtheNeg = mcParticleNeg->GetTPCTrackLength(fESD->GetMagneticField(),0.05,counterN,3.0); 

       int nPointsN =  fclRefsN->GetEntries();
       if (fRefTPC) delete fRefTPC;
       fRefTPC = new TObjArray();

       for(int iPoint=0; iPoint<nPointsN; iPoint++) {
	 AliTrackReference *ref = (AliTrackReference*)fclRefsN->At(iPoint);
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
     
       fnDecayLGeant++;

       flDecayLGeantIndex[fnDecayLGeant]=iTracks;
     
       fpiNegDecayLGeantIndex[fnDecayLGeant]=piNegLIndex;
       fpPosDecayLGeantIndex[fnDecayLGeant]=pPosLIndex;

       fpiNegDecayLGeantLength[fnDecayLGeant]=tpcTrackLengtheNeg;
       fpPosDecayLGeantLength[fnDecayLGeant]=tpcTrackLengthePos;

     
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
	   TParticle *tmpDaughter = fStack->Particle(daughterIndex);
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

       if(TMath::Abs(apNeg->Eta())> maxEta || TMath::Abs(piPosAL->Eta())> maxEta){
	 continue;
       }	
       
       if(apNeg->R()> maxVertex ){
	 continue; // cuts on distance from collision point
       }
      
      
       if( (TMath::Abs(apNeg->Vz()) * lineCutZRSlope - lineCutZValue)  > apNeg->R() ){
	 continue;               // line cut to exclude regions where we do not reconstruct
       }


     // Looking at the existance of TPC references

       TParticle* ePosTPC;
       fMCtruth->MCEvent()->GetParticleAndTR(piPosALIndex,ePosTPC,fclRefsP);

       AliMCParticle *mcParticlePos = (AliMCParticle*) (fMCtruth->MCEvent()->GetTrack(piPosALIndex));
       if(!mcParticlePos) continue;

       Int_t counter; 
       Float_t tpcTrackLengthePos = mcParticlePos->GetTPCTrackLength(fESD->GetMagneticField(),0.05,counter,3.0); 

       int nPointsP =  fclRefsP->GetEntries();
       if (fRefTPC) delete fRefTPC;
       fRefTPC = new TObjArray();

       for(int iPoint=0; iPoint<nPointsP; iPoint++) {
	 AliTrackReference *ref = (AliTrackReference*)fclRefsP->At(iPoint);
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
       fMCtruth->MCEvent()->GetParticleAndTR(apNegALIndex,eNegTPC,fclRefsN);

       AliMCParticle *mcParticleNeg = (AliMCParticle*) (fMCtruth->MCEvent()->GetTrack(apNegALIndex));
       if(!mcParticleNeg) continue;

       Int_t counterN; 
       Float_t tpcTrackLengtheNeg = mcParticleNeg->GetTPCTrackLength(fESD->GetMagneticField(),0.05,counterN,3.0); 

       int nPointsN =  fclRefsN->GetEntries();
       if (fRefTPC) delete fRefTPC;
       fRefTPC = new TObjArray();

       for(int iPoint=0; iPoint<nPointsN; iPoint++) {
	 AliTrackReference *ref = (AliTrackReference*)fclRefsN->At(iPoint);
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
     
       fnDecayALGeant++;
       falDecayALGeantIndex[fnDecayALGeant]=iTracks;
     
       fpiPosDecayALGeantIndex[fnDecayALGeant]=piPosALIndex;
       fapNegDecayALGeantIndex[fnDecayALGeant]=apNegALIndex;

       fpiPosDecayALGeantLength[fnDecayALGeant]=tpcTrackLengthePos;
       fapNegDecayALGeantLength[fnDecayALGeant]=tpcTrackLengtheNeg;

     
     }  // AntiLambda    

    } //track loop 
   }    

  }


  AliKFParticle::SetField(fESD->GetMagneticField());

  const AliESDVertex *pvertex = fESD->GetPrimaryVertex();
  Double_t xyzVtx[3];
  pvertex->GetXYZ(xyzVtx);

  AliCentrality *esdCentrality = fESD->GetCentrality();
//  Int_t centralityC = -1;
  fCentralityC =-1;
  fCentralityC = esdCentrality->GetCentralityClass10("V0M");
 
  if(fnTracksPrim>0) {

    InspectListOfChargedParticles();
    InspectListOfV0s();


    if(fnConvGamGeant>-1){
      FillHnSparseGamma();
    }

    if(fnDecayK0Geant>-1){
      FillHnSparseK0();
    }

    if(fnDecayLGeant>-1){
      FillHnSparseL();
    }

    if(fnDecayALGeant>-1){
      FillHnSparseAL();
    }

  }

 
  PostData(1, fOutputContainer );
  

}

void AliAnalysisTaskV0QA::Terminate(Option_t *) {
  // Draw some histogram at the end.

}


Int_t AliAnalysisTaskV0QA::GetTPCReference(Int_t label) {
  // Get TPC References

  int start = TMath::BinarySearch(fRefTPC->GetEntries(), fLabelsTPC, label);

  while (start >= 0) {
    AliTrackReference *ref = (AliTrackReference*)(*fRefTPC)[start];
    if (ref->GetTrack() != label) return start+1;
    start--;
  }

  return 0;
}





void AliAnalysisTaskV0QA::InspectListOfChargedParticles(){
  // Look at the list of particles for the single track reconstruction

  for(Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++){

    AliESDtrack* curTrack = fESD->GetTrack(iTracks);

    if(!curTrack){
      continue;
    }


//     if( !(curTrack->GetStatus() & AliESDtrack::kTPCrefit)){
//       continue;
//     }


    Int_t labelMC = TMath::Abs(curTrack->GetLabel());

    if ( labelMC > fStack->GetNtrack() ) continue;


    TParticle* curParticle = fStack->Particle(labelMC);
    if(curParticle->GetMother(0)==-1){
      continue;
    }

     
    if(TMath::Abs(curParticle->GetPdgCode()) == 11){ // e+/e-
      
      if( fStack->Particle(curParticle->GetMother(0))->GetPdgCode()==22 ){  // e+/e- from gamma
	if( curParticle->GetUniqueID()!=5 ){ // e+/e- from gamma conversion
	  continue;
	}

	for(Int_t iGamConv=0;iGamConv<fnConvGamGeant+1;iGamConv++ ){
	  if(curTrack->GetSign()>0){
	    if (labelMC== fePosConvGamGeantIndex[iGamConv]){
	      fePosConvGamSingleRecIndex[iGamConv]=iTracks;
	    }
	  }else{
	    if (labelMC== feNegConvGamGeantIndex[iGamConv]){
	      feNegConvGamSingleRecIndex[iGamConv]=iTracks;
	    }
	  }
	} // loop over geant converted gammas

      }
    } // condition to select reconstructed electrons
  


    if(TMath::Abs(curParticle->GetPdgCode()) == 211 || TMath::Abs(curParticle->GetPdgCode())==2212 ){ // pi+/pi-
      
      if( fStack->Particle(curParticle->GetMother(0))->GetPdgCode()==310 || 
	  fStack->Particle(curParticle->GetMother(0))->GetPdgCode()==3122 || 
	  fStack->Particle(curParticle->GetMother(0))->GetPdgCode()==-3122 ){  // pi+/proton/pi- from K0/Lambda

	for(Int_t iK0Dec=0;iK0Dec<fnDecayK0Geant+1;iK0Dec++ ){
	  if(curTrack->GetSign()>0){
	    if (labelMC== fpiPosDecayK0GeantIndex[iK0Dec]){
	      fpiPosDecayK0SingleRecIndex[iK0Dec]=iTracks;
	    }
	  }else{
	    if (labelMC== fpiNegDecayK0GeantIndex[iK0Dec]){
	      fpiNegDecayK0SingleRecIndex[iK0Dec]=iTracks;
	    }
	  }
	} // loop over geant decay K0

	for(Int_t iLDec=0;iLDec<fnDecayLGeant+1;iLDec++ ){
	  if(curTrack->GetSign()>0){
	    if (labelMC== fpPosDecayLGeantIndex[iLDec]){
	      fpPosDecayLSingleRecIndex[iLDec]=iTracks;
	    }
	  }else{
	    if (labelMC== fpiNegDecayLGeantIndex[iLDec]){
	      fpiNegDecayLSingleRecIndex[iLDec]=iTracks;
	    }
	  }
	} // loop over geant decay Lambda

	for(Int_t iALDec=0;iALDec<fnDecayALGeant+1;iALDec++ ){
	  if(curTrack->GetSign()<0){
	    if (labelMC== fapNegDecayALGeantIndex[iALDec]){
	      fapNegDecayALSingleRecIndex[iALDec]=iTracks;
	    }
	  }else{
	    if (labelMC== fpiPosDecayALGeantIndex[iALDec]){
	      fpiPosDecayALSingleRecIndex[iALDec]=iTracks;
	    }
	  }
	} // loop over geant decay antiLambda
      }
    } // condition to select reconstructed electrons
   } // all reconstructed track


}

void AliAnalysisTaskV0QA::InspectListOfV0s(){
  // Look at the list of particles for the V0 reconstruction

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
    

    if ( !fV0MIs->GetOnFlyStatus() ){ //Onfly V0 finder
//    if ( fV0MIs->GetOnFlyStatus() ){   //Offline V0 finder
	continue;
    }
    if(fnTracksPrim<=0) {
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

    if(!trackPos) return;
    if(!trackNeg) return;

    Int_t labelNeg=TMath::Abs(trackNeg->GetLabel());
    if(labelNeg > fStack->GetNtrack() ) continue;
    TParticle * particleNeg= fStack->Particle(labelNeg);

    Int_t labelPos=TMath::Abs(trackPos->GetLabel());
    if(labelPos > fStack->GetNtrack() ) continue;
    TParticle * particlePos= fStack->Particle(labelPos);


    if(particlePos->GetMother(0)>-1){
      grandMotherPos=fStack->Particle(particlePos->GetMother(0))->GetMother(0);
      motherPos=particlePos->GetMother(0);
    }
    
    if(particleNeg->GetMother(0)>-1){
      grandMotherNeg=fStack->Particle(particleNeg->GetMother(0))->GetMother(0);
      motherNeg=particleNeg->GetMother(0);
    }

    if(motherPos == motherNeg &&  motherPos!=-1 ){
      if( particlePos->GetPdgCode() ==-11  &&   particleNeg->GetPdgCode()==11 ){
	for(Int_t iGamConv=0;iGamConv<fnConvGamGeant+1;iGamConv++ ){
	  if (labelPos== fePosConvGamGeantIndex[iGamConv]){
	    fePosConvGamV0RecIndex[iGamConv]=pIndex;
	    fConvGamV0RecIndexPos[iGamConv]=iV0MI;
	  }
	  if (labelNeg== feNegConvGamGeantIndex[iGamConv]){
	    feNegConvGamV0RecIndex[iGamConv]=nIndex;
	    fConvGamV0RecIndexNeg[iGamConv]=iV0MI;
	  }

	} // loop over geant converted gammas
      }

      if( particlePos->GetPdgCode()==211  &&   particleNeg->GetPdgCode()==-211 ){
	for(Int_t iK0Dec=0;iK0Dec<fnDecayK0Geant+1;iK0Dec++ ){
	  if (labelPos== fpiPosDecayK0GeantIndex[iK0Dec]){
	    fpiPosDecayK0V0RecIndex[iK0Dec]=pIndex;
	    fDecayK0V0RecIndexPos[iK0Dec]=iV0MI;
	  }
	  if (labelNeg== fpiNegDecayK0GeantIndex[iK0Dec]){
	    fpiNegDecayK0V0RecIndex[iK0Dec]=nIndex;
	    fDecayK0V0RecIndexNeg[iK0Dec]=iV0MI;
	  }

	} // loop over geant K0
      }

      if( particlePos->GetPdgCode()==2212  && particleNeg->GetPdgCode()==-211 ){
	for(Int_t iLDec=0;iLDec<fnDecayLGeant+1;iLDec++ ){
	  if (labelPos== fpPosDecayLGeantIndex[iLDec]){
	    fpPosDecayLV0RecIndex[iLDec]=pIndex;
	    fDecayLV0RecIndexPos[iLDec]=iV0MI;
	  }
	  if (labelNeg== fpiNegDecayLGeantIndex[iLDec]){
	    fpiNegDecayLV0RecIndex[iLDec]=nIndex;
	    fDecayLV0RecIndexNeg[iLDec]=iV0MI;
	  }

	} // loop over geant Lambda
      }

      if( particleNeg->GetPdgCode()==-2212  && particlePos->GetPdgCode()==211 ){
	for(Int_t iALDec=0;iALDec<fnDecayALGeant+1;iALDec++ ){
	  if (labelNeg== fapNegDecayALGeantIndex[iALDec]){
	    fapNegDecayALV0RecIndex[iALDec]=nIndex;
	    fDecayALV0RecIndexNeg[iALDec]=iV0MI;
	  }
	  if (labelPos== fpiPosDecayALGeantIndex[iALDec]){
	    fpiPosDecayALV0RecIndex[iALDec]=pIndex;
	    fDecayALV0RecIndexPos[iALDec]=iV0MI;
	  }

	} // loop over geant antiLambda
      }


    }
    
  }
  for(Int_t iGamConv=0;iGamConv<fnConvGamGeant+1;iGamConv++ ){
    if ( fConvGamV0RecIndexNeg[iGamConv]!=  fConvGamV0RecIndexPos[iGamConv]){
      fePosConvGamV0RecIndex[iGamConv]=-1;
      feNegConvGamV0RecIndex[iGamConv]=-1; 
      fConvGamV0RecIndexNeg[iGamConv]=-1;
      fConvGamV0RecIndexPos[iGamConv]=-1;

    }
  }

  for(Int_t iLDec=0;iLDec<fnDecayLGeant+1;iLDec++ ){
    if(fDecayLV0RecIndexPos[iLDec] !=  fDecayLV0RecIndexNeg[iLDec]){
      fpiNegDecayLV0RecIndex[iLDec]=-1;
      fpPosDecayLV0RecIndex[iLDec]=-1;
      fDecayLV0RecIndexNeg[iLDec]=-1;
      fDecayLV0RecIndexPos[iLDec]=-1;
    }
  }

  for(Int_t iALDec=0;iALDec<fnDecayALGeant+1;iALDec++ ){
    if(fDecayALV0RecIndexPos[iALDec] !=  fDecayALV0RecIndexNeg[iALDec]){
      fpiPosDecayALV0RecIndex[iALDec]=-1;
      fapNegDecayALV0RecIndex[iALDec]=-1;
      fDecayALV0RecIndexNeg[iALDec]=-1;
      fDecayALV0RecIndexPos[iALDec]=-1;
    }
  }

  for(Int_t iK0Dec=0;iK0Dec<fnDecayK0Geant+1;iK0Dec++ ){
    if(fDecayK0V0RecIndexPos[iK0Dec] !=  fDecayK0V0RecIndexNeg[iK0Dec]){
      fpiNegDecayK0V0RecIndex[iK0Dec]=-1;
      fpiPosDecayK0V0RecIndex[iK0Dec]=-1;
      fDecayK0V0RecIndexNeg[iK0Dec]=-1;
      fDecayK0V0RecIndexPos[iK0Dec]=-1;
    }
  }
  

}
void AliAnalysisTaskV0QA::FillHnSparseGamma()
{
  // Fill THnSparse Gamma

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


  for(Int_t i=0;i<fnConvGamGeant+1;i++){
    TParticle* gamPart = fStack->Particle(fgConvGamGeantIndex[i]);
    TParticle* ePosPart = fStack->Particle(fePosConvGamGeantIndex[i]);
    TParticle* eNegPart = fStack->Particle(feNegConvGamGeantIndex[i]);
    if (fePosConvGamSingleRecIndex[i]!=-1){
      AliESDtrack * ePosSglTrack = fESD->GetTrack(fePosConvGamSingleRecIndex[i]);
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
    
    if (feNegConvGamSingleRecIndex[i]!=-1){
      AliESDtrack * eNegSglTrack = fESD->GetTrack(feNegConvGamSingleRecIndex[i]);
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
    
    if(fConvGamV0RecIndexPos[i]!=-1){
      AliESDv0 * fV0MIs = fESD->GetV0(fConvGamV0RecIndexPos[i]);
      AliESDtrack* trackPosTest = fESD->GetTrack(fV0MIs->GetPindex());
      AliESDtrack* trackNegTest = fESD->GetTrack(fV0MIs->GetNindex());

      if (fePosConvGamV0RecIndex[i]!=-1 ){
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
      
      if (feNegConvGamV0RecIndex[i]!=-1 ){
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
     fValueV0[8] = fePosConvGamGeantLength[i];

    fValueV0[9] = 1./TMath::Sqrt(eNegPart->Pt());
    fValueV0[10] = eNegPart->Eta();

    Double_t tmpNPhi=eNegPart->Phi();
    if( eNegPart->Phi()>TMath::Pi()){
      tmpNPhi = eNegPart->Phi()-2*TMath::Pi();
    }
    fValueV0[11] = tmpNPhi;
    fValueV0[12] = feNegConvGamGeantLength[i];    

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
    fValueV0[37] = fCentralityC;
    fSparseV0->Fill(fValueV0);
  }


}

void AliAnalysisTaskV0QA::FillHnSparseK0()
{
  // Fill THnSparse K0

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
  
  for(Int_t i=0;i<fnDecayK0Geant+1;i++){
    TParticle* k0Part = fStack->Particle(fK0DecayK0GeantIndex[i]);
    TParticle* ePosPart = fStack->Particle(fpiPosDecayK0GeantIndex[i]);
    TParticle* eNegPart = fStack->Particle(fpiNegDecayK0GeantIndex[i]);
    if (fpiPosDecayK0SingleRecIndex[i]!=-1){
      AliESDtrack * ePosSglTrack = fESD->GetTrack(fpiPosDecayK0SingleRecIndex[i]);
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

    if (fpiNegDecayK0SingleRecIndex[i]!=-1){
      AliESDtrack * eNegSglTrack = fESD->GetTrack(fpiNegDecayK0SingleRecIndex[i]);
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

    if(fDecayK0V0RecIndexPos[i]!=-1){
      AliESDv0 * fV0MIs = fESD->GetV0(fDecayK0V0RecIndexPos[i]);
      AliESDtrack* trackPosTest = fESD->GetTrack(fV0MIs->GetPindex());
      AliESDtrack* trackNegTest = fESD->GetTrack(fV0MIs->GetNindex());

      if (fpiPosDecayK0V0RecIndex[i]!=-1 ){
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
      
      if (fpiNegDecayK0V0RecIndex[i]!=-1 ){
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
    fValueK0[0] = 1./TMath::Sqrt(k0Part->Pt());
    fValueK0[1] = k0Part->Eta();

    Double_t tmpGPhi=k0Part->Phi();
    if( k0Part->Phi()>TMath::Pi()){
      tmpGPhi=k0Part->Phi()-2*TMath::Pi();
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
     fValueK0[8] = fpiPosDecayK0GeantLength[i];

    fValueK0[9] = 1./TMath::Sqrt(eNegPart->Pt());
    fValueK0[10] = eNegPart->Eta();

    Double_t tmpNPhi=eNegPart->Phi();
    if( eNegPart->Phi()>TMath::Pi()){
      tmpNPhi = eNegPart->Phi()-2*TMath::Pi();
    }
    fValueK0[11] = tmpNPhi;
    fValueK0[12] = fpiNegDecayK0GeantLength[i];    
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
    fValueK0[37] = fCentralityC;

    fSparseK0->Fill(fValueK0);
  }
  

}
void AliAnalysisTaskV0QA::FillHnSparseL()
{
  // Fill THnSparse Lambda

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

  for(Int_t i=0;i<fnDecayLGeant+1;i++){
    TParticle* lPart = fStack->Particle(flDecayLGeantIndex[i]);
    TParticle* ePosPart = fStack->Particle(fpPosDecayLGeantIndex[i]);
    TParticle* eNegPart = fStack->Particle(fpiNegDecayLGeantIndex[i]);
    if (fpPosDecayLSingleRecIndex[i]!=-1){
      AliESDtrack * ePosSglTrack = fESD->GetTrack(fpPosDecayLSingleRecIndex[i]);
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
    
    if (fpiNegDecayLSingleRecIndex[i]!=-1){
      AliESDtrack * eNegSglTrack = fESD->GetTrack(fpiNegDecayLSingleRecIndex[i]);
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
    
    if(fDecayLV0RecIndexPos[i]!=-1){
      AliESDv0 * fV0MIs = fESD->GetV0(fDecayLV0RecIndexPos[i]);
      AliESDtrack* trackPosTest = fESD->GetTrack(fV0MIs->GetPindex());
      AliESDtrack* trackNegTest = fESD->GetTrack(fV0MIs->GetNindex());

      if (fpPosDecayLV0RecIndex[i]!=-1 ){
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
      
      if (fpiNegDecayLV0RecIndex[i]!=-1 ){
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
     fValueL[8] = fpPosDecayLGeantLength[i];

    fValueL[9] = 1./TMath::Sqrt(eNegPart->Pt());
    fValueL[10] = eNegPart->Eta();

    Double_t tmpNPhi=eNegPart->Phi();
    if( eNegPart->Phi()>TMath::Pi()){
      tmpNPhi = eNegPart->Phi()-2*TMath::Pi();
    }
    fValueL[11] = tmpNPhi;
    fValueL[12] = fpiNegDecayLGeantLength[i];    
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
    fValueL[37] = fCentralityC;

    fSparseL->Fill(fValueL);
  }


}

void AliAnalysisTaskV0QA::FillHnSparseAL()
{
  // Fill THnSparse Antilambda

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


  for(Int_t i=0;i<fnDecayALGeant+1;i++){
    TParticle* alPart = fStack->Particle(falDecayALGeantIndex[i]);
    TParticle* eNegPart = fStack->Particle(fapNegDecayALGeantIndex[i]);
    TParticle* ePosPart = fStack->Particle(fpiPosDecayALGeantIndex[i]);
    if (fpiPosDecayALSingleRecIndex[i]!=-1){
      AliESDtrack * ePosSglTrack = fESD->GetTrack(fpiPosDecayALSingleRecIndex[i]);
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
    
    if (fapNegDecayALSingleRecIndex[i]!=-1){
      AliESDtrack * eNegSglTrack = fESD->GetTrack(fapNegDecayALSingleRecIndex[i]);
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
    
    if(fDecayALV0RecIndexPos[i]!=-1){
      AliESDv0 * fV0MIs = fESD->GetV0(fDecayALV0RecIndexPos[i]);
      AliESDtrack* trackPosTest = fESD->GetTrack(fV0MIs->GetPindex());
      AliESDtrack* trackNegTest = fESD->GetTrack(fV0MIs->GetNindex());

      if (fpiPosDecayALV0RecIndex[i]!=-1 ){
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
      
      if (fapNegDecayALV0RecIndex[i]!=-1 ){
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
     fValueAL[8] = fpiPosDecayALGeantLength[i];

    fValueAL[9] = 1./TMath::Sqrt(eNegPart->Pt());
    fValueAL[10] = eNegPart->Eta();

    Double_t tmpNPhi=eNegPart->Phi();
    if( eNegPart->Phi()>TMath::Pi()){
      tmpNPhi = eNegPart->Phi()-2*TMath::Pi();
    }
    fValueAL[11] = tmpNPhi;
    fValueAL[12] = fapNegDecayALGeantLength[i];    
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
    fValueAL[37] = fCentralityC;

    fSparseAL->Fill(fValueAL);
  }
}


// void AliAnalysisTaskV0QA::SetESDtrackCuts()
// {

//   fEsdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");

//   fEsdTrackCuts->SetRequireTPCRefit(kTRUE);
//   fEsdTrackCuts->SetAcceptKinkDaughters(kFALSE);



// }
