/***************************************************************          *
 * Authors : Simone Schuchmann 
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//-----------------------------------------------------------------
// AliAnalysisTaskV0ForRAA class
// This task is for analysing Lambda and K0s pt spectra in PbPb and
// pp as well as with MC. The flag for pp and MC  must be set
// accordingly, default is PbPb data.
// It works with ESD files only.
//-----------------------------------------------------------------


#define AliAnalysisTaskV0ForRAA_cxx

#include "AliAnalysisTaskV0ForRAA.h"

#include "Riostream.h"

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"

#include "AliESDVZERO.h"
#include "AliMultiplicity.h"
#include "AliCentrality.h"

#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliESDtrackCuts.h"
#include "AliESDpid.h"
#include "AliESDv0.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"



ClassImp(AliAnalysisTaskV0ForRAA)

//________________________________________________________________________
AliAnalysisTaskV0ForRAA::AliAnalysisTaskV0ForRAA(const char *name)
:AliAnalysisTaskSE(name),
  fESD(0),
  fMCev(0),
//other objects
  fESDpid(0),
  fESDTrackCuts(0),
  fESDTrackCutsCharged(0),
  fOutputContainer(0),
//event histos
  fHistNEvents(0),  
  fHistMuliplicity(0),
  fHistMuliplicityRaw(0),
  fHistMultiplicityPrimary(0),
  fHistESDVertexZ(0),
  fHistPrimVtxZESD(0),
  fHistPrimVtxZESDVSNContributors(0),
  fHistNPrim(0),
  fHistITSLayerHits(0),
  fHistOneHitWithSDD(0),
  fHistPrimVtxZESDTPC(0),
  fHistPrimVtxZESDSPD(0),  
  fHistPrimVtxZESDTPCVSNContributors(0),
  fHistPrimVtxZESDSPDVSNContributors(0),
  fHistCentBinRaw(0),
  fHistCentBin(0),
  fHistPrimVtxZESDVSNContributorsMC(0),
  fHistPrimVtxZESDTPCVSNContributorsMC(0),
  fHistPrimVtxZESDSPDVSNContributorsMC(0),
  fHistMCVertexZ(0),
//MC pdg code histos
  fHistPiPiPDGCode(0),
  fHistPiPPDGCode(0),
  fHistPiAPPDGCode(0),
//cosine of pointing angle of Xi vs pt histos
  fHistPiPCosPointAngXiVsPt(0),
  fHistPiAPCosPointAngXiVsPt(0),
// fHistUserPtShift(0),
//selection booleans and values
  fMCMode(0),
  fMCTruthMode(0),
  fSelectInjected(0),
  fUseCentrality(0),
  fUseCentralityBin(0),
  fUseCentralityRange(0),
  fAnapp(0),
  fSelSDD(0),
  fOntheFly(0),
  fVertexZCut(0),
  fVtxStatus(0),
  fUsePID(0),
  fNSigma(0),
  fPPIDcut(0),
  fMoreNclsThanRows(0),
  fMoreNclsThanFindable(0),
  fChi2PerClusterITS(0),
  fRapCutV0(0),
  fRap(0),
  fEtaCutMCDaughters(0),
  fEtaCutMCDaughtersVal(0),
  fMinPt(0),
  fAlfaCut(0),
  fQtCut(0),
  fQtCutPtLow(0),
  fArmCutK0(0),      
  fArmCutL(0),  
  fDCAToVertexK0(0),
  fDCAToVertexL(0),
  fDCAXK(0),
  fDCAYK(0),
  fDCAXL(0),
  fDCAYL(0),
  fDCAZ(0),
  fDCADaughtersL(0),
  fDCADaughtersAL(0),
  fDCADaughtersK0(0),
  fDCADaughtersToVtxLarge(0),
  fDCADaughtersToVtxSmall(0),
  fDecayRadXYMin(0),
  fDecayRadXYMax(0),
  fCosPointAngL(0),
  fCosPointAngK(0),
  fCPAPtCutK0(0),
  fCPAPtCutL(0),
  fOpengAngleDaughters(0),
  fOpAngPtCut(0),
  fDecayLengthMax(0),
  fDecayLengthMin(0),
  fCtauK0s(0),
  fCtauL(0),
  fCtauPtCutK0(0),
  fCtauPtCutL(0),
  fChiCutKf(0)				       //  fShift(0),
  // fDeltaInvP(0)
{  // Constructor.

  DefineOutput(1,TList::Class());

  // define defaults for globals
  /*
    fShift = kFALSE;                       // shift in charge/pt yes/no
    fDeltaInvP = 0.00;                     // shift value
  */
   
  fMCMode = kFALSE;
  fMCTruthMode = kFALSE;

  fUseCentrality = 0;
  fUseCentralityBin = 0;
  fUseCentralityRange =0;

  fAnapp = kFALSE;
  fSelSDD = kFALSE;

  fSelectInjected = kFALSE;

  fVertexZCut = 100000.0;
  fVtxStatus = kFALSE;
   
  fOntheFly = kTRUE;

  //----- define defaults for V0 and track cuts ----//
  fUsePID = kFALSE;
  fMoreNclsThanRows = kFALSE;
  fMoreNclsThanFindable = kFALSE;
  fChi2PerClusterITS = 100000.0;
  fNSigma = 100000.0;
  fPPIDcut = 100.0;

  fRapCutV0=kFALSE;
  fRap=1000.0;
  fRap=1000.0;

  fAlfaCut= -100.0;
  fQtCut = 100.0;
  fQtCutPtLow = -1.0;
  fArmCutK0=kFALSE;     
  fArmCutL=kFALSE;  

  fEtaCutMCDaughters = kFALSE;
  fEtaCutMCDaughtersVal = 50.0;

  fMinPt= -1.0;

  fDCAToVertexK0 = 10000.0;
  fDCAToVertexL = 10000.0;
  fDCAXK=10000.0;
  fDCAYK=10000.0;
  fDCAXL=10000.0;
  fDCAYL=10000.0;
  fDCAZ=10000.0;
   
  fDCADaughtersL=10000.0;
  fDCADaughtersAL=10000.0;
  fDCADaughtersK0=10000.0;

  fDCADaughtersToVtxLarge=-1.0;
  fDCADaughtersToVtxSmall=-1.0;

  fDecayRadXYMin=-1.0;
  fDecayRadXYMax=1000000.0;
  fDecayLengthMax=100000.0;
  fDecayLengthMin=-1.0;
   
  fCosPointAngL=-1.0;
  fCosPointAngK=-1.0;
  fCPAPtCutK0 = 1000.0;
  fCPAPtCutL =1000.0;
  fOpengAngleDaughters = -1.0;
  fOpAngPtCut = -1.0;
      
  fCtauK0s=10000.0;
  fCtauL=10000.0;
  fCtauPtCutK0=1000.0;
  fCtauPtCutL=10000.0;

  fChiCutKf=1000000.0;

   
  //---- histograms ----//
  for(Int_t j=0;j<2;j++){
    fHistV0RadiusZ[j] =NULL;
    fHistV0RadiusZVSPt[j] =NULL;
    fHistV0RadiusXY[j] =NULL;
    fHistV0RadiusXYVSY[j] =NULL;
    fHistArmenteros[j]=NULL;
 
    //K0
    fHistPiPiMass[j]=NULL;
    fHistPiPiMassVSPt[j]=NULL;
    fHistPiPiMassVSPtMCTruth[j]=NULL;
    fHistPiPiPtVSY[j]=NULL;
    // fHistPiPiMassVSAlpha[j]=NULL;
    fHistPiPiDecayLengthVsPt[j]=NULL;
    fHistPiPiMonitorCuts[j] =NULL;
      
    //Lambda
    fHistPiPMass[j]=NULL;
    fHistPiPMassVSPt[j]=NULL;
    fHistPiPMassVSPtMCTruth[j]=NULL;
    fHistPiPPtVSY[j]=NULL;
    fHistPiPDecayLengthVsPt[j]=NULL;
    fHistPiPMonitorCuts[j] =NULL;
    //ALambda
    fHistPiAPMass[j]=NULL;
    fHistPiAPMassVSPt[j]=NULL;
    fHistPiAPMassVSPtMCTruth[j]=NULL;
    fHistPiAPPtVSY[j]=NULL;
    fHistPiAPDecayLengthVsPt[j]=NULL;
    fHistPiAPMonitorCuts[j] =NULL;

    //eta daughters reco
    fHistPiPiEtaDReco[j] = NULL;
    fHistPiPEtaDReco[j] = NULL;      

    //K0s
    fHistPiPiDCADaughters[j]=NULL; 
    fHistPiPiDCADaughterPosToPrimVtxVSMass[j]=NULL;
    // fHistPiPiMassVSPtK0L[j]=NULL;
    fHistPiPiDCAVSMass[j]=NULL;
    fHistPiPiCosPointAng[j]=NULL;
    fHistPiPiDecayLengthVsMass[j]=NULL;
    fHistPiPiRadiusXY[j]=NULL;


    //Lambda
    fHistPiPDCADaughters[j]=NULL;
    fHistPiPDCADaughterPosToPrimVtxVSMass[j]=NULL;
    fHistPiPDCAVSMass[j]=NULL;
    fHistPiPCosPointAng[j]=NULL;
    fHistPiPDecayLengthVsMass[j]=NULL;
    fHistPiPRadiusXY[j]=NULL;   
    //ALambda
    fHistPiAPDCADaughters[j]=NULL;
    fHistPiAPDCADaughterPosToPrimVtxVSMass[j]=NULL;
    fHistPiAPDCAVSMass[j]=NULL;
    fHistPiAPCosPointAng[j]=NULL;
    fHistPiAPDecayLengthVsMass[j]=NULL;
    fHistPiAPRadiusXY[j]=NULL;

    //other 
    fHistDedxSecProt[j]=NULL;
    fHistDedxSecPiPlus[j]=NULL;
    fHistDedxSecAProt[j]=NULL;
    fHistDedxSecPiMinus[j]=NULL;

    fHistNclsITSPosK0[j]=NULL;
    fHistNclsITSNegK0[j]=NULL;
    fHistNclsTPCPosK0[j]=NULL;
    fHistNclsTPCNegK0[j]=NULL;
    fHistChi2PerNclsITSPosK0[j]=NULL;
    fHistChi2PerNclsITSNegK0[j]=NULL;

    fHistNclsITSPosL[j]=NULL;
    fHistNclsITSNegL[j]=NULL;
    fHistNclsTPCPosL[j]=NULL;
    fHistNclsTPCNegL[j]=NULL;
    fHistChi2PerNclsITSPosL[j]=NULL;
    fHistChi2PerNclsITSNegL[j]=NULL;

    fHistNclsITSPos[j]=NULL;
    fHistNclsITSNeg[j]=NULL;
    fHistNclsTPCPos[j]=NULL;
    fHistNclsTPCNeg[j]=NULL;
    fHistChi2PerNclsITSPos[j]=NULL;
    fHistChi2PerNclsITSNeg[j]=NULL;
    fHistNclsITS[j]=NULL;
    fHistNclsTPC[j]=NULL;
    fHistNCRowsTPCPos[j]=NULL;
    fHistNCRowsTPCNeg[j]=NULL;
      
    // eta daughters MC
    fHistPiPiEtaDMC[j] = NULL;
    fHistPiPEtaDMC[j] = NULL;      

    //Lambda sec      
    fHistPiPMassVSPtSecSigma[j]=NULL;
    fHistPiPMassVSPtSecXi[j]=NULL;
    fHistPiPMassVSYSecXi[j]=NULL;
    fHistPiPXi0PtVSLambdaPt[j]=NULL;
    fHistPiPXiMinusPtVSLambdaPt[j]=NULL;
      
    //ALambda sec
    fHistPiAPMassVSPtSecSigma[j]=NULL;
    fHistPiAPMassVSPtSecXi[j]=NULL;
    fHistPiAPMassVSYSecXi[j]=NULL;
    fHistPiAPXi0PtVSLambdaPt[j]=NULL;
    fHistPiAPXiMinusPtVSLambdaPt[j]=NULL;

  }
   
}
//_____________________________________________________
AliAnalysisTaskV0ForRAA::~AliAnalysisTaskV0ForRAA()
{
  //---- Remove all pointers ----//
  if(fOutputContainer) delete fOutputContainer;fOutputContainer=0;
  if(fESDTrackCuts) delete fESDTrackCuts;fESDTrackCuts=0;
  if(fESDTrackCutsCharged) delete fESDTrackCutsCharged;fESDTrackCutsCharged=0;

}
//________________________________________________________________________
void AliAnalysisTaskV0ForRAA::UserCreateOutputObjects(){
  //create output objects

  Int_t nbPt=800;
  Int_t nbMass=500;

   
  //-----------------  create output container -----------------//

  fOutputContainer = new TList() ;
  fOutputContainer->SetName(GetName()) ;
  fOutputContainer->SetOwner();
 

  //------------ create allways -----------------------//
  fHistNEvents = new TH1F("fHistNEvents","no of events before cuts =0, after cuts=1, after process =2",5,0.0,5.0);
  fOutputContainer->Add(fHistNEvents);
      
  fHistMuliplicity =  new TH1F("fHistMuliplicity","V0 multiplicity",3000,0.0,30000);
  fOutputContainer->Add(fHistMuliplicity);
      
  fHistMuliplicityRaw =  new TH1F("fHistMuliplicityRaw","V0 multiplicity before process",3000,0.0,30000);      
  fOutputContainer->Add(fHistMuliplicityRaw);
      
  fHistMultiplicityPrimary = new TH1F("fHistMultiplicityPrimary","number of charged tracks",5000,0.0,20000);
  fOutputContainer->Add(fHistMultiplicityPrimary);
      
  fHistESDVertexZ= new TH1F("fHistESDVertexZ"," z vertex distr in cm",500,-50,50);
  fOutputContainer->Add(fHistESDVertexZ);
   
  fHistPrimVtxZESD = new TH1F("fHistPrimVtxZESD","z vertex pos ESD",250,-50,50);
  fOutputContainer->Add(fHistPrimVtxZESD);
      
  fHistPrimVtxZESDVSNContributors = new TH2F("fHistPrimVtxZESDVSNContributors","prim vtx pos z ESD vs no. of contributers TPC",250,-50,50,500,0.0,500.0);
  fOutputContainer->Add(fHistPrimVtxZESDVSNContributors);
      
  fHistNPrim = new TH1F("fHistNPrim","Number of contributers to vertex",2500,0.0,5000);
  fOutputContainer->Add(fHistNPrim);
 
  //------------------------ pp analysis only -------------------------//
  if(fAnapp){
    fHistITSLayerHits = new TH1F("fHistITSLayerHits","SDD layer -1=0,1=1,2=2 ... 5=5,0=nothing",7,-1.5,5.5);
    fOutputContainer->Add(fHistITSLayerHits);
    fHistOneHitWithSDD = new TH1F("fHistOneHitWithSDD","min one hit in SDD",2,-0.5,1.5);
    fOutputContainer->Add(fHistOneHitWithSDD);
    fHistPrimVtxZESDTPC = new TH1F("fHistPrimVtxZESDTPC","z vertex pos TPC",250,-50,50);
    fOutputContainer->Add(fHistPrimVtxZESDTPC);
    fHistPrimVtxZESDSPD = new TH1F("fHistPrimVtxZESDSPD","z vertex pos SPD",250,-50,50);
    fOutputContainer->Add(fHistPrimVtxZESDSPD);  
    fHistPrimVtxZESDTPCVSNContributors = new TH2F("fHistPrimVtxZESDTPCVSNContributors","prim vtx pos z TPC vs no. of contributers TPC",250,-50,50,500,0.0,500.0);
    fOutputContainer->Add(fHistPrimVtxZESDTPCVSNContributors);
    fHistPrimVtxZESDSPDVSNContributors = new TH2F("fHistPrimVtxZESDSPDVSNContributors","prim vtx pos z SPD vs no. of contributers TPC",250,-50,50,500,0.0,500.0);
    fOutputContainer->Add(fHistPrimVtxZESDSPDVSNContributors);

  }
  else {
    Double_t binsCent[12]={0.0,5.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0};
    fHistCentBinRaw = new TH1F("fHistCentBinRaw","centrality bin before cent selection",11,binsCent);
    fOutputContainer->Add(fHistCentBinRaw);
    fHistCentBin = new TH1F("fHistCentBin","centrality bin",11,binsCent);
    fOutputContainer->Add(fHistCentBin);
      
  }
   

   
  fHistPrimVtxZESDVSNContributorsMC = new TH2F("fHistPrimVtxZESDVSNContributorsMC","prim vtx pos z ESD vs no. of contributers MC",250,-50,50,500,0.0,500.0);
  fOutputContainer->Add(fHistPrimVtxZESDVSNContributorsMC);
  fHistPrimVtxZESDTPCVSNContributorsMC = new TH2F("fHistPrimVtxZESDTPCVSNContributorsMC","prim vtx pos z TPC vs no. of contributers MC",250,-50,50,500,0.0,500.0);
  fOutputContainer->Add(fHistPrimVtxZESDTPCVSNContributorsMC);
  fHistPrimVtxZESDSPDVSNContributorsMC = new TH2F("fHistPrimVtxZESDSPDVSNContributorsMC","prim vtx pos z SPD vs no. of contributers MC",250,-50,50,500,0.0,500.0);
  fOutputContainer->Add(fHistPrimVtxZESDSPDVSNContributorsMC);
  fHistMCVertexZ= new TH1F("fHistMCVertexZ"," z vertex distr in cm MC",500,-50,50);
  fOutputContainer->Add(fHistMCVertexZ);
  fHistPiPiPDGCode = new TH1F("fHistPiPiPDGCode","PDG code of K0s mothers",3500,0,3500);
  fOutputContainer->Add(fHistPiPiPDGCode);
  fHistPiPPDGCode = new TH1F("fHistPiPPDGCode","PDG code of #Lambda  mothers",3500,0,3500);
  fOutputContainer->Add(fHistPiPPDGCode);
  fHistPiAPPDGCode = new TH1F("fHistPiAPPDGCode","PDG code of #bar{#Lambda} mothers",3500,0,3500);
  fOutputContainer->Add(fHistPiAPPDGCode);
  fHistPiPCosPointAngXiVsPt= new TH2F("fHistPiPCosPointAngXiVsPt","pi-p cos of pointing angle vs pt from xi",200,0.0,20.0,250,0.99,1.00);
  fOutputContainer->Add(fHistPiPCosPointAngXiVsPt);
  fHistPiAPCosPointAngXiVsPt= new TH2F("fHistPiAPCosPointAngXiVsPt","pi+p- cos of pointing angle vs pt from xi",200,0.0,20.0,250,0.98,1.00);	
  fOutputContainer->Add(fHistPiAPCosPointAngXiVsPt);    



  // ------------------- add always ---------------------------//
  fHistV0RadiusZ[0]  = new TH2F("fHistV0RadiusZ","z of decay radius vs 2D radius",100,0.0,100.0,250,-125.0,125.0);
  fOutputContainer->Add(fHistV0RadiusZ[0]);

  fHistV0RadiusZVSPt[0]  = new TH2F("fHistV0RadiusZVSPt","z of decay radius vs pt radius",200,0.0,20.0,125,0.0,125.0);
  fOutputContainer->Add( fHistV0RadiusZVSPt[0]);
  fHistV0RadiusXY[0]  = new TH2F("fHistV0RadiusXY","y vs x decay radius",250,-125.0,125.0,250,-125.0,125.0);
  fOutputContainer->Add( fHistV0RadiusXY[0]);
  fHistV0RadiusXYVSY[0]  = new TH2F("fHistV0RadiusXYVSY","2D decay radius vs rap",100,-10,10,100,0.0,100.0);
  fOutputContainer->Add( fHistV0RadiusXYVSY[0] );
  fHistArmenteros[0] = new TH2F("fHistArmenteros"," pi+pi- armenteros",nbMass,-1.,1.,500,0.,0.5);
  fOutputContainer->Add(fHistArmenteros[0]);
  //***************************************** K0s **********************************//
  fHistPiPiMass[0] = new TH1F("fHistPiPiMass"," pi+pi- InvMass distribution",2*nbMass,0.,2.);
  fOutputContainer->Add(  fHistPiPiMass[0]);
  fHistPiPiMassVSPt[0] = new TH2F("fHistPiPiMassVSPt","pi+pi- InvMass distribution",nbMass,0.25,0.75,300,0.0,30.0);
  fOutputContainer->Add( fHistPiPiMassVSPt[0]);
  fHistPiPiMassVSPtMCTruth[0] = new TH2F("fHistPiPiMassVSPtMCTruth","pi+pi- InvMass distribution vs pt MCTruth",nbMass,0.25,0.75,300,0.0,30.0);
  fOutputContainer->Add( fHistPiPiMassVSPtMCTruth[0]);
  fHistPiPiPtVSY[0] = new TH2F("fHistPiPiPtVSY","p{t} vs y",100,-1,1,100,0.0,20);
  fOutputContainer->Add(  fHistPiPiPtVSY[0]);
  fHistPiPiDecayLengthVsPt[0] = new TH2F("fHistPiPiDecayLengthVsPt","K0 decay length vs pt",200,0.0,20.0,200,0.0,100.0);
  fOutputContainer->Add(fHistPiPiDecayLengthVsPt[0]); 
  fHistPiPiMonitorCuts[0] = new TH1F("fHistPiPiMonitorCuts","K0 cut monitor",25,0.5,25.5);
  fOutputContainer->Add(  fHistPiPiMonitorCuts[0]);
  //***************************************** Lambda **********************************//
  fHistPiPMass[0] = new TH1F("fHistPiPMass"," p+pi- InvMass distribution",2*nbMass,0.,2.);
  fOutputContainer->Add(  fHistPiPMass[0]); 
  fHistPiPMassVSPt[0] = new TH2F("fHistPiPMassVSPt","p+pi- InvMass distribution",nbMass,1.05,1.25,300,0.0,30.0);
  fOutputContainer->Add( fHistPiPMassVSPt[0]);
  fHistPiPMassVSPtMCTruth[0] = new TH2F("fHistPiPMassVSPtMCTruth","p+pi- InvMass distribution vs pt MCTruth",nbMass,1.05,1.25,300,0.0,30.0);
  fOutputContainer->Add( fHistPiPMassVSPtMCTruth[0]);
  fHistPiPPtVSY[0] = new TH2F("fHistPiPPtVSY","p{t} vs y",100,-1,1,100,0.0,20);
  fOutputContainer->Add( fHistPiPPtVSY[0]);
  fHistPiPDecayLengthVsPt[0] = new TH2F("fHistPiPDecayLengthVsPt","#Lambda decay length vs pt",200,0.0,20.0,200,0.0,100.0);
  fOutputContainer->Add( fHistPiPDecayLengthVsPt[0] );
  fHistPiPMonitorCuts[0] = new TH1F("fHistPiPMonitorCuts","#Lambda cut monitor",25,0.5,25.5);
  fOutputContainer->Add(  fHistPiPMonitorCuts[0]);
  //***************************************** AntiLambda **********************************//
  fHistPiAPMass[0] = new TH1F("fHistPiAPMass"," ap-pi+ InvMass distribution",2*nbMass,0.,2.);
  fOutputContainer->Add( fHistPiAPMass[0] );
  fHistPiAPMassVSPt[0] = new TH2F("fHistPiAPMassVSPt","p-pi+ InvMass distribution",nbMass,1.05,1.25,300,0.0,30.0);
  fOutputContainer->Add( fHistPiAPMassVSPt[0]);
  fHistPiAPMassVSPtMCTruth[0] = new TH2F("fHistPiAPMassVSPtMCTruth","p-pi+ InvMass distribution vs pt MCTruth",nbMass,1.05,1.25,300,0.0,30.0);
  fOutputContainer->Add(	 fHistPiAPMassVSPtMCTruth[0]);	       
  fHistPiAPPtVSY[0] = new TH2F("fHistPiAPPtVSY","p{t} vs y",100,-1,1,100,0.0,20);
  fOutputContainer->Add( fHistPiAPPtVSY[0]);
  fHistPiAPDecayLengthVsPt[0] = new TH2F("fHistPiAPDecayLengthVsPt","#bar{#Lambda} decay length vs pt",200,0.0,20.0,200,0.0,100.0);
  fOutputContainer->Add( fHistPiAPDecayLengthVsPt[0] );
  fHistPiAPMonitorCuts[0] = new TH1F("fHistPiAPMonitorCuts","#bar{#Lambda} cut monitor",25,0.5,25.5);   
  fOutputContainer->Add(  fHistPiAPMonitorCuts[0]);
  //  if(mchist==2){// for MC reco
  fHistV0RadiusZ[1]  = new TH2F("fHistV0RadiusZSec","z of decay radius vs 2D radius",100,0.0,100.0,250,-125.0,125.0);
  fOutputContainer->Add( fHistV0RadiusZ[1]); 
  fHistV0RadiusZVSPt[1]  = new TH2F("fHistV0RadiusZVSPtSec","z of decay radius vs pt radius",200,0.0,20.0,125,0.0,125.0);
  fOutputContainer->Add( fHistV0RadiusZVSPt[1]);
  fHistV0RadiusXY[1]  = new TH2F("fHistV0RadiusXYSec","y vs x decay radius",250,-125.0,125.0,250,-125.0,125.0);
  fOutputContainer->Add(  fHistV0RadiusXY[1] );
  fHistV0RadiusXYVSY[1]  = new TH2F("fHistV0RadiusXYVSYSec","2D decay radius vs rap",100,-10,10,100,0.0,100.0);
  fOutputContainer->Add( fHistV0RadiusXYVSY[1]); 
  fHistArmenteros[1] = new TH2F("fHistArmenterosSec"," pi+pi- armenteros",nbMass,-1.,1.,500,0.,0.5);
  fOutputContainer->Add( fHistArmenteros[1]);
  //***************************************** K0s **********************************//
  fHistPiPiMass[1] = new TH1F("fHistPiPiMassSec"," pi+pi- InvMass distribution",2*nbMass,0.,2.);
  fOutputContainer->Add( fHistPiPiMass[1] );
  fHistPiPiMassVSPt[1] = new TH2F("fHistPiPiMassVSPtSec","pi+pi- InvMass distribution",nbMass,0.25,0.75,300,0.0,30.0);
  fOutputContainer->Add(fHistPiPiMassVSPt[1]);
  fHistPiPiMassVSPtMCTruth[1] = new TH2F("fHistPiPiMassVSPtMCTruthSec","pi+pi- InvMass distribution vs pt MCTruth",nbMass,0.25,0.75,300,0.0,30.0);
  fOutputContainer->Add( fHistPiPiMassVSPtMCTruth[1] );
  fHistPiPiPtVSY[1] = new TH2F("fHistPiPiPtVSYSec","p{t} vs y",100,-1,1,100,0.0,20);
  fOutputContainer->Add(fHistPiPiPtVSY[1]);
  fHistPiPiDecayLengthVsPt[1] = new TH2F("fHistPiPiDecayLengthVsPtSec","K0 decay length vs pt",200,0.0,20.0,200,0.0,100.0);
  fOutputContainer->Add(fHistPiPiDecayLengthVsPt[1]);
  fHistPiPiMonitorCuts[1] = new TH1F("fHistPiPiMonitorCutsSec","K0 cut monitor",25,0.5,25.5);
  fOutputContainer->Add(fHistPiPiMonitorCuts[1]);
  //***************************************** Lambda **********************************//
  fHistPiPMass[1] = new TH1F("fHistPiPMassSec"," p+pi- InvMass distribution",2*nbMass,0.,2.);
  fOutputContainer->Add( fHistPiPMass[1]);
  fHistPiPMassVSPt[1] = new TH2F("fHistPiPMassVSPtSec","p+pi- InvMass distribution",nbMass,1.05,1.25,300,0.0,30.0);
  fOutputContainer->Add( fHistPiPMassVSPt[1]);
  fHistPiPMassVSPtMCTruth[1] = new TH2F("fHistPiPMassVSPtMCTruthSec","p+pi- InvMass distribution vs pt MCTruth",nbMass,1.05,1.25,300,0.0,30.0);
  fOutputContainer->Add(fHistPiPMassVSPtMCTruth[1]);
  fHistPiPPtVSY[1] = new TH2F("fHistPiPPtVSYSec","p{t} vs y",100,-1,1,100,0.0,20);
  fOutputContainer->Add( fHistPiPPtVSY[1]);
  fHistPiPDecayLengthVsPt[1] = new TH2F("fHistPiPDecayLengthVsPtSec","#Lambda decay length vs pt",200,0.0,20.0,200,0.0,100.0);
  fOutputContainer->Add(  fHistPiPDecayLengthVsPt[1]); 	
  fHistPiPMonitorCuts[1] = new TH1F("fHistPiPMonitorCutsSec","#Lambda cut monitor",25,0.5,25.5);
  fOutputContainer->Add( fHistPiPMonitorCuts[1]);
  //***************************************** AntiLambda **********************************//
  fHistPiAPMass[1] = new TH1F("fHistPiAPMassSec"," ap-pi+ InvMass distribution",2*nbMass,0.,2.);
  fOutputContainer->Add( fHistPiAPMass[1]);
  fHistPiAPMassVSPt[1] = new TH2F("fHistPiAPMassVSPtSec","p-pi+ InvMass distribution",nbMass,1.05,1.25,300,0.0,30.0);
  fOutputContainer->Add(    fHistPiAPMassVSPt[1]);
  fHistPiAPMassVSPtMCTruth[1] = new TH2F("fHistPiAPMassVSPtMCTruthSec","p-pi+ InvMass distribution vs pt MCTruth",nbMass,1.05,1.25,300,0.0,30.0);
  fOutputContainer->Add(fHistPiAPMassVSPtMCTruth[1]);
  fHistPiAPPtVSY[1] = new TH2F("fHistPiAPPtVSYSec","p{t} vs y",100,-1,1,100,0.0,20);
  fOutputContainer->Add( fHistPiAPPtVSY[1]);
  fHistPiAPDecayLengthVsPt[1] = new TH2F("fHistPiAPDecayLengthVsPtSec","#bar{#Lambda} decay length vs pt",200,0.0,20.0,200,0.0,100.0);
  fOutputContainer->Add(    fHistPiAPDecayLengthVsPt[1] );
  fHistPiAPMonitorCuts[1] = new TH1F("fHistPiAPMonitorCutsSec","#bar{#Lambda} cut monitor",25,0.5,25.5);
  fOutputContainer->Add( fHistPiAPMonitorCuts[1]);
      

  //----------------- for reco or data or mc data like MC reco only -----------------//
  
    
  fHistPiPiEtaDReco[0] = new TH2F("fHistPiPiEtaDRecoRaw","K0s daughters eta raw",300,-6,6,100,0,20);
  fOutputContainer->Add(fHistPiPiEtaDReco[0]);
  fHistPiPiEtaDReco[1] = new TH2F("fHistPiPiEtaDReco","K0s daughters eta after rap V0 cut pos",300,-3,3,300,-3.00,3.0);
  fOutputContainer->Add(fHistPiPiEtaDReco[1]);	 
  fHistPiPEtaDReco[0] = new TH2F("fHistPiPEtaDRecoRaw","#Lambda daughters eta raw",300,-6,6,100,0,20);
  fOutputContainer->Add(fHistPiPEtaDReco[0]);
  fHistPiPEtaDReco[1] = new TH2F("fHistPiPEtaDReco","#Lambda daughters eta after rap V0 cut neg",300,-3,3,300,-3.00,3.0);
  fOutputContainer->Add(fHistPiPEtaDReco[1]);
	
  //********************************************** K0 *************************************************//
  
  // fHistPiPiMassVSAlpha[0] = new TH2F("fHistPiPiMassVSAlpha"," alpha armenteros vs pi+pi- InvMass distribution",nbMass,0.25,0.75,500,-1.,1.);

  fHistPiPiDCADaughters[0] = new TH2F("fHistPiPiDCADaughters","dca of K0 daughters",nbMass,0.25,0.75,250,0.0,2);
  fOutputContainer->Add( fHistPiPiDCADaughters[0]); 
  fHistPiPiDCADaughterPosToPrimVtxVSMass[0] = new TH2F("fHistPiPiDCADaughterPosToPrimVtxVSMass","pi+ DCA daughter to prim vtx vsinvmass",nbMass,0.25,0.75,250,0.0,25.0);
  fOutputContainer->Add(  fHistPiPiDCADaughterPosToPrimVtxVSMass[0] );
  fHistPiPiDCAVSMass[0] = new TH2F("fHistPiPiDCAVSMass","pi+pi- dca  vs pt",nbMass,0.25,0.75,250,0.0,5.0);
  fOutputContainer->Add(  fHistPiPiDCAVSMass[0]);
  fHistPiPiCosPointAng[0]  = new TH2F("fHistPiPiCosPointAng","K0 cosine of pointing angle vs dca to prim vtx",200,0.0,10.0,250,0.99,1.00);
  fOutputContainer->Add( fHistPiPiCosPointAng[0] );
  fHistPiPiDecayLengthVsMass[0] = new TH2F("fHistPiPiDecayLengthVsMass","K0 decay length vs mass",nbMass,0.25,0.75,200,0.0,100.0);
  fOutputContainer->Add(fHistPiPiDecayLengthVsMass[0]);
  fHistPiPiRadiusXY[0] = new TH2F("fHistPiPiRadiusXY","pi+pi- opening angle vs mass",nbMass,0.25,0.75,200,0.0,4.0);
  fOutputContainer->Add( fHistPiPiRadiusXY[0]);
  // fHistPiPiPtDaughters[0] = new TH2F("fHistPiPiPtDaughters","p_{t} pos vs p_{t} neg of daughters",400,0.0,20.0,400,0,20.0);

	 
  //***************************************** Lambda *********************************************//
  fHistPiPDCADaughters[0] = new TH2F("fHistPiPDCADaughters","dca of #Lambda daughters",nbMass,1.05,1.25,250,0.0,2.0);
  fOutputContainer->Add(fHistPiPDCADaughters[0]);
  fHistPiPDCADaughterPosToPrimVtxVSMass[0] = new TH2F("fHistPiPDCADaughterPosToPrimVtxVSMass","p DCA daughter to prim vtx vs invmass",nbMass,1.05,1.25,250,0.0,25.0);
  fOutputContainer->Add( fHistPiPDCADaughterPosToPrimVtxVSMass[0]);
  fHistPiPDCAVSMass[0] = new TH2F("fHistPiPDCAVSMass","ppi- dca  vs pt",nbMass,1.05,1.25,250,0.0,5.0);
  fOutputContainer->Add(  fHistPiPDCAVSMass[0]);
  fHistPiPCosPointAng[0]  = new TH2F("fHistPiPCosPointAng","#Lambda cosine of pointing angle vs dca to prim vtx",200,0.0,10.0,250,0.99,1.00);
  fOutputContainer->Add( fHistPiPCosPointAng[0] );
  fHistPiPDecayLengthVsMass[0] = new TH2F("fHistPiPDecayLengthVsMass","#Lambda decay length vs mass",nbMass,1.05,1.25,200,0.0,100.0);
  fOutputContainer->Add( fHistPiPDecayLengthVsMass[0]);
  fHistPiPRadiusXY[0] = new TH2F("fHistPiPRadiusXY","pi-p+ opening angle vs mass",nbMass,1.05,1.25,200,0.0,4.0);
  fOutputContainer->Add(   fHistPiPRadiusXY[0]);
  // fHistPiPPtDaughters[0] = new TH2F("fHistPiPPtDaughters","p_{t} pos vs p_{t} neg of daughters",400,0.0,20.0,400,0,20.0);
     
  //************************************** Antilambda **********************************************//
  fHistPiAPDCADaughters[0] = new TH2F("fHistPiAPDCADaughters","dca of #bar{#Lambda} daughters",nbMass,1.05,1.25,250,0.0,2.0);
  fOutputContainer->Add( fHistPiAPDCADaughters[0]);
  fHistPiAPDCADaughterPosToPrimVtxVSMass[0] = new TH2F("fHistPiAPDCADaughterPosToPrimVtxVSMass","pi+ DCA daughter to prim vtx vs invmass",nbMass,1.05,1.25,250,0.0,25.0);
  fOutputContainer->Add(fHistPiAPDCADaughterPosToPrimVtxVSMass[0]);
  fHistPiAPDCAVSMass[0] = new TH2F("fHistPiAPDCAVSMass","pi+p- dca  vs pt",nbMass,1.05,1.25,250,0.0,5.0);
  fOutputContainer->Add( fHistPiAPDCAVSMass[0]);
  fHistPiAPCosPointAng [0] = new TH2F("fHistPiAPCosPointAng","#bar{#Lambda} cosine of pointing angle vs dcs to prim vtx",200,0.0,10.0,250,0.99,1.00);
  fOutputContainer->Add( fHistPiAPCosPointAng [0] );
  fHistPiAPDecayLengthVsMass[0] = new TH2F("fHistPiAPDecayLengthVsMass","#bar{#Lambda} decay length vs mass",nbMass,1.05,1.25,200,0.0,100.0);
  fOutputContainer->Add( fHistPiAPDecayLengthVsMass[0]);
  fHistPiAPRadiusXY[0] = new TH2F("fHistPiAPRadiusXY","pi+p- opening angle vs mass",nbMass,1.05,1.25,200,0.0,4.0);
  fOutputContainer->Add( fHistPiAPRadiusXY[0]);
  // fHistPiAPPtDaughters[0] = new TH2F("fHistPiAPPtDaughters","p_{t} pos vs p_{t} neg of daughters",400,0.0,20.0,400,0,20.0);
   
  //**********************************************TPC*****************************************************//

  fHistDedxSecProt[0] = new TH2F("fHistDedxSecProt","proton", nbPt, 0, 20, 100, 0, 400);
  fOutputContainer->Add( fHistDedxSecProt[0]);
  fHistDedxSecPiPlus[0] = new TH2F("fHistDedxSecPiPlus","pi plus", nbPt, 0, 20, 100, 0, 400);
  fOutputContainer->Add( fHistDedxSecPiPlus[0]);
  fHistDedxSecAProt[0] = new TH2F("fHistDedxSecAProt","antiproton", nbPt, 0, 20, 100, 0, 400);
  fOutputContainer->Add( fHistDedxSecAProt[0]);
  fHistDedxSecPiMinus[0] = new TH2F("fHistDedxSecPiMinus","pi minus", nbPt, 0, 20, 100, 0, 400);
  fOutputContainer->Add( fHistDedxSecPiMinus[0]);
  fHistNclsITSPosK0[0] = new TH1F("fHistNclsITSPosK0","fHistNclsITSPos K0",10,-0.5,9.5);
  fOutputContainer->Add(fHistNclsITSPosK0[0]);
  fHistNclsITSNegK0[0] = new TH1F("fHistNclsITSNegK0","fHistNclsITSNeg K0",10,-0.5,9.5);
  fOutputContainer->Add( fHistNclsITSNegK0[0]);
  fHistNclsTPCPosK0[0] = new TH1F("fHistNclsTPCPosK0","fHistNclsTPCPos K0",200,-0.5,199.5);
  fOutputContainer->Add( fHistNclsTPCPosK0[0]);
  fHistNclsTPCNegK0[0] = new TH1F("fHistNclsTPCNegK0","fHistNclsTPCNeg K0",200,-0.5,199.5);
  fOutputContainer->Add( fHistNclsTPCNegK0[0] );
  fHistChi2PerNclsITSPosK0[0] = new TH1F("fHistChi2PerNclsITSPosK0","chi2 per cluster ITS pi+pi-",250,0.0,50.0);
  fOutputContainer->Add(fHistChi2PerNclsITSPosK0[0]);
  fHistChi2PerNclsITSNegK0[0] = new TH1F("fHistChi2PerNclsITSNegK0","chi2 per cluster ITS pi+pi- neg",250,0.0,50.0);
  fOutputContainer->Add(fHistChi2PerNclsITSNegK0[0]);
  fHistNclsITSPosL[0] = new TH1F("fHistNclsITSPosL","fHistNclsITSPos #Lambda",10,-0.5,9.5);
  fOutputContainer->Add( fHistNclsITSPosL[0] );
  fHistNclsITSNegL[0] = new TH1F("fHistNclsITSNegL","fHistNclsITSNeg #Lambda",10,-0.5,9.5);
  fOutputContainer->Add(fHistNclsITSNegL[0]);
  fHistNclsTPCPosL[0] = new TH1F("fHistNclsTPCPosL","fHistNclsTPCPos #Lambda",200,-0.5,199.5);
  fOutputContainer->Add(fHistNclsTPCPosL[0]);
  fHistNclsTPCNegL[0] = new TH1F("fHistNclsTPCNegL","fHistNclsTPCNeg #Lambda",200,-0.5,199.5);
  fOutputContainer->Add( fHistNclsTPCNegL[0]);
  fHistChi2PerNclsITSPosL[0] = new TH1F("fHistChi2PerNclsITSPosL","chi2 per cluster ITS pi-p+ pos",250,0.0,50.0);
  fOutputContainer->Add( fHistChi2PerNclsITSPosL[0]);
  fHistChi2PerNclsITSNegL[0] = new TH1F("fHistChi2PerNclsITSNegL","chi2 per cluster ITS pi-p+ neg",250,0.0,50.0);
  fOutputContainer->Add( fHistChi2PerNclsITSNegL[0]);
  //for 2d pt dep studies
  fHistNclsITSPos[0] = new TH2F("fHistNclsITSPos","fHistNclsITSPos  vs pt pos",200,0.0,20.0,10,-0.5,9.5);
  fOutputContainer->Add( fHistNclsITSPos[0]);
  fHistNclsITSNeg[0] = new TH2F("fHistNclsITSNeg","fHistNclsITSNeg vs pt neg",200,0.0,20.0,10,-0.5,9.5);
  fOutputContainer->Add( fHistNclsITSNeg[0]);
  fHistNclsTPCPos[0] = new TH2F("fHistNclsTPCPos","ncls vs findable ncls pos",200,0.0,200.0,200,0.0,200.0);
  fOutputContainer->Add(fHistNclsTPCPos[0] );
  fHistNclsTPCNeg[0] = new TH2F("fHistNclsTPCNeg","ncls vs findable ncls neg",200,0.0,200.0,200,0.0,200.0);
  fOutputContainer->Add(   fHistNclsTPCNeg[0]);
  fHistChi2PerNclsITSPos[0] = new TH2F("fHistChi2PerNclsITSPos","chi2 per cluster ITS pi+p-vs pt pos",200,0.0,20.0,250,0.0,50.0);
  fOutputContainer->Add(  fHistChi2PerNclsITSPos[0] );
  fHistChi2PerNclsITSNeg[0] = new TH2F("fHistChi2PerNclsITSNeg","chi2 per cluster ITS pi+p- neg vs pt neg",200,0.0,20.0,250,0.0,50.0);
  fOutputContainer->Add( fHistChi2PerNclsITSNeg[0] );
  fHistNclsITS[0] = new TH2F("fHistNclsITS","fHistNclsITS pos vs neg",10,-0.5,9.5,10,-0.5,9.5);
  fOutputContainer->Add( fHistNclsITS[0]);
  fHistNclsTPC[0] = new TH2F("fHistNclsTPC","ncls TPC neg vs crossed rows neg",200,-0.5,199.5,200,-0.5,199.5);
  fOutputContainer->Add( fHistNclsTPC[0]);
  fHistNCRowsTPCPos[0] = new TH2F("fHistNCRowsTPCPos","n crossed rows vs pt pos",200,0.0,20.0,200,0.0,200.0);
  fOutputContainer->Add(  fHistNCRowsTPCPos[0] );
  fHistNCRowsTPCNeg[0] = new TH2F("fHistNCRowsTPCNeg","n crossed rows vs pt neg",200,0.0,20.0,200,0.0,200.0);
  fOutputContainer->Add( fHistNCRowsTPCNeg[0] );

   
  // if(mchist==2){// for MC reco

  //********************************************** K0 *************************************************//
  // fHistPiPiMassVSAlpha[0] = new TH2F("fHistPiPiMassVSAlpha"," alpha armenteros vs pi+pi- InvMass distribution",nbMass,0.25,0.75,500,-1.,1.);
  fHistPiPiDCADaughters[1] = new TH2F("fHistPiPiDCADaughtersSec","dca of K0 daughters",nbMass,0.25,0.75,250,0.0,2);
  fOutputContainer->Add(   fHistPiPiDCADaughters[1] );
  fHistPiPiDCADaughterPosToPrimVtxVSMass[1] = new TH2F("fHistPiPiDCADaughterPosToPrimVtxVSMassSec","pi+ DCA daughter to prim vtx vsinvmass",nbMass,0.25,0.75,250,0.0,25.0);
  fOutputContainer->Add(  fHistPiPiDCADaughterPosToPrimVtxVSMass[1]);
  fHistPiPiDCAVSMass[1] = new TH2F("fHistPiPiDCAVSMassSec","pi+pi- dca  vs pt",nbMass,0.25,0.75,250,0.0,5.0);
  fOutputContainer->Add(fHistPiPiDCAVSMass[1]);
  fHistPiPiCosPointAng[1]  = new TH2F("fHistPiPiCosPointAngSec","K0 cosine of pointing angle vs dca to prim vtx",200,0.0,10.0,250,0.99,1.00);
  fOutputContainer->Add( fHistPiPiCosPointAng[1]);
  fHistPiPiDecayLengthVsMass[1] = new TH2F("fHistPiPiDecayLengthVsMassSec","K0 decay length vs mass",nbMass,0.25,0.75,200,0.0,100.0);
  fOutputContainer->Add( fHistPiPiDecayLengthVsMass[1]);
  fHistPiPiRadiusXY[1] = new TH2F("fHistPiPiRadiusXYSec","pi+pi- opening angle vs mass",nbMass,0.25,0.75,200,0.0,4.0);
  fOutputContainer->Add(fHistPiPiRadiusXY[1]);
  // fHistPiPiPtDaughters[0] = new TH2F("fHistPiPiPtDaughters","p_{t} pos vs p_{t} neg of daughters",400,0.0,20.0,400,0,20.0);

	 
  //********************************************** Lambda ***********************************************//
  fHistPiPDCADaughters[1] = new TH2F("fHistPiPDCADaughtersSec","dca of #Lambda daughters",nbMass,1.05,1.25,250,0.0,2.0);
  fOutputContainer->Add( fHistPiPDCADaughters[1]);
  fHistPiPDCADaughterPosToPrimVtxVSMass[1] = new TH2F("fHistPiPDCADaughterPosToPrimVtxVSMassSec","p DCA daughter to prim vtx vs invmass",nbMass,1.05,1.25,250,0.0,25.0);
  fOutputContainer->Add( fHistPiPDCADaughterPosToPrimVtxVSMass[1]);
  fHistPiPDCAVSMass[1] = new TH2F("fHistPiPDCAVSMassSec","ppi- dca  vs pt",nbMass,1.05,1.25,250,0.0,5.0);
  fOutputContainer->Add( fHistPiPDCAVSMass[1] );
  fHistPiPCosPointAng[1]  = new TH2F("fHistPiPCosPointAngSec","#Lambda cosine of pointing angle vs dca to prim vtx",200,0.0,10.0,250,0.99,1.00);
  fOutputContainer->Add( fHistPiPCosPointAng[1] );
  fHistPiPDecayLengthVsMass[1] = new TH2F("fHistPiPDecayLengthVsMassSec","#Lambda decay length vs mass",nbMass,1.05,1.25,200,0.0,100.0);
  fOutputContainer->Add( fHistPiPDecayLengthVsMass[1]);
  fHistPiPRadiusXY[1] = new TH2F("fHistPiPRadiusXYSec","pi-p+ opening angle vs mass",nbMass,1.05,1.25,200,0.0,4.0);
  // fHistPiPPtDaughters[0] = new TH2F("fHistPiPPtDaughters","p_{t} pos vs p_{t} neg of daughters",400,0.0,20.0,400,0,20.0);
  fOutputContainer->Add(  fHistPiPRadiusXY[1]);
  //******************************************* Antilambda **********************************************//
  fHistPiAPDCADaughters[1] = new TH2F("fHistPiAPDCADaughtersSec","dca of #bar{#Lambda} daughters",nbMass,1.05,1.25,250,0.0,2.0);
  fOutputContainer->Add(  fHistPiAPDCADaughters[1]);
  fHistPiAPDCADaughterPosToPrimVtxVSMass[1] = new TH2F("fHistPiAPDCADaughterPosToPrimVtxVSMassSec","pi+ DCA daughter to prim vtx vs invmass",nbMass,1.05,1.25,250,0.0,25.0);
  fOutputContainer->Add( fHistPiAPDCADaughterPosToPrimVtxVSMass[1]);
  fHistPiAPDCAVSMass[1] = new TH2F("fHistPiAPDCAVSMassSec","pi+p- dca  vs pt",nbMass,1.05,1.25,250,0.0,5.0);
  fOutputContainer->Add(    fHistPiAPDCAVSMass[1]);
  fHistPiAPCosPointAng[1] = new TH2F("fHistPiAPCosPointAngSec","#bar{#Lambda} cosine of pointing angle vs dcs to prim vtx",200,0.0,10.0,250,0.99,1.00);
  fOutputContainer->Add(  fHistPiAPCosPointAng[1]);
  fHistPiAPDecayLengthVsMass[1] = new TH2F("fHistPiAPDecayLengthVsMassSec","#bar{#Lambda} decay length vs mass",nbMass,1.05,1.25,200,0.0,100.0);
  fOutputContainer->Add( fHistPiAPDecayLengthVsMass[1]);
  fHistPiAPRadiusXY[1] = new TH2F("fHistPiAPRadiusXYSec","pi+p- opening angle vs mass",nbMass,1.05,1.25,200,0.0,4.0);
  fOutputContainer->Add( fHistPiAPRadiusXY[1] );
  // fHistPiAPPtDaughters[0] = new TH2F("fHistPiAPPtDaughters","p_{t} pos vs p_{t} neg of daughters",400,0.0,20.0,400,0,20.0);

  //******************************************* TPC ****************************************************//
	      
  fHistDedxSecProt[1] = new TH2F("fHistDedxSecProtSec","proton", nbPt, 0, 20, 100, 0, 400);
  fOutputContainer->Add( fHistDedxSecProt[1]); 
  fHistDedxSecPiPlus[1] = new TH2F("fHistDedxSecPiPlusSec","pi plus", nbPt, 0, 20, 100, 0, 400);
  fOutputContainer->Add(  fHistDedxSecPiPlus[1]);
  fHistDedxSecAProt[1] = new TH2F("fHistDedxSecAProtSec","antiproton", nbPt, 0, 20, 100, 0, 400);
  fOutputContainer->Add(fHistDedxSecAProt[1]);
  fHistDedxSecPiMinus[1] = new TH2F("fHistDedxSecPiMinusSec","pi minus", nbPt, 0, 20, 100, 0, 400);
  fOutputContainer->Add(   fHistDedxSecPiMinus[1]);
  fHistNclsITSPosK0[1] = new TH1F("fHistNclsITSPosK0Sec","fHistNclsITSPos K0 ",10,-0.5,9.5);
  fOutputContainer->Add(  fHistNclsITSPosK0[1]);
  fHistNclsITSNegK0[1] = new TH1F("fHistNclsITSNegK0Sec","fHistNclsITSNeg K0",10,-0.5,9.5);
  fOutputContainer->Add( fHistNclsITSNegK0[1]);
  fHistNclsTPCPosK0[1] = new TH1F("fHistNclsTPCPosK0Sec","fHistNclsTPCPos K0",200,-0.5,199.5);
  fOutputContainer->Add(  fHistNclsTPCPosK0[1]);
  fHistNclsTPCNegK0[1] = new TH1F("fHistNclsTPCNegK0Sec","fHistNclsTPCNeg K0",200,-0.5,199.5);
  fOutputContainer->Add(fHistNclsTPCNegK0[1] );
  fHistChi2PerNclsITSPosK0[1] = new TH1F("fHistChi2PerNclsITSPosK0Sec","chi2 per cluster ITS pi+pi- pos",250,0.0,50.0);
  fOutputContainer->Add(   fHistChi2PerNclsITSPosK0[1] );
  fHistChi2PerNclsITSNegK0[1] = new TH1F("fHistChi2PerNclsITSNegK0Sec","chi2 per cluster ITS pi+pi- neg",250,0.0,50.0);
  fOutputContainer->Add(     fHistChi2PerNclsITSNegK0[1]);
  fHistNclsITSPosL[1] = new TH1F("fHistNclsITSPosLSec","fHistNclsITSPos #Lambda",10,-0.5,9.5);
  fOutputContainer->Add( fHistNclsITSPosL[1]);
  fHistNclsITSNegL[1] = new TH1F("fHistNclsITSNegLSec","fHistNclsITSNeg #Lambda",10,-0.5,9.5);
  fOutputContainer->Add(   fHistNclsITSNegL[1]);
  fHistNclsTPCPosL[1] = new TH1F("fHistNclsTPCPosLSec","fHistNclsTPCPos #Lambda",200,-0.5,199.5);
  fOutputContainer->Add( fHistNclsTPCPosL[1]);
  fHistNclsTPCNegL[1] = new TH1F("fHistNclsTPCNegLSec","fHistNclsTPCNeg #Lambda",200,-0.5,199.5);
  fOutputContainer->Add( fHistNclsTPCNegL[1] );
  fHistChi2PerNclsITSPosL[1] = new TH1F("fHistChi2PerNclsITSPosLSec","chi2 per cluster ITS pi-p+ pos",250,0.0,50.0);
  fOutputContainer->Add( fHistChi2PerNclsITSPosL[1]);
  fHistChi2PerNclsITSNegL[1] = new TH1F("fHistChi2PerNclsITSNegLSec","chi2 per cluster ITS pi-p+ neg",250,0.0,50.0);
  fOutputContainer->Add( fHistChi2PerNclsITSNegL[1]);
  //for 2d pt dep studies
  fHistNclsITSPos[1] = new TH2F("fHistNclsITSPosSec","fHistNclsITSPos  vs pt pos",200,0.0,20.0,10,-0.5,9.5);
  fOutputContainer->Add( fHistNclsITSPos[1]);
  fHistNclsITSNeg[1] = new TH2F("fHistNclsITSNegSec","fHistNclsITSNeg vs pt neg",200,0.0,20.0,10,-0.5,9.5);
  fOutputContainer->Add(fHistNclsITSNeg[1]);
  fHistNclsTPCPos[1] = new TH2F("fHistNclsTPCPosSec","ncls vs findable ncls pos",200,0.0,200.0,200,0.0,200.0);
  fOutputContainer->Add( fHistNclsTPCPos[1] );
  fHistNclsTPCNeg[1] = new TH2F("fHistNclsTPCNegSec","ncls vs findable ncls neg",200,0.0,200.0,200,0.0,200.0);
  fOutputContainer->Add( fHistNclsTPCNeg[1]);
  fHistChi2PerNclsITSPos[1] = new TH2F("fHistChi2PerNclsITSPosSec","chi2 per cluster ITS pi+p-vs pt pos",200,0.0,20.0,250,0.0,50.0);
  fOutputContainer->Add(  fHistChi2PerNclsITSPos[1]);
  fHistChi2PerNclsITSNeg[1] = new TH2F("fHistChi2PerNclsITSNegSec","chi2 per cluster ITS pi+p- neg vs pt neg",200,0.0,20.0,250,0.0,50.0);
  fOutputContainer->Add(  fHistChi2PerNclsITSNeg[1]);
  fHistNclsITS[1] = new TH2F("fHistNclsITSSec","fHistNclsITS pos vs neg",10,-0.5,9.5,10,-0.5,9.5);
  fOutputContainer->Add(  fHistNclsITS[1]);
  fHistNclsTPC[1] = new TH2F("fHistNclsTPCSec","ncls TPC neg vs crossed rows neg",200,-0.5,199.5,200,-0.5,199.5);
  fOutputContainer->Add(fHistNclsTPC[1] );
  fHistNCRowsTPCPos[1] = new TH2F("fHistNCRowsTPCPosSec","n crossed rows vs pt pos",200,0.0,20.0,200,0.0,200.0);
  fOutputContainer->Add( fHistNCRowsTPCPos[1]);
  fHistNCRowsTPCNeg[1] = new TH2F("fHistNCRowsTPCNegSec","n crossed rows vs pt neg",200,0.0,20.0,200,0.0,200.0);
  fOutputContainer->Add( fHistNCRowsTPCNeg[1]);

  fHistPiPiEtaDMC[0] = new TH2F("fHistPiPiEtaDMCRaw","K0s daughters etaMC raw",300,-6,6,100,0,20);//
  fOutputContainer->Add(fHistPiPiEtaDMC[0]);
  fHistPiPiEtaDMC[1] = new TH2F("fHistPiPiEtaDMC","K0s daughters etaMC after rap V0 cut",300,-6,6,100,0,20);
  fOutputContainer->Add(fHistPiPiEtaDMC[1]); 
  fHistPiPEtaDMC[0] = new TH2F("fHistPiPEtaDMCRaw","#Lambda daughters etaMC raw",300,-6,6,100,0,20);
  fOutputContainer->Add(fHistPiPEtaDMC[0]); 
  fHistPiPEtaDMC[1] = new TH2F("fHistPiPEtaDMC","#Lambda daughters etaMC after rap V0 cut",300,-6,6,100,0,20);
  fOutputContainer->Add(fHistPiPEtaDMC[1]);

  //********************************************** K0 *************************************************//
   
	 
	 
  //********************************************** Lambda **********************************************//
 
  fHistPiPMassVSPtSecSigma[0] = new TH2F("fHistPiPMassVSPtSecSigmaMC"," pi-p+ InvMass distribution secondaries from sigma MC",nbMass,1.05,1.25,200,0.,20);
  fOutputContainer->Add( fHistPiPMassVSPtSecSigma[0]);
  fHistPiPMassVSPtSecSigma[1] = new TH2F("fHistPiPMassVSPtSecSigma"," pi-p+ InvMass distribution secondaries from Sigma reco",nbMass,1.05,1.25,200,0.,20);
  fOutputContainer->Add( fHistPiPMassVSPtSecSigma[1]);
  fHistPiPMassVSPtSecXi[0] = new TH2F("fHistPiPMassVSPtSecXiMC"," pi-p+ InvMass distribution secondaries from  xi MC",nbMass,1.05,1.25,200,0.,20);
  fOutputContainer->Add(  fHistPiPMassVSPtSecXi[0]);
  fHistPiPMassVSPtSecXi[1] = new TH2F("fHistPiPMassVSPtSecXi"," pi-p+ InvMass distribution secondaries from  xi  reco",nbMass,1.05,1.25,200,0.,20);
  fOutputContainer->Add( fHistPiPMassVSPtSecXi[1] );
  fHistPiPMassVSYSecXi[0] = new TH2F("fHistPiPMassVSYSecXiMC"," pi-p+ InvMass distribution secondaries from xi MC",nbMass,1.05,1.25,100,-2.,2);
  fOutputContainer->Add( fHistPiPMassVSYSecXi[0]);
  fHistPiPMassVSYSecXi[1] = new TH2F("fHistPiPMassVSYSecXi"," pi-p+ InvMass distribution secondaries from xi reco",nbMass,1.05,1.25,100,-2.,2);
  fOutputContainer->Add(fHistPiPMassVSYSecXi[1] );
  fHistPiPXi0PtVSLambdaPt[0]= new TH2F("fHistPiPXi0PtVSLambdaPtMC"," pt xi 0 vs pt lambda MC truth",200,0.0,20.0,200,0.0,20.0);
  fOutputContainer->Add(fHistPiPXi0PtVSLambdaPt[0]);
  fHistPiPXi0PtVSLambdaPt[1]= new TH2F("fHistPiPXi0PtVSLambdaPt"," pt xi 0 truth vs pt lambda reco",200,0.0,20.0,200,0.0,20.0);
  fOutputContainer->Add( fHistPiPXi0PtVSLambdaPt[1]);
  fHistPiPXiMinusPtVSLambdaPt[0]= new TH2F("fHistPiPXiMinusPtVSLambdaPtMC","pt xi- vs pt lambda MC truth",200,0.0,20.0,200,0.0,20.0);
  fOutputContainer->Add( fHistPiPXiMinusPtVSLambdaPt[0]);
  fHistPiPXiMinusPtVSLambdaPt[1]= new TH2F("fHistPiPXiMinusPtVSLambdaPt","pt xi- truth vs pt lambda reco",200,0.0,20.0,200,0.0,20.0);
  fOutputContainer->Add( fHistPiPXiMinusPtVSLambdaPt[1]);
  //******************************************* Antilambda **********************************************//

  fHistPiAPMassVSPtSecSigma[0] = new TH2F("fHistPiAPMassVSPtSecSigmaMC"," pi+p- InvMass distribution secondaries from Sigma MC",nbMass,1.05,1.25,200,0.,20);
  fOutputContainer->Add(fHistPiAPMassVSPtSecSigma[0]);
  fHistPiAPMassVSPtSecSigma[1] = new TH2F("fHistPiAPMassVSPtSecSigma"," pi+p- InvMass distribution secondaries from  Sigma  reco",nbMass,1.05,1.25,200,0.,20);
  fOutputContainer->Add( fHistPiAPMassVSPtSecSigma[1]);
  fHistPiAPMassVSPtSecXi[0] = new TH2F("fHistPiAPMassVSPtSecXiMC"," pi+p- InvMass distribution secondaries from xi MC",nbMass,1.05,1.25,200,0.,20);
  fOutputContainer->Add( fHistPiAPMassVSPtSecXi[0]);
  fHistPiAPMassVSPtSecXi[1] = new TH2F("fHistPiAPMassVSPtSecXi"," pi+p- InvMass distribution secondaries from  Xi reco",nbMass,1.05,1.25,200,0.,20);
  fOutputContainer->Add(fHistPiAPMassVSPtSecXi[1]);

  fHistPiAPMassVSYSecXi[0] = new TH2F("fHistPiAPMassVSYSecXiMC"," pi+p- InvMass distribution secondaries from  xi MC",nbMass,1.05,1.25,100,-2,2);
  fOutputContainer->Add(fHistPiAPMassVSYSecXi[0]);
  fHistPiAPMassVSYSecXi[1] = new TH2F("fHistPiAPMassVSYSecXi"," pi+p- InvMass distribution secondaries from xi reco",nbMass,1.05,1.25,100,-2.,2);
  fOutputContainer->Add(    fHistPiAPMassVSYSecXi[1]);
  fHistPiAPXi0PtVSLambdaPt[0]= new TH2F("fHistPiAPXi0PtVSLambdaPtMC"," pt xi 0 vs pt Alambda MC truth",200,0.0,20.0,200,0.0,20.0);
  fOutputContainer->Add( fHistPiAPXi0PtVSLambdaPt[0]);
  fHistPiAPXi0PtVSLambdaPt[1]= new TH2F("fHistPiAPXi0PtVSLambdaPt"," pt xi 0 truth vs pt Alambda reco",200,0.0,20.0,200,0.0,20.0);
  fOutputContainer->Add(fHistPiAPXi0PtVSLambdaPt[1]);
  fHistPiAPXiMinusPtVSLambdaPt[0]= new TH2F("fHistPiAPXiMinusPtVSLambdaPtMC","pt xi- vs pt Alambda MC truth",200,0.0,20.0,200,0.0,20.0);
  fOutputContainer->Add(fHistPiAPXiMinusPtVSLambdaPt[0]);
  fHistPiAPXiMinusPtVSLambdaPt[1]= new TH2F("fHistPiAPXiMinusPtVSLambdaPt","pt xi- truth vs pt Alambda reco",200,0.0,20.0,200,0.0,20.0);
  fOutputContainer->Add( fHistPiAPXiMinusPtVSLambdaPt[1]);
   
  PostData(1,fOutputContainer);	       
 
}

//________________________________________________________________________
void AliAnalysisTaskV0ForRAA::UserExec(Option_t *) {
  //user exec

  //-- esd handler --//
  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> 
    (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!esdH) {
    printf("ERROR: Could not get ESDInputHandler");
    return;
  } 
  fESD = esdH->GetEvent();
  if(!fESD) {
    printf("ERROR: fESD not available \n");
    return ;
  }

  //-- mc handler --//
  if(fMCMode || fMCTruthMode){
    AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*> 
      (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if(!mcH) {
      printf("ERROR: Could not get MCInputHandler");
      return;
    }
    fMCev = mcH->MCEvent();
    if (!fMCev) {
      printf("ERROR: fMCev not available \n");
      return ;
    }
  }
    

  //--  AliPIDResponse --//
  fESDpid = esdH->GetPIDResponse();
   
  //-- Count events before cuts --//
  fHistNEvents->Fill(0);

  //-- Check object existence --//
  const AliESDVertex *    vtxESD    = fESD->GetPrimaryVertexTracks();
  const AliESDVertex *    vtxESDTPC = fESD->GetPrimaryVertexTPC();  
  const AliESDVertex *    vtxESDSPD = fESD->GetPrimaryVertexSPD();  
  const AliMultiplicity * multESD   = fESD->GetMultiplicity();  

  if ( !vtxESD ){
    AliError("No Tracks Vertex");
    return;
  }

  if ( !vtxESDTPC ){
    AliError("No TPC Vertex");
    return ;
  }

  if ( !vtxESDSPD ){
    AliError("No SPD Vertex");
    return ;
  }

  if ( !multESD ){
    AliError("No Multiplicity");
    return ;
  }
   

  // ----------- MC vertex -----------------------------------//
 
  Int_t nContr =0;
  if(fMCMode && fMCTruthMode){
    Double_t vVertexPrim[3];
    fMCev->GetPrimaryVertex()->GetXYZ(vVertexPrim);
    fHistMCVertexZ->Fill(vVertexPrim[2]);
      
    if (vtxESD->GetStatus()){
      nContr=vtxESD->GetNContributors();
      fHistPrimVtxZESDVSNContributorsMC->Fill(vVertexPrim[2],nContr);
      fHistPrimVtxZESDTPCVSNContributorsMC->Fill(vVertexPrim[2],nContr);
    }
    else {
      if(vtxESDSPD->GetStatus()){
	nContr=vtxESDSPD->GetNContributors();
	fHistPrimVtxZESDTPCVSNContributorsMC->Fill(vVertexPrim[2],nContr);
	fHistPrimVtxZESDSPDVSNContributorsMC->Fill(vVertexPrim[2],nContr);
      }
      else{
	fHistPrimVtxZESDVSNContributorsMC->Fill(vVertexPrim[2],nContr);//add for correction ESD and ESDPSD!!!!
	fHistPrimVtxZESDTPCVSNContributorsMC->Fill(vVertexPrim[2],nContr);
      }
    }
  }
  
  //-- Monitor event cuts --//
  fHistNEvents->Fill(1);
  
  // -- Check fo centrality
  Bool_t process = kTRUE;
  Int_t centBin = -1;
  if(fUseCentrality) {
    centBin = CalculateCentralityBin();
    if(!fUseCentralityRange){
      if(centBin!= fUseCentralityBin) process=kFALSE;
    }
    else if(centBin < fUseCentralityBin || centBin > fUseCentralityBin+fUseCentralityRange)
      process = kFALSE;
  }

  AliESDVZERO* esdV0 = fESD->GetVZEROData();
  Float_t multV0 = esdV0->GetMTotV0A() + esdV0->GetMTotV0C();
  
  if(fAnapp){// pp Analysis
    // SDD test for 2.76TeV pp
    // select events with SDD
    //    TString trCl = fESD->GetFiredTriggerClasses();
    //if(!(trCl.Contains("ALLNOTRD")) && fSelSDD) return;
    UInt_t maskSel = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    if(maskSel& AliVEvent::kFastOnly && fSelSDD) return;  

    Int_t ntracks = fESD->GetNumberOfTracks();
    for(Int_t i=0;i<ntracks;i++){//check sdd event selection
      AliESDtrack *tr=   fESD->GetTrack(i);
      
      Bool_t sdd0 = tr->HasPointOnITSLayer(0);
      Bool_t sdd1 = tr->HasPointOnITSLayer(1);
      Bool_t sdd2 = tr->HasPointOnITSLayer(2);
      Bool_t sdd3 = tr->HasPointOnITSLayer(3);
      Bool_t sdd4 = tr->HasPointOnITSLayer(4);
      Bool_t sdd5 = tr->HasPointOnITSLayer(5);
       
      fHistITSLayerHits->Fill(Int_t(sdd0)*(-1),ntracks);
      fHistITSLayerHits->Fill(Int_t(sdd1)*1,ntracks);
      fHistITSLayerHits->Fill(Int_t(sdd2)*2,ntracks);
      fHistITSLayerHits->Fill(Int_t(sdd3)*3,ntracks);
      fHistITSLayerHits->Fill(Int_t(sdd4)*4,ntracks);
      fHistITSLayerHits->Fill(Int_t(sdd5)*5,ntracks);
    }
      
    //vertex selection
    if (vtxESD->GetStatus()){
      fHistNEvents->Fill(2);
      fHistESDVertexZ->Fill(vtxESD->GetZv());
      if(fabs(vtxESD->GetZv()) < fVertexZCut){
	fHistMuliplicityRaw->Fill(multV0);
	fHistNEvents->Fill(3);
	fHistNPrim->Fill(nContr);
	
	Process();
	
	fHistMuliplicity->Fill(multV0);
	
	nContr = vtxESD->GetNContributors();
	//  if(nContr<501){
	fHistPrimVtxZESDVSNContributors->Fill(vtxESD->GetZv(),nContr);
	fHistPrimVtxZESDTPCVSNContributors->Fill(vtxESDTPC->GetZv(),nContr);
	//fHistPrimVtxZESDSPDVSNContributorsTPC->Fill(vtxESDSPD->GetZv(),nContr);
	//   }
	fHistPrimVtxZESD->Fill(vtxESD->GetZv());
	fHistPrimVtxZESDTPC->Fill(vtxESDTPC->GetZv());
	// fHistPrimVtxZESDSPD->Fill(vtxESDSPD->GetZv());
	// -- count events after processing
	fHistNEvents->Fill(4);
      }
    }
    else{
      if(vtxESDSPD->GetStatus()){
	fHistNEvents->Fill(2);
	
	fHistESDVertexZ->Fill(vtxESDSPD->GetZv());
	if(fabs(vtxESDSPD->GetZv()) < fVertexZCut){
	  
	  fHistMuliplicityRaw->Fill(multV0);
	  fHistNEvents->Fill(3);
	  fHistNPrim->Fill(nContr);
	  
	  Process();
	  
	  fHistMuliplicity->Fill(multV0);
	  
	  nContr = vtxESDSPD->GetNContributors();
	  //  if(nContr<501){
	  //fHistPrimVtxZESDVSNContributors->Fill(vtxESD->GetZv(),nContr);
	  fHistPrimVtxZESDTPCVSNContributors->Fill(vtxESDTPC->GetZv(),nContr);
	  fHistPrimVtxZESDSPDVSNContributors->Fill(vtxESDSPD->GetZv(),nContr);
	  // }
	  // fHistPrimVtxZESD->Fill(vtxESD->GetZv());
	  fHistPrimVtxZESDTPC->Fill(vtxESDTPC->GetZv());
	  fHistPrimVtxZESDSPD->Fill(vtxESDSPD->GetZv());
	  // -- count events after processing
	  fHistNEvents->Fill(4);
	}
      }
      //else return;
    }
  }
  else{// PbPb analysis
    if(vtxESD->GetStatus()){
      Double_t vtxZ = vtxESD->GetZv();
      fHistESDVertexZ->Fill(vtxZ);
      if(process){
	fHistNEvents->Fill(2);
	if(fabs(vtxZ) < fVertexZCut){
	  nContr = vtxESD->GetNContributors();
	  fHistMuliplicityRaw->Fill(multV0);
	  fHistNEvents->Fill(3);
	  fHistNPrim->Fill(nContr);
	  Process();
	  fHistMuliplicity->Fill(multV0);
	  fHistPrimVtxZESD->Fill(vtxZ);
	  fHistPrimVtxZESDVSNContributors->Fill(vtxZ,nContr);
	  // -- count events after processing --//
	  fHistCentBin->Fill(centBin);
	  fHistNEvents->Fill(4);
	}
      }
      if(fabs(vtxZ) < fVertexZCut) fHistCentBinRaw->Fill(centBin);
    }
  }
 
    
}

//________________________________________________________________________
void AliAnalysisTaskV0ForRAA::Terminate(Option_t *) {
  //terminate
}
//________________________________________________________________________
void AliAnalysisTaskV0ForRAA::Process(){
  //run the analysis

  Int_t ntracks = fESD->GetNumberOfTracks();
  Int_t count = 0;

  //-- count number of tracks --//
   
  if(!(!fMCMode && fMCTruthMode)){
    for(Int_t i=0;i<ntracks;i++){
      AliESDtrack *track = (AliESDtrack*)fESD->GetTrack(i);
      if(!fESDTrackCuts->AcceptTrack(track)) continue;
      if( track->Eta() > fEtaCutMCDaughtersVal) continue;
      count++;
    }
    fHistMultiplicityPrimary->Fill(count);
  }
   
  //-- check number of V0s in case of data or mc data like analysis--//
  Int_t nV0 = fESD->GetNumberOfV0s();
  if(!fMCTruthMode) if(nV0 < 1) return;
   
  //-- run analysis --//
  if(fMCTruthMode)  V0MCTruthLoop();
  else  V0RecoLoop(0,0,0,0,0.0,0,0.0);

}
//________________________________________________________________________
void AliAnalysisTaskV0ForRAA::V0MCTruthLoop(){
  //loop over MC truth particles

  //-- get MC stack --//
  AliStack *stack = fMCev->Stack();

  /*
  //histo fo user defined shift in charge/pt 
  if(fShift){
  fHistUserPtShift->Fill(fDeltaInvP);
  }
  */
  /*
    AliKFVertex primVtxStart(*(fESD->GetPrimaryVertex()));
    Int_t nTracksPrim=primVtxStart.GetNContributors();
    fHistNPrim->Fill(nTracksPrim);
  */
  /*
  // MC
    
  Int_t mcPrimaries = stack->GetNprimary();
  Int_t mcParticles    = stack->GetNtrack();
    
  fHistMultiplicityPrimary->Fill(mcPrimaries);
  fHistMCMultiplicityTracks->Fill(mcParticles);
    
  // number of V0
  fHistNV0->Fill(nV0);
  if(nTracksPrim>0) {
  fHistNV0WithVertex->Fill(nV0);
  }
  */


  //-- MC truht loop for V0s --//
  for (Int_t iMc = 0; iMc < (stack->GetNtrack()); iMc++){//MC truth loop
    Int_t fillMCtruth= int(fMCTruthMode);
    if(fMCTruthMode){
      fHistPiPiMonitorCuts[0]->Fill(1*fillMCtruth);
      fHistPiPMonitorCuts[0]->Fill(1*fillMCtruth);
      fHistPiAPMonitorCuts[0]->Fill(1*fillMCtruth);
    }
    TParticle *p0 = stack->Particle(iMc);
    if(!p0) continue;

    if(fMCTruthMode){
      fHistPiPiMonitorCuts[0]->Fill(2*fillMCtruth);
      fHistPiPMonitorCuts[0]->Fill(2*fillMCtruth);
      fHistPiAPMonitorCuts[0]->Fill(2*fillMCtruth);
    }

    Int_t pdgCode = p0->GetPdgCode();

    //-------------- only K0s and Lambda ----------//
    if( (pdgCode != 310 ) && ( fabs(pdgCode) != 3122 ) ) continue;
    Int_t fillFlagK0 = (3122- fabs(pdgCode))/(3122-310)*fillMCtruth;
    Int_t fillFlagL = (fabs(pdgCode) - 310)/(3122-310)*(pdgCode+3122)/(2*3122)*fillMCtruth;
    Int_t fillFlagAL = (fabs(pdgCode) - 310)/(3122-310)*(pdgCode-3122)/(-2*3122)*fillMCtruth;
    
    fHistPiPiMonitorCuts[0]->Fill(3*fillFlagK0);
    fHistPiPMonitorCuts[0]->Fill(3*fillFlagL);
    fHistPiAPMonitorCuts[0]->Fill(3*fillFlagAL);
      
    if(p0->GetNDaughters() !=2) continue;
    fHistPiPiMonitorCuts[0]->Fill(4*fillFlagK0);
    fHistPiPMonitorCuts[0]->Fill(4*fillFlagL);
    fHistPiAPMonitorCuts[0]->Fill(4*fillFlagAL);
      
    //-------------- unique ID check-------------- //
    Int_t uniqueID =  p0->GetUniqueID();
    if(uniqueID==13) continue;
      
    fHistPiPiMonitorCuts[0]->Fill(5*fillFlagK0);
    fHistPiPMonitorCuts[0]->Fill(5*fillFlagL);
    fHistPiAPMonitorCuts[0]->Fill(5*fillFlagAL);
      
    //-------------- daughters --------------------//
    Int_t id0  = p0->GetDaughter(0);
    Int_t id1  = p0->GetDaughter(1);
    if(id0<0 || id1 <0) continue;
      
    fHistPiPiMonitorCuts[0]->Fill(6*fillFlagK0);
    fHistPiPMonitorCuts[0]->Fill(6*fillFlagL);
    fHistPiAPMonitorCuts[0]->Fill(6*fillFlagAL);
            
    Int_t pdgCodeD0 = stack->Particle(id0)->GetPdgCode();
    Int_t pdgCodeD1 = stack->Particle(id1)->GetPdgCode();
      
    if(pdgCodeD0 == pdgCodeD1) continue;
    if(pdgCodeD0*pdgCodeD1>0) continue;
      
    fHistPiPiMonitorCuts[0]->Fill(7*fillFlagK0);
    fHistPiPMonitorCuts[0]->Fill(7*fillFlagL);
    fHistPiAPMonitorCuts[0]->Fill(7*fillFlagAL);
            
    if((fabs(pdgCodeD0) != 211 ) && ( fabs(pdgCodeD0) != 2212 )) continue;
    if((fabs(pdgCodeD1) != 211 ) && ( fabs(pdgCodeD1) != 2212 )) continue;
      
    fHistPiPiMonitorCuts[0]->Fill(8*fillFlagK0);
    fHistPiPMonitorCuts[0]->Fill(8*fillFlagL);
    fHistPiAPMonitorCuts[0]->Fill(8*fillFlagAL);
      
    TParticle *p00 =stack->Particle(id0);
    TParticle *p01 =stack->Particle(id1);
    Double_t etaMC00   = p00->Eta();
    Double_t etaMC01   = p01->Eta();
    
    //----------- unique ID check daughters-------- //
    Int_t uniqueIDdaughter0 = p00->GetUniqueID();
    Int_t uniqueIDdaughter1 = p01->GetUniqueID();
    if (uniqueIDdaughter0 !=4 || uniqueIDdaughter1 !=4 ) continue;
      
    fHistPiPiMonitorCuts[0]->Fill(9*fillFlagK0);
    fHistPiPMonitorCuts[0]->Fill(9*fillFlagL);
    fHistPiAPMonitorCuts[0]->Fill(9*fillFlagAL);

    fHistPiPiMonitorCuts[1]->Fill(9*fillFlagK0);
    fHistPiPMonitorCuts[1]->Fill(9*fillFlagL);
    fHistPiAPMonitorCuts[1]->Fill(9*fillFlagAL);
      

    //-----------get geometric properties --------------//
    // DCA of mother to prim vertex = production vertex
    //-- primary and secondary vetex --//
    Double_t vVertex[3];
    fMCev->GetPrimaryVertex()->GetXYZ(vVertex);
    //Double_t x0=p0->Vx(),y0=p0->Vy(),z0=p0->Vz();//mother production vertex
    Double_t x=p00->Vx(),y=p00->Vy(),z=p00->Vz();//daughter vertex =V0 decay vertex

    //-- decay radii --//
    Double_t rMC2D  = sqrt(x*x+y*y);
    const  Double_t xyzMC[3] = {x,y,z};
    // Double_t rMC = p00->R();
      
    //-- phi --//
    Double_t pi = TMath::Pi();
    Double_t phi = p0->Phi();
    if(phi>pi) phi -=2*pi;

    //-------------------- V0 variables ----------------//
    Double_t rapidity = p0->Y();
    Double_t massV0MC = p0->GetMass();
    Double_t ptV0MC =  p0->Pt();
    Double_t pV0MC =  p0->P();

     
    //----------------- mother variables-----------------//
    Int_t indexMother1  = p0->GetMother(0);
    Int_t isSecd=0;
    Int_t pdgMother =0;
    // Int_t goodMother=1;
    Int_t uniqueIDmother=0;
    Double_t ptXiMother=0.0;
          
    //------check mother and fill mother histos---------//
      
    Bool_t isPrim= stack->IsPhysicalPrimary(iMc);
    if(!isPrim){
      isSecd=1;// is secondary V0s
      if(indexMother1 >-1){// && !isPrim){//secondary V0s
	//     isSecd=1;// is secondary V0s
	 
	//-- check for mother --//
	TParticle *mother = stack->Particle(indexMother1);
	if(!mother) {
	  Printf("no mother pointer!");continue;
	}
	pdgMother = mother->GetPdgCode();
	fHistPiPiMonitorCuts[1]->Fill(10*fillFlagK0);
	fHistPiPMonitorCuts[1]->Fill(10*fillFlagL);
	fHistPiAPMonitorCuts[1]->Fill(10*fillFlagAL);

	Bool_t isPrimMother= stack->IsPhysicalPrimary(indexMother1);
	if(!isPrimMother) continue;
	fHistPiPiMonitorCuts[1]->Fill(11*fillFlagK0);
	fHistPiPMonitorCuts[1]->Fill(11*fillFlagL);
	fHistPiAPMonitorCuts[1]->Fill(11*fillFlagAL);
      
	uniqueIDmother =  mother->GetUniqueID();

	if(uniqueIDmother==13){
	  //	goodMother =0;
	  continue;
	}
	fHistPiPiMonitorCuts[1]->Fill(12*fillFlagK0);
	fHistPiPMonitorCuts[1]->Fill(12*fillFlagL);
	fHistPiAPMonitorCuts[1]->Fill(12*fillFlagAL);

	//-- fill secondary V0s histos and pdg histos --// 
	ptXiMother=mother->Pt();
	 
	//-- K0s --//
	if(pdgCode==310){
	  if(fabs(pdgMother)==311 || fabs(pdgMother)==313 || fabs(pdgMother)==323 ) isSecd=0; // from K0L,  K0 and K* as primary
	  else fHistPiPiPDGCode->Fill(fabs(pdgMother));
	  //goodMother=1;
	}
	//-- Lambda --//
	if(pdgCode==3122){
	  fHistPiPPDGCode->Fill(fabs(pdgMother));
	  if (//sigma family
	      ( TMath::Abs(pdgMother) == 3112) || //sigma minus
	      ( TMath::Abs(pdgMother) == 3222) || //sigma plus
	      ( TMath::Abs(pdgMother) == 3224) || //sigma *plus
	      ( TMath::Abs(pdgMother) == 3114) || //sigma *minus
	      ( TMath::Abs(pdgMother) == 3214) || //sigma *0 counts as primary????
	      ( TMath::Abs(pdgMother) == 3212)    //sigma 0 counts as primary
	      )
	    {
	      isSecd=0;
	      //goodMother=1;
	    }
	   
	  if( TMath::Abs(pdgMother) == 3322) //xi0
	    {
	      if(!fRapCutV0 || fabs(rapidity)<fRap){
		if(!fEtaCutMCDaughters  ||  (fabs(etaMC00)<fEtaCutMCDaughtersVal|| fabs(etaMC01)<fEtaCutMCDaughtersVal)){
		  fHistPiPMassVSPtSecXi[0]->Fill(massV0MC,ptV0MC);
		  fHistPiPMassVSYSecXi[0]->Fill(massV0MC,ptV0MC);
		  fHistPiPXi0PtVSLambdaPt[0]->Fill(ptV0MC,ptXiMother);
		}
	      }
	      //goodMother=1;
	    }
	  if( TMath::Abs(pdgMother) == 3312) //xi minus
	    {
	      if(!fRapCutV0 || fabs(rapidity)<fRap){
		if(!fEtaCutMCDaughters  ||  (fabs(etaMC00)<fEtaCutMCDaughtersVal|| fabs(etaMC01)<fEtaCutMCDaughtersVal)){
		  fHistPiPMassVSPtSecXi[0]->Fill(massV0MC,ptV0MC);
		  fHistPiPMassVSYSecXi[0]->Fill(massV0MC,ptV0MC);
		  fHistPiPXiMinusPtVSLambdaPt[0]->Fill(ptV0MC,ptXiMother);
		}
	      }
	      //goodMother=1;

	    }	    
	  if(TMath::Abs(pdgMother) == 3334)//omega-
	    {
	      //  fHistPiPDCAtoPrimVtxOmega[0]->Fill(p0->GetMass(),dcaV0ToPrimVertex);
	      //goodMother=1;
	    }
	}
	//-- AntiLambda --//
	if(pdgCode==-3122 ){
	  fHistPiAPPDGCode->Fill(fabs(pdgMother));
	  if (//sigma family
	      ( TMath::Abs(pdgMother) == 3112) ||//sigma minus
	      ( TMath::Abs(pdgMother) == 3222) ||//sigma plus
	      ( TMath::Abs(pdgMother) == 3224) ||//sigma *plus
	      ( TMath::Abs(pdgMother) == 3114) ||//sigma *minus
	      ( TMath::Abs(pdgMother) == 3214) || //sigma *0
	      ( TMath::Abs(pdgMother) == 3212)    //sigma 0 counts as primary
	      )
	    {
	      isSecd=0;
	      //	    goodMother=1;
	    }
		   
	  if( TMath::Abs(pdgMother) == 3322) //xi0
	    {
	      if(!fRapCutV0 || fabs(rapidity)<fRap){
		if(!fEtaCutMCDaughters  ||  (fabs(etaMC00)<fEtaCutMCDaughtersVal|| fabs(etaMC01)<fEtaCutMCDaughtersVal)){
		  fHistPiAPMassVSPtSecXi[0]->Fill(massV0MC,ptV0MC);
		  fHistPiAPMassVSYSecXi[0]->Fill(massV0MC,ptV0MC);
		  fHistPiAPXi0PtVSLambdaPt[0]->Fill(ptV0MC,ptXiMother);
		}
	      }
	      // goodMother=1;
	    }
	  if( TMath::Abs(pdgMother) == 3312) //xi minus
	    {
	      if(!fRapCutV0 || fabs(rapidity)<fRap){
		if(!fEtaCutMCDaughters  ||  (fabs(etaMC00)<fEtaCutMCDaughtersVal|| fabs(etaMC01)<fEtaCutMCDaughtersVal)){
		  fHistPiAPMassVSPtSecXi[0]->Fill(massV0MC,ptV0MC);
		  fHistPiAPMassVSYSecXi[0]->Fill(massV0MC,ptV0MC);
		  fHistPiAPXiMinusPtVSLambdaPt[0]->Fill(ptV0MC,ptXiMother);
		}
	      }
	      //goodMother=1;
	    }
	  if(
	     (TMath::Abs(pdgMother) == 3334)//omega-
	     )
	    {
	      // fHistPiAPDCAtoPrimVtxOmega[0]->Fill(p0->GetMass(),dcaV0ToPrimVertex);
	      // goodMother=1;
	    }
	}
      }	
    }//end secondaries
    //else goodMother=1;
    // if(!goodMother) continue;// for (A)Lambda important
      
    //-------------- MC truth or reco mode -----------------//
    if(fMCTruthMode && !fMCMode){//MC true ana
      fHistPiPiMonitorCuts[isSecd]->Fill(13*fillFlagK0);
      fHistPiPMonitorCuts[isSecd]->Fill(13*fillFlagL);
      fHistPiAPMonitorCuts[isSecd]->Fill(13*fillFlagAL);
      
      //-- DCA daughters --//
      // values of one daugher, should be the same      
      /*
      //to primary vertex
      trackPos->GetImpactParameters(tdcaPosToVertex[0],tdcaPosToVertex[1]);
      trackNeg->GetImpactParameters(tdcaNegToVertex[0],tdcaNegToVertex[1]);
	 
      Double_t dcaPosToVertex = TMath::Sqrt(tdcaPosToVertex[0]*tdcaPosToVertex[0]+tdcaPosToVertex[1]*tdcaPosToVertex[1]);
      Double_t dcaNegToVertex = TMath::Sqrt(tdcaNegToVertex[0]*tdcaNegToVertex[0]+tdcaNegToVertex[1]*tdcaNegToVertex[1]);
      fHistDCADaughtersToPrimVtx[isSecd]->Fill(dcaPosToVertex,dcaNegToVertex);
      */
	 
	 
      //-- armenteros values --//
      TVector3 vecPip;
      TVector3 vecPin;
	
      Double_t ptPlus=0, ptMinus=0;
      Double_t pt00 = p00->Pt();
      Double_t pt01 = p01->Pt();
	 
      if(p00->GetPdgCode()<0)
	{
	  vecPip.SetXYZ(p01->Px(),p01->Py(),p01->Pz());
	  vecPin.SetXYZ(p00->Px(),p00->Py(),p00->Pz());
	  ptMinus = pt00;
	  ptPlus = pt01;
	}
      else{
	vecPin.SetXYZ(p01->Px(),p01->Py(),p01->Pz());
	vecPip.SetXYZ(p00->Px(),p00->Py(),p00->Pz());
	ptMinus = pt01;
	ptPlus = pt00;
      }
      TVector3 momTot(p0->Px(),p0->Py(),p0->Pz());
      Double_t lQlNeg = fabs(vecPin.Dot(momTot)/momTot.Mag());
      Double_t lQlPos = fabs(vecPip.Dot(momTot)/momTot.Mag());
      Double_t alfa =0.0;
      Double_t den = lQlPos + lQlNeg;
      if(den>0) alfa = (lQlPos - lQlNeg)/den;
      TVector3 qtvec= vecPin.Cross(momTot);//vecPip.Mag()*sqrt(1-pow(thetapip,2));
      Float_t qt = qtvec.Mag()/momTot.Mag();
      
      //-- check for injejcted --//
      Bool_t injected = kTRUE;
      injected = fMCev->IsFromBGEvent(iMc);

      if(fSelectInjected && injected ) continue;
      fHistPiPiMonitorCuts[isSecd]->Fill(14*fillFlagK0);
      fHistPiPMonitorCuts[isSecd]->Fill(14*fillFlagL);
      fHistPiAPMonitorCuts[isSecd]->Fill(14*fillFlagAL);
	 
      if(pdgCode == 310) {
	fHistPiPiEtaDMC[isSecd]->Fill(etaMC00,ptV0MC);
	fHistPiPiEtaDMC[isSecd]->Fill(etaMC01,ptV0MC);
      }
      if(fabs(pdgCode) == 3122) {
	fHistPiPEtaDMC[isSecd]->Fill(etaMC00,ptV0MC);
	fHistPiPEtaDMC[isSecd]->Fill(etaMC01,ptV0MC);
      }

      //-- rapidity and eta cut --//      
      if(fRapCutV0 && fabs(rapidity)>fRap) continue;
      fHistPiPiMonitorCuts[isSecd]->Fill(15*fillFlagK0);
      fHistPiPMonitorCuts[isSecd]->Fill(15*fillFlagL);
      fHistPiAPMonitorCuts[isSecd]->Fill(15*fillFlagAL);
	 
      if(fEtaCutMCDaughters) { if(fabs(etaMC00)>fEtaCutMCDaughtersVal || fabs(etaMC01)>fEtaCutMCDaughtersVal ) continue; }
      fHistPiPiMonitorCuts[isSecd]->Fill(16*fillFlagK0);
      fHistPiPMonitorCuts[isSecd]->Fill(16*fillFlagL);
      fHistPiAPMonitorCuts[isSecd]->Fill(16*fillFlagAL);

      Double_t declength=p00->Rho();      

      //-- Fill Particle histos --//
      if (pdgCode==310){//K0s
	fHistPiPiMonitorCuts[isSecd]->Fill(17);

	fHistPiPiEtaDMC[1]->Fill(etaMC00,ptV0MC);
	fHistPiPiEtaDMC[1]->Fill(etaMC01,ptV0MC);
	  
	fHistPiPiMass[isSecd]->Fill(massV0MC);
	fHistPiPiMassVSPt[isSecd]->Fill(massV0MC,ptV0MC);
	// fHistPiPiPtDaughters[isSecd]->Fill(ptMinus,ptPlus);
	fHistPiPiPtVSY[isSecd]->Fill(rapidity,ptV0MC);
	   
	Double_t ctK0s=0.0;
	if(pV0MC>0.0) ctK0s=declength*0.497614/pV0MC;
	fHistPiPiDecayLengthVsPt[isSecd]->Fill(ptV0MC,ctK0s);
	
	//all V0s histo
	fHistArmenteros[isSecd]->Fill(alfa,qt);
	fHistV0RadiusZ[isSecd]->Fill(rMC2D,xyzMC[2]);
	fHistV0RadiusXY[isSecd]->Fill(xyzMC[0],xyzMC[1]);
	fHistV0RadiusZVSPt[isSecd]->Fill(ptV0MC,xyzMC[2]);
	    
      }
      else if (pdgCode==3122){ //Lambda
	fHistPiPMonitorCuts[isSecd]->Fill(17);
	
	fHistPiPEtaDMC[1]->Fill(etaMC00,ptV0MC);
	fHistPiPEtaDMC[1]->Fill(etaMC01,ptV0MC);

	fHistPiPMassVSPt[isSecd]->Fill(massV0MC,ptV0MC);
	fHistPiPMass[isSecd]->Fill(massV0MC);  
	//  fHistPiPPtDaughters[isSecd]->Fill(ptMinus,ptPlus);
	fHistPiPPtVSY[isSecd]->Fill(rapidity,ptV0MC);
	   
	Double_t ctL=0.0;
	if(pV0MC>0.0) ctL=declength*1.115683/pV0MC;
	fHistPiPDecayLengthVsPt[isSecd]->Fill(ptV0MC,ctL);

	//all V0s hito	
	fHistArmenteros[isSecd]->Fill(alfa,qt);
	fHistV0RadiusZ[isSecd]->Fill(rMC2D,xyzMC[2]);
	fHistV0RadiusXY[isSecd]->Fill(xyzMC[0],xyzMC[1]);
	fHistV0RadiusZVSPt[isSecd]->Fill(ptV0MC,xyzMC[2]);
      }
      else if (pdgCode==-3122){ //AntiLambda
	fHistPiAPMonitorCuts[isSecd]->Fill(17);
	    
	fHistPiPEtaDMC[1]->Fill(etaMC00,ptV0MC);
	fHistPiPEtaDMC[1]->Fill(etaMC01,ptV0MC);

	fHistPiAPMassVSPt[isSecd]->Fill(massV0MC,ptV0MC);
	fHistPiAPMass[isSecd]->Fill(massV0MC);
	//  fHistPiAPPtDaughters[isSecd]->Fill(ptMinus,ptPlus);
	fHistPiAPPtVSY[isSecd]->Fill(rapidity,ptV0MC);
	  
	Double_t ctAL=0.0;
	if(pV0MC>0.0) ctAL=declength*1.115683/pV0MC;
	fHistPiAPDecayLengthVsPt[isSecd]->Fill(ptV0MC,ctAL);

	//all V0s histo	   
	fHistArmenteros[isSecd]->Fill(alfa,qt);
	fHistV0RadiusZ[isSecd]->Fill(rMC2D,xyzMC[2]);
	fHistV0RadiusXY[isSecd]->Fill(xyzMC[0],xyzMC[1]);
	fHistV0RadiusZVSPt[isSecd]->Fill(ptV0MC,xyzMC[2]);
	fHistV0RadiusXYVSY[isSecd]->Fill(rapidity,rMC2D);
      }
    }//MC true ana
    else{// V0 reco ana
      V0RecoLoop(id0,id1,isSecd,pdgCode,ptV0MC,pdgMother,ptXiMother);
    }
      
  }//end MC stack loop
   
}
//________________________________________________________________________
void AliAnalysisTaskV0ForRAA::V0RecoLoop(Int_t id0,Int_t id1,Int_t isSecd,Int_t what,Double_t ptV0MC, Int_t pdgMother,Double_t ptXiMother){
  //loop over reconstructed particles

   
  //--------------------- define variables -----------------------//
  Double_t pp[3];
  Double_t pm[3];
  Double_t xr[3];
   
  Double_t massPi=0.13957018;
  Double_t massP=0.93827203;
  
  TLorentzVector positivesMIP;
  TLorentzVector negativesMIAP;
  TLorentzVector positivesMIPi;
  TLorentzVector negativesMIPi;

  /*
    AliKFParticle::SetField(fESD->GetMagneticField());
    AliKFVertex primVtx(*(fESD->GetPrimaryVertex()));
    AliKFVertex primVtxImproved = primVtx;

    AliKFParticle* negPiKF=NULL;
    AliKFParticle* posPiKF=NULL;
    AliKFParticle* posPKF=NULL;
    AliKFParticle* negAPKF=NULL;
  */

  AliESDtrack* trackPos=NULL;
  AliESDtrack* trackNeg=NULL;
  AliESDtrack* trackPosTest = NULL;
  AliESDtrack* trackNegTest =NULL;

  Float_t nSigmaTPCtrackPosToPion;
  Float_t nSigmaTPCtrackNegToPion;
  Float_t nSigmaTPCtrackPosToProton;
  Float_t nSigmaTPCtrackNegToProton;

  Double_t primaryVtxPosition[3];
  primaryVtxPosition[0] = fESD->GetPrimaryVertex()->GetXv();
  primaryVtxPosition[1] = fESD->GetPrimaryVertex()->GetYv();
  primaryVtxPosition[2] = fESD->GetPrimaryVertex()->GetZv();
   
  Int_t nV0 = fESD->GetNumberOfV0s();
  AliESDv0 * v0MIs=NULL;
   
  //------------------------ V0 reco loop --------------------//
  for(Int_t iV0MI = 0; iV0MI < nV0; iV0MI++) {//V0 loop
    //-- get V0 info --//
    v0MIs = fESD->GetV0(iV0MI);

    fHistPiPiMonitorCuts[isSecd]->Fill(1);
    fHistPiPMonitorCuts[isSecd]->Fill(1);
    fHistPiAPMonitorCuts[isSecd]->Fill(1);

    //------------ get references of daughters --------------//
    //-- esd tracks --//
    trackPosTest = fESD->GetTrack(v0MIs->GetPindex());
    trackNegTest = fESD->GetTrack(v0MIs->GetNindex());
     
    if ( trackPosTest->GetSign() == trackNegTest->GetSign()) continue;
	 
    fHistPiPiMonitorCuts[isSecd]->Fill(2);
    fHistPiPMonitorCuts[isSecd]->Fill(2);
    fHistPiAPMonitorCuts[isSecd]->Fill(2);
      
    //-- for MC mode --//
    if(fMCMode){
      //check MC labels (and find partners for MC truth V0 daughters for fMCTruthMode=kTRUE)
      if(!GetMCTruthPartner(trackPosTest,trackNegTest,id0,id1)) continue;
    }

    fHistPiPiMonitorCuts[isSecd]->Fill(3);
    fHistPiPMonitorCuts[isSecd]->Fill(3);
    fHistPiAPMonitorCuts[isSecd]->Fill(3);

    //-- onthefly selection --//
    Bool_t onthefly = v0MIs->GetOnFlyStatus();
    if(fOntheFly!=onthefly) continue;

    fHistPiPiMonitorCuts[isSecd]->Fill(4);
    fHistPiPMonitorCuts[isSecd]->Fill(4);
    fHistPiAPMonitorCuts[isSecd]->Fill(4);
         
    //--  get eta from V0 daughters --//
    Double_t posDaughterEta=0.0;
    Double_t negDaughterEta=0.0;
      
    Double_t eta00 = trackPosTest->Eta();
    Double_t eta01 = trackNegTest->Eta();
            
    //---------- check sign assignment for daughters --------//
    Bool_t switchSign = kFALSE;
      
    if( trackPosTest->GetSign() >0){//pos
      trackPos =fESD->GetTrack(v0MIs->GetPindex());
      trackNeg =fESD->GetTrack(v0MIs->GetNindex());

      v0MIs->GetPPxPyPz(pp[0],pp[1],pp[2]);
      v0MIs->GetNPxPyPz(pm[0],pm[1],pm[2]);

      posDaughterEta =v0MIs->GetParamP()->Eta();
      negDaughterEta=v0MIs->GetParamN()->Eta();
      
      /*      
	      if (negPiKF) delete negPiKF; negPiKF=NULL;
	      if (posPiKF) delete posPiKF; posPiKF=NULL;
	      if (posPKF) delete posPKF; posPKF=NULL;
	      if (negAPKF) delete negAPKF; negAPKF=NULL;

	      negPiKF = new AliKFParticle( *(v0MIs->GetParamN()) ,-211);
	      posPiKF = new AliKFParticle( *(v0MIs->GetParamP()) ,211);
	      posPKF = new AliKFParticle( *(v0MIs->GetParamP()) ,2212);
	      negAPKF = new AliKFParticle( *(v0MIs->GetParamN()) ,-2212);
      */
    }
    if( trackPosTest->GetSign() <0){//neg
      trackPos =fESD->GetTrack(v0MIs->GetNindex());
      trackNeg =fESD->GetTrack(v0MIs->GetPindex());
     
      v0MIs->GetNPxPyPz(pp[0],pp[1],pp[2]);
      v0MIs->GetPPxPyPz(pm[0],pm[1],pm[2]);
      
      posDaughterEta = v0MIs->GetParamN()->Eta();
      negDaughterEta = v0MIs->GetParamP()->Eta();
	 
      /*
	if (negPiKF) delete negPiKF; negPiKF=NULL;
	if (posPiKF) delete posPiKF; posPiKF=NULL;
	if (posPKF) delete posPKF; posPKF=NULL;
	if (negAPKF) delete negAPKF; negAPKF=NULL;


	negPiKF = new AliKFParticle( *(v0MIs->GetParamP()) ,-211);
	posPiKF = new AliKFParticle( *(v0MIs->GetParamN()) ,211);
	posPKF = new AliKFParticle( *(v0MIs->GetParamN()) ,2212);
	negAPKF = new AliKFParticle( *(v0MIs->GetParamP()) ,-2212);
      */
      switchSign = kTRUE;
      eta01 = trackPosTest->Eta();
      eta00 = trackNegTest->Eta();

    }
    
    // ------------- calc masses and 4 vectors -------------- //
    //K0
    positivesMIPi.SetXYZM(pp[0],pp[1],pp[2],massPi);
    negativesMIPi.SetXYZM(pm[0],pm[1],pm[2],massPi);
    TLorentzVector v0K0=positivesMIPi+negativesMIPi;
    
    //Lambda
    positivesMIP.SetXYZM(pp[0],pp[1],pp[2],massP);
    TLorentzVector v0Lambda=positivesMIP+negativesMIPi;
     
    //Anitlambda
    negativesMIAP.SetXYZM(pm[0],pm[1],pm[2],massP);
    TLorentzVector v0ALambda=positivesMIPi+negativesMIAP;
     
    //---------------------AliKFParticle ---------------------//
    /*  
	Double_t chi2K0C=0.0;
	Double_t chi2LambdaC=0.0;
	Double_t chi2ALambdaC=0.0;

     
	AliKFParticle v0K0KF;
	v0K0KF +=(*negPiKF);
	v0K0KF +=(*posPiKF);
	//v0K0C.SetVtxGuess(xr[0],xr[1],xr[2]);
	v0K0KF.SetProductionVertex(primVtxImproved);
	  
	AliKFParticle v0LambdaKF;
	v0LambdaKF +=(*negPiKF);
	v0LambdaKF +=(*posPKF);
	//v0LambdaC.SetVtxGuess(xr[0],xr[1],xr[2]);
	v0LambdaKF.SetProductionVertex(primVtxImproved);
	  
	AliKFParticle v0ALambdaKF;
	v0ALambdaKF +=(*negAPKF);
	v0ALambdaKF +=(*posPiKF);
	//v0ALambdaC.SetVtxGuess(xr[0],xr[1],xr[2]);
	v0ALambdaKF.SetProductionVertex(primVtxImproved);
	
	if( v0K0KF.GetNDF() != 0) {
	chi2K0C = v0K0KF.GetChi2()/v0K0KF.GetNDF();
	}

	Double_t chi2LambdaC=100000.;
	if( v0LambdaKF.GetNDF() != 0) {
	chi2LambdaC = v0LambdaKF.GetChi2()/v0LambdaKF.GetNDF();
	}

	Double_t chi2ALambdaC=100000.;
	if( v0ALambdaKF.GetNDF() != 0) {
	chi2ALambdaC = v0ALambdaKF.GetChi2()/v0ALambdaKF.GetNDF();
	}
    */
    
    // ----------------- for MC mode ------------------------ //
    Bool_t fillK0sMC = kTRUE;
    Bool_t fillLambdaMC = kTRUE;
    Bool_t fillALambdaMC = kTRUE;

    if(fMCMode && fMCTruthMode) {
      if(what == 310) {
	fillLambdaMC = kFALSE;
	fillALambdaMC = kFALSE;
      }
      else if(what == 3122){
	fillALambdaMC = kFALSE;
	fillK0sMC = kFALSE;
      }
      else if(what == -3122){
	fillLambdaMC = kFALSE;
	fillK0sMC = kFALSE;
      }
    }
   
    //----------------- prepare for V0 ana ------------------//
    TVector3 ppTrack(pp);
    TVector3 pmTrack(pm);
      
    Double_t posDaughterPt = ppTrack.Pt();
    Double_t negDaughterPt = pmTrack.Pt();
    /*
      Double_t v0sPt=v0MIs->Pt();
      if(what == 310 || what ==0){
      fHistPiPiEtaDReco[0]->Fill(posDaughterPt,v0sPt);
      fHistPiPiEtaDReco[0]->Fill(negDaughterPt,v0sPt);
      }
      if(fabs(what) == 3122 || what == 0){
      fHistPiPEtaDReco[0]->Fill(posDaughterPt,v0sPt);
      fHistPiPEtaDReco[0]->Fill(negDaughterPt,v0sPt);
      }
    */
    //------------------- DCA daughters ---------------------//
    //-- between the daughters --//
    Double_t dcaDaughters = v0MIs->GetDcaV0Daughters();  
    
    //-- to primary vertex --
    Double_t dcaPosToVertex=0.0,dcaNegToVertex=0.0;
    Double_t dzPos=(primaryVtxPosition[0]-xr[0])*ppTrack.Y() - (primaryVtxPosition[1]-xr[1])*ppTrack.X();
    dcaPosToVertex=TMath::Sqrt(dzPos*dzPos/(pow(ppTrack.X(),2)+pow(ppTrack.Y(),2)));
    Double_t dzNeg=(primaryVtxPosition[0]-xr[0])*pmTrack.Y() - (primaryVtxPosition[1]-xr[1])*pmTrack.X();
    dcaNegToVertex=TMath::Sqrt(dzNeg*dzNeg/(pow(pmTrack.X(),2)+pow(pmTrack.Y(),2)));
     
    // Double_t dcaPosToVertex[3];dcaNegToVertex[3];
    // trackPos->GetImpactParameters(dcaPosToVertex[0],dcaPosToVertex[1]);
    // trackNeg->GetImpactParameters(dcaNegToVertex[0],dcaNegToVertex[1]);
	 
    // dcaPosToVertex = TMath::Sqrt(dcaPosToVertex[0]*dcaPosToVertex[0]+dcaPosToVertex[1]*dcaPosToVertex[1]);
    // dcaNegToVertex = TMath::Sqrt(dcaNegToVertex[0]*dcaNegToVertex[0]+dcaNegToVertex[1]*dcaNegToVertex[1]);
	 
    //dcaPosToVertex  =   posPKF->GetDistanceFromVertexXY(primaryVtxPosition);
    //dcaNegToVertex  =   negPiKF->GetDistanceFromVertexXY(primaryVtxPosition);
        
    //------------------- dca and decay radius V0 -------------//
    v0MIs->GetXYZ(xr[0],xr[1],xr[2]);
       
    Double_t dim2V0Radius      = TMath::Sqrt(xr[0]*xr[0]+xr[1]*xr[1]);

    Double_t decayLength =
      TMath::Sqrt(TMath::Power(xr[0] - primaryVtxPosition[0],2) +
		  TMath::Power(xr[1] - primaryVtxPosition[1],2) +
		  TMath::Power(xr[2] - primaryVtxPosition[2],2 ));

    Double_t  dcaV0ToPrimVertex= v0MIs->GetD(primaryVtxPosition[0],primaryVtxPosition[1]);////v0K0KF.GetDistanceFromVertexXY(tPrimaryVtxPosition);

    //-------------------- general cuts -------------------//
    //-- esd track cuts for daughters --//
    if(fESDTrackCuts){
      if(!fESDTrackCuts->AcceptTrack(trackPosTest)) continue;
      if(!fESDTrackCuts->AcceptTrack(trackNegTest)) continue;
    }
      
    fHistPiPiMonitorCuts[isSecd]->Fill(5);
    fHistPiPMonitorCuts[isSecd]->Fill(5);
    fHistPiAPMonitorCuts[isSecd]->Fill(5);
      
    //-- eta cut --//
    if( fabs(posDaughterEta) > fEtaCutMCDaughtersVal || fabs(negDaughterEta) > fEtaCutMCDaughtersVal) continue;
    fHistPiPiMonitorCuts[isSecd]->Fill(6);
    fHistPiPMonitorCuts[isSecd]->Fill(6);
    fHistPiAPMonitorCuts[isSecd]->Fill(6);

    //-- pt cut --//
    if( fabs(posDaughterPt)<fMinPt || fabs(negDaughterPt) < fMinPt ) continue;
    fHistPiPiMonitorCuts[isSecd]->Fill(7);
    fHistPiPMonitorCuts[isSecd]->Fill(7);
    fHistPiAPMonitorCuts[isSecd]->Fill(7);

      
    //-- radius xy min cut --//
    if(dim2V0Radius < fDecayRadXYMin) continue;
    //	    if(fabs(xr[1])<fDecayRadY) continue;
    fHistPiPiMonitorCuts[isSecd]->Fill(8);
    fHistPiPMonitorCuts[isSecd]->Fill(8);
    fHistPiAPMonitorCuts[isSecd]->Fill(8);

    //-- radius xy max cut --//
    if(dim2V0Radius > fDecayRadXYMax) continue;
    //	    if(fabs(xr[1])<fDecayRadY) continue;
    fHistPiPiMonitorCuts[isSecd]->Fill(9);
    fHistPiPMonitorCuts[isSecd]->Fill(9);
    fHistPiAPMonitorCuts[isSecd]->Fill(9);
      
    //-- decay length min ->ctau --//
    //  if(decayLength<fDecayLengthMin) continue;
    //       fHistPiPiMonitorCuts[isSecd]->Fill(11);
    //       fHistPiPMonitorCuts[isSecd]->Fill(11);
    //       fHistPiAPMonitorCuts[isSecd]->Fill(11);

    //-- decay length max cut --//
    if(decayLength>fDecayLengthMax) continue;
    fHistPiPiMonitorCuts[isSecd]->Fill(10);
    fHistPiPMonitorCuts[isSecd]->Fill(10);
    fHistPiAPMonitorCuts[isSecd]->Fill(10);
   
    //----------------------- V0 variables --------------------//
    //-- armenteros --//
    TVector3 momTot = ppTrack + pmTrack;
    Double_t lQlNeg = fabs(pmTrack.Dot(momTot)/momTot.Mag());
    Double_t lQlPos = fabs(ppTrack.Dot(momTot)/momTot.Mag());
    //return 1.-2./(1.+lQlNeg/lQlPos);
    Double_t alfa =0.0;
    Double_t den = lQlPos + lQlNeg;
    if(den>0) alfa = (lQlPos - lQlNeg)/den;
    TVector3 qtvec= pmTrack.Cross(momTot);//vecPip.Mag()*sqrt(1-pow(thetapip,2));
    Double_t qt = qtvec.Mag()/momTot.Mag();

    //-- momenta --//
    Double_t ptK0s = v0K0.Pt();
    Double_t ptLambda = v0Lambda.Pt();
    Double_t ptALambda = v0ALambda.Pt();
      
    Double_t pK0s = v0K0.P();
    Double_t pLambda = v0Lambda.P();
    Double_t pALambda = v0ALambda.P();
      
    //-- masses --//
    Double_t massK0s = v0K0.M();
    Double_t massLambda = v0Lambda.M();
    Double_t massALambda = v0ALambda.M();

    //-- rapidity --//
    Double_t rapK0s = v0MIs->Y(310);
    Double_t rapL   = v0MIs->Y(3122);
    Double_t rapAL  = v0MIs->Y(3122);

    //-- other variables --//
    Double_t opAng =   fabs(ppTrack.Angle(pmTrack));
    //    if( ppTrack.Angle(pmTrack)<0.001) continue;  
    //    if( ppTrack.Angle(pmTrack)<0.004) continue;   
    
    Double_t cosOPAng = v0MIs->GetV0CosineOfPointingAngle();
    
    /*     
    //introduce more histo
    Double_t errOnMassK0s = v0MIs->ChangeMassHypothesis(310);
    Double_t errOnMassLambda = 0.0;
    Double_t errOnMassALambda = 0.0;
    if(!switchSign){
    errOnMassLambda  = v0MIs->ChangeMassHypothesis(3122);
    errOnMassALambda = v0MIs->ChangeMassHypothesis(-3122);
    }
    else{
    errOnMassLambda  = v0MIs->ChangeMassHypothesis(-3122);
    errOnMassALambda = v0MIs->ChangeMassHypothesis(3122);
    }
    */

      
    //------------------------ detector values -------------------------------//
    //-- TPC ITS values pos --//
    Int_t nclsTPCPos =  trackPos->GetNcls(1);
    Int_t nclsTPCFindablePos =  trackPos->GetTPCNclsF();
    Int_t nclsITSPos =  trackPos->GetNcls(0);
    Double_t chi2PerClusterITSPos = -1.0;
    if(nclsITSPos>0) chi2PerClusterITSPos = trackPos->GetITSchi2()/Double_t(nclsITSPos);
    Double_t crossedRowsTPCPos = trackPos->GetTPCCrossedRows();
      
    //-- TPC ITS values neg --//
    Int_t nclsTPCNeg =  trackNeg->GetNcls(1);
    Int_t nclsTPCFindableNeg =  trackNeg->GetTPCNclsF();
    Int_t nclsITSNeg =  trackNeg->GetNcls(0);
    Double_t chi2PerClusterITSNeg = -1.0;
    if(nclsITSNeg>0) chi2PerClusterITSNeg =trackNeg->GetITSchi2()/Double_t(nclsITSNeg);
    Double_t crossedRowsTPCNeg = trackNeg->GetTPCCrossedRows();    

    if(fMoreNclsThanRows && (crossedRowsTPCPos < nclsTPCPos || crossedRowsTPCNeg < nclsTPCNeg  )) continue;
    fHistPiPiMonitorCuts[isSecd]->Fill(11);
    fHistPiPMonitorCuts[isSecd]->Fill(11);
    fHistPiAPMonitorCuts[isSecd]->Fill(11);
      
    if(fMoreNclsThanFindable && (nclsTPCFindablePos < nclsTPCPos || nclsTPCFindableNeg < nclsTPCNeg  )) continue;
    fHistPiPiMonitorCuts[isSecd]->Fill(12);
    fHistPiPMonitorCuts[isSecd]->Fill(12);
    fHistPiAPMonitorCuts[isSecd]->Fill(12);      

    if(chi2PerClusterITSNeg > fChi2PerClusterITS || chi2PerClusterITSPos > fChi2PerClusterITS ) continue;
    fHistPiPiMonitorCuts[isSecd]->Fill(13);
    fHistPiPMonitorCuts[isSecd]->Fill(13);
    fHistPiAPMonitorCuts[isSecd]->Fill(13);  
      
    //-- cut flags for V0 specific cuts --//
    Bool_t cutOKK0s = kTRUE;
    Bool_t cutOKLambda = kTRUE;
    Bool_t cutOKALambda = kTRUE;
      
    //-------------------------- K0 cuts -----------------------------//

    if(dcaV0ToPrimVertex > fDCAToVertexK0) continue;
    if(fabs(xr[2])>fDCAZ) continue;//like decay radius z component
    fHistPiPiMonitorCuts[isSecd]->Fill(14);
      
    Double_t ctK0 = 0.0;
    if(fabs(pK0s)>0.0)  ctK0 = decayLength*0.497614/pK0s;
    if( ctK0 > fCtauK0s &&  fabs(ptK0s) <fCtauPtCutK0) cutOKK0s = kFALSE;
    else  fHistPiPiMonitorCuts[isSecd]->Fill(15);
      
    if((cosOPAng < fCosPointAngK && fabs(ptK0s) < fCPAPtCutK0) || cosOPAng<0.99) cutOKK0s = kFALSE;
    else  fHistPiPiMonitorCuts[isSecd]->Fill(16);

    if(dcaDaughters > fDCADaughtersK0 )cutOKK0s = kFALSE;
    else  fHistPiPiMonitorCuts[isSecd]->Fill(17);
	 
    if(dcaNegToVertex < fDCADaughtersToVtxSmall || dcaPosToVertex < fDCADaughtersToVtxSmall)  cutOKK0s = kFALSE;
    else  fHistPiPiMonitorCuts[isSecd]->Fill(18);

    if(fRapCutV0 && fabs(rapK0s) > fRap) cutOKK0s = kFALSE;
    else  fHistPiPiMonitorCuts[isSecd]->Fill(19);  
	
    // if(chi2K0C > fChiCutKf) cutOKK0s = kFALSE;
    if(opAng < fOpengAngleDaughters && fabs(ptK0s) < fOpAngPtCut )  cutOKK0s = kFALSE;
    else  fHistPiPiMonitorCuts[isSecd]->Fill(20);

    Bool_t ptbinokK0s=kFALSE;
    if( ptK0s < fQtCut &&  ptK0s > fQtCutPtLow ) ptbinokK0s=kTRUE;
    Double_t qtval  = 0.2*fabs(alfa);
    if(fArmCutK0 && ptbinokK0s && qt < qtval) cutOKK0s = kFALSE;
    else  fHistPiPiMonitorCuts[isSecd]->Fill(21);
       
    //-------------------------- Lambda cuts -------------------------//
    if(dcaV0ToPrimVertex > fDCAToVertexL) continue;
    if(fabs(xr[2])>fDCAZ) continue;//like decay radius z component
    fHistPiPMonitorCuts[isSecd]->Fill(14);
         
    Double_t ctL = 0.0;
    if(fabs(pLambda)>0.0)  ctL = decayLength*1.115683/fabs(pLambda);
    if(ctL > fCtauL && fabs(ptLambda) <fCtauPtCutL)  cutOKLambda = kFALSE;
    else  fHistPiPMonitorCuts[isSecd]->Fill(15);
      
    if((cosOPAng<fCosPointAngL && fabs(ptLambda) < fCPAPtCutL)|| cosOPAng<0.99)  cutOKLambda = kFALSE;
    else fHistPiPMonitorCuts[isSecd]->Fill(16);

    if(dcaDaughters > fDCADaughtersL )cutOKLambda = kFALSE;
    else  fHistPiPMonitorCuts[isSecd]->Fill(17);
 
    if( dcaNegToVertex < fDCADaughtersToVtxSmall || dcaPosToVertex < fDCADaughtersToVtxLarge)  cutOKLambda = kFALSE;
    else  fHistPiPMonitorCuts[isSecd]->Fill(18);

    if(fRapCutV0 && fabs(rapL) > fRap) cutOKLambda = kFALSE;
    else  fHistPiPMonitorCuts[isSecd]->Fill(19);
    /*	 
	 if(chi2LambdaC > fChiCutKf) cutOKLambda = kFALSE;
	 else  fHistPiPMonitorCuts[isSecd]->Fill(20);
    */
    fHistPiPMonitorCuts[isSecd]->Fill(20);
      
    if(alfa<fAlfaCut  || (fArmCutL && qt>fQtCut)) cutOKLambda = kFALSE;
    else  fHistPiPMonitorCuts[isSecd]->Fill(21);
      
    //--------------------------- ALambda cuts --------------------------//

    if(dcaV0ToPrimVertex > fDCAToVertexL) continue;
    if(fabs(xr[2])>fDCAZ) continue;//like decay radius z component
    fHistPiAPMonitorCuts[isSecd]->Fill(14);
      
    Double_t ctAL = 0.0;
    if(fabs(pALambda)>0.0)  ctAL = decayLength*1.115683/fabs(pALambda);
    if(ctAL > fCtauL &&  fabs(ptALambda) <fCtauPtCutL)  cutOKALambda = kFALSE;
    else  fHistPiAPMonitorCuts[isSecd]->Fill(15);

    if((cosOPAng<fCosPointAngL && fabs(ptALambda) < fCPAPtCutL)|| cosOPAng<0.99)  cutOKALambda = kFALSE;
    else fHistPiAPMonitorCuts[isSecd]->Fill(16);
      
    if(dcaDaughters > fDCADaughtersAL )cutOKALambda = kFALSE;
    else  fHistPiAPMonitorCuts[isSecd]->Fill(17);
	 
    if( dcaPosToVertex < fDCADaughtersToVtxSmall || dcaNegToVertex < fDCADaughtersToVtxLarge)  cutOKALambda = kFALSE;
    else fHistPiAPMonitorCuts[isSecd]->Fill(18);
	 
    if(fRapCutV0 && fabs(rapAL) > fRap) cutOKALambda = kFALSE;
    else fHistPiAPMonitorCuts[isSecd]->Fill(19);
    /*
      if(chi2ALambdaC > fChiCutKf) cutOKALambda = kFALSE;
      else  fHistPiAPMonitorCuts[isSecd]->Fill(20);
    */
    fHistPiAPMonitorCuts[isSecd]->Fill(20);
      
    if((fArmCutL && qt>fQtCut) || alfa > -1.0*fAlfaCut) cutOKALambda = kFALSE;
    else  fHistPiAPMonitorCuts[isSecd]->Fill(21);
      
    //--------------------- PID ----------------------------//
    //-- dEdx --//  
    nSigmaTPCtrackPosToPion = fabs(fESDpid->NumberOfSigmasTPC(trackPos,AliPID::kPion));
    nSigmaTPCtrackNegToPion = fabs(fESDpid->NumberOfSigmasTPC(trackNeg,AliPID::kPion));
    nSigmaTPCtrackPosToProton = fabs(fESDpid->NumberOfSigmasTPC(trackPos,AliPID::kProton));
    nSigmaTPCtrackNegToProton = fabs(fESDpid->NumberOfSigmasTPC(trackNeg,AliPID::kProton));

    Bool_t pipidEdx=kTRUE;
    Bool_t pipdEdx =kTRUE;
    Bool_t piapdEdx=kTRUE;

    Double_t posDaughterP = ppTrack.Mag();
    Double_t negDaughterP = pmTrack.Mag();

    Double_t tpcsigPos= trackPos->GetTPCsignal();
    Double_t tpcsigNeg= trackNeg->GetTPCsignal();

    //-- dedx cut --//
    if(fUsePID){
      // 	if(fabs(posDaughterP)<fPPIDcut && nSigmaTPCtrackPosToPion > fNSigma ){
      // 	       pipidEdx =kFALSE;//k0s
      // 	}
      if(fabs(posDaughterP)<fPPIDcut && tpcsigPos < 5.0){//no zero dedx values!
	pipidEdx =kFALSE;//k0s
	piapdEdx =kFALSE;//antilambda
      }
      if(fabs(posDaughterP)<fPPIDcut && (nSigmaTPCtrackPosToProton > fNSigma || tpcsigPos < 5.0)) pipdEdx =kFALSE;//lambda
	    
      // 	if(fabs(negDaughterP)<fPPIDcut && nSigmaTPCtrackNegToPion > fNSigma ){
      // 	    pipidEdx =kFALSE;//k0s
      // 	}
      if(fabs(negDaughterP)<fPPIDcut &&  tpcsigNeg < 5.0){//no zero dedx values!
	pipidEdx =kFALSE;//k0s
	pipdEdx =kFALSE;//lambda
      }
      if(fabs(negDaughterP)<fPPIDcut && (nSigmaTPCtrackNegToProton > fNSigma || tpcsigNeg< 5.0)) piapdEdx =kFALSE;//antilambda
    }
      
     

    //-------------------- V0 ana -------------------------//
    //-- cut flags for furhter histos--//
    Bool_t k0sOK=kFALSE;
    Bool_t lambdaOK=kFALSE;
    Bool_t alambdaOK=kFALSE;

    //--  Check for K0 --//
    if( cutOKK0s  && fillK0sMC && pipidEdx){
      fHistPiPiMonitorCuts[isSecd]->Fill(22);
      k0sOK=kTRUE;		    
      if(massK0s>0.25 && massK0s<0.75 ){
	fHistPiPiMonitorCuts[isSecd]->Fill(23);
	fHistPiPiMass[isSecd]->Fill(massK0s);
	fHistPiPiMassVSPt[isSecd]->Fill(massK0s,ptK0s);
	fHistPiPiMassVSPtMCTruth[isSecd]->Fill(massK0s,ptV0MC);
	fHistPiPiRadiusXY[isSecd]->Fill(massK0s,opAng);
	fHistPiPiCosPointAng[isSecd]->Fill(dcaV0ToPrimVertex,cosOPAng);
	fHistPiPiDecayLengthVsPt[isSecd]->Fill(ptK0s,ctK0);
	fHistPiPiDecayLengthVsMass[isSecd]->Fill(massK0s,ctK0);
	fHistPiPiDCADaughters[isSecd]->Fill(massK0s,dcaDaughters);
	fHistPiPiDCADaughterPosToPrimVtxVSMass[isSecd]->Fill(massK0s,dcaPosToVertex);
	fHistPiPiDCAVSMass[isSecd]->Fill(massK0s,dcaV0ToPrimVertex);
	// fHistPiPiPtDaughters[isSecd]->Fill(posDaughterPt,negDaughterPt);
	fHistPiPiPtVSY[isSecd]->Fill(rapK0s,ptK0s);
      }
      fHistArmenteros[isSecd]->Fill(alfa,qt);
      fHistDedxSecPiPlus[isSecd]->Fill(posDaughterP,tpcsigPos);
      fHistDedxSecPiMinus[isSecd]->Fill(negDaughterP,tpcsigNeg);

      fHistV0RadiusZ[isSecd]->Fill(dim2V0Radius,xr[2]);
      fHistV0RadiusXY[isSecd]->Fill(xr[0],xr[1]);
      fHistV0RadiusZVSPt[isSecd]->Fill(ptK0s,dim2V0Radius);
      fHistV0RadiusXYVSY[isSecd]->Fill(rapK0s,dim2V0Radius);

      //-- detector values --/
      fHistNclsITSPosK0[isSecd]->Fill(nclsITSPos);
      fHistNclsITSNegK0[isSecd]->Fill(nclsITSNeg);
      fHistNclsTPCPosK0[isSecd]->Fill(nclsTPCPos);
      fHistNclsTPCNegK0[isSecd]->Fill(nclsTPCNeg);
      fHistChi2PerNclsITSPosK0[isSecd]->Fill(chi2PerClusterITSPos);
      fHistChi2PerNclsITSNegK0[isSecd]->Fill(chi2PerClusterITSNeg);
    }
    
    //--  Check for Lambda --//
    if(cutOKLambda && fillLambdaMC && pipdEdx){
      fHistPiPMonitorCuts[isSecd]->Fill(22);
      lambdaOK=kTRUE;
      if( massLambda>1.05 && massLambda<1.25 ){
	fHistPiPMonitorCuts[isSecd]->Fill(23);
	fHistPiPMass[isSecd]->Fill(massLambda);
	fHistPiPMassVSPtMCTruth[isSecd]->Fill(massLambda,ptV0MC);
	fHistPiPMassVSPt[isSecd]->Fill(massLambda,ptLambda);
	fHistPiPRadiusXY[isSecd]->Fill(massLambda,opAng);
	fHistPiPCosPointAng[isSecd]->Fill(dcaV0ToPrimVertex,cosOPAng);
	fHistPiPPtVSY[isSecd]->Fill(rapL,ptLambda);
	//   fHistPiPPtDaughters[isSecd]->Fill(posDaughterPt,negDaughterPt);
	fHistPiPDCADaughters[isSecd]->Fill(massLambda,dcaDaughters);
	fHistPiPDCAVSMass[isSecd]->Fill(massLambda,dcaV0ToPrimVertex);
	fHistPiPDCADaughterPosToPrimVtxVSMass[isSecd]->Fill(massLambda,dcaPosToVertex);
	fHistPiPDecayLengthVsPt[isSecd]->Fill(ptLambda,ctL);
	fHistPiPDecayLengthVsMass[isSecd]->Fill(massLambda,ctL);

	//-- secondaries --//
	if(isSecd==1){
	  if(fabs(pdgMother)==3112 || fabs(pdgMother)==3114 || fabs(pdgMother)==3222 || fabs(pdgMother)==3224 || fabs(pdgMother)==3214 ){
	    fHistPiPMassVSPtSecSigma[1]->Fill(massLambda,ptLambda);
		
	  }
	  //  if(pdgMother==3334){
	  // 		  fHistPiPDCAtoPrimVtxOmega[1]->Fill(massLambda,dcaV0ToPrimVertex);
	  // 	       }
	  if(pdgMother==3322){
	    fHistPiPCosPointAngXiVsPt->Fill(ptLambda,cosOPAng);
	    fHistPiPMassVSPtSecXi[1]->Fill(massLambda,ptLambda);
	    fHistPiPMassVSYSecXi[1]->Fill(massLambda,rapL);
	    fHistPiPXi0PtVSLambdaPt[1]->Fill(ptLambda,ptXiMother);
	  }
	  if(pdgMother==3312){
	    fHistPiPCosPointAngXiVsPt->Fill(ptLambda,cosOPAng);
	    fHistPiPMassVSPtSecXi[1]->Fill(massLambda,ptLambda);
	    fHistPiPMassVSYSecXi[1]->Fill(massLambda,rapL);
	    fHistPiPXiMinusPtVSLambdaPt[1]->Fill(ptLambda,ptXiMother);
	  }
	}  
      }
      if(ptLambda>0.4) fHistArmenteros[isSecd]->Fill(alfa,qt);
      fHistV0RadiusZ[isSecd]->Fill(dim2V0Radius,xr[2]);
      fHistV0RadiusXY[isSecd]->Fill(xr[0],xr[1]);
      fHistV0RadiusZVSPt[isSecd]->Fill(ptLambda,dim2V0Radius);
      fHistV0RadiusXYVSY[isSecd]->Fill(rapL,dim2V0Radius);
      fHistDedxSecProt[isSecd]->Fill(posDaughterP,tpcsigPos);
      fHistDedxSecPiMinus[isSecd]->Fill(negDaughterP,tpcsigNeg);
	 
      //-- detector values --//
      fHistNclsITSPosL[isSecd]->Fill(nclsITSPos);
      fHistNclsITSNegL[isSecd]->Fill(nclsITSNeg);
      fHistNclsTPCPosL[isSecd]->Fill(nclsTPCPos);
      fHistNclsTPCNegL[isSecd]->Fill(nclsTPCNeg);
      fHistChi2PerNclsITSPosL[isSecd]->Fill(chi2PerClusterITSPos);
      fHistChi2PerNclsITSNegL[isSecd]->Fill(chi2PerClusterITSNeg);
    }

 
    //-- Check for AntiLambda --//
    if(cutOKALambda && fillALambdaMC && piapdEdx){
      fHistPiAPMonitorCuts[isSecd]->Fill(22);
      alambdaOK=kTRUE;
      if( massALambda>1.05 && massALambda<1.25  ){
	fHistPiAPMonitorCuts[isSecd]->Fill(23);
	fHistPiAPMass[isSecd]->Fill(massALambda);
	fHistPiAPMassVSPtMCTruth[isSecd]->Fill(massALambda,ptV0MC);
	fHistPiAPMassVSPt[isSecd]->Fill(massALambda,ptALambda);
	fHistPiAPRadiusXY[isSecd]->Fill(massALambda,opAng);
	fHistPiAPCosPointAng[isSecd]->Fill(dcaV0ToPrimVertex,cosOPAng);
	fHistPiAPPtVSY[isSecd]->Fill(rapAL,ptALambda);
	//  fHistPiAPPtDaughters[isSecd]->Fill(posDaughterPt,negDaughterPt);
	fHistPiAPDCADaughters[isSecd]->Fill(massALambda,dcaDaughters);
	fHistPiAPDCAVSMass[isSecd]->Fill(massALambda,dcaV0ToPrimVertex);
	fHistPiAPDCADaughterPosToPrimVtxVSMass[isSecd]->Fill(massALambda,dcaPosToVertex);
	fHistPiAPDecayLengthVsPt[isSecd]->Fill(ptALambda,ctAL);
	fHistPiAPDecayLengthVsMass[isSecd]->Fill(massALambda,ctAL);

	//-- secondaries --//
	if(isSecd==1){
	  if(fabs(pdgMother)==3112 || fabs(pdgMother)==3114 || fabs(pdgMother)==3222 || fabs(pdgMother)==3224 || fabs(pdgMother)==3214 ){
	    fHistPiAPMassVSPtSecSigma[1]->Fill(massALambda,ptALambda);
	  }
	  //  if(fabs(pdgMother)==3334){
	  // 		  fHistPiAPDCAtoPrimVtxOmega[1]->Fill(massALambda,dcaV0ToPrimVertex);
	  // 	       }
	  if(fabs(pdgMother) == 3322){
	    fHistPiAPCosPointAngXiVsPt->Fill(ptALambda,cosOPAng);
	    fHistPiAPMassVSPtSecXi[1]->Fill(massALambda,ptALambda);
	    fHistPiAPMassVSYSecXi[1]->Fill(massALambda,rapAL);
	    fHistPiAPXi0PtVSLambdaPt[1]->Fill(ptALambda,ptXiMother);
	  }
	  if(pdgMother == -3312){
	    fHistPiAPCosPointAngXiVsPt->Fill(ptALambda,cosOPAng);
	    fHistPiAPMassVSPtSecXi[1]->Fill(massALambda,ptALambda);
	    fHistPiAPMassVSYSecXi[1]->Fill(massALambda,rapAL);
	    fHistPiAPXiMinusPtVSLambdaPt[1]->Fill(ptALambda,ptXiMother);
	  }
	}
      }
      if(ptALambda>0.4) fHistArmenteros[isSecd]->Fill(alfa,qt);
      fHistDedxSecAProt[isSecd]->Fill(negDaughterP,tpcsigNeg);
      fHistDedxSecPiPlus[isSecd]->Fill(posDaughterP,tpcsigPos);
      fHistV0RadiusZ[isSecd]->Fill(dim2V0Radius,xr[2]);
      fHistV0RadiusXY[isSecd]->Fill(xr[0],xr[1]);
      fHistV0RadiusZVSPt[isSecd]->Fill(ptALambda,dim2V0Radius);
      fHistV0RadiusXYVSY[isSecd]->Fill(rapAL,dim2V0Radius);

    }
    
    //-- fill detector histos general --//
    if(lambdaOK || alambdaOK || k0sOK){
      //-- pos --//
      fHistNclsITSPos[isSecd]->Fill(fabs(posDaughterPt),nclsITSPos);
      fHistNclsTPCPos[isSecd]->Fill(nclsTPCFindablePos,nclsTPCPos);
      fHistNCRowsTPCPos[isSecd]->Fill(fabs(posDaughterPt),crossedRowsTPCPos);
      fHistChi2PerNclsITSPos[isSecd]->Fill(fabs(posDaughterPt),chi2PerClusterITSPos);
      //--neg --//
      fHistNclsITSNeg[isSecd]->Fill(fabs(negDaughterPt),nclsITSNeg);
      fHistNclsTPCNeg[isSecd]->Fill(nclsTPCFindableNeg,nclsTPCNeg);
      fHistNCRowsTPCNeg[isSecd]->Fill(fabs(negDaughterPt),crossedRowsTPCNeg);
      fHistChi2PerNclsITSNeg[isSecd]->Fill(fabs(negDaughterPt),chi2PerClusterITSNeg);
      
      fHistNclsITS[isSecd]->Fill(nclsITSPos,nclsITSNeg);
      //if(negDaughterPt >1.0)
      fHistNclsTPC[isSecd]->Fill(crossedRowsTPCNeg,nclsTPCNeg);
      
      fHistPiPiEtaDReco[1]->Fill(posDaughterEta,eta00);
      fHistPiPEtaDReco[1]->Fill(negDaughterEta,eta01);
    }
    
    /*
    //-- AliKFParticle --//
    if (negPiKF) delete negPiKF; negPiKF=NULL;
    if (posPiKF) delete posPiKF; posPiKF=NULL;
    if (posPKF) delete posPKF; posPKF=NULL;
    if (negAPKF) delete negAPKF; negAPKF=NULL;
    */
    if(fMCMode && fMCTruthMode) break;// otherwise we would not have ended up here
  }//end V0 reco loop

}
  
//________________________________________________________________________

Int_t AliAnalysisTaskV0ForRAA::CalculateCentralityBin(){
  //find centrality bin for centrality selection

  if (fUseCentrality == 0) return -1;

  AliCentrality *esdCentrality = fESD->GetCentrality();

  Float_t centralityVZERO  = esdCentrality->GetCentralityPercentile("V0M");  
  Float_t centralitySPD    = esdCentrality->GetCentralityPercentile("CL1");

  Int_t centralityVZEROBin = -1;
  Int_t centralitySPDBin   = -1;

  //-- SPD centrality --//
  if ( fUseCentrality == 2 ){
    if      ( centralitySPD >=  0. && centralitySPD <   5.) centralitySPDBin =  0;
    else if ( centralitySPD >=  5. && centralitySPD <  10.) centralitySPDBin =  5;
    else if ( centralitySPD >= 10. && centralitySPD <  20.) centralitySPDBin = 10;
    else if ( centralitySPD >= 20. && centralitySPD <  30.) centralitySPDBin = 20;
    else if ( centralitySPD >= 30. && centralitySPD <  40.) centralitySPDBin = 30;
    else if ( centralitySPD >= 40. && centralitySPD <  50.) centralitySPDBin = 40;
    else if ( centralitySPD >= 50. && centralitySPD <  60.) centralitySPDBin = 50;
    else if ( centralitySPD >= 60. && centralitySPD <  70.) centralitySPDBin = 60;
    else if ( centralitySPD >= 70. && centralitySPD <  80.) centralitySPDBin = 70;
    else if ( centralitySPD >= 80. && centralitySPD <  90.) centralitySPDBin = 80;
    else if ( centralitySPD >= 90. && centralitySPD <  99.) centralitySPDBin = 90;
    else if ( centralitySPD >= 99. ) centralitySPDBin = 100;
    else if ( fabs(centralitySPD)< 0.0001 ) centralitySPDBin = 100;
    return centralitySPDBin;
  }
  //-- V0 centrality --//
  if ( fUseCentrality == 1 ){
    if      ( centralityVZERO >  0. && centralityVZERO <   5.) centralityVZEROBin =  0;
    else if ( centralityVZERO >=  5. && centralityVZERO <  10.) centralityVZEROBin =  5;
    else if ( centralityVZERO >= 10. && centralityVZERO <  20.) centralityVZEROBin = 10;
    else if ( centralityVZERO >= 20. && centralityVZERO <  30.) centralityVZEROBin = 20;
    else if ( centralityVZERO >= 30. && centralityVZERO <  40.) centralityVZEROBin = 30;
    else if ( centralityVZERO >= 40. && centralityVZERO <  50.) centralityVZEROBin = 40;
    else if ( centralityVZERO >= 50. && centralityVZERO <  60.) centralityVZEROBin = 50;
    else if ( centralityVZERO >= 60. && centralityVZERO <  70.) centralityVZEROBin = 60;
    else if ( centralityVZERO >= 70. && centralityVZERO <  80.) centralityVZEROBin = 70;
    else if ( centralityVZERO >= 80. && centralityVZERO <  90.) centralityVZEROBin = 80;
    else if ( centralityVZERO >= 90. && centralityVZERO <  99.) centralityVZEROBin = 90;
    else if ( centralityVZERO >= 99. ) centralityVZEROBin = 100;
    else if ( fabs(centralityVZERO)< 0.0001 ) centralityVZEROBin = 100;
    return centralityVZEROBin;
  }
  return -1;
  
}
//________________________________________________________________________
Bool_t  AliAnalysisTaskV0ForRAA::GetMCTruthPartner(AliESDtrack *pos,AliESDtrack *neg,Int_t id0,Int_t id1){
  //-- get daughter label and check it --//
  Int_t labelP = fabs(pos->GetLabel());
  Int_t labelN = fabs(neg->GetLabel());

  if (labelN==labelP)  return kFALSE;
   
  if(fMCTruthMode){
    if ((labelP!=id0) && (labelP!=id1))  return kFALSE;
    if ((labelN!=id0) && (labelN!=id1))  return kFALSE;
  }

  return kTRUE;
}
