/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
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

/* $Id$ */

/////////////////////////////////////////////////////////////
//
//  FIX SETTING OF CUT OBJECT IN CONSTRUCTOR.. also for EVENT ELECTION AND TRIGGER
//  CHECK AND FIX USAGE OF NORMALISATION COUNTER
//  FIX TRACK SELECTION METHOD
//  FIX PID SEL 
//  FIX ORDERING SCHEME (WRONG NOW)
// DELETION OF SECONDARY VTX
//  CLEAN ADD TASK AND SET CONTAINER NAMES
// SELCTION WITH CUTS TO BE FIXED FOR CASE 0: SHOULD NOT DELETE, WE MUST ADD A BIT TO THE HISTOGRAM OR DUPLICATE THE HISTO FOR INCLUDING XIC
// FIDUCIAL ACC: FOR XIC the Lc is used
// RECALCULATION OF PRIM VTX SHOULD NOT BE DEFINED BY FLAG IN TASK AS IT IS NOW
// CUT OBJECTS USED IN ANALYSIS SHOULD ALSO BE STREAMED 
// CHECK FILLING OF RECO STEP FOR EFFICIENCY CALCULATION, THERE MIGHT BE A BIAS RELATED TO THE MASS SELECTION, which is not checked to fill that histo
// DIST12 and DIST23: set hard-coded to 500 micron
// ADD MASS CUT AROUND SIGMA_C MASS: IS IT CORRECT?
// CUTS TO SIGMAVTX AND Sumd0^2 are hard coded in Init function, should be changed!!!
////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TClonesArray.h>
#include <TArrayI.h>
#include <TLorentzVector.h>
#include <TObjArray.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TDatabasePDG.h>
#include <THnSparse.h>
#include "AliVertexingHFUtils.h"
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDtrack.h"
#include "AliVertexerTracks.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoCascadeHF.h"

#include "AliRDHFCutsLctopKpi.h"
#include "AliRDHFCutsXictopKpi.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSEXicTopKpi.h"
#include "AliNormalizationCounter.h"
#include "AliPIDResponse.h"

using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSEXicTopKpi);
/// \endcond

//________________________________________________________________________
AliAnalysisTaskSEXicTopKpi::AliAnalysisTaskSEXicTopKpi():
  AliAnalysisTaskSE(),
  fCuts(0x0),
  //fCutsLc(0x0),
  fCutsXic(0x0),
  fCounter(0),
  fPidResponse(0),
  fReadMC(kFALSE),
  fAnalysisType(0),
  fRecalPrimVtx(kFALSE),
  fNentries(0),
  fhistMonitoring(0x0),
  ftrackArraySel(0),
  ftrackArraySelSoftPi(0x0),
  fnSel(),
  fnSelSoftPi(),
  ftrackSelStatusProton(),
  ftrackSelStatusKaon(),
  ftrackSelStatusPion(),
  fSys(0),
  fAODProtection(1),
  fLikeSign(0),
  fESDtrackCutsProton(0x0),
  fESDtrackCutsKaon(0x0),
  fESDtrackCutsPion(0x0),
  fESDtrackCutsSoftPion(0x0),
  fprimVtx(0x0),
  fhistInvMassCheck(0x0),
  fhistMCSpectrumAccLc(0x0),
  fhistMCSpectrumAccXic(0x0),
  fhSparseAnalysis(0x0),
  fhSparseAnalysisSigma(0x0),
  fCosPointDistrAll(0x0),
  fCosPointDistrAllFilter(0x0),
  fCosPointDistrSignal(0x0),
  fCosPointDistrSignalFilter(0x0),
  fDist12Signal(0x0),
  fDist12SignalFilter(0x0),
  fDist12All(0x0),
  fDist12AllFilter(0x0),
  fDist23Signal(0x0),
  fDist23All(0x0),
  fnSigmaPIDtofProton(0x0),
  fnSigmaPIDtofPion(0x0),
  fnSigmaPIDtofKaon(0x0),
  fnSigmaPIDtpcProton(0x0),
  fnSigmaPIDtpcPion(0x0),
  fnSigmaPIDtpcKaon(0x0),
  fOutput(0x0),
  fVertexerTracks(0x0),
  fSetTrackCutLcFilteringPP(kFALSE),
  fCutSelLevel(0),
  fApplykFirst(kFALSE),
  fMaxPtTrackkFirst(0.),
  fMaxVtxChi2Cut(0.),
  //fFillTree(0),
  fFillTree(kFALSE),
  fTreeVar(0x0)
  ,fpT_down(-1)
  ,fLowpT_down(-1)
  ,fHighpT_down(-1)
  ,fCompute_dist12_dist23(kFALSE)
  ,fExplore_PIDstdCuts(kFALSE)
{
  /// Default constructor

}

//________________________________________________________________________
AliAnalysisTaskSEXicTopKpi::AliAnalysisTaskSEXicTopKpi(const char *name,AliRDHFCutsD0toKpi* cuts):
  AliAnalysisTaskSE(name),
  fCuts(0x0),
  //fCutsLc(0x0),
  fCutsXic(0x0),
  fCounter(0),
  fPidResponse(0),
  fReadMC(kFALSE),
  fAnalysisType(0),
  fRecalPrimVtx(kFALSE),
  fNentries(0),
  fhistMonitoring(0x0),
  ftrackArraySel(0),
  ftrackArraySelSoftPi(0x0),
  fnSel(),
  fnSelSoftPi(),
  ftrackSelStatusProton(),
  ftrackSelStatusKaon(),
  ftrackSelStatusPion(),
  fSys(0),
  fAODProtection(1),
  fLikeSign(0),
  fESDtrackCutsProton(0x0),
  fESDtrackCutsKaon(0x0),
  fESDtrackCutsPion(0x0),
  fESDtrackCutsSoftPion(0x0),
  fprimVtx(0x0),
  fhistInvMassCheck(0x0),
  fhistMCSpectrumAccLc(0x0),
  fhistMCSpectrumAccXic(0x0),
  fhSparseAnalysis(0x0),
  fhSparseAnalysisSigma(0x0),
  fCosPointDistrAll(0x0),
  fCosPointDistrAllFilter(0x0),
  fCosPointDistrSignal(0x0),
  fCosPointDistrSignalFilter(0x0),
  fDist12Signal(0x0),
  fDist12SignalFilter(0x0),
  fDist12All(0x0),
  fDist12AllFilter(0x0),
  fDist23Signal(0x0),
  fDist23All(0x0),
  fnSigmaPIDtofProton(0x0),
  fnSigmaPIDtofPion(0x0),
  fnSigmaPIDtofKaon(0x0),
  fnSigmaPIDtpcProton(0x0),
  fnSigmaPIDtpcPion(0x0),
  fnSigmaPIDtpcKaon(0x0),
  fOutput(0x0),
  fVertexerTracks(0x0),
  fSetTrackCutLcFilteringPP(kFALSE),
  fCutSelLevel(0),
  fApplykFirst(kFALSE),
  fMaxPtTrackkFirst(0.),
  fMaxVtxChi2Cut(0.),
  //fFillTree(0),
  fFillTree(kFALSE),
  fTreeVar(0x0) 
  ,fpT_down(-1)
  ,fLowpT_down(-1)
  ,fHighpT_down(-1)
  ,fCompute_dist12_dist23(kFALSE)
  ,fExplore_PIDstdCuts(kFALSE)
{
  /// Default constructor


  DefineOutput(1,TH1F::Class());  //My private output
  DefineOutput(2,AliNormalizationCounter::Class());
  DefineOutput(3,TList::Class());
  fCuts=cuts;
}

//________________________________________________________________________
AliAnalysisTaskSEXicTopKpi::~AliAnalysisTaskSEXicTopKpi()
{
  if (fCuts) {
    delete fCuts;
    fCuts = 0;
  }
  //if(fCutsLc){
  //  delete fCutsLc;
  //  fCutsLc =0;
  //}
  if(fCutsXic){
    delete fCutsXic;
    fCutsXic =0;
  }
  if (fNentries){
    delete fNentries;
    fNentries = 0;
  }
  if(  fhistMonitoring){
    delete   fhistMonitoring;
    fhistMonitoring=0;
  }
  if(fCounter){
    delete fCounter;
    fCounter=0;
  }
  if(ftrackArraySel)delete ftrackArraySel;
  if(ftrackArraySelSoftPi)delete ftrackArraySelSoftPi;

  if(ftrackSelStatusProton)delete ftrackSelStatusProton;
  if(ftrackSelStatusKaon)delete ftrackSelStatusKaon;
  if(ftrackSelStatusPion)delete ftrackSelStatusPion;

  if(fESDtrackCutsProton)delete fESDtrackCutsProton;
  if(fESDtrackCutsKaon) delete fESDtrackCutsKaon;
  if(fESDtrackCutsPion) delete fESDtrackCutsPion;
  if(fESDtrackCutsSoftPion) delete fESDtrackCutsSoftPion;
  if(  fhistInvMassCheck) delete   fhistInvMassCheck;
  if(fhistMCSpectrumAccLc)delete fhistMCSpectrumAccLc;
  if(fhistMCSpectrumAccXic)delete fhistMCSpectrumAccXic;
  if(fhSparseAnalysis)delete fhSparseAnalysis;
  if(fhSparseAnalysisSigma)delete fhSparseAnalysisSigma;
  if(fCosPointDistrAll)delete fCosPointDistrAll;
  if(fCosPointDistrSignal)delete fCosPointDistrSignal;
  if(fCosPointDistrAllFilter)delete fCosPointDistrAllFilter;
  if(fCosPointDistrSignalFilter)delete fCosPointDistrSignalFilter;
  if(fDist12Signal)delete fDist12Signal;
  if(fDist12SignalFilter)delete fDist12SignalFilter;
  if(fDist12All)delete fDist12All;
  if(fDist12AllFilter)delete fDist12AllFilter;
  if(fDist23Signal)delete fDist23Signal;
  if(fDist23All)delete fDist23All;
  if(fOutput)delete fOutput;
  if(fVertexerTracks)delete fVertexerTracks;
  if(fnSigmaPIDtofProton)delete fnSigmaPIDtofProton;
  if(fnSigmaPIDtofPion)delete fnSigmaPIDtofPion;
  if(fnSigmaPIDtofKaon)delete fnSigmaPIDtofKaon;
  if(fnSigmaPIDtpcProton)delete fnSigmaPIDtpcProton;
  if(fnSigmaPIDtpcPion)delete fnSigmaPIDtpcPion;
  if(fnSigmaPIDtpcKaon)delete fnSigmaPIDtpcKaon;
  if(fTreeVar)delete fTreeVar;
}  

//________________________________________________________________________
void AliAnalysisTaskSEXicTopKpi::Init()
{
  /// Initialization

  if(fDebug > 1) printf("AnalysisTaskSEXicTopKpi::Init() \n");

  
//   AliRDHFCutsD0toKpi* copyfCuts=new AliRDHFCutsD0toKpi(*fCuts);
//   const char* nameoutput=GetOutputSlot(4)->GetContainer()->GetName();
//   copyfCuts->SetName(nameoutput);

  if(!fCutsXic){
    fCutsXic = new AliRDHFCutsXictopKpi("CutsXictopKpi");
    //Float_t cutsArrayLctopKpi[12]={0.2,0.4,0.4,0.,0.,0.01,0.06,0.02,0.,0.85,0.,10000000000.};
    //Float_t cutsArrayLctopKpi[13]={0.18,0.4,0.5,0.,0.,0.01,0.06,0.005,0.7,0.0,0.,0.05,0.4};
    Float_t cutsArrayLctopKpi[13]={0.40,0.4,0.5,0.,0.,0.00,0.06,0.000,0.0,-1.,0.,0.05,0.4}; // default filtering cuts in pp collisions with enlarged mass cut (from 0.18 to 0.400) to store also candidates in Xic mass window
    
    /*THE VARIABLES ARE: 
      "inv. mass [GeV]",
      "pTK [GeV/c]",
      "pTP [GeV/c]",
      "d0K [cm]   lower limit!",
      "d0Pi [cm]  lower limit!",
      "dist12 (cm)",
      "sigmavert (cm)",
      "dist prim-sec (cm)",
      "pM=Max{pT1,pT2,pT3} (GeV/c)",
      "cosThetaPoint",
      "Sum d0^2 (cm^2)",
      "dca cut (cm)",
      "cut on pTpion [GeV/c]" */
    fCutsXic->SetCuts(13,cutsArrayLctopKpi);
  }
  else{
    /*
    Double_t sigmaVtxMax[8]={0.09,0.09,0.05,0.035,0.035,0.03,0.03,0.025};
    Double_t sumd02[8]={0.0003,0.0003,0.0002,0.00015,0.00015,0.0001,0.,0.};
    if(fAnalysisType==0 || fAnalysisType==4){// assure mass range is large enough
      Int_t nvars=fCutsXic->GetNVars();
      Int_t nptbins=fCutsXic->GetNPtBins();
      Float_t **cutvalues=new Float_t*[nvars];
      for(Int_t iv=0;iv<nvars;iv++){
	cutvalues[iv]=new Float_t[nptbins];
      }
      fCutsXic->GetCuts(cutvalues);
      for(Int_t km=0;km<nptbins;km++){
	//cutvalues[0][km]=0.40;
	cutvalues[6][km]=sigmaVtxMax[km]; // TO BE CHANGED IN FUTURE
	cutvalues[10][km]=sumd02[km];
      }
      fCutsXic->SetCuts(nvars,nptbins,cutvalues);
      Printf("Xic Cuts modified to assure mass window is large enough, current cuts are:");*/
      Printf("\n--- Adopted cuts ---");
      fCutsXic->PrintAll();
    //}
  }
  
  if(fDebug>=0 || fSetTrackCutLcFilteringPP){// track cuts used for Lc filtering (in pp, 2018): need to set them to be sure that only tighter cuts than these are used
    AliESDtrackCuts *esdTrCuts=new AliESDtrackCuts("AliESDtrackCuts","default");
    esdTrCuts->SetRequireTPCRefit(kTRUE);
    esdTrCuts->SetMinNClustersTPC(70);
    esdTrCuts->SetRequireITSRefit(kTRUE);
    //esdTrCuts->SetMinNClustersITS(4);
    esdTrCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					   AliESDtrackCuts::kAny);
    esdTrCuts->SetMinDCAToVertexXY(0.);
    esdTrCuts->SetPtRange(0.4,1.e10);// ORIGINAL VALUE IS 0.3 BUT WE CUT AT 0.4 in ANALYSIS
    esdTrCuts->SetEtaRange(-0.8,+0.8);
    fCutsXic->AddTrackCuts(esdTrCuts);
  }

  if((fAnalysisType==0 || fAnalysisType==3) && !fESDtrackCutsSoftPion){
    fESDtrackCutsSoftPion = new AliESDtrackCuts("AliESDtrackCuts","default");
    fESDtrackCutsSoftPion->SetRequireITSRefit(kTRUE);
    fESDtrackCutsSoftPion->SetMinNClustersITS(3);
    fESDtrackCutsSoftPion->SetMaxDCAToVertexXY(0.065);
    //    fESDtrackCutsSoftPion->SetPtRange(0.1,1.e10);
    fESDtrackCutsSoftPion->SetMaxDCAToVertexZ(0.15);
    fESDtrackCutsSoftPion->SetEtaRange(-0.9,+0.9);    
  } 


  // protection against negative values for downsampling for fTreeVar filling
  if(fFillTree){
    if(fpT_down<0 || fLowpT_down<0 || fHighpT_down<0){
      if(fpT_down<0)      fpT_down     = 4;
      if(fLowpT_down<0)   fLowpT_down  = 0.005;
      if(fHighpT_down<0)  fHighpT_down = 0.05;
      printf("\n\n=== fTreeVar filling activated without downsampling parameters ===\n");
      printf("    Using downsampling rescue values.\n    Adopted values:\n");
      printf("      fpT_down = %f\n      fLowpT_down = %f\n      fHighpT_down = %f\n\n\n",fpT_down,fLowpT_down,fHighpT_down);
    }
  }

  // print wheter we require the calulcation of dist12 and dist23
  if(fCompute_dist12_dist23)  printf("\n\n===== fCompute_dist12_dist23 is kTRUE ---> dist12 and dist23 calculated =====\n\n");

  // Post the data
  //  PostData(4,copyfCuts);


  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEXicTopKpi::UserCreateOutputObjects()
{

  /// Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSEXicTopKpi::UserCreateOutputObjects() \n");
  if(fRecalPrimVtx && fSys>0)AliFatal("Cannot recalculcate prim vtx in p-Pb and Pb-Pb");
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("listOutput");


  fhistMonitoring=new TH1F("fhistMonitoring","Overview",30,-0.5,29.5);
  fhistMonitoring->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fhistMonitoring->GetXaxis()->SetBinLabel(2,"nEventsSel");
  fhistMonitoring->GetXaxis()->SetBinLabel(3,"nEvGoodVtx");
  fhistMonitoring->GetXaxis()->SetBinLabel(4,"nTracksTot");
  fhistMonitoring->GetXaxis()->SetBinLabel(5,"nTracksSel");
  fhistMonitoring->GetXaxis()->SetBinLabel(6,"nPionComp");
  fhistMonitoring->GetXaxis()->SetBinLabel(7,"nKaonComp");
  fhistMonitoring->GetXaxis()->SetBinLabel(8,"nProtonComp");
  fhistMonitoring->GetXaxis()->SetBinLabel(9,"nTriplets");
  fhistMonitoring->GetXaxis()->SetBinLabel(10,"nCandidates");
  fhistMonitoring->GetXaxis()->SetBinLabel(11,"nCandSel");

  fOutput->Add(fhistMonitoring);

  TString strnameoutput="EntryCounter";
  fNentries=new TH1F(strnameoutput.Data(), "Integral(1,2) = number of AODs *** Integral(2,3) = number of candidates selected with cuts *** Integral(3,4) = number of D0 selected with cuts *** Integral(4,5) = events with good vertex ***  Integral(5,6) = pt out of bounds", 23,-0.5,22.5);

  fNentries->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fNentries->GetXaxis()->SetBinLabel(2,"nCandSel(Cuts)");
  if(fReadMC) fNentries->GetXaxis()->SetBinLabel(3,"nD0Selected");
  else fNentries->GetXaxis()->SetBinLabel(3,"Dstar<-D0");
  fNentries->GetXaxis()->SetBinLabel(4,"nEventsGoodVtxS");
  fNentries->GetXaxis()->SetBinLabel(5,"ptbin = -1");
  fNentries->GetXaxis()->SetBinLabel(6,"no daughter");
  if(fSys==0) fNentries->GetXaxis()->SetBinLabel(7,"nCandSel(Tr)");
  //  if(fFillVarHists || fPIDCheck){
    fNentries->GetXaxis()->SetBinLabel(8,"PID=0");
    fNentries->GetXaxis()->SetBinLabel(9,"PID=1");
    fNentries->GetXaxis()->SetBinLabel(10,"PID=2");
    fNentries->GetXaxis()->SetBinLabel(11,"PID=3");
    //  }
  if(fReadMC && fSys==0){
    fNentries->GetXaxis()->SetBinLabel(12,"K");
    fNentries->GetXaxis()->SetBinLabel(13,"Lambda");
  }
  fNentries->GetXaxis()->SetBinLabel(14,"Pile-up Rej");
  fNentries->GetXaxis()->SetBinLabel(15,"N. of 0SMH");
  if(fSys==1) fNentries->GetXaxis()->SetBinLabel(16,"Nev in centr");
  //  if(fIsRejectSDDClusters)
  fNentries->GetXaxis()->SetBinLabel(17,"SDD-Cls Rej");
  fNentries->GetXaxis()->SetBinLabel(18,"Phys.Sel.Rej");
  fNentries->GetXaxis()->SetBinLabel(19,"D0 failed to be filled");
  fNentries->GetXaxis()->SetBinLabel(20,"fisFilled is 0");
  fNentries->GetXaxis()->SetBinLabel(21,"fisFilled is 1");
  fNentries->GetXaxis()->SetBinLabel(22,"AOD/dAOD mismatch");
  fNentries->GetXaxis()->SetBinLabel(23,"AOD/dAOD #events ok");
  fNentries->GetXaxis()->SetNdivisions(1,kFALSE);

  fCounter = new AliNormalizationCounter(Form("%s",GetOutputSlot(2)->GetContainer()->GetName()));
  fCounter->Init();


  ftrackArraySel=new TArrayI(10000);
  ftrackArraySelSoftPi=new TArrayI(10000);
  ftrackSelStatusProton=new TArrayI(10000);
  ftrackSelStatusKaon=new TArrayI(10000);
  ftrackSelStatusPion=new TArrayI(10000);

  
  //fhistInvMassCheck=new TH2F("fhistInvMassCheck","InvDistrCheck",1000,1.600,2.800,5,-0.5,4.5);
  fhistInvMassCheck=new TH2F("fhistInvMassCheck","InvDistrCheck",1000,1.600,2.800,5,-0.5,15.5);
  //  fhistCheckPIDTOFTPC=new TH3F("fhistCheckPIDTOFTPC","fhistCheckPIDTOFTPC",
  fhistMCSpectrumAccLc=new TH2F("fhistMCSpectrumAccLc","fhistMCSpectrumAccLc",250,0,50,15,-0.5,14.5); // 
  fhistMCSpectrumAccXic=new TH2F("fhistMCSpectrumAccXic","fhistMCSpectrumAccXic",250,0,50,15,-0.5,14.5); // 
  
  //Int_t nbinsSparse[7]={16,125,10,16,20,10,10};
  //Double_t lowEdges[7]={0,2.15,0.,0,0.8,0,-1};
  //Double_t upEdges[7]={16,2.65,0.0500,8,1.,5,9};
  //if(!fFillTree)  fhSparseAnalysis=new THnSparseF("fhSparseAnalysis","fhSparseAnalysis;pt;mass;Lxy;nLxy;cosThatPoint;normImpParXY;seleFlag",7,nbinsSparse,lowEdges,upEdges);
  //Int_t nbinsSparseSigma[9]={16,200,10,12,10,10,10,22,20};
  //Double_t lowEdgesSigma[9]={0,0.130,0.,0,0.8,0,-1,2.262,-1};
  //Double_t upEdgesSigma[9]={16,0.330,0.0500,6.,1.,5,9,2.306,1};
  //if(!fFillTree)  fhSparseAnalysisSigma=new THnSparseF("fhSparseAnalysisSigma","fhSparseAnalysis;pt;deltamass;Lxy;nLxy;cosThetaPoint;normImpParXY;seleFlag;LcMass;CosThetaStarSoftPion",9,nbinsSparseSigma,lowEdgesSigma,upEdgesSigma);
  
  // adding PID cases study
  Int_t nbinsSparse[8]={16,125,10,16,20,10,10,11};
  Double_t lowEdges[8]={0,2.15,0.,0,0.8,0,-1,-0.5};
  Double_t upEdges[8]={16,2.65,0.0500,8,1.,5,9,10.5};
  if(!fFillTree)  fhSparseAnalysis=new THnSparseF("fhSparseAnalysis","fhSparseAnalysis;pt;mass;Lxy;nLxy;cosThatPoint;normImpParXY;seleFlag;PIDcase",8,nbinsSparse,lowEdges,upEdges);
  Int_t nbinsSparseSigma[10]={16,200,10,12,10,10,10,1,22,20};
  Double_t lowEdgesSigma[10]={0,0.130,0.,0,0.8,0,-1,-0.5,2.262,-1};
  Double_t upEdgesSigma[10]={16,0.330,0.0500,6.,1.,5,9,0.5,2.306,1};
  if(!fFillTree)  fhSparseAnalysisSigma=new THnSparseF("fhSparseAnalysisSigma","fhSparseAnalysis;pt;deltamass;Lxy;nLxy;cosThetaPoint;normImpParXY;seleFlag;PIDcase;LcMass;CosThetaStarSoftPion",10,nbinsSparseSigma,lowEdgesSigma,upEdgesSigma);
  



  fCosPointDistrAll=new TH1F("fCosPointDistrAll","fCosPointDistrAll",200,-1.1,1.1);
  fCosPointDistrSignal=new TH1F("fCosPointDistrSignal","fCosPointDistrSignal",200,-1.1,1.1);
  fCosPointDistrAllFilter=new TH1F("fCosPointDistrAllFilter","fCosPointDistrAllFilter",200,-1.1,1.1);
  fCosPointDistrSignalFilter=new TH1F("fCosPointDistrSignalFilter","fCosPointDistrSignalFilter",200,-1.1,1.1);
  fDist12Signal=new TH1F("fDist12Signal","fDist12Signal",500,0.,1000.);
  fDist12SignalFilter=new TH1F("fDist12SignalFilter","fDist12SignalFilter",500,0.,1000.);
  fDist23Signal=new TH1F("fDist23Signal","fDist23Signal",500,0.,1000.);
  fDist12All=new TH1F("fDist12All","fDist12All",500,0.,1000.);
  fDist12AllFilter=new TH1F("fDist12AllFilter","fDist12AllFilter",500,0.,1000.);
  fDist23All=new TH1F("fDist23All","fDist23All",500,0.,1000.);


  fnSigmaPIDtofProton=new TH2F("fnSigmaPIDtofProton","fnSigmaPIDtofProton",100,0,20,80,-10,10);
  fnSigmaPIDtofPion=new TH2F("fnSigmaPIDtofPion","fnSigmaPIDtofPion",100,0,20,80,-10,10);
  fnSigmaPIDtofKaon=new TH2F("fnSigmaPIDtofKaon","fnSigmaPIDtofKaon",100,0,20,80,-10,10);
  fnSigmaPIDtpcProton=new TH2F("fnSigmaPIDtpcProton","fnSigmaPIDtpcProton",100,0,20,80,-10,10);
  fnSigmaPIDtpcPion=new TH2F("fnSigmaPIDtpcPion","fnSigmaPIDtpcPion",100,0,20,80,-10,10);
  fnSigmaPIDtpcKaon=new TH2F("fnSigmaPIDtpcKaon","fnSigmaPIDtpcKaon",100,0,20,80,-10,10);

/*
  //  pt vs. pointing angle, lxy, nlxy, ptP,ptK,ptPi,vtxchi2,sigmaVtx,sumd02,dca1,dca2,dca3,nd01,nd02,nd03,Lc d0
  //  Float_t pt,pAngle,lxy,nlxy,ptP,ptK,ptPi,vtxchi2,sigmaVtx,sumd02,dca1,dca2,dca3,nd01,nd02,nd03,d0Lc;
  Float_t var[20];
  TString varNames[20]={"pt","pAngle","lxy","nlxy","ptP","ptK","ptPi","vtxchi2","sigmaVtx","sumd02","dca1","dca2","dca3","nd01","nd02","nd03","d0Lc","cosThetaStar1","cosThetaStar2","flagMC"};
  fTreeVar=new TTree("T","tree with variables");
  for(Int_t k=0;k<20;k++){
    fTreeVar->Branch(varNames[k].Data(),&var[k]);
  }
  */
  // 
  //  Rescaling of reconstructed to Lc as they were Xic added in the end of the tree
  //
  Float_t var[33];
  Short_t resp;
  TString varNames[34]={"pt","pAngle","lxy","nlxy","ptP","ptK","ptPi","vtxchi2","sigmaVtx","sumd02","dca1","dca2","dca3","nd01","nd02","nd03","d0Lc","cosThetaStar1","cosThetaStar2","m_pKpi","m_piKp","flagMC","w_FromLc_toXic","nSig_TPC_prot_0","nSig_TOF_prot_0","nSig_TPC_pion_0","nSig_TOF_pion_0","nSig_TPC_kaon_1","nSig_TOF_kaon_1","nSig_TPC_prot_2","nSig_TOF_prot_2","nSig_TPC_pion_2","nSig_TOF_pion_2","massHypoFilt_respCuts_respPID"};
  fTreeVar=new TTree("T","tree with variables");
  for(Int_t k=0;k<33;k++){
    fTreeVar->Branch(varNames[k].Data(),&var[k]);
  }
  fTreeVar->Branch(varNames[33].Data(),&resp);
  fOutput->Add(fTreeVar);


  fOutput->Add(fDist12Signal);
  fOutput->Add(fDist12SignalFilter);
  fOutput->Add(fDist12All);
  fOutput->Add(fDist12AllFilter);
  fOutput->Add(fDist23Signal);
  fOutput->Add(fDist23All);
  fOutput->Add(fCosPointDistrAll);
  fOutput->Add(fCosPointDistrAllFilter);
  fOutput->Add(fCosPointDistrSignal);
  fOutput->Add(fCosPointDistrSignalFilter);
  fOutput->Add(fhistInvMassCheck);
  fOutput->Add(fhistMCSpectrumAccLc);
  fOutput->Add(fhistMCSpectrumAccXic);
  if(fhSparseAnalysis)  fOutput->Add(fhSparseAnalysis);
  if(fhSparseAnalysisSigma)  fOutput->Add(fhSparseAnalysisSigma);
  fOutput->Add(fnSigmaPIDtofProton);
  fOutput->Add(fnSigmaPIDtofPion);
  fOutput->Add(fnSigmaPIDtofKaon);
  fOutput->Add(fnSigmaPIDtpcProton);
  fOutput->Add(fnSigmaPIDtpcPion);
  fOutput->Add(fnSigmaPIDtpcKaon);


  // Post the data
  //   PostData(1,fOutputMass);
  //   PostData(2,fDistr);

  PostData(1,fNentries);
  PostData(2,fCounter);  
  PostData(3,fOutput);
//   PostData(6,fOutputMassPt);
//   PostData(7,fVariablesTree);
//   PostData(8, fDetSignal);
//   PostData(9,fOutputMassY);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEXicTopKpi::UserExec(Option_t */*option*/)

{
  if(fDebug>=0)Printf("AliAnalysisTaskSEXicTopKpi: User Exec");
  
  /// Execute analysis for current event  
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  if(!aod) {
    printf("AliAnalysisTaskSEXicTopKpi::UserExec: input event not found!\n");
    return;
  }
  fhistMonitoring->Fill(0);
  
  if(fDebug>=0 && fAODProtection>=0){
    //   Protection against different number of events in the AOD and deltaAOD
    //   In case of discrepancy the event is rejected.
    Int_t matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
    if (matchingAODdeltaAODlevel<0 || (matchingAODdeltaAODlevel==0 && fAODProtection==1)) {
      // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
      fNentries->Fill(21);
      return;
    }
    fNentries->Fill(22);
  }
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
  fPidResponse=inputHandler->GetPIDResponse();
  if(!fPidResponse)return;

  if(fCutsXic->GetIsUsePID()){
    fCutsXic->GetPidHF()->SetPidResponse(fPidResponse);
    fCutsXic->GetPidpion()->SetPidResponse(fPidResponse);
    fCutsXic->GetPidprot()->SetPidResponse(fPidResponse);
    //fUtilPid=new AliAODpidUtil(fPidResponse);
    fCutsXic->GetPidHF()->SetOldPid(kFALSE);
    fCutsXic->GetPidpion()->SetOldPid(kFALSE);
    fCutsXic->GetPidprot()->SetOldPid(kFALSE);
  }
  

  //  if(!aod && AODEvent() && IsStandardAOD()) {
  //     // In case there is an AOD handler writing a standard AOD, use the AOD 
  //     // event in memory rather than the input (ESD) event.    
  //     aod = dynamic_cast<AliAODEvent*> (AODEvent());
  //     // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
  //     // have to taken from the AOD event hold by the AliAODExtension
  //     AliAODHandler* aodHandler = (AliAODHandler*) 
  //       ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());  
  //     if(aodHandler->GetExtensions()) {
  //       AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
  //       AliAODEvent* aodFromExt = ext->GetAOD();
  //       inputArray=(TClonesArray*)aodFromExt->GetList()->FindObject(bname.Data());
  //     }
  //   } else if(aod) {
  //     inputArray=(TClonesArray*)aod->GetList()->FindObject(bname.Data());
  //   }
  
  
  // fix for temporary bug in ESDfilter
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aod->GetPrimaryVertex() || TMath::Abs(aod->GetMagneticField())<0.001){ 
    Printf("vtx failure %p or wron B %f",aod->GetPrimaryVertex(),aod->GetMagneticField());
    return;
  }
  TClonesArray *mcArray = 0;
  AliAODMCHeader *mcHeader = 0;
  TClonesArray *lcArray = 0;
  
  if(fDebug>=0){
    lcArray=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
    if(!lcArray){
      Printf("Array of filtered Lc not present, delta-file not associated?");
      return;
    }
  }

  if(fReadMC) {
    // load MC particles
    mcArray = (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!mcArray) {
      printf("AliAnalysisTaskSEXicTopKpi::UserExec: MC particles branch not found!\n");
      return;
    }
    

    // check whether lc or xic are present
    for(Int_t kmc=0;kmc<mcArray->GetEntries();kmc++){
      AliAODMCParticle *mcpart=(AliAODMCParticle*)mcArray->At(kmc);
      Int_t pdg=mcpart->GetPdgCode();
      Int_t arrayDauLab[3];
      if(TMath::Abs(pdg)==4122){
	if(AliVertexingHFUtils::CheckLcpKpiDecay(mcArray, mcpart, arrayDauLab)>=1){
	  
	  if(TMath::Abs(mcpart->Y())<0.5){
	    fhistMCSpectrumAccLc->Fill(mcpart->Pt(),kGenLimAcc);// Gen Level
	  }
	  
	  Bool_t isInAcc=kTRUE;
	  // check GenAcc level
	  if(fCutsXic){
	    if(!fCutsXic->IsInFiducialAcceptance(mcpart->Pt(),mcpart->Y())){
	      isInAcc=kFALSE;
	    }
	  }
	  else {
	    if(TMath::Abs(mcpart->Y())>0.8){
	      isInAcc=kFALSE;
	    }
	  }
	  if(isInAcc){
	    fhistMCSpectrumAccLc->Fill(mcpart->Pt(),kGenAccMother);// Gen Acc Mother
	  }
	  for(Int_t k=0;k<3;k++){
	    AliAODMCParticle *mcpartdau=(AliAODMCParticle*)mcArray->At(arrayDauLab[k]);
	    if(TMath::Abs(mcpartdau->Eta())>0.9){
	      isInAcc=kFALSE;
	    }	    
	  }
	  if(isInAcc){
	    fhistMCSpectrumAccLc->Fill(mcpart->Pt(),kGenAcc);// Gen Acc
	  }
	}
      }
      else if(TMath::Abs(pdg)==4232){
	if(CheckXicpKpiDecay(mcArray, mcpart, arrayDauLab)>=1){
	  if(TMath::Abs(mcpart->Y())<0.5){
	    fhistMCSpectrumAccLc->Fill(mcpart->Pt(),kGenLimAcc);// Gen Level
	  }	  
	  Bool_t isInAcc=kTRUE;
	  // check GenAcc level
	  if(fCutsXic){
	    if(!fCutsXic->IsInFiducialAcceptance(mcpart->Pt(),mcpart->Y())){
	      isInAcc=kFALSE;
	    }	     
	  }
	  else {
	    if(TMath::Abs(mcpart->Y())>0.8){
	      isInAcc=kFALSE;
	    }
	  }
	  if(isInAcc){
	    fhistMCSpectrumAccLc->Fill(mcpart->Pt(),kGenAccMother);// Gen Acc Mother
	  }
	  for(Int_t k=0;k<3;k++){
	    AliAODMCParticle *mcpartdau=(AliAODMCParticle*)mcArray->At(arrayDauLab[k]);
	    if(TMath::Abs(mcpartdau->Eta())>0.9){
	      isInAcc=kFALSE;
	    }	    
	  }
	  if(isInAcc){
	    fhistMCSpectrumAccLc->Fill(mcpart->Pt(),kGenAcc);// Gen Acc
	  }
	}
      }	
    }
    
    
    // load MC header
    mcHeader = (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskSEXicTopKpi::UserExec: MC header branch not found!\n");
      return;
    }
  }
  //printf("VERTEX Z %f %f\n",vtx1->GetZ(),mcHeader->GetVtxZ());

  //------- IONUT CUT  ----------
//   Int_t nTPCout=0;
//   Float_t mTotV0=0;
//   AliAODVZERO* v0data=(AliAODVZERO*)((AliAODEvent*)aod)->GetVZEROData();
//   Float_t mTotV0A=v0data->GetMTotV0A();
//   Float_t mTotV0C=v0data->GetMTotV0C();
//   mTotV0=mTotV0A+mTotV0C;
//   Int_t ntracksEv = aod->GetNumberOfTracks();
//   for(Int_t itrack=0; itrack<ntracksEv; itrack++) { // loop on tacks
//     //    ... get the track
//     AliAODTrack * track = dynamic_cast<AliAODTrack*>(aod->GetTrack(itrack));
//     if(!track) {AliFatal("Not a standard AOD");}
//     if(track->GetID()<0)continue;
//     if((track->GetFlags())&(AliESDtrack::kTPCout)) nTPCout++;
//     else continue;
//   }
//   if(fhMultVZEROTPCoutTrackCorrNoCut) fhMultVZEROTPCoutTrackCorrNoCut->Fill(nTPCout,mTotV0);
//   Float_t mV0Cut=-2200.+(2.5*nTPCout)+(0.000012*nTPCout*nTPCout);
//   if(fEnablePileupRejVZEROTPCout){
//     if(mTotV0<mV0Cut) return;	
//   }


  
  //histogram filled with 1 for every AOD
  fNentries->Fill(0);
  fCounter->StoreEvent(aod,fCuts,fReadMC); 
  //fCounter->StoreEvent(aod,fReadMC); 
  // trigger class for PbPb C0SMH-B-NOPF-ALLNOTRD, C0SMH-B-NOPF-ALL
  TString trigclass=aod->GetFiredTriggerClasses();
  if(trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD") || trigclass.Contains("C0SMH-B-NOPF-ALL")) fNentries->Fill(14);
  //   if(fReadMC && fStepMCAcc){
  //     FillMCAcceptanceHistos(mcArray, mcHeader);
  //   }
  
  if(!fCuts->IsEventSelected(aod)) {
    if(fCuts->GetWhyRejection()==1) // rejected for pileup
      fNentries->Fill(13);
    if(fSys==1 && (fCuts->GetWhyRejection()==2 || fCuts->GetWhyRejection()==3)) fNentries->Fill(15);
    if(fCuts->GetWhyRejection()==7) fNentries->Fill(17);
    if(!fReadMC){
      //      Printf("Event rejected");
      return;
    }    
  }
  fhistMonitoring->Fill(1);
  
 // Check the Nb of SDD clusters
  
  //   if (fIsRejectSDDClusters) { 
  //     Bool_t skipEvent = kFALSE;
  //     Int_t ntracks = 0;
  //     if (aod) ntracks = aod->GetNumberOfTracks();
  //     for(Int_t itrack=0; itrack<ntracks; itrack++) { // loop on tacks
  //       //    ... get the track
  //       AliAODTrack * track = dynamic_cast<AliAODTrack*>(aod->GetTrack(itrack));
  //       if(!track) AliFatal("Not a standard AOD");
  //       if(TESTBIT(track->GetITSClusterMap(),2) || TESTBIT(track->GetITSClusterMap(),3) ){
  // 	skipEvent=kTRUE;
  // 	fNentries->Fill(16);
  // 	break;
  //       }
  //     }
  //     if (skipEvent) return;
  //   }
  
  
  //   if(fhMultVZEROTPCoutTrackCorr)fhMultVZEROTPCoutTrackCorr->Fill(nTPCout,mTotV0);
  

  // AOD primary vertex
  fprimVtx = (AliAODVertex*)aod->GetPrimaryVertex();
  Bool_t isGoodVtx=kFALSE;
  
  //fprimVtx->Print();
  TString primTitle =fprimVtx->GetTitle();
  if(primTitle.Contains("VertexerTracks") && fprimVtx->GetNContributors()>0) {
    isGoodVtx=kTRUE;
    fNentries->Fill(3);
  }
  fhistMonitoring->Fill(2);

  AliAnalysisVertexingHF *vHF=new AliAnalysisVertexingHF();
  //  fEventCounter++;
  
  ftrackArraySel->Reset();
  ftrackArraySelSoftPi->Reset();
  fnSel=0;
  fnSelSoftPi=0;
  
  ftrackSelStatusProton->Reset();
  ftrackSelStatusKaon->Reset();
  ftrackSelStatusPion->Reset();

  Int_t pdgDg[3]={211,321,2212};

  if(fDebug>=0){   
    //    Printf("%d Lc filtering",lcArray->GetEntriesFast());
    for(Int_t iLcFilt=0;iLcFilt<lcArray->GetEntriesFast();iLcFilt++){
      AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)lcArray->UncheckedAt(iLcFilt);      
      if(d->GetSelectionMap()){
	if(!d->HasSelectionBit(AliRDHFCuts::kLcCuts))		continue;
	//	if(!d->HasSelectionBit(AliRDHFCuts::kLcPID))		continue;
      }     
      AliAODMCParticle* part=0x0;
      if(fReadMC){
	Int_t partind=d->MatchToMC(4122,mcArray,3,pdgDg);
	if(partind>=0){
	  part=(AliAODMCParticle*)mcArray->At(partind);
	}
      }

      if(!(vHF->FillRecoCand(aod,d))) {////Fill the data members of the candidate only if they are empty.
	continue;
      }  

      Bool_t unsetvtx=kFALSE;
      if(!d->GetOwnPrimaryVtx()){
	d->SetOwnPrimaryVtx(fprimVtx);
	unsetvtx=kTRUE;
      }
      Int_t iSel=3,iSelTrackCuts=3,iSelCuts=3,iSelPID=0;
      if(!fCutsXic->IsInFiducialAcceptance(d->Pt(),d->Y(4122)))iSel=0;
      if(!iSel){
 	if(unsetvtx)d->UnsetOwnPrimaryVtx();
 	continue;
      }
      for(Int_t itr=0;itr<3;itr++){
	AliAODTrack *track=(AliAODTrack*)d->GetDaughter(itr);
	Int_t iSelProtonCuts=-1,iSelKaonCuts=-1,iSelPionCuts=-1,iSelSoftPionCuts=-1;
	AliESDtrack *trackESD=SelectTrack(track,iSelProtonCuts,iSelKaonCuts,iSelPionCuts,iSelSoftPionCuts,fESDtrackCutsProton,fESDtrackCutsKaon,fESDtrackCutsPion,fESDtrackCutsSoftPion);
	if(!trackESD){
	  iSelTrackCuts=0;
	  break;
	}
	if(iSelProtonCuts < 0 && iSelKaonCuts < 0 && iSelPionCuts < 0 ){
	  iSelTrackCuts=0;
	  delete trackESD;
	  break;
	} 

	delete trackESD;
	Int_t iSelProton=0,iSelKaon=0,iSelPion=0,iSelSoftPion=0;
	IsSelectedPID(track,iSelPion,iSelKaon,iSelProton,iSelPionCuts,iSelKaonCuts,iSelProtonCuts,kFALSE);	
// 	if(part){
// 	  Printf("daughter %d, PID values: %d, %d, %d",itr,iSelPion,iSelKaon,iSelProton);
// 	}
	if(itr==1&&iSelKaon<=0){
// 	  if(part){
// 	    Printf("Setting PID to 0 at step 1");
// 	  }
	  iSelPID=0;
	  break;
	}
	if((itr==0||itr==2)&&(iSelProton<=0&&iSelPion<=0)){
	  //if(part)Printf("Setting PID to 0 at step 2");
	  iSelPID=0;	  
	  break;
	}
	if(itr==0){
	  if(iSelProton>0)iSelPID++;
	  if(iSelPion>0)iSelPID+=2;
	  //	  if(part)Printf(" PID at step 3 is %d",iSelPID);
	}
	if(itr==2){
	  if(iSelPion<0){if(iSelPID==1||iSelPID==3)iSelPID--;}
	  if(iSelProton<0){if(iSelPID==2||iSelPID==3)iSelPID-=2;}
	  //	  if(part)Printf(" PID at step 4 is %d",iSelPID);
	}
      }

      iSelCuts=fCutsXic->IsSelected(d,AliRDHFCuts::kCandidate,(AliAODEvent*)aod);

      if(part){
	fhistMCSpectrumAccLc->Fill(part->Pt(),kRecoCuts+2);
	if(iSelTrackCuts)fhistMCSpectrumAccLc->Fill(part->Pt(),kRecoCuts+3);
	if(iSelCuts)fhistMCSpectrumAccLc->Fill(part->Pt(),kRecoCuts+4);
	//	Printf(" PID our of loop is %d",iSelPID);
	if(iSelPID)fhistMCSpectrumAccLc->Fill(part->Pt(),kRecoCuts+5);
	if(iSelTrackCuts&&iSelPID&&iSelCuts){
	  fDist12SignalFilter->Fill(d->GetDist12toPrim()*10000.);
	  fCosPointDistrSignalFilter->Fill(d->CosPointingAngle());
	}
	//	Printf("Dist 12, filter x10000: %f",d->GetDist12toPrim()*10000.);
      }
      if(iSelTrackCuts&&iSelPID&&iSelCuts){
	fDist12AllFilter->Fill(d->GetDist12toPrim()*10000.);
	fCosPointDistrAllFilter->Fill(d->CosPointingAngle());
      }
      if(unsetvtx)d->UnsetOwnPrimaryVtx();
    }
  }
  
  Print("\nLoop on tracks\n");
  
  // SELECT TRACKS: could consider to include common tracks (dca) to reduce cpu time (call once propagatetodca)
  for(Int_t itrack=0;itrack<aod->GetNumberOfTracks();itrack++){
    fhistMonitoring->Fill(3);
    Int_t iSelProton=0,iSelKaon=0,iSelPion=0;
    Int_t iSelProtonCuts=-1,iSelKaonCuts=-1,iSelPionCuts=-1,iSelSoftPionCuts=-1;
    AliAODTrack *track = dynamic_cast<AliAODTrack*>(aod->GetTrack(itrack));
    if(!track) AliFatal("Not a standard AOD");
    //    Printf("selecting track");
    AliESDtrack *trackESD=SelectTrack(track,iSelProtonCuts,iSelKaonCuts,iSelPionCuts,iSelSoftPionCuts,fESDtrackCutsProton,fESDtrackCutsKaon,fESDtrackCutsPion,fESDtrackCutsSoftPion);
    
    if(!trackESD)continue;
    fhistMonitoring->Fill(4);
    //    Printf("good track");    

    if(iSelSoftPionCuts>=0){
      ftrackArraySelSoftPi->AddAt(itrack,fnSelSoftPi);
      fnSelSoftPi++;
    }    
    
    if(iSelProtonCuts < 0 && iSelKaonCuts < 0 && iSelPionCuts < 0 ){
      delete trackESD;
      continue;
    } 
    
    ftrackArraySel->AddAt(itrack,fnSel);
  
    // PID SELECTION
    IsSelectedPID(track,iSelPion,iSelKaon,iSelProton,iSelPionCuts,iSelKaonCuts,iSelProtonCuts,kTRUE);
    //    if(itrack%50==0)Printf("Track %d, pt: %f",itrack,track->Pt());
    ftrackSelStatusProton->AddAt(iSelProton,fnSel);
    if(iSelProton>0)fhistMonitoring->Fill(7);
    
    ftrackSelStatusKaon->AddAt(iSelKaon,fnSel);
    if(iSelKaon>0)fhistMonitoring->Fill(6);
    
    ftrackSelStatusPion->AddAt(iSelPion,fnSel);
    if(iSelPion>0)fhistMonitoring->Fill(5);
    
    fnSel++;
    delete trackESD;
  }

  /*
  // Initialize vars for TTree and set addresses if needed
  Float_t var[20];
  TString varNames[20]={"pt","pAngle","lxy","nlxy","ptP","ptK","ptPi","vtxchi2","sigmaVtx","sumd02","dca1","dca2","dca3","nd01","nd02","nd03","d0Lc","cosThetaStar1","cosThetaStar2","flagMC"};
  if(fFillTree){
    for(Int_t k=0;k<20;k++){
      fTreeVar->SetBranchAddress(varNames[k].Data(),&var[k]);
    }
  }*/
  // 
  //  Rescaling of reconstructed to Lc as they were Xic added in the end of the tree
  //
  Float_t var[33];
  Short_t resp;
  TString varNames[34]={"pt","pAngle","lxy","nlxy","ptP","ptK","ptPi","vtxchi2","sigmaVtx","sumd02","dca1","dca2","dca3","nd01","nd02","nd03","d0Lc","cosThetaStar1","cosThetaStar2","m_pKpi","m_piKp","flagMC","w_FromLc_toXic","nSig_TPC_prot_0","nSig_TOF_prot_0","nSig_TPC_pion_0","nSig_TOF_pion_0","nSig_TPC_kaon_1","nSig_TOF_kaon_1","nSig_TPC_prot_2","nSig_TOF_prot_2","nSig_TPC_pion_2","nSig_TOF_pion_2","massHypoFilt_respCuts_respPID"};
  if(fFillTree){
    for(Int_t k=0;k<33;k++){
      fTreeVar->SetBranchAddress(varNames[k].Data(),&var[k]);
    }
    fTreeVar->SetBranchAddress(varNames[33].Data(),&resp);
  }

  // NOW LOOP OVER SELECTED TRACKS
  for(Int_t itrack1=0;itrack1<fnSel;itrack1++){// First loop
    AliAODTrack *track1=(AliAODTrack*)aod->GetTrack(ftrackArraySel->At(itrack1));    
    if(!track1)continue;
    //     if(track1->Charge()<0 && !fLikeSign) continue;
    
    // HERE SHOULD CHECK DISPLACEMENT
    for(Int_t itrack2=itrack1+1;itrack2<fnSel;itrack2++){// Second loop
      AliAODTrack *track2=(AliAODTrack*)aod->GetTrack(ftrackArraySel->At(itrack2));    
      if(!track2)continue;
      //     if(track2->Charge()>0 && !fLikeSign) continue;
      
      //HERE SHOULD CHECK DISPLACEMENT and PAIR PID
      
      //      if(track1->Charge()<0 || track2->Charge()>0) continue;  // this is needed to avoid double-counting of unlike-sign
      Short_t charge=track1->Charge();
      charge+=track2->Charge();
      
      // THIRD LOOP
      for(Int_t itrackThird=itrack2+1;itrackThird<fnSel;itrackThird++){// Third loop
	if(itrackThird==itrack1 || itrackThird==itrack2)continue;// should never happen, line can be removed
	AliAODTrack *track3=(AliAODTrack*)aod->GetTrack(ftrackArraySel->At(itrackThird));    
	if(!track3)continue;       
	UShort_t trid[3];
	TObjArray trackArray(3);
	Int_t massHypothesis=0;// 0 = none, 1= p K pi only, 2 = pi, K, p only, 3=both
	if(charge!=0){// ++- or --+ : fill candidate anyway in the same way (opposite charge in the middle)
	  if(charge<0 && track3->Charge()<0)continue;
	  if(charge>0 && track3->Charge()>0)continue;
	  trackArray.AddAt(track1,0);
	  trackArray.AddAt(track3,1);
	  trackArray.AddAt(track2,2);
	  trid[0]=(UShort_t)track1->GetID();
	  trid[1]=(UShort_t)track3->GetID();
	  trid[2]=(UShort_t)track2->GetID();
	  //	  Printf("charge %d,Charged: %d, %d, %d",charge,track1->Charge(),track3->Charge(),track2->Charge());
	  // PID BASIC SELECTION HERE
	  if(ftrackSelStatusKaon->At(itrackThird)<=0)continue;// the opposite-charge particle must be compatible with the K hypo
	  if(ftrackSelStatusProton->At(itrack1)>0 && ftrackSelStatusPion->At(itrack2)>0){
	    massHypothesis+=1;// can be p K pi
	  }
	  if(ftrackSelStatusProton->At(itrack2)>0 && ftrackSelStatusPion->At(itrack1)>0){
	    massHypothesis+=2;// can be pi K p
	  }	 
	}
	else {
	  if(track3->Charge()==track1->Charge()){// +-+ or -+-
	    trackArray.AddAt(track1,0);
	    trackArray.AddAt(track2,1);
	    trackArray.AddAt(track3,2);
	    trid[0]=(UShort_t)track1->GetID();
	    trid[1]=(UShort_t)track2->GetID();
	    trid[2]=(UShort_t)track3->GetID();
	    //	    Printf("charge %d, Charged: %d, %d, %d",charge,track1->Charge(),track2->Charge(),track3->Charge());
	    // PID BASIC SELECTION HERE
	    if(ftrackSelStatusKaon->At(itrack2)<=0)continue;// the opposite-charge particle must be compatible with the K hypo
	    if(ftrackSelStatusProton->At(itrack1)>0 && ftrackSelStatusPion->At(itrackThird)>0){
	      massHypothesis+=1;// can be p K pi
	    }
	    if(ftrackSelStatusProton->At(itrackThird)>0 && ftrackSelStatusPion->At(itrack1)>0){
	      massHypothesis+=2;// can be pi K p
	    }	    
	  }
	  else {// -++ or +--
	    trackArray.AddAt(track2,0);
	    trackArray.AddAt(track1,1);
	    trackArray.AddAt(track3,2);
	    trid[0]=(UShort_t)track2->GetID();
	    trid[1]=(UShort_t)track1->GetID();
	    trid[2]=(UShort_t)track3->GetID();
	    //	    Printf("charge %d, Charged: %d, %d, %d",charge,track2->Charge(),track1->Charge(),track3->Charge());
	    // PID BASIC SELECTION HERE
	    if(ftrackSelStatusKaon->At(itrack1)<=0)continue;// the opposite-charge particle must be compatible with the K hypo
	    if(ftrackSelStatusProton->At(itrack2)>0 && ftrackSelStatusPion->At(itrackThird)>0){
	      massHypothesis+=1;// can be p K pi
	    }
	    if(ftrackSelStatusProton->At(itrackThird)>0 && ftrackSelStatusPion->At(itrack2)>0){
	      massHypothesis+=2;// can be pi K p
	    }	    
	  }
	}
	
	//      if(track3->Charge()<0 && !fLikeSign) continue;
	if(massHypothesis==0)continue; 
	fhistMonitoring->Fill(8);
	AliAODRecoDecayHF3Prong *io3Prong= new AliAODRecoDecayHF3Prong();		  
	io3Prong->SetNProngsHF(3);
	io3Prong->SetNProngs();
	io3Prong->SetProngIDs(3,trid);	
	io3Prong->SetIsFilled(0);

	//	Printf("Filling reco cand");
	if(!vHF->FillRecoCand(aod,io3Prong)){
	  //	  Printf("Filling of reco cand failed");
	  AliAODVertex *vtx3 = (AliAODVertex*)io3Prong->GetSecondaryVtx();
	  if(vtx3){delete vtx3;vtx3=0;}
	  delete io3Prong;
	  continue;
	}
	fhistMonitoring->Fill(9);
	Double_t candPt=io3Prong->Pt();
	//	  Printf("InvMass %f, decay length: %f)",io3Prong->InvMassLcpKpi(),io3Prong->CtLc());
	if(!fCutsXic->IsInFiducialAcceptance(candPt,io3Prong->Y(4122))){
	  AliAODVertex *vtx3 = (AliAODVertex*)io3Prong->GetSecondaryVtx();
	  if(vtx3){delete vtx3;vtx3=0;}
	  delete io3Prong;
	  continue;
	}
	if(io3Prong->GetReducedChi2()>fMaxVtxChi2Cut){
	  AliAODVertex *vtx3 = (AliAODVertex*)io3Prong->GetSecondaryVtx();
	  if(vtx3){delete vtx3;vtx3=0;}
	  delete io3Prong;
	  continue;	  
	}

	Bool_t recPrimVtx=kFALSE;
	AliAODVertex *origownvtx=0x0;
	if(fRecalPrimVtx && fDebug<0){
	  if(io3Prong->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*io3Prong->GetOwnPrimaryVtx());
	  if(fCutsXic->RecalcOwnPrimaryVtx(io3Prong,aod))recPrimVtx=kTRUE;
	  else fCutsXic->CleanOwnPrimaryVtx(io3Prong,aod,origownvtx);
	}

	Int_t isTrueLambdaCorXic=0;
	AliAODMCParticle *part=0x0;
	if(fReadMC){
	  Int_t partind=io3Prong->MatchToMC(4122,mcArray,3,pdgDg);
	  if(partind>=0){
	    part=(AliAODMCParticle*)mcArray->At(partind);
	    if(part){
        isTrueLambdaCorXic=1;
        Int_t pdgMother_checkQuark = AliVertexingHFUtils::CheckOrigin(mcArray,part,kTRUE);
        if(pdgMother_checkQuark==4) isTrueLambdaCorXic*=4;      // from quark c
        else if(pdgMother_checkQuark==5) isTrueLambdaCorXic*=5; // from quark b
        //
        // check if it is pKpi or piKp
        //
        AliAODTrack* trk_prong = (AliAODTrack*) io3Prong->GetDaughter(0);
        Int_t prLabel = TMath::Abs(trk_prong->GetLabel());
        if(prLabel>0){
          AliAODMCParticle* partMC_prong = (AliAODMCParticle*) mcArray->At(prLabel);
          Int_t pdg_prong = -1;
          if(partMC_prong)  pdg_prong = TMath::Abs(partMC_prong->GetPdgCode());
          if(pdg_prong==2212){      // 1st prong is a proton ---> pKpi
            isTrueLambdaCorXic*=10;
          }
          else if(pdg_prong==211){  // 1st prong is a pion ---> piKp
            isTrueLambdaCorXic*=20;
          }
        }
        else{printf("---> Lc prong label %d\n",prLabel);}
      }
	  }
	  else {
	    partind=io3Prong->MatchToMC(4232,mcArray,3,pdgDg);
	    if(partind>=0)part=(AliAODMCParticle*)mcArray->At(partind);
	    if(part){
        isTrueLambdaCorXic=3;
        Int_t pdgMother_checkQuark = AliVertexingHFUtils::CheckOrigin(mcArray,part,kTRUE);
        if(pdgMother_checkQuark==4) isTrueLambdaCorXic*=4;      // from quark c
        else if(pdgMother_checkQuark==5) isTrueLambdaCorXic*=5; // from quark b
        //
        // check if it is pKpi or piKp
        //
        AliAODTrack* trk_prong = (AliAODTrack*) io3Prong->GetDaughter(0);
        Int_t prLabel = TMath::Abs(trk_prong->GetLabel());
        if(prLabel>0){
          AliAODMCParticle* partMC_prong = (AliAODMCParticle*) mcArray->At(prLabel);
          Int_t pdg_prong = -1;
          if(partMC_prong)  pdg_prong = TMath::Abs(partMC_prong->GetPdgCode());
          if(pdg_prong==2212){      // 1st prong is a proton ---> pKpi
            isTrueLambdaCorXic*=10;
          }
          else if(pdg_prong==211){  // 1st prong is a pion ---> piKp
            isTrueLambdaCorXic*=20;
          }
        }
        else{printf("---> Xic prong label %d\n",prLabel);}
      }
	  }
	  //  static Int_t CheckLcpKpiDecay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab);
	  //  AliVertexingHFUtils::CheckLcpKpiDecay(mcArray,)
	  if(!part) {
	    if(fDebug>=2){
	      for(Int_t k=0;k<3;k++){
		AliAODTrack *trk = (AliAODTrack*)io3Prong->GetDaughter(k);
		Int_t dgLabel = trk->GetLabel();
		if(dgLabel>=0){
		  AliAODMCParticle *partMC=(AliAODMCParticle*)mcArray->At(dgLabel);		
		  Int_t labMom=-1,pdgmum=-1;
		  labMom=partMC->GetMother();
		  if(labMom>=0){
		    AliAODMCParticle *partMumMC=(AliAODMCParticle*)mcArray->At(labMom);		
		    pdgmum=partMumMC->GetPdgCode();
		  }
		  Printf("Daught %d, charge %d, label: %d, pdg: %d, mum Label: %d, pdg: %d",k,trk->Charge(),dgLabel,partMC->GetPdgCode(),labMom,pdgmum);
		  
		}
		else Printf("Daught %d, chrage %d, neg label: %d",k,trk->Charge(),dgLabel);
	      }
	    }
	  }
	}
	//if(fReadMC && isTrueLambdaCorXic==1){
    // MAYBE WE'LL CAHNGE IT !!!
	if(fReadMC && (isTrueLambdaCorXic==1 || isTrueLambdaCorXic==4 || isTrueLambdaCorXic==5)){
	  fhistMCSpectrumAccLc->Fill(part->Pt(),kReco);
	}
  //if(fReadMC && isTrueLambdaCorXic==3){
    //if(fReadMC && isTrueLambdaCorXic==1){
	if(fReadMC && (isTrueLambdaCorXic==3 || isTrueLambdaCorXic==12 || isTrueLambdaCorXic==15)){
	  fhistMCSpectrumAccXic->Fill(part->Pt(),kReco);
	}
	if(fDebug>=0){
	  FillDist12and23(io3Prong,aod->GetMagneticField());
	}
	else {
	  io3Prong->SetDist12toPrim(0.05);  //needed to pass pp filtering cuts
	  io3Prong->SetDist23toPrim(0.05);	  
	}

	Double_t pcand[3];
	io3Prong->PxPyPz(pcand);
	AliAODTrack *trPr;
	Double_t cosThetaStarP1=-2,cosThetaStarP2=-2;
	  
	if(massHypothesis==1 || massHypothesis==3){
	  trPr=(AliAODTrack*)trackArray.At(0);
	  Double_t pprot[3];
	  trPr->PxPyPz(pprot);
	  cosThetaStarP1=CosThetaStar(pcand,pprot,TDatabasePDG::Instance()->GetParticle(4122)->Mass(),TDatabasePDG::Instance()->GetParticle(2212)->Mass());	  
	}
	if(massHypothesis==2 || massHypothesis==3){
	  trPr=(AliAODTrack*)trackArray.At(2);
	  Double_t pprot[3];
	  trPr->PxPyPz(pprot);
	  cosThetaStarP2=CosThetaStar(pcand,pprot,TDatabasePDG::Instance()->GetParticle(4122)->Mass(),TDatabasePDG::Instance()->GetParticle(2212)->Mass());
	}
	var[17]=cosThetaStarP1;
	var[18]=cosThetaStarP2;
	
	
	Int_t isSeleCuts=3,pidreject;
  Int_t resp_onlyPID = 3;
	if(fCutsXic && (fAnalysisType==0 || fAnalysisType==2 || fAnalysisType==3)){// TO BE FIXED FOR CASE 0: SHOULD NOT DELETE, WE MUST ADD A BIT TO THE HISTOGRAM OR DUPLICATE THE HISTO FOR INCLUDING XIC

    if(fCompute_dist12_dist23)  FillDist12and23(io3Prong,aod->GetMagneticField());
      // here cuts + PID are considered
    isSeleCuts=fCutsXic->IsSelected(io3Prong,AliRDHFCuts::kCandidate,(AliAODEvent*)aod);
    //printf("\n==============\n");
    //printf("massHypothesis %d\n",massHypothesis);
    //printf("isSeleCuts     %d\n",isSeleCuts);
    //pidreject = fCutsXic->GetIsSelectedPID();  // 0: none; 1: case pKpi; 2: case piKp; 3: not distinguished pKpi or piKp
	  
    // store info for different selections (massHypothesis filtering, cut selection, PID selection)
    Int_t isPIDused = fCutsXic->GetIsUsePID();
    fCutsXic->SetUsePID(kFALSE);   // disable PID temporarly
    Int_t resp_onlyCuts = fCutsXic->IsSelected(io3Prong,AliRDHFCuts::kCandidate,(AliAODEvent*)aod);
    //Int_t resp_onlyPID = 3;
    if(isPIDused){  // if the PID is supposed to be used, let's restore it in the cutobject
      fCutsXic->SetUsePID(kTRUE);  // restoring PID
      resp_onlyPID = fCutsXic->IsSelected(io3Prong,AliRDHFCuts::kPID,(AliAODEvent*)aod);
    }
      //
      // fill the tree entry with a map
      //  - first two bits: massHypothesis
      //  - 3rd and 4th bit: resp_onlyCuts
      //  - last 3 bits: resp_onlyPID
      //
      resp = SetMapCutsResponse(massHypothesis,resp_onlyCuts,resp_onlyPID);


      //POSTPONED!!!
      //
    //massHypothesis=isSeleCuts&massHypothesis;
      //
      //

    //printf("massHypothesis & isSeleCuts %d\ncandPt %.2f\n==============\n",massHypothesis,candPt);

/*
	  if(fFillTree==1){
	    FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,mcArray);
	  }
	  else if(fFillTree==2){
	    if(isTrueLambdaCorXic){
	      if(candPt<4.){
		if(candPt*1000.-(Int_t)(candPt*1000)<0.3)FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,mcArray);
	      }
	      else FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,mcArray);
	    }
	    else {
	      if(candPt<4.){
		if(candPt*1000.-(Int_t)(candPt*1000)<0.003)FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,mcArray);
	      }
	      else{
		if(candPt*1000.-(Int_t)(candPt*1000)<0.02)FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,mcArray);
	      }
	    }
	  }
	  else if(fFillTree==3){
	    if(isTrueLambdaCorXic){
	      if(candPt<4.){
		if(candPt*1000.-(Int_t)(candPt*1000)<0.03)FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,mcArray);
	      }
	      else FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,mcArray);
	    }
	    else {
	      if(candPt<4.){
		if(candPt*1000.-(Int_t)(candPt*1000)<0.0003)FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,mcArray);
	      }
	      else{
		if(candPt*1000.-(Int_t)(candPt*1000)<0.002)FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,mcArray);
	      }
	    }
	  }
*/

    // fill the tree
    if(fFillTree){
      if(candPt<fpT_down){  // downsampling for low pT
        //if(candPt*1000.-(Int_t)(candPt*1000)<fLowpT_down)     FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,mcArray);
        if(candPt*1000.-(Int_t)(candPt*1000)<fLowpT_down){
          // fill only with true generated particles for MC
          if(fReadMC && isTrueLambdaCorXic>0.5)   FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,mcArray);
          else if(!fReadMC)                       FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,mcArray);
        }   
      } 
      else{   // downsampling for high pT
        //if(candPt*1000.-(Int_t)(candPt*1000)<fHighpT_down)    FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,mcArray);
        if(candPt*1000.-(Int_t)(candPt*1000)<fHighpT_down){
          // fill only with true generated particles for MC
          if(fReadMC && isTrueLambdaCorXic>0.5)   FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,mcArray);
          else if(!fReadMC)                       FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,mcArray);
        }
      }
    }

    // POSTPONED HERE !!!
    //
    //massHypothesis=isSeleCuts&massHypothesis;
    if(fExplore_PIDstdCuts)   massHypothesis=resp_onlyCuts&massHypothesis;  // we do not want the candidates to be filtered by Bayes PID
    else                      massHypothesis=isSeleCuts&massHypothesis;
    //
    //

	  if(!massHypothesis){  // true wheter massHypothesis, after PID and cuts, is null
	    if(recPrimVtx)fCutsXic->CleanOwnPrimaryVtx(io3Prong,aod,origownvtx);	    
	    AliAODVertex *vtx3 = (AliAODVertex*)io3Prong->GetSecondaryVtx();
	    if(vtx3){delete vtx3;  vtx3=0;}
	    delete io3Prong;
	    continue;
	  }
	}
	
	
	Int_t flagSel=FlagCandidateWithVariousCuts(io3Prong,aod,itrack1,itrack2,itrackThird,massHypothesis);
	
	fCosPointDistrAll->Fill(io3Prong->CosPointingAngle());
	fDist12All->Fill(io3Prong->GetDist12toPrim()*10000.);
	fDist23All->Fill(io3Prong->GetDist23toPrim()*10000.);
	
	if(fReadMC){
	  //if(isTrueLambdaCorXic==1 || isTrueLambdaCorXic==3){
        // MAYBE WE'LL CAHNGE IT !!!
	  if(isTrueLambdaCorXic==1 || isTrueLambdaCorXic==3 || isTrueLambdaCorXic==4 || isTrueLambdaCorXic==5 || isTrueLambdaCorXic==12 || isTrueLambdaCorXic==15){
	    fCosPointDistrSignal->Fill(io3Prong->CosPointingAngle());
	    fDist12Signal->Fill(io3Prong->GetDist12toPrim()*10000.);
	    //	    Printf("Dist 12, online x10000: %f",io3Prong->GetDist12toPrim()*10000.); 
	    fDist23Signal->Fill(io3Prong->GetDist23toPrim()*10000.);	    
	    //if(isTrueLambdaCorXic==1)fhistMCSpectrumAccLc->Fill(part->Pt(),kRecoCuts);
	    if(isTrueLambdaCorXic==1 || isTrueLambdaCorXic==4 || isTrueLambdaCorXic==5)fhistMCSpectrumAccLc->Fill(part->Pt(),kRecoCuts);
	    //if(isTrueLambdaCorXic==3)fhistMCSpectrumAccXic->Fill(part->Pt(),kRecoCuts);
	    if(isTrueLambdaCorXic==3 || isTrueLambdaCorXic==12 || isTrueLambdaCorXic==15)fhistMCSpectrumAccXic->Fill(part->Pt(),kRecoCuts);
	  }
	}
	
	
	Double_t diffIP[3], errdiffIP[3],normIP[3],maxIP=0;
	io3Prong->Getd0MeasMinusExpProng(0,aod->GetMagneticField(),diffIP[0],errdiffIP[0]);
	io3Prong->Getd0MeasMinusExpProng(1,aod->GetMagneticField(),diffIP[1],errdiffIP[1]);
	io3Prong->Getd0MeasMinusExpProng(2,aod->GetMagneticField(),diffIP[2],errdiffIP[2]);
	normIP[0]=diffIP[0]/errdiffIP[0];
	maxIP=TMath::Abs(  normIP[0]);
	normIP[1]=diffIP[1]/errdiffIP[1];
	if(TMath::Abs(  normIP[1])>maxIP)maxIP=TMath::Abs(  normIP[1]);
	normIP[2]=diffIP[2]/errdiffIP[2];
	if(TMath::Abs(  normIP[2])>maxIP)maxIP=TMath::Abs(  normIP[2]);
	//Double_t point[7]={candPt,0,io3Prong->DecayLengthXY(),io3Prong->NormalizedDecayLengthXY(),io3Prong->CosPointingAngle(),maxIP,(Double_t)flagSel};  // CHANGE DIMENTION!!! FROM 7 TO 8
	Double_t point[8]={candPt,0,io3Prong->DecayLengthXY(),io3Prong->NormalizedDecayLengthXY(),io3Prong->CosPointingAngle(),maxIP,(Double_t)flagSel,0};  // CHANGE DIMENTION!!! FROM 7 TO 8
	//Double_t pointSigma[9]; // CHANGE DIMENTION!!! FROM 9 TO 10
	Double_t pointSigma[10];
	//for(Int_t k=0;k<7;k++){
	for(Int_t k=0;k<8;k++){
	  pointSigma[k]=point[k];
	}
	//pointSigma[7]=0;
	pointSigma[8]=0;


  //printf("=== Candidate!!! ===\n");
	Double_t mass1=0,mass2=0;
	if(massHypothesis==1 || massHypothesis ==3) {
    //printf("massHypothesis %d\n",massHypothesis);
	  mass1=io3Prong->InvMassLcpKpi();
	  fhistInvMassCheck->Fill(mass1,isTrueLambdaCorXic);
	  point[1]=mass1;
    point[7]=0;
    // this two filling fill the sparse with PID in cut object always
    if(fhSparseAnalysis && !fExplore_PIDstdCuts)  fhSparseAnalysis->Fill(point);
    if(fExplore_PIDstdCuts){
      if(fhSparseAnalysis && ( (massHypothesis&resp_onlyPID)==1 || (massHypothesis&resp_onlyPID)==3 ) )  fhSparseAnalysis->Fill(point);
      Bool_t case_pKpi_ok, case_piKp_ok;
      for(UInt_t i=1; i<=10; i++){  // loop on PID cut combinations to be tested
        point[7] = i;
        case_pKpi_ok = kFALSE;
        case_piKp_ok = kFALSE;
        fCutsXic->ExplorePID(fPidResponse,io3Prong,point[7],case_pKpi_ok,case_piKp_ok);
        //printf("case_pKpi_ok %d, case_piKp_ok %d\n",case_pKpi_ok,case_piKp_ok);
        if(fhSparseAnalysis && case_pKpi_ok)  fhSparseAnalysis->Fill(point);
      }
    }
	  //else  {if(fhSparseAnalysis)  fhSparseAnalysis->Fill(point);}
	}
	if(massHypothesis==2 || massHypothesis ==3){
    //printf("massHypothesis %d\n",massHypothesis);
	  mass2=io3Prong->InvMassLcpiKp();
	  fhistInvMassCheck->Fill(mass2,isTrueLambdaCorXic);
	  point[1]=mass2;
    point[7]=0;
    // this two filling fill the sparse with PID in cut object always
    if(fhSparseAnalysis && !fExplore_PIDstdCuts)  fhSparseAnalysis->Fill(point);
    if(fExplore_PIDstdCuts){
      if(fhSparseAnalysis && ( (massHypothesis&resp_onlyPID)==2 || (massHypothesis&resp_onlyPID)==3 ) )  fhSparseAnalysis->Fill(point);
      Bool_t case_pKpi_ok, case_piKp_ok;
      for(UInt_t i=1; i<=10; i++){  // loop on PID cut combinations to be tested
        point[7] = i;
        case_pKpi_ok = kFALSE;
        case_piKp_ok = kFALSE;
        fCutsXic->ExplorePID(fPidResponse,io3Prong,point[7],case_pKpi_ok,case_piKp_ok);
        //printf("case_pKpi_ok %d, case_piKp_ok %d\n",case_pKpi_ok,case_piKp_ok);
        if(fhSparseAnalysis && case_piKp_ok)  fhSparseAnalysis->Fill(point);
      }
    }
	  //else  {if(fhSparseAnalysis)  fhSparseAnalysis->Fill(point);}
	}
	
	fhistMonitoring->Fill(10);
	if(fAnalysisType==0 || fAnalysisType ==3){
	  if(TMath::Abs(mass1-2.28646)<0.030 || TMath::Abs(mass2-2.28646)<0.030){//Lc mass window
	    // Loop over soft pions
	    Double_t p2=io3Prong->P2();
	    for(Int_t isoft=0;isoft<fnSelSoftPi;isoft++){
	      Int_t indsof=ftrackArraySelSoftPi->At(isoft);
	      if(indsof==itrack1 || indsof==itrack2 || indsof==itrackThird)continue;
	      AliAODTrack *tracksoft=(AliAODTrack*)aod->GetTrack(indsof);    		
	      Double_t psoft[3];
	      tracksoft->PxPyPz(psoft);
	      Double_t psigma[3]={pcand[0]+psoft[0],pcand[1]+psoft[1],pcand[2]+psoft[2]};
	      Double_t cosThetaStarSoftPi=-1.1;
	      cosThetaStarSoftPi=CosThetaStar(psigma,psoft,TDatabasePDG::Instance()->GetParticle(4222)->Mass(),TDatabasePDG::Instance()->GetParticle(211)->Mass());
	      //pointSigma[8]=cosThetaStarSoftPi;
        pointSigma[9]=cosThetaStarSoftPi;
	      Double_t e1,e2;	      
	      if((massHypothesis==1 || massHypothesis ==3)&&TMath::Abs(mass1-2.28646)<0.030){
		e1=TMath::Sqrt(mass1*mass1+p2);
		e2=TMath::Sqrt(0.019479785+tracksoft->Px()*tracksoft->Px()+tracksoft->Py()*tracksoft->Py()+tracksoft->Pz()*tracksoft->Pz());// 0.019479785 =  0.13957*0.13957
		TLorentzVector lsum(tracksoft->Px()+io3Prong->Px(),tracksoft->Py()+io3Prong->Py(),tracksoft->Pz()+io3Prong->Pz(),e1+e2);
		//pointSigma[7]=mass1;//;
    pointSigma[8]=mass1;//;
		Double_t deltaM=lsum.M()-mass1;
		pointSigma[1]=deltaM;	       
		if(fhSparseAnalysisSigma && !fExplore_PIDstdCuts)  fhSparseAnalysisSigma->Fill(pointSigma);
    if(fhSparseAnalysisSigma && fExplore_PIDstdCuts && ( (massHypothesis&resp_onlyPID)==1 || (massHypothesis&resp_onlyPID)==3 ) )  fhSparseAnalysisSigma->Fill(pointSigma);
	      }
	      if((massHypothesis==2 || massHypothesis ==3)&&TMath::Abs(mass2-2.28646)<0.030){
		e1=TMath::Sqrt(mass2*mass2+p2);
		e2=TMath::Sqrt(0.019479785+tracksoft->Px()*tracksoft->Px()+tracksoft->Py()*tracksoft->Py()+tracksoft->Pz()*tracksoft->Pz());// 0.019479785 =  0.13957*0.13957
		TLorentzVector lsum(tracksoft->Px()+io3Prong->Px(),tracksoft->Py()+io3Prong->Py(),tracksoft->Pz()+io3Prong->Pz(),e1+e2);
		//pointSigma[7]=mass2;//lsum.M();
    pointSigma[8]=mass2;//lsum.M();
		Double_t deltaM=lsum.M()-mass2;
		pointSigma[1]=deltaM;
		if(fhSparseAnalysisSigma && !fExplore_PIDstdCuts)  fhSparseAnalysisSigma->Fill(pointSigma);
    if(fhSparseAnalysisSigma && fExplore_PIDstdCuts && ( (massHypothesis&resp_onlyPID)==2 || (massHypothesis&resp_onlyPID)==3 ) )  fhSparseAnalysisSigma->Fill(pointSigma);
	      }	      
	    }
	  }
	}
	
	
	// NOW DELETE VTX AND CANDIDATE
	if(recPrimVtx)fCutsXic->CleanOwnPrimaryVtx(io3Prong,aod,origownvtx);	    
	AliAODVertex *vtx3 = (AliAODVertex*)io3Prong->GetSecondaryVtx();
	if(vtx3){delete vtx3;vtx3=0;}
	delete io3Prong;
      }      
    }
  }
  delete vHF;
  PostData(1,fNentries);
  PostData(2,fCounter);
  PostData(3,fOutput);
  
  return;
}

 




//____________________________________________________________________________
Int_t AliAnalysisTaskSEXicTopKpi::CheckXicpKpiDecay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab)const{
  /// Checks the Xic->pKpi decay channel. Returns 1 for non-resonant decays and 2 for resonant ones, -1 in other cases

  Int_t pdgD=mcPart->GetPdgCode();
  if(TMath::Abs(pdgD)!=4232) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  //Int_t labelFirstDau = mcPart->GetDaughter(0); // old
  Int_t labelFirstDau = mcPart->GetDaughterLabel(0);
  Int_t nKaons=0;
  Int_t nPions=0;
  Int_t nProtons=0;
  Double_t sumPxDau=0.;
  Double_t sumPyDau=0.;
  Double_t sumPzDau=0.;
  Int_t nFoundpKpi=0;

  Int_t codeRes=-1;
  if(nDau==3 || nDau==2){
    for(Int_t iDau=0; iDau<nDau; iDau++){
      Int_t indDau = labelFirstDau+iDau;
      if(indDau<0) return -1;
      AliAODMCParticle* dau=dynamic_cast<AliAODMCParticle*>(arrayMC->At(indDau));
      if(!dau) return -1;
      Int_t pdgdau=dau->GetPdgCode();
      if(TMath::Abs(pdgdau)==321){
	nKaons++;
	sumPxDau+=dau->Px();
	sumPyDau+=dau->Py();
	sumPzDau+=dau->Pz();
	arrayDauLab[nFoundpKpi++]=indDau;
	if(nFoundpKpi>3) return -1;
      }else if(TMath::Abs(pdgdau)==211){
	nPions++;
	sumPxDau+=dau->Px();
	sumPyDau+=dau->Py();
	sumPzDau+=dau->Pz();
	arrayDauLab[nFoundpKpi++]=indDau;
	if(nFoundpKpi>3) return -1;
      }else if(TMath::Abs(pdgdau)==2212){
	nProtons++;
	sumPxDau+=dau->Px();
	sumPyDau+=dau->Py();
	sumPzDau+=dau->Pz();
	arrayDauLab[nFoundpKpi++]=indDau;
	if(nFoundpKpi>3) return -1;
      }else if(TMath::Abs(pdgdau)==313){
	codeRes=TMath::Abs(pdgdau);
	Int_t nResDau=dau->GetNDaughters();
	if(nResDau!=2) return -1;
	//Int_t indFirstResDau=dau->GetDaughter(0); // old
	Int_t indFirstResDau=dau->GetDaughterLabel(0);
	for(Int_t resDau=0; resDau<2; resDau++){
	  Int_t indResDau=indFirstResDau+resDau;
	  if(indResDau<0) return -1;
	  AliAODMCParticle* resdau=dynamic_cast<AliAODMCParticle*>(arrayMC->At(indResDau));
	  if(!resdau) return -1;
	  Int_t pdgresdau=resdau->GetPdgCode();
	  if(TMath::Abs(pdgresdau)==321){
	    sumPxDau+=resdau->Px();
	    sumPyDau+=resdau->Py();
	    sumPzDau+=resdau->Pz();
	    nKaons++;
	    arrayDauLab[nFoundpKpi++]=indResDau;
	    if(nFoundpKpi>3) return -1;
	  }else if(TMath::Abs(pdgresdau)==211){
	    sumPxDau+=resdau->Px();
	    sumPyDau+=resdau->Py();
	    sumPzDau+=resdau->Pz();
	    nPions++;
	    arrayDauLab[nFoundpKpi++]=indResDau;
	    if(nFoundpKpi>3) return -1;
	  }
	}
      }else{
	return -1;
      }
    }
    if(nPions!=1) return -1;
    if(nKaons!=1) return -1;
    if(nProtons!=1) return -1;
    if(TMath::Abs(mcPart->Px()-sumPxDau)>0.001) return -2;
    if(TMath::Abs(mcPart->Py()-sumPyDau)>0.001) return -2;
    if(TMath::Abs(mcPart->Pz()-sumPzDau)>0.001) return -2;
    if(nDau==3) return 1;
    else if(nDau==2){
      if(codeRes==313) return 2;
    }
  }
  return -1;
  
}


///--------
AliESDtrack* AliAnalysisTaskSEXicTopKpi::SelectTrack(AliAODTrack *aodtr, Int_t &isSelProton,Int_t &isSelKaon, Int_t &isSelPion,Int_t &isSelSoftPion,AliESDtrackCuts *cutsProton, AliESDtrackCuts *cutsKaon, AliESDtrackCuts *cutsPion,AliESDtrackCuts *cutsSoftPion){
  
  isSelProton=-1;
  isSelKaon=-1;
  isSelPion=-1;
  isSelSoftPion=-1;


  if(aodtr->GetID()<0)return 0x0;
  Bool_t isFB4=kTRUE;
  if(!(aodtr->TestFilterBit(AliAODTrack::kTrkGlobalNoDCA))){
    isFB4=kFALSE;
    if(fAnalysisType!=0 && fAnalysisType !=3){
      return 0x0;
    }
  }
  
  AliESDtrack *esdTrack=new AliESDtrack(aodtr);
  // set the TPC cluster info
  esdTrack->SetTPCClusterMap(aodtr->GetTPCClusterMap());
  esdTrack->SetTPCSharedMap(aodtr->GetTPCSharedMap());
  esdTrack->SetTPCPointsF(aodtr->GetTPCNclsF());
  // needed to calculate the impact parameters
  Double_t pos[3],cov[6];
  fprimVtx->GetXYZ(pos);
  fprimVtx->GetCovarianceMatrix(cov);
  const AliESDVertex vESD(pos,cov,100.,100);

  esdTrack->RelateToVertex(&vESD,0.,3.);
  
  if(fAnalysisType==0 || fAnalysisType ==3){
    if(cutsSoftPion){
      if(cutsSoftPion->IsSelected(esdTrack)){
	isSelSoftPion=0;
      }
    }
  }
  if(!isFB4 && isSelSoftPion<0){
    delete esdTrack;
    return 0x0;   
  }
  
  //AliESDtrackCuts *esdTrCutsAll=fCuts->GetTrackCuts();
  AliESDtrackCuts *esdTrCutsAll=fCutsXic->GetTrackCuts();
  //  esdTrCutsAll->SetMinDCAToVertexXYPtDep("0.0025*TMath::Max(0.,(1-TMath::Floor(TMath::Abs(pt)/2.)))");
  //esdTrCutsAll->SetMaxDCAToVertexXY(0.15);
  //esdTrCutsAll->SetMaxDCAToVertexZ(0.25);
  AliESDtrackCuts::ITSClusterRequirement spdreq=esdTrCutsAll->GetClusterRequirementITS(AliESDtrackCuts::kSPD);

  if(fApplykFirst){
    if(aodtr->Pt()<fMaxPtTrackkFirst){
      esdTrCutsAll->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    }
  }
  if(isFB4){
    if(esdTrCutsAll->IsSelected(esdTrack)){
      isSelProton=0;
      isSelKaon=0;
      isSelPion=0;
    }    
    else if(isSelSoftPion<0){
      if(fApplykFirst)esdTrCutsAll->SetClusterRequirementITS(AliESDtrackCuts::kSPD,spdreq);
      delete esdTrack;
      esdTrack=0x0;
    }
  }
  if(fApplykFirst)esdTrCutsAll->SetClusterRequirementITS(AliESDtrackCuts::kSPD,spdreq);

  //applying ESDtrackCut
  // if(cutsProton->IsSelected(esdTrack))isSelProton=1; 
  //   if(cutsKaon){
  //     if(cutsKaon->IsSelected(esdTrack))isSelKaon=1; 
  //   }
  //   else isSelKaon=isSelProton;
  //   if(cutsPion){
  //     if(cutsPion->IsSelected(esdTrack))isSelPion=1; 
  //   }
  //   else isSelPion=isSelProton;
  
  //   if(isSelProton==0 && isSelKaon==0 && isSelPion==0){
  //     delete esdTrack; esdTrack=0x0;
  //   }
  return esdTrack;
}



void AliAnalysisTaskSEXicTopKpi::IsSelectedPID(AliAODTrack *track,Int_t &iSelPion,Int_t &iSelKaon,Int_t &iSelProton,const Int_t iSelPionCuts,const Int_t iSelKaonCuts,const Int_t iSelProtonCuts,Bool_t fillHistos){

  iSelProton=0;
  iSelKaon=0;
  iSelPion=0;
  // TOF PID SELECTION
  AliPIDResponse::EDetPidStatus status = fPidResponse->CheckPIDStatus(AliPIDResponse::kTOF,track);
  
  Double_t trpt=-1;
  if(fillHistos)trpt=track->Pt();
  
  if (status == AliPIDResponse::kDetPidOk){
    if(iSelProtonCuts>=0){
      Double_t nsigma=fPidResponse->NumberOfSigmasTOF(track,AliPID::kProton);
      if(fillHistos)fnSigmaPIDtofProton->Fill(trpt,nsigma);
      //      Printf("nsigma Proton TOF: %f",nsigma);
      if(-3.<=nsigma&&nsigma<=3)iSelProton++;
      else iSelProton--;
    }   
    if(iSelKaonCuts>=0){
      Double_t nsigma=fPidResponse->NumberOfSigmasTOF(track,AliPID::kKaon);
      if(fillHistos)fnSigmaPIDtofKaon->Fill(trpt,nsigma);
      if(-3.<=nsigma&&nsigma<=3)iSelKaon++;
      else iSelKaon--;
    }   
    if(iSelPionCuts>=0){
      Double_t nsigma=fPidResponse->NumberOfSigmasTOF(track,AliPID::kPion);
      //	Printf("nsigma Pion TOF: %f",nsigma);
      if(fillHistos)fnSigmaPIDtofPion->Fill(trpt,nsigma);
      if(-3.<=nsigma&&nsigma<=3)iSelPion++;
      else iSelPion--;
    }    
  }
  else {
    if(fillHistos && trpt>0.350 && (iSelProtonCuts>=0 || iSelKaonCuts >=0 || iSelPionCuts >=0)){
      fnSigmaPIDtofPion->Fill(trpt,-30);
      fnSigmaPIDtofKaon->Fill(trpt,-30);
      fnSigmaPIDtofProton->Fill(trpt,-30);      
    }
  }

    //    TPC PID SELECTION
  status = fPidResponse->CheckPIDStatus(AliPIDResponse::kTPC,track);
  if (status == AliPIDResponse::kDetPidOk){
    if(iSelProtonCuts>=0){
      Double_t nsigma=fPidResponse->NumberOfSigmasTPC(track,AliPID::kProton);
      if(fillHistos)fnSigmaPIDtpcProton->Fill(trpt,nsigma);
      //      Printf("nsigma Proton TPC: %f",nsigma);
      if(-3.<=nsigma&&nsigma<=3)iSelProton++;
      else iSelProton--;
    }   
    if(iSelKaonCuts>=0){
      Double_t nsigma=fPidResponse->NumberOfSigmasTPC(track,AliPID::kKaon);
      if(fillHistos)fnSigmaPIDtpcKaon->Fill(trpt,nsigma);
      if(-3.<=nsigma&&nsigma<=3)iSelKaon++;
	else iSelKaon--;
    }   
    if(iSelPionCuts>=0){
      Double_t nsigma=fPidResponse->NumberOfSigmasTPC(track,AliPID::kPion);
      if(fillHistos)fnSigmaPIDtpcPion->Fill(trpt,nsigma);
      //	Printf("nsigma Pion TPC: %f",nsigma);
      if(-3.<=nsigma&&nsigma<=3)iSelPion++;
      else iSelPion--;
    }   
  }
  else {
    if(fillHistos && trpt>0.150 && (iSelProtonCuts>=0 || iSelKaonCuts >=0 || iSelPionCuts >=0)){
      fnSigmaPIDtpcPion->Fill(trpt,-30);
      fnSigmaPIDtpcKaon->Fill(trpt,-30);
      fnSigmaPIDtpcProton->Fill(trpt,-30);
    }
  }

}


//__________________________________________________________________
void AliAnalysisTaskSEXicTopKpi::FillDist12and23(AliAODRecoDecayHF3Prong *pr,Double_t magfield){
  
  if (!fVertexerTracks) fVertexerTracks = new AliVertexerTracks(magfield);
  Double_t pos[3],cov[6];
  fprimVtx->GetXYZ(pos);
//   fprimVtx->GetCovarianceMatrix(cov);
//   const AliESDVertex vESD(pos,cov,100.,100);
  
  AliESDtrack **esdTrack=new AliESDtrack*[3];
  Double_t d0z0[2],covd0z0[3];
  for(Int_t j=0;j<3;j++){
    AliAODTrack *tr=(AliAODTrack*)pr->GetDaughter(j);
    esdTrack[j]=new AliESDtrack(tr);
    // needed to calculate the impact parameters
    esdTrack[j]->PropagateToDCA(fprimVtx,magfield,kVeryBig,d0z0,covd0z0);
  }
  TObjArray *twoTrackArray=new TObjArray(2);
  twoTrackArray->AddAt(esdTrack[0],0);
  twoTrackArray->AddAt(esdTrack[1],1);

  AliESDVertex *vertexESD = (AliESDVertex*)fVertexerTracks->VertexForSelectedESDTracks(twoTrackArray);
  Double_t dispersion;//,chi2perNDF;
  //  vertexESD->GetXYZ(pos); // position
  //  vertexESD->GetCovMatrix(cov); //covariance matrix
  //  chi2perNDF = vertexESD->GetChi2toNDF();
  //  dispersion = vertexESD->GetDispersion();
  Double_t dist12=TMath::Sqrt((vertexESD->GetX()-pos[0])*(vertexESD->GetX()-pos[0])+(vertexESD->GetY()-pos[1])*(vertexESD->GetY()-pos[1])+(vertexESD->GetZ()-pos[2])*(vertexESD->GetZ()-pos[2]));
  pr->SetDist12toPrim(dist12);
  delete vertexESD; vertexESD=NULL;
  
  esdTrack[1]->PropagateToDCA(fprimVtx,magfield,kVeryBig,d0z0,covd0z0);
  twoTrackArray->AddAt(esdTrack[2],0);
  twoTrackArray->AddAt(esdTrack[1],1);

  vertexESD = (AliESDVertex*)fVertexerTracks->VertexForSelectedESDTracks(twoTrackArray);

  //  Double_t dca[3]={dcap1n1,dcap2n1,dcap1p2};
  Double_t dist23=TMath::Sqrt((vertexESD->GetX()-pos[0])*(vertexESD->GetX()-pos[0])+(vertexESD->GetY()-pos[1])*(vertexESD->GetY()-pos[1])+(vertexESD->GetZ()-pos[2])*(vertexESD->GetZ()-pos[2]));
  pr->SetDist23toPrim(dist23);
  delete vertexESD; vertexESD=NULL;

  for(Int_t j=0;j<3;j++){
    delete esdTrack[j];
  }
  delete [] esdTrack;
  return;
}

//________________________________________________________________________
Int_t AliAnalysisTaskSEXicTopKpi::FlagCandidateWithVariousCuts(AliAODRecoDecayHF3Prong *pr,AliAODEvent *aod,Int_t itrack1,Int_t itrack2,Int_t itrack3,Int_t massHypo){
  Int_t sellevel=0;
  Double_t ptprong[3];
  ptprong[0]=pr->PtProng(0);
  ptprong[1]=pr->PtProng(1);
  ptprong[2]=pr->PtProng(2);
  
  // first set, just try harder kine selections
  if(massHypo==1||massHypo==3){
    if(ptprong[0]>1.&&ptprong[1]>0.6)sellevel|=1;
  }
  if(massHypo==2||massHypo==3){
    if(ptprong[2]>1.&&ptprong[1]>0.6)sellevel|=1;
  }


  // STRONGER PID REQUIRES DIFFERENT APPROACH
  
  // second set, improve vtx quality
  Double_t sigmaVert=pr->GetSigmaVert(aod);
  Double_t redchi2=pr->GetReducedChi2();
  Double_t maxdca2prongs=-1;
  //  Double_t maxd0prong=-1;
  Double_t sumd02=0;
  for(Int_t i=0;i<3;i++){
    Double_t d0pr=TMath::Abs(pr->Getd0Prong(i));
    sumd02+=d0pr*d0pr;
    //    if(d0pr>maxd0prong)maxd0prong=d0pr;
    Double_t dca2prongs=pr->GetDCA(i);
    if(dca2prongs>maxdca2prongs)maxdca2prongs=dca2prongs;
  }
  if(sigmaVert<0.030&&redchi2<2.&&maxdca2prongs<0.0350)sellevel|=2;
  if(sigmaVert<0.030&&sumd02>0.0001&&redchi2<1.7&&maxdca2prongs<0.0350)sellevel|=4;
  
  return sellevel;
}
//
void AliAnalysisTaskSEXicTopKpi::FillTree(AliAODRecoDecayHF3Prong *cand,Int_t massHypothesis,Float_t *varPointer, Int_t flagMC,AliAODEvent *aod, AliAODMCParticle* p=0x0, TClonesArray* array_MC=0x0){ 
  varPointer[0]=cand->Pt();
  varPointer[1]=cand->CosPointingAngle();
  varPointer[2]=cand->DecayLengthXY();
  varPointer[3]=cand->NormalizedDecayLengthXY();
    // REMOVE! put prong 0 in position 4, prong 1 in position 5, prong 2 in position 6
  /*
  UInt_t indices_prongs[3];
  if(massHypothesis<=1 || massHypothesis==3){
    varPointer[4]=cand->PtProng(0);
    varPointer[5]=cand->PtProng(1);
    varPointer[6]=cand->PtProng(2);
    indices_prongs[0]=0;
    indices_prongs[1]=1;
    indices_prongs[2]=2;
  }
  else if(massHypothesis==2){
    varPointer[4]=cand->PtProng(2);
    varPointer[5]=cand->PtProng(1);
    varPointer[6]=cand->PtProng(0);
    indices_prongs[0]=2;
    indices_prongs[1]=1;
    indices_prongs[2]=0;
  }
  */
  varPointer[4]=cand->PtProng(0);
  varPointer[5]=cand->PtProng(1);
  varPointer[6]=cand->PtProng(2);

  varPointer[7]=cand->GetReducedChi2();
  varPointer[8]=cand->GetSigmaVert(aod);
  varPointer[9]=0;
  for(Int_t i=0;i<3;i++){
    Double_t d0pr=cand->Getd0Prong(i);
    Double_t diffIP, errdiffIP;
    varPointer[10+i]=d0pr;
    varPointer[9]+=d0pr*d0pr;	    
    cand->Getd0MeasMinusExpProng(i,aod->GetMagneticField(),diffIP,errdiffIP);
    varPointer[13+i]=diffIP/errdiffIP;;
  }	     
  varPointer[16]=0.;//temporary
  varPointer[21]=flagMC;	

  // add mass (mfaggin)
    // old
  //Double_t mass_pKpi = -1.;
  //Double_t mass_piKp = -1.;
  //if(massHypothesis==3){
  //  mass_pKpi = cand->InvMassLcpKpi();
  //  mass_piKp = cand->InvMassLcpiKp();
  //}
  //else if(massHypothesis==1){ // p K pi
  //  mass_pKpi = cand->InvMassLcpKpi();
  //}
  //else if(massHypothesis==2){ // pi K p
  //  mass_piKp = cand->InvMassLcpiKp();
  //}
    // modified
  Double_t mass_pKpi = cand->InvMassLcpKpi();
  Double_t mass_piKp = cand->InvMassLcpiKp();

  varPointer[19] = mass_pKpi;
  varPointer[20] = mass_piKp;

  // add info about weight to "convert" reconstructed true LC in Xic
  //if(flagMC!=1 || !p){  // generated particle associated to this reconstructed one is not a Lc or is absent
  if(flagMC<0.5 || !p){  // generated particle associated to this reconstructed one is not a Lc or is absent
    varPointer[22]=1.;
  }
  ///////////////////////////////////////////////////////////
  // flagMC==1 means that the reconstructed particle is connected to a generated Lc
  // flagMC==4(5) means that the reconstructed particle is connected to a (non-)prompt generated Lc with found quark
  // ---
  // UPDATE
  // Information about generated pKpi (*=10) or piKp (*=20), therefore the generated Lc have
  //  - flagMC = 40: prompt Lc decaying in pKpi
  //  - flagMC = 80: prompt Lc decaying in piKp
  //  - flagMC = 50: non-prompt Lc decaying in pKpi
  //  - flagMC = 100: non-prompt Lc decaying in piKp
  // ---
  ///////////////////////////////////////////////////////////
  //else if(flagMC==1 && p && array_MC){ // flagMC==1 means that the reconstructed particle is connected to a generated Lc
  else if(flagMC>0.5 /*&& flagMC<5.5*/ && p && array_MC){ 
    //Int_t index_firstProng = p->GetDaughter(0); // old
    Int_t index_firstProng = p->GetDaughterLabel(0);
    AliAODMCParticle *mc_firstProng=(AliAODMCParticle*)array_MC->At(index_firstProng);
    varPointer[22]=Weight_fromLc_toXic(p,mc_firstProng);
  }
  else  varPointer[22]=999.;

  // store the PID variables (nSigma)
  Float_t nSigma_TPC_prot_0=99, nSigma_TOF_prot_0=99, nSigma_TPC_pion_0=99, nSigma_TOF_pion_0=99;
  Float_t nSigma_TPC_kaon_1=99, nSigma_TOF_kaon_1=99; 
  Float_t nSigma_TPC_prot_2=99, nSigma_TOF_prot_2=99, nSigma_TPC_pion_2=99, nSigma_TOF_pion_2=99;
  AliPIDResponse::EDetPidStatus status_TPC_0, status_TPC_1, status_TPC_2, status_TOF_0, status_TOF_1, status_TOF_2;
  status_TPC_0 = fPidResponse->CheckPIDStatus(AliPIDResponse::kTPC,(AliAODTrack*)cand->GetDaughter(0));
  status_TOF_0 = fPidResponse->CheckPIDStatus(AliPIDResponse::kTOF,(AliAODTrack*)cand->GetDaughter(0));
  status_TPC_1 = fPidResponse->CheckPIDStatus(AliPIDResponse::kTPC,(AliAODTrack*)cand->GetDaughter(1));
  status_TOF_1 = fPidResponse->CheckPIDStatus(AliPIDResponse::kTOF,(AliAODTrack*)cand->GetDaughter(1));
  status_TPC_2 = fPidResponse->CheckPIDStatus(AliPIDResponse::kTPC,(AliAODTrack*)cand->GetDaughter(2));
  status_TOF_2 = fPidResponse->CheckPIDStatus(AliPIDResponse::kTOF,(AliAODTrack*)cand->GetDaughter(2));
  if(status_TPC_0 == AliPIDResponse::kDetPidOk){
    nSigma_TPC_pion_0 = fPidResponse->NumberOfSigmasTPC((AliAODTrack*)cand->GetDaughter(0),AliPID::kPion);
    nSigma_TPC_prot_0 = fPidResponse->NumberOfSigmasTPC((AliAODTrack*)cand->GetDaughter(0),AliPID::kProton);
  }
  if(status_TOF_0 == AliPIDResponse::kDetPidOk){
    nSigma_TOF_pion_0 = fPidResponse->NumberOfSigmasTOF((AliAODTrack*)cand->GetDaughter(0),AliPID::kPion);
    nSigma_TOF_prot_0 = fPidResponse->NumberOfSigmasTOF((AliAODTrack*)cand->GetDaughter(0),AliPID::kProton);
  }
  if(status_TPC_1 == AliPIDResponse::kDetPidOk){
    nSigma_TPC_kaon_1 = fPidResponse->NumberOfSigmasTPC((AliAODTrack*)cand->GetDaughter(1),AliPID::kKaon);
  }
  if(status_TOF_1 == AliPIDResponse::kDetPidOk){
    nSigma_TOF_kaon_1 = fPidResponse->NumberOfSigmasTOF((AliAODTrack*)cand->GetDaughter(1),AliPID::kKaon);
  }
  if(status_TPC_2 == AliPIDResponse::kDetPidOk){
    nSigma_TPC_pion_2 = fPidResponse->NumberOfSigmasTPC((AliAODTrack*)cand->GetDaughter(2),AliPID::kPion);
    nSigma_TPC_prot_2 = fPidResponse->NumberOfSigmasTPC((AliAODTrack*)cand->GetDaughter(2),AliPID::kProton);
  }
  if(status_TOF_2 == AliPIDResponse::kDetPidOk){
    nSigma_TOF_pion_2 = fPidResponse->NumberOfSigmasTOF((AliAODTrack*)cand->GetDaughter(2),AliPID::kPion);
    nSigma_TOF_prot_2 = fPidResponse->NumberOfSigmasTOF((AliAODTrack*)cand->GetDaughter(2),AliPID::kProton);
  }
  varPointer[23] = nSigma_TPC_prot_0;
  varPointer[24] = nSigma_TOF_prot_0;
  varPointer[25] = nSigma_TPC_pion_0;
  varPointer[26] = nSigma_TOF_pion_0;
  varPointer[27] = nSigma_TPC_kaon_1;
  varPointer[28] = nSigma_TOF_kaon_1;
  varPointer[29] = nSigma_TPC_prot_2;
  varPointer[30] = nSigma_TOF_prot_2;
  varPointer[31] = nSigma_TPC_pion_2;
  varPointer[32] = nSigma_TOF_pion_2;

  fTreeVar->Fill();
}

//________________________________________________________________________
Double_t AliAnalysisTaskSEXicTopKpi::CosThetaStar(Double_t mumVector[3],Double_t daughtVector[3],Double_t massMum,Double_t massDaught){
  
  Double_t mumP2=mumVector[0]*mumVector[0]+mumVector[1]*mumVector[1]+mumVector[2]*mumVector[2];
  Double_t mumP=TMath::Sqrt(mumP2);
  Double_t eMum=TMath::Sqrt(mumP2+massMum*massMum);
  Double_t daughtP2=daughtVector[0]*daughtVector[0]+daughtVector[1]*daughtVector[1]+daughtVector[2]*daughtVector[2];
  Double_t eDaugh=TMath::Sqrt(daughtP2+massDaught*massDaught);
  Double_t plLab=(mumVector[0]*daughtVector[0]+mumVector[1]*daughtVector[1]+mumVector[2]*daughtVector[2])/mumP;
  Double_t beta = mumP/eMum;
  Double_t gamma = eMum/massMum;
  Double_t plStar=gamma*(plLab-beta*eDaugh);
  Double_t daughtpT2=daughtP2-plLab*plLab;
  return plStar/TMath::Sqrt(plStar*plStar+daughtpT2);  
  
}

//________________________________________________________________________
void AliAnalysisTaskSEXicTopKpi::Terminate(Option_t */*option*/)
{
  return;
}

//________________________________________________________________________
Double_t AliAnalysisTaskSEXicTopKpi::Weight_fromLc_toXic(AliAODMCParticle* p, AliAODMCParticle* prong)
{
  // Function to calculate weight to treat reco true Lc as Xic (mfaggin)
  Double_t cTau_Lc  = 59.9E-04;  // [cTau]=cm
  Double_t cTau_Xic = 132E-04;  // [cTau]=cm  
  Double_t mass_Lc  = 2.28646;   // [m]=GeV/c^2
  Double_t mass_Xic = 2.46787;  // [m]=GeV/c^2

  TLorentzVector vecLc;
  TLorentzVector vecXic;
  vecLc.SetXYZM(p->Px(),p->Py(),p->Pz(),mass_Lc);   // [p]=GeV/c
  vecXic.SetXYZM(p->Px(),p->Py(),p->Pz(),mass_Xic);  // [p]=GeV/c

  /*
  // calc. decay length
  Double_t X_prodVtx_Lc = p->Xv();  // coordinates of vertex where the Lc is produced (primary vertex for prompt)
  Double_t Y_prodVtx_Lc = p->Yv();
  Double_t Z_prodVtx_Lc = p->Zv();
  Double_t X_prodVtx_prong = prong->Xv(); // coordinates of space point where the Lc decays
  Double_t Y_prodVtx_prong = prong->Yv();
  Double_t Z_prodVtx_prong = prong->Zv();
  Double_t decLength = TMath::Sqrt( pow(X_prodVtx_Lc-X_prodVtx_prong,2) + pow(Y_prodVtx_Lc-Y_prodVtx_prong,2) + pow(Z_prodVtx_Lc-Z_prodVtx_prong,2) );

  cout << "****** X_prodVtx_Lc: " << X_prodVtx_Lc << "   Y_prodVtx_Lc: " << Y_prodVtx_Lc << "   Z_prodVtx_Lc: " << Z_prodVtx_Lc << endl;
  cout << "****** X_prodVtx_prong: " << X_prodVtx_prong << "   Y_prodVtx_prong: " << Y_prodVtx_prong << "   Z_prodVtx_prong: " << Z_prodVtx_prong << endl; 
  Double_t ct = decLength*vecLc.E()/p->P();  // [decaylength]=cm, [m/p] = /

  // ct
  cout << "****** prong->Tv(): " << prong->Tv() << "     p->Tv(): " << p->Tv() << endl;
  */
  Double_t ct1 = 2.99792458E+10*(prong->Tv()-p->Tv());  // [c] = cm/s

  //cout << "************** ct: " << ct << "    ct1 " << ct1 << endl;

  Double_t gamma_Lc = vecLc.Gamma();
  Double_t gamma_Xic = vecXic.Gamma();

  //Double_t exp_Xic = TMath::Exp(-ct/(gamma_Xic*cTau_Xic));
  //Double_t exp_Lc  = TMath::Exp(-ct/(gamma_Lc*cTau_Lc));  
  Double_t exp_Xic = TMath::Exp(-ct1/(gamma_Xic*cTau_Xic));
  Double_t exp_Lc  = TMath::Exp(-ct1/(gamma_Lc*cTau_Lc));

  return exp_Xic/exp_Lc; 
}
//________________________________________________________________________
Short_t AliAnalysisTaskSEXicTopKpi::SetMapCutsResponse(Int_t massHypo_filtering, Int_t response_onlyCuts, Int_t response_onlyPID)
{
  //
  // Create the following map:
  //    - massHypothesis in 0 and 1 bits (00: massHypothesis = 0; 01: massHypothesis = 1; 10: massHypothesis = 2; 11: massHypothesis = 3)
  //    - resp_onlyCuts  in 2 and 3 bits (same logic, but two bits left)
  //    - resp_onlyPID   in 4 and 5 bits (same logic, but two bits left further) 
  //

  Short_t returnValue = 0;
  // massHypthesis
  if(massHypo_filtering==1)       SETBIT(returnValue,0);
  else if(massHypo_filtering==2)  SETBIT(returnValue,1);
  else if(massHypo_filtering==3){
    SETBIT(returnValue,0);
    SETBIT(returnValue,1);
  }
  // resp_onlyCuts
  if(response_onlyCuts==1)        SETBIT(returnValue,2);
  else if(response_onlyCuts==2)   SETBIT(returnValue,3);
  else if (response_onlyCuts==3)
  {
    SETBIT(returnValue,2);
    SETBIT(returnValue,3);
  }
  //resp_onlyPID
  if(response_onlyPID==1)         SETBIT(returnValue,4);
  else if(response_onlyPID==2)    SETBIT(returnValue,5);
  else if(response_onlyPID==3){
    SETBIT(returnValue,4);
    SETBIT(returnValue,5);
  }
  
  return returnValue;
}
