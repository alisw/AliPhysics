#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "TObject.h"
#include "AliVEvent.h"
#include "AliAnalysisTaskdNdEtapp13.h"
#include "AliAnalysisManager.h"
#include "TProof.h"
#include "TROOT.h"
#endif


const Double_t centBinsMultV0MSPD[] = {1e-10, 1., 5., 10., 15., 20., 30., 40., 50., 70., 100};
const Double_t centBinsMultV0M[] =   {0., 1., 5., 10., 15., 20., 30., 40., 50., 70., 100};
const Double_t centBinsMultV0M_Ridge[] = {0., 0.001, 0.01, 0.1, 0.5, 1., 3., 5., 10., 20., 30., 40., 50., 90., 100};
const Double_t centBinsMultRef[] = {0., 1., 5., 10., 15., 20., 30., 40., 50., 70., 100};
const Double_t centBinsMB[] = {0., 100.};




//__________________________________________________________________

AliAnalysisTaskdNdEtapp13 *
AddAnalysisTaskdNdEtaPP13(const Char_t *outfilename = "AnalysisResults.root",
const Char_t *listname = "clist",
//
Float_t etaMin     =-5,          // min eta range to fill in histos
Float_t etaMax     = 5,          // max eta range to fill in histos
Float_t zMin       = -10,         // process events with Z vertex min
Float_t zMax       =  10,         //                     max positions
const char* useCentVar = "V0M",          // centrality variable to use
//
Float_t cutSigNStd  = 1.5,        // cut on weighed distance used to extract signal
Float_t cutSigDPhiS = -1,        // cut on dPhi-phiBent used to extract signal (if negative -> dphi*sqrt(cutSigNStd)
Bool_t  useMC  = kTRUE,          // fill MC info
//
Bool_t doRec  = kFALSE,          // fill data histos from new reco
Bool_t doInj  = kFALSE,          // create Inj. bg
Bool_t doRot  = kFALSE,          // create Rot. bg
//
// specific parameters for reconstruction
float  phiRot      = 3.14159e+00, // angle for bg. generation with rotation
float  injScale    = 1.,//0.7,    // inject injScale*Ncl(Lr1/Lr2) hits
Bool_t scaleDTheta = kTRUE,       // scale dTheta by 1/sin^2(theta) in trackleting
float  nStdDev     = 25.,         // number of st.dev. for tracklet cut to keep
float  dphi        = 0.08,        // dphi window (sigma of tracklet cut)
float  dtht        = 0.025,       // dtheta .... (if negative, abs will be used with additional cut on |dthetaX|, apart from w.distance
  float  phishift    = 0.0045,      // bending shift
  Bool_t remOvl      = kTRUE,
  float  ovlPhiCut   = 0.005,
  float  ovlZetaCut  = 0.05,
  Bool_t checkReconstructables = kFALSE,
  UInt_t trigSel = AliVEvent::kINT7,//kTRUE, // fill histos for reconstructable (needs useMC and doRec)
  Bool_t ridgeBins      = kFALSE,             // VOM percentiles with ridge binning
  Bool_t useBCmod      = kFALSE,              // set Bunch crossing mode 4
  Int_t BCmod4      = 2,              // set Bunch crossing mode 4
  Bool_t phicuts      = kFALSE,       // set cut on affected phi regions
  TString calibfile = "$ALICE_PHYSICS/PWGLF/ppVsMult/dNdEtaPP13/task/V0M_bins_LHC15g3a3.root"

)
{

  if (cutSigDPhiS<0) cutSigDPhiS = TMath::Sqrt(cutSigNStd)*dphi;

  /* init analysis name */
  TString analysisName = "dNdEtapp13";

  /* check analysis manager */
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("", "cannot get analysis manager");
    return NULL;
  }

  /* check input event handler */
  if (!mgr->GetInputEventHandler()) {
    Error("", "cannot get input event handler");
    return NULL;
  }

  /* get common input data container */
  AliAnalysisDataContainer *inputc = mgr->GetCommonInputContainer();
  if (!inputc) {
    Error("", "cannot get common input container");
    return NULL;
  }

  /* create output data container */
  AliAnalysisDataContainer *outputc1 = mgr->CreateContainer(listname, TList::Class(), AliAnalysisManager::kOutputContainer, outfilename);

  if (!outputc1) {
    Error("", "cannot create output container \"clist\"");
    return NULL;
  }


  /*  create task and connect input/output */
  AliAnalysisTaskdNdEtapp13 *task = new AliAnalysisTaskdNdEtapp13("AliAnalysisTaskdNdEtapp13");
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, inputc);
  mgr->ConnectOutput(task, 1, outputc1);
  //__________________________________________________________________
  //__________________________________________________________________
  //__________________________________________________________________

  /* configure task */

  //__________________________________________________________________
  //__________________________________________________________________
  //__________________________________________________________________

  //

  TFile *fcalib = TFile::Open(calibfile);
  TH1D * hMCcalib2 = (TH1D*)fcalib->Get("h3") ;

  task->SetCalibfilePath(calibfile);
  task->SetCalibHisto(hMCcalib2);

  task->SetUseCentralityVar(useCentVar);

  const Double_t *centBins = NULL;
  //const Float_t *centBins = NULL;

  Int_t nCentBins = 0;
  TString strCentr = useCentVar;
  if (strCentr == "MB"){
    centBins = centBinsMB;
    nCentBins = sizeof(centBinsMB) / 8 - 1;
  }
  else if (strCentr == "V0M" && ridgeBins){
    centBins = centBinsMultV0M_Ridge;
    nCentBins = sizeof(centBinsMultV0M_Ridge) / 8 - 1;
  }
  else if (strCentr == "V0M" || strCentr == "V0av"){
    centBins = centBinsMultV0M;
    nCentBins = sizeof(centBinsMultV0M) / 8 - 1;
  }
  else if(strCentr == "RefMult08"){
    centBins = centBinsMultRef;
    nCentBins = sizeof(centBinsMultRef) / 8 - 1;
  }
  else if(strCentr == "SPDClusters1"){
    centBins = centBinsMultV0MSPD;
    nCentBins = sizeof(centBinsMultV0MSPD) / 8 - 1;
  }
  else if(strCentr == "CL0"){
    centBins = centBinsMultV0M;
    nCentBins = sizeof(centBinsMultV0M) / 8 - 1;
  }
  else if(strCentr == "SPDTracklets08"){
    centBins = centBinsMultV0M;
    nCentBins = sizeof(centBinsMultV0M) / 8 - 1;
  }
  else if(strCentr == "SPDTracklets08to15"){
    centBins = centBinsMultV0M;
    nCentBins = sizeof(centBinsMultV0M) / 8 - 1;
  }




  task->SetCentPercentiles(centBins, nCentBins);
  task->SetTriggerSelection(trigSel);
  //
  task->SetDoNormalReco(doRec);
  task->SetDoInjection(doInj);
  task->SetDoRotation(doRot);
  //
  task->SetUseMC(useMC);
  task->SetCheckReconstructables(checkReconstructables);
  //
  task->SetEtaMin(etaMin);
  task->SetEtaMax(etaMax);
  task->SetZVertexMin(zMin);
  task->SetZVertexMax(zMax);
  //
  task->SetDPhiSCut(cutSigDPhiS);
  task->SetNStdCut(cutSigNStd);
  //
  task->SetScaleDThetaBySin2T(scaleDTheta);
  task->SetNStdDev(nStdDev);
  task->SetPhiWindow(dphi);
  task->SetThetaWindow(dtht);
  task->SetPhiShift(phishift);
  task->SetPhiOverlapCut(ovlPhiCut);
  task->SetZetaOverlapCut(ovlZetaCut);
  task->SetPhiRot(phiRot);
  task->SetInjScale(injScale);
  task->SetRemoveOverlaps(remOvl);
  task->SetUseBCMod(useBCmod);
  task->SetBCMod(BCmod4);
  task->SetCutOnPhi(phicuts);



  //
  //    task->Dump();
  return task;
}
