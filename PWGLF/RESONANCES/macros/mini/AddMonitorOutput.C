/***************************************************************************
              fbellini@cern.ch - last modified on 09/04/2013

 Macro to add monitoring histograms for track and PID cuts to the rsn task 
 Tuned for monitoring TOF KStar analysis of PbPb 2010 data

Options ("opt" argument):
- dim1 --> use TH1 only (no p or pt dependence)
- dim3 --> use TH3 (p dependence) for TOF vs TPC pid histos
- NoTOFSIGMA  --> disable the nsigma_TOF vs p histos
- NoTPCSIGMA  --> disable the nsigma_TPC vs p_TPC histos
- NoSIGN --> disable splitting for single track charge
- NoTrackQ  --> disable track quality monitoring (DCA, nCls) histos
***************************************************************************/

#ifndef ALIRSNADDMONITOROUTPUT_C
#define ALIRSNADDMONITOROUTPUT_C

#if !defined (__CINT__) || defined (__CLING__)
#include "AliRsnValueEvent.h"
#include "AliRsnValueDaughter.h"
#include "AliRsnListOutput.h"
#include "AliRsnLoopDaughter.h"
#endif

void AddMonitorOutput(Bool_t useMCMon, TObjArray *mon=0,TString opt="NoSIGN",AliRsnLoopDaughter *lm=0)
{
  //Set options
  Bool_t useTH1 = opt.Contains("dim1");
  Bool_t useTH3 = opt.Contains("dim3");
  Bool_t monitorSign = !opt.Contains("NoSIGN");
  Bool_t monitorTrackQ = !opt.Contains("NoTrackQ");
  Bool_t monitorTPCpid = !opt.Contains("NoTPCSIGMA");
  Bool_t monitorTOFpid = !opt.Contains("NoTOFSIGMA");
  
  // Multiplicity/centrality
  AliRsnValueEvent *multi = new AliRsnValueEvent("multi",AliRsnValueEvent::kMult);
  multi->SetBins(0.0, 100.0, 1);
  // Momentum
  AliRsnValueDaughter *axisMomTPC = new AliRsnValueDaughter("pTPC", AliRsnValueDaughter::kPtpc);
  axisMomTPC->SetBins(0.0, 10.0, 0.02);
  AliRsnValueDaughter *axisMomP = new AliRsnValueDaughter("p", AliRsnValueDaughter::kP);
  axisMomP->SetBins(0.0, 10.0, 0.02);
  // Momentum Pt
  AliRsnValueDaughter *axisMomPt = new AliRsnValueDaughter("pt", AliRsnValueDaughter::kPt);
  axisMomPt->SetBins(0.0,10.0,0.02);
  // Eta
  AliRsnValueDaughter *axisMomEta = new AliRsnValueDaughter("eta", AliRsnValueDaughter::kEta);
  axisMomEta->SetBins(-1.0, 1.0, 0.1);
  //ITS clusters
  AliRsnValueDaughter *axisITScls = new AliRsnValueDaughter("ITScls", AliRsnValueDaughter::kNITSclusters);
  axisITScls->SetBins(0.0, 12.0, 1.0);
  //TPC clusters
  AliRsnValueDaughter *axisTPCcls = new AliRsnValueDaughter("TPCcls", AliRsnValueDaughter::kNTPCclusters);
  axisTPCcls->SetBins(0.0, 160.0, 1.0);
  //TPC crossed rows
  AliRsnValueDaughter *axisTPCcrossedRows = new AliRsnValueDaughter("TPCcrossedRows", AliRsnValueDaughter::kNTPCcrossedRows);
  axisTPCcrossedRows->SetBins(0.0, 160.0, 1.0);
  //TPC crossed rows / findable clusters
  AliRsnValueDaughter *axisTPCcrossedRows2Fcls = new AliRsnValueDaughter("TPCcrossedRows2Fcls", AliRsnValueDaughter::kNTPCcrossedRowsFclusters);
  axisTPCcrossedRows2Fcls->SetBins(0.0, 2.0, 0.01);
  //ITS chi2
  AliRsnValueDaughter *axisITSchi2 = new AliRsnValueDaughter("ITSchi2", AliRsnValueDaughter::kITSchi2);
  axisITSchi2->SetBins(0.0, 12.0, 0.1);
  //TPC chi2
  AliRsnValueDaughter *axisTPCchi2 = new AliRsnValueDaughter("TPCchi2", AliRsnValueDaughter::kTPCchi2);
  axisTPCchi2->SetBins(0.0, 10.0, 0.1);
  //DCA xy
  AliRsnValueDaughter *axisDCAxy = new AliRsnValueDaughter("DCAxy", AliRsnValueDaughter::kDCAXY);
  axisDCAxy->SetBins(-2.5, 2.5, 0.001);
  //DCA z
  AliRsnValueDaughter *axisDCAz = new AliRsnValueDaughter("DCAz", AliRsnValueDaughter::kDCAZ);
  axisDCAz->SetBins(-5.0, 5.0, 0.01);
  //Charge
  AliRsnValueDaughter *axisCharge = new AliRsnValueDaughter("charge",AliRsnValueDaughter::kCharge);
  axisCharge->SetBins(-1.5, 1.5, 1.0);
  //Phi
  AliRsnValueDaughter *axisPhi = new AliRsnValueDaughter("phi", AliRsnValueDaughter::kPhi);
  axisPhi->SetBins(0.0, 360.0, 1.0);
  AliRsnValueDaughter *axisPhiOuterTPC = new AliRsnValueDaughter("phiOuterTPC", AliRsnValueDaughter::kPhiOuterTPC);
  axisPhiOuterTPC->SetBins(0.0, 360.0, 1.0);
  
  /* dEdx tpc */
  AliRsnValueDaughter *axisSigTPC = new AliRsnValueDaughter("sTPC", AliRsnValueDaughter::kTPCsignal);
  axisSigTPC->SetBins(0.0, 500.0, 2.0);
  // kTPCnsigmaPi
  AliRsnValueDaughter *axisTPCnsigmaPi = new AliRsnValueDaughter("pi", AliRsnValueDaughter::kTPCnsigmaPi);
  axisTPCnsigmaPi->SetBins(-10.,10., 0.1);
  // kTPCnsigmaK
  AliRsnValueDaughter *axisTPCnsigmaK = new AliRsnValueDaughter("K", AliRsnValueDaughter::kTPCnsigmaK);
  axisTPCnsigmaK->SetBins(-10.,10., 0.1);
  // kTPCnsigmaP
  AliRsnValueDaughter *axisTPCnsigmaP = new AliRsnValueDaughter("p", AliRsnValueDaughter::kTPCnsigmaP);
  axisTPCnsigmaP->SetBins(-10.,10., 0.1);

  /* tof signal - time */
  AliRsnValueDaughter *axisSigTOF = new AliRsnValueDaughter("sTOF", AliRsnValueDaughter::kTOFsignal);
  axisSigTOF->SetBins(10000.0, 50000.0, 250);//in ps
  
  // kTOFdeltaPi
  AliRsnValueDaughter *axisTOFdeltaPi = new AliRsnValueDaughter("Dpi", AliRsnValueDaughter::kTOFdeltaPi);
  axisTOFdeltaPi->SetBins(-3000.,3000., 10.);
  // kTOFdeltaK
  AliRsnValueDaughter *axisTOFdeltaK = new AliRsnValueDaughter("DK", AliRsnValueDaughter::kTOFdeltaK);
  axisTOFdeltaK->SetBins(-3000.,3000., 10.);
  // kTOFdeltaP
  AliRsnValueDaughter *axisTOFdeltaP = new AliRsnValueDaughter("Dp", AliRsnValueDaughter::kTOFdeltaP);
  axisTOFdeltaP->SetBins(-3000.,3000., 10.);

  // kTOFnsigmaPi
  AliRsnValueDaughter *axisTOFnsigmaPi = new AliRsnValueDaughter("pi", AliRsnValueDaughter::kTOFnsigmaPi);
  axisTOFnsigmaPi->SetBins(-10.,10., 0.1);
  // kTOFnsigmaK
  AliRsnValueDaughter *axisTOFnsigmaK = new AliRsnValueDaughter("K", AliRsnValueDaughter::kTOFnsigmaK);
  axisTOFnsigmaK->SetBins(-10.,10., 0.1);
  // kTOFnsigmaP
  AliRsnValueDaughter *axisTOFnsigmaP = new AliRsnValueDaughter("p", AliRsnValueDaughter::kTOFnsigmaP);
  axisTOFnsigmaP->SetBins(-10.,10., 0.1);
   
  /****************************************************************/
  /***************       MONITOR AOB           ********************/
  /****************************************************************/
  AliRsnListOutput *outMonitorPTvsMult = new AliRsnListOutput("PtVsMult",AliRsnListOutput::kHistoDefault);
  outMonitorPTvsMult->AddValue(axisMomPt);
  outMonitorPTvsMult->AddValue(multi);
  if (mon) mon->Add(outMonitorPTvsMult);
  if (lm) lm->AddOutput(outMonitorPTvsMult);

  //    if (lm) lm->SetTrueMC(kTRUE);
  /****************************************************************/
  /***************       MONITOR kinematics    ********************/
  /****************************************************************/
  // output: TH1D for momentum
  AliRsnListOutput *outMonitorP = new AliRsnListOutput("P", AliRsnListOutput::kHistoDefault);
  outMonitorP->AddValue(axisMomP);
  if (monitorSign) outMonitorP->AddValue(axisCharge);
  if (mon) mon->Add(outMonitorP);
  if (lm) lm->AddOutput(outMonitorP);
  
  // output:  TH1D for pt
  AliRsnListOutput *outMonitorPt = new AliRsnListOutput("Pt", AliRsnListOutput::kHistoDefault);
  outMonitorPt->AddValue(axisMomPt);
  if (monitorSign) outMonitorPt->AddValue(axisCharge);
  if (mon) mon->Add(outMonitorPt);
  if (lm) lm->AddOutput(outMonitorPt);
  
  // output: TH1D for pseudorapidity
  AliRsnListOutput *outMonitorEta = new AliRsnListOutput("Eta", AliRsnListOutput::kHistoDefault);
  outMonitorEta->AddValue(axisMomEta);
  if (monitorSign) outMonitorEta->AddValue(axisCharge);
  if (mon) mon->Add(outMonitorEta);
  if (lm) lm->AddOutput(outMonitorEta);
  
  // output:  TH1D for phi at vertex
  AliRsnListOutput *outMonitorPhi = new AliRsnListOutput("Phi", AliRsnListOutput::kHistoDefault);
  outMonitorPhi->AddValue(axisPhi);
  if (monitorSign) outMonitorPhi->AddValue(axisCharge);
  if (mon) mon->Add(outMonitorPhi);
  if (lm) lm->AddOutput(outMonitorPhi);
  
  // output:  TH1D for phiOuterTPC at TPC outer radius
  AliRsnListOutput *outMonitorPhiOuterTPC = new AliRsnListOutput("PhiOuterTPC", AliRsnListOutput::kHistoDefault);
  outMonitorPhiOuterTPC->AddValue(axisPhiOuterTPC);
  if (monitorSign) outMonitorPhiOuterTPC->AddValue(axisCharge);
  if (mon) mon->Add(outMonitorPhiOuterTPC);
  if (lm) lm->AddOutput(outMonitorPhiOuterTPC);
  
  // output:  TH2D for phi vs pt
  AliRsnListOutput *outMonitorPhiVsPt = new AliRsnListOutput("PhiVsPt", AliRsnListOutput::kHistoDefault);
  outMonitorPhiVsPt->AddValue(axisMomPt);
  outMonitorPhiVsPt->AddValue(axisPhi);
  if (monitorSign) outMonitorPhiVsPt->AddValue(axisCharge);
  if (mon) mon->Add(outMonitorPhiVsPt);
  if (lm) lm->AddOutput(outMonitorPhiVsPt);

  /****************************************************************/
  /***************      MONITOR TRACK QUALITY  ********************/
  /****************************************************************/
  if (monitorTrackQ) {
    // output: 2D histogram of DCAxy vs pt
    AliRsnListOutput *outMonitorDCAxy = new AliRsnListOutput("DCAxyVsPt", AliRsnListOutput::kHistoDefault);
    if (!useTH1) outMonitorDCAxy->AddValue(axisMomPt);
    outMonitorDCAxy->AddValue(axisDCAxy);
    if (mon) mon->Add(outMonitorDCAxy);
    if (lm) lm->AddOutput(outMonitorDCAxy);

    // output: 2D histogram of DCAz vs P
    AliRsnListOutput *outMonitorDCAz = new AliRsnListOutput("DCAzVsP", AliRsnListOutput::kHistoDefault);
    if (!useTH1) outMonitorDCAz->AddValue(axisMomP);
    outMonitorDCAz->AddValue(axisDCAz);
    if (mon) mon->Add(outMonitorDCAz);
    if (lm) lm->AddOutput(outMonitorDCAz);

    // output: 2D histogram of ITS cls vs pt
    AliRsnListOutput *outMonitorITScls = new AliRsnListOutput("ITSclsVsPt", AliRsnListOutput::kHistoDefault);
    if (!useTH1) outMonitorITScls->AddValue(axisMomPt);
    outMonitorITScls->AddValue(axisITScls);
    if (mon) mon->Add(outMonitorITScls);
    if (lm) lm->AddOutput(outMonitorITScls);

    // output: 2D histogram of TPC cls vs. pt
    AliRsnListOutput *outMonitorTPCcls = new AliRsnListOutput("TPCclsVsPt", AliRsnListOutput::kHistoDefault);
    if (!useTH1) outMonitorTPCcls->AddValue(axisMomPt);
    outMonitorTPCcls->AddValue(axisTPCcls);
    if (mon) mon->Add(outMonitorTPCcls);
    if (lm) lm->AddOutput(outMonitorTPCcls);

    // output: 2D histogram of TPC cls vs. TPC momentum
    AliRsnListOutput *outMonitorTPCclsVsPtpc = new AliRsnListOutput("TPCclsVsPtpc", AliRsnListOutput::kHistoDefault);
    if (!useTH1) outMonitorTPCclsVsPtpc->AddValue(axisMomTPC);
    outMonitorTPCclsVsPtpc->AddValue(axisTPCcls);
    if (mon) mon->Add(outMonitorTPCclsVsPtpc);
    if (lm) lm->AddOutput(outMonitorTPCclsVsPtpc);

    // output: 2D histogram of TPC crossed rows vs. TPC momentum
    AliRsnListOutput *outMonitorTPCcrossedRowsVsPtpc = new AliRsnListOutput("TPCcrossedRowsVsPtpc", AliRsnListOutput::kHistoDefault);
    if (!useTH1) outMonitorTPCcrossedRowsVsPtpc->AddValue(axisMomTPC);
    outMonitorTPCcrossedRowsVsPtpc->AddValue(axisTPCcrossedRows);
    if (mon) mon->Add(outMonitorTPCcrossedRowsVsPtpc);
    if (lm) lm->AddOutput(outMonitorTPCcrossedRowsVsPtpc);

    // output: 2D histogram of TPC crossed rows/Findable cls vs. TPC momentum
    AliRsnListOutput *outMonitorTPCcrossedRows2FclsVsPtpc = new AliRsnListOutput("TPCcrossedRows2FclsVsPtpc", AliRsnListOutput::kHistoDefault);
    if (!useTH1) outMonitorTPCcrossedRows2FclsVsPtpc->AddValue(axisMomTPC);
    outMonitorTPCcrossedRows2FclsVsPtpc->AddValue(axisTPCcrossedRows2Fcls);
    if (mon) mon->Add(outMonitorTPCcrossedRows2FclsVsPtpc);
    if (lm) lm->AddOutput(outMonitorTPCcrossedRows2FclsVsPtpc);

    // output: 2D histogram of ITS chi2 vs pt
    AliRsnListOutput *outMonitorITSchi2 = new AliRsnListOutput("ITSchi2VsPt", AliRsnListOutput::kHistoDefault);
    if (!useTH1) outMonitorITSchi2->AddValue(axisMomPt);
    outMonitorITSchi2->AddValue(axisITSchi2);
    if (mon) mon->Add(outMonitorITSchi2);
    if (lm) lm->AddOutput(outMonitorITSchi2);

    // output: 2D histogram of TPC chi2 vs. pt
    AliRsnListOutput *outMonitorTPCchi2 = new AliRsnListOutput("TPCchi2VsPt", AliRsnListOutput::kHistoDefault);
    if (!useTH1) outMonitorTPCchi2->AddValue(axisMomPt);
    outMonitorTPCchi2->AddValue(axisTPCchi2);
    if (mon) mon->Add(outMonitorTPCchi2);
    if (lm) lm->AddOutput(outMonitorTPCchi2);

    // output: 2D histogram of TPC chi2 vs. TPC momentum
    AliRsnListOutput *outMonitorTPCchi2VsPtpc = new AliRsnListOutput("TPCchi2VsPtpc", AliRsnListOutput::kHistoDefault);
   if (!useTH1) outMonitorTPCchi2VsPtpc->AddValue(axisMomTPC);
    outMonitorTPCchi2VsPtpc->AddValue(axisTPCchi2);
    if (mon) mon->Add(outMonitorTPCchi2VsPtpc);
    if (lm) lm->AddOutput(outMonitorTPCchi2VsPtpc);
  }
  /****************************************************************/
  /***************       MONITOR TPC           ********************/
  /****************************************************************/
  if (monitorTPCpid) {
    // output: 2D histogram of TPC signal vs. TPC momentum
    AliRsnListOutput *outMonitordEdxTPC = new AliRsnListOutput("dEdx_VsPtpc", AliRsnListOutput::kHistoDefault);
    outMonitordEdxTPC->AddValue(axisMomTPC);
    outMonitordEdxTPC->AddValue(axisSigTPC);
    if (mon) mon->Add(outMonitordEdxTPC);
    if (lm) lm->AddOutput(outMonitordEdxTPC);

    // output: 2D histogram of TPC nsigma pi vs. TPC momentum
    AliRsnListOutput *outMonitorTPCnsigmaPi = new AliRsnListOutput("TPC_nsigmaPi_VsPtpc", AliRsnListOutput::kHistoDefault);
    outMonitorTPCnsigmaPi->AddValue(axisMomTPC);
    outMonitorTPCnsigmaPi->AddValue(axisTPCnsigmaPi);
    if (mon) mon->Add(outMonitorTPCnsigmaPi);
    if (lm) lm->AddOutput(outMonitorTPCnsigmaPi);

    // output: 2D histogram of TPC nsigma K vs. TPC momentum
    AliRsnListOutput *outMonitorTPCnsigmaK = new AliRsnListOutput("TPC_nsigmaK_VsPtpc", AliRsnListOutput::kHistoDefault);
    outMonitorTPCnsigmaK->AddValue(axisMomTPC);
    outMonitorTPCnsigmaK->AddValue(axisTPCnsigmaK);
    if (mon) mon->Add(outMonitorTPCnsigmaK);
    if (lm) lm->AddOutput(outMonitorTPCnsigmaK);

    // output: 2D histogram of TPC nsigma pro vs. TPC momentum
    AliRsnListOutput *outMonitorTPCnsigmaP = new AliRsnListOutput("TPC_nsigmaPro_VsPtpc", AliRsnListOutput::kHistoDefault);
    outMonitorTPCnsigmaP->AddValue(axisMomTPC);
    outMonitorTPCnsigmaP->AddValue(axisTPCnsigmaP);
    if (mon) mon->Add(outMonitorTPCnsigmaP);
    if (lm) lm->AddOutput(outMonitorTPCnsigmaP);
  }
  /****************************************************************/
  /***************       MONITOR TOF           ********************/
  /****************************************************************/
  // output:2D histogram of TOF Nsigma pi vs. TPC Nsigma pi vs momentum
  AliRsnListOutput *outMonitorTOFvsTPCnsigmaPi = new AliRsnListOutput("TOFnsigmaPi_TPCnsigmaPi", AliRsnListOutput::kHistoDefault);
  outMonitorTOFvsTPCnsigmaPi->AddValue(axisTOFnsigmaPi);
  outMonitorTOFvsTPCnsigmaPi->AddValue(axisTPCnsigmaPi);
  if (useTH3) outMonitorTOFvsTPCnsigmaPi->AddValue(axisMomP);
  if (mon) mon->Add(outMonitorTOFvsTPCnsigmaPi);
  if (lm) lm->AddOutput(outMonitorTOFvsTPCnsigmaPi);
  
  // output:2D histogram of TOF Nsigma pi vs. TPC Nsigma pi vs momentum
  AliRsnListOutput *outMonitorTOFvsTPCnsigmaK = new AliRsnListOutput("TOFnsigmaK_TPCnsigmaK", AliRsnListOutput::kHistoDefault);
  outMonitorTOFvsTPCnsigmaK->AddValue(axisTOFnsigmaK);
  outMonitorTOFvsTPCnsigmaK->AddValue(axisTPCnsigmaK);
  if (useTH3) outMonitorTOFvsTPCnsigmaK->AddValue(axisMomP);
  if (mon) mon->Add(outMonitorTOFvsTPCnsigmaK);
  if (lm) lm->AddOutput(outMonitorTOFvsTPCnsigmaK);

  AliRsnListOutput *outMonitorTOFvsTPCnsigmaP = new AliRsnListOutput("TOFnsigmaP_TPCnsigmaP", AliRsnListOutput::kHistoDefault);
  outMonitorTOFvsTPCnsigmaP->AddValue(axisTOFnsigmaP);
  outMonitorTOFvsTPCnsigmaP->AddValue(axisTPCnsigmaP);
  if (useTH3) outMonitorTOFvsTPCnsigmaP->AddValue(axisMomP);
  if (mon) mon->Add(outMonitorTOFvsTPCnsigmaP);
  if (lm) lm->AddOutput(outMonitorTOFvsTPCnsigmaP);

  // // output: 2D histogram of TOF signal vs. momentum
  // AliRsnListOutput *outMonitorTimeTOF = new AliRsnListOutput("time_VsP", AliRsnListOutput::kHistoDefault);
  // outMonitorTimeTOF->AddValue(axisMomP);
  // outMonitorTimeTOF->AddValue(axisSigTOF);
  // if (mon) mon->Add(outMonitorTimeTOF);
  // if (lm) lm->AddOutput(outMonitorTimeTOF);

  // // output: 2D histogram of TOF signal vs. pt
  // AliRsnListOutput *outMonitorTimeTOFPt = new AliRsnListOutput("time_VsPt", AliRsnListOutput::kHistoDefault);
  // outMonitorTimeTOFPt->AddValue(axisMomPt);
  // outMonitorTimeTOFPt->AddValue(axisSigTOF);
  // if (mon) mon->Add(outMonitorTimeTOFPt);
  // if (lm) lm->AddOutput(outMonitorTimeTOFPt);

  if (monitorTOFpid) {
    // output: 2D histogram of TOF Nsigma pi vs. TPC momentum
    AliRsnListOutput *outMonitorTOFnsigmaPi = new AliRsnListOutput("TOF_nsigmaPi_vsP", AliRsnListOutput::kHistoDefault);
    outMonitorTOFnsigmaPi->AddValue(axisMomP);
    outMonitorTOFnsigmaPi->AddValue(axisTOFnsigmaPi);
    if (mon) mon->Add(outMonitorTOFnsigmaPi);
    if (lm) lm->AddOutput(outMonitorTOFnsigmaPi);
     
    // output: 2D histogram of TOF signal vs. TOF momentum
    AliRsnListOutput *outMonitorTOFnsigmaK = new AliRsnListOutput("TOF_nsigmaK_vsP", AliRsnListOutput::kHistoDefault);
    outMonitorTOFnsigmaK->AddValue(axisMomP);
    outMonitorTOFnsigmaK->AddValue(axisTOFnsigmaK);
    if (mon) mon->Add(outMonitorTOFnsigmaK);
    if (lm) lm->AddOutput(outMonitorTOFnsigmaK);
      
    // output: 2D histogram of TOF signal vs. TOF momentum
    AliRsnListOutput *outMonitorTOFnsigmaP = new AliRsnListOutput("TOF_nsigmaPro_vsP", AliRsnListOutput::kHistoDefault);
    outMonitorTOFnsigmaP->AddValue(axisMomP);
    outMonitorTOFnsigmaP->AddValue(axisTOFnsigmaP);
    if (mon) mon->Add(outMonitorTOFnsigmaP);
    if (lm) lm->AddOutput(outMonitorTOFnsigmaP);

   // output: 2D histogram of TOF Delta pi vs. TPC momentum
    AliRsnListOutput *outMonitorTOFdeltaPi = new AliRsnListOutput("TOF_deltaPi_vsP", AliRsnListOutput::kHistoDefault);
    outMonitorTOFdeltaPi->AddValue(axisMomP);
    outMonitorTOFdeltaPi->AddValue(axisTOFdeltaPi);
    if (mon) mon->Add(outMonitorTOFdeltaPi);
    if (lm) lm->AddOutput(outMonitorTOFdeltaPi);
     
    // output: 2D histogram of TOF signal vs. TOF momentum
    AliRsnListOutput *outMonitorTOFdeltaK = new AliRsnListOutput("TOF_deltaK_vsP", AliRsnListOutput::kHistoDefault);
    outMonitorTOFdeltaK->AddValue(axisMomP);
    outMonitorTOFdeltaK->AddValue(axisTOFdeltaK);
    if (mon) mon->Add(outMonitorTOFdeltaK);
    if (lm) lm->AddOutput(outMonitorTOFdeltaK);
      
    // output: 2D histogram of TOF signal vs. TOF momentum
    AliRsnListOutput *outMonitorTOFdeltaP = new AliRsnListOutput("TOF_deltaPro_vsP", AliRsnListOutput::kHistoDefault);
    outMonitorTOFdeltaP->AddValue(axisMomP);
    outMonitorTOFdeltaP->AddValue(axisTOFdeltaP);
    if (mon) mon->Add(outMonitorTOFdeltaP);
    if (lm) lm->AddOutput(outMonitorTOFdeltaP);
  }

  /****************************************************************/
  /***************       MONITOR MC            ********************/
  /****************************************************************/
  if (useMCMon) {
    //Momentum 
    AliRsnValueDaughter *axisMomPMC = new AliRsnValueDaughter("pMC", AliRsnValueDaughter::kP);
    axisMomPMC->SetUseMCInfo(kTRUE);
    axisMomPMC->SetBins(0.0,10.0,0.01);
    // Momentum Pt
    AliRsnValueDaughter *axisMomPtMC = new AliRsnValueDaughter("ptMC", AliRsnValueDaughter::kPt);
    axisMomPtMC->SetUseMCInfo(kTRUE);
    axisMomPtMC->SetBins(0.0,10.0,0.01);
    // Eta
    AliRsnValueDaughter *axisEtaMC = new AliRsnValueDaughter("etaMC", AliRsnValueDaughter::kEta);
    axisEtaMC->SetUseMCInfo(kTRUE);
    axisEtaMC->SetBins(-1.0,1.0,0.1);

    // output: 2D histo for kine acceptance
    AliRsnListOutput *outMonitorPtVsEtaMC = new AliRsnListOutput("Pt_VsEtaMC", AliRsnListOutput::kHistoDefault);
    outMonitorPtVsEtaMC->AddValue(axisMomPtMC);
    outMonitorPtVsEtaMC->AddValue(axisEtaMC);
    if (mon) mon->Add(outMonitorPtVsEtaMC);
    if (lm) lm->AddOutput(outMonitorPtVsEtaMC);

    // output: 1D histo pt from MC
    AliRsnListOutput *outMonitorPtMC = new AliRsnListOutput("PtMC", AliRsnListOutput::kHistoDefault);
    outMonitorPtMC->AddValue(axisMomPtMC);
    if (mon) mon->Add(outMonitorPtMC);
    if (lm) lm->AddOutput(outMonitorPtMC);
 
    // output: 1D histo eta from MC
    AliRsnListOutput *outMonitorEtaMC = new AliRsnListOutput("EtaMC", AliRsnListOutput::kHistoDefault);
    outMonitorEtaMC->AddValue(axisEtaMC);
    if (mon) mon->Add(outMonitorEtaMC);
    if (lm) lm->AddOutput(outMonitorEtaMC);
  }
 
}

#endif
