/***************************************************************************
              fbellini@cern.ch - last modified on 09/04/2013

 Macro to add monitoring histograms for track and PID cuts to the rsn task 
 Tuned for monitoring TOF KStar analysis of PbPb 2010 data

Options ("opt" argument):
- NoTOFSIGMA  --> disable the nsigma_TOF vs p histos
- NoTPCSIGMA  --> disable the nsigma_TPC vs p_TPC histos
- NoSIGN --> disable splitting for single track charge
- NoTrackQ  --> disable track quality monitoring (DCA, nCls) histos
/***************************************************************************/

void AddMonitorOutput(Bool_t useMCMon = 0, TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *lm=0)
{
  // Multiplicity/centrality
  AliRsnValueEvent *multi = new AliRsnValueEvent("multi",AliRsnValueEvent::kMult);
  multi->SetBins(0.0, 400.0, 1);
  // Momentum
  AliRsnValueDaughter *axisMomTPC = new AliRsnValueDaughter("pTPC", AliRsnValueDaughter::kPtpc);
  axisMomTPC->SetBins(0.0, 10.0, 0.01);
  AliRsnValueDaughter *axisMomP = new AliRsnValueDaughter("p", AliRsnValueDaughter::kP);
  axisMomP->SetBins(0.0, 10.0, 0.01);
  // Momentum Pt
  AliRsnValueDaughter *axisMomPt = new AliRsnValueDaughter("pt", AliRsnValueDaughter::kPt);
  axisMomPt->SetBins(0.0,10.0,0.01);
  // Eta
  AliRsnValueDaughter *axisMomEta = new AliRsnValueDaughter("eta", AliRsnValueDaughter::kEta);
  axisMomEta->SetBins(-1.0, 1.0, 0.1);
  //ITS clusters
  AliRsnValueDaughter *axisITScls = new AliRsnValueDaughter("ITScls", AliRsnValueDaughter::kNITSclusters);
  axisITScls->SetBins(0.0, 10.0, 1.0);
  //TPC clusters
  AliRsnValueDaughter *axisTPCcls = new AliRsnValueDaughter("TPCcls", AliRsnValueDaughter::kNTPCclusters);
  axisTPCcls->SetBins(0.0, 300.0, 1.0);
  //ITS chi2
  AliRsnValueDaughter *axisITSchi2 = new AliRsnValueDaughter("ITSchi2", AliRsnValueDaughter::kITSchi2);
  axisITSchi2->SetBins(0.0, 10.0, 0.1);
  //TPC chi2
  AliRsnValueDaughter *axisTPCchi2 = new AliRsnValueDaughter("TPCchi2", AliRsnValueDaughter::kTPCchi2);
  axisTPCchi2->SetBins(0.0, 10.0, 0.1);
  //DCA xy
  AliRsnValueDaughter *axisDCAxy = new AliRsnValueDaughter("DCAxy", AliRsnValueDaughter::kDCAXY);
  axisDCAxy->SetBins(-2.0, 2.0, 0.1);
  //DCA z
  AliRsnValueDaughter *axisDCAz = new AliRsnValueDaughter("DCAz", AliRsnValueDaughter::kDCAZ);
  axisDCAz->SetBins(-10.0, 10.0, 0.1);
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
  if (!opt.Contains("NoSIGN")) {
    outMonitorP->AddValue(axisCharge);
  }
  if (mon) mon->Add(outMonitorP);
  if (lm) lm->AddOutput(outMonitorP);
  
  // output:  TH1D for pt
  AliRsnListOutput *outMonitorPt = new AliRsnListOutput("Pt", AliRsnListOutput::kHistoDefault);
  outMonitorPt->AddValue(axisMomPt);
  if (!opt.Contains("NoSIGN")) {
    outMonitorPt->AddValue(axisCharge);
  }
  if (mon) mon->Add(outMonitorPt);
  if (lm) lm->AddOutput(outMonitorPt);
  
  // output: TH1D for pseudorapidity
  AliRsnListOutput *outMonitorEta = new AliRsnListOutput("Eta", AliRsnListOutput::kHistoDefault);
  outMonitorEta->AddValue(axisMomEta);
  if (!opt.Contains("NoSIGN")) {
    outMonitorEta->AddValue(axisCharge);
  }
  if (mon) mon->Add(outMonitorEta);
  if (lm) lm->AddOutput(outMonitorEta);
  
  // output:  TH2D for phi at vertex
  AliRsnListOutput *outMonitorPhi = new AliRsnListOutput("Phi", AliRsnListOutput::kHistoDefault);
  //outMonitorPhi->AddValue(axisMomPt);
  outMonitorPhi->AddValue(axisPhi);
  if (!opt.Contains("NoSIGN")) {
    outMonitorPhi->AddValue(axisCharge);
  }
  if (mon) mon->Add(outMonitorPhi);
  if (lm) lm->AddOutput(outMonitorPhi);
  
  // output:  TH2D for phiOuterTPC at TPC outer radius
  AliRsnListOutput *outMonitorPhiOuterTPC = new AliRsnListOutput("PhiOuterTPC", AliRsnListOutput::kHistoDefault);
  //outMonitorPhiOuterTPC->AddValue(axisMomPt);
  outMonitorPhiOuterTPC->AddValue(axisPhiOuterTPC);
  if (!opt.Contains("NoSIGN")) {
    outMonitorPhiOuterTPC->AddValue(axisCharge);
  }
  if (mon) mon->Add(outMonitorPhiOuterTPC);
  if (lm) lm->AddOutput(outMonitorPhiOuterTPC);

  // output:  TH2D for phi vs pt
  AliRsnListOutput *outMonitorPhiVsPt = new AliRsnListOutput("PhiVsPt", AliRsnListOutput::kHistoDefault);
  outMonitorPhiVsPt->AddValue(axisMomPt);
  outMonitorPhiVsPt->AddValue(axisPhi);
  if (!opt.Contains("NoSIGN")) {
    outMonitorPhiVsPt->AddValue(axisCharge);
  }
  if (mon) mon->Add(outMonitorPhiVsPt);
  if (lm) lm->AddOutput(outMonitorPhiVsPt);

  /****************************************************************/
  /***************      MONITOR TRACK QUALITY  ********************/
  /****************************************************************/
  if (!opt.Contains("NoTrackQ")) {
    // output: 2D histogram of DCAxy vs pt
    AliRsnListOutput *outMonitorDCAxy = new AliRsnListOutput("DCAxyVsPt", AliRsnListOutput::kHistoDefault);
    outMonitorDCAxy->AddValue(axisMomPt);
    outMonitorDCAxy->AddValue(axisDCAxy);
    if (mon) mon->Add(outMonitorDCAxy);
    if (lm) lm->AddOutput(outMonitorDCAxy);

    // output: 2D histogram of DCAz vs P
    AliRsnListOutput *outMonitorDCAz = new AliRsnListOutput("DCAzVsP", AliRsnListOutput::kHistoDefault);
    outMonitorDCAz->AddValue(axisMomP);
    outMonitorDCAz->AddValue(axisDCAz);
    if (mon) mon->Add(outMonitorDCAz);
    if (lm) lm->AddOutput(outMonitorDCAz);

    // output: 2D histogram of ITS cls vs pt
    AliRsnListOutput *outMonitorITScls = new AliRsnListOutput("ITSclsVsPt", AliRsnListOutput::kHistoDefault);
    outMonitorITScls->AddValue(axisMomPt);
    outMonitorITScls->AddValue(axisITScls);
    if (mon) mon->Add(outMonitorITScls);
    if (lm) lm->AddOutput(outMonitorITScls);

    // output: 2D histogram of TPC cls vs. pt
    AliRsnListOutput *outMonitorTPCcls = new AliRsnListOutput("TPCclsVsPt", AliRsnListOutput::kHistoDefault);
    outMonitorTPCcls->AddValue(axisMomPt);
    outMonitorTPCcls->AddValue(axisTPCcls);
    if (mon) mon->Add(outMonitorTPCcls);
    if (lm) lm->AddOutput(outMonitorTPCcls);

    // output: 2D histogram of TPC cls vs. TPC momentum
    AliRsnListOutput *outMonitorTPCclsVsPtpc = new AliRsnListOutput("TPCclsVsPtpc", AliRsnListOutput::kHistoDefault);
    outMonitorTPCclsVsPtpc->AddValue(axisMomTPC);
    outMonitorTPCclsVsPtpc->AddValue(axisTPCcls);
    if (mon) mon->Add(outMonitorTPCclsVsPtpc);
    if (lm) lm->AddOutput(outMonitorTPCclsVsPtpc);

    // output: 2D histogram of ITS chi2 vs pt
    AliRsnListOutput *outMonitorITSchi2 = new AliRsnListOutput("ITSchi2VsPt", AliRsnListOutput::kHistoDefault);
    outMonitorITSchi2->AddValue(axisMomPt);
    outMonitorITSchi2->AddValue(axisITSchi2);
    if (mon) mon->Add(outMonitorITSchi2);
    if (lm) lm->AddOutput(outMonitorITSchi2);

    // output: 2D histogram of TPC chi2 vs. pt
    AliRsnListOutput *outMonitorTPCchi2 = new AliRsnListOutput("TPCchi2VsPt", AliRsnListOutput::kHistoDefault);
    outMonitorTPCchi2->AddValue(axisMomPt);
    outMonitorTPCchi2->AddValue(axisTPCchi2);
    if (mon) mon->Add(outMonitorTPCchi2);
    if (lm) lm->AddOutput(outMonitorTPCchi2);

    // output: 2D histogram of TPC chi2 vs. TPC momentum
    AliRsnListOutput *outMonitorTPCchi2VsPtpc = new AliRsnListOutput("TPCchi2VsPtpc", AliRsnListOutput::kHistoDefault);
    outMonitorTPCchi2VsPtpc->AddValue(axisMomTPC);
    outMonitorTPCchi2VsPtpc->AddValue(axisTPCchi2);
    if (mon) mon->Add(outMonitorTPCchi2VsPtpc);
    if (lm) lm->AddOutput(outMonitorTPCchi2VsPtpc);
  }
  /****************************************************************/
  /***************       MONITOR TPC           ********************/
  /****************************************************************/
  if (!opt.Contains("NoTPCSIGMA")) {
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

  if (!opt.Contains("NoTOFSIGMA")) {
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

    // output: 1D histo pt from MC
    AliRsnListOutput *outMonitorPtMC = new AliRsnListOutput("PtMC", AliRsnListOutput::kHistoDefault);
    outMonitorPtMC->AddValue(axisMomPtMC);
    if (mon) mon->Add(outMonitorPtMC);
    if (lm) lm->AddOutput(outMonitorPtMC);
  }
 
}
