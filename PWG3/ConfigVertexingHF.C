AliAnalysisVertexingHF* ConfigVertexingHF() {

  printf("Call to AliAnalysisVertexingHF parameters setting :\n");
  vHF = new AliAnalysisVertexingHF();
 
  //vHF->SetJPSItoEleOff();
  vHF->Set3ProngOff();
  vHF->Set4ProngOff();
  vHF->SetITSrefitRequired();
  vHF->SetBothSPDNotRequired();
  vHF->SetMinITSCls(5);
  vHF->SetMinPtCut(0.5);
  vHF->SetMind0Cut(0.);
  vHF->SetD0toKpiCuts(0.2,999999.,1.1,0.,0.,999999.,999999.,999999.,0.3);
  vHF->SetBtoJPSICuts(0.350);
  vHF->SetDplusCuts(0.2,0.,0.,0.,0.,0.01,0.06,0.,0.,0.8);
 
  return vHF;
}



