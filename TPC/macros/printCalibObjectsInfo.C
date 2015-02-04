/// \file printCalibObjectsInfo.C

printCalibObjectsInfo(const char* filename="CalibObjects.root")
{
  gROOT->Macro("$ALICE_ROOT/PWGPP/CalibMacros/CPass0/LoadLibraries.C");
  TFile f(filename,"read");
 
  TDirectory* tpcCalib = dynamic_cast<TDirectory*>(f.Get("TPCCalib"));
  tpcCalib->ls();

  AliTPCcalibCalib* calibCalib = tpcCalib->Get("calibTPC");
  AliTPCcalibTimeGain* calibTimeGain = tpcCalib->Get("calibTimeGain");
  AliTPCcalibGainMult* calibGainMult = tpcCalib->Get("calibGainMult");
  AliTPCcalibTime* calibTime = tpcCalib->Get("calibTime");
  
  
  if (!calibCalib || !calibTimeGain || !calibGainMult || !calibTime)
  {
    printf("file empty\n");
    return;
  }

  printf("\ncalibTimeGain->GetHistGainTime()->GetEntries() = %10i, size: %.2f MB\n", calibTimeGain->GetHistGainTime()->GetEntries(), 
      (AliSysInfo::EstimateObjectSize(calibTimeGain->GetHistGainTime()))/1024./1024);
  printf("\n");
  
  printf("calibGainMult->GetHistGainSector()->GetEntries() = %10i, size: %.2f MB\n", calibGainMult->GetHistGainSector()->GetEntries(), 
      (AliSysInfo::EstimateObjectSize(calibGainMult->GetHistGainSector())/1024./1024));
  printf("calibGainMult->GetHistPadEqual()->GetEntries()   = %10i, size: %.2f MB\n", calibGainMult->GetHistPadEqual()->GetEntries(), 
      (AliSysInfo::EstimateObjectSize(calibGainMult->GetHistPadEqual())/1024./1024));
  printf("calibGainMult->GetHistGainMult()->GetEntries()   = %10i, size: %.2f MB\n", calibGainMult->GetHistGainMult()->GetEntries(),
      (AliSysInfo::EstimateObjectSize(calibGainMult->GetHistGainMult())/1024./1024));
  printf("calibGainMult->GetHistdEdxMap()->GetEntries()    = %10i, size: %.2f MB\n", calibGainMult->GetHistdEdxMap()->GetEntries(),
      (AliSysInfo::EstimateObjectSize(calibGainMult->GetHistdEdxMap())/1024./1024));
  printf("calibGainMult->GetHistdEdxMax()->GetEntries()    = %10i, size: %.2f MB\n", calibGainMult->GetHistdEdxMax()->GetEntries(),
      (AliSysInfo::EstimateObjectSize(calibGainMult->GetHistdEdxMax())/1024./1024));
  printf("calibGainMult->GetHistdEdxTot()->GetEntries()    = %10i, size: %.2f MB\n", calibGainMult->GetHistdEdxTot()->GetEntries(),
      (AliSysInfo::EstimateObjectSize(calibGainMult->GetHistdEdxTot())/1024./1024));
  printf("\n");

  for (int n=0; n<3; n++)
  {
    printf("calibTime->GetHistVdriftLaserA(%i)->GetEntries()          = %10i, size: %.2f MB\n", n, calibTime->GetHistVdriftLaserA(n)->GetEntries(),
        (AliSysInfo::EstimateObjectSize(calibTime->GetHistVdriftLaserA(n))/1024./1024));
    printf("calibTime->GetHistVdriftLaserC(%i)->GetEntries()          = %10i, size: %.2f MB\n", n, calibTime->GetHistVdriftLaserC(n)->GetEntries(),
        (AliSysInfo::EstimateObjectSize(calibTime->GetHistVdriftLaserC(n))/1024./1024));
  }

  for (int n=0; n<12; n++)
  {
    printf("calibTime->GetTPCVertexHisto(%i)->GetEntries()            = %10i, size: %.2f MB\n", n, calibTime->GetTPCVertexHisto(n)->GetEntries(),
        (AliSysInfo::EstimateObjectSize(calibTime->GetTPCVertexHisto(n))/1024./1024));
  }

  for (int n=0; n<5; n++)
  {
    printf("calibTime->GetTPCVertexHistoCorrelation(%i)->GetEntries() = %10i, size: %.2f MB\n", n, calibTime->GetTPCVertexHistoCorrelation(n)->GetEntries(),
        (AliSysInfo::EstimateObjectSize(calibTime->GetTPCVertexHistoCorrelation(n))/1024./1024));
    printf("calibTime->GetResHistoTPCCE(%i)->GetEntries()             = %10i, size: %.2f MB\n", n, calibTime->GetResHistoTPCCE(n)->GetEntries(),
        (AliSysInfo::EstimateObjectSize(calibTime->GetResHistoTPCCE(n))/1024./1024));
    printf("calibTime->GetResHistoTPCITS(%i)->GetEntries()            = %10i, size: %.2f MB\n", n, calibTime->GetResHistoTPCITS(n)->GetEntries(),
        (AliSysInfo::EstimateObjectSize(calibTime->GetResHistoTPCITS(n))/1024./1024));
    printf("calibTime->GetResHistoTPCvertex(%i)->GetEntries()         = %10i, size: %.2f MB\n", n, calibTime->GetResHistoTPCvertex(n)->GetEntries(),
        (AliSysInfo::EstimateObjectSize(calibTime->GetResHistoTPCvertex(n))/1024./1024));
    printf("calibTime->GetResHistoTPCTRD(%i)->GetEntries()            = %10i, size: %.2f MB\n", n, calibTime->GetResHistoTPCTRD(n)->GetEntries(),
        (AliSysInfo::EstimateObjectSize(calibTime->GetResHistoTPCTRD(n))/1024./1024));
    printf("calibTime->GetResHistoTPCTOF(%i)->GetEntries()            = %10i, size: %.2f MB\n", n, calibTime->GetResHistoTPCTOF(n)->GetEntries(),
        (AliSysInfo::EstimateObjectSize(calibTime->GetResHistoTPCTOF(n))/1024./1024));
  }

  TObjArray* TPCCluster = f.Get("TPCCluster");
  if (TPCCluster)
  {
    TPCCluster->Print();
    AliTPCcalibTracks* calibTracks = TPCCluster->FindObject("calibTracks");
    printf("calibTracks->fHisDeltaY->GetEntries() = %10i, size: %.2f MB\n", n, (calibTracks->fHisDeltaY)->GetEntries(),
        (AliSysInfo::EstimateObjectSize(calibTracks->fHisDeltaY)/1024./1024.));
    printf("calibTracks->fHisDeltaZ->GetEntries() = %10i, size: %.2f MB\n", n, (calibTracks->fHisDeltaZ)->GetEntries(),
        (AliSysInfo::EstimateObjectSize(calibTracks->fHisDeltaZ)/1024./1024.));
    printf("calibTracks->fHisRMSY->GetEntries() = %10i, size: %.2f MB\n", n, (calibTracks->fHisRMSY)->GetEntries(),
        (AliSysInfo::EstimateObjectSize(calibTracks->fHisRMSY)/1024./1024.));
    printf("calibTracks->fHisRMSZ->GetEntries() = %10i, size: %.2f MB\n", n, (calibTracks->fHisRMSZ)->GetEntries(),
        (AliSysInfo::EstimateObjectSize(calibTracks->fHisRMSZ)/1024./1024.));
    printf("calibTracks->fHisQmax->GetEntries() = %10i, size: %.2f MB\n", n, (calibTracks->fHisQmax)->GetEntries(),
        (AliSysInfo::EstimateObjectSize(calibTracks->fHisQmax)/1024./1024.));
    printf("calibTracks->fHisQtot->GetEntries() = %10i, size: %.2f MB\n", n, (calibTracks->fHisQtot)->GetEntries(),
        (AliSysInfo::EstimateObjectSize(calibTracks->fHisQtot)/1024./1024.));
  }


}
