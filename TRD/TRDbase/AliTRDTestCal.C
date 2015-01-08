void Dump(Char_t* comment, Float_t* f)
{
  cout << f[0] << " " << f[1] << " " << f[2] << endl;
}

void Dump(Char_t* comment, Float_t f)
{
  cout << comment << ": " << f << endl;
}

void AliTRDTestCal()
{
  AliTRDcalibDB* calib = AliTRDcalibDB::Instance();
  if (!calib)
  {
    cerr << "calibDB singleton has already been terminated." << endl;
    return;
  }

  Float_t f[3];

  calib->GetChamberPos(0, f); Dump("chamberpos", f);
  calib->GetChamberRot(0, f); Dump("chamberrot", f);

  calib->GetStackPos(0, 0, f); Dump("stackpos", f);
  calib->GetStackRot(0, 0, f); Dump("stackrot", f);

  calib->GetSuperModulePos(0, f); Dump("smpos", f);
  calib->GetSuperModuleRot(0, f); Dump("smrot", f);

  Dump("vdrift", calib->GetVdrift(0, 0, 0));
  Dump("vdrift-av", calib->GetVdriftAverage(0));

  Dump("t0", calib->GetT0(0, 0, 0));
  Dump("t0-av", calib->GetT0Average(0));

  Dump("gain", calib->GetGainFactor(0, 0, 0));
  Dump("gain-av", calib->GetGainFactorAverage(0));

  Dump("prf", calib->GetPRFWidth(0, 0, 0));

  Dump("sf", calib->GetSamplingFrequency());
  Dump("timebins", calib->GetNumberOfTimeBins());

  Dump("padstatus", calib->GetPadStatus(0, 0, 0));
  Dump("mcmstatus", calib->GetMCMStatus(0, 0, 0));
  Dump("chamberstatus", calib->GetChamberStatus(0));
  Dump("smstatus", calib->GetSuperModuleStatus(0));

  AliTRDCalMonitoring* mon = calib->GetMonitoringObject();
  Dump("monitoring", (Float_t) mon);
  AliTRDCalPIDLQ* pid = calib->GetPIDLQObject();
  Dump("pid", (Float_t) pid);

}