/*

  Simple macro to read info from OADB file TOFPIDParams.root

  @P. Antonioli / INFN - BO

*/
readTOFparams(Int_t run, char* pass,char *oadbfile="$ALICE_PHYSICS/OADB/COMMON/PID/data/TOFPIDParams.root")

{
  TString fname = Form("%s",oadbfile);

  Printf("-----------------------------------------------------------------------------------------------");
  Printf(" Reading OADB TOF Response params from: %s",fname.Data());
  
  TString passName = Form("%s",pass);
  TString stMethod[4]={"kFILL_T0","kTOF_T0","kT0_T0","kBest_T0"};
  TFile *oadbf = new TFile(fname);
  AliOADBContainer *oadbc = oadbf->Get("TOFoadb");
  AliTOFPIDParams *tofPar = dynamic_cast<AliTOFPIDParams*>oadbc->GetObject(run,"TOFparams",passName);
  oadbf->Close();
  delete oadbf;
  delete oadbc;

  if (tofPar) {
    Printf("-----------------------------------------------------------------------------------------------");
    Printf(" TOF PID params for run # %d [entry tag: %s]",run,tofPar->GetOADBentryTag());
    Printf("   Intrinsic resolution (MRPC+electronics+clock+calibration): %6.2f ps",tofPar->GetTOFresolution());
    Printf("   Start Time method: %s",stMethod[tofPar->GetStartTimeMethod()].Data());
    for (Int_t i=0;i<4;i++) Printf("   TOF PID (p) params (%d): %6.4f",i,tofPar->GetSigParams(i));
    Printf("   Fraction of tracks within gaussian behaviour: %6.4f",tofPar->GetTOFtail());
    Printf("   MC: Fraction of tracks (percentage) to cut to fit matching in data: %6.2f%%",tofPar->GetTOFmatchingLossMC());
    Printf("   MC: Fraction of tracks (percentage) to add to fit mismatch in data: %6.2f%%",tofPar->GetTOFadditionalMismForMC());
    Printf("   Time Offset (intercalibration adjust with T0): %6.2f ps",tofPar->GetTOFtimeOffset());
  }


}
