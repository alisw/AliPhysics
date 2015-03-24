void MakeOCDBConfigPreprocessor()
{
  // create the object and set some values:
  AliCDBEntry *fCDBEntry = new AliCDBEntry();
  TEnv *fConfEnv = new TEnv();

  fConfEnv->SetValue("Pedestal","DAQ");
  fConfEnv->SetValue("Signal","DAQ");
  fConfEnv->SetValue("Temperature","ON");
  fConfEnv->SetValue("ErrorHandling","ON");

  fCDBEntry->SetObject(fConfEnv);

  // done; now save it..; add some metadata business etc.
  Int_t firstRun   =  0;
  Int_t lastRun    =  999999999;
  Int_t version = 0;
  Int_t subversion = 0;
  Int_t beamPeriod =  1;

  char filename[200];
  sprintf(filename, "Preprocessor/Run%d_%d_v%d_s%d.root",
	  firstRun, lastRun, version, subversion);

  AliCDBMetaData md;
  md.SetBeamPeriod(beamPeriod);
  md.SetResponsible("David Silvermyr");
  
  AliCDBId id("EMCAL/Config/Preprocessor", firstRun, lastRun, version, subversion);

  fCDBEntry->SetId(id);
  fCDBEntry->SetMetaData(&md);

  // ok, write the file
  TFile f(filename, "recreate");
  if (!f.IsZombie()) {
    f.cd();
    fCDBEntry->Write("AliCDBEntry");
    f.Close();
  }

}
