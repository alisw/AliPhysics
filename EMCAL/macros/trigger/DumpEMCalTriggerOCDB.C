// Simple macro to access the information of the trigger
// contained in the OCDB.
// Input is the patch to an OCDB file in alien.

void DumpEMCalTriggerOCDB(const char *ocdb_file = "")
{
  TGrid::Connect("alien://");
  f = TFile::Open(Form("alien://%s",ocdb_file));
  e = (AliCDBEntry*)f->Get("AliCDBEntry");
  d = (AliEMCALTriggerDCSConfig*)e->GetObject();
  c = (AliEMCALTriggerSTUDCSConfig*)d->GetSTUDCSConfig();

  cout << "L1 fw version: " << c->GetFw() << endl;
  cout << "L1-jet patch size: " << 2 + (c->GetFw() >> 16) << " sub-regions (4x4 FOR)" << endl;

}

