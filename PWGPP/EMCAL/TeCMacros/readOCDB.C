#include "readOCDB_Temperature.C"
#include "readOCDB_LED.C"

Int_t readTemp(Int_t *runs, Int_t nruns=1, const char *pername="unspec") 
{
  Int_t ret = 0;

  TObjArray arr;
  arr.SetOwner(1);

  for (Int_t i=0;i<nruns;++i) {
    Int_t rn = runs[i];
    cout << "*** Working on " <<  i+1 << "/" << nruns << " with run number " << rn << " ***" << endl;
    TInfo *ti = readOCDB_Temperature(rn,0);
    if (!ti) 
      continue;
    ti->Print();
    arr.Add(ti);
    ++ret;
  }

  TFile *outf = TFile::Open("tempinfo.root","update");
  arr.Write(Form("temperatures_%s",pername),TObject::kSingleKey);
  outf->ls();
  outf->Close();
  delete outf;
  Double_t frac = 100.*ret/nruns;
  cout << "Finished temperature objects with " << frac << " percent!" << endl;
  return ret;
}

Int_t readLed(Int_t *runs, Int_t nruns=1, const char *pername="unspec") 
{
  Int_t ret = 0;

  TObjArray arr;
  arr.SetOwner(1);

  for (Int_t i=0;i<nruns;++i) {
    Int_t rn = runs[i];
    cout << "*** Working on " <<  i+1 << "/" << nruns << " with run number " << rn << " ***" << endl;
    LInfo *ti = readOCDB_LED(rn,0);
    if (!ti) 
      continue;
    ti->Print();
    arr.Add(ti);
    ++ret;
  }

  TFile *outf = TFile::Open("ledinfo.root","update");
  arr.Write(Form("led_%s",pername),TObject::kSingleKey);
  outf->ls();
  outf->Close();
  delete outf;
  Double_t frac = 100.*ret/nruns;
  cout << "Finished LED objects with " << frac << " percent!" << endl;
  return ret;
}
 

void read_LHC18d(Bool_t test=0) 
{
  Int_t runs[] = {285978,285979,285980,286014,286018,286025,286026,286027,286030,286064,286124,286127,286129,286130,286154,286157,286159,286198,286201,286202,286203,286229,286230,286231,286254,286255,286256,286257,286258,286261,286263,286282,286284,286287,286288,286289,286308,286309,286310,286311,286312,286313,286314,286336,286337,286340,286341,286345,286348,286349,286350};
  
  Int_t nruns = sizeof(runs)/sizeof(Int_t);
  if (test)
    nruns=3;

  //readTemp(runs,nruns,"lhc18d");
  readLed(runs,nruns,"lhc18d");
}
