void Config(Int_t version)
{
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libTPC");
  gROOT->LoadMacro("SetTPCParam.C");
  AliTPCParam *param = SetTPCParam();

  AliTPC* TPC = 0;
  switch (version) {
    case 0: TPC  = new AliTPCv0("TPC","TPCv0 detector"); break;
    case 1: TPC  = new AliTPCv1("TPC","normal TPC");     break;
    case 2: TPC  = new AliTPCv2("TPC","TPCv2 detector"); break;
    case 3: TPC  = new AliTPCv3("TPC","TPCv3 detector"); break;
  }   

//============================ TPC parameters ================================
// --- This allows the user to specify sectors for the SLOW (TPC geometry 2)
// --- Simulator. SecAL (SecAU) <0 means that ALL lower (upper)
// --- sectors are specified, any value other than that requires at least one 
// --- sector (lower or upper)to be specified!
// --- Reminder: sectors 1-24 are lower sectors (1-12 -> z>0, 13-24 -> z<0)
// ---           sectors 25-72 are the upper ones (25-48 -> z>0, 49-72 -> z<0)
// --- SecLows - number of lower sectors specified (up to 6)
// --- SecUps - number of upper sectors specified (up to 12)
// --- Sens - sensitive strips for the Slow Simulator !!!
// --- This does NOT work if all S or L-sectors are specified, i.e.
// --- if SecAL or SecAU < 0
//
//
//-----------------------------------------------------------------------------

  TPC->SetParam(param); // pass the parameter object to the TPC

// set gas mixture

TPC->SetGasMixt(2,20,10,-1,0.9,0.1,0.);
TPC->SetSecAL(4);
TPC->SetSecAU(4);
TPC->SetSecLows(1,  2,  3, 19, 20, 21);
TPC->SetSecUps(37, 38, 39, 37+18, 38+18, 39+18, -1, -1, -1, -1, -1, -1);
TPC->SetSens(1);

if (TPC->IsVersion()==1) param->Write(param->GetTitle());
}
