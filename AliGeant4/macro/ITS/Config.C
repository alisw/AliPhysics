void Config(Int_t version)
{
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libITS");

  AliITS* ITS = 0;
  switch (version) {
    case 1: ITS  = new AliITSv1("ITS","ITSv1 detector"); break;
    case 3: ITS  = new AliITSv3("ITS","ITSv3 detector"); break; 
    case 5: ITS  = new AliITSv5("ITS","normal ITS");       break;
  }  

//=================== ITS parameters ============================
//
// EUCLID is a flag to output (=1) both geometry and media to two ASCII files 
// (called by default ITSgeometry.euc and ITSgeometry.tme) in a format
// understandable to the CAD system EUCLID. The default (=0) means that you 
// dont want to use this facility.
//
ITS->SetEUCLID(0);
}
