void Config(Int_t version)
{
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libITS");

  AliITS* ITS = 0;
  switch (version) {
    case 1: ITS  = new AliITSv1("ITS","Old ITS coarse version as of the ALICE TP"); break;
    case 3: ITS  = new AliITSv3("ITS","Old ITS detailed version as of the ALICE TP"); break; 
    case 5: ITS  = new AliITSv5("ITS","Current ITS detailed version used for the ITS TDR"); break;
    case 6: ITS  = new AliITSv5symm( "ITS","Updated ITS TDR detailed version with symmetric services");break;
    case 7: ITS  = new AliITSv5asymm("ITS","Updates ITS TDR detailed version with asymmetric services");break;
    case 8: ITS  = new AliITSvPPRcoarseasymm("ITS","New ITS coarse version with asymmetric services"); break;
    case 9: ITS  = new AliITSvPPRcoarsesymm( "ITS","New ITS coarse version with symmetric services"); break;
    // 
  }  

    //=================== ITS parameters ============================
    //
    // As the innermost detector in ALICE, the Inner Tracking System "impacts" on
    // almost all other detectors. This involves the fact that the ITS geometry 
    // still has several options to be followed in parallel in order to determine 
    // the best set-up which minimizes the induced background. All the geometries
    // available to date are described in the following. Read carefully the comments 
    // and use the default version (the only one uncommented) unless you are making
    // comparisons and you know what you are doing. In this case just uncomment the
    // ITS geometry you want to use and run Aliroot. 
    //
    // Detailed geometries:
    // ====================
    //AliITS *ITS  = new AliITSv3("ITS","Old ITS detailed version as of the ALICE TP");
    //AliITS *ITS  = new AliITSv5("ITS","Current ITS detailed version used for the ITS TDR");
    //AliITS *ITS  = new AliITSv5symm("ITS","Updated ITS TDR detailed version with symmetric services");
    //AliITS *ITS  = new AliITSv5asymm("ITS","Updates ITS TDR detailed version with asymmetric services");
    //
    //
    // Coarse geometries (warning: no hits are produced with these coarse geometries and they unuseful for reconstruction !):
    // ======================================================================================================================
    //AliITS *ITS  = new AliITSv1("ITS","Old ITS coarse version as of the ALICE TP");
    //AliITS *ITS  = new AliITSvPPRcoarseasymm("ITS","New ITS coarse version with asymmetric services");
    //AliITS *ITS  = new AliITSvPPRcoarsesymm("ITS","New ITS coarse version with symmetric services");
    //
    //
    // Geant3 <-> EUCLID conversion
    // ============================
    //
    // SetEUCLID is a flag to output (=1) or not to output (=0) both geometry and 
    // media to two ASCII files (called by default ITSgeometry.euc and 
    // ITSgeometry.tme) in a format understandable to the CAD system EUCLID. 
    // The default (=0) means that you dont want to use this facility.
    //
    ITS->SetEUCLID(0);
}
