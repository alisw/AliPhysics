void Config(Int_t version)
{
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libITS");

  AliITS* ITS = 0;
  switch (version) {
    case 1: ITS  = new AliITSvPPRcoarseasymm("ITS","New ITS coarse version with asymmetric services"); break;
    case 2: ITS  = new AliITSvPPRcoarsesymm("ITS","New ITS coarse version with symmetric services"); break;
    case 3: ITS  = new AliITSvPPRasymm("ITS","New ITS PPR detailed version with asymmetric services"); break;
    case 4: ITS  = new AliITSvPPRsymm("ITS","New ITS PPR detailed version with symmetric services"); break;
    case 5: ITS  = new AliITSv5asymm("ITS","Updates ITS TDR detailed version with asymmetric services");break;
    case 6: ITS  = new AliITSv5symm("ITS","Updated ITS TDR detailed version with symmetric services");break;
    // 
  }  

//=================== ITS parameters ============================
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
    //
    //
    //AliITS *ITS  = new AliITSv5symm("ITS","Updated ITS TDR detailed version with symmetric services");
    //
    //AliITS *ITS  = new AliITSv5asymm("ITS","Updates ITS TDR detailed version with asymmetric services");
    //
    //AliITSvPPRasymm *ITS  = new AliITSvPPRasymm("ITS","New ITS PPR detailed version with asymmetric services");
/*
    ITS->SetMinorVersion(2);					 // don't touch this parameter if you're not an ITS developer
    ITS->SetReadDet(kFALSE);					 // don't touch this parameter if you're not an ITS developer
    //    ITS->SetWriteDet("$ALICE_ROOT/ITS/ITSgeometry_vPPRasymm2.det");  // don't touch this parameter if you're not an ITS developer
    ITS->SetThicknessDet1(200.);   // detector thickness on layer 1 must be in the range [100,300]
    ITS->SetThicknessDet2(200.);   // detector thickness on layer 2 must be in the range [100,300]
    ITS->SetThicknessChip1(200.);  // chip thickness on layer 1 must be in the range [150,300]
    ITS->SetThicknessChip2(200.);  // chip thickness on layer 2 must be in the range [150,300]
    ITS->SetRails(1);	     // 1 --> rails in ; 0 --> rails out
    ITS->SetCoolingFluid(1);   // 1 --> water ; 0 --> freon
*/
    //
    //AliITSvPPRsymm *ITS  = new AliITSvPPRsymm("ITS","New ITS PPR detailed version with symmetric services");
    //ITS->SetMinorVersion(2);                                       // don't touch this parameter if you're not an ITS developer
    //ITS->SetReadDet(kFALSE);                                       // don't touch this parameter if you're not an ITS developer
    //ITS->SetWriteDet("$ALICE_ROOT/ITS/ITSgeometry_vPPRsymm2.det"); // don't touch this parameter if you're not an ITS developer
    //ITS->SetThicknessDet1(200.);   // detector thickness on layer 1 must be in the range [100,300]
    //ITS->SetThicknessDet2(200.);   // detector thickness on layer 2 must be in the range [100,300]
    //ITS->SetThicknessChip1(200.);  // chip thickness on layer 1 must be in the range [150,300]
    //ITS->SetThicknessChip2(200.);  // chip thickness on layer 2 must be in the range [150,300]
    //ITS->SetRails(1);              // 1 --> rails in ; 0 --> rails out
    //ITS->SetCoolingFluid(1);       // 1 --> water ; 0 --> freon
    //
    //
    // Coarse geometries (warning: no hits are produced with these coarse geometries and they unuseful 
    // for reconstruction !):
    //                                                     
    //
    //AliITSvPPRcoarseasymm *ITS  = new AliITSvPPRcoarseasymm("ITS","New ITS PPR coarse version with asymmetric services");
    //ITS->SetRails(1);                // 1 --> rails in ; 0 --> rails out
    //ITS->SetSupportMaterial(0);      // 0 --> Copper ; 1 --> Aluminum ; 2 --> Carbon
    //
    //AliITS *ITS  = new AliITSvPPRcoarsesymm("ITS","New ITS PPR coarse version with symmetric services");
    //ITS->SetRails(1);                // 1 --> rails in ; 0 --> rails out
    //ITS->SetSupportMaterial(0);      // 0 --> Copper ; 1 --> Aluminum ; 2 --> Carbon
    //                      
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
