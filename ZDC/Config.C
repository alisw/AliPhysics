 void Config()
 {

 new AliGeant3("C++ Interface to Geant3");

 //=======================================================================
 //  Create the output file

 TFile *rootfile = new TFile("galice.root","recreate");
 rootfile->SetCompressionLevel(2);
 TGeant3 *geant3 = (TGeant3*)gMC;
 //
 // Set External decayer
  AliDecayer* decayer = new AliDecayerPythia();
  decayer->SetForceDecay(all);
  decayer->Init();
  gMC->SetExternalDecayer(decayer);
 //
 //
 //=======================================================================
 // ******* GEANT STEERING parameters FOR ALICE SIMULATION *******
 geant3->SetTRIG(1); //Number of events to be processed
 geant3->SetSWIT(4,10);
 geant3->SetDEBU(0,0,1);
 //geant3->SetSWIT(2,2);
 //geant3->SetSWIT(2,3); //for drawing
 geant3->SetDCAY(1);
 geant3->SetPAIR(1);
 geant3->SetCOMP(1);
 geant3->SetPHOT(1);
 geant3->SetPFIS(0);
 geant3->SetDRAY(0);
 geant3->SetANNI(1);
 geant3->SetBREM(1);
 geant3->SetMUNU(1);
 geant3->SetCKOV(1);
 geant3->SetHADR(1); //Select pure GEANH (HADR 1) or GEANH/NUCRIN (HADR 3)
 geant3->SetLOSS(2);
 geant3->SetMULS(1);
 geant3->SetRAYL(1);
 geant3->SetAUTO(1); //Select automatic STMIN etc... calc. (AUTO 1) or manual (AUTO 0)
 geant3->SetABAN(0); //Restore 3.16 behaviour for abandoned tracks
 geant3->SetOPTI(2); //Select optimisation level for GEANT geometry searches (0,1,2)
 geant3->SetERAN(5.e-7);

 Float_t cut	= 1.e-3; // 1MeV cut by default
 Float_t tofmax = 1.e10;
 //		GAM ELEC NHAD CHAD MUON EBREM MUHAB EDEL MUDEL MUPA TOFMAX
 geant3->SetCUTS(cut,cut, cut, cut, cut, cut,  cut,  cut, cut,  cut, tofmax);
 //
 //=======================================================================
 // ************* STEERING parameters FOR ALICE SIMULATION **************
 // --- Specify event type to be tracked through the ALICE setup
 // --- All positions are in cm, angles in degrees, and P and E in GeV
//
// ####  AliGenZDC generation
//
AliGenZDC *gener = new AliGenZDC(1);
gener->SetDirection(0,0,0,1);
gener->SetFermi(1);--------------------//Nucleon in ZP or ZN
gener->SetDiv(0.000032,0.0001,2);
gener->SetOrigin(0.,0.,0.);
//gener->SetParticle(kNeutron);
gener->SetParticle(kProton);
gener->SetMomentum(2760.);
//gener->SetFermi(0); -----------------// Gamma in ZEM
//gener->SetDiv(0,0,0);
//gener->SetOrigin(0.,5.8,11400.);
//gener->SetParticle(kGamma);
//gener->SetMomentum(100.);	
//
gener->SetTrackingFlag(1);
gener->Init();
 //
 // Activate this line if you want the vertex smearing to happen
 // track by track
 //
 //gener->SetVertexSmear(perTrack);

 gAlice->SetField(-999,2);    //Specify maximum magnetic field in Tesla (neg. ==> default field)

 Int_t iABSO=0;
 Int_t iCASTOR=0;
 Int_t iDIPO=1;
 Int_t iFMD=0;
 Int_t iFRAME=0;
 Int_t iHALL=0;
 Int_t iITS=0;
 Int_t iMAG=0;
 Int_t iMUON=0;
 Int_t iPHOS=0;
 Int_t iPIPE=0;
 Int_t iPMD=0;
 Int_t iHMPID=0;
 Int_t iSHIL=0;
 Int_t iSTART=0;
 Int_t iTOF=0;
 Int_t iTPC=0;
 Int_t iTRD=0;
 Int_t iZDC=1;

 //=================== Alice BODY parameters =============================
 AliBODY *BODY = new AliBODY("BODY","Alice envelop");

 if(iMAG) {
 //=================== MAG parameters ============================
 // --- Start with Magnet since detector layouts may be depending ---
 // --- on the selected Magnet dimensions ---
 AliMAG *MAG  = new AliMAG("MAG","Magnet");
 }

 if(iABSO) {
 //=================== ABSO parameters ============================
 AliABSO *ABSO  = new AliABSOv0("ABSO","Muon Absorber");
 }

 if(iDIPO) {
 //=================== DIPO parameters ============================

 AliDIPO *DIPO  = new AliDIPOv2("DIPO","Dipole version 2");
 }

 if(iHALL) {
 //=================== HALL parameters ============================

 AliHALL *HALL  = new AliHALL("HALL","Alice Hall");
 }

 if(iFRAME) {
 //=================== FRAME parameters ============================

 AliFRAME *FRAME  = new AliFRAMEv1("FRAME","Space Frame");

 }

 if(iSHIL) {
 //=================== SHIL parameters ============================

 AliSHIL *SHIL  = new AliSHILv0("SHIL","Shielding");
 }

 if(iPIPE) {
 //=================== PIPE parameters ============================

 AliPIPE *PIPE  = new AliPIPEv0("PIPE","Beam Pipe");
 // PIPE->Dump();
 }

 if(iITS) {
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
 //
 //
 //AliITS *ITS  = new AliITSv3("ITS","Old ITS detailed version as of the ALICE TP");
 //
 //AliITS *ITS  = new AliITSv5("ITS","Current ITS detailed version used for the ITS TDR");
 //
 //AliITS *ITS  = new AliITSv5symm("ITS","Updated ITS TDR detailed version with symmetric services");
 //
 AliITS *ITS  = new AliITSv5asymm("ITS","Updates ITS TDR detailed version with asymmetric services");
 //
 //
 // Coarse geometries (warning: no hits are produced with these coarse geometries and they unuseful for reconstruction !):
 // ======================================================================================================================
 //
 //
 //AliITS *ITS  = new AliITSv1("ITS","Old ITS coarse version as of the ALICE TP");
 //
 //AliITS *ITS  = new AliITSvPPRcoarseasymm("ITS","New ITS coarse version with asymmetric services");
 //
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

 if(iTPC) {
 //============================ TPC parameters ================================
 // --- This allows the user to specify sectors for the SLOW (TPC geometry 2)
 // --- Simulator. SecAL (SecAU) <0 means that ALL lower (upper)
 // --- sectors are specified, any value other than that requires at least one
 // --- sector (lower or upper)to be specified!
 // --- Reminder: sectors 1-24 are lower sectors (1-12 -> z>0, 13-24 -> z<0)
 // --- 	  sectors 25-72 are the upper ones (25-48 -> z>0, 49-72 -> z<0)
 // --- SecLows - number of lower sectors specified (up to 6)
 // --- SecUps - number of upper sectors specified (up to 12)
 // --- Sens - sensitive strips for the Slow Simulator !!!
 // --- This does NOT work if all S or L-sectors are specified, i.e.
 // --- if SecAL or SecAU < 0
 //
 //
 //-----------------------------------------------------------------------------

   //  gROOT->LoadMacro("SetTPCParam.C");
   //  AliTPCParam *param = SetTPCParam();
   AliTPC *TPC  = new AliTPCv2("TPC","Default"); //v1 is default
   //  TPC->SetParam(param); // pass the parameter object to the TPC

 // set gas mixture

   //TPC->SetGasMixt(2,20,10,-1,0.9,0.1,0.);
   TPC->SetSecAL(-1);
   TPC->SetSecAU(-1);
   //  TPC->Dump();
   //TPC->SetSecLows(1,  2,  3, 19, 20, 21);
   //TPC->SetSecUps(37, 38, 39, 37+18, 38+18, 39+18, -1, -1, -1, -1, -1, -1);
   //TPC->SetSens(1);

   //if (TPC->IsVersion()==1) param->Write(param->GetTitle());
 }

 if(iTOF) {
 //=================== TOF parameters ============================
   AliTOF *TOF  = new AliTOFv2("TOF","normal TOF");
 }

 if(iHMPID) {
 //=================== HMPID parameters ===========================
   AliHMPID *HMPID  = new AliHMPIDv1("HMPID","normal HMPID");

 }

 if(iZDC) {
 //=================== ZDC parameters ============================

   AliZDC *ZDC  = new AliZDCv1("ZDC","normal ZDC");
   AliZDCv1 *ZDCv1 = (AliZDC*)ZDC;
   //ZDCv1->NoShower();
   ZDCv1->Shower();
 }

 if(iCASTOR) {
 //=================== CASTOR parameters ============================

   AliCASTOR *CASTOR  = new AliCASTORv1("CASTOR","normal CASTOR");
 }

 if(iTRD) {
 //=================== TRD parameters ============================

   AliTRD *TRD  = new AliTRDv1("TRD","TRD slow simulator");

   // Select the gas mixture (0: 97% Xe + 3% isobutane, 1: 90% Xe + 10% CO2)
   TRD->SetGasMix(1);

   // With hole in front of PHOS
   TRD->SetPHOShole();
   // With hole in front of HMPID
   TRD->SetHMPIDhole();
 }

 if(iFMD) {
 //=================== FMD parameters ============================

   AliFMD *FMD  = new AliFMDv1("FMD","normal FMD");
 }

 if(iMUON) {
 //=================== MUON parameters ===========================

   AliMUON *MUON  = new AliMUONv1("MUON","normal MUON");
   MUON->SetIshunt(0);
   MUON->SetMaxStepGas(0.1);
   MUON->SetMaxStepAlu(0.1);
 //
 // Version 0
 //
 // First define the number of planes that are segmented (1 or 2) by a call
 // to SetNsec.
 // Then chose for each chamber (chamber plane) the segmentation
 // and response model.
 // They should be equal for the two chambers of each station. In a future
 // version this will be enforced.
 //
 //
   Int_t chamber;
   Int_t station;
 // Default response
   AliMUONResponseV0* response0 = new AliMUONResponseV0;
   response0->SetSqrtKx3(0.7131);
   response0->SetKx2(1.0107);
   response0->SetKx4(0.4036);
   response0->SetSqrtKy3(0.7642);
   response0->SetKy2(0.9706);
   response0->SetKy4(0.3831);
   response0->SetPitch(0.25);
   response0->SetSigmaIntegration(10.);
   response0->SetChargeSlope(50);
   response0->SetChargeSpread(0.18, 0.18);
   response0->SetMaxAdc(4096);
   response0->SetZeroSuppression(6);
 //--------------------------------------------------------
 // Configuration for Chamber TC1/2  (Station 1) ----------
 //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   Float_t rseg1[4]={17.5, 55.2, 71.3, 95.5};
   Int_t   nseg1[4]={4, 4, 2, 1};
   //
   chamber=1;
   //^^^^^^^^^
   MUON->SetNsec(chamber-1,2);
   //
   AliMUONSegmentationV01 *seg11=new AliMUONSegmentationV01;

   seg11->SetSegRadii(rseg1);
   seg11->SetPadSize(3, 0.5);
   seg11->SetDAnod(3.0/3./4);
   seg11->SetPadDivision(nseg1);

   MUON->SetSegmentationModel(chamber-1, 1, seg11);
 //
   AliMUONSegmentationV02 *seg12=new AliMUONSegmentationV02;
   seg12->SetSegRadii(rseg1);
   seg12->SetPadSize(0.75, 2.0);
   seg12->SetDAnod(3.0/3./4);
   seg12->SetPadDivision(nseg1);

   MUON->SetSegmentationModel(chamber-1, 2, seg12);

   MUON->SetResponseModel(chamber-1, response0);

   chamber=2;
 //^^^^^^^^^
 //
   MUON->SetNsec(chamber-1,2);
 //
   AliMUONSegmentationV01 *seg21=new AliMUONSegmentationV01;
   seg21->SetSegRadii(rseg1);
   seg21->SetPadSize(3, 0.5);
   seg21->SetDAnod(3.0/3./4);
   seg21->SetPadDivision(nseg1);
   MUON->SetSegmentationModel(chamber-1, 1, seg21);
   //
   AliMUONSegmentationV02 *seg22=new AliMUONSegmentationV02;
   seg22->SetSegRadii(rseg1);
   seg22->SetPadSize(0.75, 2.);
   seg22->SetDAnod(3.0/3./4);
   seg22->SetPadDivision(nseg1);
   MUON->SetSegmentationModel(chamber-1, 2, seg22);

   MUON->SetResponseModel(chamber-1, response0);
   //
   //--------------------------------------------------------
   // Configuration for Chamber TC3/4 -----------------------
   ///^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   // Float_t rseg2[4]={23.5, 87.7, 122.4, 122.5};
   Float_t rseg2[4]={23.5, 47.1, 87.7, 122.5};
   Int_t   nseg2[4]={4, 4, 2, 1};
   //
   chamber=3;
   //^^^^^^^^^
   MUON->SetNsec(chamber-1,2);
   //
   AliMUONSegmentationV01 *seg31=new AliMUONSegmentationV01;
   seg31->SetSegRadii(rseg2);
   seg31->SetPadSize(6, 0.5);
   seg31->SetDAnod(3.0/3./4);
   seg31->SetPadDivision(nseg2);
   MUON->SetSegmentationModel(chamber-1, 1, seg31);
   //
   AliMUONSegmentationV02 *seg32=new AliMUONSegmentationV02;
   seg32->SetSegRadii(rseg2);
   seg32->SetPadSize(0.75, 4.);
   seg32->SetPadDivision(nseg2);
   seg32->SetDAnod(3.0/3./4);

   MUON->SetSegmentationModel(chamber-1, 2, seg32);

   MUON->SetResponseModel(chamber-1, response0);

   chamber=4;
   //^^^^^^^^^
   //
   MUON->SetNsec(chamber-1,2);
   //
   AliMUONSegmentationV01 *seg41=new AliMUONSegmentationV01;
   seg41->SetSegRadii(rseg2);
   seg41->SetPadSize(6, 0.5);
   seg41->SetDAnod(3.0/3./4);
   seg41->SetPadDivision(nseg2);
   MUON->SetSegmentationModel(chamber-1, 1, seg41);
   //
   AliMUONSegmentationV02 *seg42=new AliMUONSegmentationV02;
   seg42->SetSegRadii(rseg2);
   seg42->SetPadSize(0.75, 4.);
   seg42->SetPadDivision(nseg2);
   seg42->SetDAnod(3.0/3./4);

   MUON->SetSegmentationModel(chamber-1, 2, seg42);

   MUON->SetResponseModel(chamber-1, response0);

 //--------------------------------------------------------
 // Configuration for Chamber TC5/6 -----------------------
 //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   seg5 =  new AliMUONSegmentationV1;
   AliMUONResponseV0* response5 =  new AliMUONResponseV0;
   // K3 = 0.62
   response5->SetSqrtKx3(0.78740079);
   response5->SetKx2(0.95237319); //  0.5 * kPI * (1- 0.5*sqrtky3 )
   response5->SetKx4(0.37480633); //  0.25/TMath::ATan(sqrtkx3)
   // K3 = 0.55
   response5->SetSqrtKy3(0.74161985);
   response5->SetKy2(0.98832946);
   response5->SetKy4(0.39177817);
   response5->SetPitch(0.325);
   response5->SetSigmaIntegration(10.);
   response5->SetChargeSlope(50);
   response5->SetChargeSpread(0.4, 0.4);
   response5->SetMaxAdc(4096);
   response5->SetZeroSuppression(6);


   chamber=5;
   MUON->SetNsec(chamber-1,1);
   MUON->SetSegmentationModel(chamber-1, 1, seg5);
   MUON->SetResponseModel(chamber-1, response5);

   chamber=6;
   MUON->SetNsec(chamber-1,1);
   MUON->SetSegmentationModel(chamber-1, 1, seg5);
   MUON->SetResponseModel(chamber-1, response5);
   //
   // Station 3
   station=3;
   MUON->SetPadSize(station, 1, 0.975, 0.55);

 //--------------------------------------------------------
 // Configuration for Chamber TC7/8  (Station 4) ----------
 //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   Int_t   nseg4[4]={4, 4, 2, 1};

   chamber=7;
   //^^^^^^^^^
   MUON->SetNsec(chamber-1,2);
   //
   AliMUONSegmentationV04 *seg71=new AliMUONSegmentationV04;
   seg71->SetPadSize(10.,0.5);
   seg71->SetDAnod(0.25);
   seg71->SetPadDivision(nseg4);
   MUON->SetSegmentationModel(chamber-1, 1, seg71);
   AliMUONSegmentationV05 *seg72=new AliMUONSegmentationV05;
   seg72->SetPadSize(1,10);
   seg72->SetDAnod(0.25);
   seg72->SetPadDivision(nseg4);
   MUON->SetSegmentationModel(chamber-1, 2, seg72);

   MUON->SetResponseModel(chamber-1, response0);

   chamber=8;
   //^^^^^^^^^
   MUON->SetNsec(chamber-1,2);
   AliMUONSegmentationV04 *seg81=new AliMUONSegmentationV04;
   seg81->SetPadSize(10., 0.5);
   seg81->SetPadDivision(nseg4);
   seg81->SetDAnod(0.25);
   MUON->SetSegmentationModel(chamber-1, 1, seg81);

   AliMUONSegmentationV05 *seg82=new AliMUONSegmentationV05;
   seg82->SetPadSize(1, 10);
   seg82->SetPadDivision(nseg4);
   seg82->SetDAnod(0.25);
   MUON->SetSegmentationModel(chamber-1, 2, seg82);

   MUON->SetResponseModel(chamber-1, response0);
   //--------------------------------------------------------
   // Configuration for Chamber TC9/10  (Station 5) ---------
   //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   chamber=9;
   //^^^^^^^^^
   MUON->SetNsec(chamber-1,2);
   //
   AliMUONSegmentationV04 *seg91=new AliMUONSegmentationV04;
   seg91->SetPadSize(10.,0.5);
   seg91->SetDAnod(0.25);
   seg91->SetPadDivision(nseg4);
   MUON->SetSegmentationModel(chamber-1, 1, seg91);

   AliMUONSegmentationV05 *seg92=new AliMUONSegmentationV05;
   seg92->SetPadSize(1,10);
   seg92->SetDAnod(0.25);
   seg92->SetPadDivision(nseg4);

   MUON->SetSegmentationModel(chamber-1, 2, seg92);

   MUON->SetResponseModel(chamber-1, response0);

   chamber=10;
   //^^^^^^^^^
   MUON->SetNsec(chamber-1,2);
   AliMUONSegmentationV04 *seg101=new AliMUONSegmentationV04;
   seg101->SetPadSize(10., 0.5);
   seg101->SetPadDivision(nseg4);
   seg101->SetDAnod(0.25);
   MUON->SetSegmentationModel(chamber-1, 1, seg101);

   AliMUONSegmentationV05 *seg102=new AliMUONSegmentationV05;
   seg102->SetPadSize(1,10);
   seg102->SetPadDivision(nseg4);
   seg102->SetDAnod(0.25);
   MUON->SetSegmentationModel(chamber-1, 2, seg102);

   MUON->SetResponseModel(chamber-1, response0);

 //--------------------------------------------------------
 // Configuration for Trigger staions ---------------------
 //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   AliMUONResponseTrigger* responseTrigger0 =  new AliMUONResponseTrigger;

   chamber=11;
   MUON->SetNsec(chamber-1,2);
   AliMUONSegmentationTriggerX *seg111=new AliMUONSegmentationTriggerX;
   MUON->SetSegmentationModel(chamber-1, 1, seg111);
   AliMUONSegmentationTriggerY *seg112=new AliMUONSegmentationTriggerY;
   MUON->SetSegmentationModel(chamber-1, 2, seg112);

   MUON->SetResponseModel(chamber-1, responseTrigger0);

   chamber=12;
   MUON->SetNsec(chamber-1,2);
   AliMUONSegmentationTriggerX *seg121=new AliMUONSegmentationTriggerX;
   MUON->SetSegmentationModel(chamber-1, 1, seg121);
   AliMUONSegmentationTriggerY *seg122=new AliMUONSegmentationTriggerY;
   MUON->SetSegmentationModel(chamber-1, 2, seg122);

   MUON->SetResponseModel(chamber-1, responseTrigger0);

   chamber=13;
   MUON->SetNsec(chamber-1,2);
   AliMUONSegmentationTriggerX *seg131=new AliMUONSegmentationTriggerX;
   MUON->SetSegmentationModel(chamber-1, 1, seg131);
   AliMUONSegmentationTriggerY *seg132=new AliMUONSegmentationTriggerY;
   MUON->SetSegmentationModel(chamber-1, 2, seg132);
   MUON->SetResponseModel(chamber-1, responseTrigger0);

   chamber=14;
   MUON->SetNsec(chamber-1,2);
   AliMUONSegmentationTriggerX *seg141=new AliMUONSegmentationTriggerX;
   MUON->SetSegmentationModel(chamber-1, 1, seg141);
   AliMUONSegmentationTriggerY *seg142=new AliMUONSegmentationTriggerY;
   MUON->SetSegmentationModel(chamber-1, 2, seg142);

   MUON->SetResponseModel(chamber-1, responseTrigger0);
 }

 //=================== PHOS parameters ===========================

 if(iPHOS) {
   AliPHOS *PHOS  = new AliPHOSv1("PHOS","GPS2");
 }

 if(iPMD) {
 //=================== PMD parameters ============================

   AliPMD *PMD  = new AliPMDv1("PMD","normal PMD");
   PMD->SetPAR(1., 1., 0.8, 0.02);
   PMD->SetIN(6., 18., -580., 27., 27.);
   PMD->SetGEO(0.0, 0.2, 4.);
   PMD->SetPadSize(0.8, 1.0, 1.0, 1.5);

 }

 if(iSTART) {
 //=================== START parameters ============================
   AliSTART *START  = new AliSTARTv1("START","START Detector");
 }


 }
