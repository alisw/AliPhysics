enum gentype_t {hijing, gun, box, pythia, param, cocktail, fluka, halo, ntuple, scan, doublescan};

gentype_t gentype=param;
//Int_t ntracks=2828*20;
Int_t ntracks=200;
void Config()
{

new AliGeant3("C++ Interface to Geant3");

//=======================================================================
//  Create the output file
   
TFile *rootfile = new TFile("galice.root","recreate");
rootfile->SetCompressionLevel(2);
TGeant3 *geant3 = (TGeant3*)gMC;
 
 
//=======================================================================
// ******* GEANT STEERING parameters FOR ALICE SIMULATION *******
 geant3->fGctrak->maxnst=1000000;
 
geant3->SetTRIG(1);          //Number of events to be processed 
geant3->SetSWIT(4,100);
geant3->SetDEBU(0,0,1);
geant3->SetDCAY(1);
geant3->SetPAIR(1);
geant3->SetCOMP(1);
geant3->SetPHOT(1);
geant3->SetPFIS(0);
geant3->SetDRAY(1);
geant3->SetANNI(1);
geant3->SetBREM(1);
geant3->SetMUNU(1);
geant3->SetCKOV(0);
geant3->SetHADR(4); //Select pure GEANH (HADR 1) or GEANH/NUCRIN (HADR 3)
geant3->SetLOSS(1);
geant3->SetMULS(1);
geant3->SetRAYL(0);
geant3->SetAUTO(0); //Select automatic STMIN etc... calc. (AUTO 1) or manual (AUTO 0)
geant3->SetABAN(0); //Restore 3.16 behaviour for abandoned tracks
geant3->SetOPTI(2); //Select optimisation level for GEANT geometry searches (0,1,2)
Float_t cut    = 1.e-1; // 100MeV cut by default
Float_t tofmax = 1.e10;
//              GAM    ELEC   NHAD   CHAD   MUON  EBREM  MUHAB EDEL MUDEL MUPA TOFMAX
geant3->SetCUTS(1.e-4, 1.e-4, 1.e-3, 1.e-4, 1.e-4, cut,  cut,  cut, cut,  cut, 1.e-5);

 gAlice->TrackingLimits( 700, 2200);
 
//
//=======================================================================
// ************* STEERING parameters FOR ALICE SIMULATION **************
// --- Specify event type to be tracked through the ALICE setup
// --- All positions are in cm, angles in degrees, and P and E in GeV

 switch(gentype)
 {
 case gun:
//*********************************************
// Example for Fixed Particle Gun             *
//*********************************************
     AliGenFixed *gener = new AliGenFixed(ntracks);
     gener->SetMomentum(3);
     gener->SetPhiRange(0);
     gener->SetThetaRange(0);
     gener->SetOrigin(0,0,0);          //vertex position
     gener->SetPart(kProton);          //GEANT particle type
     break;
 case box:  
//*********************************************
// Example for Moving Particle Gun            *
//*********************************************
     AliGenBox *gener = new AliGenBox(ntracks);
     gener->SetMomentumRange(3,4);
     gener->SetPhiRange(-360,360);
     gener->SetThetaRange(2., 10.);
     gener->SetOrigin(25,25,510.5);   
     gener->SetVertexSmear(kPerTrack); 
     //vertex position
     gener->SetSigma(1.8, 1.8,0);           //Sigma in (X,Y,Z) (cm) on IP position
     gener->SetPart(kProton);                    //GEANT particle type
     break;
 case scan:  
//*********************************************
// Scanning on a grid                         *
//*********************************************
     AliGenScan *gener = new AliGenScan(-1);
     gener->SetMomentumRange(20,20);
     gener->SetPhiRange(90,90);
     gener->SetThetaRange(0,0);
     //vertex position
     gener->SetSigma(3,3,0);           //Sigma in (X,Y,Z) (cm) on IP position
     gener->SetPart(kMuonMinus); 
     gener->SetRange(20, -100, 100, 20, -100, 100, 1, 500, 500);
     break;
 case doublescan:  
//*********************************************
// Scanning on a grid                         *
//*********************************************
     AliGenDoubleScan *gener = new AliGenDoubleScan(-1);
     gener->SetMomentumRange(4,4);
     gener->SetPhiRange(0,360);
     gener->SetThetaRange(0,0);
     //vertex position
     gener->SetSigma(3,3,0);           //Sigma in (X,Y,Z) (cm) on IP position
     gener->SetPart(8); 
     gener->SetRange(20, -100, 100, 20, -100, 100, 1, 500, 500);
     gener->SetDistance(1);
     
     break;
     
 case hijing:
     AliGenHIJINGpara *gener = new AliGenHIJINGpara(ntracks);
     gener->SetMomentumRange(0,999);
     gener->SetPhiRange(0,360);
     gener->SetThetaRange(0.104,33.52);
     gener->SetOrigin(0,0,0);          //vertex position
     gener->SetSigma(0,0,5.6);         //Sigma in (X,Y,Z) (cm) on IP position
     break;

 case pythia:
//********************************************
// Example for Charm  Production with Pythia *
//********************************************

     AliGenPythia *gener = new AliGenPythia(ntracks);
     gener->SetMomentumRange(0,999);
     gener->SetPhiRange(0,360);
     gener->SetThetaRange(0., 180.);
     gener->SetYRange(-10,10);
     gener->SetPtRange(0,100);
     gener->SetOrigin(0,0,0);          // vertex position
     gener->SetVertexSmear(kPerEvent); 
     gener->SetSigma(0,0,5.6);         // Sigma in (X,Y,Z) (cm) on IP position
//     gener->SetStrucFunc(DO_Set_1);
     gener->SetProcess(mb); 
     gener->SetEnergyCMS(5500.);
     break;
     
 case param:
//*******************************************************
// Example for J/psi  Production from  Parameterisation *
//*******************************************************
     AliGenParam *gener = new AliGenParam(ntracks,upsilon_p,);
     gener->SetMomentumRange(0,999);
     gener->SetPhiRange(0,360);
     gener->SetYRange(2.5,4);
     gener->SetThetaRange(2,9);
     gener->SetPtRange(0,10);
     gener->SetOrigin(0,0,0);          //vertex position
     gener->SetSigma(0,0,0);           //Sigma in (X,Y,Z) (cm) on IP position
     gener->SetForceDecay(dimuon);
     gener->SetCutOnChild(0);
     gener->SetTrackingFlag(0);

     break;
     
 case fluka:
//*******************************************************
// Example for a FLUKA Boundary Source                  *
//*******************************************************
     AliGenFLUKAsource *gener = new AliGenFLUKAsource(-1);
     gener->AddFile("$(ALICE_ROOT)/data/all32.root"); 
     rootfile->cd();
     gener->SetPartFlag(7);
     gener->SetMomentumRange(0,999);
     gener->SetPhiRange(0,360);
     gener->SetThetaRange(0., 180.); 
     gener->SetAgeMax(1.e-5);
     
//  31.7 events     
     gener->SetFraction(0.0315);     
     break;

 case ntuple:
//*******************************************************
// Example for reading from a external file                  *
//*******************************************************
     AliGenExtFile *gener = new AliGenExtFile(-1); 
     gener->SetFileName("$(ALICE_ROOT)/data/pbpb.root");
     gener->SetThetaRange(0.104,33.52);
     gener->SetOrigin(0,0,0);          //vertex position
     gener->SetSigma(0,0,5.6);         //Sigma in (X,Y,Z) (cm) on IP position
     gener->SetVertexSmear(kPerTrack); 
     gener->SetTrackingFlag(1);
     break;

 case halo:
//*******************************************************
// Example for Tunnel Halo Source                       *
//*******************************************************
     AliGenHalo *gener = new AliGenHalo(ntracks); 
     gener->SetFileName("/h1/morsch/marsip/marsip5.mu");
     break;
     
 case cocktail:
//*******************************************************
// Example for a Cocktail                               *
//*******************************************************
     AliGenCocktail *gener = new AliGenCocktail();
     gener->SetMomentumRange(0,10);
     gener->SetPhiRange(0,360);
     gener->SetThetaRange(45.,135);
    
     pions   = new AliGenParam(100, pion_p);
//     kaons   = new AliGenParam(10 , kaon_p);
//     protons = new AliGenParam(10 , proton_p);
     gener->AddGenerator(pions  , "Pions"  , 100);
//     gener->AddGenerator(kaons  , "Kaons"  , 10);
//     gener->AddGenerator(protons, "Protons", 10);
	 
     break;
 }
 
// Activate this line if you want the vertex smearing to happen
// track by track
//

gener->Init();
gAlice->SetField(-2,1);    //Specify maximum magnetic field in Tesla (neg. ==> default field)
Int_t iFRAME  =0;
Int_t iMAG    =1;
Int_t iITS    =0;
Int_t iTPC    =0;
Int_t iTOF    =0;
Int_t iRICH   =0;
Int_t iZDC    =0;
Int_t iCASTOR =0;
Int_t iTRD    =0;
Int_t iABSO   =0;
Int_t iDIPO   =1;
Int_t iHALL   =1;
Int_t iSHIL   =0;
Int_t iPIPE   =1;
Int_t iFMD    =0;
Int_t iMUON   =1;
Int_t iPHOS   =0;
Int_t iPMD    =0;

//=================== Alice BODY parameters =============================
AliBODY *BODY = new AliBODY("BODY","Alice envelop");

if(iFRAME) {
//=================== FRAME parameters ============================
AliFRAME *FRAME  = new AliFRAMEv0("FRAME","Space Frame");
}

if(iMAG) {
//=================== MAG parameters ============================
// --- Start with Magnet since detector layouts may be depending ---
// --- on the selected Magnet dimensions ---
AliMAG *MAG  = new AliMAG("MAG","Magnet");
}

if(iITS) {
//=================== ITS parameters ============================
//
// EUCLID is a flag to output (=1) both geometry and media to two ASCII files 
// (called by default ITSgeometry.euc and ITSgeometry.tme) in a format
// understandable to the CAD system EUCLID. The default (=0) means that you 
// dont want to use this facility.
//
    AliITS *ITS  = new AliITSv5("ITS","normal ITS");
    ITS->SetEUCLID(0);
}

if(iTPC) {
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

    AliTPC *TPC  = new AliTPCv3("TPC","Normal TPC");
    AliTPCD *paramd = TPC->GetDigParam();
    AliTPCParam *param = &(paramd->GetParam());
    
// Set geometrical parameters

    param->SetSectorAngles(40.,0.,20.,10.);
    param->SetInnerRadiusLow(83.7);
    param->SetInnerRadiusUp(132.9);
    param->SetOuterRadiusLow(146.9);
    param->SetOuterRadiusUp(249.4);
    param->SetInSecLowEdge(81.6);
    param->SetInSecUpEdge(135.);
    param->SetOuSecLowEdge(144.2);
    param->SetOuSecUpEdge(252.1);
    param->SetEdge(1.5);
    param->SetDeadZone(1.15);
    param->Update();
    
// set gas mixture

    TPC->SetGasMixt(2,20,10,-1,0.9,0.1,0.);
    TPC->SetSecAL(1);
    TPC->SetSecAU(1);
    TPC->SetSecLows(0, -1, -1, -1, -1, -1);
    TPC->SetSecUps(18, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
    TPC->SetSens(-1);
}

if(iTOF) {
//=================== TOF parameters ============================
    AliTOF *TOF  = new AliTOFv2("TOF","normal TOF");
}

if(iRICH) {
//=================== RICH parameters ===========================
    AliRICH *RICH  = new AliRICHv0("RICH","normal RICH");
    
    RICH->SetSMAXAR(0.03);
    RICH->SetSMAXAL(-1);
//
// Version 0
// Default Segmentation
    AliRICHsegmentationV0* RsegV0 = new AliRICHsegmentationV0;
    RsegV0->SetPadSize(.8, .8);
    RsegV0->SetDAnod(0.8/3);
// Default Response
    AliRICHresponseV0* Rresponse0 = new AliRICHresponseV0;
    AliRICHresponseCkv* RresponseCkv = new AliRICHresponseCkv;
    
//------------------------Chambers 0-6 ----------------------------
  for (Int_t i=0; i<7; i++) {
      RICH->SetSegmentationModel(i, 1, RsegV0);
      RICH->SetResponseModel(i, mip     , Rresponse0);
      RICH->SetResponseModel(i, cerenkov, RresponseCkv);
      RICH->Chamber(i).SetRSIGM(5.);
      RICH->Chamber(i).SetMUCHSP(20.);
      RICH->Chamber(i).SetMUSIGM(0.18, 0.18);
      RICH->Chamber(i).SetMAXADC( 1024);
      RICH->Chamber(i).SetSqrtKx3(0.77459667);
      RICH->Chamber(i).SetKx2(0.962);
      RICH->Chamber(i).SetKx4(0.379);
      RICH->Chamber(i).SetSqrtKy3(0.77459667);
      RICH->Chamber(i).SetKy2(0.962);
      RICH->Chamber(i).SetKy4(0.379);
      RICH->Chamber(i).SetPitch(0.25);
      RICH->SetNsec(i,1);
  }
}
 
if(iZDC) {
//=================== ZDC parameters ============================
    
    AliZDC *ZDC  = new AliZDCv1("ZDC","normal ZDC");
}

if(iCASTOR) {
//=================== CASTOR parameters ============================

    AliCASTOR *CASTOR  = new AliCASTORv1("CASTOR","normal CASTOR");
}

if(iTRD) {
//=================== TRD parameters ============================
    
    AliTRD *TRD  = new AliTRDv2("TRD","TRD version 2");
}


if(iABSO) {
//=================== ABSO parameters ============================
    AliABSO *ABSO  = new AliABSOv1("ABSO","Muon Absorber");
}

if(iDIPO) {
//=================== DIPO parameters ============================

    AliDIPO *DIPO  = new AliDIPOv2("DIPO","Dipole version 2");
}

if(iHALL) {
//=================== HALL parameters ============================
    AliHALL *HALL  = new AliHALL("HALL","Alice Hall");
}



if(iSHIL) {
//=================== SHIL parameters ============================
    AliSHIL *SHIL  = new AliSHILv1("SHIL","Shielding");
}


if(iPIPE) {
//=================== PIPE parameters ============================
    AliPIPE *PIPE  = new AliPIPEv0("PIPE","Beam Pipe");
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
 

if(iPHOS) {
//=================== PHOS parameters ===========================

AliPHOS *PHOS  = new AliPHOSv1("PHOS","normal PHOS");
// * PHOSflags:    YES: X<>0   NO: X=0
// * PHOSflags(1) : -----X  Create branch for TObjArray of AliPHOSCradle
// *                ----X-  Create file (ftn03 on HP-UX) with list of SHAKER particles (7Mb/event)
// *                
PHOS->SetFlags(000001);
PHOS->SetRadius(460); //Distance from beam to PHOS crystals.
// (crystal_side_size,crystal_length,wrap_thikness,air_thikness,PIN_size,PIN length)
PHOS->SetCell(2.2,          18.,         0.01,        0.01,        1.,      0.1);
PHOS->SetCradleSize(104, 88, 4); // Nz (along beam), Nphi, Ncradles
PHOS->SetCradleA(0);   //Angle between Cradles
PHOS->SetCPV(1., 2.); //CPV thikness, CPV-PHOS distance
// *  ===============
// * PHOS extra parameters (contact Maxim Volkov volkov@mail.cern.ch)
// * 1. STE_THICK         Steel cover thickness
// * 2. SUP_Y             Crystal support height
// * 3. FTIU_THICK        Thermo Insulating outer cover Upper plate thickness
// * 4. UFP_Y             Upper Polystyrene Foam plate thickness
// * 5. TCB_THICK         Thermo insulating Crystal Block wall thickness
// * 6. UCP_Y             Upper Cooling Plate thickness
// * 7. ASP_Y             Al Support Plate thickness
// * 8. TIP_Y             Lower Thermo Insulating Plate thickness
// * 9. TXP_Y             Lower Textolit Plate thickness
PHOS->SetExtra(0.001, 6.95, 4., 5., 2., 0.06, 10., 3., 1.);   
PHOS->SetTextolitWall(209., 71., 250.);    //Textolit Wall box dimentions
PHOS->SetInnerAir(206.,    66.,     244.); //Inner AIR volume dimensions
// *  ===============================
// * 1. FTI_X             Foam Thermo Insulating outer cover dimensions
// * 2. FTI_Y             ==//==
// * 3. FTI_Z             ==//==
// * 4. FTI_R             Distance from IP to Foam Thermo Insulating top plate
PHOS->SetFoam(214.6,  80.,  260., 467.); 
//    =================================
// *******************************************************************************
// * KINE 700  - SHAKER generator
// * KINE 700 x y z NDNDY YLIM PTLIM ChargeFlag
// *     JWEAK=0
// *     JPI0=JETA=1
// *     JPIC=JPRO=JKAC=JKA0=JRHO=JOME=JPHI=JPSI=JDRY=ChargeFlag
// *     Int_t               JWEI;           // Unweighted generation
// *     Int_t               NDNDY;          // Density of charged particles
// *     Float_t             YLIM;           // Rapidity Limit
// *     Float_t             PTLIM;          // Pt limit in GeV/c
// *     Int_t               JWEAK;          // Disable weak decays
// *     Int_t               JPI0;           // pi0 generation
// *     Int_t               JETA;           // eta generation
// *     Int_t               JPIC;           // pi+/- generation
// *     Int_t               JPRO;           // proton generation
// *     Int_t               JKAC;           // K+/- generation
// *     Int_t               JKA0;           // K0 generation
// *     Int_t               JRHO;           // rho generation
// *     Int_t               JOME;           // omega generation
// *     Int_t               JPHI;           // phi generation
// *     Int_t               JPSI;           // J/psi generation
// *     Int_t               JDRY;           // Drell-Yan generation
// * KINE  700     5.    175.    0.          800. 1.5 5. 1.
// *******************************************************************************
}

if(iPMD) {
//=================== PMD parameters ============================

//         Must be defined AFTER PHOS
AliPMD *PMD  = new AliPMDv1("PMD","normal PMD");
PMD->SetPAR(1., 1., 0.8, 0.02);
PMD->SetIN(6., 20., 600., 27., 27.);
PMD->SetGEO(0.0, 0.2, 4.);
PMD->SetPadSize(0.8, 1.0, 1.2, 1.5);
}
}
		



