void Config()
{

new AliGEODB("Test with the GEO database");

gSystem->Load("/home/dcollado/AliRoot/pro/lib/libGEODB");
//=======================================================================
//  Create the output file
   
//TFile *rootfile = new TFile("GeoDB.root","recreate");
//rootfile->SetCompressionLevel(2);

AliGEODB* geant3 = (AliGEODB*) gMC;

Int_t iMAG=0;
Int_t iITS=0;
Int_t iTPC=1;
Int_t iTOF=0;
Int_t iRICH=0;
Int_t iZDC=0;
Int_t iCASTOR=0;
Int_t iTRD=0;
Int_t iABSO=0;
Int_t iDIPO=0;
Int_t iHALL=0;
Int_t iFRAME=0;
Int_t iSHIL=0;
Int_t iPIPE=0;
Int_t iFMD=0;
Int_t iMUON=0;
Int_t iPHOS=0;
Int_t iPMD=0;

//=================== Alice BODY parameters =============================

//AliBODY *BODY = new AliBODY("BODY","Alice envelop");


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
AliITS *ITS  = new AliITSv3("ITS","normal ITS");
ITS->SetEUCLID(1);
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

    AliTPC *TPC  = new AliTPCv1("TPC","Normal TPC");
    TPC->SetSecAL(1);
    TPC->SetSecAU(1);
    TPC->SetSecLows(1, -1, -1, -1, -1, -1);
    TPC->SetSecUps(25, 26, 48, -1, -1, -1, -1, -1, -1, -1, -1, -1);
    TPC->SetSens(1);
}

if(iTOF) {
//=================== TOF parameters ============================
AliTOF *TOF  = new AliTOFv2("TOF","normal TOF");
}

if(iRICH) {
//=================== RICH parameters ===========================


AliRICH *RICH  = new AliRICHv1("RICH","normal RICH");
RICH->SetSP(40);
RICH->SetFEED(0.04);
RICH->SetSIGM(0.18);
RICH->SetTRIG(0);
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
AliABSO *ABSO  = new AliABSO("ABSO","Muon Absorber");
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
AliFRAME *FRAME  = new AliFRAMEv0("FRAME","Space Frame");
}

if(iSHIL) {
//=================== SHIL parameters ============================

AliSHIL *SHIL  = new AliSHIL("SHIL","Shielding");
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

AliMUON *MUON  = new AliMUONv0("MUON","normal MUON");

MUON->SetSMAXAR(0.03);
MUON->SetSMAXAL(-1);
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
// Default Segmentation
 AliMUONsegmentationV0* segV0 = new AliMUONsegmentationV0;
// Default response
 AliMUONresponseV0* response0 = new AliMUONresponseV0;
 response0->SetSqrtKx3(0.761577);
 response0->SetKx2(0.972655);
 response0->SetKx4(0.3841);
 response0->SetSqrtKy3(0.714143);
 response0->SetKy2(1.0099);
 response0->SetKy4(0.403);
 response0->SetPitch(0.25);
 response0->SetRSIGM(10.);
 response0->SetMUCHSP(5.);
 response0->SetMUSIGM(0.18, 0.18);
 response0->SetMAXADC( 1024);
//--------------------------------------------------------
// Configuration for Chamber TC1/2  (Station 1) ----------           
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 Float_t rseg[4]={17.5, 55.2, 71.3, 95.5};
 Int_t   nseg[4]={4, 4, 2, 1};

 chamber=1;
//^^^^^^^^^
 MUON->SetNsec(chamber-1,2);
//
 AliMUONsegmentationV01 *seg11=new AliMUONsegmentationV01;
 seg11->SetSegRadii(rseg);
 seg11->SetPADSIZ(3.048, 0.508);
 seg11->SetPadDivision(nseg);
 MUON->SetSegmentationModel(chamber-1, 1, seg11);
//
 AliMUONsegmentationV01 *seg12=new AliMUONsegmentationV01;
 seg12->SetSegRadii(rseg); 
 seg12->SetPADSIZ(2.032, 0.762);
 seg12->SetPadDivision(nseg);

 MUON->SetSegmentationModel(chamber-1, 2, seg12);

 chamber=2;
//^^^^^^^^^
 MUON->SetNsec(chamber-1,2);
 MUON->SetSegmentationModel(chamber-1, 1, seg11);
 MUON->SetSegmentationModel(chamber-1, 2, seg12);

 station=1;
//^^^^^^^^^ 
 MUON->SetResponseModel(0, response0);	    
 MUON->SetResponseModel(1, response0);	    
//
//--------------------------------------------------------
// Configuration for Chamber TC3/4 -----------------------
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

 chamber=3;
 MUON->SetNsec(chamber-1,1);
 AliMUONsegmentationV0 *seg34=new AliMUONsegmentationV0;
 seg34->SetDAnod(0.51/3.);
 
 MUON->SetSegmentationModel(chamber-1, 1, seg34);
 MUON->SetResponseModel(chamber-1, response0);	    

 chamber=4;
 MUON->SetNsec(chamber-1,1);
 MUON->SetSegmentationModel(chamber-1, 1, seg34);
 MUON->SetResponseModel(chamber-1, response0);	    
//
// Station 2
 station=2;
 MUON->SetPADSIZ(station, 1, 0.75, 0.51);
 MUON->SetMUCHSP(station, 5.);
 MUON->SetMUSIGM(station, 0.18, 0.18);
 MUON->SetRSIGM(station, 10.);
 MUON->SetMAXADC(station, 1024);

//
//--------------------------------------------------------
// Configuration for Chamber TC5/6 -----------------------
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

 seg5 =  new AliMUONsegmentationV1;
 AliMUONresponseV0* response5 =  new AliMUONresponseV0;
 // K3 = 0.62
 response5->SetSqrtKx3(0.78740079);
 response5->SetKx2(0.95237319); //  0.5 * kPI * (1- 0.5*sqrtky3 )
 response5->SetKx4(0.37480633); // 0.25/TMath::ATan(sqrtkx3)
 // K3 = 0.55
 response5->SetSqrtKy3(0.74161985);
 response5->SetKy2(0.98832946);
 response5->SetKy4(0.39177817);
 response5->SetPitch(0.325);
 response5->SetRSIGM(10.);
 response5->SetMUCHSP(5.);
 response5->SetMUSIGM( 0.4, 0.4);
 response5->SetMAXADC( 1024);

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
 MUON->SetPADSIZ(station, 1, 0.975, 0.55);

//
//--------------------------------------------------------
// Configuration for Chamber TC7/8/9/10-------------------
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

 chamber=7;
 MUON->SetNsec(chamber-1,1);
 AliMUONsegmentationV0 *seg78=new AliMUONsegmentationV0;
 seg78->SetDAnod(0.51/3.);

 MUON->SetSegmentationModel(chamber-1, 1, seg78);
 MUON->SetResponseModel(chamber-1, response0);	    

 chamber=8;
 MUON->SetNsec(chamber-1,1);
 MUON->SetSegmentationModel(chamber-1, 1, seg78);
 MUON->SetResponseModel(chamber-1, response0);	    
//
// Station 4
 station=4;
 MUON->SetPADSIZ(station, 1, 0.75, 0.5);

 chamber=9;
 MUON->SetNsec(chamber-1,1);
 AliMUONsegmentationV0 *seg910=new AliMUONsegmentationV0;
 seg910->SetDAnod(0.51/3.);

 MUON->SetSegmentationModel(chamber-1, 1, seg910);
 MUON->SetResponseModel(chamber-1, response0);	    

 chamber=10;
 MUON->SetNsec(chamber-1,1);
 MUON->SetSegmentationModel(chamber-1, 1, seg910);
 MUON->SetResponseModel(chamber-1, response0);	    
//
// Station 5
 station=5;
 MUON->SetPADSIZ(station, 1, 0.75, 0.5);

 chamber=11;
 MUON->SetNsec(chamber-1,1);
 AliMUONsegmentationV0 *seg1112=new AliMUONsegmentationV0;
 seg1112->SetDAnod(0.51/3.);

 MUON->SetSegmentationModel(chamber-1, 1, seg1112);
 MUON->SetResponseModel(chamber-1, response0);	    

 chamber=12;
 MUON->SetNsec(chamber-1,1);
 MUON->SetSegmentationModel(chamber-1, 1, seg1112);
 MUON->SetResponseModel(chamber-1, response0);	    
//
// Trigger Station 1
 station=6;
 MUON->SetPADSIZ(station, 1, 0.75, 0.5);

 chamber=13;
 MUON->SetNsec(chamber-1,1);
 AliMUONsegmentationV0 *seg1314=new AliMUONsegmentationV0;
 seg1314->SetDAnod(0.51/3.);

 MUON->SetSegmentationModel(chamber-1, 1, seg1314);
 MUON->SetResponseModel(chamber-1, response0);	    

 chamber=14;
 MUON->SetNsec(chamber-1,1);
 MUON->SetSegmentationModel(chamber-1, 1, seg1314);
 MUON->SetResponseModel(chamber-1, response0);	    
//
// Trigger Station 2
 station=7;
 MUON->SetPADSIZ(station, 1, 0.75, 0.5);
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
