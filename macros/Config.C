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
geant3->SetTRIG(1); //Number of events to be processed 
geant3->SetSWIT(4,10);
geant3->SetDEBU(0,0,1);
//geant3->SetSWIT(2,2);
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

Float_t cut    = 1.e-3; // 1MeV cut by default
Float_t tofmax = 1.e10;
//             GAM ELEC NHAD CHAD MUON EBREM MUHAB EDEL MUDEL MUPA TOFMAX
geant3->SetCUTS(cut,cut, cut, cut, cut, cut,  cut,  cut, cut,  cut, tofmax);
//
//=======================================================================
// ************* STEERING parameters FOR ALICE SIMULATION **************
// --- Specify event type to be tracked through the ALICE setup
// --- All positions are in cm, angles in degrees, and P and E in GeV
AliGenHIJINGpara *gener = new AliGenHIJINGpara(50);
gener->SetMomentumRange(0,999);
gener->SetPhiRange(0,360);
gener->SetThetaRange(10,170);
gener->SetOrigin(0,0,0);        //vertex position
gener->SetSigma(0,0,0);         //Sigma in (X,Y,Z) (cm) on IP position
gener->Init();
// 
// Activate this line if you want the vertex smearing to happen
// track by track
//
//gener->SetVertexSmear(perTrack); 

gAlice->SetField(-999,2);    //Specify maximum magnetic field in Tesla (neg. ==> default field)

Int_t iMAG=1;
Int_t iITS=1;
Int_t iTPC=1;
Int_t iTOF=1;
Int_t iRICH=1;
Int_t iZDC=0;
Int_t iCASTOR=1;
Int_t iTRD=1;
Int_t iABSO=1;
Int_t iDIPO=1;
Int_t iHALL=1;
Int_t iFRAME=1;
Int_t iSHIL=1;
Int_t iPIPE=1;
Int_t iFMD=0;
Int_t iMUON=1;
Int_t iPHOS=1;
Int_t iPMD=1;
Int_t iSTART=1;

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

AliFRAME *FRAME  = new AliFRAMEv0("FRAME","Space Frame");

}

if(iSHIL) {
//=================== SHIL parameters ============================

AliSHIL *SHIL  = new AliSHILv0("SHIL","Shielding");
}


if(iPIPE) {
//=================== PIPE parameters ============================

AliPIPE *PIPE  = new AliPIPEv0("PIPE","Beam Pipe");
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

  gROOT->LoadMacro("SetTPCParam.C");
  AliTPCParam *param = SetTPCParam();
  AliTPC *TPC  = new AliTPCv1("TPC","Normal TPC"); //v1 is default
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

if(iTOF) {
//=================== TOF parameters ============================
AliTOF *TOF  = new AliTOFv1("TOF","normal TOF");
}

if(iRICH) {
//=================== RICH parameters ===========================
    AliRICH *RICH  = new AliRICHv0("RICH","normal RICH");
    
//
// Version 0
// Default Segmentation
    AliRICHSegmentationV0* SegmentationV0 = new AliRICHSegmentationV0;
//
//  Segmentation parameters
    SegmentationV0->SetPadSize(0.84,0.80);
    SegmentationV0->SetDAnod(0.84/2);
//
//  Geometry parameters
    AliRICHGeometry* GeometryV0 = new AliRICHGeometryV0;
    GeometryV0->SetGapThickness(7.6);
    GeometryV0->SetProximityGapThickness(.4);
    GeometryV0->SetQuartzLength(131);
    GeometryV0->SetQuartzWidth(126.2);
    GeometryV0->SetQuartzThickness(.5);
    GeometryV0->SetOuterFreonLength(131);
    GeometryV0->SetOuterFreonWidth(40.3);
    GeometryV0->SetInnerFreonLength(131);
    GeometryV0->SetInnerFreonWidth(40.3);
    GeometryV0->SetFreonThickness(1);
//
//  Response parameters
    AliRICHResponseV0*  Rresponse0   = new AliRICHResponseV0;
    Rresponse0->SetSigmaIntegration(5.);
    Rresponse0->SetChargeSlope(41.);
    Rresponse0->SetChargeSpread(0.18, 0.18);
    Rresponse0->SetMaxAdc(1024);
    Rresponse0->SetAlphaFeedback(0.05);
    Rresponse0->SetEIonisation(26.e-9);
    Rresponse0->SetSqrtKx3(0.77459667);
    Rresponse0->SetKx2(0.962);
    Rresponse0->SetKx4(0.379);
    Rresponse0->SetSqrtKy3(0.77459667);
    Rresponse0->SetKy2(0.962);
    Rresponse0->SetKy4(0.379);
    Rresponse0->SetPitch(0.25);
//
//      
  for (Int_t i=0; i<7; i++) {
    RICH->SetGeometryModel(i,GeometryV0);
    RICH->SetSegmentationModel(i, SegmentationV0);
    RICH->SetResponseModel(i, Rresponse0);
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
  
  AliTRD *TRD  = new AliTRDv0("TRD","TRD fast simulator");
  //TRD->SetHits();
  
  //AliTRD *TRD  = new AliTRDv1("TRD","TRD slow simulator");
  //TRD->SetSensPlane(0);
  //TRD->SetSensChamber(2);
  //TRD->SetSensSector(17);
  
  // Select the gas mixture (0: 97% Xe + 3% isobutane, 1: 90% Xe + 10% CO2)
  TRD->SetGasMix(1);
  
  // With hole in front of PHOS
  TRD->SetPHOShole();
  // With hole in front of RICH
  TRD->SetRICHhole();
}

if(iFMD) {
//=================== FMD parameters ============================

AliFMD *FMD  = new AliFMDv1("FMD","normal FMD");
}

if(iMUON) {
//=================== MUON parameters ===========================

AliMUON *MUON  = new AliMUONv0("MUON","normal MUON");

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
 AliMUONresponseV0* response0 = new AliMUONresponseV0;
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
//--------------------------------------------------------
// Configuration for Chamber TC1/2  (Station 1) ----------           
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// Float_t rseg1[4]={17.5, 55.2, 71.3, 95.5};
 Float_t rseg1[4]={15.5, 55.2, 71.3, 95.5};
 Int_t   nseg1[4]={4, 4, 2, 1};
//
 chamber=1;
//^^^^^^^^^
 MUON->SetNsec(chamber-1,2);
//
 AliMUONsegmentationV01 *seg11=new AliMUONsegmentationV01;
 
 seg11->SetSegRadii(rseg1);
 seg11->SetPADSIZ(3, 0.5);
 seg11->SetDAnod(3.0/3./4);
 seg11->SetPadDivision(nseg1);
 
 MUON->SetSegmentationModel(chamber-1, 1, seg11);
//
 AliMUONsegmentationV02 *seg12=new AliMUONsegmentationV02;
 seg12->SetSegRadii(rseg1); 
 seg12->SetPADSIZ(0.75, 2.0);
 seg12->SetDAnod(3.0/3./4);
 seg12->SetPadDivision(nseg1);

 MUON->SetSegmentationModel(chamber-1, 2, seg12);

 MUON->SetResponseModel(chamber-1, response0);	    

 chamber=2;
//^^^^^^^^^
//
 MUON->SetNsec(chamber-1,2);
//
 AliMUONsegmentationV01 *seg21=new AliMUONsegmentationV01;
 seg21->SetSegRadii(rseg1);
 seg21->SetPADSIZ(3, 0.5);
 seg21->SetDAnod(3.0/3./4);
 seg21->SetPadDivision(nseg1);
 MUON->SetSegmentationModel(chamber-1, 1, seg21);
//
 AliMUONsegmentationV02 *seg22=new AliMUONsegmentationV02;
 seg22->SetSegRadii(rseg1); 
 seg22->SetPADSIZ(0.75, 2.);
 seg22->SetDAnod(3.0/3./4);
 seg22->SetPadDivision(nseg1);
 MUON->SetSegmentationModel(chamber-1, 2, seg22);

 MUON->SetResponseModel(chamber-1, response0);	    
//
//--------------------------------------------------------
// Configuration for Chamber TC3/4 -----------------------
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// Float_t rseg2[4]={23.5, 47.1, 87.7, 122.5};
 Float_t rseg2[4]={21.5, 47.1, 87.7, 122.5};
 Int_t   nseg2[4]={4, 4, 2, 1};
//
 chamber=3;
//^^^^^^^^^
 MUON->SetNsec(chamber-1,2);
//
 AliMUONsegmentationV01 *seg31=new AliMUONsegmentationV01;
 seg31->SetSegRadii(rseg2);
 seg31->SetPADSIZ(3, 0.5);
 seg31->SetDAnod(3.0/3./4);
 seg31->SetPadDivision(nseg2);
 MUON->SetSegmentationModel(chamber-1, 1, seg31);
//
 AliMUONsegmentationV02 *seg32=new AliMUONsegmentationV02;
 seg32->SetSegRadii(rseg2); 
 seg32->SetPADSIZ(0.75, 2.);
 seg32->SetPadDivision(nseg2);
 seg32->SetDAnod(3.0/3./4);

 MUON->SetSegmentationModel(chamber-1, 2, seg32);

 MUON->SetResponseModel(chamber-1, response0);	    

 chamber=4;
//^^^^^^^^^
//
 MUON->SetNsec(chamber-1,2);
//
 AliMUONsegmentationV01 *seg41=new AliMUONsegmentationV01;
 seg41->SetSegRadii(rseg2);
 seg41->SetPADSIZ(3, 0.5);
 seg41->SetDAnod(3.0/3./4);
 seg41->SetPadDivision(nseg2);
 MUON->SetSegmentationModel(chamber-1, 1, seg41);
//
 AliMUONsegmentationV02 *seg42=new AliMUONsegmentationV02;
 seg42->SetSegRadii(rseg2); 
 seg42->SetPADSIZ(0.75, 2.);
 seg42->SetPadDivision(nseg2);
 seg42->SetDAnod(3.0/3./4);

 MUON->SetSegmentationModel(chamber-1, 2, seg42);

 MUON->SetResponseModel(chamber-1, response0);	    


//--------------------------------------------------------
// Configuration for Chamber TC5/6 -----------------------
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
/*
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
 response5->SetSigmaIntegration(10.);
 response5->SetChargeSlope(50);
 response5->SetChargeSpread(0.4, 0.4);
 response5->SetMaxAdc(4096);

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
*/

 chamber=5;
//^^^^^^^^^
 MUON->SetNsec(chamber-1,2);
//
 AliMUONsegmentationV0 *seg51=new AliMUONsegmentationV0;
 seg51->SetPADSIZ(0.75, 0.5);
 seg51->SetDAnod(3.0/3./4);
 MUON->SetSegmentationModel(chamber-1, 1, seg51);
//
 AliMUONsegmentationV0 *seg52=new AliMUONsegmentationV0;
 seg52->SetPADSIZ(0.5,0.75);
 seg52->SetDAnod(3.0/3./4);
 MUON->SetSegmentationModel(chamber-1, 2, seg52);

 MUON->SetResponseModel(chamber-1, response0);	    

 chamber=6;
//^^^^^^^^^
 MUON->SetNsec(chamber-1,2);
//
 AliMUONsegmentationV0 *seg61=new AliMUONsegmentationV0;
 seg61->SetPADSIZ(0.75, 0.5);
 seg61->SetDAnod(3.0/3./4);
 MUON->SetSegmentationModel(chamber-1, 1, seg61);
//
 AliMUONsegmentationV0 *seg62=new AliMUONsegmentationV0;
 seg62->SetPADSIZ(0.5,0.75);
 seg62->SetDAnod(3.0/3./4);
 MUON->SetSegmentationModel(chamber-1, 2, seg62);

 MUON->SetResponseModel(chamber-1, response0);	  

//--------------------------------------------------------
// Configuration for Chamber TC7/8  (Station 4) ----------           
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

 Int_t   nseg4[4]={4, 4, 2, 1};

 chamber=7;
//^^^^^^^^^
 MUON->SetNsec(chamber-1,2);
//
 AliMUONsegmentationV04 *seg71=new AliMUONsegmentationV04;
 seg71->SetPADSIZ(10.,0.5);
 seg71->SetDAnod(0.25);
 seg71->SetPadDivision(nseg4);
 MUON->SetSegmentationModel(chamber-1, 1, seg71);

 AliMUONsegmentationV05 *seg72=new AliMUONsegmentationV05;
 seg72->SetPADSIZ(1,10);
 seg72->SetDAnod(0.25);
 seg72->SetPadDivision(nseg4);
 MUON->SetSegmentationModel(chamber-1, 2, seg72);

 MUON->SetResponseModel(chamber-1, response0);	    

 chamber=8;
//^^^^^^^^^
 MUON->SetNsec(chamber-1,2);
 AliMUONsegmentationV04 *seg81=new AliMUONsegmentationV04;
 seg81->SetPADSIZ(10., 0.5);
 seg81->SetPadDivision(nseg4);
 seg81->SetDAnod(0.25);
 MUON->SetSegmentationModel(chamber-1, 1, seg81);

 AliMUONsegmentationV05 *seg82=new AliMUONsegmentationV05;
 seg82->SetPADSIZ(1, 10);
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
 AliMUONsegmentationV04 *seg91=new AliMUONsegmentationV04;
 seg91->SetPADSIZ(10.,0.5);
 seg91->SetDAnod(0.25);
 seg91->SetPadDivision(nseg4);
 MUON->SetSegmentationModel(chamber-1, 1, seg91);

 AliMUONsegmentationV05 *seg92=new AliMUONsegmentationV05;
 seg92->SetPADSIZ(1,10);
 seg92->SetDAnod(0.25);
 seg92->SetPadDivision(nseg4);

 MUON->SetSegmentationModel(chamber-1, 2, seg92);

 MUON->SetResponseModel(chamber-1, response0);	    

 chamber=10;
//^^^^^^^^^
 MUON->SetNsec(chamber-1,2);
 AliMUONsegmentationV04 *seg101=new AliMUONsegmentationV04;
 seg101->SetPADSIZ(10., 0.5);
 seg101->SetPadDivision(nseg4);
 seg101->SetDAnod(0.25);
 MUON->SetSegmentationModel(chamber-1, 1, seg101);

 AliMUONsegmentationV05 *seg102=new AliMUONsegmentationV05;
 seg102->SetPADSIZ(1,10);
 seg102->SetPadDivision(nseg4);
 seg102->SetDAnod(0.25);
 MUON->SetSegmentationModel(chamber-1, 2, seg102);

 MUON->SetResponseModel(chamber-1, response0);	    
//--------------------------------------------------------
// Configuration for Trigger staions --------------------- 
// (not yet used/implemented) ----------------------------          
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
 
//=================== PHOS parameters ===========================

if(iPHOS) {
  AliPHOS *PHOS  = new AliPHOSv0("PHOS","GPS2");
}


if(iPMD) {
//=================== PMD parameters ============================

AliPMD *PMD  = new AliPMDv0("PMD","normal PMD");
PMD->SetPAR(1., 1., 0.8, 0.02);
PMD->SetIN(6., 18., -580., 27., 27.);
PMD->SetGEO(0.0, 0.2, 4.);
PMD->SetPadSize(0.8, 1.0, 1.0, 1.5);

}

if(iSTART) {
//=================== START parameters ============================
AliSTART *START  = new AliSTARTv0("START","START Detector");
}

         
}
