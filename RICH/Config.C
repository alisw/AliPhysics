enum gentype_t {hijing, gun, box, pythia, param, cocktail, fluka, halo, ntuple, scan, doublescan};

gentype_t gentype=gun;
ntracks=1;

void Config()
{

new AliGeant3("C++ Interface to Geant3");

//=======================================================================
//  Create the output file
   
TFile *rootfile = new TFile("galice.root","recreate");
rootfile->SetCompressionLevel(2);
TGeant3 *geant3 = (TGeant3*)gMC;

//geant3->Grndmq(0,0,55," ");
 
//=======================================================================
// ******* GEANT STEERING parameters FOR ALICE SIMULATION *******
geant3->SetTRIG(1); //Number of events to be processed 
geant3->SetSWIT(4,100);
geant3->SetDEBU(0,0,1);
//geant3->SetSWIT(2,2);
geant3->SetDCAY(1);
geant3->SetPAIR(1);
geant3->SetCOMP(1);
geant3->SetPHOT(1);
geant3->SetPFIS(0);
geant3->SetDRAY(1);
geant3->SetANNI(1);
geant3->SetBREM(1);
geant3->SetMUNU(1); 
geant3->SetCKOV(1);
geant3->SetHADR(3); //Select pure GEANH (HADR 1) or GEANH/NUCRIN (HADR 3)
geant3->SetLOSS(1);
geant3->SetMULS(1);
geant3->SetRAYL(0);
geant3->SetAUTO(1); //Select automatic STMIN etc... calc. (AUTO 1) or manual (AUTO 0)
geant3->SetABAN(0); //Restore 3.16 behaviour for abandoned tracks
geant3->SetOPTI(2); //Select optimisation level for GEANT geometry searches (0,1,2)
Float_t cut    = 1.e-1; // 100MeV cut by default
Float_t tofmax = 1.e10;
//             GAM ELEC NHAD CHAD MUON EBREM MUHAB EDEL MUDEL MUPA TOFMAX
geant3->SetCUTS(1.e-5,5.e-5, 1.e-3, 1.e-4, cut, cut,  cut,  cut, cut,  cut, tofmax);

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
     gener->SetPhiRange(92);
     gener->SetThetaRange(90);
     gener->SetOrigin(0,0,0);                 //vertex position
     gener->SetPart(kPiPlus);                 //GEANT particle type
     break;
 case box:  
//*********************************************
// Example for Moving Particle Gun            *
//*********************************************
     AliGenBox *gener = new AliGenBox(ntracks);
     gener->SetMomentumRange(3,3);
     gener->SetPhiRange(80,100);
     gener->SetThetaRange(80,100);
     gener->SetOrigin(0,0,0);   
     gener->SetVertexSmear(kPerTrack); 
     //vertex position
     gener->SetSigma(1.8, 1.8,0);           //Sigma in (X,Y,Z) (cm) on IP position
     gener->SetPart(kPiPlus);                    //GEANT particle type
     break;
 case scan:  
//*********************************************
// Scanning on a grid                         *
//*********************************************
     AliGenScan *gener = new AliGenScan(-1);
     gener->SetMomentumRange(3,3);
     gener->SetPhiRange(90,90);
     gener->SetThetaRange(90,90);
     //vertex position
     gener->SetSigma(0,0,0);           //Sigma in (X,Y,Z) (cm) on IP position
     gener->SetPart(kPiPlus); 
     gener->SetRange(5, -80, 60, 1, 480, 480, 5, -80, 60);
     break;
 case doublescan:  
//*********************************************
// Scanning on a grid                         *
//*********************************************
     AliGenDoubleScan *gener = new AliGenDoubleScan(-1);
     gener->SetMomentumRange(3,3);
     gener->SetPhiRange(0,360);
     gener->SetThetaRange(0,0);
     //vertex position
     gener->SetSigma(0,0,0);           //Sigma in (X,Y,Z) (cm) on IP position
     gener->SetPart(kPiPlus); 
     gener->SetRange(20, -60, 60, 1, 480, 480, 20, -60, 60);
     gener->SetDistance(1);
     
     break;
     
 case hijing:
     AliGenHIJINGpara *gener = new AliGenHIJINGpara(ntracks);
     gener->SetMomentumRange(0,999);
     gener->SetPhiRange(0,360);
     gener->SetThetaRange(.77,179.23);
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
     gener->SetVertexSmear(perEvent); 
     gener->SetSigma(0,0,5.6);         // Sigma in (X,Y,Z) (cm) on IP position
//     gener->SetStrucFunc(DO_Set_1);
     gener->SetProcess(mb); 
     gener->SetEnergyCMS(5500.);
     break;
     
 case param:
//*******************************************************
// Example for J/psi  Production from  Parameterisation *
//*******************************************************
     AliGenParam *gener = new AliGenParam(178,Eta,
					     AliGenPHOSlib::GetPt(Eta),
					     AliGenPHOSlib::GetY(Eta),
					     AliGenPHOSlib::GetIp(Eta) );

     gener->SetMomentumRange(0,999);
     gener->SetPhiRange(0,360);
     gener->SetYRange(2.5,4);
     gener->SetThetaRange(2,9);
     gener->SetPtRange(0,10);
     gener->SetOrigin(0,0,0);      //vertex position
     gener->SetSigma(0,0,0);       //Sigma in (X,Y,Z) (cm) on IP position
     gener->SetCutOnChild(1);
     break;
     
 case fluka:
//*******************************************************
// Example for a FLUKA Boundary Source                  *
//*******************************************************
     AliGenFLUKAsource *gener = new AliGenFLUKAsource(-1);
     gener->AddFile("$(ALICE_ROOT)/data/all32.root"); 
     rootfile->cd();
     gener->SetPartFlag(9);
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
     gener->SetFileName("$(ALICE_ROOT)/data/dtujet93.root");
     gener->SetVertexSmear(perEvent); 
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
//
//   test 
//
     
     Float_t   p2(Float_t);
     Float_t  (*f1)(Float_t);
     Double_t (*f2)(Double_t);

     
     
     
     Float_t p2(Float_t x) 
	 {
	     return x*x;
	 }
     f1=p2;
     Float_t x = TMath::Sqrt(2);
     
     f1=TMath::Sqrt;
     
     printf("\n Result %f %f \n", (*f1)(2.), TMath::Sqrt(2));
     
	 
     break;
 }
 
// Activate this line if you want the vertex smearing to happen
// track by track
//
gener->SetVertexSmear(kPerTrack); 
gener->Init();
gAlice->SetField(0,2);    //Specify maximum magnetic field in Tesla (neg. ==> default field)

Int_t iMAG=0;
Int_t iITS=0;
Int_t iTPC=0;
Int_t iTOF=0;
Int_t iRICH=1;
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
Int_t iSTART=0; 

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
}


if(iITS) {
//=================== ITS parameters ============================
//
// EUCLID is a flag to output (=1) both geometry and media to two ASCII files 
// (called by default ITSgeometry.euc and ITSgeometry.tme) in a format
// understandable to the CAD system EUCLID. The default (=0) means that you 
// dont want to use this facility.
//
AliITS *ITS  = new AliITSv1("ITS","normal ITS");
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
  AliTPC *TPC  = new AliTPCv0("TPC","Normal TPC"); //v1 is default
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
AliTOF *TOF  = new AliTOFv0("TOF","normal TOF");
}

if(iRICH) {
//=================== RICH parameters ===========================
    AliRICH *RICH  = new AliRICHv2("RICH","normal RICH");
    
//
// Version 0
// Default Segmentation
    AliRICHSegmentationV1* SegmentationV0 = new AliRICHSegmentationV1;
//
//  Segmentation parameters
    SegmentationV0->SetPadSize(0.8,0.84);
    SegmentationV0->SetDAnod(0.84/2);

//  Geometry parameters
    AliRICHGeometry* GeometryV0 = new AliRICHGeometry;
    GeometryV0->SetGapThickness(8);
    GeometryV0->SetProximityGapThickness(.4);
    GeometryV0->SetQuartzLength(133);
    GeometryV0->SetQuartzWidth(127.9);
    GeometryV0->SetQuartzThickness(.5);
    GeometryV0->SetOuterFreonLength(133);
    GeometryV0->SetOuterFreonWidth(41.3);
    GeometryV0->SetInnerFreonLength(133);
    GeometryV0->SetInnerFreonWidth(41.3);
    GeometryV0->SetFreonThickness(1.5);

//  Response parameters
    AliRICHResponseV0*  Rresponse0   = new AliRICHResponseV0;
    Rresponse0->SetSigmaIntegration(5.);
    Rresponse0->SetChargeSlope(27.);
    Rresponse0->SetChargeSpread(0.18, 0.18);
    Rresponse0->SetMaxAdc(4096);
    Rresponse0->SetAlphaFeedback(0.036);
    Rresponse0->SetEIonisation(26.e-9);
    Rresponse0->SetSqrtKx3(0.77459667);
    Rresponse0->SetKx2(0.962);
    Rresponse0->SetKx4(0.379);
    Rresponse0->SetSqrtKy3(0.77459667);
    Rresponse0->SetKy2(0.962);
    Rresponse0->SetKy4(0.379);
    Rresponse0->SetPitch(0.25);

      
  for (Int_t i=0; i<7; i++) {
    RICH->SetGeometryModel(i,GeometryV0);
    RICH->SetSegmentationModel(i, SegmentationV0);
    RICH->SetResponseModel(i, Rresponse0);
    RICH->SetNsec(i,1);
  }  
  RICH->SetDebugLevel(0);
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

AliTRD *TRD  = new AliTRDv1("TRD","TRD version 0");
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
