void Config()
{

new TGeant3("C++ Interface to Geant3");

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
geant3->SetDCAY(0);
geant3->SetPAIR(0);
geant3->SetCOMP(0);
geant3->SetPHOT(0);
geant3->SetPFIS(0);
geant3->SetDRAY(0);
geant3->SetANNI(0);
geant3->SetBREM(0);
geant3->SetMUNU(1);
geant3->SetCKOV(1);
geant3->SetHADR(0); //Select pure GEANH (HADR 1) or GEANH/NUCRIN (HADR 3)
geant3->SetLOSS(1);
geant3->SetMULS(1);
geant3->SetRAYL(1);
geant3->SetAUTO(1); //Select automatic STMIN etc... calc. (AUTO 1) or manual (AUTO 0)
geant3->SetABAN(0); //Restore 3.16 behaviour for abandoned tracks
geant3->SetOPTI(2); //Select optimisation level for GEANT geometry searches (0,1,2)
Float_t cut    = 1.e-3; // 1MeV cut by default
Float_t tofmax = 1.e10;
//             GAM ELEC NHAD CHAD MUON EBREM MUHAB EDEL MUDEL MUPA TOFMAX
geant3->SetCUTS(cut,cut, cut, cut, cut, cut,  cut,  cut, cut,  cut, tofmax);
//
//=======================================================================
// ************* STEERING parameters FOR ALICE SIMULATION **************
// --- Specify event type to be tracked through the ALICE setup
// --- All positions are in cm, angles in degrees, and P and E in GeV
//
//*********************************************
// Example for Fixed Particle Gun             *
//*********************************************
//AliGenFixed *gener = new AliGenFixed(200);
//gener->SetMomentumRange(0,999);
//gener->SetPhiRange(0,0);
//gener->SetThetaRange(5., 5.);
//gener->SetOrigin(0,0,0);          //vertex position
//gener->SetPart(14)                //GEANT particle type

//*********************************************
// Example for Moving Particle Gun            *
//*********************************************
/*
AliGenBox *gener = new AliGenBox(500);
gener->SetMomentumRange(0,10);
gener->SetPhiRange(0,360);
gener->SetThetaRange(2., 10.);
gener->SetOrigin(0,0,0);   
       //vertex position
gener->SetSigma(0,0,5.6);           //Sigma in (X,Y,Z) (cm) on IP position
gener->SetPart(14)                //GEANT particle type
 */
//**************************************
// Example for HIJING Parameterisation *
//**************************************
/*
AliGenHIJINGpara *gener = new AliGenHIJINGpara(84210);
gener->SetMomentumRange(0,999);
gener->SetPhiRange(0,360);
gener->SetThetaRange(0.28,179.72);
gener->SetOrigin(0,0,0);                //vertex position
gener->SetSigma(0,0,0);     //Sigma in (X,Y,Z) (cm) on IP position
/*
//********************************************
// Example for Charm  Production with Pythia *
//********************************************
*/
/*
AliGenPythia *gener = new AliGenPythia(200);
gener->SetMomentumRange(0,999);
gener->SetPhiRange(0,360);
gener->SetThetaRange(0., 180.);
gener->SetYRange(2,5);
gener->SetOrigin(0,0,0);          // vertex position
gener->SetSigma(0,0,5.6);         // Sigma in (X,Y,Z) (cm) on IP position
gener->SetProcess(AliGenPythia::jpsi);       
gener->ForceDecay(AliGenPythia::dimuon);       

//*******************************************************
// Example for J/psi  Production from  Parameterisation *
//*******************************************************
/*
 AliGenParam *gener = new AliGenParam(1000, 443);
 gener->SetMomentumRange(0,999);
 gener->SetPhiRange(0,360);
 gener->SetYRange(2,4);
 gener->SetPtRange(1,10);
 gener->SetOrigin(0,0,0);          //vertex position
 gener->SetSigma(0,0,5.6);         //Sigma in (X,Y,Z) (cm) on IP position

//*******************************************************
// Example for a FLUKA Boundary Source                  *
//*******************************************************
/*
AliGenFLUKAsource *gener = new AliGenFLUKAsource(1000); 
gener->SetPartFlag(9);
gener->SetMomentumRange(0,999);
gener->SetPhiRange(0,360);
gener->SetThetaRange(0., 180.);      
 */
//*******************************************************
// Example for a Cocktail                               *
//*******************************************************

AliGenCocktail *gener = new AliGenCocktail(); 

gener->SetMomentumRange(0,999);
gener->SetPhiRange(0,360);
gener->SetYRange(-4,4);
gener->SetPtRange(0,10);
gener->SetOrigin(0,0,0);          //vertex position
gener->SetSigma(0,0,5.6);         //Sigma in (X,Y,Z) (cm) on IP position
//
 AliGenPythia *jpsi = new AliGenPythia(200);
 
 jpsi->SetProcess(AliGenPythia::jpsi);       
 jpsi->ForceDecay(AliGenPythia::dimuon);

 AliGenPythia *beauty = new AliGenPythia(200);
 beauty->SetProcess(AliGenPythia::beauty_unforced);       
 beauty->ForceDecay(AliGenPythia::semielectronic);

 AliGenPythia *charm = new AliGenPythia(200);
 charm->SetProcess(AliGenPythia::charm_unforced);       
 charm->ForceDecay(AliGenPythia::semimuonic);
 charm->SetPtHard(5,10);

 AliGenParam *jpsi_to_muons = new AliGenParam(100,443);
 jpsi_to_muons->ForceDecay(AliGenParam::dimuon);

 AliGenParam *jpsi_to_electrons = new AliGenParam(100,443);
 jpsi_to_electrons->ForceDecay(AliGenParam::dielectron);

 AliGenParam *phi_to_electrons = new AliGenParam(100,333);
 phi_to_electrons->ForceDecay(AliGenParam::dielectron);

// gener->AddGenerator(jpsi,"Jpsi",1.);
// gener->AddGenerator(beauty,"Beauty",1.);
// gener->AddGenerator(charm,"Charm",1.);
// gener->AddGenerator(jpsi_to_muons,"Jpsi_to_Muons",1.);
 gener->AddGenerator(jpsi_to_electrons,"Jpsi_to_Electrons",1.);
 // gener->AddGenerator(phi_to_electrons,"Phi_to_Electrons",1.);
//
gener->Init();
//**************************************************************************
// Specify maximum magnetic field in Tesla (neg. ==> default field)
gAlice->SetField(-999,2);    
// gAlice->TrackingLimits(2000.,200);

//=================== Alice BODY parameters =============================

AliBODY *BODY = new AliBODY("BODY","Alice envelop");
/*
AliFRAME *FRAME = new AliFRAMEv0("FRAME", "Space Frame");
 
/*
//=================== ABSO parameters ============================

AliABSO *ABSO  = new AliABSO("ABSO","Muon Absorber");

//=================== DIPO parameters ============================

AliDIPO *DIPO  = new AliDIPOv2("DIPO","Dipole version 2");

//=================== SHIL parameters ============================

AliSHIL *SHIL  = new AliSHIL("SHIL","Shielding");

//=================== PIPE parameters ============================
*/
// AliPIPE *PIPE  = new AliPIPEv0("PIPE","Beam Pipe");
/*
*/
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
 AliMUONSegmentationV0* segV0 = new AliMUONSegmentationV0;
// Default response
 AliMUONResponseV0* response0 = new AliMUONResponseV0;
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
 AliMUONSegmentationV01 *seg11=new AliMUONSegmentationV01;
 seg11->SetSegRadii(rseg);
 seg11->SetPADSIZ(3.048, 0.508);
 seg11->SetPadDivision(nseg);
 MUON->SetSegmentationModel(chamber-1, 1, seg11);
//
 AliMUONSegmentationV01 *seg12=new AliMUONSegmentationV01;
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
 AliMUONSegmentationV0 *seg34=new AliMUONSegmentationV0;
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

 seg5 =  new AliMUONSegmentationV1;
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
 AliMUONSegmentationV0 *seg78=new AliMUONSegmentationV0;
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
 AliMUONSegmentationV0 *seg910=new AliMUONSegmentationV0;
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
 AliMUONSegmentationV0 *seg1112=new AliMUONSegmentationV0;
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
 AliMUONSegmentationV0 *seg1314=new AliMUONSegmentationV0;
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





















