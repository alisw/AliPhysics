enum gentype_t {hijing, gun, box, pythia, param, cocktail, fluka, halo, ntuple, scan, doublescan};

gentype_t gentype=param;
Int_t ntracks=1000;

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
     gener->SetVertexSmear(perTrack); 
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
     AliGenParam *gener = new AliGenParam(ntracks,upsilon_p);
//     AliGenParam *gener = new AliGenParam(ntracks, jpsi_p);
     gener->SetMomentumRange(0,999);
     gener->SetPhiRange(0,360);
     gener->SetYRange(2.5,4);
     gener->SetPtRange(0,50);
     gener->SetOrigin(0,0,0);          //vertex position
     gener->SetSigma(0,0,0);           //Sigma in (X,Y,Z) (cm) on IP position
     gener->SetForceDecay(dimuon);
     gener->SetCutOnChild(1);
     gener->SetChildThetaRange(2,9);
     gener->SetTrackingFlag(1);

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
     gener->SetVertexSmear(perTrack); 
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
Int_t iMAG    =0;
Int_t iABSO   =0;
Int_t iDIPO   =0;
Int_t iSHIL   =0;
Int_t iPIPE   =0;
Int_t iMUON   =1;

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
    AliABSO *ABSO  = new AliABSOv1("ABSO","Muon Absorber");
}

if(iDIPO) {
//=================== DIPO parameters ============================

    AliDIPO *DIPO  = new AliDIPOv2("DIPO","Dipole version 2");
}


if(iSHIL) {
//=================== SHIL parameters ============================
    AliSHIL *SHIL  = new AliSHILv1("SHIL","Shielding");
}


if(iPIPE) {
//=================== PIPE parameters ============================
    AliPIPE *PIPE  = new AliPIPEv0("PIPE","Beam Pipe");
}


if(iMUON) {
//=================== MUON parameters ===========================

AliMUON *MUON  = new AliMUONv0("MUON","normal MUON");
 MUON->SetAcceptance(kTRUE, 0.1, 12.);
 MUON->SetIshunt(0);
 MUON->SetMaxStepGas(0.1);
 MUON->SetMaxStepAlu(0.1);
}
}
