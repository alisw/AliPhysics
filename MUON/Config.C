enum gentype_t {hijing, gun, box, pythia, param, cocktail, fluka, halo, ntuple, scan, doublescan, hijing_g};

gentype_t gentype=param;
Int_t ntracks=1;

void Config()

{
 gSystem->Load("libgeant321");
 new     TGeant3("C++ Interface to Geant3");

 

//=======================================================================
//  Create the output file
   
 TFile *rootfile = new TFile("galice.root","recreate");
 rootfile->SetCompressionLevel(2);
 AliDecayer* decayer = new AliDecayerPythia();
 decayer->SetForceDecay(kAll);
 decayer->Init();
 gMC->SetExternalDecayer(decayer);

 TGeant3 *geant3 = (TGeant3*)gMC;
//=======================================================================
// ******* GEANT STEERING parameters FOR ALICE SIMULATION *******
geant3->SetTRIG(1);          //Number of events to be processed 
geant3->SetSWIT(4,100);
geant3->SetDEBU(0,0,1);
geant3->SetDCAY(1);
geant3->SetPAIR(1);
geant3->SetCOMP(1);
geant3->SetPHOT(1);
geant3->SetPFIS(0);
geant3->SetDRAY(0);
geant3->SetANNI(1);
geant3->SetBREM(1);
geant3->SetMUNU(1);
geant3->SetCKOV(0);
geant3->SetHADR(4); //Select pure GEANH (HADR 1) or GEANH/NUCRIN (HADR 3)
geant3->SetLOSS(1);
geant3->SetMULS(1);
geant3->SetRAYL(0);
geant3->SetAUTO(1); //Select automatic STMIN etc... calc. (AUTO 1) or manual (AUTO 0)
geant3->SetABAN(1); //Restore 3.16 behaviour for abandoned tracks
geant3->SetOPTI(2); //Select optimisation level for GEANT geometry searches (0,1,2)
Float_t cut    = 1.e-4; // 100MeV cut by default
Float_t tofmax = 1.e10;
//              GAM    ELEC   NHAD   CHAD   MUON  EBREM  MUHAB EDEL MUDEL MUPA TOFMAX
geant3->SetCUTS(1.e-4, 1.e-4, 1.e-3, 1.e-4, 1.e-3, cut,  cut,  cut, cut,  cut, 1.e-5);

gAlice->TrackingLimits(700, 2000);
 
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
     gener->SetMomentum(10);
     gener->SetPhiRange(0);
     gener->SetThetaRange(0.);
     gener->SetOrigin(30,30,1200);          //vertex position
     gener->SetPart(kProton);          //GEANT particle type
     break;
 case box:  
//*********************************************
// Example for Moving Particle Gun            *
//*********************************************
     AliGenBox *gener = new AliGenBox(ntracks);
     gener->SetMomentumRange(33,34);
     gener->SetPhiRange(-180,180);
     gener->SetThetaRange(2., 9.);
     gener->SetOrigin(0,0,0);   
     gener->SetVertexSmear(kPerTrack); 
     //vertex position
     gener->SetSigma(0, 0, 0);   // Sigma in (X,Y,Z) (cm) on IP position
     gener->SetPart(kMuonPlus);    // GEANT particle type
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
     gener->SetSigma(1,1,0);           //Sigma in (X,Y,Z) (cm) on IP position
     gener->SetPart(kMuonMinus); 
     gener->SetRange(100, -300., 300., 100, -300., 300., 1, 900, 900);
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
     gener->SetPtRange(0,999);
     gener->SetPhiRange(0,360);
//     gener->SetThetaRange(0.104,33.52);
     gener->SetThetaRange(0.104,90.0);
//     gener->SetThetaRange(2.,9.);
     gener->SetOrigin(0., 0.0 ,0);          // vertex position
     gener->SetSigma(0,0,5.3);              // Sigma in (X,Y,Z) (cm) on IP position
     gener->SetVertexSmear(kPerTrack); 
     gener->SetTrackingFlag(0);
     break;
 case hijing_g:
     AliGenHijing *gener = new AliGenHijing(-1);

     gener->SetEnergyCMS(5600.);
     gener->SetReferenceFrame("CMS");
     gener->SetProjectile("A", 208, 82);
     gener->SetTarget    ("A", 208, 82);
     gener->SetImpactParameterRange(0, 5.);
     gener->SetEvaluate(0);
     gener->KeepFullEvent();
     gener->SetJetQuenching(1);
     gener->SetShadowing(1);
     gener->SetDecaysOff(1);
     gener->SetTrigger(0);
     gener->SetSelectAll(1);
     gener->SetMomentumRange(0,9999);
     gener->SetPhiRange(-180,180);
     gener->SetThetaRange(0.104,90.0);
//     gener->SetFlavor(4);
     gener->SetOrigin(0., 0.0 ,0);
     gener->SetSigma(0,0,5.3);
     gener->SetVertexSmear(kPerEvent); 
     gener->SetTrackingFlag(0);
     
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
     //gener->SetOrigin(0,0,0);          // vertex position
     //gener->SetVertexSmear(kPerEvent);
     //gener->SetSigma(0,0,5.6);         // Sigma in (X,Y,Z) (cm) on IP position
     gener->SetStrucFunc(kDOSet1);
     gener->SetProcess(kPyCharm);
     gener->SetEnergyCMS(5500.);
     break;              
/*
     AliGenPythia *gener = new AliGenPythia(ntracks);
     gener->SetMomentumRange(0,999);
     gener->SetPhiRange(0,360);
     gener->SetThetaRange(0., 180.);
     gener->SetYRange(-10,10);
     gener->SetPtRange(0,100);
     gener->SetOrigin(0,0,0);          // vertex position
     gener->SetVertexSmear(kPerEvent); 
     gener->SetSigma(0,0,5.6);         // Sigma in (X,Y,Z) (cm) on IP position
     gener->SetStrucFunc(DO_Set_1);   
     gener->SetProcess(charm);     
     gener->SetForceDecay(dimuon);
     gener->SetEnergyCMS(5500.);
     gener->SetTrackingFlag(0);
     
     break;
 */   
 case param:
//*******************************************************
// Example for J/psi or Upsilon Production from  Parameterisation *
//*******************************************************
     AliGenParam *gener = new AliGenParam(ntracks, AliGenMUONlib::kUpsilon);
     gener->SetMomentumRange(0,999);
     gener->SetPtRange(0,100.);
     gener->SetPhiRange(-180, 180);
     gener->SetYRange(2.5,4);
     gener->SetCutOnChild(1);
     gener->SetChildThetaRange(2.0,9);
     gener->SetOrigin(0,0,0);          //vertex position
     gener->SetSigma(0,0,0);           //Sigma in (X,Y,Z) (cm) on IP position
     gener->SetForceDecay(kDiMuon);
     gener->SetTrackingFlag(1);

     break;
     
 case fluka:
//*******************************************************
// Example for a FLUKA Boundary Source                  *
//*******************************************************
     AliGenFLUKAsource *gener = new AliGenFLUKAsource(-1);
     gener->AddFile("$(ALICE_ROOT)/data/alice.root"); 
     rootfile->cd();
     gener->SetPartFlag(7);
     gener->SetMomentumRange(0,999);
     gener->SetPhiRange(0,360);
     gener->SetThetaRange(0., 180.); 
     gener->SetAgeMax(1.e-5);
     
//  31.7 events     
     gener->SetFraction(1.);     
     break;

 case ntuple:
//*******************************************************
// Example for reading from a external file                  *
//*******************************************************
     AliGenExtFileCH *gener = new AliGenExtFileCH(-1); 
     gener->SetFileName("$(ALICE_ROOT)/data/pbpb.root");
     gener->SetThetaRange(0.104,90.);
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
     gener->SetThetaRange(2.,9.);
     gener->SetTrackingFlag(0);
     AliGenParam *Pi0 = new AliGenParam(100, new AliGenPMDlib(), AliGenPMDlib::kPion);
     AliGenParam *Eta = new AliGenParam( 10, new AliGenPMDlib(), AliGenPMDlib::kEta);
     gener->AddGenerator(Pi0, "neutral pions"  , 1.);
     gener->AddGenerator(Eta, "neutral etas"  ,  1.);
     break;
 }
 

gener->Init();
gAlice->SetField(2,1);    //Specify maximum magnetic field in Tesla (neg. ==> default field)

Int_t iFRAME  =0;
Int_t iMAG    =0;
Int_t iITS    =0;
Int_t iABSO   =1;
Int_t iDIPO   =0;
Int_t iHALL   =0;
Int_t iSHIL   =1;
Int_t iPIPE   =1;
Int_t iFMD    =0;
Int_t iMUON   =1;

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



if(iSHIL) {
//=================== SHIL parameters ============================
//    AliSHIL *SHIL  = new AliSAROV("SHIL","Shielding");
    AliSHILvF *SHIL  = new AliSHILvF("SHIL","Shielding");
    SHIL->SetPbCone(1);
//    AliSAROV *SHIL  = new AliSAROV("SHIL","Shielding");
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

    AliMUON *MUON  = new AliMUONv1("MUON","default");
    
}

}

		



