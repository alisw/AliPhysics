enum gentype_t {hijing, hijingParam, gun, box, pythia, 
		param1, param2, param3, param4, 
		cocktail, fluka, halo, ntuple, scan, doublescan};

gentype_t gentype=param4;

ntracks=1;

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
 decayer->SetForceDecay(kAll);
 decayer->Init();
 gMC->SetExternalDecayer(decayer);

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
Float_t cut    = 1.e-3; // 1MeV cut by default
Float_t tofmax = 1.e10;
//             GAM ELEC NHAD CHAD MUON EBREM MUHAB EDEL MUDEL MUPA TOFMAX
geant3->SetCUTS(cut,cut, cut, cut, cut, cut,  cut,  cut, cut,  cut, tofmax);
//
//=======================================================================
// ************* STEERING parameters FOR ALICE SIMULATION **************
// --- Specify event type to be tracked through the ALICE setup
// --- All positions are in cm, angles in degrees, and P and E in GeV

 switch(gentype)
 {
 case gun:
//*********************************************
// Example for Fixed Particle Gun             
//*********************************************
     AliGenFixed *gener = new AliGenFixed(ntracks);
     gener->SetMomentum(50);
     gener->SetPhi(180.);
     gener->SetTheta(5.);
     gener->SetOrigin(0,0,0);          //vertex position
     gener->SetPart(13);                //GEANT particle type
     break;
 case box:  
//*********************************************
// Example for Moving Particle Gun            *
//*********************************************
     AliGenBox *gener = new AliGenBox(ntracks);
     gener->SetMomentumRange(3,4);
     gener->SetPhiRange(0,360);
     gener->SetThetaRange(90, 180. );
     gener->SetOrigin(0,0,0);   
     //vertex position
     gener->SetSigma(0,0,0);           //Sigma in (X,Y,Z) (cm) on IP position
     gener->SetPart(5);                //GEANT particle type
     break;
 case scan:  
//*********************************************
// Scanning on a grid                         *
//*********************************************
     AliGenScan *gener = new AliGenScan(-1);
     gener->SetMomentumRange(4,4);
     gener->SetPhiRange(0,360);
     gener->SetThetaRange(9,9);
     //vertex position
     gener->SetSigma(6,6,0);           //Sigma in (X,Y,Z) (cm) on IP position
     gener->SetPart(5); 
     gener->SetRange(20, -100, 100, 20, -100, 100, 1, 500, 500);
     break;
     
 case hijingParam:
     AliGenHIJINGpara *gener = new AliGenHIJINGpara(ntracks);
     gener->SetMomentumRange(0,999);
     gener->SetPhiRange(0,360);
     gener->SetThetaRange(2,10);
     gener->SetOrigin(0,0,0);        //vertex position
     gener->SetSigma(0,0,0);         //Sigma in (X,Y,Z) (cm) on IP position
     break;
 case hijing:
     AliGenHijing *gener = new AliGenHijing(-1);
// centre of mass energy 
     gener->SetEnergyCMS(5500);
// reference frame
     gener->SetReferenceFrame("CMS     ");
// projectile
     gener->SetProjectile("A       ", 208, 82);
     gener->SetTarget    ("A       ", 208, 82);
// impact parameter range
     gener->SetImpactParameterRange(0, 3.);
// evaluate cross section before run
     gener->SetEvaluate(1);
// tell hijing to keep the full parent child chain
     gener->KeepFullEvent();
// enable jet quenching
     gener->SetJetQuenching(1);
// enable shadowing
     gener->SetShadowing(1);
// neutral pion and heavy particle decays switched off
     gener->SetDecaysOff(1);
// trigger
     gener->SetTrigger(0);
// kinematic selection
     gener->SetSelectAll(0);
// momentum range
     gener->SetMomentumRange(0,999);
// phi range
     gener->SetPhiRange(-180,180);
// theta range 
     gener->SetThetaRange(0,180.);
// select flavor (0: no, 4: charm+beauty, 5:beauty)
     gener->SetFlavor(4);
//     
     gener->SetOrigin(0., 0.0 ,0);
     gener->SetSigma(0,0,5.3);
     gener->SetVertexSmear(kPerEvent); 
// no tracking
     gener->SetTrackingFlag(0);
     break;
     
 case pythia:
//********************************************
// Example for Charm  Production with Pythia *
//********************************************
     AliGenPythia *gener = new AliGenPythia(-1);
//   final state kinematic cuts
     gener->SetMomentumRange(0,999);
     gener->SetPhiRange(-180,180);
     gener->SetThetaRange(0., 180.);
     gener->SetYRange(-10,10);
     gener->SetPtRange(0,100);
//   vertex position and smearing 
     gener->SetOrigin(0,0,0);         // vertex position
     gener->SetVertexSmear(kPerEvent);
     gener->SetSigma(0,0,5.6);        // Sigma in (X,Y,Z) (cm) on IP position
//   Structure function
// DO_Set_1
// GRV_LO
// GRV_HO
// MRS_D_minus
// MRS_D0
// MRS_G
// CTEQ_2pM
// CTEQ_4M
     gener->SetStrucFunc(kGRV_HO);
// Select corection for nuclear structure functions
//     gener->SetNuclei(208,208);
//
//   Process type
//   charm, beauty, charm_unforced, beauty_unforced, jpsi, jpsi_chi, mb
     gener->SetProcess(kPyBeauty);
//   
//   Pt transfer of the hard scattering
     gener->SetPtHard(0.,5.);
//   Decay type (semielectronic, semimuonic, nodecay)
     gener->SetForceDecay(kSemiElectronic);
//   Centre of mass energy 
     gener->SetEnergyCMS(5500.);
//   No Tracking 
     gener->SetTrackingFlag(0);
     break;              

 case param1:
//*******************************************************
// Example for J/psi  Production from  Parameterisation 
// using default library (AliMUONlib)                                       
//*******************************************************
     AliGenParam *gener =
	 new AliGenParam(ntracks, AliGenMUONlib::kUpsilon);
     gener->SetMomentumRange(0,999);
     gener->SetPtRange(0,999);     
     gener->SetPhiRange(-180, 180);
     gener->SetYRange(2.5,4);
     gener->SetCutOnChild(1);
     gener->SetChildThetaRange(2,9);
     gener->SetOrigin(0,0,0);          //vertex position
     gener->SetSigma(0,0,5.3);         //Sigma in (X,Y,Z) (cm) on IP position
     gener->SetForceDecay(kDiMuon);
     gener->SetTrackingFlag(0);
     break;

 case param2:
//*******************************************************
// Example for Omega  Production from  Parameterisation 
// specifying library.                                       
//*******************************************************
     AliGenParam *gener = new AliGenParam(1000,new AliGenPHOSlib(), 
					  AliGenPHOSlib::kOmega);
     gener->SetWeighting(kNonAnalog);
     gener->SetForceDecay(kNoDecay);
     gener->SetPtRange(0,100);
     gener->SetThetaRange(45,135);
     gener->SetTrackingFlag(0);
     break;

 case param3:
//*******************************************************
// Example for Upsilon  Production from  Parameterisation 
// specifying library.                                       
// GSI style
//*******************************************************
     AliGenParam *gener = new AliGenParam(1000,new AliGenGSIlib(), 
					  AliGenGSIlib::kUpsilon, "MUON");
     gener->SetMomentumRange(0,999);
     gener->SetPtRange(0,999);     
     gener->SetPhiRange(-180, 180);
     gener->SetYRange(2.5,4);
     gener->SetCutOnChild(1);
     gener->SetChildThetaRange(2,9);
     gener->SetOrigin(0,0,0);          //vertex position
     gener->SetSigma(0,0,5.3);         //Sigma in (X,Y,Z) (cm) on IP position
     gener->SetForceDecay(kDiMuon);
     gener->SetTrackingFlag(0);
     break;
     
 case param4:
//*******************************************************
// Example for Omega  Production from  Parameterisation 
// specifying library.
// The alternative way.                                       
//*******************************************************
     AliGenLib* Lib=new AliGenPHOSlib();
     Int_t iOmega = AliGenPHOSlib::kOmega;
     AliGenParam *gener = new AliGenParam(50, iOmega,            
					  Lib->GetPt(iOmega, ""),
					  Lib->GetY (iOmega, ""),
					  Lib->GetIp(iOmega, ""));
     gener->SetWeighting(kNonAnalog);
     gener->SetForceDecay(kNoDecay);
     gener->SetTrackingFlag(0);
     break;
     
 case fluka:
//*******************************************************
// Example for a FLUKA Boundary Source                  *
//*******************************************************
     AliGenFLUKAsource *gener = new AliGenFLUKAsource(-1);
     gener->SetFileName("$(ALICE_ROOT)/data/all32.root"); 
     gener->SetPartFlag(9);
     gener->SetAgeMax(1.e-5);
//  31.7 events     
     gener->SetFraction(0.0315);     
//     gener->SetFraction(0.75*0.0315);     
     rootfile->cd();
//     gener->SetPartFlag(10);
     gener->SetMomentumRange(0,999);
     gener->SetPhiRange(0,360);
     gener->SetThetaRange(0., 180.); 
     gener->SetAgeMax(1.e-5);
     
//  31.7 events     
//     gener->SetFraction(0.0315);     
     break;

 case ntuple:
//*******************************************************
// Example for reading from a external file                  *
//*******************************************************
     AliGenExtFile *gener = new AliGenExtFile(-1); 
     gener->SetFileName("$(ALICE_ROOT)/data/dtujet93.root");
     gener->SetVertexSmear(kPerEvent); 
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

     gener->SetPhiRange(0,360);
     gener->SetYRange(2.5,4);
     gener->SetThetaRange(2,9);
     gener->SetPtRange(0,10);
     gener->SetOrigin(0,0,0);          //vertex position
     gener->SetSigma(0,0,0);           //Sigma in (X,Y,Z) (cm) on IP position
     gener->SetMomentumRange(0,999);

     AliGenParam *jpsi = new AliGenParam(1,jpsi_p);
     jpsi->SetForceDecay(dimuon);
     jpsi->SetCutOnChild(1);

     
     AliGenFLUKAsource *bg = new AliGenFLUKAsource(-1);
     bg->AddFile("$(ALICE_ROOT)/data/all32.root"); 
     rootfile->cd();
     bg->SetPartFlag(9);
     bg->SetAgeMax(1.e-5);
//  31.7 events     
//     gener->SetFraction(0.0315);     
     bg->SetFraction(0.01*0.0315);     
      
     gener->AddGenerator(jpsi,"J/Psi", 1);
     gener->AddGenerator(bg,"Background",1);

     break;
 }
 
// Activate this line if you want the vertex smearing to happen
// track by track
//
// gener->SetVertexSmear(kPerTrack); 

gener->Init();

gAlice->SetField(-999,2);    //Specify maximum magnetic field in Tesla (neg. ==> default field)

Int_t iMAG=1;

//=================== Alice BODY parameters =============================
AliBODY *BODY = new AliBODY("BODY","Alice envelop");


if(iMAG) {
//=================== MAG parameters ============================
// --- Start with Magnet since detector layouts may be depending ---
// --- on the selected Magnet dimensions ---
AliMAG *MAG  = new AliMAG("MAG","Magnet");
}
}
