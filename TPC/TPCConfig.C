void Config()
{
//=======================================================================
//  Create the output file
   
TFile *rootfile = new TFile("galice.root","recreate");
rootfile->SetCompressionLevel(2);

//=======================================================================
// ******* GEANT STEERING parameters FOR ALICE SIMULATION *******
geant3->SetTRIG(1); //Number of events to be processed 
geant3->SetSWIT(4,1);
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
AliGenHIJINGpara *gener = new AliGenHIJINGpara(1000);
gener->SetMomentumRange(0,100);
gener->SetPhiRange(0,360);
gener->SetThetaRange(20,160);
gener->SetOrigin(0,0,0);        //vertex position
gener->SetSigma(0,0,0);         //Sigma in (X,Y,Z) (cm) on IP position
gener->Init();

gAlice->SetField(-999,2);    //Specify maximum magnetic field in Tesla (neg. ==> default field)

//=================== Alice BODY parameters =============================
AliBODY *BODY = new AliBODY("BODY","Alice envelop");



//============================ TPC parameters ================================
// --- This allows the user to specify sectors for the SLOW (TPC geometry 2)
// --- Simulator. TPCSECL (TPSECU)  <0 means that ALL lower (upper)
// --- sectors are specified, any value other than that requires at least one 
// --- sector (lower or upper)to be specified!
// --- Reminder: sectors 1-24 are lower sectors (1-12 -> z>0, 13-24 -> z<0)
// ---           sectors 25-72 are the upper ones (25-48 -> z>0, 49-72 -> z<0)
// --- TPCLOWS - number of lower sectors specified (up to 6)
// --- TPCUPS - number of upper sectors specified (up to 12)
// --- TPCSENS - sensitive strips for the Slow Simulator !!!
// --- This does NOT work if all S or L-sectors are specified, i.e.
// --- if TPCSECL or TPCSECU < 0
//-----------------------------------------------------------------------------

AliTPC *TPC  = new AliTPCv1("TPC","Normal TPC");
TPC->SetSecAL(1);
TPC->SetSecAU(1);
TPC->SetSecLows(1, -1, -1, -1, -1, -1);
TPC->SetSecUps(25, 26, 48, -1, -1, -1, -1, -1, -1, -1, -1, -1);
TPC->SetSens(1);



         
}
