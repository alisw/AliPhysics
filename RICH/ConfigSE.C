enum EGeneratorTypes {kHijing, kGun, kBox, kPythia, kParam, kCcocktail, kFluka, kHalo, kNtuple, kScan,
                      kDoubleScan};

EGeneratorTypes gentype=kScan;
Int_t ntracks=1;

void Config()
{
   cout<<"RICH private Config.C> Start\n";
new AliGeant3("C++ Interface to Geant3");

//=======================================================================
//  Create the output file
   
TFile *rootfile = new TFile("galice.root","recreate");
rootfile->SetCompressionLevel(2);
TGeant3 *geant3 = (TGeant3*)gMC;

//
// Set Random Number seed
 gRandom->SetSeed(10);
 
//=======================================================================
// ******* GEANT STEERING parameters FOR ALICE SIMULATION *******
geant3->SetTRIG(1); //Number of events to be processed 
geant3->SetSWIT(4,100);
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
geant3->SetMUNU(0); 
geant3->SetCKOV(1);
geant3->SetHADR(0); //Select pure GEANH (HADR 1) or GEANH/NUCRIN (HADR 3)
geant3->SetLOSS(1);
geant3->SetMULS(0);
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

   switch(gentype){
   case kGun:
     AliGenFixed *gener = new AliGenFixed(ntracks);
     gener->SetMomentum(3);
     gener->SetPhiRange(90);
     gener->SetThetaRange(101);
     gener->SetOrigin(0,480,-20);                 //vertex position
     gener->SetPart(kPiPlus);                 //GEANT particle type
   break;
   case kBox:  
     AliGenBox *gener = new AliGenBox(ntracks);
     gener->SetMomentumRange(2,2);
     //gener->SetPhiRange(30,30);                //for inclined HMPID
     gener->SetPhiRange(82,98);
     gener->SetThetaRange(82,98);
     gener->SetOrigin(0,0,0);   
     gener->SetVertexSmear(kPerTrack); 
     //vertex position
     gener->SetSigma(1.8, 1.8,0);           //Sigma in (X,Y,Z) (cm) on IP position
     gener->SetPart(kPiPlus);                    //GEANT particle type
   break;
   case kScan:  
     AliGenScan *gener = new AliGenScan(-1);
     gener->SetMomentumRange(2,2);
     //gener->SetPhiRange(30,30);           //for inclined HMPID
     gener->SetPhiRange(90,90);             //for normal HMPID
     gener->SetThetaRange(90,90);
     //vertex position
     gener->SetSigma(0,0,0);           //Sigma in (X,Y,Z) (cm) on IP position
     gener->SetPart(kPiPlus); 
     //gener->SetRange(1, 415, 415, 1, 245, 245, 1, -20, -20);   //for inclined HMPID
     gener->SetRange(1, 0, 0, 1, 480, 480, 1, -20, -20);         //for normal HMPID

   break;
   case kSoubleScan:  
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
   case kHijing:
     AliGenHIJINGpara *gener = new AliGenHIJINGpara(ntracks);
     gener->SetMomentumRange(0,999);
     gener->SetPhiRange(0,360);
     gener->SetThetaRange(.77,179.23);
     gener->SetOrigin(0,0,0);          //vertex position
     gener->SetSigma(0,0,5.6);         //Sigma in (X,Y,Z) (cm) on IP position
   break;
   case kPythia:
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
   case kParam:
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
   case kFluka:
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
   case kNtuple:
     AliGenExtFile *gener = new AliGenExtFile(-1); 
     gener->SetFileName("$(ALICE_ROOT)/data/dtujet93.root");
     gener->SetVertexSmear(perEvent); 
     gener->SetTrackingFlag(1);
   break;
   case kHalo:
     AliGenHalo *gener = new AliGenHalo(ntracks); 
     gener->SetFileName("/h1/morsch/marsip/marsip5.mu");
   break;     
   case kCocktail:
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
   }//switch
 
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
  //=================== ITS parameters ===========================
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
  //AliITS *ITS  = new AliITSv5asymm("ITS","Updates ITS TDR detailed version with asymmetric services");
  //
  //
  // Coarse geometries (warning: no hits are produced with these coarse geometries and they unuseful for reconstruction !):
  // ======================================================================================================================
  //
  //
  //AliITS *ITS  = new AliITSv1("ITS","Old ITS coarse version as of the ALICE TP");
  //
  AliITS *ITS  = new AliITSvPPRcoarseasymm("ITS","New ITS coarse version with asymmetric services");
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
// ---           sectors 25-72 are the upper ones (25-48 -> z>0, 49-72 -> z<0)
// --- SecLows - number of lower sectors specified (up to 6)
// --- SecUps - number of upper sectors specified (up to 12)
// --- Sens - sensitive strips for the Slow Simulator !!!
// --- This does NOT work if all S or L-sectors are specified, i.e.
// --- if SecAL or SecAU < 0
//
//
//-----------------------------------------------------------------------------

 /* gROOT->LoadMacro("SetTPCParam.C");
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

  if (TPC->IsVersion()==1) param->Write(param->GetTitle());*/

   AliTPC *TPC  = new AliTPCv0("TPC","Default");
   // All sectors included 
   TPC->SetSecAL(-1);
   TPC->SetSecAU(-1);

}

if(iTOF) {
//=================== TOF parameters ============================
AliTOF *TOF  = new AliTOFv2("TOF","normal TOF");
}
   

 gAlice->SetDebug(1);


if(iRICH){
//====================== RICH parameters =========================
AliRICH *RICH  = new AliRICHv3("RICH","normal RICH");
}//if(iRICH)


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

   cout<<"RICH private Config.C> End\n";         
}//void Config()
