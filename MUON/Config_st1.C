enum gentype_t {hijing, gun, box, pythia, param, cocktail, fluka, halo, ntuple, scan, doublescan, hijing_g};

gentype_t gentype=param;
// Int_t ntracks=6407;
// Int_t ntracks=12000;
// Int_t ntracks=28380;
// Int_t ntracks=19900;
Int_t ntracks=1;

void Config()

{
// Load geant321 library
 gSystem->Load("libgeant321");

 new TGeant3("C++ Interface to Geant3");

//=======================================================================
//  Create the output file
   
 TFile *rootfile = new TFile("galice.root","recreate");
 rootfile->SetCompressionLevel(2);
 TGeant3 *geant3 = (TGeant3*)gMC;
 AliDecayer* decayer = new AliDecayerPythia();
 decayer->SetForceDecay(kAll);
 decayer->Init();
 gMC->SetExternalDecayer(decayer);

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
     gener->SetMomentum(20);
     gener->SetPhiRange(0);
     gener->SetThetaRange(0.);
     gener->SetOrigin(30,30,500); //vertex position       
     gener->SetPart(kMuonMinus);  //GEANT particle type       
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
     gener->SetPhiRange(0,0);
     gener->SetThetaRange(0,0);
     //vertex position
//     gener->SetSigma(1,1,0);           //Sigma in (X,Y,Z) (cm) on IP position
     gener->SetPart(kMuonPlus); 
     gener->SetRange(30, -100., 100., 30, -100., 100., 1, 500, 500);
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
     //gener->SetSigma(0,0,5.6);         // Sigma in (X,Y,Z) (cm) on IP
position
     gener->SetStrucFunc(kDO_Set_1);
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
     gener->SetPtRange(0,999);
     gener->SetPhiRange(-180, 180);
     gener->SetYRange(2.5,9);
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
if (gentype==param) {
  gAlice->SetField(2,1); //Specify maximum magnetic field in Tesla (neg. ==> default field)
} 
else {
  gAlice->SetField(0,2); // No magnetic field
}


Int_t iFRAME  =0;
Int_t iMAG    =0;
Int_t iITS    =0;
Int_t iABSO   =1;
Int_t iDIPO   =1;
Int_t iHALL   =0;
Int_t iSHIL   =1;
Int_t iPIPE   =0;
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
//
// Parameters for selection MUON station 1 configuration:
//   detectorVersion:  
//     0:  AliMUONv1.cxx
//     2:  AliMUONv2.cxx  (detailed geometry for station 1)
//
//   responseVersion:  
//     0:  AliMUONResponseV0
//     2:  AliMUONSt1Response
//
//   segmentationVersion: 
//     0:  AliMUONSegmentationV01
//     2:  AliMUONSt1Segmentation

  Int_t detectorVersion = 2;        
  Int_t responseVersion = 2;    
  Int_t segmentationVersion = 2;

  // MUON detector
  //
  AliMUON *MUON = 0;
  if (detectorVersion == 2) {
    MUON= new AliMUONv2("MUON","normal MUON");
  } 
  else {
    MUON= new AliMUONv1("MUON","normal MUON");
  }
  MUON->SetIshunt(0);
  MUON->SetMaxStepGas(0.1);
  MUON->SetMaxStepAlu(0.1);

  // Response
  //
  AliMUONResponse* responseChamber1;
  AliMUONResponse* responseChamber2;

  if (responseVersion == 2) {
    AliMUONSt1Response* responseCh1 = new AliMUONSt1Response(1);
    responseCh1->SetSqrtKx3AndDeriveKx2Kx4(0.7000); // sqrt(0.4900)
    responseCh1->SetSqrtKy3AndDeriveKy2Ky4(0.7550); // sqrt(0.5700)
    responseCh1->SetPitch(0.20); // anode-cathode distance
    responseCh1->SetSigmaIntegration(10.);
       // Mathieson parameters from L.Kharmandarian's thesis, page 190
    responseCh1->SetChargeSlope(62.5);    //(62.5);
       // ChargeSlope larger to compensate for the smaller anode-cathode distance
       // and keep the same most probable ADC channel for mip's
    responseCh1->SetChargeSpread(0.144, 0.144);
       // assumed proportionality to anode-cathode distance for ChargeSpread
    responseCh1->SetMaxAdc(4095); 
    responseCh1->SetZeroSuppression(3);
    responseChamber1 = responseCh1;

    AliMUONSt1Response* responseCh2 = new AliMUONSt1Response(2);
    responseCh2->SetSqrtKx3AndDeriveKx2Kx4(0.7000); // sqrt(0.4900)
    responseCh2->SetSqrtKy3AndDeriveKy2Ky4(0.7550); // sqrt(0.5700)
    responseCh2->SetPitch(0.20); // anode-cathode distance
    responseCh2->SetSigmaIntegration(10.);
       // Mathieson parameters from L.Kharmandarian's thesis, page 190
    responseCh2->SetChargeSlope(62.5);    //(62.5);
       // ChargeSlope larger to compensate for the smaller anode-cathode distance
       // and keep the same most probable ADC channel for mip's
    responseCh2->SetChargeSpread(0.144, 0.144);
    responseCh2->SetMaxAdc(4095);
       // assumed proportionality to anode-cathode distance for ChargeSpread 
    responseCh2->SetZeroSuppression(3);
    responseChamber2 = responseCh2;
  }
  else {
    // Default response: 4 mm of gas
    AliMUONResponseV0* response0 = new AliMUONResponseV0;
    response0->SetSqrtKx3AndDeriveKx2Kx4(0.7131); // sqrt(0.5085)
    response0->SetSqrtKy3AndDeriveKy2Ky4(0.7642); // sqrt(0.5840)
    response0->SetPitch(0.2); // anode-cathode distance
    response0->SetSigmaIntegration(10.);
    response0->SetChargeSlope(62.5);
    response0->SetChargeSpread(0.144, 0.144);
    response0->SetMaxAdc(4095);
    response0->SetZeroSuppression(0);
  
    responseChamber1 = response0;
    responseChamber2 = response0;
  }   

  // Segmentation
  //
  AliSegmentation* segmentation11;
  AliSegmentation* segmentation12;
  AliSegmentation* segmentation21;
  AliSegmentation* segmentation22;
 
  if (segmentationVersion == 2) { 
    AliMUONSt1Segmentation* seg11 = new AliMUONSt1Segmentation(kBendingPlane);
    seg11->SetDAnod(0.20); // smaller distance between anod
    segmentation11 = seg11;

    AliMUONSt1Segmentation* seg12 = new AliMUONSt1Segmentation(kNonBendingPlane);
    seg12->SetDAnod(0.20); // smaller distance between anod
    segmentation12 = seg12;

    AliMUONSt1Segmentation* seg21 = new AliMUONSt1Segmentation(kBendingPlane);
    seg21->SetDAnod(0.20); // smaller distance between anod
    segmentation21 = seg21;

    AliMUONSt1Segmentation* seg22 = new AliMUONSt1Segmentation(kNonBendingPlane);
    seg22->SetDAnod(0.20); // smaller distance between anod
    segmentation22 = seg22;
  }
  else {	
    Float_t rseg1[4]={17.5, 55.2, 71.3, 95.5};
    Int_t   nseg1[4]={4, 4, 2, 1};

    AliMUONSegmentationV01* seg11=new AliMUONSegmentationV01(4);
    seg11->SetSegRadii(rseg1);
    seg11->SetPadSize(0.42, 0.63); // smaller pad size .....modif Marion 3/12/2
    seg11->SetDAnod(0.20); // smaller distance between anode wires
    seg11->SetPadDivision(nseg1);
    segmentation11 = seg11;

    AliMUONSegmentationV02* seg12=new AliMUONSegmentationV02(4);
    seg12->SetSegRadii(rseg1); 
    seg12->SetPadSize(0.42, 0.63); // smaller pad size.....modif Marion 3/12/2
    seg12->SetDAnod(0.20); // smaller distance between anode wires
    seg12->SetPadDivision(nseg1);
    segmentation12 = seg12;

    AliMUONSegmentationV01* seg21=new AliMUONSegmentationV01(4);
    seg21->SetSegRadii(rseg1);
    seg21->SetPadSize(2.4, 0.4); // smaller pad size
    seg21->SetDAnod(0.20); // smaller distance between anode wires
    seg21->SetPadDivision(nseg1);
    segmentation21 = seg21;

    AliMUONSegmentationV02* seg22=new AliMUONSegmentationV02(4);
    seg22->SetSegRadii(rseg1); 
    seg22->SetPadSize(0.6, 1.6); // smaller pad size
    seg22->SetDAnod(0.20); // smaller distance between anode wires
    seg22->SetPadDivision(nseg1);
    segmentation22 = seg22;
  } 
  
  // Configure station 1
  //
  Int_t chamber=1;
  MUON->SetNsec(chamber-1,2);
  MUON->SetSegmentationModel(chamber-1, 1, segmentation11);
  MUON->SetSegmentationModel(chamber-1, 2, segmentation12);
  MUON->SetResponseModel(chamber-1, responseChamber1);
  MUON->Chamber(chamber-1).SetChargeCorrel(0.11); // 11% charge spread
    
  chamber=2;
  MUON->SetNsec(chamber-1,2);
  MUON->SetSegmentationModel(chamber-1, 1, segmentation21);
  MUON->SetSegmentationModel(chamber-1, 2, segmentation22);
  MUON->SetResponseModel(chamber-1, responseChamber2);
  MUON->Chamber(chamber-1).SetChargeCorrel(0.11); // 11% charge spread

  //--------------------------------------------------------
  // Other stations 2, 3, 4, 5, 6 (Trigger) - default setting
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  AliMUONFactory factory;
  for (Int_t i=2; i<7; i++) factory.BuildStation(MUON, i);

     // Was: response0->SetMaxAdc(4095);
     // In AliMUONFactory is: fResponse0->SetMaxAdc(4096);

     // Station2:
     // Was: Float_t rseg2[4]={23.5, 47.1, 87.7, 122.5};
     // In In AliMUONFactory: Float_t rseg2[4]={23.5, 53.5, 90.5, 122.5};       

     // in macros: SHILv2
     // in Config_slat: SHILvF
}
}
