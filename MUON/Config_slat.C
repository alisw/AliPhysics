enum gentype_t {hijing, gun, box, pythia, param, cocktail, fluka, halo, ntuple, scan, doublescan, hijing_g};

gentype_t gentype=scan;
// Int_t ntracks=6407;
// Int_t ntracks=12000;
// Int_t ntracks=28380;
// Int_t ntracks=19900;
Int_t ntracks=1;

void Config()

{
new AliGeant3("C++ Interface to Geant3");

//=======================================================================
//  Create the output file
   
 TFile *rootfile = new TFile("galice.root","recreate");
 rootfile->SetCompressionLevel(2);
 TGeant3 *geant3 = (TGeant3*)gMC;
 AliDecayer* decayer = new AliDecayerPythia();
 decayer->SetForceDecay(all);
 decayer->Init();
 gMC->SetExternalDecayer(decayer);

//=======================================================================
// ******* GEANT STEERING parameters FOR ALICE SIMULATION *******
 geant3->fGctrak->maxnst=1000000;
 geant3->fGcflag->nrndm[0]=10;
 geant3->fGcflag->nrndm[1]=11; 
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

 gAlice->TrackingLimits( 700, 2000);
 
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
     gener->SetRange(60, -300, 300, 60, -300., 300., 1, 900, 900);
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
     gener->SetStrucFunc(DO_Set_1);
     gener->SetProcess(charm);
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
// Example for J/psi  Production from  Parameterisation *
//*******************************************************
     AliGenParam *gener = new AliGenParam(ntracks,upsilon_p);
     gener->SetMomentumRange(0,999);
     gener->SetPhiRange(-180, 180);
     gener->SetYRange(2.5,4);
     gener->SetCutOnChild(1);
     gener->SetChildThetaRange(2,9);
     gener->SetOrigin(0,0,0);          //vertex position
     gener->SetSigma(0,0,0);           //Sigma in (X,Y,Z) (cm) on IP position
     gener->SetForceDecay(dimuon);
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
     AliGenParam *Pi0 = new AliGenParam(100, new AliGenPMDlib(), Pion);
     AliGenParam *Eta = new AliGenParam( 10, new AliGenPMDlib(), Eta);
     gener->AddGenerator(Pi0, "neutral pions"  , 1.);
     gener->AddGenerator(Eta, "neutral etas"  ,  1.);
     break;
 }
 

gener->Init();
gAlice->SetField(2,1);    //Specify maximum magnetic field in Tesla (neg. ==> default field)

Int_t iFRAME  =0;
Int_t iMAG    =0;
Int_t iITS    =0;
Int_t iABSO   =0;
Int_t iDIPO   =0;
Int_t iHALL   =0;
Int_t iSHIL   =0;
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
    AliSHILv0 *SHIL  = new AliSHILv0("SHIL","Shielding");
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

AliMUON *MUON  = new AliMUONv1("MUON","normal MUON");
 MUON->SetIshunt(1);
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
 // Default response: 5 mm of gas
 AliMUONResponseV0* response0 = new AliMUONResponseV0;
 response0->SetSqrtKx3AndDeriveKx2Kx4(0.7131); // sqrt(0.5085)
 response0->SetSqrtKy3AndDeriveKy2Ky4(0.7642); // sqrt(0.5840)
 response0->SetPitch(0.25); // anode-cathode distance
 response0->SetSigmaIntegration(10.);
 response0->SetChargeSlope(50);
 response0->SetChargeSpread(0.18, 0.18);
 response0->SetMaxAdc(4096);
 response0->SetZeroSuppression(6);

 // Response for 4 mm of gas (station 1)
 // automatic consistency with width of sensitive medium in CreateGeometry ????
 AliMUONResponseV0* responseSt1 = new AliMUONResponseV0;
 // Mathieson parameters from L.Kharmandarian's thesis, page 190
 responseSt1->SetSqrtKx3AndDeriveKx2Kx4(0.7000); // sqrt(0.4900)
 responseSt1->SetSqrtKy3AndDeriveKy2Ky4(0.7550); // sqrt(0.5700)
 responseSt1->SetPitch(0.20); // anode-cathode distance
 responseSt1->SetSigmaIntegration(10.);
 // ChargeSlope larger to compensate for the smaller anode-cathode distance
 // and keep the same most probable ADC channel for mip's
 responseSt1->SetChargeSlope(62.5); 
 // assumed proportionality to anode-cathode distance for ChargeSpread
 responseSt1->SetChargeSpread(0.144, 0.144);
 responseSt1->SetMaxAdc(4096);
 responseSt1->SetZeroSuppression(6);

//--------------------------------------------------------
// Configuration for Chamber TC1/2  (Station 1) ----------           
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 Float_t rseg1[4]={17.5, 55.2, 71.3, 95.5};
 Int_t   nseg1[4]={4, 4, 2, 1};
//
 chamber=1;
//^^^^^^^^^
 MUON->SetNsec(chamber-1,2);
//
 AliMUONSegmentationV01 *seg11=new AliMUONSegmentationV01;
 
 seg11->SetSegRadii(rseg1);
//  seg11->SetPadSize(3, 0.5);
 seg11->SetPadSize(2.4, 0.4); // smaller pad size
//  seg11->SetDAnod(3.0/3./4);
 seg11->SetDAnod(0.20); // smaller distance between anode wires
 seg11->SetPadDivision(nseg1);
 
 MUON->SetSegmentationModel(chamber-1, 1, seg11);
//
 AliMUONSegmentationV02 *seg12=new AliMUONSegmentationV02;
 seg12->SetSegRadii(rseg1); 
//  seg12->SetPadSize(0.75, 2.0);
 seg12->SetPadSize(0.6, 1.6); // smaller pad size
//  seg12->SetDAnod(3.0/3./4);
 seg12->SetDAnod(0.20); // smaller distance between anode wires
 seg12->SetPadDivision(nseg1);

 MUON->SetSegmentationModel(chamber-1, 2, seg12);

//  MUON->SetResponseModel(chamber-1, response0);	    
 MUON->SetResponseModel(chamber-1, responseSt1); // special response	    

 chamber=2;
//^^^^^^^^^
//
 MUON->SetNsec(chamber-1,2);
//
 AliMUONSegmentationV01 *seg21=new AliMUONSegmentationV01;
 seg21->SetSegRadii(rseg1);
//  seg21->SetPadSize(3, 0.5);
 seg21->SetPadSize(2.4, 0.4); // smaller pad size
//  seg21->SetDAnod(3.0/3./4);
 seg21->SetDAnod(0.20); // smaller distance between anode wires
 seg21->SetPadDivision(nseg1);
 MUON->SetSegmentationModel(chamber-1, 1, seg21);
//
 AliMUONSegmentationV02 *seg22=new AliMUONSegmentationV02;
 seg22->SetSegRadii(rseg1); 
//  seg22->SetPadSize(0.75, 2.);
 seg22->SetPadSize(0.6, 1.6); // smaller pad size
//  seg22->SetDAnod(3.0/3./4);
 seg22->SetDAnod(0.20); // smaller distance between anode wires
 seg22->SetPadDivision(nseg1);
 MUON->SetSegmentationModel(chamber-1, 2, seg22);

//  MUON->SetResponseModel(chamber-1, response0);	    
 MUON->SetResponseModel(chamber-1, responseSt1); // special response
	    
//
//--------------------------------------------------------
// Configuration for Chamber TC3/4 (Station 2) -----------
///^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// Float_t rseg2[4]={23.5, 87.7, 122.4, 122.5};
     Float_t rseg2[4]={23.5, 47.1, 87.7, 122.5};
     Int_t   nseg2[4]={4, 4, 2, 1};
//
 chamber=3;
//^^^^^^^^^
 MUON->SetNsec(chamber-1,2);
//
 AliMUONSegmentationV01 *seg31=new AliMUONSegmentationV01;
 seg31->SetSegRadii(rseg2);
 seg31->SetPadSize(3.0, 0.5);
 seg31->SetDAnod(3.0/3./4);
 seg31->SetPadDivision(nseg2);
 MUON->SetSegmentationModel(chamber-1, 1, seg31);
//
 AliMUONSegmentationV02 *seg32=new AliMUONSegmentationV02;
 seg32->SetSegRadii(rseg2); 
 seg32->SetPadSize(0.75, 2.0);
 seg32->SetPadDivision(nseg2);
 seg32->SetDAnod(3.0/3./4);

 MUON->SetSegmentationModel(chamber-1, 2, seg32);

 MUON->SetResponseModel(chamber-1, response0);	    

 chamber=4;
//^^^^^^^^^
//
 MUON->SetNsec(chamber-1,2);
//
 AliMUONSegmentationV01 *seg41=new AliMUONSegmentationV01;
 seg41->SetSegRadii(rseg2);
 seg41->SetPadSize(3.0, 0.5);
 seg41->SetDAnod(3.0/3./4);
 seg41->SetPadDivision(nseg2);
 MUON->SetSegmentationModel(chamber-1, 1, seg41);
//
 AliMUONSegmentationV02 *seg42=new AliMUONSegmentationV02;
 seg42->SetSegRadii(rseg2); 
 seg42->SetPadSize(0.75, 2.0);
 seg42->SetPadDivision(nseg2);
 seg42->SetDAnod(3.0/3./4);

 MUON->SetSegmentationModel(chamber-1, 2, seg42);

 MUON->SetResponseModel(chamber-1, response0);	    


//--------------------------------------------------------
// Configuration for Chamber TC5/6  (Station 3) ----------          
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 Int_t   nseg3[4]={4, 4, 2, 1};
 Int_t   npcb5[36] = {0,0,2,0,
		      0,0,3,0,
		      0,1,3,0,
		      0,2,2,0,
		      0,1,2,0, 
		      0,2,2,0, 
		      0,1,3,0, 
		      0,0,3,0,
                      0,0,2,0};

 Float_t shift = 1.5/2.;
 // Float_t xpos5[8]    = {2., 2., 2., 42., 42., 2., 2., 2.};
 Float_t xpos5[9]    = {2., 2., 2., 2.,32., 2., 2., 2., 2.};
 Float_t ypos5       = -(20.+4.*(40.-2.*shift));

 chamber=5;
 MUON->SetNsec(chamber-1,2);
 AliMUONSegmentationSlat *seg51=new AliMUONSegmentationSlat;
 seg51->SetNSlats(9); 
 seg51->SetShift(shift);  
 seg51->SetNPCBperSector(npcb5); 
 seg51->SetSlatXPositions(xpos5);
 seg51->SetSlatYPosition(ypos5);
 seg51->SetPadSize(10.,0.5);
 seg51->SetDAnod(0.25);
 seg51->SetPadDivision(nseg3);
 MUON->SetSegmentationModel(chamber-1, 1, seg51);

 AliMUONSegmentationSlatN *seg52=new AliMUONSegmentationSlatN;
 seg52->SetNSlats(9); 
 seg52->SetShift(shift);  
 seg52->SetNPCBperSector(npcb5); 
 seg52->SetSlatXPositions(xpos5);
 seg52->SetSlatYPosition(ypos5);
 seg52->SetPadSize(1., 10.); // DeltaX(non bending) = 2 * DeltaY(bending)
 seg52->SetDAnod(0.25);
 seg52->SetPadDivision(nseg3);
 MUON->SetSegmentationModel(chamber-1, 2, seg52);
 MUON->SetResponseModel(chamber-1, response0);	    

 chamber=6;
 MUON->SetNsec(chamber-1,2);
 AliMUONSegmentationSlat *seg61=new AliMUONSegmentationSlat;
 seg61->SetNSlats(9); 
 seg61->SetShift(shift);  
 seg61->SetNPCBperSector(npcb5); 
 seg61->SetSlatXPositions(xpos5);
 seg61->SetSlatYPosition(ypos5);
 seg61->SetPadSize(10.,0.5);
 seg61->SetDAnod(0.25);
 seg61->SetPadDivision(nseg3);
 MUON->SetSegmentationModel(chamber-1, 1, seg61);

 AliMUONSegmentationSlatN *seg62=new AliMUONSegmentationSlatN;
 seg62->SetNSlats(9); 
 seg62->SetShift(shift);  
 seg62->SetNPCBperSector(npcb5); 
 seg62->SetSlatXPositions(xpos5);
 seg62->SetSlatYPosition(ypos5);
 seg62->SetPadSize(1., 10.); // DeltaX(non bending) = 2 * DeltaY(bending)
 seg62->SetDAnod(0.25);
 seg62->SetPadDivision(nseg3);
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
 AliMUONSegmentationSlat *seg71=new AliMUONSegmentationSlat;
 Int_t npcb7[44] = {0,0,0,3,
		    0,0,2,2,
		    0,0,3,2,
		    0,2,2,1,
		    0,2,2,1,
                    0,1,2,1, 
		    0,2,2,1, 
		    0,2,2,1, 
		    0,0,3,2, 
		    0,0,2,2, 
		    0,0,0,3};
 Float_t xpos7[11]   = {2., 2., 2., 2., 2., 39.5, 2., 2., 2., 2., 2.};
 Float_t ypos7       = -(20.+5.*(40.-2.*shift));
 
 seg71->SetNSlats(11);  
 seg71->SetShift(shift);  
 seg71->SetNPCBperSector(npcb7); 
 seg71->SetSlatXPositions(xpos7);
 seg71->SetSlatYPosition(ypos7);
 
 seg71->SetPadSize(10.,0.5);
 seg71->SetDAnod(0.25);
 seg71->SetPadDivision(nseg4);
 MUON->SetSegmentationModel(chamber-1, 1, seg71);

 AliMUONSegmentationSlatN *seg72=new AliMUONSegmentationSlatN;

 MUON->SetSegmentationModel(chamber-1, 2, seg72);
 seg72->SetNSlats(11);  
 seg72->SetShift(shift);   
 seg72->SetNPCBperSector(npcb7); 
 seg72->SetSlatXPositions(xpos7);
 seg72->SetSlatYPosition(ypos7);
 seg72->SetPadSize(1., 10.); // DeltaX(non bending) = 2 * DeltaY(bending)
 seg72->SetDAnod(0.25);
 seg72->SetPadDivision(nseg4);

 MUON->SetResponseModel(chamber-1, response0);	    

 chamber=8;
//^^^^^^^^^
 MUON->SetNsec(chamber-1,2);
//
 AliMUONSegmentationSlat *seg81=new AliMUONSegmentationSlat;

 seg81->SetNSlats(11);  
 seg81->SetShift(shift);  
 seg81->SetNPCBperSector(npcb7); 
 seg81->SetSlatXPositions(xpos7);
 seg81->SetSlatYPosition(ypos7);
 seg81->SetPadSize(10.,0.5);
 seg81->SetDAnod(0.25);
 seg81->SetPadDivision(nseg4);
 MUON->SetSegmentationModel(chamber-1, 1, seg81);

 AliMUONSegmentationSlat *seg82=new AliMUONSegmentationSlatN;

 MUON->SetSegmentationModel(chamber-1, 2, seg82);
 seg82->SetNSlats(11);  
 seg82->SetShift(shift);  
 seg82->SetNPCBperSector(npcb7); 
 seg82->SetSlatXPositions(xpos7);
 seg82->SetSlatYPosition(ypos7);
 seg82->SetPadSize(1., 10.); // DeltaX(non bending) = 2 * DeltaY(bending)
 seg82->SetDAnod(0.25);
 seg82->SetPadDivision(nseg4);

 MUON->SetResponseModel(chamber-1, response0);	    


//--------------------------------------------------------
// Configuration for Chamber TC9/10  (Station 5) ---------           
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 chamber=9;
//^^^^^^^^^

 MUON->SetNsec(chamber-1,2);
//
 AliMUONSegmentationSlat *seg91=new AliMUONSegmentationSlat;
 Int_t   npcb9[52] = {0,0,0,3,
		      0,0,0,4,
		      0,0,2,3,
		      0,0,3,3,
		      0,2,2,2,
		      0,2,2,2,
                      0,1,2,2, 
		      0,2,2,2, 
		      0,2,2,2, 
		      0,0,3,3, 
		      0,0,2,3, 
		      0,0,0,4, 
		      0,0,0,3};   

 // Float_t xpos9[13]   = {2., 2., 2., 2., 2., 2., 39.5 , 2., 2., 2., 2., 2., 2.};
 Float_t xpos9[13]   = {2., 2., 2., 2., 2., 2., 39.5, 2., 2., 2., 2., 2., 2.};
 Float_t ypos9       = -(20.+6.*(40.-2.*shift));

 seg91->SetNSlats(13);  
 seg91->SetShift(shift);  
 seg91->SetNPCBperSector(npcb9); 
 seg91->SetSlatXPositions(xpos9);
 seg91->SetSlatYPosition(ypos9);
 seg91->SetPadSize(10.,0.5);
 seg91->SetDAnod(0.25);
 seg91->SetPadDivision(nseg4);
 MUON->SetSegmentationModel(chamber-1, 1, seg91);

 AliMUONSegmentationSlatN *seg92=new AliMUONSegmentationSlatN;

 MUON->SetSegmentationModel(chamber-1, 2, seg92);
 seg92->SetNSlats(13);  
 seg92->SetShift(shift);   
 seg92->SetNPCBperSector(npcb9); 
 seg92->SetSlatXPositions(xpos9);
 seg92->SetSlatYPosition(ypos9);
 seg92->SetPadSize(1., 10.); // DeltaX(non bending) = 2 * DeltaY(bending)
 seg92->SetDAnod(0.25);
 seg92->SetPadDivision(nseg4);

 MUON->SetResponseModel(chamber-1, response0);	    

 chamber=10;
//^^^^^^^^^
 MUON->SetNsec(chamber-1,2);
//
 AliMUONSegmentationSlat *seg101=new AliMUONSegmentationSlat;
 
 seg101->SetNSlats(13);  
 seg101->SetShift(shift);  
 seg101->SetNPCBperSector(npcb9); 
 seg101->SetSlatXPositions(xpos9);
 seg101->SetSlatYPosition(ypos9);
 seg101->SetPadSize(10.,0.5);
 seg101->SetDAnod(0.25);
 seg101->SetPadDivision(nseg4);
 MUON->SetSegmentationModel(chamber-1, 1, seg101);

 AliMUONSegmentationSlatN *seg102=new AliMUONSegmentationSlatN;

 MUON->SetSegmentationModel(chamber-1, 2, seg102);
 seg102->SetNSlats(13);  
 seg102->SetShift(shift);   
 seg102->SetNPCBperSector(npcb9); 
 seg102->SetSlatXPositions(xpos9);
 seg102->SetSlatYPosition(ypos9);
 seg102->SetPadSize(1., 10.); // DeltaX(non bending) = 2 * DeltaY(bending)
 seg102->SetDAnod(0.25);
 seg102->SetPadDivision(nseg4);

 MUON->SetResponseModel(chamber-1, response0);	    

//--------------------------------------------------------
// Configuration for Trigger Stations -------------------- 
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// Cluster-size off
 AliMUONResponseTrigger* responseTrigger0 =  new AliMUONResponseTrigger;
// Cluster-size on  
// AliMUONResponseTriggerV1* responseTrigger0 =  new AliMUONResponseTriggerV1;
 
 chamber=11;
 MUON->SetNsec(chamber-1,2);
 AliMUONSegmentationTriggerX *seg111=new AliMUONSegmentationTriggerX;
 MUON->SetSegmentationModel(chamber-1, 1, seg111);
 AliMUONSegmentationTriggerY *seg112=new AliMUONSegmentationTriggerY;
 MUON->SetSegmentationModel(chamber-1, 2, seg112);

 MUON->SetResponseModel(chamber-1, responseTrigger0);      

 
 chamber=12;
 MUON->SetNsec(chamber-1,2);
 AliMUONSegmentationTriggerX *seg121=new AliMUONSegmentationTriggerX;
 MUON->SetSegmentationModel(chamber-1, 1, seg121);
 AliMUONSegmentationTriggerY *seg122=new AliMUONSegmentationTriggerY;
 MUON->SetSegmentationModel(chamber-1, 2, seg122);

 MUON->SetResponseModel(chamber-1, responseTrigger0);      
  chamber=13;
 MUON->SetNsec(chamber-1,2);
 AliMUONSegmentationTriggerX *seg131=new AliMUONSegmentationTriggerX;
 MUON->SetSegmentationModel(chamber-1, 1, seg131);
 AliMUONSegmentationTriggerY *seg132=new AliMUONSegmentationTriggerY;
 MUON->SetSegmentationModel(chamber-1, 2, seg132);
 MUON->SetResponseModel(chamber-1, responseTrigger0);      

 chamber=14;
 MUON->SetNsec(chamber-1,2);
 AliMUONSegmentationTriggerX *seg141=new AliMUONSegmentationTriggerX;
 MUON->SetSegmentationModel(chamber-1, 1, seg141);
 AliMUONSegmentationTriggerY *seg142=new AliMUONSegmentationTriggerY;
 MUON->SetSegmentationModel(chamber-1, 2, seg142);

 MUON->SetResponseModel(chamber-1, responseTrigger0); 

}
 
}

		



