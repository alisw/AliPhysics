void Config(Int_t version)
{
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libminicern");
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libMUON");

  AliMUON *MUON = 0;
  switch (version) {
    case 0: MUON  = new AliMUONv0("MUON","normal MUON"); break;
  }  

//=================== MUON parameters ===========================
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
