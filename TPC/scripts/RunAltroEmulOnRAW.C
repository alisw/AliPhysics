/// \file RunAltroEmulOnRAW.C



// +++++++++++++++++++++++++++++++++++++++++++++++++
void RunAltroEmulOnRAW(Int_t mode, const char *uri,
                           const char *outputRootFileName ="rawdata.root"){

  // modes

  // 0: MAF off / TCF off - ONLY ZS

  // 1: MAF on / TCF off
  // 2: MAF on / TCF on - params from Rehak (baseline)
  // 3: MAF on / TCF on - params from Rehak (baseline & pulse shortening)
  // 4: MAF on / TCF on - params from Rehak (baseline, used in P2)
  // 5: MAF on - presample 1, postsample 1
  // 6: MAF on - presample 1, postsample 3
  // 7: MAF on - presample 1, postsample 5
 

  if (mode<0 || mode>7) {
    printf(" ERROR: chosen mode not avalailable ...\n");
    return;
  }


  TStopwatch timer;
  timer.Start();

  // READ RAW DATA +++++++++++

  AliRawReader *reader=AliRawReader::Create(uri);
  reader->NextEvent();
  Int_t runNumber=reader->GetRunNumber();
  reader->RewindEvents();

  // ALTRO EMULATOR +++++++++++

  Int_t baseline = 50; // standard in run 137228

  Int_t *K0= new Int_t[2];  Int_t *K1= new Int_t[2];  Int_t *K2= new Int_t[2];
  Int_t *L0= new Int_t[2];  Int_t *L1= new Int_t[2];  Int_t *L2= new Int_t[2];
 
  if ( mode == 2 ) { // "middle aggressive"   
    // params from Rehak (baseline)
    K0[0]=64386;  K1[0]=65184;  K2[0]=   77;
    L0[0]=64675;  L1[0]=64950;  L2[0]=    0;
    K0[1]=63675;  K1[1]=65355;  K2[1]=  218;
    L0[1]=63917;  L1[1]=65263;  L2[1]=    0;

  } else  if ( mode ==3 ) { // "most aggressive"
    // params from Rehak (baseline & pulse shortening)
    K0[0]=58684;  K1[0]=65273;  K2[0]= 1361;   // IROC
    L0[0]=59754;  L1[0]=65187;  L2[0]=    0;
    K0[1]=59688;  K1[1]=65347;  K2[1]= 1349;   // OROC
    L0[1]=60704;  L1[1]=65273;  L2[1]=    0;

  } else  if ( mode == 4 ) {  // a bit "aggressive"
    // params from Rehak (baseline, used in P2);
    K0[0]=62055;  K1[0]=65300;  K2[0]=  277;  // IROC
    L0[0]=62339;  L1[0]=65208;  L2[0]=    0; 
    K0[1]=63151;  K1[1]=65435;  K2[1]=  266;  // OROC  
    L0[1]=63387;  L1[1]=65371;  L2[1]=    0;     
  }


  AliTPCAltroEmulator *a1 = new AliTPCAltroEmulator();

  a1->ConfigBaselineCorrection1(4, baseline, 0, 0);
  a1->ConfigZerosuppression(3,2,0,0);
   
  if (mode ==0) { // do nothing, just ZS
 
    a1->ConfigAltro(1,0,0,1,1,0); 

  } else if (mode==1) {  // MAF

    a1->ConfigBaselineCorrection2(3,3,0,0,0);
    a1->ConfigAltro(1,0,1,1,1,0);
    
  } else if (mode==2 || mode ==3 || mode==4) { // MAF and different TCF  

    a1->ConfigTailCancellationFilterForRAWfiles(K0,K1,K2,L0,L1,L2);
    a1->ConfigBaselineCorrection2(3,3,0,0,0);
    a1->ConfigAltro(1,1,1,1,1,0); 

  } else if (mode == 5) { // different MAF
    
    a1->ConfigBaselineCorrection2(3,3,0,1,1);  // 1 presample, 1 postsamples
    a1->ConfigAltro(1,0,1,1,1,0);

  } else if (mode == 6) { // different MAF
    
    a1->ConfigBaselineCorrection2(3,3,0,1,3);  // 1 presample, 3 postsamples
    a1->ConfigAltro(1,0,1,1,1,0);
    
  } else if (mode == 7) { // different MAF
    
    a1->ConfigBaselineCorrection2(3,3,0,1,5);  // 1 presample, 5 postsamples
    a1->ConfigAltro(1,0,1,1,1,0);

  } 


  a1->SetOutputRootFileName(outputRootFileName);
  a1->RunEmulationOnRAWdata(reader);
  
  timer.Stop(); timer.Print();
}


