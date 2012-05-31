// inserted you find the pixel routine: From a macro it is called like:

//    Int_t spdTriggerBit = CalculateFOTrigger(itsLoader,geom,hDigits,nFO);


// spdTriggerBit is a bit array which can be queried like:

// if (spdTriggerBit & (1 << 1)) globalFO = kTRUE;

//The arguments hDigits and NFO are things that you don't need. The routine
//CalculateFOTrigger uses anouther subroutine RequireZ10cm() which
//calculates the VERTEX coincidence requirement.

//Don't hesitate to contact me in case of questions. My suggestion to you is
//that you first only implement a simple version of the trigger (for example
//the OR of all signals).


//Best regards,

//Jan



Int_t CalculateFOTrigger( AliITSLoader *itsl, AliITSgeom* geom, TH1F* hDigits, Float_t& nFO)
{

  TObjArray *digDet = 0;
  digDet = new TObjArray(3);

  Int_t startSPD = geom->GetStartSPD();
  Int_t lastSPD  = geom->GetLastSPD();

  // Cut on Signal In the Pixel Detector
  treeD = itsl->TreeD();
  br = treeD->GetBranch("ITSDigitsSPD");
  br->SetAddress(&((*digDet)[0]));
  ((TClonesArray*)(digDet->At(0)))->Clear();

  Int_t ndig = 0;
  Int_t ndigfo = 0;

  Int_t singhitthreshold = 1; // single hit threshold
  Int_t threshold = 1;

  const Int_t nlay =2;
  const Int_t nlad =240;
  const Int_t nstave=40;
  const Int_t ndet =4;
  const Int_t nchip =5;
  const Int_t ntotal = 1200;

  Int_t ndigA[5];
  Int_t FOperlayer[2];
  Int_t FOperladder[240];
  Int_t FOperstave[40][2];
  Int_t FOperchip[ntotal];
  Int_t FOperChipinStave[20][40][2];

  for (Int_t m=startSPD;m<lastSPD;m++) {
    FOperladder[m] = 0;
  }

  for (Int_t k = 0;k<2;k++){
    FOperlayer[k] = 0;
    for (Int_t o=0;o<40;o++) {
      FOperstave[o][k] = 0;
      for (Int_t ich=0;ich<20;ich++) {
        FOperChipinStave[ich][o][k] = 0;
      }
    }
  }


  nFO=0.0;

  Int_t mInStaveCounter = 0.;

  Int_t checkStave = 0;


  for (m=startSPD;m<lastSPD;m++) {

    treeD->GetEvent(m);
    digits = (TClonesArray*) (digDet->At(0)); // SPD only.


    Int_t lay,stav,det;  geom->GetModuleId(m,lay,stav,det);

    ndig = digits->GetEntriesFast();

    for(Int_t l = 0;l<5;l++){
      ndigA[l] =0 ;
    }

    //
   
    for(Int_t dig=0; dig<ndig; dig++){
      dp = (AliITSdigitSPD*) digits->At(dig);
      Int_t column = dp->GetCoord1();
      //   Int_t row = dp->GetCoord2();
      Int_t chip = Int_t(column/32.);
      ndigA[chip]++;
    }

    if (checkStave !=stav){
      mInStaveCounter = 0;
    } else{
      mInStaveCounter += 1;
    }

    // m 0,.., 239
    // stav 1,..,40
    // mInStave 0,..,3
    // chipInStave 0,..,19
  
    //cout << "m " << m << " stav "  << stav  << " mInStave " <<
mInStaveCounter << " " <<lay << endl;

    for (Int_t ichip=0;ichip<5;ichip++){
      Int_t seq = (m*5+ichip);
      Int_t chipInStave = (mInStaveCounter *5) + ichip;

      if (ndigA[ichip] >= 1) {
        FOperladder[m]++;
        FOperlayer[lay-1]++;
        FOperstave[stav-1][lay-1]++;
        //FOperHstave[hstav-1][lay-1]++;
        FOperChipinStave[chipInStave][stav-1][lay-1]++;
        nFO++;
     }
     
    }
    // SIMPLE FO ---> ANY HIT TRIGGERS
    ndigfo += ndig;   
    checkStave = stav;

  }

  //cout << 2 << endl;

   Int_t bit1 = 0;
   Int_t upper_cut = 120;
 
   hDigits->Fill(ndigfo);

   
   //   nFO  = ndigfo;


   //cout << nFO << endl;

  if ( ndigfo >= singhitthreshold) {bit1 |= (1 << 1);}
 
  //  if ( ndigfo <= upper_cut) {bit1 |= (1 << 10);}



  if (FOperlayer[0] >= threshold && FOperlayer[1] >=threshold) {
     bit1 |= (1 << 2);
     if ( ndigfo <= upper_cut) {bit1 |= (1 << 10);}
  }
 

  // Sector coincidence

  Int_t nsec = 0;
  Int_t finstav = 0;

  // staves layer 1  1-20
  // staves layer 2: 0-39

  //cout << 3 << endl;


   for (Int_t istav=1;istav<21;istav++){
     for (Int_t jstav = finstav; jstav<finstav+4;jstav++) {
       if ((FOperstave[istav-1][0] >= threshold) &&
           (FOperstave[jstav][1] >= threshold)) {
         bit1 |= (1 << 3);
        
          
     if (RequireZ10cm(FOperChipinStave,istav-1,jstav) == kTRUE) {
             //cout << (RequireZ(FOperChipinStave,istav-1,probe_stav) <<
endl;
              bit1 |= (1 << 7);
     }
 
       //   cout << " " << istav << " " << jstav << " " << endl;
       }   
  
     }
     if (TMath::Even(istav)) {
       finstav = jstav;
     }
    
   }
   //cout << 4 << endl;

   // half sector coincidence

   Int_t finstav = 0;
    for (Int_t istav=1;istav<21;istav++){
     for (Int_t jstav = finstav; jstav<finstav+2;jstav++)
     {
       if (FOperstave[istav-1][0] >= threshold && FOperstave[jstav][1] >=
threshold)
       {
         bit1 |= (1 << 4);

           if (RequireZ10cm(FOperChipinStave,istav-1,jstav) == kTRUE) {
           
             bit1 |= (1 << 8);
           }
       }
     }
   
      finstav = jstav;

    }
  

 
  
    Int_t finstav = 0;
  
    for (Int_t istav=1;istav<21;istav++){
    
      for (Int_t jstav = finstav-1; jstav<finstav+3;jstav++) {
      
       Int_t probe_stav = jstav;

       if (jstav == -1) probe_stav = 39;
       if (jstav == 40) probe_stav = 0;
      
       if (FOperstave[istav-1][0] >= threshold &&
FOperstave[probe_stav][1] >= threshold) {
         bit1 |= (1 << 5);
       }
     
      }

     finstav = jstav-1;
    
    }



 
  // sliding window coincidence (symmetric): 1 (layer 1), 5 (layer 2)

   Int_t finstav = 0;

   for (Int_t istav=1;istav<21;istav++) {
     for (Int_t jstav = finstav-2; jstav<finstav+3;jstav++) {
      
       // cout << 7 << endl;

       probe_stav = jstav;
       if (jstav == -2) probe_stav = 38;
       if (jstav == -1) probe_stav = 39;
       if (jstav == 40) probe_stav = 0;
       if (jstav == 41) probe_stav = 1;
      
       if ((FOperstave[istav-1][0] >= threshold) &&
           (FOperstave[probe_stav][1] >= threshold)) {
            bit1 |= (1 << 6);
            if (RequireZ10cm(FOperChipinStave,istav-1,probe_stav) ==
kTRUE) {
                     bit1 |= (1 << 9);
            }
     
       

       }
     }
     finstav = jstav-1;
   }

  return bit1;
}



Bool_t RequireZ10cm(Int_t FOperChipinStave[][40][2], Int_t stave1, Int_t
stave2){

// z  sliding window

  Bool_t zFlag = kFALSE;

  Int_t threshold = 1;
  Int_t start1 = 1;
  Int_t start2 = 3;
  Int_t start3 = 6;
  Int_t start4 = 7;
  Int_t i = 1;

   for (Int_t ic=0;ic<=19;ic++) {
    

       if(ic <= 5) { 


         for (Int_t jc=0;jc<=ic-1;jc++) {
          if (FOperChipinStave[ic][stave1][0] >=threshold) {
           if (FOperChipinStave[jc][stave2][1] >= threshold){
             zFlag = kTRUE;
           }
           }
         }
       }
      
       if(ic >=6 && ic <=8){
         for (jc=(2*start1-1);jc<=(2*start1-1)+5;jc++) {
           if (FOperChipinStave[ic][stave1][0] >=threshold) {
           if  (FOperChipinStave[jc][stave2][1] >= threshold){
             zFlag = kTRUE;
           }
           }

         }
         start1++;
       }

       if(ic >=9 && ic <=11){

         for (jc=(2*start2);jc<=(2*start2+5);jc++) {
           if (FOperChipinStave[ic][stave1][0] >=threshold) {
           if  (FOperChipinStave[jc][stave2][1] >= threshold){
             zFlag = kTRUE;
           }
           }
         }
         start2++;
       }
      
       if(ic >=12 && ic <=13){

         for (jc=(2*start3-1);jc<=(2*start3-1)+5;jc++) {
           if (FOperChipinStave[ic][stave1][0] >=threshold) {
           if  (FOperChipinStave[jc][stave2][1] >= threshold){
             zFlag = kTRUE;
           }
           }
         }
         start3++;
       }


        if(ic >=14){
         for (jc=(2*start4);jc<=19;jc++) {
           if (FOperChipinStave[ic][stave1][0] >=threshold) {
           if  (FOperChipinStave[jc][stave2][1] >= threshold){
             zFlag = kTRUE;
           }
           }
         }
         start4++;

        }

   }

   return zFlag;
  
}