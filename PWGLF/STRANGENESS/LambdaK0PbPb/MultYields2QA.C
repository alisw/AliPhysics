#include "AliMassFitControl.h"

// Two functions now.
// MultYields3 takes a 3D histogram and calls the 2D function MultYields2
// MultYields2 could be called directly if we only have a 2D histogram

void MultYields3(TH3F *PtMassMult, Int_t particleMode, Int_t MultBin, Char_t* label){

  // Leave open possibility to choose different values depending on MultBin - LSB
  Float_t MultLo[1] = {0};
  Float_t MultHi[1] = {300};

  //Make 2D projections from the 3D histogram
  PtMassMult->GetZaxis()->SetRange(MultLo[MultBin],MultHi[MultBin]);
  TH2F* hParMass = (TH2F*)PtMassMult->Project3D("XY");// FIX:MF
  hParMass->SetTitle("PtMass");
  
  //hParMass->Draw(); // Drawing invokes default c1 canvas making it inaccessible in MultYields2 - TO FIX
  MultYields2(hParMass, particleMode,MultBin,label);
  
}

void MultYields2QA(TH2F *hParMass, Int_t particleMode, Int_t ihist,Int_t Nev = 1, Int_t MultBin, Char_t* label){

  hParMass->Draw();
  ////////////////////////// minhist is needed by AliMassFitControl.h (in case the minimum on the parameter's axis is arbitrary value); 
  Double_t minhist = 0;
  minhist = hParMass->GetBinLowEdge(1);
  cout<<" value of the first bin of parameter : " <<minhist << endl;

  // Do .L on command line first
  //gROOT->LoadMacro("macros/PtMassAna2.C");
  /* Modifications to produce a single uncorrected spectrum in a particular mult. bin for K0 or Lambda
     Old ratio code preserved in MultYieldsRatio.C
     Dec 2003 */
  Char_t* part; //for name of particle
  Float_t BR; //branching ratio - check them
  if (particleMode==0){
    part = "K0";
    BR=0.686;
  } else if (particleMode==1){
    part = "Lambda";
    BR=0.639;
  } else if (particleMode==2){
    part = "AntiLambda";
    BR=0.639;
  } else if (particleMode==3){
    part = "Lambda+antiLambda";
    BR=0.639;
  } else if (particleMode==4){
    part = "Xi";
    BR=1.0; // Should be Lam ->p pi Br
  } else if (particleMode ==6){
    part = "Omega";
    BR=1.0; // Should be Lam->p pi  * Om -> K Lam
  }
  


  TString title[1]={"MinimumBias"}; // Not used?
  //Make 2D projections from the 3D histogram
  //Minbias (i.e. everything)
  
  // 1st argument is initial size but array will expand as necessary
  // 2nd arg means can count from 1!
  TObjArray *controllerArray = new TObjArray(10,1); 
  
  //Here probably need switch-case, depending on mult bin and particle
  // LoPt, HiPt, polynomial order, rebinning factor

  /////  LAMBDA and LAMBDA+ANTILAMBDA (Combination)



  if(particleMode == 1 || particleMode == 3 || particleMode == 2){ // Lambda or Lambda+Anti-Lambda
    if(ihist ==0){
      //           controllerArray->AddLast(new AliMassFitControl(0.2,0.3, 2,2, 1.095,1.17)); //1
      //	     controllerArray->AddLast(new AliMassFitControl(0.3,0.4, 2,2, 1.095,1.17)); //2
      //      	     controllerArray->AddLast(new AliMassFitControl(0.4,0.5, 2,2, 1.095,1.17)); //3
      //	     controllerArray->AddLast(new AliMassFitControl(0.5,0.6, 2,2, 1.095,1.17)); //4
      controllerArray->AddLast(new AliMassFitControl(0.6,0.7, 2,2, 1.105,1.14)); //5
      controllerArray->AddLast(new AliMassFitControl(0.7,0.8, 2,2, 1.095,1.14)); //6
      controllerArray->AddLast(new AliMassFitControl(0.8,0.9, 2,2, 1.095,1.14)); //7
      controllerArray->AddLast(new AliMassFitControl(0.9,1.0, 2,2, 1.095,1.14)); //1
      controllerArray->AddLast(new AliMassFitControl(1.0,1.1, 2,2, 1.099,1.14)); //2
      controllerArray->AddLast(new AliMassFitControl(1.1,1.2, 2,2, 1.099,1.14)); //3
      controllerArray->AddLast(new AliMassFitControl(1.2,1.3, 2,2, 1.099,1.14)); //4 
      controllerArray->AddLast(new AliMassFitControl(1.3,1.4, 2,2, 1.099,1.14)); //5
      controllerArray->AddLast(new AliMassFitControl(1.4,1.5, 2,2, 1.095,1.14)); //6
      controllerArray->AddLast(new AliMassFitControl(1.5,1.6, 2,2, 1.095,1.15)); //7
      controllerArray->AddLast(new AliMassFitControl(1.6,1.7, 2,2, 1.095,1.145)); //8
      controllerArray->AddLast(new AliMassFitControl(1.7,1.8, 2,2, 1.095,1.15)); //9
      controllerArray->AddLast(new AliMassFitControl(1.8,1.9, 2,2, 1.095,1.14)); //10
      controllerArray->AddLast(new AliMassFitControl(1.9,2.0, 2,2, 1.095,1.14)); //11
      controllerArray->AddLast(new AliMassFitControl(2.0,2.2, 2,2, 1.095,1.17)); //12
      controllerArray->AddLast(new AliMassFitControl(2.2,2.4, 2,2, 1.095,1.17)); //13
      controllerArray->AddLast(new AliMassFitControl(2.4,2.6, 2,2, 1.095,1.17)); //14
      controllerArray->AddLast(new AliMassFitControl(2.6,2.8, 2,2, 1.095,1.17)); //15
      controllerArray->AddLast(new AliMassFitControl(2.8,3.0, 2,2, 1.095,1.16)); //16
      controllerArray->AddLast(new AliMassFitControl(3.0,3.2, 2,2, 1.095,1.17)); //17
      controllerArray->AddLast(new AliMassFitControl(3.2,3.4, 2,2, 1.095,1.16)); //18
      controllerArray->AddLast(new AliMassFitControl(3.4,3.6, 2,2, 1.095,1.17)); //19
      controllerArray->AddLast(new AliMassFitControl(3.6,3.8, 2,2, 1.095,1.17)); //20
      controllerArray->AddLast(new AliMassFitControl(3.8,4.0, 1,2, 1.090,1.16)); //21
      controllerArray->AddLast(new AliMassFitControl(4.0,4.5, 1,2, 1.095,1.17)); //22
      controllerArray->AddLast(new AliMassFitControl(4.5,5.0, 1,2, 1.083,1.17)); //23
      controllerArray->AddLast(new AliMassFitControl(5.0,5.5, 1,2, 1.083,1.17)); //24  bin05
      controllerArray->AddLast(new AliMassFitControl(5.5,6.5, 1,2, 1.083,1.17)); //25
      controllerArray->AddLast(new AliMassFitControl(6.5,8.0, 0,2, 1.095,1.17)); //33
      //controllerArray->AddLast(new AliMassFitControl(8.0,12.0, 1,2, 1.096,1.17));//34
      

    }
    if(ihist == 1){
      cout << "histogram : " <<1<<endl;
      for(int i = 0; i<49; i++)
	//     for(int i = 0; i<33; i++) //for 05 centrality
	//	controllerArray->AddLast(new AliMassFitControl(minhist,i,i+1.0, 2,1, 0.44,0.56));
	controllerArray->AddLast(new AliMassFitControl(minhist,0.1+i*0.2,0.1+(i+1.0)*0.2, 2,1, 1.095,1.15));
    }
    if(ihist == 2){
      cout << "histogram : " <<2<<endl;
      for(int i = 0; i<49; i++)
	//      for(int i = 0; i<33; i++) //for 05 centrality
	//	controllerArray->AddLast(new AliMassFitControl(minhist,i,i+1.0, 2,1, 0.44,0.56));
	controllerArray->AddLast(new AliMassFitControl(minhist,0.1+i*0.2,0.1+(i+1.0)*0.2, 2,1, 1.095,1.17));
    }
    if(ihist == 3){
      cout << "histogram : " <<3<<endl;
      for(int i = 0; i<=40; i++)
	//	controllerArray->AddLast(new AliMassFitControl(minhist,i,i+1.0, 2,1, 0.44,0.56));
	controllerArray->AddLast(new AliMassFitControl(minhist,i,i+1.0, 2,1, 1.095,1.17));
    }
    if(ihist == 4){
      cout << "histogram : " <<4<<endl;
      for(int i = 0; i<=40; i++)
	//	controllerArray->AddLast(new AliMassFitControl(minhist,i,i+1.0, 2,1, 0.44,0.56));
	controllerArray->AddLast(new AliMassFitControl(minhist,i,i+1.0, 2,1, 1.095,1.17));
    }
    if(ihist == 5){
      cout << "histogram : " <<5<<endl;
      for(int i = 0; i<=33; i++)
	//	if(i==11 || i==15 || i==16)
	//        controllerArray->AddLast(new AliMassFitControl(minhist,0+i*0.03,(i+1)*0.03, 2,1, 0.43,0.57));
	//        else
	//controllerArray->AddLast(new AliMassFitControl(minhist,0+i*0.03,(i+1)*0.03, 2,1, 0.44,0.56));
	controllerArray->AddLast(new AliMassFitControl(minhist,0+i*0.03,(i+1)*0.03, 2,1, 1.095,1.17));   //00
    }
    if(ihist == 6){
      cout << "histogram : " <<6<<endl;
      for(int i = 1; i<40; i++){
	//   if(i==23) 
	controllerArray->AddLast(new AliMassFitControl(minhist,0.998 +i*0.00005,0.998 +(i+1)*0.00005, 2,1, 1.095,1.17));
	/*        if(  i==17 || i==18 || i==23) 
		  controllerArray->AddLast(new AliMassFitControl(minhist,0.998 +i*0.00005,0.998 +(i+1)*0.00005, 2,1, 0.446,0.55));
		  else
		  controllerArray->AddLast(new AliMassFitControl(minhist,0.998 +i*0.00005,0.998 +(i+1)*0.00005, 2,1, 0.44,0.56));
	*/        }
    }
  }
  /// ANTI LAMBDA ---->
//  else if (particleMode == 2){ // Anti-Lambdas
//    controllerArray->AddLast(new AliMassFitControl(5.0,5.5, 0,4, 1.095,1.17));
//    controllerArray->AddLast(new AliMassFitControl(5.5,6.0, 0,6, 1.085,1.17));
//  } // end if anti-Lambda
  else if (particleMode == 0){ // K0s case
    if(ihist == 0){
      //     controllerArray->AddLast(new AliMassFitControl(0.2,0.3, 2,2, 0.45,0.56)); //1
      //controllerArray->AddLast(new AliMassFitControl(0.3,0.4, 2,2, 0.45,0.56)); //2
     controllerArray->AddLast(new AliMassFitControl(0.4,0.5, 1,2, 0.45,0.6)); //3
     controllerArray->AddLast(new AliMassFitControl(0.5,0.6, 1,2, 0.45,0.59)); //4
     controllerArray->AddLast(new AliMassFitControl(0.6,0.7, 2,2, 0.45,0.55)); //5
     controllerArray->AddLast(new AliMassFitControl(0.7,0.8, 2,2, 0.45,0.55)); //6
      controllerArray->AddLast(new AliMassFitControl(0.8,0.9, 2,2, 0.44,0.55)); //7
      controllerArray->AddLast(new AliMassFitControl(0.9,1.0, 2,2, 0.443,0.55)); //8
      controllerArray->AddLast(new AliMassFitControl(1.0,1.1, 2,2, 0.443,0.55)); //9
      controllerArray->AddLast(new AliMassFitControl(1.1,1.2, 2,2, 0.443,0.56)); //10
      controllerArray->AddLast(new AliMassFitControl(1.2,1.3, 2,2, 0.44,0.55)); //11 
      controllerArray->AddLast(new AliMassFitControl(1.3,1.4, 2,2, 0.44,0.56)); //12
      controllerArray->AddLast(new AliMassFitControl(1.4,1.5, 2,2, 0.44,0.56)); //13
      controllerArray->AddLast(new AliMassFitControl(1.5,1.6, 2,2, 0.44,0.55)); //14
      controllerArray->AddLast(new AliMassFitControl(1.6,1.7, 2,2, 0.44,0.55)); //15
      controllerArray->AddLast(new AliMassFitControl(1.7,1.8, 2,2, 0.44,0.55)); //16
      controllerArray->AddLast(new AliMassFitControl(1.8,1.9, 2,2, 0.45,0.55)); //17
      controllerArray->AddLast(new AliMassFitControl(1.9,2.0, 2,2, 0.44,0.55)); //18
      controllerArray->AddLast(new AliMassFitControl(2.0,2.2, 2,2, 0.44,0.55)); //19
      controllerArray->AddLast(new AliMassFitControl(2.2,2.4, 2,2, 0.44,0.55)); //20
      controllerArray->AddLast(new AliMassFitControl(2.4,2.6, 1,2, 0.44,0.54)); //21
      controllerArray->AddLast(new AliMassFitControl(2.6,2.8, 1,2, 0.44,0.54)); //22
      controllerArray->AddLast(new AliMassFitControl(2.8,3.0, 2,2, 0.44,0.54)); //23
      controllerArray->AddLast(new AliMassFitControl(3.0,3.2, 2,2, 0.443,0.54)); //24
      controllerArray->AddLast(new AliMassFitControl(3.2,3.4, 2,2, 0.443,0.54)); //25
      controllerArray->AddLast(new AliMassFitControl(3.4,3.6, 2,2, 0.44,0.56)); //26
      controllerArray->AddLast(new AliMassFitControl(3.6,3.8, 2,2, 0.444,0.56)); //27
      controllerArray->AddLast(new AliMassFitControl(3.8,4.0, 2,2, 0.44,0.56)); //28
      controllerArray->AddLast(new AliMassFitControl(4.0,4.5, 2,2, 0.44,0.56)); //29
      controllerArray->AddLast(new AliMassFitControl(4.5,5.0, 1,2, 0.44,0.56)); //30
      controllerArray->AddLast(new AliMassFitControl(5.0,5.5, 1,2, 0.44,0.55)); //31
      controllerArray->AddLast(new AliMassFitControl(5.5,6.5, 1,2, 0.44,0.56)); //32
              controllerArray->AddLast(new AliMassFitControl(6.5,8.0, 1,2, 0.43,0.56)); //33
        controllerArray->AddLast(new AliMassFitControl(8.0,12.0, 0,4, 0.42,0.57));//34
      
    }
    if(ihist == 1){
      cout << "histogram : " <<1<<endl;
      //      for(int i = 0; i<49; i++)
      // for(int i = 0; i<=20; i++) //for pt3
	     for(int i = 0; i<33; i++) //for 05 centrality
	//	controllerArray->AddLast(new AliMassFitControl(minhist,i,i+1.0, 2,1, 0.44,0.56));
        if(i==0)
        controllerArray->AddLast(new AliMassFitControl(minhist,0.1+i*0.2,0.1+(i+1.0)*0.2, 2,1, 0.443,0.54));
        else
        if(i>=25)
        controllerArray->AddLast(new AliMassFitControl(minhist,0.1+i*0.2,0.1+(i+1.0)*0.2, 1,1, 0.443,0.53));
        else
        if(i==24)
        controllerArray->AddLast(new AliMassFitControl(minhist,0.1+i*0.2,0.1+(i+1.0)*0.2, 1,1, 0.45,0.57));
        
        else
	controllerArray->AddLast(new AliMassFitControl(minhist,0.1+i*0.2,0.1+(i+1.0)*0.2, 2,1, 0.45,0.56));
    }
    if(ihist == 2){
      cout << "histogram : " <<2<<endl;
      //      for(int i = 0; i<=20; i++){ //for pt3
	//      for(int i = 0; i<42; i++){
      for(int i = 0; i<33; i++){ //for 05 centrality
	//	controllerArray->AddLast(new AliMassFitControl(minhist,i,i+1.0, 2,1, 0.44,0.56));
        if(i<=5)controllerArray->AddLast(new AliMassFitControl(minhist,0.1+i*0.2,0.1+(i+1.0)*0.2, 2,1, 0.445,0.55));
        else
        if(i>25)
        controllerArray->AddLast(new AliMassFitControl(minhist,0.1+i*0.2,0.1+(i+1.0)*0.2, 1,1, 0.444,0.56));
        else
	  controllerArray->AddLast(new AliMassFitControl(minhist,0.1+i*0.2,0.1+(i+1.0)*0.2, 2,1, 0.445,0.56));
      }
    }
    if(ihist == 3){
      cout << "histogram : " <<3<<endl;
      //for(int i = 0; i<=28; i++) //bin05
      for(int i = 1; i<=20; i++) //for pt3
	// for(int i = 0; i<=40; i++)
	//	controllerArray->AddLast(new AliMassFitControl(minhist,i,i+1.0, 2,1, 0.44,0.56));
        if(i==2 || i==4 || i==5)
        controllerArray->AddLast(new AliMassFitControl(minhist,i,i+1.0, 2,1, 0.442,0.535));
        else if(i<2)
	controllerArray->AddLast(new AliMassFitControl(minhist,i,i+1.0, 2,1, 0.442,0.53));
        else if(i>=3)
	controllerArray->AddLast(new AliMassFitControl(minhist,i,i+1.0, 2,1, 0.443,0.55));
         
        else
	controllerArray->AddLast(new AliMassFitControl(minhist,i,i+1.0, 2,1, 0.443,0.54));
    }
    if(ihist == 4){
      cout << "histogram : " <<4<<endl;
      for(int i = 0; i<=40; i++)
	//	controllerArray->AddLast(new AliMassFitControl(minhist,i,i+1.0, 2,1, 0.44,0.56));
	if(i<4)controllerArray->AddLast(new AliMassFitControl(minhist,i,i+1.0, 1,1, 0.443,0.56));
        else if (i >=4 && i<=6)
        controllerArray->AddLast(new AliMassFitControl(minhist,i,i+1.0, 1,1, 0.443,0.547));
         else if (i==10 || i==11)
        controllerArray->AddLast(new AliMassFitControl(minhist,i,i+1.0, 1,1, 0.443,0.544));
         else if (i==26)
        controllerArray->AddLast(new AliMassFitControl(minhist,i,i+1.0, 1,1, 0.443,0.544));
         else if (i==29 || i == 34)
        controllerArray->AddLast(new AliMassFitControl(minhist,i,i+1.0, 1,1, 0.454,0.544));
        else
        controllerArray->AddLast(new AliMassFitControl(minhist,i,i+1.0, 1,1, 0.445,0.55));
    }
    if(ihist == 5){
      cout << "histogram : " <<5<<endl;
      for(int i = 0; i<25; i++) //for pt3
	//      for(int i = 0; i<=32; i++)
	//	if(i==11 || i==15 || i==16)
	//        controllerArray->AddLast(new AliMassFitControl(minhist,0+i*0.03,(i+1)*0.03, 2,1, 0.43,0.57));
	//        else
	//controllerArray->AddLast(new AliMassFitControl(minhist,0+i*0.03,(i+1)*0.03, 2,1, 0.44,0.56));
      if(i<=5) controllerArray->AddLast(new AliMassFitControl(minhist,0+i*0.03,(i+1)*0.03, 2,1, 0.443,0.55));   //00
      else
      if(i>=6 && i<=9) controllerArray->AddLast(new AliMassFitControl(minhist,0+i*0.03,(i+1)*0.03, 2,1, 0.44,0.55));   //00
      else
      if(i>=21) controllerArray->AddLast(new AliMassFitControl(minhist,0+i*0.03,(i+1)*0.03, 1,1, 0.44,0.55));   //00
      else
	controllerArray->AddLast(new AliMassFitControl(minhist,0+i*0.03,(i+1)*0.03, 2,1, 0.43,0.56));

    }
    if(ihist == 6){
      cout << "histogram : " <<6<<endl;
            for(int i = 18; i<40; i++){ //for pt3
      //  for(int i = 1; i<40; i++){
	//   if(i==23) 
	        if(  i==15 || i==0) 
		  controllerArray->AddLast(new AliMassFitControl(minhist,0.998 +i*0.00005,0.998 +(i+1)*0.00005, 1,1, 0.444,0.56));
                else if(  i==10 || i==11) 
		  controllerArray->AddLast(new AliMassFitControl(minhist,0.998 +i*0.00005,0.998 +(i+1)*0.00005, 1,1, 0.444,0.56));
	        else if(  i==20 || i==21) 
		  controllerArray->AddLast(new AliMassFitControl(minhist,0.998 +i*0.00005,0.998 +(i+1)*0.00005, 1,1, 0.445,0.56));
		else if (i<=22)
		  controllerArray->AddLast(new AliMassFitControl(minhist,0.998 +i*0.00005,0.998 +(i+1)*0.00005, 1,1, 0.445,0.56));
		else if (i>22 && i<=24)
		  controllerArray->AddLast(new AliMassFitControl(minhist,0.998 +i*0.00005,0.998 +(i+1)*0.00005, 1,1, 0.445,0.555));
		else if (i>24 && i<29)
		  controllerArray->AddLast(new AliMassFitControl(minhist,0.998 +i*0.00005,0.998 +(i+1)*0.00005, 1,1, 0.44,0.53));
		else if (i>29)
		  controllerArray->AddLast(new AliMassFitControl(minhist,0.998 +i*0.00005,0.998 +(i+1)*0.00005, 1,1, 0.445,0.54));
			        else if( i>31 ) 
		controllerArray->AddLast(new AliMassFitControl(minhist,0.998 +i*0.00005,0.998 +(i+1)*0.00005, 2,1, 0.44,0.54));

                else
                  	controllerArray->AddLast(new AliMassFitControl(minhist,0.998 +i*0.00005,0.998 +(i+1)*0.00005, 2,1, 0.445,0.55));
			 }
    }
  }  else if (particleMode == 4) { //Xi case
    //controllerArray->AddLast(new AliMassFitControl(0.5,0.7, 1,1, 1.28,1.45)); //signal not visible with 
    controllerArray->AddLast(new AliMassFitControl(0.7,0.9, 1,1, 1.285,1.45));
    controllerArray->AddLast(new AliMassFitControl(0.9,1.1, 1,1, 1.285,1.45));
    controllerArray->AddLast(new AliMassFitControl(1.1,1.3, 1,1, 1.27,1.45));
    controllerArray->AddLast(new AliMassFitControl(1.3,1.5, 1,1, 1.27,1.45));
    controllerArray->AddLast(new AliMassFitControl(1.5,1.7, 1,1, 1.27,1.45));
    controllerArray->AddLast(new AliMassFitControl(1.7,2.0, 1,1, 1.27,1.45));
    controllerArray->AddLast(new AliMassFitControl(2.0,2.5, 1,1, 1.27,1.44));
    controllerArray->AddLast(new AliMassFitControl(2.5,3.0, 1,1, 1.27,1.45));
    controllerArray->AddLast(new AliMassFitControl(3.0,3.5, 1,1, 1.27,1.44));
    controllerArray->AddLast(new AliMassFitControl(3.5,4.0, 1,1, 1.27,1.45));
    controllerArray->AddLast(new AliMassFitControl(4.0,4.5, 2,1, 1.27,1.44));
    controllerArray->AddLast(new AliMassFitControl(4.5,5.0, 2,1, 1.27,1.45));
  }
  controllerArray->Sort();

  // Make the proper label to pass in for use in labelling the diagnostics
  // saved canvas files and histos
  Char_t fulllabel[80];
  //sprintf(fulllabel,"%s%s",title[MultBin].Data(),label);
  sprintf(fulllabel,"%s",label);

  //Slice up projection into mass histograms to extract yield
  TH1F* hYield = (TH1F*)PtMassAna2(hParMass,particleMode,ihist,controllerArray->GetEntries(),controllerArray,fulllabel);


  // CORRECTIONS Nev -comming from "FitSpectrum.C" (events with Zvertex < |10cm|) 
  hYield->Scale(1.0/Nev); //Divide by the number of events

  //hYield->Scale(Veff[MultBin]); //Multiply by the vertex efficiency effectively increaing number of events 
  //(since Veff<1) therefore decreases yield
  //hYield->Scale(1.0/BR);  //Divide by branching ratio (again increases yield since BR<1)
  //hYield->Scale(1.0/(2*TMath::Pi())); // Always plot 1/2pi ...

  Char_t yieldTitle[80];
  //sprintf(yieldTitle,"Uncorrected %s yield: %s",part,title[MultBin].Data());
  sprintf(yieldTitle,"Uncorrected %s yield",part);
  hYield->SetTitle(yieldTitle);
  if(ihist == 0)hYield->SetXTitle("p_{t} / [GeV/c]");
  if(ihist == 1)hYield->SetXTitle("DCA / [cm]");
  if(ihist == 2)hYield->SetXTitle("DCA / [cm]");
  if(ihist == 3)hYield->SetXTitle("Radius / [cm]");
  if(ihist == 4)hYield->SetXTitle("Decay Length / [cm]");
  if(ihist == 5)hYield->SetXTitle("V0 Daughters / [cm]");
  if(ihist == 6)hYield->SetXTitle("Cos of pointing angle");
  hYield->SetYTitle("1/Nev.dN/dp_{t}");

  // Create plots

  Char_t fileNameBase[80];
  if(ihist == 0)sprintf(fileNameBase,"Masses%s",part);
  else sprintf(fileNameBase,"Masses%s",label);
  Char_t fileNamePng[80];
  sprintf(fileNamePng,"%s.png",fileNameBase);
  Char_t fileNameEps[80];
  sprintf(fileNameEps,"%s.eps",fileNameBase);
  Char_t fileNamePdf[80];
  sprintf(fileNamePdf,"%s.pdf",fileNameBase);

  c1->SaveAs(fileNamePng);
  c1->SaveAs(fileNameEps);
  c1->SaveAs(fileNamePdf);

  //c1->Clear();

  //c1->SetLogy();
  TCanvas *cYield = new TCanvas("Yield","Corrected Yield",600,400);
  cYield->cd();
  //cYield->SetLogy();
  hYield->SetStats(kFALSE);
  hYield->Draw();
  // cYield->SetLogy();
  cYield->Update();

  //  hRC_MB->SetMarkerStyle(20);
  //  hRC_MB->SetMarkerColor(4);
  //  hRC_MB->Scale(NBinMB/NBin3);
  Char_t fnametext[80];
  if(ihist == 0)sprintf(fnametext,"Yield%s",part);
  else sprintf(fnametext,"Yield%s",label);
  Char_t fnamePng[80];
  sprintf(fnamePng,"%s.png",fnametext);
  c1->SaveAs(fnamePng);
  Char_t fnameEps[80];
  sprintf(fnameEps,"%s.eps",fnametext);
  c1->SaveAs(fnameEps);

  // This section for yield scaled by number of binary collisions.
  // Could add array of values and do scaling according to 'MultBin' index
  TH1F* hScYield = hYield->Clone("ScYield");
  Char_t scYieldTitle[80];
  //sprintf(scYieldTitle,"<N_{bin}> Scaled %s",hYield->GetTitle());
  //hScYield->SetTitle(scYieldTitle);
  //SCALING for scaled yield only  divide by mean Nbin (scaled yield is therefore smaller)
  //hScYield->Scale(1/NBin[MultBin]);

  Char_t fnameRoot[80];
  sprintf(fnameRoot,"%s.root",fnametext);
  TFile *YieldFile = new TFile(fnameRoot,"RECREATE");
  hYield->Write();
  hScYield->Write();
  YieldFile->Close();

  controllerArray->Delete();
}
