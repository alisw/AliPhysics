//Macros to calculate PHOS alignment matrixes from survay
//To run select module number (mod=1,2,3,4) and uncomment correcponding 
//pairs of matrixes: survay points in crystal plane and in ALICE global system
//These data one takes from syrvay results.
// Note that for Mod.4 there was now survay repers in crystal surface frame.
// So we choose module with closest ditances between repers (mod2==mod4) 
// and assume that mod4 has similr reper coordinates wrt crytal surface.
//

TGeoHMatrix *mPHOS[6] ;
Int_t mod=1 ; //Which module minimize


//Module 1, points 500x survey
  const Float_t srv5000[3]={-217.61,-386.44,-77.98} ;
  const Float_t srv5001[3]={ -79.32,-439.05,-77.54} ;
  const Float_t srv5002[3]={ -79.40,-439.29, 80.80} ;
  const Float_t srv5003[3]={-217.64,-386.49, 80.04} ;
  
//Module 1, points 500x wrt crystals (copied from Mod3)
  const Float_t p5000[3]={ -3.08,28.97,-16.43} ;
  const Float_t p5001[3]={145.94,28.90,-16.03} ;
  const Float_t p5002[3]={145.58,29.01,142.17} ;
  const Float_t p5003[3]={ -3.58,28.88,141.88} ;

  
  


/*
//Module 2, points 500x survey
  const Float_t srv5000[3]={-74.67,-429.97,-90.95} ;
  const Float_t srv5001[3]={ 72.39,-430.96,-91.25} ;
  const Float_t srv5002[3]={ 64.43,-431.06, 92.01} ;
  const Float_t srv5003[3]={-70.56,-430.27, 92.43} ;
*/ 
/*  
//Module 3, points 500x survey
  const Float_t srv5000[3]={ 77.95,-430.77,-78.85} ;
  const Float_t srv5001[3]={217.78,-380.47,-78.76} ;
  const Float_t srv5002[3]={217.86,-380.25, 79.28} ;
  const Float_t srv5003[3]={ 77.57,-430.93, 79.58} ;
*/
/*  
//Module 4, points 500x survey
  const Float_t srv5000[3]={222.58,-377.81,-78.96} ;
  const Float_t srv5001[3]={333.99,-284.63,-78.21} ;
  const Float_t srv5002[3]={335.38,-283.62, 80.39} ;
  const Float_t srv5003[3]={220.94,-379.38, 79.89} ;
*/  
  

/*
  //Module 2 reference points on crystal surface (in cm)
  const Float_t p5000[3]={0.0000,0.0000,  0.00} ;
  const Float_t p5001[3]={0.0000,0.0000,122.62} ;
  const Float_t p5002[3]={142.42,0.0000,  1.95} ;
  const Float_t p5003[3]={142.49,0.0001,124.69} ;
  
  //Module 2 survey: position of these points in global ALICE system
  const Float_t srv5000[3]={-71.06,-459.00,-61.82} ;
  const Float_t srv5001[3]={-71.16,-459.15, 61.06} ;
  const Float_t srv5002[3]={ 71.65,-459.99,-59.76} ;
  const Float_t srv5003[3]={ 71.63,-460.13, 63.24} ;
*/
/*  
  //Module 3 reference points on crystal surface
  const Float_t p5000[3]={0.0000,0.0000,0.0000 } ;
  const Float_t p5001[3]={0.0000,0.0000,122.63 } ;
  const Float_t p5002[3]={142.53,0.0000,  2.37 } ;
  const Float_t p5003[3]={142.28, -0.11,125.06 } ;
  
  //Module 3 survey 
  const Float_t srv5000[3]={ 90.51,-457.01,-62.35 } ;
  const Float_t srv5001[3]={ 90.77,-456.86, 60.19 } ;
  const Float_t srv5002[3]={224.51,-408.73,-60.33 } ;
  const Float_t srv5003[3]={224.57,-408.77, 62.27 } ;

*/
/*
  //Module 4 reference points on crystal surface
  const Float_t p5000[3]={0.0000,0.0000,0.0000 } ;
  const Float_t p5001[3]={0.0000,0.0000,122.98 } ;
  const Float_t p5002[3]={142.39,0.0000,  2.05 } ;
  const Float_t p5003[3]={142.41,0.0019,124.82 } ;

  //Module 4 survey 
  const Float_t srv5000[3]={244.57,-399.74,-62.67 } ;
  const Float_t srv5001[3]={244.10,-400.79, 60.38 } ;
  const Float_t srv5002[3]={353.95,-308.47,-59.42 } ;
  const Float_t srv5003[3]={353.38,-309.36, 63.42 } ;
*/
  
  
  
  

Double_t DistV(TVector3 &a, Float_t * b){
  
 return sqrt((a[0]-b[0])*(a[0]-b[0])+ (a[1]-b[1])*(a[1]-b[1])+ (a[2]-b[2])*(a[2]-b[2])   ) ; 
  
}

Double_t Dist(Float_t *a, Float_t *b){
  
 return sqrt((a[0]-b[0])*(a[0]-b[0])+ (a[1]-b[1])*(a[1]-b[1])+ (a[2]-b[2])*(a[2]-b[2])   ) ; 
  
}
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  //Function used in Minuit minimization:
  //Minimize distances between measred reper coordinates
  //and calculate after PHOS rotaion/shifts.
  //Uses dedicated function AliPHOSGeoUtils::TestSurvey()
  
   Double_t chisq = 0;

  //Create current rotation matrix   
  TGeoRotation r1;
  r1.SetAngles(par[0],par[1],par[2]) ;
  mPHOS[mod]->SetRotation(r1.GetRotationMatrix()) ;
  mPHOS[mod]->SetDx(par[3]) ;
  mPHOS[mod]->SetDy(par[4]) ;
  mPHOS[mod]->SetDz(par[5]) ;
  
  AliPHOSGeoUtils * geom = new AliPHOSGeoUtils("PHOS") ;     
  geom->SetMisalMatrix(mPHOS[mod],mod-1) ;

  TVector3 globaPos ;
  //Compare ditances for all 4 points  
  geom->TestSurvey(mod,p5000,globaPos) ; 
  chisq+=DistV(globaPos,srv5000) ;
 
  geom->TestSurvey(mod,p5001,globaPos) ; 
  chisq+=DistV(globaPos,srv5001) ;

  geom->TestSurvey(mod,p5002,globaPos) ; 
  chisq+=DistV(globaPos,srv5002) ;

  geom->TestSurvey(mod,p5003,globaPos) ; 
  chisq+=DistV(globaPos,srv5003) ;
 
  delete geom ;
  f = chisq;
}


void CheckPHOSAlignment(){
  //Checks if PHOS alignment agree with results of photogrammetry
  
  
  mPHOS[mod] = new TGeoHMatrix() ;
  
  //Initial conditions for minimization, close to ideal PHOS positions
  TGeoRotation r1;
  r1.SetAngles(0,90+(mod-2)*20,0) ;
  r1.Print() ;
  
  mPHOS[mod]->SetRotation(r1.GetRotationMatrix()) ;
  mPHOS[mod]->SetDx(-480.*TMath::Cos((210+mod*20)*TMath::DegToRad())) ;
  mPHOS[mod]->SetDy(-480.*TMath::Sin((210+mod*20)*TMath::DegToRad())) ;
  mPHOS[mod]->SetDz(0.) ;
  
  
  
  //Each module has 4 reference points
  //Define them w.r.t. center of top right crystal (see more details in EDH comenents
  // 992344, 1002387, 1012391.

  //Print distances between repers in current module
//printf("%6.2f,\t %6.2f,\t  %6.2f,\t %6.2f\n",Dist(srv5000,srv5001),Dist(srv5000,srv5002),Dist(srv5000,srv5003),Dist(srv5001,srv5003)) ;
//return ;

  
  
  //Module 3
  const Float_t m3p5000[3]={-30.76,  191.15, -164.32} ;
  const Float_t m3p5001[3]={1459.44, 190.47, -160.25} ;
  const Float_t m3p5002[3]={1455.82, 191.57, 1421.70} ;
  const Float_t m3p5003[3]={-35.83,  190.32, 1418.78} ;
  
  //Module 3 survay (in meters!)  
  const Float_t m3srv5000[3]={0.7795,-4.3077,-0.7885} ;
  const Float_t m3srv5001[3]={2.1778,-3.8047,-0.7876} ;
  const Float_t m3srv5002[3]={2.1786,-3.8025,0.7928} ;
  const Float_t m3srv5003[3]={0.7757,-4.3093,0.7958} ;

  
  //Module 4
  const Float_t m4p5000[3]={-34.59,  210.66, -165.23} ;
  const Float_t m4p5001[3]={1456.73, 210.71, -164.91} ;
  const Float_t m4p5002[3]={1454.76, 213.62, 1418.01} ;
  const Float_t m4p5003[3]={-36.60,  210.14, 1419.84} ;
  
  
  const Float_t m4srv5000[3]={2.2258,-3.7781,-0.7896} ;
  const Float_t m4srv5001[3]={3.3399,-2.8463,-0.7821} ;
  const Float_t m4srv5002[3]={3.3538,-2.8362, 0.8039} ;
  const Float_t m4srv5003[3]={2.2094,-3.7938, 0.7989} ;

  
  
  
  // prepare minuit
  const Double_t nPar = 6; // number of fit parameters 
  TMinuit* gMinuit = new TMinuit(nPar);
  gMinuit->SetFCN(fcn);
  gMinuit->SetPrintLevel(1); // -1 quiet, 0 normal, 1 verbose    

  Double_t arglist[10];
  Int_t ierflg = 0;

  // start minimization
  Double_t chi2 = 0;

  ierflg = 0;
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
  gMinuit->SetMaxIterations(5000);

  // parameters: parameter no., name, start value, step size, parameter range min., parameter range max., status
  gMinuit->mnparm(0,"A1",0., 1., 0, 0, ierflg);
  gMinuit->mnparm(1,"A2",90., 1., 0, 0, ierflg);
  gMinuit->mnparm(2,"A3",0., 1., 0, 0, ierflg);
  gMinuit->mnparm(3,"dx",0., 1., 0, 0, ierflg);
  gMinuit->mnparm(4,"dy",-460., 1., 0, 0, ierflg);
  gMinuit->mnparm(5,"dz",0., 1., 0, 0, ierflg);

  // Now ready for minimization step
  arglist[0] = 5000;
  arglist[1] = 1.;
  gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  
  AliPHOSGeometry *  geom = AliPHOSGeometry::GetInstance("PHOS_mod1234") ;
  geom->SetMisalMatrix(mPHOS[mod],mod-1) ;
  
  
  
  TVector3 globaPos ;
  

  geom->TestSurvey(mod,p5000,globaPos) ; 
  printf("P5000: (%f, %f, %f) \n",globaPos.X(),globaPos.Y(),globaPos.Z()) ;
  printf("  srv: (%f, %f, %f) \n",srv5000[0],srv5000[1],srv5000[2]) ;
  geom->TestSurvey(mod,p5001,globaPos) ; 
  printf("P5001: (%f, %f, %f) \n",globaPos.X(),globaPos.Y(),globaPos.Z()) ;
  printf("  srv: (%f, %f, %f) \n",srv5001[0],srv5001[1],srv5001[2]) ;
  geom->TestSurvey(mod,p5002,globaPos) ; 
  printf("P5002: (%f, %f, %f) \n",globaPos.X(),globaPos.Y(),globaPos.Z()) ;
  printf("  srv: (%f, %f, %f) \n",srv5002[0],srv5002[1],srv5002[2]) ;
  geom->TestSurvey(mod,p5003,globaPos) ; 
  printf("P5003: (%f, %f, %f) \n",globaPos.X(),globaPos.Y(),globaPos.Z()) ;
  printf("  srv: (%f, %f, %f) \n",srv5003[0],srv5003[1],srv5003[2]) ;
 
  mPHOS[mod]->Print() ;
}



