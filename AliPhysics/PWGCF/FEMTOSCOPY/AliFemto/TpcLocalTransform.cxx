
#include "AliFmThreeVectorD.h"
// - not needed (which is good because we want no STAR-dependence! 25apr2006) - #include "StEvent/StTpcHit.h"
#include "TMath.h"


//__________________________________



int TpcLocalTransform(AliFmThreeVectorD& aPoint, int& aSector, int& aRow, 
		      float& aU, double& aPhi){
  static int tNPadAtRow[45]={
  88,96,104,112,118,126,134,142,150,158,166,174,182,
  98,100,102,104,106,106,108,110,112,112,114,116,118,120,122,122,
  124,126,128,128,130,132,134,136,138,138,140,142,144,144,144,144};
  static double tSectToPhi[24]={2.,1.,0.,11.,10.,9.,8. ,7. ,6.,5.,4.,3.,
				4.,5.,6., 7., 8.,9.,10.,11.,0.,1.,2.,3.};
  //static double tPhiToSect[24]={2.,1.,0.,11.,10.,9.,8. ,7. ,6.,5.,4.,3.,
	//			4.,5.,6., 7., 8.,9.,10.,11.,0.,1.,2.,3.};
  static double tPadWidthInner = 0.335;
  static double tPadWidthOuter = 0.67;

  static double tPi = TMath::Pi();
  // --- find sector number
  aPhi = aPoint.Phi();
  if(aPhi<0.) aPhi+=(2*tPi);
  aPhi += tPi/12.;
  if(aPhi>2*tPi) aPhi-=2*tPi;
  int tiPhi = (int) (aPhi/tPi*6.);
  if(aPoint.z()<0) {
    aSector = (tiPhi<3)? 3-tiPhi : 15-tiPhi;
  }
  else{
    aSector = (tiPhi<4)? 21+tiPhi : 9+tiPhi;
  }
  aPhi = tSectToPhi[aSector-1]*tPi/6.;
  //if((fabs(aPhi-aPoint.phi())>(tPi/12)){
  //cout << "Sector missmatch " << aPhi << " " << aPoint.phi() << " "
  // << aSector << endl;
  //}

  // --- calculate local coordinate
  float tR = aPoint.x()*cos(aPhi)+aPoint.y()*sin(aPhi);
  aU =      -aPoint.x()*sin(aPhi)+aPoint.y()*cos(aPhi);

  // --- find pad row 
  if(tR<57.6) {
    aRow = 0;
    return 1;
  }
  float radmax = 62.4;
  float spacing= 4.8;
  aRow=1;
  while(tR>radmax && aRow<46){
    aRow++;
    if(aRow==8){
      radmax = 96.2;
      spacing = 5.2;
    }
    else{
      if (aRow==13){
	radmax = 126.195; // lots of stuf in row 13!
	spacing = 2.0;
      }
      else{
	radmax+=spacing;
      }
    }
  }
  if(aRow>45){
    //cout << "No pad row " << tR << endl;
    return 2;
  }
  
  // --- Check if u (=aU) inbound
  double tPadWidth = aRow<14? tPadWidthInner : tPadWidthOuter;
  if(fabs(aU) > tNPadAtRow[aRow-1]*tPadWidth/2.){
    return 3;
  }

  return 0;
}






