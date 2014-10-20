#include "PhotosUtilities.h"
#include <cstdlib>
#include <cstdio>
using std::max;

namespace Photospp
{

namespace PhotosUtilities
{


void fill_val(int beg, int end, double* array, double value) 
{
  for (int i = beg; i < end; i++)
    array[i] = value;
}


//----------------------------------------------------------------------
//
//    PHOEPS:   PHOeps vector product (normalized to unity)
//
//    Purpose:  calculates vector product, then normalizes its length.
//              used to generate orthogonal vectors, i.e. to
//              generate polarimetric vectors for photons.
//
//    Input Parameters:  VEC1,VEC2 - input 4-vectors
//                                          
//    Output Parameters: EPS - normalized 4-vector, orthogonal to
//                             VEC1 and VEC2
//
//    Author(s):  Z. Was, P.Golonka               Created at:  19/01/05
//                                                Last Update: 10/06/13
//
//----------------------------------------------------------------------

void PHOEPS(double vec1[4], double vec2[4], double eps[4]){
  double xn;
  int j=1;  // convention of indices of Riemann space must be preserved.

  eps[1-j]=vec1[2-j]*vec2[3-j] - vec1[3-j]*vec2[2-j];
  eps[2-j]=vec1[3-j]*vec2[1-j] - vec1[1-j]*vec2[3-j];      
  eps[3-j]=vec1[1-j]*vec2[2-j] - vec1[2-j]*vec2[1-j];
  eps[4-j]=0.0;
      
  xn=sqrt( eps[1-j]*eps[1-j] + eps[2-j]*eps[2-j] + eps[3-j]*eps[3-j]);
      
  eps[1-j]=eps[1-j]/xn;
  eps[2-j]=eps[2-j]/xn;
  eps[3-j]=eps[3-j]/xn;

}



//----------------------------------------------------------------------
//
//    PHOTOS:   PHOton radiation  in decays function for SPIn determina-
//              tion
//
//    Purpose:  Calculate  the spin  of particle  with  code IDHEP.  The
//              code  of the particle  is  defined  by the Particle Data
//              Group in Phys. Lett. B204 (1988) 1.
//
//    Input Parameter:   IDHEP
//
//    Output Parameter:  Funtion  value = spin  of  particle  with  code
//                       IDHEP
//
//    Author(s):  E. Barberio and B. van Eijk     Created at:  29/11/89
//                                                Last update: 10/06/13
//
//----------------------------------------------------------------------
double PHOSPI(int idhep){
  static double SPIN[100] = { 0 };
  static int j=0;  
  //--
  //--   Array 'SPIN' contains the spin  of  the first 100 particles accor-
  //--   ding to the PDG particle code...

  if(j==0) // initialization
    {   
      j=1;
      fill_val(0 ,  8, SPIN, 0.5);
      fill_val(8 ,  9, SPIN, 1.0);
      fill_val(9 , 10, SPIN, 0.0);
      fill_val(10, 18, SPIN, 0.5);
      fill_val(18, 20, SPIN, 0.0);
      fill_val(20, 24, SPIN, 1.0);
      fill_val(24,100, SPIN, 0.0);
    }

  int idabs=abs(idhep);
  //--
  //--   Spin of quark, lepton, boson etc....
  if (idabs-1<100) return SPIN[idabs-1];

  //--   ...other particles, however...
  double xx=((idabs % 10)-1.0)/2.0;
  //--
  //--   ...K_short and K_long are special !!
  xx=max(xx,0.0);
  return xx;
}

//----------------------------------------------------------------------
//
//    PHOTOS:   PHOton radiation in decays CHArge determination
//
//    Purpose:  Calculate the charge  of particle  with code IDHEP.  The
//              code  of the  particle  is  defined by the Particle Data
//              Group in Phys. Lett. B204 (1988) 1.
//
//    Input Parameter:   IDHEP
//
//    Output Parameter:  Funtion value = charge  of  particle  with code
//                       IDHEP
//
//    Author(s):  E. Barberio and B. van Eijk     Created at:  29/11/89
//                                                Last update: 11/06/13
//
//----------------------------------------------------------------------
double PHOCHA(int idhep){
  static double CHARGE[101] = { 0 };
  static int j=0;  
  //--
  //--   Array 'SPIN' contains the spin  of  the first 100 particles accor-
  //--   ding to the PDG particle code...

  if(j==0) // initialization
    {   
      j=1;
      fill_val(0 ,  1, CHARGE, 0.0         );
      fill_val(1 ,  2, CHARGE,-0.3333333333);
      fill_val(2 ,  3, CHARGE, 0.6666666667);
      fill_val(3 ,  4, CHARGE,-0.3333333333);
      fill_val(4 ,  5, CHARGE, 0.6666666667);
      fill_val(5 ,  6, CHARGE,-0.3333333333);
      fill_val(6 ,  7, CHARGE, 0.6666666667);
      fill_val(7 ,  8, CHARGE,-0.3333333333);
      fill_val(8 ,  9, CHARGE, 0.6666666667);
      fill_val(9 , 11, CHARGE, 0.0         );
      fill_val(11 ,12, CHARGE,-1.0         );
      fill_val(12 ,13, CHARGE, 0.0         );
      fill_val(13 ,14, CHARGE,-1.0         );
      fill_val(14, 15, CHARGE, 0.0         );
      fill_val(15 ,16, CHARGE,-1.0         );
      fill_val(16, 17, CHARGE, 0.0         );
      fill_val(17 ,18, CHARGE,-1.0         );
      fill_val(18, 24, CHARGE, 0.0         );
      fill_val(24, 25, CHARGE, 1.0         );
      fill_val(25, 37, CHARGE, 0.0         );
      fill_val(37, 38, CHARGE, 1.0         );
      fill_val(38,101, CHARGE, 0.0         );
    }

  int idabs=abs(idhep);
  double phoch=0.0;

  //--
  //--   Charge of quark, lepton, boson etc....
  if (idabs<=100) phoch=CHARGE[idabs];
  else {
    int Q3= idabs/1000 % 10;
    int Q2= idabs/100  % 10;
    int Q1= idabs/10   % 10;
    if (Q3==0){
      //--
      //-- ...meson...
      if(Q2 % 2==0) phoch=CHARGE[Q2]-CHARGE[Q1];
      else          phoch=CHARGE[Q1]-CHARGE[Q2];
    }
    else{
      //--
      //--   ...diquarks or baryon.
      phoch=CHARGE[Q1]+CHARGE[Q2]+CHARGE[Q3];
    }
  }
  //--
  //--   Find the sign of the charge...
  if (idhep<0.0) phoch=-phoch;
  if (phoch*phoch<0.000001) phoch=0.0;
  
  return phoch;
}




//----------------------------------------------------------------------
//
//    PHOTOS:   PHOton radiation in decays calculation of TRIangle fie
//
//    Purpose:  Calculation of triangle function for phase space.
//
//    Input Parameters:  A, B, C (Virtual) particle masses.
//
//    Output Parameter:  Function value =
//                       SQRT(LAMBDA(A**2,B**2,C**2))/(2*A)
//
//    Author(s):  B. van Eijk                     Created at:  15/11/89
//                                                Last Update: 12/06/13
//
//----------------------------------------------------------------------
double PHOTRI(double A,double B,double C){
  double DA,DB,DC,DAPB,DAMB,DTRIAN;
  DA=A;
  DB=B;
  DC=C;
  DAPB=DA+DB;
  DAMB=DA-DB;
  DTRIAN=sqrt((DAMB-DC)*(DAPB+DC)*(DAMB+DC)*(DAPB-DC));
  return DTRIAN/(DA+DA);
}
//----------------------------------------------------------------------
//
//    PHOTOS:   PHOton radiation in decays calculation of ANgle '1'
//
//    Purpose:  Calculate angle from X and Y
//
//    Input Parameters:  X, Y
//
//    Output Parameter:  Function value
//
//    Author(s):  S. Jadach                       Created at:  01/01/89
//                B. van Eijk                     Last Update: 12/06/13
//
//----------------------------------------------------------------------
double PHOAN1(double X,double Y){

  double phoan1 = 0.0;
  
  static double PI=3.14159265358979324, TWOPI=6.28318530717958648;
 
  if (fabs(Y)<fabs(X)){
    phoan1=atan(fabs(Y/X));
    if (X<0.0) phoan1=PI-phoan1;
  }
  else phoan1=acos(X/sqrt(X*X+Y*Y));
  //
  if (Y<0.0) phoan1=TWOPI-phoan1;
  return phoan1;
 
}

//----------------------------------------------------------------------
//
//    PHOTOS:   PHOton radiation in decays calculation of ANgle '2'
//
//    Purpose:  Calculate angle from X and Y
//
//    Input Parameters:  X, Y
//
//    Output Parameter:  Function value
//
//    Author(s):  S. Jadach                       Created at:  01/01/89
//                B. van Eijk                     Last Update: 12/06/13
//
//----------------------------------------------------------------------
double PHOAN2(double X,double Y){

  double phoan2 = 0.0;

  static double PI=3.14159265358979324; //, TWOPI=6.28318530717958648;

  if (fabs(Y)<fabs(X)){
    phoan2=atan(fabs(Y/X));
    if (X<0.0) phoan2=PI-phoan2;
  }
  else phoan2=acos(X/sqrt(X*X+Y*Y));
  return phoan2;
}

//----------------------------------------------------------------------
//
//    PHOTOS:   PHOton radiation in decays ROtation routine '2'
//
//    Purpose:  Rotate  x and z components  of vector PVEC  around angle
//              'ANGLE'.
//
//    Input Parameters:  ANGLE, PVEC
//
//    Output Parameter:  PVEC
//
//    Author(s):  S. Jadach                       Created at:  01/01/89
//                B. van Eijk                     Last Update: 12/06/13
//
//----------------------------------------------------------------------
void PHORO2(double ANGLE,double PVEC[4]){
  int j=1;  // convention of indices of Riemann space must be preserved.

  double CS,SN;
  CS= cos(ANGLE)*PVEC[1-j]+sin(ANGLE)*PVEC[3-j];
  SN=-sin(ANGLE)*PVEC[1-j]+cos(ANGLE)*PVEC[3-j];
  PVEC[1-j]=CS;
  PVEC[3-j]=SN;
}

//----------------------------------------------------------------------
//
//    PHOTOS:   PHOton radiation in decays ROtation routine '3'
//
//    Purpose:  Rotate  x and y components  of vector PVEC  around angle
//              'ANGLE'.
//
//    Input Parameters:  ANGLE, PVEC
//
//    Output Parameter:  PVEC
//
//    Author(s):  S. Jadach     RO                 Created at:  01/01/89
//                B. van Eijk                     Last Update: 12/06/13
//
//----------------------------------------------------------------------
void PHORO3(double ANGLE,double PVEC[4]){
  int j=1;  // convention of indices of Riemann space must be preserved.
  double CS,SN;
  CS=cos(ANGLE)*PVEC[1-j]-sin(ANGLE)*PVEC[2-j];
  SN=sin(ANGLE)*PVEC[1-j]+cos(ANGLE)*PVEC[2-j];
  PVEC[1-j]=CS;
  PVEC[2-j]=SN;
}

//----------------------------------------------------------------------
//
//
//    PHOB:     PHotosBoost
//
//    Purpose:  Boosts VEC to (MODE=1)  rest frame of PBOOS1;  
//              or back (MODE=1)
//
//    Input Parameters:   MODE,PBOOS1,VEC
//
//    Output Parameters:  VEC
//
//    Author(s):                                  Created at:  08/12/05
//                Z. Was                          Last Update: 13/06/13
//
//----------------------------------------------------------------------

void PHOB(int MODE,double PBOOS1[4],double vec[4]){
  double BET1[3],GAM1,PB;
  static int j0=1;
  int J;


  PB=sqrt(PBOOS1[4-j0]*PBOOS1[4-j0]-PBOOS1[3-j0]*PBOOS1[3-j0]-PBOOS1[2-j0]*PBOOS1[2-j0]-PBOOS1[1-j0]*PBOOS1[1-j0]);
  for( J=1; J<4;J++){
    if (MODE==1) BET1[J-j0]=-PBOOS1[J-j0]/PB;
    else BET1[J-j0]= PBOOS1[J-j0]/PB;
  }

  GAM1=PBOOS1[4-j0]/PB;

  //--
  //--   Boost vector 

  PB=BET1[1-j0]*vec[1-j0]+BET1[2-j0]*vec[2-j0]+BET1[3-j0]*vec[3-j0];

  for( J=1; J<4;J++) vec[J-j0]=vec[J-j0]+BET1[J-j0]*(vec[4-j0]+PB/(GAM1+1.0));
  vec[4-j0]=GAM1*vec[4-j0]+PB;
  //--
}


//     *******************************
// Boost along arbitrary axis (as implemented by Ronald Kleiss).
// The method is described in book of Bjorken and Drell
// p boosted into r  from actual frame to rest frame of q
// forth (mode = 1) or back (mode = -1).
// q must be a timelike, p may be arbitrary.
void bostdq(int mode,double qq[4],double pp[4],double r[4]){
  double q[4],p[4],amq,fac;
  static int i=1;
  int k;

  for(k=1;k<=4;k++){
    p[k-i]=pp[k-i];
    q[k-i]=qq[k-i];
  }
  amq =sqrt(q[4-i]*q[4-i]-q[1-i]*q[1-i]-q[2-i]*q[2-i]-q[3-i]*q[3-i]);

  if    (mode == -1){
    r[4-i] = (p[1-i]*q[1-i]+p[2-i]*q[2-i]+p[3-i]*q[3-i]+p[4-i]*q[4-i])/amq;
    fac  = (r[4-i]+p[4-i])/(q[4-i]+amq);
  }
  else if(mode ==  1){
    r[4-i] =(-p[1-i]*q[1-i]-p[2-i]*q[2-i]-p[3-i]*q[3-i]+p[4-i]*q[4-i])/amq;
    fac  =-(r[4-i]+p[4-i])/(q[4-i]+amq);
  }
  else{
    cout << " ++++++++ wrong mode in boostdq " << endl;
    exit(-1);
  }
  r[1-i]=p[1-i]+fac*q[1-i];
  r[2-i]=p[2-i]+fac*q[2-i];
  r[3-i]=p[3-i]+fac*q[3-i];
}


//----------------------------------------------------------------------
//
//    PHOTOS:   PHOton radiation in decays BOost routine '3'
//
//    Purpose:  Boost  vector PVEC  along z-axis where ANGLE = EXP(ETA),
//              ETA is the hyperbolic velocity.
//
//    Input Parameters:  ANGLE, PVEC
//
//    Output Parameter:  PVEC
//
//    Author(s):  S. Jadach                       Created at:  01/01/89
//                B. van Eijk                     Last Update: 12/06/13
//
//----------------------------------------------------------------------
void PHOBO3(double ANGLE,double PVEC[4]){
  int j=1;  // convention of indices of Riemann space must be preserved.
  double QPL,QMI;
  QPL=(PVEC[4-j]+PVEC[3-j])*ANGLE;
  QMI=(PVEC[4-j]-PVEC[3-j])/ANGLE;
  PVEC[3-j]=(QPL-QMI)/2.0;
  PVEC[4-j]=(QPL+QMI)/2.0;
}

} // namespace PhotosUtilities
	
} // namespace Photospp

