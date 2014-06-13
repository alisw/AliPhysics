#include "forZ-MEc.h"
#include "Photos.h"
#include "PhotosUtilities.h"
#include "PH_HEPEVT_Interface.h"
#include "f_Init.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using namespace Photospp;
using namespace PhotosUtilities;

namespace Photospp
{

// from photosC.cxx

extern void PHODMP();
extern double PHINT(int idumm);
// ----------------------------------------------------------------------
// PROVIDES ELECTRIC CHARGE AND WEAK IZOSPIN OF A FAMILY FERMION
// IDFERM=1,2,3,4 DENOTES NEUTRINO, LEPTON, UP AND DOWN QUARK
// NEGATIVE IDFERM=-1,-2,-3,-4, DENOTES ANTIPARTICLE
// IHELIC=+1,-1 DENOTES RIGHT AND LEFT HANDEDNES ( CHIRALITY)
// SIZO3 IS THIRD PROJECTION OF WEAK IZOSPIN (PLUS MINUS HALF)
// AND CHARGE IS ELECTRIC CHARGE IN UNITS OF ELECTRON CHARGE
// KOLOR IS A QCD COLOUR, 1 FOR LEPTON, 3 FOR QUARKS
//
//     called by : EVENTE, EVENTM, FUNTIH, .....
// ----------------------------------------------------------------------

void PhotosMEforZ::GIVIZO(int IDFERM,int IHELIC,double *SIZO3,double *CHARGE,int *KOLOR) {
  //
  int IH, IDTYPE, IC, LEPQUA, IUPDOW; 
  if (IDFERM==0 || abs(IDFERM)>4 || abs(IHELIC)!=1){
    cout << "STOP IN GIVIZO: WRONG PARAMS" << endl;
    exit(-1);
   }

  IH  =IHELIC;
  IDTYPE =abs(IDFERM);
  IC  =IDFERM/IDTYPE;
  LEPQUA=(int)(IDTYPE*0.4999999);
  IUPDOW=IDTYPE-2*LEPQUA-1;
  *CHARGE  =(-IUPDOW+2.0/3.0*LEPQUA)*IC;
  *SIZO3   =0.25*(IC-IH)*(1-2*IUPDOW);
  *KOLOR=1+2*LEPQUA;
  //** NOTE THAT CONVENTIONALY Z0 COUPLING IS
  //** XOUPZ=(SIZO3-CHARGE*SWSQ)/SQRT(SWSQ*(1-SWSQ))
  return;
}


////////////////////////////////////////////////////////////////////////////
///                                                                       //
/// This routine provides unsophisticated Born differential cross section //
/// at the crude x-section level, with Z and gamma s-chanel exchange.     //
///////////////////////////////////////////////////////////////////////////
double PhotosMEforZ::PHBORNM(double svar,double costhe,double T3e,double qe,double T3f,double qf,int NCf){

  double   s,Sw2,MZ,MZ2,GammZ,AlfInv,GFermi;  // t,MW,MW2,
  double   Ve,Ae,thresh;          //  sum,deno,
  double   xe,yf,xf,ye,ff0,ff1,amx2,amfin,Vf,Af;
  double   ReChiZ,SqChiZ,RaZ;     //,RaW,ReChiW,SqChiW;
  double   Born;                  //, BornS;
  //  int  KeyZet,HadMin,KFbeam;
  //  int  i,ke,KFfin,kf,IsGenerated,iKF;
  int  KeyWidFix;
 
  AlfInv= 137.0359895;
  GFermi=1.16639e-5;

  //--------------------------------------------------------------------
  s = svar;
  //------------------------------
  //     EW paratemetrs taken from BornV
  MZ=91.187;
  GammZ=2.50072032;
  Sw2=.22276773;
  //------------------------------
  // Z and gamma couplings to beams (electrons)
  // Z and gamma couplings to final fermions
  // Loop over all flavours defined in m_xpar(400+i)


  //------ incoming fermion
  Ve=  2*T3e -4*qe*Sw2;
  Ae=  2*T3e;
  //------ final fermion couplings
  amfin = 0.000511; //  m_xpar(kf+6)
  Vf =  2*T3f -4*qf*Sw2;
  Af =  2*T3f;
  if(fabs(costhe) > 1.0){
    cout << "+++++STOP in PHBORN: costhe>0 =" << costhe << endl;
    exit(-1);
  }
  MZ2  = MZ*MZ;
  RaZ  = (GFermi *MZ2 *AlfInv  )/( sqrt(2.0) *8.0 *PI); //
  RaZ  = 1/(16.0*Sw2*(1.0-Sw2));
  KeyWidFix = 1;       // fixed width
  KeyWidFix = 0;       // variable width
  if( KeyWidFix == 0 ){
    ReChiZ=(s-MZ2)*s/((s-MZ2)*(s-MZ2)+(GammZ*s/MZ)*(GammZ*s/MZ)) *RaZ;     // variable width
    SqChiZ=      s*s/((s-MZ2)*(s-MZ2)+(GammZ*s/MZ)*(GammZ*s/MZ)) *RaZ*RaZ; // variable width
  }
  else{
      ReChiZ=(s-MZ2)*s/((s-MZ2)*(s-MZ2)+(GammZ*MZ)*(GammZ*MZ)) *RaZ;     // fixed width
      SqChiZ=      s*s/((s-MZ2)*(s-MZ2)+(GammZ*MZ)*(GammZ*MZ)) *RaZ*RaZ; // fixed width
  }
  xe= Ve*Ve +Ae*Ae;
  xf= Vf*Vf +Af*Af;
  ye= 2*Ve*Ae;
  yf= 2*Vf*Af;
  ff0= qe*qe*qf*qf +2*ReChiZ*qe*qf*Ve*Vf +SqChiZ*xe*xf;
  ff1=             +2*ReChiZ*qe*qf*Ae*Af +SqChiZ*ye*yf;
  Born    = (1.0+ costhe*costhe)*ff0 +2.0*costhe*ff1;
  // Colour factor
  Born = NCf*Born;
  // Crude method of correcting threshold, cos(theta) depencence incorrect!!!
  if(    svar <  4.0*amfin*amfin){ 
    thresh=0.0;
  }
  else if(svar < 16.0*amfin*amfin){
    amx2=4.0*amfin*amfin/svar;
    thresh=sqrt(1.0-amx2)*(1.0+amx2/2.0);
  }
  else{
    thresh=1.0;
  }

  Born= Born*thresh;
  return Born;
}


// ----------------------------------------------------------------------
// THIS ROUTINE CALCULATES  BORN ASYMMETRY.
// IT EXPLOITS THE FACT THAT BORN X. SECTION = A + B*C + D*C**2
//
//     called by : EVENTM
// ----------------------------------------------------------------------
//
double PhotosMEforZ::AFBCALC(double SVAR,int IDEE,int IDFF){
  int KOLOR,KOLOR1;
  double T3e,qe,T3f,qf,A,B;
  GIVIZO(IDEE,-1,&T3e,&qe,&KOLOR);
  GIVIZO(IDFF,-1,&T3f,&qf,&KOLOR1);

  A=PHBORNM(SVAR,0.5,T3e,qe,T3f,qf,KOLOR*KOLOR1);
  B=PHBORNM(SVAR,-0.5,T3e,qe,T3f,qf,KOLOR*KOLOR1);
  return (A-B)/(A+B)*5.0/2.0 *3.0/8.0;
}


int PhotosMEforZ::GETIDEE(int IDE){

  int IDEE;
  IDEE=-555;
  if((IDE==11)       || (IDE== 13) || (IDE== 15)){
    IDEE=2;
  }
  else if((IDE==-11) || (IDE==-13) || (IDE==-15)){
    IDEE=-2;
  }
  else if((IDE== 12) || (IDE== 14) || (IDE== 16)){
    IDEE=1;
  }
  else if((IDE==-12) || (IDE==-14) || (IDE==-16)){
    IDEE=-1;
  }
  else if((IDE==  1) || (IDE==  3) || (IDE==  5)){
    IDEE=4;
  }
  else if((IDE== -1) || (IDE== -3) || (IDE== -5)){
    IDEE=-4;
  }
  else if((IDE==  2) || (IDE==  4) || (IDE==  6)){
    IDEE=3;
  }
  else if((IDE==- 2) || (IDE== -4) || (IDE== -6)){
    IDEE=-3;
  }
  if(IDEE==-555) {cout << " ERROR IN GETIDEE of PHOTS Z-ME: I3= &4i"<<IDEE<<endl;}
  return IDEE;
}




//----------------------------------------------------------------------
//
//    PHASYZ:   PHotosASYmmetry of Z
//
//    Purpose:  Calculates born level asymmetry for Z
//              between distributions (1-C)**2 and (1+C)**2
//              At present dummy, requrires effective Z and gamma 
//              Couplings and also spin polarization states
//              For initial and final states.
//              To be correct this function need to be tuned
//              to host generator. Axes orientation polarisation
//              conventions etc etc. 
//              Modularity of PHOTOS would break. 
//
//    Input Parameters:   SVAR
//
//    Output Parameters:  Function value
//
//    Author(s):  Z. Was                          Created at:  10/12/05
//                                                Last Update: 19/06/13
//
//----------------------------------------------------------------------
double PhotosMEforZ::PHASYZ(double SVAR,int IDE, int IDF){

  double AFB;
  int IDEE,IDFF;

  IDEE=abs(GETIDEE(IDE));
  IDFF=abs(GETIDEE(IDF));
  AFB= -AFBCALC(SVAR,IDEE,IDFF);
  //      AFB=0
  return 4.0/3.0*AFB;
  //      write(*,*) 'IDE=',IDE,'  IDF=',IDF,'  SVAR=',SVAR,'AFB=',AFB
}

//----------------------------------------------------------------------
//
//    PHWTNLO:   PHotosWTatNLO
//
//    Purpose:  calculates instead of interference weight
//              complete NLO weight for vector boson decays
//              of pure vector (or pseudovector) couplings
//              Proper orientation of beams required.
//              This is not standard in PHOTOS.
//              At NLO more tuning than in standard is needed.
//               
//              
//
//    Input Parameters:   as in function declaration
//
//    Output Parameters:  Function value
//
//    Author(s):  Z. Was                          Created at:  08/12/05
//                                                Last Update: 20/06/13
//
//----------------------------------------------------------------------
double PhotosMEforZ::Zphwtnlo(double svar,double xk,int IDHEP3,int IREP,double qp[4],double qm[4],double ph[4],double pp[4],double pm[4],double COSTHG,double BETA,double th1,int IDE,int IDF){
  double C,s,xkaM,xkaP,t,u,t1,u1,BT,BU;
  double waga,wagan2;
  static int i=1;
  int IBREM;


  // IBREM is spurious but it numbers branches of MUSTRAAL
  IBREM=1;
  if (IREP==1)  IBREM=-1;

  // we calculate C and S, note that TH1 exists in MUSTRAAL as well. 

  C=cos(th1); // this parameter is calculated outside of the class

  // from off line application we had:
  if(IBREM==-1) C=-C;
  // ... we may want to re-check it. 
  s=sqrt(1.0-C*C);

  if (IBREM==1){
    xkaM=(qp[4-i]*ph[4-i]-qp[3-i]*ph[3-i]-qp[2-i]*ph[2-i]-qp[1-i]*ph[1-i])/xk;
    xkaP=(qm[4-i]*ph[4-i]-qm[3-i]*ph[3-i]-qm[2-i]*ph[2-i]-qm[1-i]*ph[1-i])/xk;
  }
  else{
    xkaP=(qp[4-i]*ph[4-i]-qp[3-i]*ph[3-i]-qp[2-i]*ph[2-i]-qp[1-i]*ph[1-i])/xk;
    xkaM=(qm[4-i]*ph[4-i]-qm[3-i]*ph[3-i]-qm[2-i]*ph[2-i]-qm[1-i]*ph[1-i])/xk;
  }   

  //        XK=2*PHEP(4,nhep)/PHEP(4,1)/xphmax   ! it is not used becuse here
  //                                             ! order of emissions is meaningless
  //
  //        DELTA=2*PHEP(5,4)**2/svar/(1+(1-XK)**2)*(xKAP/xKAM+xKAM/xKAP)
  //        waga=svar/4./xkap
  //        waga=waga*(1.D0-COSTHG*BETA) ! sprawdzone 1= svar/xKAp/4   * (1.D0-COSTHG*BETA)
  //        waga=waga*(1-delta) /wt2 ! sprawdzone ze to jest =2/(1.D0+COSTHG*BETA)
  //                                 ! czyli ubija de-interferencje
 

  // this is true only for intermediate resonances with afb=0!
  t =2*(qp[4-i]*pp[4-i]-qp[3-i]*pp[3-i]-qp[2-i]*pp[2-i]-qp[1-i]*pp[1-i]);
  u =2*(qm[4-i]*pp[4-i]-qm[3-i]*pp[3-i]-qm[2-i]*pp[2-i]-qm[1-i]*pp[1-i]);
  u1=2*(qp[4-i]*pm[4-i]-qp[3-i]*pm[3-i]-qp[2-i]*pm[2-i]-qp[1-i]*pm[1-i]);
  t1=2*(qm[4-i]*pm[4-i]-qm[3-i]*pm[3-i]-qm[2-i]*pm[2-i]-qm[1-i]*pm[1-i]);

  // basically irrelevant lines  ...
  t =t - (qp[4-i]*qp[4-i]-qp[3-i]*qp[3-i]-qp[2-i]*qp[2-i]-qp[1-i]*qp[1-i]);
  u =u - (qm[4-i]*qm[4-i]-qm[3-i]*qm[3-i]-qm[2-i]*qm[2-i]-qm[1-i]*qm[1-i]);
  u1=u1- (qp[4-i]*qp[4-i]-qp[3-i]*qp[3-i]-qp[2-i]*qp[2-i]-qp[1-i]*qp[1-i]);
  t1=t1- (qm[4-i]*qm[4-i]-qm[3-i]*qm[3-i]-qm[2-i]*qm[2-i]-qm[1-i]*qm[1-i]);

  // we adjust to what is f-st,s-nd beam flavour 
  if (IDE*IDHEP3>0){
    BT=1.0+PHASYZ(svar,IDE,IDF);
    BU=1.0-PHASYZ(svar,IDE,IDF);
  }
  else{
    BT=1.0-PHASYZ(svar,IDE,IDF);
    BU=1.0+PHASYZ(svar,IDE,IDF);
  }  
  wagan2=2*(BT*t*t+BU*u*u+BT*t1*t1+BU*u1*u1)
    /(1+(1-xk)*(1-xk))* 2.0/(BT*(1-C)*(1-C)+BU*(1+C)*(1+C))/svar/svar;

  //!        waga=waga*wagan2
  //!        waga=waga*(1-delta) /wt2 ! sprawdzone ze to jest =2/(1.D0+COSTHG*BETA)
  waga=2/(1.0+COSTHG*BETA)*wagan2;  
  //!     %       * svar/4./xkap*(1.D0-COSTHG*BETA)*sqrt(1.0-xk)

  if(wagan2<=3.8) return waga;

  // 
  // exceptional case  wagan2>3.8
  // it should correspond to extremely high bremssthahlung in multiphot conf.  
  //
  FILE *PHLUN = stdout;


  //         fprintf(PHLUN,"") 'phwtnlo= ',phwtnlo
  //         fprintf(PHLUN,"") 'idhepy= ',IDHEP[1-i],IDHEP[2-i],IDHEP[3-i],IDHEP[4-i],IDHEP(5)
  fprintf(PHLUN," IDE= %i  IDF= %i",IDE,IDF);
  fprintf(PHLUN,"bt,bu,bt+bu= %f %f %f",BT,BU,BT+BU);
  PHODMP();  // we will activate this once PHODMP(); is re-written

  fprintf(PHLUN," "); 
  fprintf(PHLUN,"%i %i <-- IREP,IBREM", IREP,IBREM);
  //!        fprintf(PHLUN,"%f %f %f %f %f") 'pneutr= ',phomom_.pneutr[0],phomom_.pneutr[1],phomom_.pneutr[2],phomom_.pneutr[3],phomom_.pneutr[4];
  fprintf(PHLUN,"%f %f %f %f  qp    = ",qp[0],qp[1],qp[2],qp[3]);
  fprintf(PHLUN,"%f %f %f %f  qm    = ",qm[0],qm[1],qm[2],qm[3]);
  fprintf(PHLUN," ");
  fprintf(PHLUN,"%f %f %f %f  ph    = ",ph[0],ph[1],ph[2],ph[3]);
  //        fprintf(PHLUN,"") 'p1= ',PHEP(1,1),PHEP(2,1),PHEP(3,1),PHEP(4,1)
  //        fprintf(PHLUN,"") 'p2= ',PHEP(1,2),PHEP(2,2),PHEP(3,2),PHEP(4,2)
  //        fprintf(PHLUN,"") 'p3= ',PHEP(1,3),PHEP(2,3),PHEP(3,3),PHEP(4,3)
  //        fprintf(PHLUN,"") 'p4= ',PHEP(1,4),PHEP(2,4),PHEP(3,4),PHEP(4,4)
  //        fprintf(PHLUN,"") 'p5= ',PHEP(1,5),PHEP(2,5),PHEP(3,5),PHEP(4,5)

  fprintf(PHLUN," c= %f theta= %f",C,th1);
  //         fprintf(PHLUN,"")  'photos waga daje ... IBREM=',IBREM,' waga=',waga
  //         fprintf(PHLUN,"") 'xk,COSTHG,c',xk,COSTHG,c
  //         fprintf(PHLUN,"") svar/4./xkap*(1.D0-COSTHG*BETA), 
  //     $   (1-delta) /wt2 *(1.D0+COSTHG*BETA)/2, wagan2
  //         fprintf(PHLUN,"") ' delta, wt2,beta',  delta, wt2,beta
  fprintf(PHLUN,"   -  ");
  fprintf(PHLUN,"t,u       = %f %f",t,u);
  fprintf(PHLUN,"t1,u1     = %f %f",t1,u1);
  fprintf(PHLUN,"sredniaki = %f %f",svar*(1-C)/2,svar*(1+C)/2);
  //	   !         fprintf(PHLUN,"") 'xk= %f c= %f COSTHG=  %f' ,xk,c,COSTHG
  fprintf(PHLUN,"PHASYZ(svar)=',%f,' svar= %f',' waga= %f",PHASYZ(svar,IDE,IDF),svar,waga);
  fprintf(PHLUN,"  -  ");
  fprintf(PHLUN,"BT-part= %f BU-part= %f",
                 2*(BT*t*t+BT*t1*t1)
                   /(1+(1-xk)*(1-xk))* 2.0/(BT*(1-C)*(1-C))/svar/svar,
                 2*(BU*u*u+BU*u1*u1)
	           /(1+(1-xk)*(1-xk))* 2.0/(BU*(1+C)*(1+C))/svar/svar);
  fprintf(PHLUN,"BT-part*BU-part= %f wagan2= %f",
                 2*(BT*t*t+BT*t1*t1)
                   /(1+(1-xk)*(1-xk))* 2.0/(BT*(1-C)*(1-C))/svar/svar
                *2*(BU*u*u+BU*u1*u1)
	           /(1+(1-xk)*(1-xk))* 2.0/(BU*(1+C)*(1+C))/svar/svar,  wagan2);

  fprintf(PHLUN,"wagan2= %f",wagan2);
  fprintf(PHLUN," ###################  ");


  wagan2=3.8; //  ! overwrite 
  waga=2/(1.0+COSTHG*BETA)*wagan2 ; 
	   //     %       * svar/4./xkap*(1.D0-COSTHG*BETA)*sqrt(1.0-xk)

  return waga;

}




//----------------------------------------------------------------------
//
//    PHWTNLO:   PHotosWTatNLO
//
//    Purpose:  calculates instead of interference weight
//              complete NLO weight for vector boson decays
//              of pure vector (or pseudovector) couplings
//              Proper orientation of beams required.
//              Uses Zphwtnlo encapsulating actual matrix element
//              At NLO more tuning on kinematical conf.
//              than in standard is needed.
//              Works with KORALZ and KKM// 
//              Note some commented out commons from MUSTAAL, KORALZ
//
//    Input Parameters:   Common /PHOEVT/ /PHOPS/ /PHOREST/ /PHOPRO/
//
//    Output Parameters:  Function value
//
//    Author(s):  Z. Was                          Created at:  08/12/05
//                                                Last Update: 23/06/13
//
//----------------------------------------------------------------------

double PhotosMEforZ::phwtnlo(){
  // fi3 orientation of photon, fi1,th1 orientation of neutral

      //      COMMON/PHOPHS/XPHMAX,XPHOTO,COSTHG,SINTHG

      //      COMMON /PHOREST/ FI3,fi1,th1
      //      COMMON /PHWT/ BETA,WT1,WT2,WT3
      //      COMMON/PHOPRO/PROBH,CORWT,XF,IREP
  //      COMMON/PHOMOM/MCHSQR,MNESQR,PNEUTR(5)
  //  static double PI=3.141592653589793238462643;
  static int i=1;
  int K,L,IDHEP3,IDUM=0;
  int IDE,IDF;
  double  QP[4],QM[4],PH[4],QQ[4],PP[4],PM[4],QQS[4];
  double XK,ENE,svar;

	//      REAL*8 s,c,svar,xkaM,xkaP,xk,phwtnlo,xdumm,PHINT
	//      REAL*8 ENE,a,t,u,t1,u1,wagan2,waga,PHASYZ,BT,BU,ENEB
	//      INTEGER IBREM,K,L,IREP,IDUM,IDHEP3
	//      integer icont,ide,idf
	//      REAL*8 delta

/////////////////////
//         phlupa(299500);


/////////////////////
//        phlupa(299500);

  XK=2.0*pho.phep[pho.nhep-i][4-i]/pho.phep[1-i][4-i];

//  XK=2.0*pho.phep[pho.nhep-i][4-i]/pho.phep[1-i][4-i]/phophs_.xphmax;  // it is not used becuse here
                                                               //order of emissions is meaningless
  if(pho.nhep<=4) XK=0.0;
  // the mother must be Z or gamma*  !!!!
      
  if (XK>1.0e-10 &&(pho.idhep[1-i]==22 || pho.idhep[1-i]==23)){

    //        write(*,*) 'nhep=',nhep
    //      DO K=1,3 ENDDO
    //      IF (K.EQ.1) IBREM= 1
    //      IF (K.EQ.2) IBREM=-1
    //      ICONT=ICONT+1
    //      IBREM=IBREX        ! that will be input parameter.
    //      IBREM=IBREY        ! that IS now   input parameter.

    // We initialize twice 4-vectors, here and again later after boost 
    // must be the same way. Important is how the reduction procedure will work.
    // It seems at present that the beams must be translated to be back to back.
    // this may be done after initialising, thus on 4-vectors.

    for( K=1;K<5;K++){
      PP[K-i]=pho.phep[1-i][K-i];
      PM[K-i]=pho.phep[2-i][K-i];
      QP[K-i]=pho.phep[3-i][K-i];
      QM[K-i]=pho.phep[4-i][K-i];
      PH[K-i]=pho.phep[pho.nhep-i][K-i];
      QQ[K-i]=0.0;
      QQS[K-i]=QP[K-i]+QM[K-i];
    }


    PP[4-i]=(pho.phep[1-i][4-i]+pho.phep[2-i][4-i])/2.0;
    PM[4-i]=(pho.phep[1-i][4-i]+pho.phep[2-i][4-i])/2.0;
    PP[3-i]= PP[4-i];
    PM[3-i]=-PP[4-i];
        
    for(L=5;L<=pho.nhep-1;L++){
      for( K=1;K<5;K++){      
	QQ [K-i]=QQ [K-i]+ pho.phep[L-i][K-i];
	QQS[K-i]=QQS[K-i]+ pho.phep[L-i][K-i];
      }
    }       

    // go to the restframe of 3        
    PHOB(1,QQS,QP);
    PHOB(1,QQS,QM);
    PHOB(1,QQS,QQ);
    ENE=(QP[4-i]+QM[4-i]+QQ[4-i])/2;

    // preserve direction of emitting particle and wipeout QQ 
    if (phopro_.irep==1){
    double  a=sqrt(ENE*ENE-pho.phep[3-i][5-i]*pho.phep[3-i][5-i])/sqrt(QM[4-i]*QM[4-i]-pho.phep[3-i][5-i]*pho.phep[3-i][5-i]);
      QM[1-i]= QM[1-i]*a;
      QM[2-i]= QM[2-i]*a;
      QM[3-i]= QM[3-i]*a;
      QP[1-i]=-QM[1-i];
      QP[2-i]=-QM[2-i];
      QP[3-i]=-QM[3-i];
    }
    else{
    double  a=sqrt(ENE*ENE-pho.phep[3-i][5-i]*pho.phep[3-i][5-i])/sqrt(QP[4-i]*QP[4-i]-pho.phep[3-i][5-i]*pho.phep[3-i][5-i]);
      QP[1-i]= QP[1-i]*a;
      QP[2-i]= QP[2-i]*a;
      QP[3-i]= QP[3-i]*a;
      QM[1-i]=-QP[1-i];
      QM[2-i]=-QP[2-i];
      QM[3-i]=-QP[3-i];
    }
    QP[4-i]=ENE;
    QM[4-i]=ENE;
    // go back to reaction frame (QQ eliminated) 
    PHOB(-1,QQS,QP);
    PHOB(-1,QQS,QM);
    PHOB(-1,QQS,QQ);

    svar=pho.phep[1-i][4-i]*pho.phep[1-i][4-i];

    IDE=hep.idhep[1-i];
    IDF=hep.idhep[4-i];
    if(abs(hep.idhep[4-i])==abs(hep.idhep[3-i])) IDF=hep.idhep[3-i];

    IDHEP3=pho.idhep[3-i];
    return Zphwtnlo(svar,XK,IDHEP3,phopro_.irep,QP,QM,PH,PP,PM,phophs_.costhg,phwt_.beta,phorest_.th1,IDE,IDF);
  }
  else{
      // in other cases we just use default setups.
    return PHINT(IDUM);
  }
}

} // namespace Photospp

