#include "forW-MEc.h"
#include "Photos.h"
#include "f_Init.h"
#include "PH_HEPEVT_Interface.h"
#include <cstdlib>
#include<iostream>
using std::cout;
using std::endl;

namespace Photospp
{

// COMMON /Kleiss_Stirling/spV,bet
double PhotosMEforW::spV[4],PhotosMEforW::bet[4];

// COMMON /mc_parameters/pi,sw,cw,alphaI,qb,mb,mf1,mf2,qf1,qf2,vf,af,mcLUN
double PhotosMEforW::pi,PhotosMEforW::sw,PhotosMEforW::cw,PhotosMEforW::alphaI,PhotosMEforW::qb,PhotosMEforW::mb,PhotosMEforW::mf1,PhotosMEforW::mf2,PhotosMEforW::qf1,PhotosMEforW::qf2,PhotosMEforW::vf,PhotosMEforW::af,PhotosMEforW::mcLUN;

//////////////////////////////////////////////////////////////////
//         small s_{+,-}(p1,p2) for massless case:              //
//                 p1^2 = p2^2 = 0                              // 
//                                                              //
//     k0(0) = 1.d0                                             //
//     k0(1) = 1.d0                                             //
//     k0(2) = 0.d0  Kleisse_Stirling k0 points to X-axis       // 
//     k0(3) = 0.d0                                             //
//                                                              //
//////////////////////////////////////////////////////////////////
complex<double> PhotosMEforW::InProd_zero(double p1[4],int l1,double p2[4],int l2){


  double  forSqrt1,forSqrt2,sqrt1,sqrt2;
  complex<double>    Dcmplx;
  static complex<double>    i_= complex<double>(0.0,1.0);
  bool           equal;



  equal = true;  
  for (int i = 0; i < 4; i++){
 
    if (p1[i]!=p2[i])  equal = equal && false ;
  }               


  if ( (l1==l2) || equal ) return complex<double>(0.0,0.0);

 
  else if ( (l1==+1) && (l2==-1) ){

    forSqrt1 = (p1[0]-p1[1])/(p2[0]-p2[1]);
    forSqrt2 = 1.0/forSqrt1;
    sqrt1    = sqrt(forSqrt2);
    sqrt2    = sqrt(forSqrt1);

    return (p1[2]+i_*p1[3])*sqrt1 -
	   (p2[2]+i_*p2[3])*sqrt2 ;
  }
  else if ( (l1==-1) && (l2==+1) ){

    forSqrt1 = (p1[0]-p1[1])/(p2[0]-p2[1]);
    forSqrt2 = 1.0/forSqrt1;
    sqrt1    = sqrt(forSqrt2);
    sqrt2    = sqrt(forSqrt1);

    return (p2[2]-i_*p2[3])*sqrt2 -
           (p1[2]-i_*p1[3])*sqrt1 ;
  }
  else{
                 

    cout << " "<<endl;             
    cout << " ERROR IN InProd_zero:"<<endl;
    cout << "   WRONG VALUES FOR l1,l2: l1,l2 = -1,+1 "<<endl;
    cout << " "  <<endl;           
    exit(-1);
  }
}

double PhotosMEforW::InSqrt(double p[4],double q[4]){
            
  return sqrt( (p[0]-p[1]) / (q[0]-q[1]) );
}
    
//////////////////////////////////////////////////////////////////
//                                                              //
//  Inner product for massive spinors: Ub(p1,m1,l1)*U(p2,m2,l2) //
//                                                              //
//////////////////////////////////////////////////////////////////

complex<double> PhotosMEforW::InProd_mass(double p1[4],double m1,int l1,double p2[4],double m2,int l2){
  double sqrt1,sqrt2,forSqrt1;


  if ((l1==+1)&&(l2==+1)) {               
    forSqrt1    = (p1[0]-p1[1])/(p2[0]-p2[1]);
    sqrt1       = sqrt(forSqrt1);
    sqrt2       = 1.0/sqrt1;
    return complex<double>(m1*sqrt2+m2*sqrt1,0.0);
  }
  else if  ((l1==+1)&&(l2==-1))                              
    return InProd_zero(p1,+1,p2,-1);

  else if ((l1==-1)&&(l2==+1))                         
    return  InProd_zero(p1,-1,p2,+1);               

  else if ((l1==-1)&&(l2==-1)){                             
    forSqrt1    = (p1[0]-p1[1])/(p2[0]-p2[1]);
    sqrt1       = sqrt(forSqrt1);
    sqrt2       = 1.0/sqrt1;
    return complex<double>(m1*sqrt2+m2*sqrt1,0.0);
  }
  else {        
    cout <<" " <<endl;            
    cout <<" ERROR IN InProd_mass.."<<endl;
    cout <<"       WRONG VALUES FOR l1,l2"<<endl;
    cout <<" " <<endl;            
    exit(-1);
  }
}

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  this is small B_{s}(k,p) function when TrMartix is diaagonal!! //
//                                                                 //
/////////////////////////////////////////////////////////////////////
complex<double> PhotosMEforW::BsFactor(int s,double k[4],double p[4],double m){
    double forSqrt1,sqrt1;
    complex<double>  inPr1;

  if ( s==1 ){ 

    inPr1    = InProd_zero(k,+1,p,-1);
    forSqrt1 = (p[0]-p[1])/(k[0]-k[1]);
    sqrt1    = sqrt(2.0*forSqrt1);  
    //BsFactor = 
    return inPr1*sqrt1;
  }

  else if ( s==-1 ){

    inPr1    = InProd_zero(k,-1,p,+1);
    forSqrt1 = (p[0]-p[1])/(k[0]-k[1]);
    sqrt1    = sqrt(2.0*forSqrt1); 
    //BsFactor = 
    return inPr1*sqrt1;
  }
  else{

    cout << " "<<endl;             
    cout << " ERROR IN BsFactor: "<<endl;
    cout << "       WRONG VALUES FOR s : s = -1,+1"<<endl;
    cout << " "  <<endl;           
    exit(-1);
  }
}




//====================================================================== 
//     
// Eikonal factor of decay W->l_1+l_2+\gamma in terms of K&S objects !
// 
//   EikFactor = q1*eps.p1/k.p1 + q2*eps.p2/k.p2 - q3*eps.p3/k.p3
//
//   indices 1,2 are for charged decay products
//   index 3 is for W
//   
//   q - charge
//    
//======================================================================
complex<double> PhotosMEforW::WDecayEikonalKS_1ph(double p3[4],double p1[4],double p2[4],double k[4],int s){

  double scalProd1,scalProd2,scalProd3;
  complex<double> wdecayeikonalks_1ph,BSoft1,BSoft2;  

  scalProd1 = p1[0]*k[0]-p1[1]*k[1]-p1[2]*k[2]-p1[3]*k[3];
  scalProd2 = p2[0]*k[0]-p2[1]*k[1]-p2[2]*k[2]-p2[3]*k[3];
  scalProd3 = p3[0]*k[0]-p3[1]*k[1]-p3[2]*k[2]-p3[3]*k[3];


  BSoft1  = BsFactor(s,k,p1,mf1);
  BSoft2  = BsFactor(s,k,p2,mf2);
 
  //WDecayEikonalKS_1ph =   
   return sqrt(pi/alphaI)*(-(qf1/scalProd1+qb/scalProd3)*BSoft1   
                           +(qf2/scalProd2-qb/scalProd3)*BSoft2);

}

//======================================================================
//
//       Gauge invariant soft factor for decay!!
//       Gmass2 -- photon mass square       
// 
//======================================================================
complex<double>  PhotosMEforW::SoftFactor(int s,double k[4],double p1[4],double m1,double p2[4],double m2,double Gmass2){

  double ScalProd1,ScalProd2;
  complex<double>  BsFactor2,BsFactor1;
           

  ScalProd1 = k[0]*p1[0]-k[1]*p1[1]-k[2]*p1[2]-k[3]*p1[3];
  ScalProd2 = k[0]*p2[0]-k[1]*p2[1]-k[2]*p2[2]-k[3]*p2[3];
          
  BsFactor1 = BsFactor(s,k,p1,m1);
  BsFactor2 = BsFactor(s,k,p2,m2);

  return + BsFactor2/2.0/(ScalProd2-Gmass2)
	 - BsFactor1/2.0/(ScalProd1-Gmass2);
}

//############################################################################# 
//                                                                            #
//                         \ eps(k,0,s)                                       # 
//                         /                                                  #   
//                        _\                                                  # 
//                         /\                                                 #
//                         \                                                  #
//                         /                                                  #
//           ---<----------\-------------<---                                 #
//       Ub(p1,m1,l1)                  U(p2,m2,l2)                            #
//                                                                            #
//                                                                            #
//             definition of arbitrary light-like vector beta!!               #
//                                                                            #
//              bet[0] = 1.d0                                                 #
//              bet[1] = 1.d0                                                 #
//              bet[2] = 0.d0      <==> bet == k0  expression becomes easy!!  #
//              bet[3] = 0.d0                                                 #
//#############################################################################

complex<double> PhotosMEforW::TrMatrix_zero(double p1[4],double m1,int l1,double k[4],int s,double p2[4],double m2,int l2){

  double forSqrt1,forSqrt2;
  //                            double p1_1[4],p2_1[4];
  double sqrt1,sqrt2;        //       ,scalProd1,scalProd2;
  complex<double>   inPr1,inPr2,inPr3;
  bool          equal;

  equal = true;    
  for (int i = 0; i < 4; i++) 
    if (p1[i] != p2[i])  equal = equal&&false;

                    

  if ( (m1==m2)&&(equal) ){
    //..          
    //..             when:  p1=p2=p <=> m1=m2 TrMatrix_zero is diagonal
    //..               
    if ( (l1==+1)&&(l2==+1) ){ 

      inPr1    = InProd_zero(k,+s,p1,-s);
      forSqrt1 = (p1[0]-p1[1])/(k[0]-k[1]); 
      sqrt1    = sqrt(2.0*forSqrt1);

      return sqrt1*inPr1;
    }  
 
    else if ( (l1==+1)&&(l2==-1) ){                

      return complex<double>(0.0,0.0);}
                     

    else if ( (l1==-1)&&(l2==+1) ){               

      return complex<double>(0.0,0.0);
    } 

    else if ( (l1==-1)&&(l2==-1) ){                

      inPr1    = InProd_zero(k,+s,p1,-s);
      forSqrt1 = (p1[0]-p1[1])/(k[0]-k[1]); 
      sqrt1    = sqrt(2.0*forSqrt1);

      return sqrt1*inPr1;
    }  
          
    else{ 
        
      cout << ""  <<endl;           
      cout << " ERROR IN  TrMatrix_zero: " <<endl;
      cout << "       WRONG VALUES FOR l1,l2,s" <<endl; 
      cout <<  "" <<endl;             
      exit(-1);

    }       

  }

  if ( (l1==+1)&&(l2==+1)&&(s==+1) ){

    inPr1    = InProd_zero(k,+1,p1,-1);
    forSqrt1 = (p2[0]-p2[1])/(k[0]-k[1]);
    sqrt1    = sqrt(2.0*forSqrt1);                   
 
    return sqrt1*inPr1;
  }
  else if ( (l1==+1)&&(l2==-1)&&(s==+1) ) {
 
    return complex<double>(0.0,0.0);
  }

  else if( (l1==-1)&&(l2==+1)&&(s==+1) ){
  
    forSqrt1 = (p1[0]-p1[1])/(p2[0]-p2[1]);             
    forSqrt2 = 1.0/forSqrt1;
    sqrt1    = sqrt(2.0*forSqrt1);                   
    sqrt2    = sqrt(2.0*forSqrt2);                   
                     
    return complex<double>(m2*sqrt1-m1*sqrt2,0.0);
  }
  else if ( (l1==-1)&&(l2==-1)&&(s==+1) ){ 

    inPr1    = InProd_zero(k,+1,p2,-1);
    forSqrt1 = (p1[0]-p1[1])/(k[0]-k[1]);
    sqrt1    = sqrt(2.0*forSqrt1);                   
  
    return inPr1*sqrt1;
  }

  else if ( (l1==+1)&&(l2==+1)&&(s==-1) ){
 
    inPr1    = -InProd_zero(k,-1,p2,+1);
    forSqrt1 = (p1[0]-p1[1])/(k[0]-k[1]);
    sqrt1    = sqrt(2.0*forSqrt1);                   
 
    return   -sqrt1*inPr1;
  }

  else if ( (l1==+1)&&(l2==-1)&&(s==-1) ){ 
           
    forSqrt1 = (p1[0]-p1[1])/(p2[0]-p2[1]);     
    forSqrt2 = 1.0/forSqrt1;
    sqrt1    = sqrt(2.0*forSqrt1);                   
    sqrt2    = sqrt(2.0*forSqrt2);                   
                     
    return complex<double>(m2*sqrt1-m1*sqrt2,0.0);
  }

  else if ( (l1==-1)&&(l2==+1)&&(s==-1) ){ 

    return complex<double>(0.0,0.0);
  }

  else if( (l1==-1)&&(l2==-1)&&(s==-1) ){ 

    inPr1    = -InProd_zero(k,-1,p1,+1);
    forSqrt1 = (p2[0]-p2[1])/(k[0]-k[1]);
    sqrt1    = sqrt(2.0*forSqrt1);                   
  
    return -inPr1*sqrt1;
  }
  else {     

    cout << "" << endl;
    cout << " ERROR IN TrMatrix_zero: " << endl;
    cout << "    WRONG VALUES FOR l1,l2,s" << endl;
    cout << "" << endl;             
    exit(-1);
  }

}



////////////////////////////////////////////////////////////////
//          transition matrix for massive boson               //
//                                                            // 
//                                                            //
//                         \ eps(k,m,s)                       //
//                         /                                  // 
//                        _\                                  //
//                         /\ k                               // 
//                         \                                  //
//             <-- p1      /         <-- p2                   //                       
//           ---<----------\----------<---                    //
//       Ub(p1,m1,l1)                  U(p2,m2,l2)            //
//                                                            // 
////////////////////////////////////////////////////////////////                         
complex<double> PhotosMEforW::TrMatrix_mass(double p1[4],double m1,int l1,double k[4],double m,int s,double p2[4],double m2,int l2){


  double forSqrt1,forSqrt2;
  double k_1[4],k_2[4];
  double forSqrt3,forSqrt4,sqrt3,sqrt1,sqrt2,sqrt4;
  complex<double>   inPr1,inPr2,inPr3,inPr4;

  for (int i = 0; i < 4; i++) {
    k_1[i] = 1.0/2.0*(k[i] - m*spV[i]);
    k_2[i] = 1.0/2.0*(k[i] + m*spV[i]);                                
  }

  if ( (l1==+1)&&(l2==+1)&&(s==0) ){ 
                
    inPr1 = InProd_zero(p1,+1,k_2,-1);
    inPr2 = InProd_zero(p2,-1,k_2,+1);
    inPr3 = InProd_zero(p1,+1,k_1,-1);
    inPr4 = InProd_zero(p2,-1,k_1,+1);
    sqrt1 = sqrt(p1[0]-p1[1]);
    sqrt2 = sqrt(p2[0]-p2[1]);
    sqrt3 = m1*m2/sqrt1/sqrt2;

              return                 
                            (inPr1*inPr2-inPr3*inPr4)*(vf+af)/m 
		+ (k_1[0]-k_2[0]-k_1[1]+k_2[1])*sqrt3*(vf-af)/m; 
  }       
                 
  else if ( (l1==+1)&&(l2==-1)&&(s==0) ){

    inPr1 = InProd_zero(p1,+1,k_1,-1);
    inPr2 = InProd_zero(p1,+1,k_2,-1);
    inPr3 = InProd_zero(p2,+1,k_2,-1);
    inPr4 = InProd_zero(p2,+1,k_1,-1);

    forSqrt1 = (k_1[0]-k_1[1])/(p2[0]-p2[1]);
    forSqrt2 = (k_2[0]-k_2[1])/(p2[0]-p2[1]);
    forSqrt3 = (k_2[0]-k_2[1])/(p1[0]-p1[1]);
    forSqrt4 = (k_1[0]-k_1[1])/(p1[0]-p1[1]);
    sqrt1 = sqrt(forSqrt1);
    sqrt2 = sqrt(forSqrt2);
    sqrt3 = sqrt(forSqrt3);
    sqrt4 = sqrt(forSqrt4);     

              return 
                  (inPr1*sqrt1 - inPr2*sqrt2)*(vf+af)*m2/m
		+ (inPr3*sqrt3 - inPr4*sqrt4)*(vf-af)*m1/m;
  }
  else if ( (l1==-1)&&(l2==+1)&&(s==0) ){ 

    inPr1 = InProd_zero(p1,-1,k_1,+1);
    inPr2 = InProd_zero(p1,-1,k_2,+1);
    inPr3 = InProd_zero(p2,-1,k_2,+1);
    inPr4 = InProd_zero(p2,-1,k_1,+1);

    forSqrt1 = (k_1[0]-k_1[1])/(p2[0]-p2[1]);
    forSqrt2 = (k_2[0]-k_2[1])/(p2[0]-p2[1]);
    forSqrt3 = (k_2[0]-k_2[1])/(p1[0]-p1[1]);
    forSqrt4 = (k_1[0]-k_1[1])/(p1[0]-p1[1]);
    sqrt1 = sqrt(forSqrt1);
    sqrt2 = sqrt(forSqrt2);
    sqrt3 = sqrt(forSqrt3);
    sqrt4 = sqrt(forSqrt4);     
        
              return 
                  (inPr1*sqrt1 - inPr2*sqrt2)*(vf-af)*m2/m
		+ (inPr3*sqrt3 - inPr4*sqrt4)*(vf+af)*m1/m;
  }
  else if  ( (l1==-1)&&(l2==-1)&&(s==0) ){ 

    inPr1 = InProd_zero(p2,+1,k_2,-1);
    inPr2 = InProd_zero(p1,-1,k_2,+1);
    inPr3 = InProd_zero(p2,+1,k_1,-1);
    inPr4 = InProd_zero(p1,-1,k_1,+1);
    sqrt1 = sqrt(p1[0]-p1[1]);
    sqrt2 = sqrt(p2[0]-p2[1]);
    sqrt3 = m1*m2/sqrt1/sqrt2;

             return                    
                         (inPr1*inPr2 - inPr3*inPr4)*(vf-af)/m  
	       + (k_1[0]-k_2[0]-k_1[1]+k_2[1])*sqrt3*(vf+af)/m;
  }
  else if ( (l1==+1)&&(l2==+1)&&(s==+1) ){ 

    inPr1 = InProd_zero(p1,+1,k_1,-1);
    inPr2 = InProd_zero(k_2,-1,p2,+1);
    inPr3 = inPr1*inPr2;

    forSqrt1 = (k_1[0]-k_1[1])/(p1[0]-p1[1]);                       
    forSqrt2 = (k_2[0]-k_2[1])/(p2[0]-p2[1]);  
    sqrt1 = sqrt(forSqrt1);                   
    sqrt2 = sqrt(forSqrt2);                   
    sqrt3 = m1*m2*sqrt1*sqrt2;

             return
	       sqrt(2.0)/m*(inPr3*(vf+af)+sqrt3*(vf-af));
  }

  else if ( (l1==+1)&&(l2==-1)&&(s==+1) ){

    inPr1 = InProd_zero(p1,+1,k_1,-1);
    inPr2 = InProd_zero(p2,+1,k_1,-1); 

    forSqrt1 = (k_2[0]-k_2[1])/(p2[0]-p2[1]);                      
    forSqrt2 = (k_2[0]-k_2[1])/(p1[0]-p1[1]);                       
    sqrt1 = m2*sqrt(forSqrt1);                   
    sqrt2 = m1*sqrt(forSqrt2);                                     
                     
              return
                      sqrt(2.0)/m*( + inPr1*sqrt1*(vf+af)
                                    - inPr2*sqrt2*(vf-af)
				  );
  }
  else if  ( (l1==-1)&&(l2==+1)&&(s==+1) ){

    inPr1 = InProd_zero(k_2,-1,p2,+1);
    inPr2 = InProd_zero(k_2,-1,p1,+1);

    forSqrt1 = (k_1[0]-k_1[1])/(p1[0]-p1[1]);                       
    forSqrt2 = (k_1[0]-k_1[1])/(p2[0]-p2[1]);                       
    sqrt1 = m1*sqrt(forSqrt1);                   
    sqrt2 = m2*sqrt(forSqrt2);                                     
                     
              return
                      sqrt(2.0)/m*( + inPr1*sqrt1*(vf+af)
                                    - inPr2*sqrt2*(vf-af)
				  );
  }
  else if ( (l1==-1)&&(l2==-1)&&(s==+1) ){ 

    inPr1 = InProd_zero(p2,+1,k_1,-1);
    inPr2 = InProd_zero(k_2,-1,p1,+1);
    inPr3 = inPr1*inPr2;

    forSqrt1 = (k_1[0]-k_1[1])/(p1[0]-p1[1]);                       
    forSqrt2 = (k_2[0]-k_2[1])/(p2[0]-p2[1]);  
    sqrt1 = sqrt(forSqrt1);                  
    sqrt2 = sqrt(forSqrt2);                   
    sqrt3 = m1*m2*sqrt1*sqrt2;

              return 
		sqrt(2.0)/m*(inPr3*(vf-af)+sqrt3*(vf+af));
  }

  else if ( (l1==+1)&&(l2==+1)&&(s==-1) ){ 

    inPr1 = InProd_zero(p2,-1,k_1,+1);
    inPr2 = InProd_zero(k_2,+1,p1,-1);
    inPr3 = inPr1*inPr2;

    forSqrt1 = (k_1[0]-k_1[1])/(p1[0]-p1[1]);                       
    forSqrt2 = (k_2[0]-k_2[1])/(p2[0]-p2[1]);  
    sqrt1 = sqrt(forSqrt1);                   
    sqrt2 = sqrt(forSqrt2);                   
    sqrt3 = m1*m2*sqrt1*sqrt2;

             return               
	       sqrt(2.0)/m*(inPr3*(vf+af)+sqrt3*(vf-af));
  }
  else if ( (l1==+1)&&(l2==-1)&&(s==-1) ){ 

    inPr1 = InProd_zero(k_2,+1,p2,-1);
    inPr2 = InProd_zero(k_2,+1,p1,-1);

    forSqrt1 = (k_1[0]-k_1[1])/(p1[0]-p1[1]);                       
    forSqrt2 = (k_1[0]-k_1[1])/(p2[0]-p2[1]);                       
    sqrt1 = m1*sqrt(forSqrt1);                   
    sqrt2 = m2*sqrt(forSqrt2);                                     
                     
              return
                      sqrt(2.0)/m*(+ inPr1*sqrt1*(vf-af)
                                   - inPr2*sqrt2*(vf+af)
				  );
  }
  else if ( (l1==-1)&&(l2==+1)&&(s==-1) ){

    inPr1 = InProd_zero(p1,-1,k_1,+1);
    inPr2 = InProd_zero(p2,-1,k_1,+1);

    forSqrt1 = (k_2[0]-k_2[1])/(p2[0]-p2[1]);                       
    forSqrt2 = (k_2[0]-k_2[1])/(p1[0]-p1[1]);                       
    sqrt1 = m2*sqrt(forSqrt1);                   
    sqrt2 = m1*sqrt(forSqrt2);                                     
                     
              return
                      sqrt(2.0)/m*(+ inPr1*sqrt1*(vf-af)
                                   - inPr2*sqrt2*(vf+af) 
				  );
  }
  else if ( (l1==-1)&&(l2==-1)&&(s==-1) ){ 

    inPr1 = InProd_zero(p1,-1,k_1,+1);
    inPr2 = InProd_zero(k_2,+1,p2,-1);
    inPr3 = inPr1*inPr2;

    forSqrt1 = (k_1[0]-k_1[1])/(p1[0]-p1[1]);                       
    forSqrt2 = (k_2[0]-k_2[1])/(p2[0]-p2[1]);  
    sqrt1 = sqrt(forSqrt1);                   
    sqrt2 = sqrt(forSqrt2);                   
    sqrt3 = m1*m2*sqrt1*sqrt2;

             return 
	       sqrt(2.0)/m*(inPr3*(vf-af)+sqrt3*(vf+af));
  }

  else{ 

    cout << " "<< endl;             
    cout << " TrMatrix_mass: Wrong values for l1,l2,s:"<< endl;
    cout << "          l1,l2 = -1,+1; s = -1,0,1 "<< endl;
    cout << " "<< endl;             
    exit(-1);

  } 
         
}


              
//======================================================================                 
//                                                                     =
//                            p1,mf1,l1                                =
//                           /                                         =
//                         \/_                                         = 
//                         /                                           =
//        p3,mb,l3        /                                            =
//              \/\/\/\/\/\      ------> g_(mu,1)*(1+g5_(1))           =
//                         \                                           =
//                         _\/                                         = 
//                           \                                         =
//                            p2,mf2,l2                                =
// INPUT : p1,m1,l1; p2,m2,l2; p3,m3,l3  -- momenta,mass and helicity  =
//                                                                     =
// OUTPUT: value of functions            -- decay amplitude            =
//                                                                     =
//======================================================================
complex<double> PhotosMEforW::WDecayBornAmpKS_1ph(double p3[4],int l3,double p1[4],int l1,double p2[4],int l2){

  double coeff;
 

  coeff = sqrt(pi/alphaI/2.0)/sw;      // vertex: g/2/sqrt(2)

  return coeff*TrMatrix_mass(p2,mf2,l2,p3,mb,l3,p1,-mf1,-l1);
}        


//======================================================================                 
//                k,0,l                                                =
//                   \        p1,mf1,l1                                =
//                   /       /                                         =
//                   \     \/_                                         = 
//                   /     /                                           =
//        p3,mb,l3   \    /                                            =
//              \/\/\/\/\/\      ------> g_(mu,1)*(1+g5_(1))           =
//                         \                                           =
//                         _\/                                         = 
//                           \                                         =
//                            p2,mf2,l2                                =
//           { + }                                                     =
//                            p1,mf1,l1                                =
//                           /                                         = 
//                         \/_~~~~~~~ k,0,s                            = 
//                         /                                           =
//        p3,mb,l3        /                                            =
//              \/\/\/\/\/\      ------> g_(mu,1)*(1+g5_(1))           =
//                         \                                           =
//                         _\/                                         = 
//                           \                                         =
//                            p2,mf2,l2                                =
//           { + }                                                     =
//                            p1,mf1,l1                                =
//                           /                                         =
//                         \/_                                         =
//                         /                                           =
//        p3,mb,l3        /                                            =
//              \/\/\/\/\/\      ------> g_(mu,1)*(1+g5_(1))           =
//                         \                                           =
//                         _\/ ~~~~~~~ k,0,s                           =
//                           \                                         =
//                             p2,mf2,l2                               =
//                                                                     =
//                   all momentas, exept k are incoming !!!            =
//                                                                     =
// This function culculates The W-ff\gamma decay amplitude into permion=
// pair and one photon using Kleisse&Stirling method for helicity      =
// amplitudes, which includes three above feynman diagramms..          = 
//                                                                     =
// INPUT : p1,m1,l1; p2,m2,l2; p3,m3,l3  -- momenta,mass and helicity  =
//                                                                     =
// OUTPUT: value of functions            -- decay amplitude            =
//                                                                     =
//======================================================================

complex<double> PhotosMEforW::WDecayAmplitudeKS_1ph(double p3[4],int l3,double p1[4],int l1,double p2[4],int l2,double k[4],int s){
 
  double scalProd1,scalProd2,scalProd3,coeff;  //,theta3,ph3;
  complex<double>  bornAmp,TrMx1,TrMx2;
  complex<double>  BSoft1,BSoft2;  

  coeff = sqrt(2.0)*pi/sw/alphaI;      // vertex: g/2/sqrt[2] * e

  scalProd1 = p1[0]*k[0]-p1[1]*k[1]-p1[2]*k[2]-p1[3]*k[3];
  scalProd2 = p2[0]*k[0]-p2[1]*k[1]-p2[2]*k[2]-p2[3]*k[3];
  scalProd3 = p3[0]*k[0]-p3[1]*k[1]-p3[2]*k[2]-p3[3]*k[3];

  BSoft1  = BsFactor(s,k,p1,mf1);
  BSoft2  = BsFactor(s,k,p2,mf2);
  bornAmp = TrMatrix_mass(p2,mf2,l2,p3,mb,l3,p1,-mf1,-l1);
  TrMx1   = complex<double>(0.0,0.0);  
  TrMx2   = complex<double>(0.0,0.0);  
 
  for (int la1 = -1; la1< 3 ; la1+=2) {
    //      DO la1=-1,1,2            
    TrMx1 = TrMx1 + TrMatrix_zero(k,0.0,-la1,k,s,p1,-mf1,-l1)*
	            TrMatrix_mass(p2,mf2,l2,p3,mb,l3,k,0.0,-la1);
    TrMx2 = TrMx2 + TrMatrix_zero(p2,mf2,l2,k,s,k,0.0,la1)*
	            TrMatrix_mass(k,0.0,la1,p3,mb,l3,p1,-mf1,-l1);
  }

    return coeff * (        
	     + (-(qf1/scalProd1+qb/scalProd3)*BSoft1              // IR-divergent part of amplitude      
                +(qf2/scalProd2-qb/scalProd3)*BSoft2)/2.0*bornAmp
	     //
	     - (qf1/scalProd1+qb/scalProd3)*TrMx1/2.0             // IR-finite part of amplitude            
             + (qf2/scalProd2-qb/scalProd3)*TrMx2/2.0
		    ); 
}



//========================================================
//        The squared eikonal factor for W decay         =
//        into fermion pair and one photon               =
//  INPUT :                                              =
//                                                       =
//  OUTPUT:                                              =
//========================================================       

double PhotosMEforW::WDecayEikonalSqrKS_1ph(double p3[4],double p1[4],double p2[4],double k[4]){
  double spinSumAvrg;
  complex<double>  wDecAmp;

  spinSumAvrg = 0.0;
  for (int s = -1; s< 3 ; s+=2) {
    wDecAmp     = WDecayEikonalKS_1ph(p3,p1,p2,k,s);
    spinSumAvrg = spinSumAvrg + real(wDecAmp*conj(wDecAmp)); 
  }
  return spinSumAvrg;
}
             
//========================================================
//        The squared eikonal factor for W decay         =
//        into fermion pair and one photon               =
//  INPUT :                                              =
//                                                       =
//  OUTPUT:                                              =
//========================================================       

double PhotosMEforW::WDecayBornAmpSqrKS_1ph(double p3[4],double p1[4],double p2[4]){
  double spinSumAvrg;
  complex<double> wDecAmp;

  spinSumAvrg = 0.0;
  for (int l3 = -1; l3< 2 ; l3++) {
    for (int l1 = -1; l1< 3 ; l1+=2) {
      for (int l2 = -1; l2< 3 ; l2+=2) {
	wDecAmp     = WDecayBornAmpKS_1ph(p3,l3,p1,l1,p2,l2);
	spinSumAvrg = spinSumAvrg + real(wDecAmp*conj(wDecAmp)); 
      }
    }
  }
  return spinSumAvrg;
}



//========================================================
//        The squared amplitude for W decay              =
//        into fermion pair and one photon               =
//  INPUT :                                              =
//                                                       =
//  OUTPUT:                                              =
//========================================================       

double PhotosMEforW::WDecayAmplitudeSqrKS_1ph(double p3[4],double p1[4],double p2[4], double k[4]){

  double spinSumAvrg;
  complex<double> wDecAmp;

  spinSumAvrg = 0.0;
  for (int l3 = -1; l3< 2 ; l3++) {
    for (int l1 = -1; l1< 3 ; l1+=2) {
      for (int l2 = -1; l2< 3 ; l2+=2) {
	for (int s  = -1; s < 3 ;  s+=2) {
	  wDecAmp     = WDecayAmplitudeKS_1ph(p3,l3,p1,l1,p2,l2,k,s);
	  spinSumAvrg = spinSumAvrg + real(wDecAmp*conj(wDecAmp)); 
	}
      }
    }
  }
  return spinSumAvrg;
      


//$$$
//$$$
//$$$
//$$$
//$$$
//$$$
//$$$
//$$$
//$$$
//$$$
//$$$
//$$$
//$$$
//$$$
//$$$
//$$$   WffGammaME.f  ends above: 
//$$$
//$$$
//$$$
//$$$
}



//C========================================================== ==
//C========================================================== ==
//C these will be public for PHOTOS functions of W_ME class   ==
//C========================================================== ==
//C========================================================== ==

double PhotosMEforW::SANC_WT(double PW[4],double PNE[4],double PMU[4],double PPHOT[4],double B_PW[4],double B_PNE[4],double B_PMU[4]){


  //..        Exact amplitude square      
  double AMPSQR=WDecayAmplitudeSqrKS_1ph(PW,PNE,PMU,PPHOT);

  double EIKONALFACTOR=WDecayBornAmpSqrKS_1ph(B_PW,B_PNE,B_PMU)
                       *WDecayEikonalSqrKS_1ph(PW,PNE,PMU,PPHOT);
      
    //..        New weight

    //           cout << 'B_pne=',B_PNE  << endl;
    //           cout << 'B_PMU=',B_PMU  << endl;
    //           cout << 'bornie=',WDecayBornAmpSqrKS_1ph(B_PW,B_PNE,B_PMU)  << endl;

    //           cout << ' '  << endl;
    //           cout << '  pne=',pne  << endl;
    //           cout << '  pmu=',pmu  << endl;
    //           cout << 'pphot=',pphot  << endl;
    //           cout << ' '  << endl;
    //           cout << '  b_pw=',B_PW  << endl;
    //           cout << '  b_pne=',B_PNE  << endl;
    //           cout << 'b_pmu=',B_PMU  << endl;
 
    //          cout << 'cori=',AMPSQR/EIKONALFACTOR,AMPSQR,EIKONALFACTOR  << endl;
 
  return AMPSQR/EIKONALFACTOR;
    //           
    //          return (1-8*EMU*XPH*(1-COSTHG*BETA)*     
    //                 (MCHREN+2*XPH*SQRT(MPASQR))/
    //                 MPASQR**2/(1-MCHREN/MPASQR)/(4-MCHREN/MPASQR)) 
}


void PhotosMEforW::SANC_INIT1(double QB0,double QF20,double MF10,double MF20,double MB0){
  qb =QB0;
  qf2=QF20;
  mf1=MF10;
  mf2=MF20;
  mb =MB0;
}

void PhotosMEforW::SANC_INIT(double ALPHA,int PHLUN){


  static int SANC_MC_INIT=-123456789;

  //...       Initialization of the W->l\nu\gamma 
  //...       decay Matrix Element parameters 
  if (SANC_MC_INIT==-123456789){
    SANC_MC_INIT=1;

    pi=4*atan(1.0);
    qf1=0.0;                           // neutrino charge
    mf1=1.0e-10;                       // newutrino mass
    vf=1.0;                            // V&A couplings
    af=1.0;
    alphaI=1.0/ALPHA;
    cw=0.881731727;                    // Weak Weinberg angle
    sw=0.471751166;
           

    //...          An auxilary K&S vectors
    bet[0]= 1.0;
    bet[1]= 0.0722794881816159;
    bet[2]=-0.994200045099866;
    bet[3]= 0.0796363353729248; 

    spV[0]= 0.0; 
    spV[1]= 7.22794881816159e-2;
    spV[2]=-0.994200045099866;     
    spV[3]= 7.96363353729248e-2;

    mcLUN = PHLUN;
  } 
}
//----------------------------------------------------------------------
//
//    PHOTOS:   PHOtos Boson W correction weight
//
//    Purpose:  calculates correction weight due to amplitudes of 
//              emission from W boson. It is ecact, but not verified
//              for exponentiation yet.
//              
//              
//              
//
//    Input Parameters:  Common /PHOEVT/, with photon added.
//                       wt  to be corrected
//                       
//                       
//                       
//    Output Parameters: wt
//
//    Author(s):  G. Nanava, Z. Was               Created at:  13/03/03
//                                                Last Update: 22/06/13
//
//----------------------------------------------------------------------
void PhotosMEforW::PHOBWnlo(double *WT){
  //  FILE *PHLUN = stdout;  // printouts from matrix element calculations
  //                            directed with phlun still
  int phlun=6;
  double EMU,MCHREN,BETA,COSTHG,MPASQR,XPH;
  double  PW[4],PMU[4],PPHOT[4],PNE[4];
  double  B_PW[4],B_PNE[4],B_PMU[4]; //,AMPSQR;
  static int i=1;
  int I,IJ,I3,I4,JJ;
  double MB,MF1,MF2,QB,QF2;
  //  double  pi,sw,cw,alphaI,qb,mb,mf1,mf2,qf1,qf2,vf,af;


  //!      write(*,*) 'IDPHOs=',IDPHO(1),IDPHO(2),IDPHO(3),IDPHO(4),IDPHO(5)
  //!      write(*,*) 'IDPHOs=',pho.jdahep[1-i][1-i],npho
  //!      write(*,*) 'hep.IDPHOs=',hep.IDhep(1),hep.IDhep(2),hep.IDhep(3),hep.IDhep(4),hep.IDhep(5)

  //--
        if(abs(pho.idhep[1-i])==24&&
           abs(pho.idhep[pho.jdahep[1-i][1-i]-i  ])>=11&&
           abs(pho.idhep[pho.jdahep[1-i][1-i]-i  ])<=16&&
           abs(pho.idhep[pho.jdahep[1-i][1-i]-i+1])>=11&&
           abs(pho.idhep[pho.jdahep[1-i][1-i]-i+1])<=16     ){

	  if(
            abs(pho.idhep[pho.jdahep[1-i][1-i]-i  ])==11||
            abs(pho.idhep[pho.jdahep[1-i][1-i]-i  ])==13||
            abs(pho.idhep[pho.jdahep[1-i][1-i]-i  ])==15    ){
	    I=pho.jdahep[1-i][1-i];
	  }
	  else{
	    I=pho.jdahep[1-i][1-i]+1;
	  }
	  //..        muon energy   
	  EMU=pho.phep[I-i][4-i];
	  //..        muon mass square
	  MCHREN=fabs(pho.phep[I-i][4-i]*pho.phep[I-i][4-i]-pho.phep[I-i][3-i]*pho.phep[I-i][3-i]
		     -pho.phep[I-i][2-i]*pho.phep[I-i][2-i]-pho.phep[I-i][1-i]*pho.phep[I-i][1-i]);
	  BETA=sqrt(1- MCHREN/ pho.phep[I-i][4-i]*pho.phep[I-i][4-i]);
          COSTHG=((pho.phep[I-i][3-i]*pho.phep[pho.nhep-i][3-i]+pho.phep[I-i][2-i]*pho.phep[pho.nhep-i][2-i]
	          +pho.phep[I-i][1-i]*pho.phep[pho.nhep-i][1-i])/
                  sqrt(pho.phep[I-i][3-i]*pho.phep[I-i][3-i]+pho.phep[I-i][2-i]*pho.phep[I-i][2-i]+pho.phep[I-i][1-i]*pho.phep[I-i][1-i])   /
		  sqrt(pho.phep[pho.nhep-i][3-i]*pho.phep[pho.nhep-i][3-i]+pho.phep[pho.nhep-i][2-i]*pho.phep[pho.nhep-i][2-i]+pho.phep[pho.nhep-i][1-i]*pho.phep[pho.nhep-i][1-i]));
	  MPASQR=pho.phep[1-i][4-i]*pho.phep[1-i][4-i];
	  XPH=pho.phep[pho.nhep-i][4-i];

	  //...       Initialization of the W->l\nu\gamma 
	  //...       decay Matrix Element parameters 
	  SANC_INIT(phocop_.alpha,phlun);


	  MB=pho.phep[1-i][4-i];//                      ! W boson mass
	  MF2=sqrt(MCHREN);//                 ! muon mass
	  I3=-1;
	  for(IJ=1;IJ<=hep.nhep;IJ++){
            if(abs(hep.idhep[IJ-i])==24){ I3=IJ;} //! position of W 
	  }
          if(I3==-1) {cout << " ERROR IN PHOBWnlo of PHOTS W-ME: I3= &2i"<<I3<<endl;}
           if(
              abs(hep.idhep[hep.jdahep[I3-i][1-i]-i  ])==11||
              abs(hep.idhep[hep.jdahep[I3-i][1-i]-i  ])==13||
              abs(hep.idhep[hep.jdahep[I3-i][1-i]-i  ])==15    ){ 
	     I4=hep.jdahep[I3-i][1-i];} //              ! position of lepton
           else{
	     I4=hep.jdahep[I3-i][1-i]+1 ;  //         ! position of lepton
	   }

	   if (hep.idhep[I3-i]==-24) QB=-1.0;//  ! W boson charge
	   if (hep.idhep[I3-i]==+24) QB=+1.0;//   
	   if (hep.idhep[I4-i]>0.0) QF2=-1.0; // ! lepton charge
	   if (hep.idhep[I4-i]<0.0) QF2=+1.0;


	   //...          Particle momenta before foton radiation; effective Born level
	   for( JJ=1; JJ<=4;JJ++){
	     B_PW [(JJ % 4)]=hep.phep[I3-i][JJ-i];//  ! W boson
	     B_PNE[(JJ % 4)]=hep.phep[I3-i][JJ-i]-hep.phep[I4-i][JJ-i];// ! neutrino
	     B_PMU[(JJ % 4)]=hep.phep[I4-i][JJ-i]; // ! muon
	   }

	   //..        Particle monenta after photon radiation
	   for( JJ=1; JJ<=4;JJ++){
             PW   [(JJ % 4)]=pho.phep[1-i][JJ-i];
             PMU  [(JJ % 4)]=pho.phep[I-i][JJ-i];          
             PPHOT[(JJ % 4)]=pho.phep[pho.nhep-i][JJ-i];
             PNE  [(JJ % 4)]=pho.phep[1-i][JJ-i]-pho.phep[I-i][JJ-i]-pho.phep[pho.nhep-i][JJ-i];
           }

	   // two options of calculating neutrino (spectator) mass
           MF1=sqrt(fabs(B_PNE[0]*B_PNE[0]*-B_PNE[1]*B_PNE[1]-B_PNE[2]*B_PNE[2]-B_PNE[3]*B_PNE[3]));
           MF1=sqrt(fabs(  PNE[0]*PNE[0]-  PNE[1]*PNE[1]-  PNE[2]*PNE[2]-  PNE[3]*PNE[3]));

	   SANC_INIT1(QB,QF2,MF1,MF2,MB);
	   *WT=(*WT)*SANC_WT(PW,PNE,PMU,PPHOT,B_PW,B_PNE,B_PMU);
        }
	//      write(*,*)   'AMPSQR/EIKONALFACTOR= ',   AMPSQR/EIKONALFACTOR
}

} // namespace Photospp

