//This class introduce the weights calculation according with Lednicky's algorithm.
//The detailed description of the algorithm can be found in comments to fortran code:
//fsiw.f, fsiini.f  
#ifndef ALIHBTLLWEIGHTS_H
#define ALIHBTLLWEIGHTS_H

#include <TObject.h>
#include "WLedCOMMONS.h"

class AliHBTPair;
class AliHBTLLWeights: public TObject
 {
   public:
     virtual ~AliHBTLLWeights(){;}
     static AliHBTLLWeights* Instance();
     
     void Init(); //put the initial values in fortran commons fsiini, led_bldata
     Double_t GetWeight(const AliHBTPair* partpair); //get weight calculated by Lednicky's algorithm

     void SetTest(Bool_t rtest = kTRUE){ftest = rtest;} //if ftest=0: 
     //physical values of the following  parameters are put automatically                       
     //            in FSIINI (their values are not required)          
     // ftest=1: any values of the following parameters are allowed,    
     //the following parameters are required:                           

     void SetColoumb(Bool_t col = kTRUE){fColoumbSwitch = col;}//: (ICH in fortran code) Coulomb interaction between the two particles ON (OFF)
     void SetQuantumStatistics(Bool_t qss = kTRUE){fQuantStatSwitch = qss;}//IQS: quantum statistics for the two particles ON (OFF) //if non-identical particles automatically off
     void SetStrongInterSwitch(Bool_t sis = kTRUE){fStrongInterSwitch = sis;}//ISI: strong interaction between the two particles ON (OFF)
     void SetColWithResidNuclSwitch(Bool_t crn = kTRUE){fColWithResidNuclSwitch = crn;}//I3C: Coulomb interaction with residual nucleus ON (OFF)  
     void SetApproxModel(Int_t ap){approximationModel=ap;}//NS in Fortran code, 
     //   NS=1  Square well potential,                                                             
     //   NS=3  not used                                                                           
     //   NS=4  scattered wave approximated by the spherical wave,                                 
     //   NS=2  same as NS=4 but the approx. of equal emission times in PRF                        
     //         not required (t=0 approx. used in all other cases).      

     
     void SetRandomPosition(Bool_t rp = kTRUE){fRandomPosition = rp;} //ON=kTRUE(OFF=kFALSE)
     // ON -- calculation of the Gauss source radii if the generator don't allows the source generation (for example MeVSim)
     //if ON the following parameters are requested:
     void SetR1dw(Double_t R){fRadius=R;}   //spherical source model radii                                                                           		                                                                                                        
     void SetLambdaw(Double_t la){flambda=la;}  //lambda=haoticity                                                   

     
     void SetParticlesTypes(Int_t pid1, Int_t pid2){fPID1 = pid1; fPID2 = pid2;} //set AliRoot particles types   
    
     void SetNucleusCharge(Double_t ch){fNuclCharge=ch;} // not used now  (see comments in fortran code)
     void SetNucleusMass(Double_t mass){fNuclMass=mass;} // (see comments in fortran code)


   protected:
     
     Bool_t ftest; 
     Bool_t fColoumbSwitch; 
     Bool_t fQuantStatSwitch; 
     Bool_t fStrongInterSwitch;//Switches strong interactions TRUE=ON
     Bool_t fColWithResidNuclSwitch;//Switches couloumb interaction 
                                    //with residual nucleus TRUE=ON          
     Double_t fNuclMass; //mass 
     Double_t fNuclCharge; //charge	
     
     Bool_t  fRandomPosition;
     Double_t fRadius;
     Double_t flambda;

     
     Double_t wein;

     Int_t approximationModel; //approximation used to calculate Bethe-Salpeter amplitude
                             //   ==1  Square well potential,
                             //   ==3  not used
                             //   ==4  scattered wave approximated by the spherical wave,
                             //   ==2  same as NS=4 but the approx. of equal emission times in PRF
                             //         not required (t=0 approx. used in all other cases).
                             //  Note: if ==2,4, the B-S amplitude diverges at zero distance r* in
                             //         the two-particle c.m.s.; user can specify a cutoff AA in
                             //         SUBROUTINE FSIINI, for example:
                             //         IF(NS.EQ.2.OR.NS.EQ.4)AA=5.D0 !! in 1/GeV --> AA=1. fm
     
     Int_t fPID1;
     Int_t fPID2;

     static  AliHBTLLWeights *fgLLWeights;// pointer to wrapper of Fortran Lednicky code


     static Int_t GetPairCode(Int_t pid1,Int_t pid2);
     static Int_t GetPairCode(const AliHBTPair* partpair);//calculate automatically internal FSIW 
	  //     C----------------------------------------------------------------------               
	  //     C-   LL       1  2  3  4   5    6   7  8  9 10  11  12  13  14 15 16 17               
	  //     C-   part. 1: n  p  n alfa pi+ pi0 pi+ n  p pi+ pi+ pi+ pi- K+ K+ K+ K-               
	  //     C-   part. 2: n  p  p alfa pi- pi0 pi+ d  d  K-  K+  p   p  K- K+ p  p                
	  //     C   NS=1 y/n: +  +  +  +   +    -   -  -  -  -   -   -   -  -  -  -  -                
	  //     C----------------------------------------------------------------------               
	  //     C-   LL       18  19 20  21  22 23 24 25 26    27     28                              
	  //     C-   part. 1:  d  d   t  t   K0 K0  d p  p      p      n                              
	  //     C-   part. 2:  d alfa t alfa K0 K0b t t alfa lambda lambda                            
	  //     C   NS=1 y/n:  -  -   -  -   -  -   - -  -      +      +                              
	  //     C----------------------------------------------------------------------               
	                                                                           

     Double_t fsigma; //constants for spherical source model 
     Double_t fwcons; //
     
   private:
   AliHBTLLWeights();

   public:
     ClassDef(AliHBTLLWeights,1)
 };

#endif
