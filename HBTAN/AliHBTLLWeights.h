/* $Id$ */

// This class introduces the weight's calculation 
// according to the Lednicky's algorithm.
// The detailed description of the algorithm can be found 
// in comments to fortran code:
// fsiw.f, fsiini.f  
#ifndef ALIHBTLLWEIGHTS_H
#define ALIHBTLLWEIGHTS_H

#include <TObject.h>

class AliHBTPair;
class AliHBTLLWeights: public TObject
{
 public:
  virtual ~AliHBTLLWeights(){;}
  AliHBTLLWeights(const AliHBTLLWeights &source) {
    //Copy ctor needed by the coding conventions but not used
    Fatal("AliHBTLLWeights","copy ctor not implemented");
  }
  AliHBTLLWeights & operator=(const AliHBTLLWeights &source) {
    //Assignment operator needed by the coding conventions but not used
    Fatal("AliHBTLLWeights","assignment operator not implemented");
    return * this;
  }
  static AliHBTLLWeights* Instance();
  
  void Init(); //put the initial values in fortran commons fsiini, led_bldata
  Double_t GetWeight(const AliHBTPair* partpair); //get weight calculated by Lednicky's algorithm

  void SetTest(Bool_t rtest = kTRUE){ftest = rtest;} //if ftest=0: 
  //physical values of the following  parameters are put automatically                       
  //            in FSIINI (their values are not required)          
  // ftest=1: any values of the following parameters are allowed,    
  //the following parameters are required:                           

  void SetColoumb(Bool_t col = kTRUE){
    //(ICH in fortran code) Coulomb interaction between 
    //the two particles ON (OFF)
    fColoumbSwitch = col;
  }
  void SetQuantumStatistics(Bool_t qss = kTRUE){
    //IQS: quantum statistics for the two particles ON (OFF)
    //if non-identical particles automatically off
    fQuantStatSwitch = qss;
  }
  void SetStrongInterSwitch(Bool_t sis = kTRUE){
    //ISI: strong interaction between the two particles ON (OFF)
    fStrongInterSwitch = sis;
  }
  void SetColWithResidNuclSwitch(Bool_t crn = kTRUE){
    //I3C: Coulomb interaction with residual nucleus ON (OFF)
    fColWithResidNuclSwitch = crn;
  }
  void SetApproxModel(Int_t ap){
    //NS in Fortran code,
    fApproximationModel=ap;
  }
  //   NS=1  Square well potential,
  //   NS=3  not used
  //   NS=4  scattered wave approximated by the spherical wave,
  //   NS=2  same as NS=4 but the approx. of equal emission times in PRF
  //         not required (t=0 approx. used in all other cases).    

     
  void SetRandomPosition(Bool_t rp = kTRUE){
    //ON=kTRUE(OFF=kFALSE)
    fRandomPosition = rp;
  } 
  // ON -- calculation of the Gauss source radii if the generator 
  //       don't allows the source generation (for example MeVSim)
  //if ON the following parameters are requested:
  void SetR1dw(Double_t R){fRadius=R;}   //spherical source model radii
  void SetLambdaw(Double_t la){flambda=la;}  //lambda=haoticity
  void SetParticlesTypes(Int_t pid1, Int_t pid2){
    //set AliRoot particles types   
    fPID1 = pid1; fPID2 = pid2;
  }
    
  void SetNucleusCharge(Double_t ch){
    // not used now  (see comments in fortran code)
    fNuclCharge=ch;
  }
  void SetNucleusMass(Double_t mass){
    // (see comments in fortran code)
    fNuclMass=mass;
  }
  
  
 protected:
  
  Bool_t ftest; // Switch for automatic setting of all parameters
  Bool_t fColoumbSwitch; // Swith for Couloumb interaction in the pair
  Bool_t fQuantStatSwitch; //Switch for quantum statistics
  Bool_t fStrongInterSwitch;//Switches strong interactions TRUE=ON
  Bool_t fColWithResidNuclSwitch;//Switches couloumb interaction 
  //with residual nucleus TRUE=ON          
  Double_t fNuclMass; //mass 
  Double_t fNuclCharge; //charge	
  
  Bool_t  fRandomPosition; // Radius of Gaussian source
  Double_t fRadius; // Raduis of spheric source
  Double_t flambda; // Chaoticity
  
  
  //  Double_t fWein;
  
  Int_t fApproximationModel; //approximation used to calculate 
  //  Bethe-Salpeter amplitude
  //   ==1  Square well potential,
  //   ==3  not used
  //   ==4  scattered wave approximated by the spherical wave,
  //   ==2  same as NS=4 but the approx. of equal emission times in PRF
  //         not required (t=0 approx. used in all other cases).
  //  Note: if ==2,4, the B-S amplitude diverges at zero distance r* in
  //         the two-particle c.m.s.; user can specify a cutoff AA in
  //         SUBROUTINE FSIINI, for example:
  //         IF(NS.EQ.2.OR.NS.EQ.4)AA=5.D0 !! in 1/GeV --> AA=1. fm
  
  Int_t fPID1; // Type of the first particle
  Int_t fPID2; // Type of the second particle
  
  static  AliHBTLLWeights *fgLLWeights;// pointer to wrapper of Fortran Lednicky code
  
  
  static Int_t GetPairCode(Int_t pid1,Int_t pid2);
  static Int_t GetPairCode(const AliHBTPair* partpair);//calculate automatically internal FSIW
  //----------------------------------------------------------------------
  // LL        1  2  3  4   5    6   7  8  9 10  11  12  13  14 15 16 17
  // part. 1:  n  p  n alfa pi+ pi0 pi+ n  p pi+ pi+ pi+ pi- K+ K+ K+ K-
  // part. 2:  n  p  p alfa pi- pi0 pi+ d  d  K-  K+  p   p  K- K+ p  p
  // NS=1 y/n: +  +  +  +   +    -   -  -  -  -   -   -   -  -  -  -  -
  //----------------------------------------------------------------------
  // LL       18  19 20  21  22 23 24 25 26    27     28
  // part. 1:  d  d   t  t   K0 K0  d p  p      p      n
  // part. 2:  d alfa t alfa K0 K0b t t alfa lambda lambda
  // NS=1 y/n:  -  -   -  -   -  -   - -  -      +      +
  //----------------------------------------------------------------------
  
  Double_t fsigma; //constants for spherical source model 
  Double_t fwcons; //weight of final state interaction?
  
 private:
  AliHBTLLWeights();
  
  ClassDef(AliHBTLLWeights,1)
 };

#endif
