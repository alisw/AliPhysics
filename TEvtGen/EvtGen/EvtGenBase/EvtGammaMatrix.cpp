//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtGammaMatrix.cc
//
// Description: Make gamma matrices availible for the calc. of amplitudes, etc.
//
// Modification history:
//
//    DJL/RYD     September 25, 1996         Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <iostream>
#include <math.h>
#include <assert.h>
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtGammaMatrix.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include <stdlib.h>
using std::endl;
using std::ostream;

EvtGammaMatrix::EvtGammaMatrix(){
  int i,j;

  static EvtComplex zero(0.0,0.0);

  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      _gamma[i][j]=zero;
    }
  }
}

EvtGammaMatrix operator*(const EvtGammaMatrix& g, const EvtComplex& c)
{
    return c*g;
}


EvtGammaMatrix operator*(const EvtComplex& c,const EvtGammaMatrix& g){
  int i,j;

  EvtGammaMatrix temp;

  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      temp._gamma[i][j]=g._gamma[i][j]*c;
    }
  }

  return temp;

}


ostream& operator<<(ostream& s, const EvtGammaMatrix& g){


  s<<"["<<g._gamma[0][0]<<","<<g._gamma[0][1]<<","<<g._gamma[0][2]<<","<<g._gamma[0][3]<<"]"<<endl;
  s<<"["<<g._gamma[1][0]<<","<<g._gamma[1][1]<<","<<g._gamma[1][2]<<","<<g._gamma[1][3]<<"]"<<endl;
  s<<"["<<g._gamma[2][0]<<","<<g._gamma[2][1]<<","<<g._gamma[2][2]<<","<<g._gamma[2][3]<<"]"<<endl;
  s<<"["<<g._gamma[3][0]<<","<<g._gamma[3][1]<<","<<g._gamma[3][2]<<","<<g._gamma[3][3]<<"]"<<endl;

  return s;

}



EvtGammaMatrix::EvtGammaMatrix(const EvtGammaMatrix& gm){
  int i,j;
  
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      _gamma[i][j]=gm._gamma[i][j];
    }
  }
}

EvtGammaMatrix::~EvtGammaMatrix() {}

EvtGammaMatrix& EvtGammaMatrix::operator=(const EvtGammaMatrix& gm){
  int i,j;
  
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      _gamma[i][j]=gm._gamma[i][j];
    }
  }
  return *this;
}

void EvtGammaMatrix::init(){
  int i,j;

  static EvtComplex zero(0.0,0.0);

  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      _gamma[i][j]=zero;
    }
  }
}

const EvtGammaMatrix& EvtGammaMatrix::va0(){

  static EvtGammaMatrix g;
  static int first=1;

  if (first){
    first = 0;
    g._gamma[0][0]=EvtComplex(1.0,0.0);
    g._gamma[0][1]=EvtComplex(0.0,0.0);
    g._gamma[0][2]=EvtComplex(-1.0,0.0);
    g._gamma[0][3]=EvtComplex(0.0,0.0);
    g._gamma[1][0]=EvtComplex(0.0,0.0);
    g._gamma[1][1]=EvtComplex(1.0,0.0);
    g._gamma[1][2]=EvtComplex(0.0,0.0);
    g._gamma[1][3]=EvtComplex(-1.0,0.0);
    g._gamma[2][0]=EvtComplex(-1.0,0.0);
    g._gamma[2][1]=EvtComplex(0.0,0.0);
    g._gamma[2][2]=EvtComplex(1.0,0.0);
    g._gamma[2][3]=EvtComplex(0.0,0.0);
    g._gamma[3][0]=EvtComplex(0.0,0.0);
    g._gamma[3][1]=EvtComplex(-1.0,0.0);
    g._gamma[3][2]=EvtComplex(0.0,0.0);
    g._gamma[3][3]=EvtComplex(1.0,0.0);
  }

  return g;

}


const EvtGammaMatrix& EvtGammaMatrix::va1(){

  static EvtGammaMatrix g;
  static int first=1;

  if (first){
    first = 0;
    g._gamma[0][0]=EvtComplex(0.0,0.0);
    g._gamma[0][1]=EvtComplex(-1.0,0.0);
    g._gamma[0][2]=EvtComplex(0.0,0.0);
    g._gamma[0][3]=EvtComplex(1.0,0.0);
    g._gamma[1][0]=EvtComplex(-1.0,0.0);
    g._gamma[1][1]=EvtComplex(0.0,0.0);
    g._gamma[1][2]=EvtComplex(1.0,0.0);
    g._gamma[1][3]=EvtComplex(0.0,0.0);
    g._gamma[2][0]=EvtComplex(0.0,0.0);
    g._gamma[2][1]=EvtComplex(1.0,0.0);
    g._gamma[2][2]=EvtComplex(0.0,0.0);
    g._gamma[2][3]=EvtComplex(-1.0,0.0);
    g._gamma[3][0]=EvtComplex(1.0,0.0);
    g._gamma[3][1]=EvtComplex(0.0,0.0);
    g._gamma[3][2]=EvtComplex(-1.0,0.0);
    g._gamma[3][3]=EvtComplex(0.0,0.0);
  }

  return g;

}



const EvtGammaMatrix& EvtGammaMatrix::va2(){

  static EvtGammaMatrix g;
  static int first=1;

  if (first){
    first = 0;
    g._gamma[0][0]=EvtComplex(0.0,0.0);
    g._gamma[0][1]=EvtComplex(0.0,1.0);
    g._gamma[0][2]=EvtComplex(0.0,0.0);
    g._gamma[0][3]=EvtComplex(0.0,-1.0);
    g._gamma[1][0]=EvtComplex(0.0,-1.0);
    g._gamma[1][1]=EvtComplex(0.0,0.0);
    g._gamma[1][2]=EvtComplex(0.0,1.0);
    g._gamma[1][3]=EvtComplex(0.0,0.0);
    g._gamma[2][0]=EvtComplex(0.0,0.0);
    g._gamma[2][1]=EvtComplex(0.0,-1.0);
    g._gamma[2][2]=EvtComplex(0.0,0.0);
    g._gamma[2][3]=EvtComplex(0.0,1.0);
    g._gamma[3][0]=EvtComplex(0.0,1.0);
    g._gamma[3][1]=EvtComplex(0.0,0.0);
    g._gamma[3][2]=EvtComplex(0.0,-1.0);
    g._gamma[3][3]=EvtComplex(0.0,0.0);
  }

  return g;

}




const EvtGammaMatrix& EvtGammaMatrix::va3(){

  static EvtGammaMatrix g;
  static int first=1;

  if (first){
    first = 0;
    g._gamma[0][0]=EvtComplex(-1.0,0.0);
    g._gamma[0][1]=EvtComplex(0.0,0.0);
    g._gamma[0][2]=EvtComplex(1.0,0.0);
    g._gamma[0][3]=EvtComplex(0.0,0.0);
    g._gamma[1][0]=EvtComplex(0.0,0.0);
    g._gamma[1][1]=EvtComplex(1.0,0.0);
    g._gamma[1][2]=EvtComplex(0.0,0.0);
    g._gamma[1][3]=EvtComplex(-1.0,0.0);
    g._gamma[2][0]=EvtComplex(1.0,0.0);
    g._gamma[2][1]=EvtComplex(0.0,0.0);
    g._gamma[2][2]=EvtComplex(-1.0,0.0);
    g._gamma[2][3]=EvtComplex(0.0,0.0);
    g._gamma[3][0]=EvtComplex(0.0,0.0);
    g._gamma[3][1]=EvtComplex(-1.0,0.0);
    g._gamma[3][2]=EvtComplex(0.0,0.0);
    g._gamma[3][3]=EvtComplex(1.0,0.0);
  }

  return g;
  
}





const EvtGammaMatrix& EvtGammaMatrix::g0(){

  static EvtGammaMatrix g;
  static int first=1;

  if (first){

    first=0;

    int i,j;
  
    for(i=0;i<4;i++){
      for(j=0;j<4;j++){
	g._gamma[i][j]=EvtComplex(0.0,0.0);
      }
    }
    
    g._gamma[0][0]=EvtComplex(1.0,0.0);
    g._gamma[1][1]=EvtComplex(1.0,0.0);
    g._gamma[2][2]=EvtComplex(-1.0,0.0);
    g._gamma[3][3]=EvtComplex(-1.0,0.0);
  }

  return g;

}




const EvtGammaMatrix& EvtGammaMatrix::g1(){

  static EvtGammaMatrix g;
  static int first=1;

  if (first){
    first=0;
    int i,j;
    
    for(i=0;i<4;i++){
      for(j=0;j<4;j++){
	g._gamma[i][j]=EvtComplex(0.0,0.0);
      }
    }
    
    g._gamma[0][3]=EvtComplex(1.0,0.0);
    g._gamma[1][2]=EvtComplex(1.0,0.0);
    g._gamma[2][1]=EvtComplex(-1.0,0.0);
    g._gamma[3][0]=EvtComplex(-1.0,0.0);
  }
  
  return g;

}




const EvtGammaMatrix& EvtGammaMatrix::g2(){

  static EvtGammaMatrix g;
  static int first=1;

  if (first){
    first=0;
    int i,j;
    
    for(i=0;i<4;i++){
      for(j=0;j<4;j++){
	g._gamma[i][j]=EvtComplex(0.0,0.0);
      }
    }
    
    g._gamma[0][3]=EvtComplex(0.0,-1.0);
    g._gamma[1][2]=EvtComplex(0.0,1.0);
    g._gamma[2][1]=EvtComplex(0.0,1.0);
    g._gamma[3][0]=EvtComplex(0.0,-1.0);
  }
  
  return g;

}





const EvtGammaMatrix& EvtGammaMatrix::g3(){

  static EvtGammaMatrix g;
  static int first=1;

  if (first){
    first=0;
    int i,j;
    
    for(i=0;i<4;i++){
      for(j=0;j<4;j++){
	g._gamma[i][j]=EvtComplex(0.0,0.0);
      }
    }
    
    g._gamma[0][2]=EvtComplex(1.0,0.0);
    g._gamma[1][3]=EvtComplex(-1.0,0.0);
    g._gamma[2][0]=EvtComplex(-1.0,0.0);
    g._gamma[3][1]=EvtComplex(1.0,0.0);
  }

  return g;

}




const EvtGammaMatrix& EvtGammaMatrix::g5(){

  static EvtGammaMatrix g;
  static int first=1;

  if (first){
    first = 0;
    int i,j;
    
    for(i=0;i<4;i++){
      for(j=0;j<4;j++){
	g._gamma[i][j]=EvtComplex(0.0,0.0);
      }
    }
    
    g._gamma[0][2]=EvtComplex(1.0,0.0);
    g._gamma[1][3]=EvtComplex(1.0,0.0);
    g._gamma[2][0]=EvtComplex(1.0,0.0);
    g._gamma[3][1]=EvtComplex(1.0,0.0);
  }

  return g;

}




const EvtGammaMatrix& EvtGammaMatrix::g(int index) {
  switch (index) {
  case 0:
    return g0();
  case 1:
    return g1();
  case 2:
    return g2();
  case 3:
    return g3();
  case 5:
    return g5();
  default:
    report(Severity::Error, "EvtGen") << "Invalid index for four vector: " << index << endl;
    exit(-2);
  }
}



const EvtGammaMatrix& EvtGammaMatrix::v0(){

  static EvtGammaMatrix g;
  static int first=1;

  if (first){
    first = 0;
    int i,j;
  
    for(i=0;i<4;i++){
      for(j=0;j<4;j++){
	g._gamma[i][j]=EvtComplex(0.0,0.0);
      }
    }
    
    g._gamma[0][0]=EvtComplex(1.0,0.0);
    g._gamma[1][1]=EvtComplex(1.0,0.0);
    g._gamma[2][2]=EvtComplex(1.0,0.0);
    g._gamma[3][3]=EvtComplex(1.0,0.0);
  }

  return g;

}





const EvtGammaMatrix& EvtGammaMatrix::v1(){

  static EvtGammaMatrix g;
  static int first=1;

  if (first){
    first = 0;
    int i,j;
    
    for(i=0;i<4;i++){
      for(j=0;j<4;j++){
	g._gamma[i][j]=EvtComplex(0.0,0.0);
      }
    }
    
    g._gamma[0][3]=EvtComplex(1.0,0.0);
    g._gamma[1][2]=EvtComplex(1.0,0.0);
    g._gamma[2][1]=EvtComplex(1.0,0.0);
    g._gamma[3][0]=EvtComplex(1.0,0.0);
  }

  return g;

}




const EvtGammaMatrix& EvtGammaMatrix::v2(){

  static EvtGammaMatrix g;
  static int first=1;

  if (first){
    first = 0;
    int i,j;

    for(i=0;i<4;i++){
      for(j=0;j<4;j++){
	g._gamma[i][j]=EvtComplex(0.0,0.0);
      }
    }
    
    g._gamma[0][3]=EvtComplex(0.0,-1.0);
    g._gamma[1][2]=EvtComplex(0.0,1.0);
    g._gamma[2][1]=EvtComplex(0.0,-1.0);
    g._gamma[3][0]=EvtComplex(0.0,1.0);
  }

  return g;

}




const EvtGammaMatrix& EvtGammaMatrix::v3(){

  static EvtGammaMatrix g;
  static int first=1;

  if (first){
    first = 0;
    int i,j;
  
    for(i=0;i<4;i++){
      for(j=0;j<4;j++){
	g._gamma[i][j]=EvtComplex(0.0,0.0);
      }
    }
    
    g._gamma[0][2]=EvtComplex(1.0,0.0);
    g._gamma[1][3]=EvtComplex(-1.0,0.0);
    g._gamma[2][0]=EvtComplex(1.0,0.0);
    g._gamma[3][1]=EvtComplex(-1.0,0.0);
  }

  return g;

}





const EvtGammaMatrix& EvtGammaMatrix::id(){

  static EvtGammaMatrix g;
  static int first=1;

  if (first){
    first = 0;
    int i,j;
    
    for(i=0;i<4;i++){
      for(j=0;j<4;j++){
	g._gamma[i][j]=EvtComplex(0.0,0.0);
      }
    }
    
    g._gamma[0][0]=EvtComplex(1.0,0.0);
    g._gamma[1][1]=EvtComplex(1.0,0.0);
    g._gamma[2][2]=EvtComplex(1.0,0.0);
    g._gamma[3][3]=EvtComplex(1.0,0.0);
  }

  return g;

}




EvtGammaMatrix& EvtGammaMatrix::operator+=(const EvtGammaMatrix &g){

  int i,j;
  
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      _gamma[i][j]+=g._gamma[i][j];
    }
  }
  return *this;
}





EvtGammaMatrix& EvtGammaMatrix::operator-=(const EvtGammaMatrix &g){

  int i,j;
  
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      _gamma[i][j]-=g._gamma[i][j];
    }
  }
  return *this;
}



EvtGammaMatrix& EvtGammaMatrix::operator*=(const EvtGammaMatrix &g){

  int i,j,k;
  EvtGammaMatrix temp;

  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      temp._gamma[i][j]=EvtComplex(0.0,0.0);
      for(k=0;k<4;k++){
	temp._gamma[i][j]+=_gamma[i][k]*g._gamma[k][j];
      }
    }
  }

  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
       _gamma[i][j]=temp._gamma[i][j];
    }
  }

  return *this;
}


EvtDiracSpinor operator*(const EvtGammaMatrix& g,const EvtDiracSpinor& d){

  int i,j;
  EvtDiracSpinor temp;
  
   for(i=0;i<4;i++){
     temp.set_spinor(i,EvtComplex(0.0,0.0));
     for(j=0;j<4;j++){
       temp.set_spinor(i,temp.get_spinor(i)+g._gamma[i][j]*d.get_spinor(j));
     }
   }
   
   return temp;
}

// upper index
const EvtGammaMatrix& EvtGammaMatrix::sigmaUpper(unsigned int mu, unsigned int nu)
{
    EvtGammaMatrix a, b;
    static const EvtTensor4C eta = EvtTensor4C::g(); //metric
    static EvtGammaMatrix sigma[4][4];
    static bool hasBeenCalled = false;
    if (!hasBeenCalled)
    {
        EvtComplex I(0, 1);
        for (int i=0; i<4; ++i)
            sigma[i][i].init(); // set to 0
        
        EvtGammaMatrix s01 = I/2 * (g0()*g1() - g1()*g0());
        EvtGammaMatrix s02 = I/2 * (g0()*g2() - g2()*g0());
        EvtGammaMatrix s03 = I/2 * (g0()*g3() - g3()*g0());
        EvtGammaMatrix s12 = I/2 * (g1()*g2() - g2()*g1());
        EvtGammaMatrix s13 = I/2 * (g1()*g3() - g3()*g1());
        EvtGammaMatrix s23 = I/2 * (g2()*g3() - g3()*g2());
        sigma[0][1] = s01;
        sigma[1][0] = -1*s01;
        sigma[0][2] = s02;
        sigma[2][0] = -1*s02;
        sigma[0][3] = s03;
        sigma[3][0] = -1*s03;
        sigma[1][2] = s12;
        sigma[2][1] = -1*s12;
        sigma[1][3] = s13;
        sigma[3][1] = -1*s13;
        sigma[2][3] = s23;
        sigma[3][2] = -1*s23;
    }
    hasBeenCalled = true;
        
    if (mu > 3 || nu > 3)
    {
        report(Severity::Error, "EvtSigmaTensor") << "Expected index between 0 and 3, but found " << nu << "!" << endl;
        assert(0);
    }
    return sigma[mu][nu];
    
}

const EvtGammaMatrix& EvtGammaMatrix::sigmaLower(unsigned int mu, unsigned int nu)
{
    const EvtComplex I(0, 1);
    EvtGammaMatrix a, b;
    static EvtGammaMatrix sigma[4][4];
    static bool hasBeenCalled = false;
    static const EvtTensor4C eta = EvtTensor4C::g();
    
    if (!hasBeenCalled) // has to be initialized only at the first call
    {
        // lower index
        for (int i=0; i<4; ++i)
        {
            a = eta.get(i, 0)*g0() + eta.get(i, 1)*g1() + eta.get(i, 2)*g2() + eta.get(i, 3)*g3();
            for (int j=0; j<4; ++j)
            {
                b = eta.get(j, 0)*g0() + eta.get(j, 1)*g1() + eta.get(j, 2)*g2() + eta.get(j, 3)*g3();
                sigma[i][j] = I/2 * (a*b - b*a);
            }
        }
    }
    return sigma[mu][nu];    
}


EvtGammaMatrix EvtGenFunctions::slash(const EvtVector4C& p)
{
    return EvtGammaMatrix::g0()*p.get(0) + 
      EvtGammaMatrix::g1()*p.get(1) + 
      EvtGammaMatrix::g2()*p.get(2) + 
      EvtGammaMatrix::g3()*p.get(3);
}

EvtGammaMatrix EvtGenFunctions::slash(const EvtVector4R& p)
{
  return EvtGammaMatrix::g0()*p.get(0) + 
    EvtGammaMatrix::g1()*p.get(1) + 
    EvtGammaMatrix::g2()*p.get(2) + 
    EvtGammaMatrix::g3()*p.get(3);
}
