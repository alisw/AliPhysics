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
// Module: EvtSpinDensity.cc
//
// Description: Class to reperesent spindensity matrices.
//
// Modification history:
//
//    RYD       May 29,1997       Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <assert.h>
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtSpinDensity.hh"
#include "EvtGenBase/EvtReport.hh"
using std::endl;
using std::ostream;


EvtSpinDensity::EvtSpinDensity(const EvtSpinDensity& density){
  dim=0;
  rho=0;

  int i,j;
  setDim(density.dim);

  for(i=0;i<dim;i++){
    for(j=0;j<dim;j++){
      rho[i][j]=density.rho[i][j];
    }
  }
}

EvtSpinDensity& EvtSpinDensity::operator=(const EvtSpinDensity& density){
  int i,j;
  setDim(density.dim);

  for(i=0;i<dim;i++){
    for(j=0;j<dim;j++){
      rho[i][j]=density.rho[i][j];
    }
  }

  return *this;

}

EvtSpinDensity::~EvtSpinDensity(){

  if (dim!=0){
    int i;
    for(i=0;i<dim;i++) delete [] rho[i];
  }
  
  delete [] rho;

}

EvtSpinDensity::EvtSpinDensity(){
  dim=0;
  rho=0;
}

void EvtSpinDensity::setDim(int n){
  if (dim==n) return;
  if (dim!=0){
    int i;
    for(i=0;i<dim;i++) delete [] rho[i];
    delete [] rho;
    rho=0;
    dim=0;
  }
  if (n==0) return;
  dim=n;
  rho=new EvtComplexPtr[n];
  int i;
  for(i=0;i<n;i++){
    rho[i]=new EvtComplex[n];
  }


}

int EvtSpinDensity::getDim() const {
  return dim;
}

void EvtSpinDensity::set(int i,int j,const EvtComplex& rhoij){
  assert(i<dim&&j<dim);
  rho[i][j]=rhoij;
}

const EvtComplex& EvtSpinDensity::get(int i,int j)const{
  assert(i<dim&&j<dim);
  return rho[i][j];
}

void EvtSpinDensity::setDiag(int n){
  setDim(n);
  int i,j;

  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      rho[i][j]=EvtComplex(0.0);
    }
    rho[i][i]=EvtComplex(1.0);
  }
}

double EvtSpinDensity::normalizedProb(const EvtSpinDensity& d){

  int i,j;
  EvtComplex prob(0.0,0.0);
  double norm=0.0;

  if (dim!=d.dim) {
    report(Severity::Error,"EvtGen")<<"Not matching dimensions in NormalizedProb"<<endl;
    ::abort();
  }

  for(i=0;i<dim;i++){
    norm+=real(rho[i][i]);
    for(j=0;j<dim;j++){
      prob+=rho[i][j]*d.rho[i][j];
    }
  }

  if (imag(prob)>0.00000001*real(prob)) {
    report(Severity::Error,"EvtGen")<<"Imaginary probability:"<<prob<<" "<<norm<<endl;
  }
  if (real(prob)<0.0) {
    report(Severity::Error,"EvtGen")<<"Negative probability:"<<prob<<" "<<norm<<endl;
  }

  return real(prob)/norm;

}

int EvtSpinDensity::check(){

  if (dim<1) {
    report(Severity::Error,"EvtGen")<<"dim="<<dim<<"in SpinDensity::Check"<<endl;
  }

  int i,j;

  double trace(0.0);

  for (i=0;i<dim;i++) {
    trace += abs(rho[i][i]);
  }

  for(i=0;i<dim;i++){

    if (real(rho[i][i])<0.0) return 0;
    if (imag(rho[i][i])*1000000.0>trace) {
      report(Severity::Info,"EvtGen") << *this << endl;
      report(Severity::Info,"EvtGen") << trace << endl;
      report(Severity::Info,"EvtGen") << "Failing 1"<<endl;
      return 0;
    }
  }

  for(i=0;i<dim;i++){
    for(j=i+1;j<dim;j++){
      if (fabs(real(rho[i][j]-rho[j][i]))>
	  0.00000001*(abs(rho[i][i])+abs(rho[j][j]))) {
	report(Severity::Info,"EvtGen") << "Failing 2"<<endl;
	return 0;
      }
      if (fabs(imag(rho[i][j]+rho[j][i]))>
	  0.00000001*(abs(rho[i][i])+abs(rho[j][j]))) {
	report(Severity::Info,"EvtGen") << "Failing 3"<<endl;
	return 0;
      }
    }
  }

  return 1;
}

ostream& operator<<(ostream& s,const EvtSpinDensity& d){

  int i,j;

  s << endl;
  s << "Dimension:"<<d.dim<<endl;

  for (i=0;i<d.dim;i++){
    for (j=0;j<d.dim;j++){
     s << d.rho[i][j]<<" ";
    }
    s <<endl;
  }

  return s;

}



