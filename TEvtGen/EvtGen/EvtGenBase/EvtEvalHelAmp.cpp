//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2000      Caltech, UCSB
//
// Module: EvtHelAmp.cc
//
// Description: Decay model for implementation of generic 2 body
//              decay specified by the partial wave amplitudes
//
//
// Modification history:
//
//    fkw        February 2, 2001     changes to satisfy KCC
//    RYD       September 7, 2000       Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtEvalHelAmp.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtdFunction.hh"
#include "EvtGenBase/EvtAmp.hh"
using std::endl;


EvtEvalHelAmp::~EvtEvalHelAmp() {

  //deallocate memory
  delete [] _lambdaA2;
  delete [] _lambdaB2;
  delete [] _lambdaC2;

  int ia,ib,ic;
  for(ib=0;ib<_nB;ib++){
    delete [] _HBC[ib];
  }

  delete [] _HBC;


  for(ia=0;ia<_nA;ia++){
    delete [] _RA[ia];
  }
  delete [] _RA;

  for(ib=0;ib<_nB;ib++){
    delete [] _RB[ib];
  }
  delete [] _RB;

  for(ic=0;ic<_nC;ic++){
    delete [] _RC[ic];
  }
  delete [] _RC;

 
  for(ia=0;ia<_nA;ia++){
    for(ib=0;ib<_nB;ib++){
      delete [] _amp[ia][ib];
      delete [] _amp1[ia][ib];
      delete [] _amp3[ia][ib];
    }
    delete [] _amp[ia];
    delete [] _amp1[ia];
    delete [] _amp3[ia];
  }

  delete [] _amp;
  delete [] _amp1;
  delete [] _amp3;

}


EvtEvalHelAmp::EvtEvalHelAmp(EvtId idA,
			     EvtId idB,
			     EvtId idC,
			     EvtComplexPtrPtr HBC){


  EvtSpinType::spintype typeA=EvtPDL::getSpinType(idA);
  EvtSpinType::spintype typeB=EvtPDL::getSpinType(idB);
  EvtSpinType::spintype typeC=EvtPDL::getSpinType(idC);

  //find out how many states each particle have
  _nA=EvtSpinType::getSpinStates(typeA);
  _nB=EvtSpinType::getSpinStates(typeB);
  _nC=EvtSpinType::getSpinStates(typeC);

  //find out what 2 times the spin is
  _JA2=EvtSpinType::getSpin2(typeA);
  _JB2=EvtSpinType::getSpin2(typeB);
  _JC2=EvtSpinType::getSpin2(typeC);


  //allocate memory
  _lambdaA2=new int[_nA];
  _lambdaB2=new int[_nB];
  _lambdaC2=new int[_nC];

  _HBC=new EvtComplexPtr[_nB];
  int ia,ib,ic;
  for(ib=0;ib<_nB;ib++){
    _HBC[ib]=new EvtComplex[_nC];
  }


  _RA=new EvtComplexPtr[_nA];
  for(ia=0;ia<_nA;ia++){
    _RA[ia]=new EvtComplex[_nA];
  }
  _RB=new EvtComplexPtr[_nB];
  for(ib=0;ib<_nB;ib++){
    _RB[ib]=new EvtComplex[_nB];
  }
  _RC=new EvtComplexPtr[_nC];
  for(ic=0;ic<_nC;ic++){
    _RC[ic]=new EvtComplex[_nC];
  }
  
  _amp=new EvtComplexPtrPtr[_nA];
  _amp1=new EvtComplexPtrPtr[_nA];
  _amp3=new EvtComplexPtrPtr[_nA];
  for(ia=0;ia<_nA;ia++){
    _amp[ia]=new EvtComplexPtr[_nB];
    _amp1[ia]=new EvtComplexPtr[_nB];
    _amp3[ia]=new EvtComplexPtr[_nB];
    for(ib=0;ib<_nB;ib++){
      _amp[ia][ib]=new EvtComplex[_nC];
      _amp1[ia][ib]=new EvtComplex[_nC];
      _amp3[ia][ib]=new EvtComplex[_nC];
    }
  }

  //find the allowed helicities (actually 2*times the helicity!)

  fillHelicity(_lambdaA2,_nA,_JA2,idA);
  fillHelicity(_lambdaB2,_nB,_JB2,idB);
  fillHelicity(_lambdaC2,_nC,_JC2,idC);

  for(ib=0;ib<_nB;ib++){
    for(ic=0;ic<_nC;ic++){
      _HBC[ib][ic]=HBC[ib][ic];
    }
  }
}






double EvtEvalHelAmp::probMax(){

  double c=1.0/sqrt(4*EvtConst::pi/(_JA2+1));

  int ia,ib,ic;


  double theta;
  int itheta;

  double maxprob=0.0;

  for(itheta=-10;itheta<=10;itheta++){
    theta=acos(0.099999*itheta);
    for(ia=0;ia<_nA;ia++){
      double prob=0.0;
      for(ib=0;ib<_nB;ib++){
	for(ic=0;ic<_nC;ic++){
	  _amp[ia][ib][ic]=0.0;
	  if (abs(_lambdaB2[ib]-_lambdaC2[ic])<=_JA2) {
	    _amp[ia][ib][ic]=c*_HBC[ib][ic]*
	      EvtdFunction::d(_JA2,_lambdaA2[ia],
			      _lambdaB2[ib]-_lambdaC2[ic],theta);
	    prob+=real(_amp[ia][ib][ic]*conj(_amp[ia][ib][ic]));
	  }
	}
      }
      
      prob*=sqrt(1.0*_nA);
      
      if (prob>maxprob) maxprob=prob;

    }
  }

  return maxprob;

}


void EvtEvalHelAmp::evalAmp( EvtParticle *p, EvtAmp& amp){

  //find theta and phi of the first daughter
  
  EvtVector4R pB=p->getDaug(0)->getP4();

  double theta=acos(pB.get(3)/pB.d3mag());
  double phi=atan2(pB.get(2),pB.get(1));

  double c=sqrt((_JA2+1)/(4*EvtConst::pi));

  int ia,ib,ic;

  double prob1=0.0;

  for(ia=0;ia<_nA;ia++){
    for(ib=0;ib<_nB;ib++){
      for(ic=0;ic<_nC;ic++){
	_amp[ia][ib][ic]=0.0;
	if (abs(_lambdaB2[ib]-_lambdaC2[ic])<=_JA2) {
	  double dfun=EvtdFunction::d(_JA2,_lambdaA2[ia],
				      _lambdaB2[ib]-_lambdaC2[ic],theta);

	  _amp[ia][ib][ic]=c*_HBC[ib][ic]*
	    exp(EvtComplex(0.0,phi*0.5*(_lambdaA2[ia]-_lambdaB2[ib]+
	  				_lambdaC2[ic])))*dfun;
	}
	prob1+=real(_amp[ia][ib][ic]*conj(_amp[ia][ib][ic]));
      }
    }
  }

  setUpRotationMatrices(p,theta,phi);

  applyRotationMatrices();

  double prob2=0.0;

  for(ia=0;ia<_nA;ia++){
    for(ib=0;ib<_nB;ib++){
      for(ic=0;ic<_nC;ic++){
	prob2+=real(_amp[ia][ib][ic]*conj(_amp[ia][ib][ic]));
	if (_nA==1){
	  if (_nB==1){
	    if (_nC==1){
	      amp.vertex(_amp[ia][ib][ic]);
	    }
	    else{
	      amp.vertex(ic,_amp[ia][ib][ic]);
	    }
	  }
	  else{
	    if (_nC==1){
	      amp.vertex(ib,_amp[ia][ib][ic]);
	    }
	    else{
	      amp.vertex(ib,ic,_amp[ia][ib][ic]);
	    }
	  }
	}else{
	  if (_nB==1){
	    if (_nC==1){
	      amp.vertex(ia,_amp[ia][ib][ic]);
	    }
	    else{
	      amp.vertex(ia,ic,_amp[ia][ib][ic]);
	    }
	  }
	  else{
	    if (_nC==1){
	      amp.vertex(ia,ib,_amp[ia][ib][ic]);
	    }
	    else{
	      amp.vertex(ia,ib,ic,_amp[ia][ib][ic]);
	    }
	  }
	}
      }
    }
  }

  if (fabs(prob1-prob2)>0.000001*prob1){
    report(Severity::Info,"EvtGen") << "prob1,prob2:"<<prob1<<" "<<prob2<<endl;
    ::abort();
  }
    
  return ;

}


void EvtEvalHelAmp::fillHelicity(int* lambda2,int n,int J2, EvtId id){
  
  int i;
  
  //photon is special case!
  if (n==2&&J2==2) {
    lambda2[0]=2;
    lambda2[1]=-2;
    return;
  }
  
  //and so is the neutrino!
  if (n==1&&J2==1) {
    if (EvtPDL::getStdHep(id)>0){
	//particle i.e. lefthanded
        lambda2[0]=-1;
    }else{
	//anti particle i.e. righthanded
        lambda2[0]=1;
    }
    return;
  }


  assert(n==J2+1);

  for(i=0;i<n;i++){
    lambda2[i]=n-i*2-1;
  }

  return;

}


void EvtEvalHelAmp::setUpRotationMatrices(EvtParticle* p,double theta, double phi){

  switch(_JA2){

  case 0: case 1: case 2: case 3: case 4: case 5: case 6: case 7: case 8:

    {

      EvtSpinDensity R=p->rotateToHelicityBasis();

      
      int i,j,n;
      
      n=R.getDim();
      
      assert(n==_nA);
	
      
      for(i=0;i<n;i++){
	for(j=0;j<n;j++){
	  _RA[i][j]=R.get(i,j);
	}
      }

    }

    break;

  default:
    report(Severity::Error,"EvtGen") << "Spin2(_JA2)="<<_JA2<<" not supported!"<<endl;
    ::abort();
  }
  
  
  switch(_JB2){


  case 0: case 1: case 2: case 3: case 4: case 5: case 6: case 7: case 8:

    {
      
      int i,j,n;

      EvtSpinDensity R=p->getDaug(0)->rotateToHelicityBasis(phi,theta,-phi);
      
      n=R.getDim();
      
      assert(n==_nB);
	
      for(i=0;i<n;i++){
	for(j=0;j<n;j++){
	  _RB[i][j]=conj(R.get(i,j));
	}
      }

    }

    break;

  default:
    report(Severity::Error,"EvtGen") << "Spin2(_JB2)="<<_JB2<<" not supported!"<<endl;
    ::abort();
  }
  
  switch(_JC2){

  case 0: case 1: case 2: case 3: case 4: case 5: case 6: case 7: case 8:

    {

      int i,j,n;

      EvtSpinDensity R=p->getDaug(1)->rotateToHelicityBasis(phi,EvtConst::pi+theta,phi-EvtConst::pi);
            
      n=R.getDim();

      assert(n==_nC);

      for(i=0;i<n;i++){
	for(j=0;j<n;j++){
	  _RC[i][j]=conj(R.get(i,j));
	}
      }

    }

    break;

  default:
    report(Severity::Error,"EvtGen") << "Spin2(_JC2)="<<_JC2<<" not supported!"<<endl;
    ::abort();
  }
  
  

}


void EvtEvalHelAmp::applyRotationMatrices(){

  int ia,ib,ic,i;
  
  EvtComplex temp;



  for(ia=0;ia<_nA;ia++){
    for(ib=0;ib<_nB;ib++){
      for(ic=0;ic<_nC;ic++){
	temp=0;
	for(i=0;i<_nC;i++){
	  temp+=_RC[i][ic]*_amp[ia][ib][i];
	}
	_amp1[ia][ib][ic]=temp;
      }
    }
  }



  for(ia=0;ia<_nA;ia++){
    for(ic=0;ic<_nC;ic++){
      for(ib=0;ib<_nB;ib++){
  	temp=0;
  	for(i=0;i<_nB;i++){
  	  temp+=_RB[i][ib]*_amp1[ia][i][ic];
  	}
  	_amp3[ia][ib][ic]=temp;
      }
    }
  }
  


  for(ib=0;ib<_nB;ib++){
    for(ic=0;ic<_nC;ic++){
      for(ia=0;ia<_nA;ia++){
	temp=0;
	for(i=0;i<_nA;i++){
	  temp+=_RA[i][ia]*_amp3[i][ib][ic];
	}
	_amp[ia][ib][ic]=temp;
      }
    }
  }


}












