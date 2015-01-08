//-----------------------------------------------------------------------
// File and Version Information: 
//      $Id: EvtPto3PAmpFactory.cpp,v 1.3 2009-03-16 15:44:04 robbep Exp $
// 
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations. If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information:
//      Copyright (C) 1998 Caltech, UCSB
//
// Module creator:
//      Alexei Dvoretskii, Caltech, 2001-2002.
//-----------------------------------------------------------------------
#include "EvtGenBase/EvtPatches.hh"

// AmpFactory for building a P -> 3P decay
// (pseudoscalar to three pseudoscalars)

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtCyclic3.hh"
#include "EvtGenBase/EvtSpinType.hh"
#include "EvtGenBase/EvtPto3PAmp.hh"
#include "EvtGenBase/EvtNonresonantAmp.hh"
#include "EvtGenBase/EvtFlatAmp.hh"
#include "EvtGenBase/EvtLASSAmp.hh"
#include "EvtGenBase/EvtPto3PAmpFactory.hh"
#include "EvtGenBase/EvtPropBreitWigner.hh"
#include "EvtGenBase/EvtPropFlatte.hh"
#include "EvtGenBase/EvtPropBreitWignerRel.hh"
#include "EvtGenBase/EvtDalitzResPdf.hh"
#include "EvtGenBase/EvtDalitzFlatPdf.hh"

using namespace EvtCyclic3;
#include <iostream>

void EvtPto3PAmpFactory::processAmp(EvtComplex c, std::vector<std::string> vv, bool conj)
{
  if(_verbose) {
    
    printf("Make %samplitude\n",conj ? "CP conjugate" : "");
    unsigned i;
    for(i=0;i<vv.size();i++) printf("%s\n",vv[i].c_str());
    printf("\n");
  }
  
  EvtAmplitude<EvtDalitzPoint>* amp = 0;
  EvtPdf<EvtDalitzPoint>* pdf = 0;
  std::string name;
  Pair pairRes=AB;

  size_t i;
  /*
         Experimental amplitudes
  */
  if(vv[0] == "PHASESPACE") {
    
    pdf = new EvtDalitzFlatPdf(_dp);
    amp = new EvtFlatAmp<EvtDalitzPoint>();
    name = "NR";
  }
  else if (!vv[0].find("NONRES")) {
    double alpha=0;
    EvtPto3PAmp::NumType typeNRes=EvtPto3PAmp::NONRES;
    if      (vv[0]=="NONRES_LIN") {
      typeNRes=EvtPto3PAmp::NONRES_LIN;
      pairRes=strToPair(vv[1].c_str());
    }
    else if (vv[0]=="NONRES_EXP") {
      typeNRes=EvtPto3PAmp::NONRES_EXP;
      pairRes = strToPair(vv[1].c_str());
      alpha   = strtod(vv[2].c_str(),0);
    }
    else assert(0);
    pdf = new EvtDalitzFlatPdf(_dp);
    amp = new EvtNonresonantAmp( &_dp, typeNRes, pairRes, alpha);  
  }
  else if (vv[0]=="LASS" || vv[0]=="LASS_ELASTIC" || vv[0]=="LASS_RESONANT") {
    pairRes = strToPair(vv[1].c_str());
    double m0 = strtod(vv[2].c_str(),0);
    double g0 = strtod(vv[3].c_str(),0);
    double a  = strtod(vv[4].c_str(),0);
    double r  = strtod(vv[5].c_str(),0);
    double cutoff  = strtod(vv[6].c_str(),0);
    pdf = new EvtDalitzResPdf(_dp,m0,g0,pairRes);
    amp = new EvtLASSAmp( &_dp, pairRes, m0, g0, a, r, cutoff, vv[0]);
  }

  /*
      Resonant amplitudes
  */
  else if(vv[0] == "RESONANCE") {
    EvtPto3PAmp* partAmp = 0;
      
    // RESONANCE stanza
    
    pairRes = strToPair(vv[1].c_str());
    EvtSpinType::spintype spinR = EvtSpinType::SCALAR;
    double mR, gR;      
    name = vv[2];
    EvtId resId = EvtPDL::getId(vv[2]);
    if(_verbose) printf("Particles %s form %sresonance %s\n",
			vv[1].c_str(),vv[2].c_str(), conj ? "(conj) " : "");

    // If no valid particle name is given, assume that 
    // it is the spin, the mass and the width of the particle.
      
    if(resId.getId() == -1) {
	
      switch(atoi(vv[2].c_str())) {
	
      case 0: { spinR = EvtSpinType::SCALAR; break; }
      case 1: { spinR = EvtSpinType::VECTOR; break; }
      case 2: { spinR = EvtSpinType::TENSOR; break; }
      case 3: { spinR = EvtSpinType::SPIN3; break; }
      case 4: { spinR = EvtSpinType::SPIN4; break; }
      default: { assert(0); break; }
      }
	
      mR = strtod(vv[3].c_str(),0);
      gR = strtod(vv[4].c_str(),0);
      i = 4;
    }
    else {
      
      // For a valid particle get spin, mass and width
      
      spinR = EvtPDL::getSpinType(resId);
      mR = EvtPDL::getMeanMass(resId);
      gR = EvtPDL::getWidth(resId);
      i = 2;
      
      // It's possible to specify mass and width of a particle 
      // explicitly
      
      if(vv[3] != "ANGULAR") {
	
	if(_verbose) 
	  printf("Setting m(%s)=%s g(%s)=%s\n",
		 vv[2].c_str(),vv[3].c_str(),vv[2].c_str(),vv[4].c_str());

	mR = strtod(vv[3].c_str(),0);
	gR = strtod(vv[4].c_str(),0);
	i = 4;
      }
    }
    
    // ANGULAR stanza
    
    if(vv[++i] != "ANGULAR") {

      printf("%s instead of ANGULAR\n",vv[i].c_str());
      exit(0);
    }
    Pair pairAng = strToPair(vv[++i].c_str());
    if(_verbose) printf("Angle is measured between particles %s\n",vv[i].c_str());
      
    // TYPE stanza
    
    std::string typeName = vv[++i];
    assert(typeName == "TYPE");
    std::string type = vv[++i];
    if(_verbose) printf("Propagator type %s\n",vv[i].c_str());
    
    if(type == "NBW") {      

      EvtPropBreitWigner prop(mR,gR);
      partAmp = new EvtPto3PAmp(_dp,pairAng,pairRes,spinR,prop,EvtPto3PAmp::NBW);
    }
    else if(type == "RBW_ZEMACH") {
      
      EvtPropBreitWignerRel prop(mR,gR);
      partAmp = new EvtPto3PAmp(_dp,pairAng,pairRes,spinR,prop,EvtPto3PAmp::RBW_ZEMACH);
    }
    else if(type == "RBW_KUEHN") {
      
      EvtPropBreitWignerRel prop(mR,gR);
      partAmp = new EvtPto3PAmp(_dp,pairAng,pairRes,spinR,prop,EvtPto3PAmp::RBW_KUEHN);
    }
    else if(type == "RBW_CLEO") {
      
      EvtPropBreitWignerRel prop(mR,gR);
      partAmp = new EvtPto3PAmp(_dp,pairAng,pairRes,spinR,prop,EvtPto3PAmp::RBW_CLEO);
    }     
    else if(type == "FLATTE") {
      
      double m1a = _dp.m( first(pairRes) );
      double m1b = _dp.m( second(pairRes) );    
      // 2nd channel
      double g2  = strtod(vv[++i].c_str(),0);
      double m2a = strtod(vv[++i].c_str(),0);
      double m2b = strtod(vv[++i].c_str(),0);
      EvtPropFlatte  prop( mR, gR, m1a, m1b, g2, m2a, m2b );
      partAmp = new EvtPto3PAmp(_dp,pairAng,pairRes,spinR,prop,EvtPto3PAmp::FLATTE);
    }
    else assert(0);
      
    // Optional DVFF, BVFF stanzas
    
    if(i < vv.size() - 1) {
      if(vv[i+1] == "DVFF") {	
	i++;
	if(vv[++i] == "BLATTWEISSKOPF") {
	  
	  double R = strtod(vv[++i].c_str(),0);
	  partAmp->set_fd(R);
	}
	else assert(0);
      }
    }
      
    if(i < vv.size() - 1) {
      if(vv[i+1] == "BVFF") {	
	i++;
	if(vv[++i] == "BLATTWEISSKOPF") {

	  if(_verbose) printf("BVFF=%s\n",vv[i].c_str());
	  double R = strtod(vv[++i].c_str(),0);
	  partAmp->set_fb(R);
	}
	else assert(0);
      }
    }

    const int minwidths=5;
    //Optional resonance minimum and maximum
    if(i < vv.size() - 1) {
      if(vv[i+1] == "CUTOFF") {	
	i++;
	if(vv[i+1] == "MIN") {
	  i++;
	  double min = strtod(vv[++i].c_str(),0);
	  if(_verbose) std::cout<<"CUTOFF MIN = "<<min<<" "<<minwidths<<std::endl;
	  //ensure against cutting off too close to the resonance
	  assert( min<(mR-minwidths*gR) );
	  partAmp->setmin(min);
	}
	else if (vv[i+1] == "MAX") {
	  i++;
	  double max = strtod(vv[++i].c_str(),0);
	  if(_verbose) std::cout<<"CUTOFF MAX = "<<max<<" "<<minwidths<<std::endl;
	  //ensure against cutting off too close to the resonance
	  assert( max>(mR+minwidths*gR) );
	  partAmp->setmax(max);
	}
	else assert(0);
      }
    }

    //2nd iteration in case min and max are both specified
    if(i < vv.size() - 1) {
      if(vv[i+1] == "CUTOFF") {	
	i++;
	if(vv[i+1] == "MIN") {
	  i++;
	  double min = strtod(vv[++i].c_str(),0);
	  if(_verbose) std::cout<<"CUTOFF MIN = "<<min<<std::endl;
	  //ensure against cutting off too close to the resonance
	  assert( min<(mR-minwidths*gR) );
	  partAmp->setmin(min);
	}
	else if (vv[i+1] == "MAX") {
	  i++;
	  double max = strtod(vv[++i].c_str(),0);
	  if(_verbose) std::cout<<"CUTOFF MAX = "<<max<<std::endl;
	  //ensure against cutting off too close to the resonance
	  assert( max>(mR+minwidths*gR) );
	  partAmp->setmax(max);
	}
	else assert(0);
      }
    }


    i++;
    
    pdf = new EvtDalitzResPdf(_dp,mR,gR,pairRes);
    amp = partAmp;
  }

  assert(amp);
  assert(pdf);

  if(!conj) {
    _amp->addOwnedTerm(c,amp);
  }
  else {
    _ampConj->addOwnedTerm(c,amp);
  }

  double scale = matchIsobarCoef(_amp, pdf, pairRes);
  _pc->addOwnedTerm(abs2(c)*scale,pdf);

  _names.push_back(name);
}
  
double EvtPto3PAmpFactory::matchIsobarCoef(EvtAmplitude<EvtDalitzPoint>* amp,
					   EvtPdf<EvtDalitzPoint>* pdf, 
					   EvtCyclic3::Pair ipair) {

  // account for differences in the definition of amplitudes by matching 
  //        Integral( c'*pdf ) = Integral( c*|A|^2 ) 
  // to improve generation efficiency ...

  double Ipdf  = pdf->compute_integral(10000).value();
  double Iamp2 = 0;


  EvtCyclic3::Pair jpair = EvtCyclic3::next(ipair);
  EvtCyclic3::Pair kpair = EvtCyclic3::next(jpair);

  // Trapezoidal integral
  int N=10000;
  
  double di = (_dp.qAbsMax(ipair) - _dp.qAbsMin(ipair))/((double) N);
  
  double siMin = _dp.qAbsMin(ipair);
  
  double s[3]; // playing with fire
  for(int i=1; i<N; i++) {
    
    s[ipair] = siMin + di*i;
    s[jpair] = _dp.q(jpair, 0.9999, ipair, s[ipair]);    
    s[kpair] = _dp.bigM()*_dp.bigM() - s[ipair] - s[jpair]
      + _dp.mA()*_dp.mA() + _dp.mB()*_dp.mB() + _dp.mC()*_dp.mC();
    
    EvtDalitzPoint point( _dp.mA(), _dp.mB(), _dp.mC(), 
			  s[EvtCyclic3::AB], s[EvtCyclic3::BC], s[EvtCyclic3::CA]);
    
    if (!point.isValid()) continue;
    
    double p = point.p(other(ipair), ipair);
    double q = point.p(first(ipair), ipair);
    
    double itg = abs2( amp->evaluate(point) )*di*4*q*p;
    Iamp2 += itg;
    
  }
  if (_verbose) std::cout << "integral = " << Iamp2 << "  pdf="<<Ipdf << std::endl;
  
  assert(Ipdf>0 && Iamp2>0);
  
  return Iamp2/Ipdf;
}
