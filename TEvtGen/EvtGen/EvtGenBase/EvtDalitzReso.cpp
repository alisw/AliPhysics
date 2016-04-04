#include "EvtGenBase/EvtPatches.hh"
/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtDalitzReso.cpp,v 1.1 2009-03-16 16:47:51 robbep Exp $
 *
 * Description:
 *   Class to compute Dalitz amplitudes based on many models that cannot be
 *     handled with EvtResonance.
 *
 * Modification history:
 *   Jordi Garra Ticó     2008/07/03         File created
 *****************************************************************************/


#include <assert.h>
#include <cmath>
#include <iostream>

#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtMatrix.hh"
#include "EvtGenBase/EvtDalitzReso.hh"

#include "EvtGenBase/EvtdFunction.hh"
#include "EvtGenBase/EvtCyclic3.hh"

#define PRECISION ( 1.e-3 )

using EvtCyclic3::Index;
using EvtCyclic3::Pair;


// single Breit-Wigner
EvtDalitzReso::EvtDalitzReso(const EvtDalitzPlot& dp, Pair pairAng, Pair pairRes, 
			     EvtSpinType::spintype spin, double m0, double g0, NumType typeN, double f_b, double f_d) 
  : _dp(dp),
    _pairAng(pairAng),
    _pairRes(pairRes),
    _spin(spin),
    _typeN(typeN),
    _m0(m0),_g0(g0),
    _massFirst(dp.m(first(pairRes))),_massSecond(dp.m(second(pairRes))),
    _m0_mix(-1.),_g0_mix(0.),_delta_mix(0.),_amp_mix(0.,0.),
    _g1(-1.),_g2(-1.),_coupling2(Undefined),
    _f_b(f_b), _f_d(f_d),
    _kmatrix_index(-1),_fr12prod(0.,0.),_fr13prod(0.,0.),_fr14prod(0.,0.),_fr15prod(0.,0.),_s0prod(0.),
    _a(0.),_r(0.),_Blass(0.),_phiB(0.),_R(0.),_phiR(0.),_cutoff(-1.), _scaleByMOverQ(false),
    _alpha(0.)
{
  _vb = EvtTwoBodyVertex(_m0,_dp.m(EvtCyclic3::other(_pairRes)),_dp.bigM(),_spin); 
  _vd = EvtTwoBodyVertex(_massFirst,_massSecond,_m0,_spin);
  _vb.set_f( _f_b ); // Default values for Blatt-Weisskopf factors are 0.0 and 1.5.
  _vd.set_f( _f_d );
  assert(_typeN != K_MATRIX && _typeN != K_MATRIX_I && _typeN != K_MATRIX_II);  // single BW cannot be K-matrix
}


// Breit-Wigner with electromagnetic mass mixing
EvtDalitzReso::EvtDalitzReso(const EvtDalitzPlot& dp, Pair pairAng, Pair pairRes, 
			     EvtSpinType::spintype spin, double m0, double g0, NumType typeN,
			     double m0_mix, double g0_mix, double delta_mix, EvtComplex amp_mix) 
  : _dp(dp),
    _pairAng(pairAng),
    _pairRes(pairRes),
    _spin(spin),
    _typeN(typeN),
    _m0(m0),_g0(g0),
    _massFirst(dp.m(first(pairRes))),_massSecond(dp.m(second(pairRes))),
    _m0_mix(m0_mix),_g0_mix(g0_mix),_delta_mix(delta_mix),_amp_mix(amp_mix),
    _g1(-1.),_g2(-1.),_coupling2(Undefined),
    _f_b(0.0), _f_d(1.5),
    _kmatrix_index(-1),_fr12prod(0.,0.),_fr13prod(0.,0.),_fr14prod(0.,0.),_fr15prod(0.,0.),_s0prod(0.),
    _a(0.),_r(0.),_Blass(0.),_phiB(0.),_R(0.),_phiR(0.),_cutoff(-1.), _scaleByMOverQ(false),
    _alpha(0.)
{
  _vb = EvtTwoBodyVertex(_m0,_dp.m(EvtCyclic3::other(_pairRes)),_dp.bigM(),_spin); 
  _vd = EvtTwoBodyVertex(_massFirst,_massSecond,_m0,_spin);
  _vb.set_f( 0.0 ); // Default values for Blatt-Weisskopf factors.
  _vd.set_f( 1.5 );
  // single BW (with electromagnetic mixing) cannot be K-matrix
  assert(_typeN != K_MATRIX && _typeN != K_MATRIX_I && _typeN != K_MATRIX_II);
}

// coupled Breit-Wigner
EvtDalitzReso::EvtDalitzReso(const EvtDalitzPlot& dp, Pair pairAng, Pair pairRes, 
			     EvtSpinType::spintype spin, double m0, NumType typeN, double g1, double g2, CouplingType coupling2)
  : _dp(dp),
    _pairAng(pairAng),
    _pairRes(pairRes),
    _spin(spin),
    _typeN(typeN),
    _m0(m0),_g0(-1.),
    _massFirst(dp.m(first(pairRes))),_massSecond(dp.m(second(pairRes))),
    _m0_mix(-1.),_g0_mix(0.),_delta_mix(0.),_amp_mix(0.,0.),
    _g1(g1),_g2(g2),_coupling2(coupling2),
    _f_b(0.0), _f_d(1.5),
    _kmatrix_index(-1),_fr12prod(0.,0.),_fr13prod(0.,0.),_fr14prod(0.,0.),_fr15prod(0.,0.),_s0prod(0.),
    _a(0.),_r(0.),_Blass(0.),_phiB(0.),_R(0.),_phiR(0.),_cutoff(-1.), _scaleByMOverQ(false),
    _alpha(0.)
{
  _vb = EvtTwoBodyVertex(_m0,_dp.m(EvtCyclic3::other(_pairRes)),_dp.bigM(),_spin);   
  _vd = EvtTwoBodyVertex(_massFirst,_massSecond,_m0,_spin);
  _vb.set_f( 0.0 ); // Default values for Blatt-Weisskopf factors.
  _vd.set_f( 1.5 );
  assert(_coupling2 != Undefined);
  assert(_typeN != K_MATRIX && _typeN != K_MATRIX_I && _typeN != K_MATRIX_II); // coupled BW cannot be K-matrix
  assert(_typeN != LASS);     // coupled BW cannot be LASS
  assert(_typeN != NBW);      // for coupled BW, only relativistic BW 
}


// K-Matrix (A&S)
EvtDalitzReso::EvtDalitzReso(const EvtDalitzPlot& dp, Pair pairRes, std::string nameIndex, NumType typeN,
			     EvtComplex fr12prod, EvtComplex fr13prod, EvtComplex fr14prod, EvtComplex fr15prod, double s0prod) 
  : _dp(dp),
    _pairRes(pairRes),
    _typeN(typeN),
    _m0(0.),_g0(0.),
    _massFirst(dp.m(first(pairRes))),_massSecond(dp.m(second(pairRes))),
    _m0_mix(-1.),_g0_mix(0.),_delta_mix(0.),_amp_mix(0.,0.),
    _g1(-1.),_g2(-1.),_coupling2(Undefined),
    _f_b(0.), _f_d(0.),
    _kmatrix_index(-1),_fr12prod(fr12prod),_fr13prod(fr13prod),_fr14prod(fr14prod),_fr15prod(fr15prod),_s0prod(s0prod),
    _a(0.),_r(0.),_Blass(0.),_phiB(0.),_R(0.),_phiR(0.),_cutoff(-1.), _scaleByMOverQ(false),
    _alpha(0.)
{
  assert(_typeN==K_MATRIX || _typeN==K_MATRIX_I || _typeN==K_MATRIX_II);
  _spin=EvtSpinType::SCALAR;
  if (nameIndex=="Pole1") _kmatrix_index=1;
  else if (nameIndex=="Pole2") _kmatrix_index=2;
  else if (nameIndex=="Pole3") _kmatrix_index=3;
  else if (nameIndex=="Pole4") _kmatrix_index=4;
  else if (nameIndex=="Pole5") _kmatrix_index=5;
  else if (nameIndex=="f11prod") _kmatrix_index=6;
  else assert(0);
}


// LASS parameterization
EvtDalitzReso::EvtDalitzReso(const EvtDalitzPlot& dp, Pair pairRes, 
			     double m0, double g0, double a, double r, double B, double phiB, double R, double phiR, double cutoff, bool scaleByMOverQ) 
  : _dp(dp),
    _pairRes(pairRes),
    _typeN(LASS),
    _m0(m0),_g0(g0),
    _massFirst(dp.m(first(pairRes))),_massSecond(dp.m(second(pairRes))),
    _m0_mix(-1.),_g0_mix(0.),_delta_mix(0.),_amp_mix(0.,0.),
    _g1(-1.),_g2(-1.),_coupling2(Undefined),
    _f_b(0.0), _f_d(1.5),
    _kmatrix_index(-1),_fr12prod(0.,0.),_fr13prod(0.,0.),_fr14prod(0.,0.),_fr15prod(0.,0.),_s0prod(0.),
    _a(a),_r(r),_Blass(B),_phiB(phiB),_R(R),_phiR(phiR), _cutoff(cutoff), _scaleByMOverQ(scaleByMOverQ),
    _alpha(0.)
{
  _spin=EvtSpinType::SCALAR;
  _vd = EvtTwoBodyVertex(_massFirst,_massSecond,_m0,_spin);
  _vd.set_f( 1.5 ); // Default values for Blatt-Weisskopf factors.
}


//Flatte
EvtDalitzReso::EvtDalitzReso(const EvtDalitzPlot& dp, EvtCyclic3::Pair pairRes, double m0)
  : _dp(dp),
    _pairRes(pairRes),
    _typeN(FLATTE),
    _m0(m0), _g0(0.),
    _massFirst(dp.m(first(pairRes))),_massSecond(dp.m(second(pairRes))),
    _m0_mix(-1.),_g0_mix(0.),_delta_mix(0.),_amp_mix(0.,0.),
    _g1(-1.),_g2(-1.),_coupling2(Undefined),
    _f_b(0.), _f_d(0.),
    _kmatrix_index(-1),_fr12prod(0.,0.),_fr13prod(0.,0.),_fr14prod(0.,0.),_fr15prod(0.,0.),_s0prod(0.),
    _a(0.),_r(0.),_Blass(0.),_phiB(0.),_R(0.),_phiR(0.),_cutoff(-1.), _scaleByMOverQ(false),
    _alpha(0.)
{
  _spin=EvtSpinType::SCALAR;
}


EvtDalitzReso::EvtDalitzReso(const EvtDalitzReso& other) 
  : _dp(other._dp),
    _pairAng(other._pairAng),
    _pairRes(other._pairRes),
    _spin(other._spin),
    _typeN(other._typeN),
    _m0(other._m0),_g0(other._g0),
    _vb(other._vb),_vd(other._vd),
    _massFirst(other._massFirst),_massSecond(other._massSecond),
    _m0_mix(other._m0_mix),_g0_mix(other._g0_mix),_delta_mix(other._delta_mix),_amp_mix(other._amp_mix),
    _g1(other._g1),_g2(other._g2),_coupling2(other._coupling2),
    _f_b(other._f_b), _f_d(other._f_d),
    _kmatrix_index(other._kmatrix_index),
    _fr12prod(other._fr12prod),_fr13prod(other._fr13prod),_fr14prod(other._fr14prod),_fr15prod(other._fr15prod),
    _s0prod(other._s0prod),
    _a(other._a),_r(other._r),_Blass(other._Blass),_phiB(other._phiB),_R(other._R),_phiR(other._phiR),_cutoff(other._cutoff), _scaleByMOverQ(other._scaleByMOverQ),
    _alpha(other._alpha),
    _flatteParams(other._flatteParams)
{}


EvtDalitzReso::~EvtDalitzReso()
{}


EvtComplex EvtDalitzReso::evaluate(const EvtDalitzPoint& x) 
{
  double m = sqrt(x.q(_pairRes));

  if (_typeN==NON_RES) 
    return EvtComplex(1.0,0.0);

  if (_typeN==NON_RES_LIN)
    return m*m;

  if (_typeN==NON_RES_EXP)
    return exp(-_alpha*m*m);

  // do use always hash table (speed up fitting)
  if (_typeN==K_MATRIX || _typeN==K_MATRIX_I || _typeN==K_MATRIX_II)
    return Fvector( m*m, _kmatrix_index );

  if (_typeN==LASS)
    return lass(m*m);

  if (_typeN==FLATTE)
    return flatte(m);

  EvtComplex amp(1.0,0.0);

  if (fabs(_dp.bigM() - x.bigM()) > 0.000001) {
    _vb = EvtTwoBodyVertex(_m0,_dp.m(EvtCyclic3::other(_pairRes)),x.bigM(),_spin);
    _vb.set_f(_f_b);
  }
  EvtTwoBodyKine vb(m,x.m(EvtCyclic3::other(_pairRes)),x.bigM());
  EvtTwoBodyKine vd(_massFirst,_massSecond,m);   

  EvtComplex prop(0,0);
  if (_typeN==NBW) {
    prop = propBreitWigner(_m0,_g0,m);
  } else if (_typeN==GAUSS_CLEO || _typeN==GAUSS_CLEO_ZEMACH) {
    prop = propGauss(_m0,_g0,m);
  } else {
    if (_coupling2==Undefined) {  
      // single BW
      double g = (_g0<=0. || _vd.pD()<=0.)? -_g0 : _g0*_vd.widthFactor(vd);  // running width
      if (_typeN==GS_CLEO || _typeN==GS_CLEO_ZEMACH) {
	// Gounaris-Sakurai (GS)
	assert(_massFirst==_massSecond);
	prop = propGounarisSakurai(_m0,fabs(_g0),_vd.pD(),m,g,vd.p());
      } else {
	// standard relativistic BW
	prop = propBreitWignerRel(_m0,g,m);
      }
    } else {    
      // coupled width BW
      EvtComplex G1,G2;
      switch (_coupling2) { 
      case PicPic: {
	G1 = _g1*_g1*psFactor(_massFirst,_massSecond,m);
	static double mPic = EvtPDL::getMass( EvtPDL::getId( "pi+" ) );
	G2 = _g2*_g2*psFactor(mPic,mPic,m);
	break;
      }
      case PizPiz: {
	G1 = _g1*_g1*psFactor(_massFirst,_massSecond,m);
	static double mPiz = EvtPDL::getMass( EvtPDL::getId( "pi0" ) );
	G2 = _g2*_g2*psFactor(mPiz,mPiz,m);
	break;
      }
      case PiPi: {
	G1 = _g1*_g1*psFactor(_massFirst,_massSecond,m);
	static double mPic = EvtPDL::getMass( EvtPDL::getId( "pi+" ) );
	static double mPiz = EvtPDL::getMass( EvtPDL::getId( "pi0" ) );
	G2 = _g2*_g2*psFactor(mPic,mPic,mPiz,mPiz,m);
	break;
      }
      case KcKc: {
	G1 = _g1*_g1*psFactor(_massFirst,_massSecond,m);
	static double mKc = EvtPDL::getMass( EvtPDL::getId( "K+" ) );
	G2 = _g2*_g2*psFactor(mKc,mKc,m);
	break;
      }
      case KzKz: {
	G1 = _g1*_g1*psFactor(_massFirst,_massSecond,m);
	static double mKz = EvtPDL::getMass( EvtPDL::getId( "K0" ) );
	G2 = _g2*_g2*psFactor(mKz,mKz,m);
	break;
      }
      case KK: {
	G1 = _g1*_g1*psFactor(_massFirst,_massSecond,m);
	static double mKc = EvtPDL::getMass( EvtPDL::getId( "K+" ) );
	static double mKz = EvtPDL::getMass( EvtPDL::getId( "K0" ) );
	G2 = _g2*_g2*psFactor(mKc,mKc,mKz,mKz,m);
	break;
      }
      case EtaPic: {
	G1 = _g1*_g1*psFactor(_massFirst,_massSecond,m);
	static double mEta = EvtPDL::getMass( EvtPDL::getId( "eta" ) );
	static double mPic = EvtPDL::getMass( EvtPDL::getId( "pi+" ) );
	G2 = _g2*_g2*psFactor(mEta,mPic,m);
	break;
      }
      case EtaPiz: {
	G1 = _g1*_g1*psFactor(_massFirst,_massSecond,m);
	static double mEta = EvtPDL::getMass( EvtPDL::getId( "eta" ) );
	static double mPiz = EvtPDL::getMass( EvtPDL::getId( "pi0" ) );
	G2 = _g2*_g2*psFactor(mEta,mPiz,m);
	break;
      }
      case PicPicKK: {
	static double mPic = EvtPDL::getMass( EvtPDL::getId( "pi+" ) );
	//G1 = _g1*_g1*psFactor(mPic,mPic,m);
	G1 = _g1*psFactor(mPic,mPic,m);
	static double mKc = EvtPDL::getMass( EvtPDL::getId( "K+" ) );
	static double mKz = EvtPDL::getMass( EvtPDL::getId( "K0" ) );
	//G2 = _g2*_g2*psFactor(mKc,mKc,mKz,mKz,m);
	G2 = _g2*psFactor(mKc,mKc,mKz,mKz,m);
	break;
      }
      default:
	std::cout << "EvtDalitzReso:evaluate(): PANIC, wrong coupling2 state." << std::endl;
	assert(0);
	break;
      }
      // calculate standard couple BW propagator
      if (_coupling2 != WA76)
	prop = _g1*propBreitWignerRelCoupled(_m0,G1,G2,m);
    } 
  }
  amp *= prop;

  // Compute form-factors (Blatt-Weisskopf penetration factor)
  amp *= _vb.formFactor(vb);  
  amp *= _vd.formFactor(vd);  

  // Compute numerator (angular distribution)
  amp *= numerator(x,vb,vd);  

  // Compute electromagnetic mass mixing factor
  if (_m0_mix>0.) {
    EvtComplex prop_mix;
    if (_typeN==NBW) {
      prop_mix = propBreitWigner(_m0_mix,_g0_mix,m);
    } else {
      assert(_g1<0.); // running width only
      double g_mix = _g0_mix*_vd.widthFactor(vd);
      prop_mix = propBreitWignerRel(_m0_mix,g_mix,m);
    }
    amp *= mixFactor(prop,prop_mix);
  }

  return amp;
}


EvtComplex EvtDalitzReso::psFactor(double & ma, double & mb, double& m)
{
  if (m>(ma+mb)) {
    EvtTwoBodyKine vd(ma,mb,m);
    return EvtComplex(0,2*vd.p()/m);
  } else { 
    // analytical continuation
    double s = m*m;
    double phaseFactor_analyticalCont = -0.5*(sqrt(4*ma*ma/s-1)+sqrt(4*mb*mb/s-1)); 
    return EvtComplex(phaseFactor_analyticalCont,0);
  }
}


EvtComplex EvtDalitzReso::psFactor(double & ma1,double & mb1, double & ma2, double & mb2, double& m)
{
  return 0.5*(psFactor(ma1,mb1,m)+psFactor(ma2,mb2,m));
}


EvtComplex EvtDalitzReso::propGauss(const double& m0, const double& s0, const double& m) 
{
  // Gaussian
  double gauss = 1./sqrt(EvtConst::twoPi)/s0*exp(-(m-m0)*(m-m0)/2./(s0*s0));
  return EvtComplex(gauss,0.);
}


EvtComplex EvtDalitzReso::propBreitWigner(const double& m0, const double& g0, const double& m) 
{
  // non-relativistic BW
  return sqrt(g0/EvtConst::twoPi)/(m-m0-EvtComplex(0.0,g0/2.));
}


EvtComplex EvtDalitzReso::propBreitWignerRel(const double& m0, const double& g0, const double& m) 
{
  // relativistic BW with real width
  return 1./(m0*m0-m*m-EvtComplex(0.,m0*g0));
}



EvtComplex EvtDalitzReso::propBreitWignerRel(const double& m0, const EvtComplex& g0, const double& m) 
{
  // relativistic BW with complex width
  return 1./(m0*m0-m*m-EvtComplex(0.,m0)*g0);
}


EvtComplex EvtDalitzReso::propBreitWignerRelCoupled(const double& m0, const EvtComplex& g1, const EvtComplex& g2, const double& m)
{
  // relativistic coupled BW
  return 1./(m0*m0-m*m-(g1+g2));
}

EvtComplex EvtDalitzReso::propGounarisSakurai(const double& m0, const double& g0, const double& k0,
					    const double& m, const double& g, const double& k) 
{
  // Gounaris-Sakurai parameterization of pi+pi- P wave. PRD, Vol61, 112002. PRL, Vol21, 244.
  // Expressions taken from BAD637v4, after fixing the imaginary part of the BW denominator: i M_R Gamma_R(s) --> i sqrt(s) Gamma_R(s) 
  return (1.+GS_d(m0,k0)*g0/m0)/(m0*m0-m*m-EvtComplex(0.,m*g)+GS_f(m0,g0,k0,m,k));
}


inline double EvtDalitzReso::GS_f(const double& m0, const double& g0, const double& k0, const double& m, const double& k) 
{
  // m: sqrt(s)
  // m0: nominal resonance mass
  // k: momentum of pion in resonance rest frame (at m)
  // k0: momentum of pion in resonance rest frame (at nominal resonance mass)
  return g0*m0*m0/(k0*k0*k0)*( k*k*(GS_h(m,k)-GS_h(m0,k0)) + (m0*m0-m*m)*k0*k0*GS_dhods(m0,k0) );
}

inline double EvtDalitzReso::GS_h(const double& m, const double& k) 
{return 2./EvtConst::pi*k/m*log((m+2.*k)/(2.*_massFirst)) ;}

inline double EvtDalitzReso::GS_dhods(const double& m0, const double& k0)  
{return GS_h(m0,k0)*( 0.125/(k0*k0) - 0.5/(m0*m0) ) + 0.5/(EvtConst::pi*m0*m0) ;}

inline double EvtDalitzReso::GS_d(const double& m0, const double& k0) 
{return 3./EvtConst::pi*_massFirst*_massFirst/(k0*k0)*log((m0+2.*k0)/(2.*_massFirst)) + 
   m0/(2.*EvtConst::pi*k0) - _massFirst*_massFirst*m0/(EvtConst::pi*k0*k0*k0) ;}


EvtComplex EvtDalitzReso::numerator(const EvtDalitzPoint& x, const EvtTwoBodyKine& vb, const EvtTwoBodyKine& vd) 
{
  EvtComplex ret(0.,0.);

  // Non-relativistic Breit-Wigner
  if(NBW == _typeN) {
    ret = angDep(x);
  }

  // Standard relativistic Zemach propagator
  else if(RBW_ZEMACH == _typeN) {
    ret = _vd.phaseSpaceFactor(vd,EvtTwoBodyKine::AB)*angDep(x);
  }

  // Standard relativistic Zemach propagator
  else if(RBW_ZEMACH2 == _typeN) {
    ret = _vd.phaseSpaceFactor(vd,EvtTwoBodyKine::AB)*_vb.phaseSpaceFactor(vb,EvtTwoBodyKine::AB)*angDep(x);
    if(_spin == EvtSpinType::VECTOR) {
      ret *= -4.;
    } else if(_spin == EvtSpinType::TENSOR) {
      ret *= 16./3.;
    } else if(_spin != EvtSpinType::SCALAR)
      assert(0);
  }

  // Kuehn-Santamaria normalization:
  else if(RBW_KUEHN == _typeN) {
    ret = _m0*_m0 * angDep(x);
  }  

  // CLEO amplitude 
  else if( ( RBW_CLEO        == _typeN ) || ( GS_CLEO           == _typeN ) ||
	   ( RBW_CLEO_ZEMACH == _typeN ) || ( GS_CLEO_ZEMACH    == _typeN ) ||
	   ( GAUSS_CLEO      == _typeN ) || ( GAUSS_CLEO_ZEMACH == _typeN)) {

    Index iA = other(_pairAng);           // A = other(BC)
    Index iB = common(_pairRes,_pairAng); // B = common(AB,BC)
    Index iC = other(_pairRes);           // C = other(AB)
    
    double M = x.bigM();
    double mA = x.m(iA);
    double mB = x.m(iB);
    double mC = x.m(iC);
    double qAB = x.q(combine(iA,iB));
    double qBC = x.q(combine(iB,iC));
    double qCA = x.q(combine(iC,iA));

    double M2 = M*M;
    double m02 = ((RBW_CLEO_ZEMACH == _typeN)||(GS_CLEO_ZEMACH == _typeN)||(GAUSS_CLEO_ZEMACH == _typeN))?  qAB : _m0*_m0;
    double mA2 = mA*mA;
    double mB2 = mB*mB;
    double mC2 = mC*mC;
    
    if (_spin == EvtSpinType::SCALAR) ret = EvtComplex(1.,0.);
    else if(_spin == EvtSpinType::VECTOR) {
      ret = qCA - qBC + (M2 - mC2)*(mB2 - mA2)/m02;
    } else if(_spin == EvtSpinType::TENSOR) {
      double x1 = qBC - qCA + (M2 - mC2)*(mA2 - mB2)/m02;       
      double x2 = M2 - mC2;      
      double x3 = qAB - 2*M2 - 2*mC2 + x2*x2/m02;      
      double x4 = mA2 - mB2;
      double x5 = qAB - 2*mB2 - 2*mA2 + x4*x4/m02;
      ret = x1*x1 - x3*x5/3.;
    } else assert(0);
  }
  
  return ret;
}


double EvtDalitzReso::angDep(const EvtDalitzPoint& x)  
{ 
  // Angular dependece for factorizable amplitudes  
  // unphysical cosines indicate we are in big trouble
  double cosTh = x.cosTh(_pairAng,_pairRes);  // angle between common(reso,ang) and other(reso)
  if(fabs(cosTh) > 1.) {
    report(Severity::Info,"EvtGen") << "cosTh " << cosTh << std::endl; 
    assert(0);
  }
  
  // in units of half-spin
  return EvtdFunction::d(EvtSpinType::getSpin2(_spin),2*0,2*0,acos(cosTh));
}


EvtComplex EvtDalitzReso::mixFactor(EvtComplex prop, EvtComplex prop_mix) 
{
  double Delta = _delta_mix*(_m0+_m0_mix);
  return 1/(1-Delta*Delta*prop*prop_mix)*(1+_amp_mix*Delta*prop_mix);
}



EvtComplex EvtDalitzReso::Fvector( double s, int index )
{
  assert(index>=1 && index<=6);

  //Define the complex coupling constant
  //The convection is as follow
  //i=0 --> pi+ pi-
  //i=1 --> KK
  //i=2 --> 4pi
  //i=3 --> eta eta
  //i=4 --> eta eta'
  //The first index is the resonace-pole index
      
  double g[5][5]; // Coupling constants. The first index is the pole index. The second index is the decay channel
  double ma[5];   // Pole masses. The unit is in GeV

  int solution = (_typeN==K_MATRIX)? 3 : (    (_typeN==K_MATRIX_I)? 1 : ( (_typeN==K_MATRIX_II)? 2 : 0 )    ) ;
  if (solution==0) { std::cout << "EvtDalitzReso::Fvector() error. Kmatrix solution incorrectly chosen ! " << std::endl; abort(); } 

  if (solution == 3 ) {

    // coupling constants
    //pi+pi- channel
    g[0][0]=0.22889;
    g[1][0]=0.94128;
    g[2][0]=0.36856;
    g[3][0]=0.33650;
    g[4][0]=0.18171;
    //K+K- channel
    g[0][1]=-0.55377;
    g[1][1]=0.55095;
    g[2][1]=0.23888;
    g[3][1]=0.40907;
    g[4][1]=-0.17558;
    //4pi channel
    g[0][2]=0;
    g[1][2]=0;
    g[2][2]=0.55639;
    g[3][2]=0.85679;
    g[4][2]=-0.79658;
    //eta eta channel
    g[0][3]=-0.39899;
    g[1][3]=0.39065;
    g[2][3]=0.18340;
    g[3][3]=0.19906;
    g[4][3]=-0.00355;
    //eta eta' channel
    g[0][4]=-0.34639;
    g[1][4]=0.31503;
    g[2][4]=0.18681;
    g[3][4]=-0.00984;
    g[4][4]=0.22358;

    // Pole masses
    ma[0]=0.651;      
    ma[1]=1.20360;
    ma[2]=1.55817;
    ma[3]=1.21000;
    ma[4]=1.82206;

  } else if (solution == 1) { // solnI.txt 
    
    // coupling constants
    //pi+pi- channel
    g[0][0]=0.31896;
    g[1][0]=0.85963;
    g[2][0]=0.47993;
    g[3][0]=0.45121;
    g[4][0]=0.39391;
    //K+K- channel
    g[0][1]=-0.49998;
    g[1][1]=0.52402;
    g[2][1]=0.40254;
    g[3][1]=0.42769;
    g[4][1]=-0.30860;
    //4pi channel
    g[0][2]=0;
    g[1][2]=0;
    g[2][2]=1.0;
    g[3][2]=1.15088;
    g[4][2]=0.33999;
    //eta eta channel
    g[0][3]=-0.21554;
    g[1][3]=0.38093;
    g[2][3]=0.21811;
    g[3][3]=0.22925;
    g[4][3]=0.06919;
    //eta eta' channel
    g[0][4]=-0.18294;
    g[1][4]=0.23788;
    g[2][4]=0.05454;
    g[3][4]=0.06444;
    g[4][4]=0.32620;

    // Pole masses
    ma[0]=0.7369;
    ma[1]=1.24347;
    ma[2]=1.62681;
    ma[3]=1.21900;
    ma[4]=1.74932;

  } else if (solution == 2) { // solnIIa.txt 
    
    // coupling constants
    //pi+pi- channel
    g[0][0]=0.26014;
    g[1][0]=0.95289;
    g[2][0]=0.46244;
    g[3][0]=0.41848;
    g[4][0]=0.01804;
    //K+K- channel
    g[0][1]=-0.57849;
    g[1][1]=0.55887;
    g[2][1]=0.31712;
    g[3][1]=0.49910;
    g[4][1]=-0.28430;
    //4pi channel
    g[0][2]=0;
    g[1][2]=0;
    g[2][2]=0.70340;
    g[3][2]=0.96819;
    g[4][2]=-0.90100;
    //eta eta channel
    g[0][3]=-0.32936;
    g[1][3]=0.39910;
    g[2][3]=0.22963;
    g[3][3]=0.24415;
    g[4][3]=-0.07252;
    //eta eta' channel
    g[0][4]=-0.30906;
    g[1][4]=0.31143;
    g[2][4]=0.19802;
    g[3][4]=-0.00522;
    g[4][4]=0.17097;

    // Pole masses
    ma[0]=0.67460;
    ma[1]=1.21094;
    ma[2]=1.57896;
    ma[3]=1.21900;
    ma[4]=1.86602;
  } 

  //Now define the K-matrix pole
  double  rho1sq,rho2sq,rho4sq,rho5sq;
  EvtComplex rho[5];
  double f[5][5];

  //Initalize the mass of the resonance
  double mpi=0.13957;
  double mK=0.493677;     //using charged K value
  double meta=0.54775;    //using PDG value
  double metap=0.95778;   //using PDG value
    
  //Initialize the matrix to value zero
  EvtComplex K[5][5];
  for(int i=0;i<5;i++) { 
    for(int j=0;j<5;j++) {
      K[i][j]=EvtComplex(0,0);
      f[i][j]=0;
    }
  }

  //Input the _f[i][j] scattering data
  double s_scatt=0.0 ; 
  if (solution == 3) 
    s_scatt=-3.92637; 
  else if (solution == 1) 
    s_scatt= -5.0 ;
  else if (solution == 2) 
    s_scatt= -5.0 ; 
  double sa=1.0;
  double sa_0=-0.15;
  if (solution == 3) {
    f[0][0]=0.23399;  // f^scatt
    f[0][1]=0.15044;
    f[0][2]=-0.20545;
    f[0][3]=0.32825;
    f[0][4]=0.35412;
  }else if (solution == 1) {
    f[0][0]=0.04214;  // f^scatt
    f[0][1]=0.19865;
    f[0][2]=-0.63764;
    f[0][3]=0.44063;
    f[0][4]=0.36717;
  }else if (solution == 2) {
    f[0][0]=0.26447;  // f^scatt
    f[0][1]=0.10400;
    f[0][2]=-0.35445;
    f[0][3]=0.31596;
    f[0][4]=0.42483;
  }   
  f[1][0]=f[0][1];
  f[2][0]=f[0][2];
  f[3][0]=f[0][3];
  f[4][0]=f[0][4];

  //Now construct the phase-space factor
  //For eta-eta' there is no difference term
  rho1sq = 1. - pow( mpi + mpi, 2 ) / s;   //pi+ pi- phase factor
  if( rho1sq >= 0 )
    rho[ 0 ] = EvtComplex( sqrt( rho1sq ), 0 );
  else
    rho[ 0 ] = EvtComplex( 0, sqrt( -rho1sq ) );  

  rho2sq = 1. - pow( mK + mK, 2 ) / s;
  if( rho2sq >= 0 )
    rho[ 1 ] = EvtComplex( sqrt( rho2sq ), 0 );
  else
    rho[ 1 ] = EvtComplex( 0, sqrt( -rho2sq ) );

  //using the A&S 4pi phase space Factor:
  //Shit, not continue
  if( s <= 1 )
    {
      double real   = 1.2274 + .00370909 / ( s * s ) - .111203 / s - 6.39017 * s + 16.8358*s*s - 21.8845*s*s*s + 11.3153*s*s*s*s;
      double cont32 = sqrt(1.0-(16.0*mpi*mpi));
      rho[ 2 ] = EvtComplex( cont32 * real, 0 );
    }
  else
    rho[ 2 ] = EvtComplex( sqrt( 1. - 16. * mpi * mpi / s ), 0 );

  rho4sq = 1. - pow( meta + meta, 2 ) / s;
  if( rho4sq >= 0 )
    rho[ 3 ] = EvtComplex( sqrt( rho4sq ), 0 );
  else
    rho[ 3 ] = EvtComplex( 0, sqrt( -rho4sq ) );

  rho5sq = 1. - pow( meta + metap, 2 ) / s;
  if( rho5sq >= 0 )
    rho[ 4 ] = EvtComplex( sqrt( rho5sq ), 0 );
  else
    rho[ 4 ] = EvtComplex( 0, sqrt( -rho5sq ) );

  double smallTerm = 1; // Factor to prevent divergences.

  // Check if some pole may arise problems.
  for ( int pole = 0; pole < 5; pole++ )
    if ( fabs( pow( ma[ pole ], 2 ) - s ) < PRECISION )
      smallTerm = pow( ma[ pole ], 2 ) - s;

  //now sum all the pole
  //equation (3) in the E791 K-matrix paper
  for(int i=0;i<5;i++) { 
    for(int j=0;j<5;j++) {  
      for (int pole_index=0;pole_index<5;pole_index++) {
	double A=g[pole_index][i]*g[pole_index][j];
	double B=ma[pole_index]*ma[pole_index]-s;

	if ( fabs( B ) < PRECISION )
	  K[ i ][ j ] += EvtComplex( A    , 0 );
	else
	  K[ i ][ j ] += EvtComplex( A / B, 0 ) * smallTerm;
      }
    }
  }

  //now add the SVT part
  for(int i=0;i<5;i++) { 
    for(int j=0;j<5;j++) {
      double C=f[i][j]*(1.0-s_scatt);
      double D=(s-s_scatt);
      K[ i ][ j ] += EvtComplex( C / D, 0 ) * smallTerm;
    }
  }

  //Fix the bug in the FOCUS paper
  //Include the Alder zero term:
  for(int i=0;i<5;i++) { 
    for(int j=0;j<5;j++) {
      double E=(s-(sa*mpi*mpi*0.5))*(1.0-sa_0);
      double F=(s-sa_0);    
      K[ i ][ j ] *= EvtComplex(E/F,0);
    }
  }

  //This is not correct!
  //(1-ipK) != (1-iKp)
  static EvtMatrix< EvtComplex > mat;
  mat.setRange( 5 ); // Try to do in only the first time. DEFINE ALLOCATION IN CONSTRUCTOR.

  for ( int row = 0; row < 5; row++ )
    for ( int col = 0; col < 5; col++ )
      mat( row, col ) = ( row == col ) * smallTerm - EvtComplex( 0., 1. ) * K[ row ][ col ] * rho[ col ];


  EvtMatrix< EvtComplex >* matInverse = mat.inverse();  //The 1st row of the inverse matrix. This matrix is {(I-iKp)^-1}_0j
  vector< EvtComplex > U1j;
  for ( int j = 0; j < 5; j++ )
    U1j.push_back( (*matInverse)[ 0 ][ j ] );

  delete matInverse;

  //this calculates final F0 factor
  EvtComplex value( 0, 0 );
  if (index<=5) {
    //this calculates the beta_idx Factors
    for(int j=0;j<5;j++) {        // sum for 5 channel
      EvtComplex top    = U1j[j]*g[index-1][j];
      double     bottom = ma[index-1]*ma[index-1]-s;

      if ( fabs( bottom ) < PRECISION )
	value += top;
      else
	value += top / bottom * smallTerm;
    }
  } else {
    //this calculates fprod Factors
    value += U1j[0];
    value += U1j[1]*_fr12prod;
    value += U1j[2]*_fr13prod;
    value += U1j[3]*_fr14prod;
    value += U1j[4]*_fr15prod;

    value *= (1-_s0prod)/(s-_s0prod) * smallTerm;
  }

  return value;
}


//replace Breit-Wigner with LASS
EvtComplex EvtDalitzReso::lass(double s)
{
  EvtTwoBodyKine vd(_massFirst,_massSecond, sqrt(s));
  double q = vd.p();
  double GammaM = _g0*_vd.widthFactor(vd);  // running width;

  //calculate the background phase motion
  double cot_deltaB = 1.0/(_a*q) + 0.5*_r*q;
  double deltaB = atan( 1.0/cot_deltaB);
  double totalB = deltaB + _phiB ;

  //calculate the resonant phase motion
  double deltaR = atan((_m0*GammaM/(_m0*_m0 - s)));
  double totalR = deltaR + _phiR ;

  //sum them up
  EvtComplex  bkgB,resT;
  bkgB = EvtComplex(_Blass*sin(totalB),0)*EvtComplex(cos(totalB),sin(totalB));
  resT = EvtComplex(_R*sin(deltaR),0)*EvtComplex(cos(totalR),sin(totalR))*EvtComplex(cos(2*totalB),sin(2*totalB));

  EvtComplex T;
  if(_cutoff>0 && sqrt(s)>_cutoff) T = resT;
  else T = bkgB + resT;

  if(_scaleByMOverQ) T*=(sqrt(s)/q);

  return T;
}


EvtComplex EvtDalitzReso::flatte(const double& m) {

  EvtComplex w;

  for (vector<EvtFlatteParam>::const_iterator param = _flatteParams.begin();
       param != _flatteParams.end();
       ++param) {
    double m1 = (*param).m1(); double m2 = (*param).m2();
    double g = (*param).g();
    w += (g*g*sqrtCplx((1-((m1-m2)*(m1-m2))/(m*m))*(1-((m1+m2)*(m1+m2))/(m*m))));
  }
  
  EvtComplex denom = _m0*_m0 - m*m - EvtComplex(0,1)*w;

  return EvtComplex(1.0,0.0)/denom;
}
