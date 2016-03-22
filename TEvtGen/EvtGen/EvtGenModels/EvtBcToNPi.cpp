//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package. If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//
// Module: EvtGenModels/EvtBcToNPi.hh
//
// Description: General decay model for Bc -> V + npi and Bc -> P + npi
//
// Modification history:
//
//    A.Berezhnoy, A.Likhoded, A.Luchinsky  April 2011   Module created
//
//------------------------------------------------------------------------

#include "EvtGenBase/EvtPatches.hh"

#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtSpinType.hh"
#include "EvtGenBase/EvtScalarParticle.hh"

#include "EvtGenModels/EvtBcToNPi.hh"

#include <iostream>
using std::endl;

EvtBcToNPi::EvtBcToNPi(bool printAuthorInfo) {
  nCall=0; maxAmp2=0;
  if (printAuthorInfo == true) {this->printAuthorInfo();}
}

EvtBcToNPi::~EvtBcToNPi() { 
}

std::string EvtBcToNPi::getName(){

  return "EvtBcToNPi";     

}

EvtDecayBase* EvtBcToNPi::clone(){

  return new EvtBcToNPi;

}


void EvtBcToNPi::init(){

	// check spins
	checkSpinParent(EvtSpinType::SCALAR);


	// the others are scalar
	for (int i=1; i<=(getNDaug()-1);i++) {
		checkSpinDaughter(i,EvtSpinType::SCALAR);
	};

	_beta=-0.108; _mRho=0.775; _gammaRho=0.149;
	_mRhopr=1.364; _gammaRhopr=0.400; _mA1=1.23; _gammaA1=0.4;

	// read arguments
	if( EvtPDL::getSpinType(getDaug(0)) == EvtSpinType::VECTOR) {
		checkNArg(10);
		int n=0;
		_maxProb=getArg(n++);
		FA0_N=getArg(n++);
		FA0_c1=getArg(n++);
		FA0_c2=getArg(n++);
		FAp_N=getArg(n++);
		FAp_c1=getArg(n++);
		FAp_c2=getArg(n++);
		FV_N=getArg(n++);
		FV_c1=getArg(n++);
		FV_c2=getArg(n++);
		FAm_N=0;
		FAm_c1=0;
		FAm_c2=0;
	}
	else if( EvtPDL::getSpinType(getDaug(0)) == EvtSpinType::SCALAR) {
		checkNArg(4);
		int n=0;
		_maxProb=getArg(n++);
		Fp_N=getArg(n++);
		Fp_c1=getArg(n++);
		Fp_c2=getArg(n++);
		Fm_N=0;
		Fm_c1=0;
		Fm_c2=0;
	}
	else {
		report(Severity::Error,"EvtGen") << "Have not yet implemented this final state in BCPSINPI model" << endl;
		report(Severity::Error,"EvtGen") << "Ndaug="<<getNDaug() << endl;
		for ( int id=0; id<(getNDaug()-1); id++ ) 
			report(Severity::Error,"EvtGen") << "Daug " << id << " "<<EvtPDL::name(getDaug(id)).c_str() << endl;
		return;

	};

	if(getNDaug()<2 || getNDaug()>4) {
		report(Severity::Error,"EvtGen") << "Have not yet implemented this final state in BCPSINPI model" << endl;
		report(Severity::Error,"EvtGen") << "Ndaug="<<getNDaug() << endl;
		for ( int id=0; id<(getNDaug()-1); id++ ) 
			report(Severity::Error,"EvtGen") << "Daug " << id << " "<<EvtPDL::name(getDaug(id)).c_str() << endl;
		return;
	}

}

double EvtBcToNPi::_ee(double M, double m1, double m2) {
        return (M*M+m1*m1-m2*m2)/(2*M);
};

double EvtBcToNPi::_pp(double M, double m1, double m2) {
        double __ee=_ee(M,m1,m2);
        return sqrt(__ee*__ee-m1*m1);
};


void EvtBcToNPi::initProbMax(){
	if(_maxProb>0.)	setProbMax(_maxProb);
	else {
                EvtId id=getParentId();
                EvtScalarParticle *p=new EvtScalarParticle();
                p->init(id, EvtPDL::getMass(id),0., 0., 0.);
                p->setDiagonalSpinDensity();
                // add daughters
                p->makeDaughters(getNDaug(), getDaugs() );

                // fill the momenta
                if(getNDaug()==2) {
                        double M=EvtPDL::getMass(id), m1=EvtPDL::getMass(getDaug(0)), m2=EvtPDL::getMass(getDaug(1));
                        double __pp=_pp(M,m1,m2);
                        p->getDaug(0)->setP4( EvtVector4R( _ee(M,m1,m2), 0., 0., __pp) );
                        p->getDaug(1)->setP4( EvtVector4R( _ee(M,m2,m1), 0., 0., -__pp) );
                }
                else if( getNDaug()==3) {
                        double M=EvtPDL::getMass(id), 
                                m1=EvtPDL::getMass(getDaug(0)), 
                                m2=EvtPDL::getMass(getDaug(1)),
                                m3=EvtPDL::getMass(getDaug(2));
                        double __ppRho=_pp(M,m1,_mRho), __ppPi=_pp(_mRho,m2,m3);
                        p->getDaug(0)->setP4( EvtVector4R( _ee(M,m1,_mRho), 0., 0., __ppRho) );
                        EvtVector4R _pRho( _ee(M,_mRho,m1), 0., 0., -__ppRho);
                        EvtVector4R _p2( _ee(_mRho, m2, m3), 0., 0., __ppPi); _p2.applyBoostTo(_pRho);
                        EvtVector4R _p3( _ee(_mRho, m2, m3), 0., 0., -__ppPi); _p3.applyBoostTo(_pRho);
                        p->getDaug(1)->setP4(_p2);
                        p->getDaug(2)->setP4(_p3);

                }
                else if( getNDaug()==4) {
                        double M=EvtPDL::getMass(id), 
                                m1=EvtPDL::getMass(getDaug(0)), 
                                m2=EvtPDL::getMass(getDaug(1)),
                                m3=EvtPDL::getMass(getDaug(2)),
                                m4=EvtPDL::getMass(getDaug(3));
                        if(M<m1+_mA1) return;
                        double   __ppA1=_pp(M,m1,_mA1),
                                 __ppRho=_pp(_mA1,_mRho,m4),
                                 __ppPi=_pp(_mRho, m2, m3);
                        p->getDaug(0)->setP4( EvtVector4R( _ee(M,m1,_mRho), 0., 0., __ppA1) );
                        EvtVector4R _pA1( _ee(M,_mA1,m1), 0., 0., -__ppA1);
                        EvtVector4R _pRho(_ee(_mA1, _mRho, m4), 0, 0, __ppRho);
                        _pRho.applyBoostTo(_pA1);
                        EvtVector4R _p4( _ee(_mA1, m4, _mRho), 0, 0, -__ppRho); _p4.applyBoostTo(_pA1);
                        p->getDaug(3)->setP4(_p4);
                        EvtVector4R _p2( _ee(_mRho, m2, m3), 0, 0, __ppPi); _p2.applyBoostTo(_pRho);
                        p->getDaug(1)->setP4(_p2);
                        EvtVector4R _p3( _ee(_mRho, m2, m3), 0, 0, -__ppPi); _p2.applyBoostTo(_pRho);
                        p->getDaug(2)->setP4(_p3);
                };


                _amp2.init(p->getId(),getNDaug(),getDaugs());

		decay(p); 	

                EvtSpinDensity rho=_amp2.getSpinDensity();

                double prob=p->getSpinDensityForward().normalizedProb(rho);

                if(prob>0) setProbMax(0.9*prob);


	};
}

void EvtBcToNPi::decay( EvtParticle *root_particle ){
	++nCall;

	EvtIdSet thePis("pi+","pi-","pi0");
	EvtComplex I=EvtComplex(0.0, 1.0);


  	root_particle->initializePhaseSpace(getNDaug(),getDaugs());

	EvtVector4R
		p(root_particle->mass(), 0., 0., 0.),                  // Bc momentum
		k=root_particle->getDaug(0)->getP4(),     		   // J/psi momenta
		Q=p-k;
	
	double Q2=Q.mass2();


// check pi-mesons and calculate hadronic current
	EvtVector4C hardCur;
	bool foundHadCurr=false;

	if ( getNDaug() == 2 ) // Bc -> psi pi+
	{
		hardCur=Q;
		foundHadCurr=true;
	}
	else if ( getNDaug() == 3 ) // Bc -> psi pi+ pi0
	{
		EvtVector4R p1,p2;
		p1=root_particle->getDaug(1)->getP4(),     // pi+ momenta
		p2=root_particle->getDaug(2)->getP4(),    // pi0 momentum
		hardCur=Fpi(p1,p2)*(p1-p2);
		foundHadCurr=true;
	}
	else if( getNDaug()==4 ) // Bc -> psi pi+ pi pi
	{
		int diffPi(0),samePi1(0),samePi2(0);
		if ( getDaug( 1) == getDaug( 2) ) {diffPi= 3; samePi1= 1; samePi2= 2;}
		if ( getDaug( 1) == getDaug( 3) ) {diffPi= 2; samePi1= 1; samePi2= 3;}
		if ( getDaug( 2) == getDaug( 3) ) {diffPi= 1; samePi1= 2; samePi2= 3;}

		EvtVector4R p1=root_particle->getDaug(samePi1)->getP4();
		EvtVector4R p2=root_particle->getDaug(samePi2)->getP4();
		EvtVector4R p3=root_particle->getDaug(diffPi)->getP4();
		
		EvtComplex BA1;
		double GA1=_gammaA1*pi3G(Q2,samePi1)/pi3G(_mA1*_mA1,samePi1);
		EvtComplex denBA1(_mA1*_mA1 - Q.mass2(),-1.*_mA1*GA1);
		BA1 = _mA1*_mA1 / denBA1;

		hardCur = BA1*( (p1-p3) - (Q*(Q*(p1-p3))/Q2)*Fpi(p2,p3) +
				(p2-p3) - (Q*(Q*(p2-p3))/Q2)*Fpi(p1,p3) ); 
		foundHadCurr=true;
	}

	if ( !foundHadCurr ) {
		report(Severity::Error,"EvtGen") << "Have not yet implemented this final state in BCNPI model" << endl;
		report(Severity::Error,"EvtGen") << "Ndaug="<<getNDaug() << endl;
		int id;
		for ( id=0; id<(getNDaug()-1); id++ ) 
		report(Severity::Error,"EvtGen") << "Daug " << id << " "<<EvtPDL::name(getDaug(id)).c_str() << endl;
		::abort();
	};

 	EvtTensor4C H;
	double amp2=0.;
	if( root_particle->getDaug(0)->getSpinType() == EvtSpinType::VECTOR) {
		double FA0=FA0_N*exp(FA0_c1*Q2 + FA0_c2*Q2*Q2);
		double FAp=FAp_N*exp(FAp_c1*Q2 + FAp_c2*Q2*Q2);
		double FAm=FAm_N*exp(FAm_c1*Q2 + FAm_c2*Q2*Q2);
		double FV= FV_N* exp( FV_c1*Q2 + FV_c2*Q2*Q2 );
		H=-FA0*EvtTensor4C::g()
			-FAp*EvtGenFunctions::directProd(p,p+k)
			+FAm*EvtGenFunctions::directProd(p,p-k)
			+2*I*FV*dual(EvtGenFunctions::directProd(p,k));
		EvtVector4C  Heps=H.cont2(hardCur);

		for(int i=0; i<4; i++) {
			EvtVector4C  eps=root_particle->getDaug(0)->epsParent(i).conj(); // psi-meson polarization vector
			EvtComplex amp=eps*Heps;
			vertex(i,amp);
			amp2+=pow( abs(amp),2);
		}
	}
	else if( root_particle->getDaug(0)->getSpinType() == EvtSpinType::SCALAR) {
		double Fp=Fp_N*exp(Fp_c1*Q2 + Fp_c2*Q2*Q2);
		double Fm=Fm_N*exp(Fm_c1*Q2 + Fm_c2*Q2*Q2);
		EvtVector4C H=Fp*(p+k)+Fm*(p-k);
		EvtComplex amp=H*hardCur;
		vertex(amp);
		amp2+=pow( abs(amp),2);
	};
	if(amp2>maxAmp2) maxAmp2=amp2;

  return ;
}

EvtComplex EvtBcToNPi::Fpi( EvtVector4R q1, EvtVector4R q2) {
	double m1=q1.mass();
	double m2=q2.mass();
	
	EvtVector4R Q = q1 + q2;
	double mQ2= Q*Q;
	
	// momenta in the rho->pipi decay
	double dRho= _mRho*_mRho - m1*m1 - m2*m2;
	double pPiRho = (1.0/_mRho)*sqrt((dRho*dRho)/4.0 - m1*m1*m2*m2);
	
	double dRhopr= _mRhopr*_mRhopr - m1*m1 - m2*m2;
	double pPiRhopr = (1.0/_mRhopr)*sqrt((dRhopr*dRhopr)/4.0 - m1*m1*m2*m2);
	
	double dQ= mQ2 - m1*m1 - m2*m2;
	double pPiQ = (1.0/sqrt(mQ2))*sqrt((dQ*dQ)/4.0 - m1*m1*m2*m2);
	
	
	double gammaRho = _gammaRho*_mRho/sqrt(mQ2)*pow((pPiQ/pPiRho),3);
	EvtComplex BRhoDem(_mRho*_mRho - mQ2,-1.0*_mRho*gammaRho);
	EvtComplex BRho= _mRho*_mRho / BRhoDem;
	
	double gammaRhopr = _gammaRhopr*_mRhopr/sqrt(mQ2)*pow((pPiQ/pPiRhopr),3);
	EvtComplex BRhoprDem(_mRhopr*_mRhopr - mQ2,-1.0*_mRho*gammaRhopr);
	EvtComplex BRhopr= _mRhopr*_mRhopr / BRhoprDem;
	
	return (BRho + _beta*BRhopr)/(1+_beta);
}

double EvtBcToNPi::pi3G(double m2,int dupD) {
	double mPi= EvtPDL::getMeanMass(getDaug(dupD));
	if ( m2 > (_mRho+mPi) ) {
	return m2*(1.623 + 10.38/m2 - 9.32/(m2*m2) + 0.65/(m2*m2*m2));
	}
	else {
	double t1=m2-9.0*mPi*mPi;
	return 4.1*pow(t1,3.0)*(1.0 - 3.3*t1+5.8*t1*t1);
	}
}

void EvtBcToNPi::printAuthorInfo() {
  
  report(Severity::Info,"EvtGen")<<"Defining EvtBcToNPi model: Bc -> V + npi and Bc -> P + npi decays\n"
		       <<"from A.V. Berezhnoy, A.K. Likhoded, A.V. Luchinsky: "
		       <<"Phys.Rev.D 82, 014012 (2010) and arXiV:1104.0808."<<endl;

}

