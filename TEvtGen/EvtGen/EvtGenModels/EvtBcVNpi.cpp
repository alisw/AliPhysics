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
// Module: EvtBcVNpi.cc
//
// Description: Module to implement Bc -> psi + (n pi) decays
//
// Modification history:
//
//    AVL     July 6, 2012        Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenModels/EvtTauHadnu.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtParser.hh"

#include "EvtGenModels/EvtBcVNpi.hh"
#include "EvtGenModels/EvtWnPi.hh"



EvtBcVNpi::~EvtBcVNpi() {
//   cout<<"BcVNpi::destructor : nCall = "<<nCall<<" getProbMax(-1) = "<<getProbMax(-1)<<endl;
  
}

std::string EvtBcVNpi::getName(){ return "BC_VNPI";}

EvtDecayBase* EvtBcVNpi::clone() {  return new EvtBcVNpi;}

//======================================================
void EvtBcVNpi::init(){
    //cout<<"BcVNpi::init()"<<endl;
    
    checkNArg(1);
    checkSpinParent(EvtSpinType::SCALAR);
    checkSpinDaughter(0,EvtSpinType::VECTOR);
    for (int i=1; i<=(getNDaug()-1);i++) {
      checkSpinDaughter(i,EvtSpinType::SCALAR);
    };

    if(getNDaug()<2 || getNDaug()>6) {
      report(Severity::Error,"EvtGen") << "Have not yet implemented this final state in BcVNpi model" << endl;
      report(Severity::Error,"EvtGen") << "Ndaug="<<getNDaug() << endl;
      for ( int id=0; id<(getNDaug()-1); id++ ) 
	report(Severity::Error,"EvtGen") << "Daug " << id << " "<<EvtPDL::name(getDaug(id)).c_str() << endl;
      return;
    }

  
//     for(int i=0; i<getNDaug(); i++)
//       cout<<"BcVNpi::init \t\t daughter "<<i<<" : "<<getDaug(i).getId()<<"   "<<EvtPDL::name(getDaug(i)).c_str()<<endl;

   idVector = getDaug(0).getId();
    whichfit = int(getArg(0)+0.1);
//     cout<<"BcVNpi: whichfit ="<<whichfit<<"  idVector="<<idVector<<endl;
    ffmodel = new EvtBCVFF(idVector,whichfit);
    
    wcurr = new EvtWnPi();
    
    nCall = 0;
}

//======================================================
void EvtBcVNpi::initProbMax() {
//     cout<<"BcVNpi::initProbMax()"<<endl;
    if(idVector == EvtPDL::getId("J/psi").getId() && whichfit == 1 && getNDaug()==6) setProbMax(720000.);
    else if(idVector == EvtPDL::getId("J/psi").getId() && whichfit == 2 && getNDaug()==6) setProbMax(471817.);
    else if(idVector == EvtPDL::getId("J/psi").getId() && whichfit == 1 && getNDaug()==4) setProbMax(42000.);
    else if(idVector == EvtPDL::getId("J/psi").getId() && whichfit == 2 && getNDaug()==4) setProbMax(16000.);
    
    else if(idVector == EvtPDL::getId("psi(2S)").getId() && whichfit == 1 && getNDaug()==4) setProbMax(1200.);
    else if(idVector == EvtPDL::getId("psi(2S)").getId() && whichfit == 2 && getNDaug()==4) setProbMax(2600.);
    else if(idVector == EvtPDL::getId("psi(2S)").getId() && whichfit == 1 && getNDaug()==6) setProbMax(40000.);
    else if(idVector == EvtPDL::getId("psi(2S)").getId() && whichfit == 2 && getNDaug()==6) setProbMax(30000.);
}

//======================================================
void EvtBcVNpi::decay( EvtParticle *root_particle ) {
   ++nCall;
//     cout<<"BcVNpi::decay()"<<endl;
    root_particle->initializePhaseSpace(getNDaug(),getDaugs());

    EvtVector4R
	    p4b(root_particle->mass(), 0., 0., 0.),                  // Bc momentum
	    p4meson=root_particle->getDaug(0)->getP4(),     		   // J/psi momenta
	    Q=p4b-p4meson;
    double Q2=Q.mass2();


// check pi-mesons and calculate hadronic current
    EvtVector4C hardCur;
//     bool foundHadCurr=false;
    if( getNDaug() == 2) {
      hardCur = wcurr->WCurrent( root_particle->getDaug(1)->getP4() );
//       foundHadCurr=true;
    }
    else if( getNDaug() == 3) {
      hardCur = wcurr->WCurrent( root_particle->getDaug(1)->getP4() , 
				 root_particle->getDaug(2)->getP4() 
			       );
//       foundHadCurr=true;   
    }
    else if( getNDaug() == 4) {
      hardCur = wcurr->WCurrent( root_particle->getDaug(1)->getP4() , 
				 root_particle->getDaug(2)->getP4(), 
				 root_particle->getDaug(3)->getP4() 
			       );
//       foundHadCurr=true;         
    }
    else if( getNDaug() == 6) // Bc -> psi pi+ pi+ pi- pi- pi+ from [Kuhn, Was, hep-ph/0602162
    {

		hardCur = wcurr->WCurrent(root_particle->getDaug(1)->getP4(),
					  root_particle->getDaug(2)->getP4(),
					  root_particle->getDaug(3)->getP4(),
					  root_particle->getDaug(4)->getP4(),
					  root_particle->getDaug(5)->getP4()
				 );
// 		foundHadCurr=true;
    }	
    else {
	    report(Severity::Error,"EvtGen") << "Have not yet implemented this final state in BCNPI model" << endl;
	    report(Severity::Error,"EvtGen") << "Ndaug="<<getNDaug() << endl;
	    int id;
	    for ( id=0; id<(getNDaug()-1); id++ ) 
	    report(Severity::Error,"EvtGen") << "Daug " << id << " "<<EvtPDL::name(getDaug(id)).c_str() << endl;
	    ::abort();
    };  

// calculate Bc -> V W form-factors
	double a1f, a2f, vf, a0f;
	double m_meson = root_particle->getDaug(0)->mass();
	double m_b = root_particle->mass();
	ffmodel->getvectorff(root_particle->getId(),
				root_particle->getDaug(0)->getId(),
				Q2,
				m_meson,
				&a1f, 
				&a2f, 
				&vf, 
				&a0f);
	double a3f = ((m_b+m_meson)/(2.0*m_meson))*a1f -
	      ((m_b-m_meson)/(2.0*m_meson))*a2f;

// calculate Bc -> V W current
	EvtTensor4C H;
	H = a1f*(m_b+m_meson)*EvtTensor4C::g();
	H.addDirProd((-a2f/(m_b+m_meson))*p4b,p4b+p4meson);
	H+=EvtComplex(0.0,vf/(m_b+m_meson))*dual(EvtGenFunctions::directProd(p4meson+p4b,p4b-p4meson));
	H.addDirProd((a0f-a3f)*2.0*(m_meson/Q2)*p4b,p4b-p4meson);
	EvtVector4C  Heps=H.cont2(hardCur);
	
	for(int i=0; i<4; i++) {
		EvtVector4C  eps=root_particle->getDaug(0)->epsParent(i).conj(); // psi-meson polarization vector
		EvtComplex amp=eps*Heps;
		vertex(i,amp);
	};

}

