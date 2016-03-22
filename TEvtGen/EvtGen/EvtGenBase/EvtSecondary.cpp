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
// Module: EvtSecondary.cc
//
// Description: Class to store the decays of the secondary particles.
//
// Modification history:
//
//    RYD       March 12, 1998       Module created
//
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPatches.hh"
#include <iostream>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtSecondary.hh"
#include "EvtGenBase/EvtReport.hh"
using std::endl;
using std::ostream;


void EvtSecondary::init(){
  _npart=0;
}
  
int EvtSecondary::getNPart(){
  return _npart;
}

void EvtSecondary::createSecondary(int stdhepindex,EvtParticle* prnt){

  _stdhepindex[_npart]=stdhepindex;
  if (prnt->getNDaug()==0){
    _id1[_npart]=0;
    _id2[_npart]=0;
    _id3[_npart]=0;
    _npart++;
    return;
  }
  if (prnt->getNDaug()==1){
    _id1[_npart]=EvtPDL::getStdHep(prnt->getDaug(0)->getId());
    _id2[_npart]=0;
    _id3[_npart]=0;
    _npart++;
    return;
  }
  if (prnt->getNDaug()==2){
    _id1[_npart]=EvtPDL::getStdHep(prnt->getDaug(0)->getId());
    _id2[_npart]=EvtPDL::getStdHep(prnt->getDaug(1)->getId());
    _id3[_npart]=0;
    _npart++;
    return;
  }
  if (prnt->getNDaug()==3){
    _id1[_npart]=EvtPDL::getStdHep(prnt->getDaug(0)->getId());
    _id2[_npart]=EvtPDL::getStdHep(prnt->getDaug(1)->getId());
    _id3[_npart]=EvtPDL::getStdHep(prnt->getDaug(2)->getId());
    _npart++;
    return;
  }
  
  report(Severity::Error,"EvtGen") << 
    "More than 3 decay products in a secondary particle!"<<endl;


}
 

ostream& operator<<(ostream& s, const EvtSecondary& secondary){

  s <<endl;
  s << "Secondary decays:"<<endl;

  int i;
  for(i=0;i<secondary._npart;i++){

    report(Severity::Info,"EvtGen") <<i<<" "
	 <<secondary._stdhepindex[i]<<" "
	 <<secondary._id1[i]<<" "
	 <<secondary._id2[i]<<" "
	 <<secondary._id3[i]<<endl;

  }
  
  s<<endl;
  
  return s;

}  

