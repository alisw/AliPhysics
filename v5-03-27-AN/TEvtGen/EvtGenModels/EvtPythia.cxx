//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See BelEvtGen/COPYRIGHT
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtJetSet.cc
//
// Description: Routine to use JetSet for decaying particles.
//
// Modification history:
//
//    RYD     July 24, 1997        Module created
//    RS      October 28, 2002        copied from JETSET module
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtStringParticle.hh"
#include "EvtGenBase/EvtDecayTable.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenModels/EvtPythia.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtId.hh"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
using std::endl;
using std::fstream;
using std::ios;
using std::ofstream;
using std::resetiosflags;
using std::setiosflags;
using std::setw;

using std::string;

int EvtPythia::njetsetdecays=0;
  EvtDecayBasePtr* EvtPythia::jetsetdecays=0; 
int EvtPythia::ntable=0;

int EvtPythia::ncommand=0;
int EvtPythia::lcommand=0;
std::string* EvtPythia::commands=0;


extern "C" {
  extern void pycontinuum_(double *,int *, int *, 
                            double *,double *,double *,double *);
}

extern "C" {
  extern void evtpythiainit_(const char* fname, int len);
}

extern "C" {
  extern void init_cont_();
}

extern "C" {
  extern void pythiadec_(int *,double *,int *,int *,int *,
			  double *,double *,double *,double *);
}

extern "C" {
  extern void initpythia_(int *);
}

extern "C" {
  extern void pygive_(const char *cnfgstr,int length);
}

extern "C" {
  extern int pycomp_(int* kf);
}

extern "C" {
  extern void pylist_(int &);
}


EvtPythia::EvtPythia(){}

EvtPythia::~EvtPythia(){


  int i;


  //the deletion of commands is really uggly!

  if (njetsetdecays==0) {
    delete [] commands;
    commands=0;
    return;
  }

  for(i=0;i<njetsetdecays;i++){
    if (jetsetdecays[i]==this){
      jetsetdecays[i]=jetsetdecays[njetsetdecays-1];
      njetsetdecays--;
      if (njetsetdecays==0) {
	delete [] commands;
	commands=0;
      }
      return;
    }
  }

  report(ERROR,"EvtGen") << "Error in destroying Pythia model!"<<endl;
 
}


std::string EvtPythia::getName(){

  return "PYTHIA";     

}

EvtDecayBase* EvtPythia::clone(){

  return new EvtPythia;

}


void EvtPythia::initProbMax(){

  noProbMax();

}


void EvtPythia::init(){

  checkNArg(1);


  if (getParentId().isAlias()){

    report(ERROR,"EvtGen") << "EvtPythia finds that you are decaying the"<<endl
                           << " aliased particle "
			   << EvtPDL::name(getParentId()).c_str()
			   << " with the Pythia model"<<endl
			   << " this does not work, please modify decay table."
			   << endl;
    report(ERROR,"EvtGen") << "Will terminate execution!"<<endl;
    ::abort();

  }

  store(this);

}


std::string EvtPythia::commandName(){

  return std::string("JetSetPar");
  
}


void EvtPythia::command(std::string cmd){

  if (ncommand==lcommand){

    lcommand=10+2*lcommand;

    std::string* newcommands=new std::string[lcommand];
    
    int i;

    for(i=0;i<ncommand;i++){
      newcommands[i]=commands[i];
    }
    
    delete [] commands;

    commands=newcommands;

  }

  commands[ncommand]=cmd;

  ncommand++;

}

void EvtPythia::pythiacont(double *energy, int *ndaugjs, int *kf,
			   double *px, double *py, double *pz, double *e)
{
  pycontinuum_(energy,ndaugjs,kf,px,py,pz,e);
}



void EvtPythia::decay( EvtParticle *p){
  
  
  //added by Lange Jan4,2000
  static EvtId STRNG=EvtPDL::getId("string");
  
  int istdheppar=EvtPDL::getStdHep(p->getId());
  
  if (pycomp_(&istdheppar)==0){
    report(ERROR,"EvtGen") << "Pythia can not decay:"
      <<EvtPDL::name(p->getId()).c_str()<<endl;
    return;
  }
  
  double mp=p->mass();
  
  EvtVector4R p4[20];
  
  int i,more;
  int ip=EvtPDL::getStdHep(p->getId());
  int ndaugjs;
  int kf[100];
  EvtId evtnumstable[100],evtnumparton[100];
  int stableindex[100],partonindex[100];
  int numstable;
  int numparton;
  int km[100];
  EvtId type[MAX_DAUG];
  
  pythiaInit(0);
  
  double px[100],py[100],pz[100],e[100];
  if ( p->getNDaug() != 0 ) { p->deleteDaughters(true);}
  
  int count=0;
  
  do{
    
    pythiadec_(&ip,&mp,&ndaugjs,kf,km,px,py,pz,e);
    
    
    numstable=0;
    numparton=0;
    
    for(i=0;i<ndaugjs;i++){
      
      if (EvtPDL::evtIdFromStdHep(kf[i])==EvtId(-1,-1)) {
	report(ERROR,"EvtGen") << "Pythia returned particle:"<<kf[i]<<endl;
	report(ERROR,"EvtGen") << "This can not be translated to evt number"<<endl;
	report(ERROR,"EvtGen") << "and the decay will be rejected!"<<endl;
	report(ERROR,"EvtGen") << "The decay was of particle:"<<ip<<endl;
	int i=1;
	pylist_(i);
      }

      //sort out the partons
      if (abs(kf[i])<=6||kf[i]==21){
	partonindex[numparton]=i;
	evtnumparton[numparton]=EvtPDL::evtIdFromStdHep(kf[i]);
	numparton++;
      }
      else{
	stableindex[numstable]=i;
	evtnumstable[numstable]=EvtPDL::evtIdFromStdHep(kf[i]); 
	numstable++;
      }
      
      
      // have to protect against negative mass^2 for massless particles
      // i.e. neutrinos and photons.
      // this is uggly but I need to fix it right now....
      
      if (px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]>=e[i]*e[i]){
	
        e[i]=sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i])+0.0000000000001;
	
      }
      
      p4[i].set(e[i],px[i],py[i],pz[i]);
      
      
    }
    
    int channel=EvtDecayTable::inChannelList(p->getId(),numstable,evtnumstable);
    
    
   more=(channel!=-1);
   
   
   
   
   count++;
   
  }while( more && (count<10000) );
  
  if (count>9999) {
    report(INFO,"EvtGen") << "Too many loops in EvtPythia!!!"<<endl;
    report(INFO,"EvtGen") << "Parent:"<<EvtPDL::name(getParentId()).c_str()<<endl;
    for(i=0;i<numstable;i++){
      report(INFO,"EvtGen") << "Daug("<<i<<")"<<EvtPDL::name(evtnumstable[i]).c_str()<<endl;
    }
    
  }
  
  
  
  if (numparton==0){
    
    p->makeDaughters(numstable,evtnumstable);
    
    for(i=0;i<numstable;i++){
      p->getDaug(i)->init(evtnumstable[i],p4[stableindex[i]]);
    }
    
    fixPolarizations(p);
    
    return ;
   
  }
  else{
    
    //have partons in JETSET
    
    EvtVector4R p4string(0.0,0.0,0.0,0.0);
    
    for(i=0;i<numparton;i++){
      p4string+=p4[partonindex[i]];
    }
    
    int nprimary=1;
    type[0]=STRNG;
    for(i=0;i<numstable;i++){
      if (km[stableindex[i]]==0){
	type[nprimary++]=evtnumstable[i];
      }
    }
    
    p->makeDaughters(nprimary,type);
    
    p->getDaug(0)->init(STRNG,p4string);
    
    EvtVector4R p4partons[10];
    
    for(i=0;i<numparton;i++){
      p4partons[i]=p4[partonindex[i]];
    }
    
    ((EvtStringParticle*)p->getDaug(0))->initPartons(numparton,p4partons,evtnumparton);
    
    
    
    nprimary=1;
    
    for(i=0;i<numstable;i++){
      
      if (km[stableindex[i]]==0){
	p->getDaug(nprimary++)->init(evtnumstable[i],p4[stableindex[i]]);
      }
    }
    
    
    int nsecond=0;
    for(i=0;i<numstable;i++){
      if (km[stableindex[i]]!=0){
	type[nsecond++]=evtnumstable[i];
      }
    }
    
    
    p->getDaug(0)->makeDaughters(nsecond,type);
    
    nsecond=0;
    for(i=0;i<numstable;i++){
      if (km[stableindex[i]]!=0){
	p4[stableindex[i]]=boostTo(p4[stableindex[i]],p4string);
	p->getDaug(0)->getDaug(nsecond)->init(evtnumstable[i],p4[stableindex[i]]);
	p->getDaug(0)->getDaug(nsecond)->setDiagonalSpinDensity();
	p->getDaug(0)->getDaug(nsecond)->decay();
	nsecond++;
      }
    }
    
    fixPolarizations(p);
    
    return ;
    
  }

}

void EvtPythia::fixPolarizations(EvtParticle *p){

  //special case for now to handle the J/psi polarization

  int ndaug=p->getNDaug();
  
  int i;

  static EvtId Jpsi=EvtPDL::getId("J/psi");

  for(i=0;i<ndaug;i++){
    if(p->getDaug(i)->getId()==Jpsi){
  
      EvtSpinDensity rho;
      
      rho.setDim(3);
      rho.set(0,0,0.5);
      rho.set(0,1,0.0);
      rho.set(0,2,0.0);

      rho.set(1,0,0.0);
      rho.set(1,1,1.0);
      rho.set(1,2,0.0);

      rho.set(2,0,0.0);
      rho.set(2,1,0.0);
      rho.set(2,2,0.5);

      EvtVector4R p4Psi=p->getDaug(i)->getP4();

      double alpha=atan2(p4Psi.get(2),p4Psi.get(1));
      double beta=acos(p4Psi.get(3)/p4Psi.d3mag());


      p->getDaug(i)->setSpinDensityForwardHelicityBasis(rho,alpha,beta,0.0);
      setDaughterSpinDensity(i);

    }
  }

}

void EvtPythia::store(EvtDecayBase* jsdecay){

  if (njetsetdecays==ntable){

    EvtDecayBasePtr* newjetsetdecays=new EvtDecayBasePtr[2*ntable+10];
    int i;
    for(i=0;i<ntable;i++){
      newjetsetdecays[i]=jetsetdecays[i];
    }
    ntable=2*ntable+10;
    delete [] jetsetdecays;
    jetsetdecays=newjetsetdecays;
  }

  jetsetdecays[njetsetdecays++]=jsdecay;



}


void EvtPythia::WritePythiaEntryHeader(ofstream &outdec, int lundkc,
			       EvtId evtnum,std::string name,
			       int chg, int cchg, int spin2,double mass,
			       double width, double maxwidth,double ctau,
			       int stable,double rawbrfrsum){

  char sname[100];
  char ccsname[100];
  
  // RS changed to 16, new PYTHIA standard
  int namelength=16;
  
  int i,j;
  int temp;
  temp = spin2;

  if (ctau>1000000.0) ctau=0.0;
  
  strcpy(sname,name.c_str());
  
  i=0;

  while (sname[i]!=0){
    i++;
  }

  // strip up to two + or -
 
  if(evtnum.getId()>=0) {
    if (sname[i-1]=='+'||sname[i-1]=='-'){ 
      sname[i-1]=0;
      i--;
    }
    if (sname[i-1]=='+'||sname[i-1]=='-'){ 
      sname[i-1]=0;
      i--;
    }
    // strip 0 except for _0 and chi...0
    if (sname[i-1]=='0' && sname[i-2]!='_' && !(sname[0]=='c' && sname[1]=='h')){
      sname[i-1]=0;
      i--;
    }
  }
  
  if (i>namelength) {
    for(j=1;j<namelength;j++){
      sname[j]=sname[j+i-namelength];
    }
    sname[namelength]=0;
  }
  
  // RS: copy name for cc particle
  for(j=0;j<=namelength;j++)
    ccsname[j]=sname[j];
  i=0;
  while (ccsname[i]!=' '){
    i++;
    if(ccsname[i]==0) break;
  }
  if(i<namelength)
    {
      ccsname[i]='b';
      ccsname[i+1]=0;
    }
  
  //cchg=0;
  
  if(evtnum.getId()>=0) {
    if (abs(EvtPDL::getStdHep(evtnum))==21) cchg=2;
    if (abs(EvtPDL::getStdHep(evtnum))==90) cchg=-1;
    if ((abs(EvtPDL::getStdHep(evtnum))<=8)&&
	(abs(EvtPDL::getStdHep(evtnum))!=0)) cchg=1;

  }
  
  // RS output format changed to new PYTHIA style
  outdec << " " << setw(9) << lundkc;
  outdec << "  ";
  outdec.width(namelength);
  outdec << setiosflags(ios::left)
	 << sname;
  // RS: name for cc paricle
  if ((evtnum.getId()>=0) && (EvtPDL::chargeConj(evtnum)!=evtnum))
    {
      outdec << "  ";
      outdec.width(namelength);
      outdec << ccsname;
    }else{
      // 2+16 spaces
      outdec << "                  ";
    }
  
  outdec << resetiosflags(ios::left);
  outdec << setw(3) << chg;
  outdec << setw(3) << cchg;
  outdec.width(3);
  if (evtnum.getId()>=0) {
    if (EvtPDL::chargeConj(evtnum)==evtnum) {
      outdec << 0;
    }
    else{
      outdec << 1;
    }
  }
  else{
    outdec << 0;
  }
  outdec.setf(ios::fixed,ios::floatfield);
  outdec.precision(5);
  outdec << setw(12) << mass;
  outdec.setf(ios::fixed,ios::floatfield);
  outdec << setw(12) << width;
  outdec.setf(ios::fixed,ios::floatfield);
  outdec.width(12);
  if (fabs(width)<0.0000000001) {
    outdec << 0.0 ;
  }
  else{
    outdec << maxwidth;
  }
  // scientific notation ...  outdec << setw(14) << ctau;
  outdec.setf(ios::scientific,ios::floatfield);
  outdec << "  ";
  outdec << ctau;
  outdec.width(3);
  if (evtnum.getId()>=0) {
    if (ctau>1.0 || rawbrfrsum<0.000001) {  
      stable=0;
    }
  }
  //resonance width treatment
  outdec.width(3);
  outdec << 0;
  outdec.width(3);
  outdec << stable;
  outdec << endl;
  outdec.width(0);
  //outdec.setf(0,0);

}

void EvtPythia::WritePythiaParticle(ofstream &outdec,EvtId ipar,
				    EvtId iparname,int &first){

  int ijetset;

  double br_sum=0.0;

  for(ijetset=0;ijetset<njetsetdecays;ijetset++){
   
    if (jetsetdecays[ijetset]->getParentId()==ipar){
      br_sum+=jetsetdecays[ijetset]->getBranchingFraction();
    }
    if (jetsetdecays[ijetset]->getParentId()!=
	EvtPDL::chargeConj(jetsetdecays[ijetset]->getParentId())&&
	EvtPDL::chargeConj(jetsetdecays[ijetset]->getParentId())==ipar){
      br_sum+=jetsetdecays[ijetset]->getBranchingFraction();
    }


  }

  double br_sum_true=br_sum;

  if (br_sum<0.000001) br_sum=1.0;
  
  for(ijetset=0;ijetset<njetsetdecays;ijetset++){
    if (jetsetdecays[ijetset]->getParentId()==ipar){
      
      double br=jetsetdecays[ijetset]->getBranchingFraction();
      
      int i,daugs[5];
      EvtId cdaugs[5];
      
      for(i=0;i<5;i++){
	
	if(i<jetsetdecays[ijetset]->getNDaug()){
	  daugs[i]=EvtPDL::getStdHep(
			 jetsetdecays[ijetset]->getDaugs()[i]);
	  cdaugs[i]=EvtPDL::chargeConj(jetsetdecays[ijetset]->getDaugs()[i]);
	}
	else{
	  daugs[i]=0;
	}
      }
      
      int channel;
      
      channel=EvtDecayTable::findChannel(EvtPDL::chargeConj(ipar),
			     jetsetdecays[ijetset]->getModelName(),
			     jetsetdecays[ijetset]->getNDaug(),
			     cdaugs,
			     jetsetdecays[ijetset]->getNArg(),
			     jetsetdecays[ijetset]->getArgsStr());     
      
      if (jetsetdecays[ijetset]->getModelName()=="PYTHIA"){
	
	if (first) {
	  first=0;      
	  WritePythiaEntryHeader(outdec,
				 //EvtPDL::getLundKC(iparname),
				 EvtPDL::getStdHep(iparname),
				 iparname,
				 EvtPDL::name(iparname), 
				 EvtPDL::chg3(iparname),
				 0,0,EvtPDL::getMeanMass(ipar),
				 EvtPDL::getWidth(ipar),
				 EvtPDL::getMeanMass(ipar)-EvtPDL::getMinMass(ipar),
				 EvtPDL::getctau(ipar),1,br_sum_true);
	}
	
	int dflag=2;
	
	if (EvtPDL::getStdHep(ipar)<0) {
	  dflag=3;
	  for(i=0;i<jetsetdecays[ijetset]->getNDaug();i++){
	    daugs[i]=EvtPDL::getStdHep(cdaugs[i]);
	  }
	  
	}
	
	/* RS
	  PYTHIA allows to introduce new particles via a call to PYUPDA
	  so no need for this check any more
	  
	  //now lets check to make sure that jetset, lucomp, knows
	  //about all particles!
	  int unknown=0;
	  for(i=0;i<jetsetdecays[ijetset]->getNDaug();i++){
	  if (pycomp_(&daugs[i])==0) {
	  unknown=1;
	  report(ERROR,"EvtGen") << "Pythia (pycomp) does not "
	  << "know the particle:"<<
	  EvtPDL::name(jetsetdecays[ijetset]->getDaugs()[i])<<endl;
	  }
	  }
	  
	  int istdheppar=EvtPDL::getStdHep(ipar);
	  
	  if (pycomp_(&istdheppar)==0){
	  unknown=1;
	  report(ERROR,"EvtGen") << "Pythia (pycomp) does not "
	  << "know the particle:"<<
	  EvtPDL::name(ipar)<<endl;
	  }
	  
	  
	  
	  if (unknown){
	  report(ERROR,"EvtGen") << "Therfore the decay:"<<endl;
	  report(ERROR,"EvtGen") << EvtPDL::name(jetsetdecays[ijetset]->getParentId())<<" -> ";
	  for(i=0;i<jetsetdecays[ijetset]->getNDaug();i++){
	  report(ERROR,"") << EvtPDL::name(jetsetdecays[ijetset]->getDaugs()[i])<<" ";
	  }
	  report(ERROR,"")<<endl;
	  report(ERROR,"EvtGen")<<"Will not be generated."<<endl;
	  return;
	  }
	  */
	
	if (EvtPDL::chargeConj(ipar)==ipar) {
	  dflag=1;
	  //report(INFO,"EvtGen") << EvtPDL::name(iparname) << " dflag=1 because C(ipar)=ipar!"<<endl;
	}
	
	
	//if (channel>=0) {
	//  dflag=1;
	//report(INFO,"EvtGen") << EvtPDL::name(iparname) << " dflag=1 because channel>=0"<<endl;
	//}
	
	//	if (!(EvtPDL::getStdHep(ipar)<0&&channel>=0)){
	if (1){
	  
	  // RS changed format to new PYTHIA one
	  outdec << "          ";
	  outdec.width(5);
	  outdec <<dflag;
	  outdec.width(5);
	  outdec <<(int)jetsetdecays[ijetset]->getArgs()[0];
	  outdec.width(12);
	  if (fabs(br)<0.000000001) {
	    outdec <<"0.00000";
	  }
	  else{
	    outdec <<br/br_sum;
	  }
	  outdec.width(10);
	  outdec <<daugs[0];
	  outdec.width(10);
	  outdec <<daugs[1];
	  outdec.width(10);
	  outdec <<daugs[2];
	  outdec.width(10);
	  outdec <<daugs[3];
	  outdec.width(10);
	  outdec <<daugs[4];
	  outdec<<endl;
	  outdec.width(0);
	}
      }
    }
  }
}

bool
EvtPythia::diquark(int ID)
{
  switch(ID)
    {
    case 1103:
    case 2101:
    case 2103:
    case 2203:
    case 3101:
    case 3103:
    case 3201:
    case 3203:
    case 3303:
    case 4101:
    case 4103:
    case 4201:
    case 4203:
    case 4301:
    case 4303:
    case 4403:
    case 5101:
    case 5103:
    case 5201:
    case 5203:
    case 5301:
    case 5303:
    case 5401:
    case 5403:
    case 5503:
      return true;
      break;
    default:
      return false;
      break;
    }
}

double
EvtPythia::NominalMass(int ID)
{
  // return default mass in PYTHIA
  switch(ID)
    {
    case 1103:
      return 0.77133;
    case 2101:
      return 0.57933;
    case 2103:
      return 0.77133;
    case 2203:
      return 0.77133;
    case 3101:
      return 0.80473;
    case 3103:
      return 0.92953;
    case 3201:
      return 0.80473;
    case 3203:
      return 0.92953;
    case 3303:
      return 1.09361;
    case 4101:
      return 1.96908;
    case 4103:
      return 2.00808;
    case 4201:
      return 1.96908;
    case 4203:
      return 2.00808;
    case 4301:
      return 2.15432;
    case 4303:
      return 2.17967;
    case 4403:
      return 3.27531;
    case 5101:
      return 5.38897;
    case 5103:
      return 5.40145;
    case 5201:
      return 5.38897;
    case 5203:
      return 5.40145;
    case 5301:
      return 5.56725;
    case 5303:
      return 5.57536;
    case 5401:
      return 6.67143;
    case 5403:
      return 6.67397;
    case 5503:
      return 10.07354;
      break;
    default:
      return 0.0;
      break;
    }
}

int
NominalCharge(int ID)
{
  // return default mass in PYTHIA
  switch(ID)
    {
    case 1103:
      return -2;
    case 2101:
      return  1;
    case 2103:
      return  1;
    case 2203:
      return  4;
    case 3101:
      return -2;
    case 3103:
      return -2;
    case 3201:
      return  1;
    case 3203:
      return  1;
    case 3303:
      return -2;
    case 4101:
      return  1;
    case 4103:
      return  1;
    case 4201:
      return  4;
    case 4203:
      return  4;
    case 4301:
      return  1;
    case 4303:
      return  1;
    case 4403:
      return  4;
    case 5101:
      return -2;
    case 5103:
      return -2;
    case 5201:
      return  1;
    case 5203:
      return  1;
    case 5301:
      return -2;
    case 5303:
      return -2;
    case 5401:
      return  1;
    case 5403:
      return  1;
    case 5503:
      return -2;
      break;
    default:
      return 0;
      break;
    }
}

void EvtPythia::MakePythiaFile(char* fname){
  
  EvtId ipar;
  int lundkc;
  
  //int part_list[MAX_PART];
  
  ofstream outdec;
  
  outdec.open(fname);
  
  //outdec << "ERROR;"<<endl;
  //outdec << ";"<<endl;
  //outdec << ";This decayfile has been automatically created by"<<endl;
  //outdec << ";EvtGen from the DECAY.DEC file"<<endl;
  //outdec << ";"<<endl;
  
  int nokcentry;
  
  for(lundkc=1;lundkc<500;lundkc++){
    
    nokcentry=1;
    
    for(size_t iipar=0;iipar<EvtPDL::entries();iipar++){
      
      ipar=EvtId(iipar,iipar);
      //no aliased particles!
      std::string tempStr = EvtPDL::name(ipar);
      EvtId realId = EvtPDL::getId(tempStr);
      if ( realId.isAlias() != 0 ) continue;
      
      if(!(
	       EvtPDL::getStdHep(ipar)==21 ||
	       EvtPDL::getStdHep(ipar)==22 ||
	       EvtPDL::getStdHep(ipar)==23))
	{
	  
	  if (lundkc==EvtPDL::getLundKC(ipar)){
	    
	    nokcentry=0;
	    
	    int first=1;
	    
	    WritePythiaParticle(outdec,ipar,ipar,first);
	    
	    
	    EvtId ipar2=EvtPDL::chargeConj(ipar);
	    
	    
	    if (ipar2!=ipar){
	      WritePythiaParticle(outdec,ipar2,ipar,first);
	    }
	    
	    if (first){
	      WritePythiaEntryHeader(outdec, 
				     //EvtPDL::getLundKC(ipar),
				     EvtPDL::getStdHep(ipar),
				     ipar,
				     EvtPDL::name(ipar),
				     EvtPDL::chg3(ipar),
				     0,0,EvtPDL::getMeanMass(ipar),
				     EvtPDL::getWidth(ipar),
				     EvtPDL::getMeanMass(ipar)-EvtPDL::getMinMass(ipar),
				     EvtPDL::getctau(ipar),0,0.0);
	      
	    }
	  }
	}
    }
    if (lundkc==99999) // Write out diquarks after quarks, but only once
      for(size_t iipar=0;iipar<EvtPDL::entries();iipar++){
	
	ipar=EvtId(iipar,iipar);
	
	if (diquark(EvtPDL::getStdHep(ipar))){
	  
	  nokcentry=0;
	  
	  int first=1;
	  
	  WritePythiaParticle(outdec,ipar,ipar,first);
	  
	  
	  EvtId ipar2=EvtPDL::chargeConj(ipar);
	  
	  
	  if (ipar2!=ipar){
	    WritePythiaParticle(outdec,ipar2,ipar,first);
	  }
	  
	  if (first){
	    WritePythiaEntryHeader(outdec, 
				   EvtPDL::getStdHep(ipar),
				   ipar,
				   EvtPDL::name(ipar),
				   NominalCharge(EvtPDL::getStdHep(ipar)),
				   -1,0,
				   NominalMass(EvtPDL::getStdHep(ipar)),
				   0, 0, 0, 0,0.0);
	    
	  }
	}
      }
    /* if (nokcentry){
       
       WritePythiaEntryHeader(outdec, 
       lundkc,EvtId(-1,-1),"  ",
       0,0,0,EvtPDL::getNominalMass(ipar),0.0,0.0,
       EvtPDL::getctau(ipar),0,0.0);
       
       } */
  }
  outdec.close();
}

void EvtPythia::pythiaInit(int dummy){
  
  //static int first=1; 
  static int first=0; //if first=0 Pythia is not reinitialize   
  if (first){
    
    first=0;
    
    report(INFO,"EvtGen") << "Will initialize Pythia."<<endl;
    for(int i=0;i<ncommand;i++)
      pygive_(commands[i].c_str(),strlen(commands[i].c_str()));
    
    char fname[200];
    
    char hostBuffer[100];
    
    if ( gethostname( hostBuffer, 100 ) != 0 ){
      report(ERROR,"EvtGen") << " couldn't get hostname." << endl;
      strncpy( hostBuffer, "hostnameNotFound", 100 );
    }
    
    char pid[100];
    
    int thePid=getpid();
    
    if ( sprintf( pid, "%d", thePid ) == 0 ){
      report(ERROR,"EvtGen") << " couldn't get process ID." << endl;
      strncpy( pid, "666", 100 );
    }
    
    strcpy(fname,"jet.d-");
    strcat(fname,hostBuffer);
    strcat(fname,"-");
    strcat(fname,pid);
    
    MakePythiaFile(fname);
    evtpythiainit_(fname,strlen(fname));
    initpythia_(&dummy);
    
    if (0==getenv("EVTSAVEJETD")){
      char delcmd[300];
      strcpy(delcmd,"rm -f ");
      strcat(delcmd,fname);
      system(delcmd);
    }
    
    report(INFO,"EvtGen") << "Done initializing Pythia."<<endl;
    
  }

}
