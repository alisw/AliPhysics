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
// Module: EvtDecayTable.cc
//
// Description:
//
// Modification history:
//
//    DJL/RYD     September 25, 1996         Module created
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
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtDecayTable.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtSymTable.hh"
#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtModel.hh"
#include "EvtGenBase/EvtParser.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtModelAlias.hh"
#include "EvtGenBase/EvtRadCorr.hh"
using std::endl;
using std::fstream;
using std::ifstream;

std::vector<EvtParticleDecayList> EvtDecayTable::_decaytable;

int EvtDecayTable::getNMode(int ipar){
   return _decaytable[ipar].getNMode();
} 

EvtDecayBase* EvtDecayTable::getDecay(int ipar, int imode){
  return _decaytable[ipar].getDecayModel(imode);
}

void EvtDecayTable::printSummary(){
  
  for(size_t i=0;i<EvtPDL::entries();i++){
    _decaytable[i].printSummary();
  }

}

EvtDecayBase* EvtDecayTable::getDecayFunc(EvtParticle *p){
  int partnum;
  
  partnum=p->getId().getAlias();

  if ( _decaytable[partnum].getNMode()==0 ) return 0;
  return _decaytable[partnum].getDecayModel(p);

}

void EvtDecayTable::readDecayFile(const std::string dec_name, bool verbose){

  if ( _decaytable.size() < EvtPDL::entries() ) _decaytable.resize(EvtPDL::entries());
  EvtModel &modelist=EvtModel::instance();
  int i;

  report(INFO,"EvtGen") << "In readDecayFile, reading:"<<dec_name.c_str()<<endl;
  
  ifstream fin;
  
  fin.open(dec_name.c_str());
  if (!fin) {
    report(ERROR,"EvtGen") << "Could not open "<<dec_name.c_str()<<endl;
  }
  fin.close();

  EvtParser parser;
  parser.read(dec_name);

  int itok;

  int hasend=0;

  std::string token;

  for(itok=0;itok<parser.getNToken();itok++){

    token=parser.getToken(itok);
    
    if (token=="End") hasend=1;

  }

  if (!hasend){
    report(ERROR,"EvtGen") << "Could not find an 'End' in "<<dec_name.c_str()<<endl;
    report(ERROR,"EvtGen") << "Will terminate execution."<<endl;
    ::abort();
  }



  std::string model,parent,sdaug;  

  EvtId ipar;

  int n_daugh;
  EvtId daught[MAX_DAUG];
  double brfr;

  int itoken=0;

  std::vector<EvtModelAlias> modelAliasList;

  
  do{

    token=parser.getToken(itoken++);

    //Easy way to turn off photos... Lange September 5, 2000
    if (token=="noPhotos"){ 
      EvtRadCorr::setNeverRadCorr();
      if ( verbose )
	report(INFO,"EvtGen") 
	  << "As requested, PHOTOS will be turned off."<<endl; 
    }
    else if (token=="yesPhotos"){ 
      EvtRadCorr::setAlwaysRadCorr();
      if ( verbose) 
	report(INFO,"EvtGen") 
	  << "As requested, PHOTOS will be turned on for all decays."<<endl; 
    }
    else if (token=="normalPhotos"){ 
      EvtRadCorr::setNormalRadCorr();
      if ( verbose) 
	report(INFO,"EvtGen") 
	  << "As requested, PHOTOS will be turned on only when requested."<<endl; 
    }
    else if (token=="Alias"){

      std::string newname;
      std::string oldname;

      newname=parser.getToken(itoken++);
      oldname=parser.getToken(itoken++);

      EvtId id=EvtPDL::getId(oldname);

      if (id==EvtId(-1,-1)) {
	report(ERROR,"EvtGen") <<"Unknown particle name:"<<oldname.c_str()
			       <<" on line "<<parser.getLineofToken(itoken)<<endl;
	report(ERROR,"EvtGen") <<"Will terminate execution!"<<endl;
	::abort();
      }

      EvtPDL::alias(id,newname);
      if ( _decaytable.size() < EvtPDL::entries() ) _decaytable.resize(EvtPDL::entries());

    } else if (token=="ModelAlias"){
      std::vector<std::string> modelArgList;

      std::string aliasName=parser.getToken(itoken++);
      std::string modelName=parser.getToken(itoken++);

      std::string nameTemp;
      do{
	nameTemp=parser.getToken(itoken++);
	if (nameTemp!=";") {
	  modelArgList.push_back(nameTemp);
	}
      }while(nameTemp!=";");
      EvtModelAlias newAlias(aliasName,modelName,modelArgList);
      modelAliasList.push_back(newAlias);
    } else if (token=="ChargeConj"){

      std::string aname;
      std::string abarname;

      aname=parser.getToken(itoken++);
      abarname=parser.getToken(itoken++);

      EvtId a=EvtPDL::getId(aname);
      EvtId abar=EvtPDL::getId(abarname);

      if (a==EvtId(-1,-1)) {
	report(ERROR,"EvtGen") <<"Unknown particle name:"<<aname.c_str()
			       <<" on line "<<parser.getLineofToken(itoken)<<endl;
	report(ERROR,"EvtGen") <<"Will terminate execution!"<<endl;
	::abort();
      }

      if (abar==EvtId(-1,-1)) {
	report(ERROR,"EvtGen") <<"Unknown particle name:"<<abarname.c_str()
			       <<" on line "<<parser.getLineofToken(itoken)<<endl;
	report(ERROR,"EvtGen") <<"Will terminate execution!"<<endl;
	::abort();
      }


      EvtPDL::aliasChgConj(a,abar);

    } else if (modelist.isCommand(token)){

      std::string cnfgstr;

      cnfgstr=parser.getToken(itoken++);

      modelist.storeCommand(token,cnfgstr);

    } else if (token=="CDecay"){

      std::string name;

      name=parser.getToken(itoken++);
      ipar=EvtPDL::getId(name);

      if (ipar==EvtId(-1,-1)) {
	report(ERROR,"EvtGen") <<"Unknown particle name:"<<name.c_str()
			       <<" on line "
			       <<parser.getLineofToken(itoken-1)<<endl;
	report(ERROR,"EvtGen") <<"Will terminate execution!"<<endl;
	::abort();
      }

      EvtId cipar=EvtPDL::chargeConj(ipar);

      if (_decaytable[ipar.getAlias()].getNMode()!=0) {
	if ( verbose )
	  report(DEBUG,"EvtGen") << 
	    "Redefined decay of "<<name.c_str()<<" in CDecay"<<endl;

	_decaytable[ipar.getAlias()].removeDecay();
      }

      //take contents of cipar and conjugate and store in ipar
      _decaytable[ipar.getAlias()].makeChargeConj(&_decaytable[cipar.getAlias()]);

    } else if (token=="Define"){

      std::string name;

      name=parser.getToken(itoken++);
      //      value=atof(parser.getToken(itoken++).c_str());

      EvtSymTable::define(name,parser.getToken(itoken++));

      //New code Lange April 10, 2001 - allow the user
      //to change particle definitions of EXISTING
      //particles on the fly
    } else if (token=="Particle"){

      std::string pname;
      pname=parser.getToken(itoken++);
      if ( verbose )
	report(INFO,"EvtGen") << pname.c_str() << endl;
      //There should be at least the mass 
      double newMass=atof(parser.getToken(itoken++).c_str());
      EvtId thisPart = EvtPDL::getId(pname);
      double newWidth=EvtPDL::getMeanMass(thisPart);
      if ( parser.getNToken() > 3 ) newWidth=atof(parser.getToken(itoken++).c_str());

      //Now make the change!
      EvtPDL::reSetMass(thisPart, newMass);
      EvtPDL::reSetWidth(thisPart, newWidth);

      if (verbose )
	report(INFO,"EvtGen") << "Changing particle properties of " <<
	  pname.c_str() << " Mass=" << newMass << " Width="<<newWidth<<endl;

    } else if ( token=="ChangeMassMin") {
      std::string pname;
      pname=parser.getToken(itoken++);
      double tmass=atof(parser.getToken(itoken++).c_str());

      EvtId thisPart = EvtPDL::getId(pname);
      EvtPDL::reSetMassMin(thisPart,tmass);
      if ( verbose )
	report(DEBUG,"EvtGen") <<"Refined minimum mass for " << EvtPDL::name(thisPart).c_str() << " to be " << tmass << endl;

    } else if ( token=="ChangeMassMax") {
      std::string pname;
      pname=parser.getToken(itoken++);
      double tmass=atof(parser.getToken(itoken++).c_str());
      EvtId thisPart = EvtPDL::getId(pname);
      EvtPDL::reSetMassMax(thisPart,tmass);
      if ( verbose )
	report(DEBUG,"EvtGen") <<"Refined maximum mass for " << EvtPDL::name(thisPart).c_str() << " to be " << tmass << endl;

    } else if ( token=="IncludeBirthFactor") {
      std::string pname;
      pname=parser.getToken(itoken++);
      bool yesno=false;
      if ( parser.getToken(itoken++).c_str()=="yes") yesno=true;
      EvtId thisPart = EvtPDL::getId(pname);
      EvtPDL::includeBirthFactor(thisPart,yesno);
      if ( verbose ) {
	if ( yesno ) report(DEBUG,"EvtGen") <<"Include birth factor for " << EvtPDL::name(thisPart).c_str() <<endl;
	if ( !yesno ) report(DEBUG,"EvtGen") <<"No longer include birth factor for " << EvtPDL::name(thisPart).c_str() <<endl;
      }

    } else if ( token=="IncludeDecayFactor") {
      std::string pname;
      pname=parser.getToken(itoken++);
      bool yesno=false;
      if ( parser.getToken(itoken++).c_str()=="yes") yesno=true;
      EvtId thisPart = EvtPDL::getId(pname);
      EvtPDL::includeDecayFactor(thisPart,yesno);
      if ( verbose ) {
	if ( yesno ) report(DEBUG,"EvtGen") <<"Include decay factor for " << EvtPDL::name(thisPart).c_str() <<endl;
	if ( !yesno ) report(DEBUG,"EvtGen") <<"No longer include decay factor for " << EvtPDL::name(thisPart).c_str() <<endl;
      }
    } else if ( token=="LSNONRELBW") {
      std::string pname;
      pname=parser.getToken(itoken++);
      EvtId thisPart = EvtPDL::getId(pname);
      std::string tstr="NONRELBW";
      EvtPDL::changeLS(thisPart,tstr);
      if ( verbose )
	report(DEBUG,"EvtGen") <<"Change lineshape to non-rel BW for " << EvtPDL::name(thisPart).c_str() <<endl;
    } else if ( token=="SP8LSFIX") {
      //this was a bug, but preserve functionality as not to confuse people...
      std::string pname;
      pname=parser.getToken(itoken++);
      EvtId thisPart = EvtPDL::getId(pname);
      EvtPDL::fixLSForSP8(thisPart);
      if ( verbose )
	report(DEBUG,"EvtGen") <<"Fixed lineshape for SP8 --from D.Lange,J.Smith " << EvtPDL::name(thisPart).c_str() <<endl;

    } else if ( token=="SP6LSFIX") {
      std::string pname;
      pname=parser.getToken(itoken++);
      EvtId thisPart = EvtPDL::getId(pname);
      EvtPDL::fixLSForSP8(thisPart);
      if ( verbose ) 
	report(DEBUG,"EvtGen") <<"Fixed lineshape for SP8 --from D.Lange,J.Smith " << EvtPDL::name(thisPart).c_str() <<endl;

    } else if ( token=="LSFLAT") {
      std::string pname;
      pname=parser.getToken(itoken++);
      EvtId thisPart = EvtPDL::getId(pname);
      std::string tstr="FLAT";
      EvtPDL::changeLS(thisPart,tstr);
      if (verbose) 
	report(DEBUG,"EvtGen") <<"Change lineshape to flat for " << EvtPDL::name(thisPart).c_str() <<endl;
    } else if ( token=="LSMANYDELTAFUNC") {
      std::string pname;
      pname=parser.getToken(itoken++);
      EvtId thisPart = EvtPDL::getId(pname);
      std::string tstr="MANYDELTAFUNC";
      EvtPDL::changeLS(thisPart,tstr);
      if ( verbose )
	report(DEBUG,"EvtGen") <<"Change lineshape to spikes for " << EvtPDL::name(thisPart).c_str() <<endl;

    } else if ( token=="BlattWeisskopf") {
      std::string pname;
      pname=parser.getToken(itoken++);
      double tnum=atof(parser.getToken(itoken++).c_str());
      EvtId thisPart = EvtPDL::getId(pname);
      EvtPDL::reSetBlatt(thisPart,tnum);
      if ( verbose )
	report(DEBUG,"EvtGen") <<"Redefined Blatt-Weisskopf factor " << EvtPDL::name(thisPart).c_str() << " to be " << tnum << endl;
    } else if ( token=="SetLineshapePW") {
      std::string pname;
      pname=parser.getToken(itoken++);
      EvtId thisPart = EvtPDL::getId(pname);
      std::string pnameD1=parser.getToken(itoken++);
      EvtId thisD1 = EvtPDL::getId(pnameD1);
      std::string pnameD2=parser.getToken(itoken++);
      EvtId thisD2 = EvtPDL::getId(pnameD2);
      int pw=atoi(parser.getToken(itoken++).c_str());
      if ( verbose ) 
	report(DEBUG,"EvtGen") <<"Redefined Partial wave for " << pname.c_str() << " to " << pnameD1.c_str() << " " << pnameD2.c_str() << " ("<<pw<<")"<<endl;
      EvtPDL::setPWForDecay(thisPart,pw,thisD1,thisD2);
      EvtPDL::setPWForBirthL(thisD1,pw,thisPart,thisD2);
      EvtPDL::setPWForBirthL(thisD2,pw,thisPart,thisD1);


    } else if (token=="Decay") {

      std::string temp_fcn_new_model;

      EvtDecayBase* temp_fcn_new;
      
      double brfrsum=0.0;

  

      parent=parser.getToken(itoken++);
      ipar=EvtPDL::getId(parent);

      if (ipar==EvtId(-1,-1)) {
	report(ERROR,"EvtGen") <<"Unknown particle name:"<<parent.c_str()
			       <<" on line "
			       <<parser.getLineofToken(itoken-1)<<endl;
	report(ERROR,"EvtGen") <<"Will terminate execution!"<<endl;
	::abort();
      }

      if (_decaytable[ipar.getAlias()].getNMode()!=0) {
	report(DEBUG,"EvtGen") <<"Redefined decay of "
			       <<parent.c_str()<<endl;
	_decaytable[ipar.getAlias()].removeDecay();
      }


      do{

        token=parser.getToken(itoken++);

        if (token!="Enddecay"){

	  i=0;
	  while (token.c_str()[i++]!=0){
	    if (isalpha(token.c_str()[i])){
	      report(ERROR,"EvtGen") << 
		"Expected to find a branching fraction or Enddecay "<<
		"but found:"<<token.c_str()<<" on line "<<
		parser.getLineofToken(itoken-1)<<endl;
	      report(ERROR,"EvtGen") << "Possibly to few arguments to model "<<
		"on previous line!"<<endl;
	      report(ERROR,"EvtGen") << "Will terminate execution!"<<endl;
	      ::abort();
	    }
	  }
	  
          brfr=atof(token.c_str());

	  int isname=EvtPDL::getId(parser.getToken(itoken)).getId()>=0;
          int ismodel=modelist.isModel(parser.getToken(itoken));

          if (!(isname||ismodel)){
	    //see if this is an aliased model
	    for(size_t iAlias=0;iAlias<modelAliasList.size();iAlias++){
	      if ( modelAliasList[iAlias].matchAlias(parser.getToken(itoken)) ) {
		ismodel=2;
		break;
	      }
	    }
	  }

          if (!(isname||ismodel)){

	    report(INFO,"EvtGen") << parser.getToken(itoken).c_str()
	     << " is neither a particle name nor "
	     << "the name of a model. "<<endl;
	    report(INFO,"EvtGen") << "It was encountered on line "<<
	      parser.getLineofToken(itoken)<<" of the decay file."<<endl;
	    report(INFO,"EvtGen") << "Please fix it. Thank you."<<endl;
	    report(INFO,"EvtGen") << "Be sure to check that the "
	     << "correct case has been used. \n";
	    report(INFO,"EvtGen") << "Terminating execution. \n";
	    ::abort();

	    itoken++;
	  }

          n_daugh=0;

	  while(EvtPDL::getId(parser.getToken(itoken)).getId()>=0){
            sdaug=parser.getToken(itoken++);
            daught[n_daugh++]=EvtPDL::getId(sdaug);
	    if (daught[n_daugh-1]==EvtId(-1,-1)) {
	      report(ERROR,"EvtGen") <<"Unknown particle name:"<<sdaug.c_str()
				     <<" on line "<<parser.getLineofToken(itoken)<<endl;
	      report(ERROR,"EvtGen") <<"Will terminate execution!"<<endl;
	      ::abort();
	    }
          }

	  
          model=parser.getToken(itoken++);


          int photos=0;
          int verbose=0;
          int summary=0;
	  
	  do{
	    if (model=="PHOTOS"){
	      photos=1;
	      model=parser.getToken(itoken++);
	    }
	    if (model=="VERBOSE"){
	      verbose=1;
	      model=parser.getToken(itoken++);
	    }
	    if (model=="SUMMARY"){
	      summary=1;
	      model=parser.getToken(itoken++);
	    }
	  }while(model=="PHOTOS"||
		 model=="VERBOSE"||
		 model=="SUMMARY");

	  //see if this is an aliased model
	  int foundAnAlias=-1;
	  for(size_t iAlias=0;iAlias<modelAliasList.size();iAlias++){
	    if ( modelAliasList[iAlias].matchAlias(model) ) {
	      foundAnAlias=iAlias;
	      break;
	    }
	  }

	  if ( foundAnAlias==-1 ) {
	    if(!modelist.isModel(model)){
	      report(ERROR,"EvtGen") << 
		"Expected to find a model name,"<<
		"found:"<<model.c_str()<<" on line "<<
		parser.getLineofToken(itoken)<<endl;
	      report(ERROR,"EvtGen") << "Will terminate execution!"<<endl;
	      ::abort();
	    }
	  }
	  else{
	    model=modelAliasList[foundAnAlias].getName();
	  }

	  temp_fcn_new_model=model;
	  temp_fcn_new=modelist.getFcn(model);


          if (photos){
	    temp_fcn_new->setPHOTOS();
	  }
          if (verbose){
	    temp_fcn_new->setVerbose();
	  }
          if (summary){
	    temp_fcn_new->setSummary();
	  }
	  

	  std::vector<std::string> temp_fcn_new_args;

	  std::string name;
	  int ierr;

	  if ( foundAnAlias==-1 ) {
	    do{
	      name=parser.getToken(itoken++);
	      if (name!=";") {
		temp_fcn_new_args.push_back(EvtSymTable::get(name,ierr));
		if (ierr) {
		  report(ERROR,"EvtGen")
		    <<"Reading arguments and found:"<<
		    name.c_str()<<" on line:"<<
		    parser.getLineofToken(itoken-1)<<endl;
		  report(ERROR,"EvtGen") 
		    << "Will terminate execution!"<<endl;
		  ::abort();
		}
	      }
	      //int isname=EvtPDL::getId(name).getId()>=0;
	      int ismodel=modelist.isModel(name);
	      if (ismodel) {
		report(ERROR,"EvtGen")
		  <<"Expected ';' but found:"<<
		  name.c_str()<<" on line:"<<
		  parser.getLineofToken(itoken-1)<<endl;
		report(ERROR,"EvtGen") 
		  << "Most probable error is omitted ';'."<<endl;
		report(ERROR,"EvtGen") 
		  << "Will terminate execution!"<<endl;
		::abort();
	      }
	    }while(name!=";");
	  }
	  else{
	    std::vector<std::string> copyMe=modelAliasList[foundAnAlias].getArgList();
	    temp_fcn_new_args=copyMe;
	    itoken++;
	  }
	  //Found one decay.

	  brfrsum+=brfr;

	  temp_fcn_new->saveDecayInfo(ipar,n_daugh,
				      daught,
				      temp_fcn_new_args.size(),
				      temp_fcn_new_args,
				      temp_fcn_new_model,
				      brfr);

	  double massmin=0.0;

	  //          for (i=0;i<n_daugh;i++){
          for (i=0;i<temp_fcn_new->nRealDaughters();i++){
	    if ( EvtPDL::getMinMass(daught[i])>0.0001 ){
              massmin+=EvtPDL::getMinMass(daught[i]);
	    } else {
              massmin+=EvtPDL::getMeanMass(daught[i]);
	    }  
	  } 
	  
	  _decaytable[ipar.getAlias()].addMode(temp_fcn_new,brfrsum,massmin);
	  

	}
      } while(token!="Enddecay");      

      _decaytable[ipar.getAlias()].finalize();

    }
    // Allow copying of decays from one particle to another; useful
    // in combination with RemoveDecay
    else if (token=="CopyDecay") {
      std::string newname;
      std::string oldname;
      
      newname=parser.getToken(itoken++);
      oldname=parser.getToken(itoken++);
      
      EvtId newipar=EvtPDL::getId(newname);
      EvtId oldipar=EvtPDL::getId(oldname);
      
      if (oldipar==EvtId(-1,-1)) {
	report(ERROR,"EvtGen") <<"Unknown particle name:"<<oldname.c_str()
			       <<" on line "<<parser.getLineofToken(itoken)<<endl;
	report(ERROR,"EvtGen") <<"Will terminate execution!"<<endl;
	::abort();
      }
      if (newipar==EvtId(-1,-1)) {
	report(ERROR,"EvtGen") <<"Unknown particle name:"<<newname.c_str()
			       <<" on line "<<parser.getLineofToken(itoken)<<endl;
	report(ERROR,"EvtGen") <<"Will terminate execution!"<<endl;
	::abort();
      }
      if (_decaytable[newipar.getAlias()].getNMode()!=0) {
	report(DEBUG,"EvtGen") <<"Redefining decay of "
			       <<newname<<endl;
	_decaytable[newipar.getAlias()].removeDecay();
      }
      _decaytable[newipar.getAlias()] = _decaytable[oldipar.getAlias()];
    }
    // Enable decay deletion; intended primarily for aliases
    // Peter Onyisi, March 2008
    else if (token=="RemoveDecay") {
      parent = parser.getToken(itoken++);
      ipar = EvtPDL::getId(parent);
      
      if (ipar==EvtId(-1,-1)) {
        report(ERROR,"EvtGen") <<"Unknown particle name:"<<parent.c_str()
                               <<" on line "
                               <<parser.getLineofToken(itoken-1)<<endl;
        report(ERROR,"EvtGen") <<"Will terminate execution!"<<endl;
        ::abort();
      }
       
      if (_decaytable[ipar.getAlias()].getNMode()==0) {
        report(DEBUG,"EvtGen") << "No decays to delete for "
                               << parent.c_str() << endl;
      } else {
        report(DEBUG,"EvtGen") <<"Deleting selected decays of "
                               <<parent.c_str()<<endl;
      }
      
      do {
        token = parser.getToken(itoken);
        
        if (token != "Enddecay") {
	  n_daugh = 0;
	  while (EvtPDL::getId(parser.getToken(itoken)).getId() >= 0) {
	    sdaug = parser.getToken(itoken++);
	    daught[n_daugh++] = EvtPDL::getId(sdaug);
	    if (daught[n_daugh-1]==EvtId(-1,-1)) {
	      report(ERROR,"EvtGen") <<"Unknown particle name:"<<sdaug.c_str()
				     <<" on line "<<parser.getLineofToken(itoken)<<endl;
	      report(ERROR,"EvtGen") <<"Will terminate execution!"<<endl;
	      ::abort();
	    }
	  }
	  token = parser.getToken(itoken);
	  if (token != ";") {
	    report(ERROR,"EvtGen")
	      <<"Expected ';' but found:"<<
	      token <<" on line:"<<
	      parser.getLineofToken(itoken-1)<<endl;
	    report(ERROR,"EvtGen") 
	      << "Most probable error is omitted ';'."<<endl;
	    report(ERROR,"EvtGen") 
	      << "Will terminate execution!"<<endl;
	    ::abort();
	  }
	  token = parser.getToken(itoken++);
	  EvtDecayBase* temp_fcn_new = modelist.getFcn("PHSP");
	  std::vector<std::string> temp_fcn_new_args;
	  std::string temp_fcn_new_model("PHSP");
	  temp_fcn_new->saveDecayInfo(ipar, n_daugh,
                                      daught,
                                      0,
                                      temp_fcn_new_args,
                                      temp_fcn_new_model,
                                      0.);
	  _decaytable[ipar.getAlias()].removeMode(temp_fcn_new);
        }
      } while (token != "Enddecay");
      itoken++;
    }
    else if (token!="End"){

      report(ERROR,"EvtGen") << "Found unknown command:'"<<token.c_str()<<"' on line "
			     <<parser.getLineofToken(itoken)<<endl;
      report(ERROR,"EvtGen") << "Will terminate execution!"<<endl;
      ::abort();

    }

  } while ((token!="End")&&itoken!=parser.getNToken());

  //Now we may need to reset the minimum mass for some particles????

  for (size_t ii=0; ii<EvtPDL::entries(); ii++){
    EvtId temp(ii,ii);
    int nModTot=getNMode(ii);
    //no decay modes
    if ( nModTot == 0 ) continue;
    //0 width?
    if ( EvtPDL::getWidth(temp) < 0.0000001 ) continue;
    int jj;
    double minMass=EvtPDL::getMaxMass(temp);
    for (jj=0; jj<nModTot; jj++) {
      double tmass=_decaytable[ii].getDecay(jj).getMassMin();
      if ( tmass< minMass) minMass=tmass;
    }
    if ( minMass > EvtPDL::getMinMass(temp) ) {
      if ( verbose )
	report(INFO,"EvtGen") << "Given allowed decays, resetting minMass " << EvtPDL::name(temp).c_str() << " " 
			      << EvtPDL::getMinMass(temp) << " to " << minMass << endl;
      EvtPDL::reSetMassMin(temp,minMass);
    }
  }
}

int  EvtDecayTable::findChannel(EvtId parent, std::string model, 
				int ndaug, EvtId *daugs, 
				int narg, std::string *args){

   int i,j,right;
   EvtId daugs_scratch[50];
   int nmatch,k;

   for(i=0;i<_decaytable[parent.getAlias()].getNMode();i++){

     right=1;

     right=right&&model==_decaytable[parent.getAlias()].
		    getDecay(i).getDecayModel()->getModelName();
     right=right&&(ndaug==_decaytable[parent.getAlias()].
	     getDecay(i).getDecayModel()->getNDaug());
     right=right&&(narg==_decaytable[parent.getAlias()].
	     getDecay(i).getDecayModel()->getNArg());

     if ( right ){

       

       for(j=0;j<ndaug;j++){
	 daugs_scratch[j]=daugs[j];
       }

       nmatch=0;

       for(j=0;j<_decaytable[parent.getAlias()].
	     getDecay(i).getDecayModel()->getNDaug();j++){

         for(k=0;k<ndaug;k++){
	   if (daugs_scratch[k]==_decaytable[parent.getAlias()].
	       getDecay(i).getDecayModel()->getDaug(j)){
             daugs_scratch[k]=EvtId(-1,-1);
             nmatch++;
	     break;
	   }
	 }
       } 

       right=right&&(nmatch==ndaug);

       for(j=0;j<_decaytable[parent.getAlias()].
	     getDecay(i).getDecayModel()->getNArg();j++){
	 right=right&&(args[j]==_decaytable[parent.getAlias()].
		 getDecay(i).getDecayModel()->getArgStr(j));
       } 
     }
     if (right) return i;
   }
   return -1;
}

int  EvtDecayTable::inChannelList(EvtId parent, int ndaug, EvtId *daugs){

   int i,j,k;
   EvtId daugs_scratch[MAX_DAUG];

   int dsum=0;
   for(i=0;i<ndaug;i++){
     dsum+=daugs[i].getAlias();
   }

   int nmatch;

   int ipar=parent.getAlias();

   int nmode=_decaytable[ipar].getNMode();

   for(i=0;i<nmode;i++){

     EvtDecayBase* thedecaymodel=_decaytable[ipar].getDecay(i).getDecayModel();

     if (thedecaymodel->getDSum()==dsum){

       int nd=thedecaymodel->getNDaug();

       if (ndaug==nd){
	 for(j=0;j<ndaug;j++){
	   daugs_scratch[j]=daugs[j];
	 }
	 nmatch=0;
	 for(j=0;j<nd;j++){
	   for(k=0;k<ndaug;k++){
	     if (EvtId(daugs_scratch[k])==thedecaymodel->getDaug(j)){
	       daugs_scratch[k]=EvtId(-1,-1);
	       nmatch++;
	       break;
	     }
	   }
	 } 
         if ((nmatch==ndaug)&&
             (!
              ((thedecaymodel->getModelName()=="JETSET")||
               (thedecaymodel->getModelName()=="PYTHIA")))){
           return i;
         }
       }
     }
   }

   return -1;
}
   
      
