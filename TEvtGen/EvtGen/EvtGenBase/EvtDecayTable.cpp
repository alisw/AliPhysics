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
#include <sstream>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtDecayTable.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtSymTable.hh"
#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtModel.hh"
#include "EvtGenBase/EvtParser.hh"
#include "EvtGenBase/EvtParserXml.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtModelAlias.hh"
#include "EvtGenBase/EvtRadCorr.hh"
#include "EvtGenBase/EvtExtGeneratorCommandsTable.hh"

using std::endl;
using std::fstream;
using std::ifstream;

EvtDecayTable::EvtDecayTable() {
  _decaytable.clear();
}

EvtDecayTable::~EvtDecayTable() {
  _decaytable.clear();
}

EvtDecayTable* EvtDecayTable::getInstance() {

  static EvtDecayTable* theDecayTable = 0;

  if (theDecayTable == 0) {
    theDecayTable = new EvtDecayTable();
  }

  return theDecayTable;

}

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
  EvtExtGeneratorCommandsTable* extGenCommands = EvtExtGeneratorCommandsTable::getInstance();
  std::string colon(":"), equals("=");

  int i;

  report(Severity::Info,"EvtGen") << "In readDecayFile, reading:"<<dec_name.c_str()<<endl;
  
  ifstream fin;
  
  fin.open(dec_name.c_str());
  if (!fin) {
    report(Severity::Error,"EvtGen") << "Could not open "<<dec_name.c_str()<<endl;
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
    report(Severity::Error,"EvtGen") << "Could not find an 'End' in "<<dec_name.c_str()<<endl;
    report(Severity::Error,"EvtGen") << "Will terminate execution."<<endl;
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
	report(Severity::Info,"EvtGen") 
	  << "As requested, PHOTOS will be turned off."<<endl; 
    }
    else if (token=="yesPhotos"){ 
      EvtRadCorr::setAlwaysRadCorr();
      if ( verbose) 
	report(Severity::Info,"EvtGen") 
	  << "As requested, PHOTOS will be turned on for all decays."<<endl; 
    }
    else if (token=="normalPhotos"){ 
      EvtRadCorr::setNormalRadCorr();
      if ( verbose) 
	report(Severity::Info,"EvtGen") 
	  << "As requested, PHOTOS will be turned on only when requested."<<endl; 
    }
    else if (token=="Alias"){

      std::string newname;
      std::string oldname;

      newname=parser.getToken(itoken++);
      oldname=parser.getToken(itoken++);

      EvtId id=EvtPDL::getId(oldname);

      if (id==EvtId(-1,-1)) {
	report(Severity::Error,"EvtGen") <<"Unknown particle name:"<<oldname.c_str()
			       <<" on line "<<parser.getLineofToken(itoken)<<endl;
	report(Severity::Error,"EvtGen") <<"Will terminate execution!"<<endl;
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
	report(Severity::Error,"EvtGen") <<"Unknown particle name:"<<aname.c_str()
			       <<" on line "<<parser.getLineofToken(itoken)<<endl;
	report(Severity::Error,"EvtGen") <<"Will terminate execution!"<<endl;
	::abort();
      }

      if (abar==EvtId(-1,-1)) {
	report(Severity::Error,"EvtGen") <<"Unknown particle name:"<<abarname.c_str()
			       <<" on line "<<parser.getLineofToken(itoken)<<endl;
	report(Severity::Error,"EvtGen") <<"Will terminate execution!"<<endl;
	::abort();
      }


      EvtPDL::aliasChgConj(a,abar);

    } else if (token == "JetSetPar") {

      // Check if any old Pythia 6 commands are present
      std::string pythiaCommand = parser.getToken(itoken++);

      Command command;

      // The old command format is NAME(INT)=VALUE
      int i1 = pythiaCommand.find_first_of("(");
      int i2 = pythiaCommand.find_first_of(")");
      int i3 = pythiaCommand.find_first_of("=");

      std::string pythiaModule = pythiaCommand.substr(0, i1);
      std::string pythiaParam = pythiaCommand.substr(i1+1, i2-i1-1);
      std::string pythiaValue = pythiaCommand.substr(i3+1);

      command["MODULE"] = pythiaModule;
      command["PARAM"] = pythiaParam;
      command["VALUE"] = pythiaValue;

      command["GENERATOR"] = "Both";
      command["VERSION"] = "PYTHIA6";

      extGenCommands->addCommand("PYTHIA", command);

    } else if (modelist.isCommand(token)){

      std::string cnfgstr;

      cnfgstr=parser.getToken(itoken++);

      modelist.storeCommand(token,cnfgstr);

    } else if (token == "PythiaGenericParam" || token == "PythiaAliasParam" ||
	       token == "PythiaBothParam") {

      // Read in any Pythia 8 commands, which will be of the form 
      // pythia<type>Param module:param=value, with no spaces in the parameter 
      // string! Here, <type> specifies whether the command is for generic 
      // decays, alias decays, or both.

      // Pythia 6 commands will be defined by the old JetSetPar command
      // name, which is handled by the modelist.isCommand() statement above.

      std::string pythiaCommand = parser.getToken(itoken++);
      std::string pythiaModule(""), pythiaParam(""), pythiaValue("");

      // Separate out the string into the 3 sections using the delimiters
      // ":" and "=".
      
      std::vector<std::string> pComVect1 = this->splitString(pythiaCommand, colon);

      if (pComVect1.size() == 2) {

	pythiaModule = pComVect1[0];

	std::string pCom2 = pComVect1[1];

	std::vector<std::string> pComVect2 = this->splitString(pCom2, equals);

	if (pComVect2.size() == 2) {

	  pythiaParam = pComVect2[0];
	  pythiaValue = pComVect2[1];

	}

      }

      // Define the Pythia 8 command and pass it to the external generator
      // command list.
      Command command;
      if (token == "PythiaGenericParam") {
	command["GENERATOR"] = "Generic";
      } else if (token == "PythiaAliasParam") {
	command["GENERATOR"] = "Alias";
      } else {
	command["GENERATOR"] = "Both";
      }
      
      command["MODULE"]  = pythiaModule;
      command["PARAM"]   = pythiaParam;
      command["VALUE"]   = pythiaValue;

      command["VERSION"] = "PYTHIA8";
      extGenCommands->addCommand("PYTHIA", command);
      
    } else if (token=="CDecay"){

      std::string name;

      name=parser.getToken(itoken++);
      ipar=EvtPDL::getId(name);

      if (ipar==EvtId(-1,-1)) {
	report(Severity::Error,"EvtGen") <<"Unknown particle name:"<<name.c_str()
			       <<" on line "
			       <<parser.getLineofToken(itoken-1)<<endl;
	report(Severity::Error,"EvtGen") <<"Will terminate execution!"<<endl;
	::abort();
      }

      EvtId cipar=EvtPDL::chargeConj(ipar);

      if (_decaytable[ipar.getAlias()].getNMode()!=0) {
	if ( verbose )
	  report(Severity::Debug,"EvtGen") << 
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
	report(Severity::Info,"EvtGen") << pname.c_str() << endl;
      //There should be at least the mass 
      double newMass=atof(parser.getToken(itoken++).c_str());
      EvtId thisPart = EvtPDL::getId(pname);
      double newWidth=EvtPDL::getMeanMass(thisPart);
      if ( parser.getNToken() > 3 ) newWidth=atof(parser.getToken(itoken++).c_str());

      //Now make the change!
      EvtPDL::reSetMass(thisPart, newMass);
      EvtPDL::reSetWidth(thisPart, newWidth);

      if (verbose )
	report(Severity::Info,"EvtGen") << "Changing particle properties of " <<
	  pname.c_str() << " Mass=" << newMass << " Width="<<newWidth<<endl;

    } else if ( token=="ChangeMassMin") {
      std::string pname;
      pname=parser.getToken(itoken++);
      double tmass=atof(parser.getToken(itoken++).c_str());

      EvtId thisPart = EvtPDL::getId(pname);
      EvtPDL::reSetMassMin(thisPart,tmass);
      if ( verbose )
	report(Severity::Debug,"EvtGen") <<"Refined minimum mass for " << EvtPDL::name(thisPart).c_str() << " to be " << tmass << endl;

    } else if ( token=="ChangeMassMax") {
      std::string pname;
      pname=parser.getToken(itoken++);
      double tmass=atof(parser.getToken(itoken++).c_str());
      EvtId thisPart = EvtPDL::getId(pname);
      EvtPDL::reSetMassMax(thisPart,tmass);
      if ( verbose )
	report(Severity::Debug,"EvtGen") <<"Refined maximum mass for " << EvtPDL::name(thisPart).c_str() << " to be " << tmass << endl;

    } else if ( token=="IncludeBirthFactor") {
      std::string pname;
      pname=parser.getToken(itoken++);
      bool yesno=false;
      if ( parser.getToken(itoken++)=="yes") yesno=true;
      EvtId thisPart = EvtPDL::getId(pname);
      EvtPDL::includeBirthFactor(thisPart,yesno);
      if ( verbose ) {
	if ( yesno ) report(Severity::Debug,"EvtGen") <<"Include birth factor for " << EvtPDL::name(thisPart).c_str() <<endl;
	if ( !yesno ) report(Severity::Debug,"EvtGen") <<"No longer include birth factor for " << EvtPDL::name(thisPart).c_str() <<endl;
      }

    } else if ( token=="IncludeDecayFactor") {
      std::string pname;
      pname=parser.getToken(itoken++);
      bool yesno=false;
      if ( parser.getToken(itoken++)=="yes") yesno=true;
      EvtId thisPart = EvtPDL::getId(pname);
      EvtPDL::includeDecayFactor(thisPart,yesno);
      if ( verbose ) {
	if ( yesno ) report(Severity::Debug,"EvtGen") <<"Include decay factor for " << EvtPDL::name(thisPart).c_str() <<endl;
	if ( !yesno ) report(Severity::Debug,"EvtGen") <<"No longer include decay factor for " << EvtPDL::name(thisPart).c_str() <<endl;
      }
    } else if ( token=="LSNONRELBW") {
      std::string pname;
      pname=parser.getToken(itoken++);
      EvtId thisPart = EvtPDL::getId(pname);
      std::string tstr="NONRELBW";
      EvtPDL::changeLS(thisPart,tstr);
      if ( verbose )
	report(Severity::Debug,"EvtGen") <<"Change lineshape to non-rel BW for " << EvtPDL::name(thisPart).c_str() <<endl;
    } else if ( token=="LSFLAT") {
      std::string pname;
      pname=parser.getToken(itoken++);
      EvtId thisPart = EvtPDL::getId(pname);
      std::string tstr="FLAT";
      EvtPDL::changeLS(thisPart,tstr);
      if (verbose) 
	report(Severity::Debug,"EvtGen") <<"Change lineshape to flat for " << EvtPDL::name(thisPart).c_str() <<endl;
    } else if ( token=="LSMANYDELTAFUNC") {
      std::string pname;
      pname=parser.getToken(itoken++);
      EvtId thisPart = EvtPDL::getId(pname);
      std::string tstr="MANYDELTAFUNC";
      EvtPDL::changeLS(thisPart,tstr);
      if ( verbose )
	report(Severity::Debug,"EvtGen") <<"Change lineshape to spikes for " << EvtPDL::name(thisPart).c_str() <<endl;

    } else if ( token=="BlattWeisskopf") {
      std::string pname;
      pname=parser.getToken(itoken++);
      double tnum=atof(parser.getToken(itoken++).c_str());
      EvtId thisPart = EvtPDL::getId(pname);
      EvtPDL::reSetBlatt(thisPart,tnum);
      if ( verbose )
	report(Severity::Debug,"EvtGen") <<"Redefined Blatt-Weisskopf factor " << EvtPDL::name(thisPart).c_str() << " to be " << tnum << endl;
    } else if ( token=="BlattWeisskopfBirth") {
      std::string pname;
      pname=parser.getToken(itoken++);
      double tnum=atof(parser.getToken(itoken++).c_str());
      EvtId thisPart = EvtPDL::getId(pname);
      EvtPDL::reSetBlattBirth(thisPart,tnum);
      if ( verbose )
	report(Severity::Debug,"EvtGen") <<"Redefined Blatt-Weisskopf birth factor " << EvtPDL::name(thisPart).c_str() << " to be " << tnum << endl;
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
	report(Severity::Debug,"EvtGen") <<"Redefined Partial wave for " << pname.c_str() << " to " << pnameD1.c_str() << " " << pnameD2.c_str() << " ("<<pw<<")"<<endl;
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
        report(Severity::Error,"EvtGen") <<"Unknown particle name:"<<parent.c_str()
                               <<" on line "
                               <<parser.getLineofToken(itoken-1)<<endl;
        report(Severity::Error,"EvtGen") <<"Will terminate execution!"<<endl;
        ::abort();
      }

      if (_decaytable[ipar.getAlias()].getNMode()!=0) {
        report(Severity::Debug,"EvtGen") <<"Redefined decay of "
                               <<parent.c_str()<<endl;
        _decaytable[ipar.getAlias()].removeDecay();
      }


      do{
        
        token=parser.getToken(itoken++);
        
        if (token!="Enddecay"){
          
          i=0;
          while (token.c_str()[i++]!=0){
            if (isalpha(token.c_str()[i])){
              report(Severity::Error,"EvtGen") << 
                "Expected to find a branching fraction or Enddecay "<<
                "but found:"<<token.c_str()<<" on line "<<
                parser.getLineofToken(itoken-1)<<endl;
              report(Severity::Error,"EvtGen") << "Possibly to few arguments to model "<<
                "on previous line!"<<endl;
              report(Severity::Error,"EvtGen") << "Will terminate execution!"<<endl;
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
            
            report(Severity::Info,"EvtGen") << parser.getToken(itoken).c_str()
                                  << " is neither a particle name nor "
                                  << "the name of a model. "<<endl;
            report(Severity::Info,"EvtGen") << "It was encountered on line "<<
              parser.getLineofToken(itoken)<<" of the decay file."<<endl;
            report(Severity::Info,"EvtGen") << "Please fix it. Thank you."<<endl;
            report(Severity::Info,"EvtGen") << "Be sure to check that the "
                                  << "correct case has been used. \n";
            report(Severity::Info,"EvtGen") << "Terminating execution. \n";
            ::abort();
            
            itoken++;
          }
          
          n_daugh=0;
          
          while(EvtPDL::getId(parser.getToken(itoken)).getId()>=0){
            sdaug=parser.getToken(itoken++);
            daught[n_daugh++]=EvtPDL::getId(sdaug);
            if (daught[n_daugh-1]==EvtId(-1,-1)) {
              report(Severity::Error,"EvtGen") <<"Unknown particle name:"<<sdaug.c_str()
                                     <<" on line "<<parser.getLineofToken(itoken)<<endl;
              report(Severity::Error,"EvtGen") <<"Will terminate execution!"<<endl;
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
              report(Severity::Error,"EvtGen") << 
                "Expected to find a model name,"<<
                "found:"<<model.c_str()<<" on line "<<
                parser.getLineofToken(itoken)<<endl;
              report(Severity::Error,"EvtGen") << "Will terminate execution!"<<endl;
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
		  report(Severity::Error,"EvtGen")
		    <<"Reading arguments and found:"<<
		    name.c_str()<<" on line:"<<
		    parser.getLineofToken(itoken-1)<<endl;
		  report(Severity::Error,"EvtGen") 
		    << "Will terminate execution!"<<endl;
		  ::abort();
		}
	      }
	      //int isname=EvtPDL::getId(name).getId()>=0;
	      int ismodel=modelist.isModel(name);
	      if (ismodel) {
		report(Severity::Error,"EvtGen")
		  <<"Expected ';' but found:"<<
		  name.c_str()<<" on line:"<<
		  parser.getLineofToken(itoken-1)<<endl;
		report(Severity::Error,"EvtGen") 
		  << "Most probable error is omitted ';'."<<endl;
		report(Severity::Error,"EvtGen") 
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
	report(Severity::Error,"EvtGen") <<"Unknown particle name:"<<oldname.c_str()
			       <<" on line "<<parser.getLineofToken(itoken)<<endl;
	report(Severity::Error,"EvtGen") <<"Will terminate execution!"<<endl;
	::abort();
      }
      if (newipar==EvtId(-1,-1)) {
	report(Severity::Error,"EvtGen") <<"Unknown particle name:"<<newname.c_str()
			       <<" on line "<<parser.getLineofToken(itoken)<<endl;
	report(Severity::Error,"EvtGen") <<"Will terminate execution!"<<endl;
	::abort();
      }
      if (_decaytable[newipar.getAlias()].getNMode()!=0) {
	report(Severity::Debug,"EvtGen") <<"Redefining decay of "
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
        report(Severity::Error,"EvtGen") <<"Unknown particle name:"<<parent.c_str()
                               <<" on line "
                               <<parser.getLineofToken(itoken-1)<<endl;
        report(Severity::Error,"EvtGen") <<"Will terminate execution!"<<endl;
        ::abort();
      }
       
      if (_decaytable[ipar.getAlias()].getNMode()==0) {
        report(Severity::Debug,"EvtGen") << "No decays to delete for "
                               << parent.c_str() << endl;
      } else {
        report(Severity::Debug,"EvtGen") <<"Deleting selected decays of "
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
	      report(Severity::Error,"EvtGen") <<"Unknown particle name:"<<sdaug.c_str()
				     <<" on line "<<parser.getLineofToken(itoken)<<endl;
	      report(Severity::Error,"EvtGen") <<"Will terminate execution!"<<endl;
	      ::abort();
	    }
	  }
	  token = parser.getToken(itoken);
	  if (token != ";") {
	    report(Severity::Error,"EvtGen")
	      <<"Expected ';' but found:"<<
	      token <<" on line:"<<
	      parser.getLineofToken(itoken-1)<<endl;
	    report(Severity::Error,"EvtGen") 
	      << "Most probable error is omitted ';'."<<endl;
	    report(Severity::Error,"EvtGen") 
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

      report(Severity::Error,"EvtGen") << "Found unknown command:'"<<token.c_str()<<"' on line "
			     <<parser.getLineofToken(itoken)<<endl;
      report(Severity::Error,"EvtGen") << "Will terminate execution!"<<endl;
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
	report(Severity::Info,"EvtGen") << "Given allowed decays, resetting minMass " << EvtPDL::name(temp).c_str() << " " 
			      << EvtPDL::getMinMass(temp) << " to " << minMass << endl;
      EvtPDL::reSetMassMin(temp,minMass);
    }
  }
}

void EvtDecayTable::readXMLDecayFile(const std::string dec_name, bool verbose){
  if ( _decaytable.size() < EvtPDL::entries() ) _decaytable.resize(EvtPDL::entries());
  EvtModel &modelist=EvtModel::instance();
  EvtExtGeneratorCommandsTable* extGenCommands = EvtExtGeneratorCommandsTable::getInstance();

  EvtParserXml parser;
  parser.open(dec_name);

  EvtId ipar;
  std::string decayParent = "";
  double brfrSum = 0.;
  std::vector<EvtModelAlias> modelAliasList;
  bool endReached = false;

  while(parser.readNextTag()) {
      //TAGS FOUND UNDER DATA
      if(parser.getParentTagTitle() == "data") {
        if(parser.getTagTitle() == "photos") {
          std::string usage = parser.readAttribute("usage");
          if(usage == "always") {
            EvtRadCorr::setAlwaysRadCorr();
            if ( verbose )
              report(Severity::Info,"EvtGen")
                << "As requested, PHOTOS will be turned on for all decays."<<endl;
          } else if(usage == "never") {
            EvtRadCorr::setNeverRadCorr();
            if ( verbose )
              report(Severity::Info,"EvtGen")
                << "As requested, PHOTOS will be turned off."<<endl;
          } else {
            EvtRadCorr::setNormalRadCorr();
            if ( verbose )
              report(Severity::Info,"EvtGen")
                << "As requested, PHOTOS will be turned on only when requested."<<endl;
          }

        } else if(parser.getTagTitle() == "alias") {
          std::string alias = parser.readAttribute("name");
          std::string particle = parser.readAttribute("particle");
          checkParticle(particle);
          EvtId id=EvtPDL::getId(particle);

          EvtPDL::alias(id,alias);
          if ( _decaytable.size() < EvtPDL::entries() ) _decaytable.resize(EvtPDL::entries());

        } else if(parser.getTagTitle() == "modelAlias") {
          std::vector<std::string> modelArgList;

          std::string alias = parser.readAttribute("name");
          std::string model = parser.readAttribute("model");
          std::string paramStr = parser.readAttribute("params");
          std::istringstream paramStream(paramStr);

          std::string param;

          if(paramStr=="") {
            EvtDecayBase* fcn = modelist.getFcn(model);
            int i(0);
            std::string paramName = fcn->getParamName(0);
            while(paramName!="") {
              param = parser.readAttribute(paramName,fcn->getParamDefault(i));
              if(param=="") break;
              modelArgList.push_back(param);
              ++i;
              paramName = fcn->getParamName(i);
            }
          } else {
            while(std::getline(paramStream, param, ' ')) {
              modelArgList.push_back(param);
            }
          }
          EvtModelAlias newAlias(alias,model,modelArgList);
          modelAliasList.push_back(newAlias);

        } else if(parser.getTagTitle() == "chargeConj") {
          std::string particle = parser.readAttribute("particle");
          std::string conjugate = parser.readAttribute("conjugate");

          EvtId a=EvtPDL::getId(particle);
          EvtId abar=EvtPDL::getId(conjugate);

          checkParticle(particle);
          checkParticle(conjugate);

          EvtPDL::aliasChgConj(a,abar);

        } else if(parser.getTagTitle() == "conjDecay") {
          std::string particle = parser.readAttribute("particle");

          EvtId a=EvtPDL::getId(particle);
          EvtId abar=EvtPDL::chargeConj(a);

          checkParticle(particle);
          checkParticle(abar.getName());

          if (_decaytable[a.getAlias()].getNMode()!=0) {
            if ( verbose )
              report(Severity::Debug,"EvtGen") <<
                "Redefined decay of "<<particle.c_str()<<" in ConjDecay"<<endl;

            _decaytable[a.getAlias()].removeDecay();
          }

          //take contents of abar and conjugate and store in a
          _decaytable[a.getAlias()].makeChargeConj(&_decaytable[abar.getAlias()]);
          
        } else if(parser.getTagTitle() == "define") {
          std::string name = parser.readAttribute("name");
          std::string value = parser.readAttribute("value");
          EvtSymTable::define(name,value);

        } else if(parser.getTagTitle() == "particle") {
          std::string name = parser.readAttribute("name");
          double mass = parser.readAttributeDouble("mass");
          double width = parser.readAttributeDouble("width");
          double minMass = parser.readAttributeDouble("massMin");
          double maxMass = parser.readAttributeDouble("massMax");
          std::string birthFactor = parser.readAttribute("includeBirthFactor");
          std::string decayFactor = parser.readAttribute("includeDecayFactor");
          std::string lineShape = parser.readAttribute("lineShape");
          double blattWeisskopfD = parser.readAttributeDouble("blattWeisskopfFactor");
          double blattWeisskopfB = parser.readAttributeDouble("blattWeisskopfBirth");

          EvtId thisPart = EvtPDL::getId(name);
          checkParticle(name);

          if(mass != -1) {
            EvtPDL::reSetMass(thisPart, mass);
            report(Severity::Debug,"EvtGen") <<"Refined mass for " << EvtPDL::name(thisPart).c_str() << " to be " << mass << endl;
          }
          if(width != -1) {
            EvtPDL::reSetWidth(thisPart, width);
            report(Severity::Debug,"EvtGen") <<"Refined width for " << EvtPDL::name(thisPart).c_str() << " to be " << width << endl;
          }
          if(minMass != -1) {
            EvtPDL::reSetMassMin(thisPart,minMass);
            report(Severity::Debug,"EvtGen") <<"Refined minimum mass for " << EvtPDL::name(thisPart).c_str() << " to be " << minMass << endl;
          }
          if(maxMass != -1) {
            EvtPDL::reSetMassMax(thisPart,maxMass);
            report(Severity::Debug,"EvtGen") <<"Refined maximum mass for " << EvtPDL::name(thisPart).c_str() << " to be " << maxMass << endl;
          }
          if(!birthFactor.empty()) {
            EvtPDL::includeBirthFactor(thisPart,stringToBoolean(birthFactor));
            if(verbose) {
              if(stringToBoolean(birthFactor)) {
                report(Severity::Debug,"EvtGen") <<"Include birth factor for " << EvtPDL::name(thisPart).c_str() <<endl;
              } else {
                report(Severity::Debug,"EvtGen") <<"No longer include birth factor for " << EvtPDL::name(thisPart).c_str() <<endl;
              }
            }
          }
          if(!decayFactor.empty()) {
            EvtPDL::includeDecayFactor(thisPart,stringToBoolean(decayFactor));
            if(verbose) {
              if(stringToBoolean(decayFactor)) {
                report(Severity::Debug,"EvtGen") <<"Include decay factor for " << EvtPDL::name(thisPart).c_str() <<endl;
              } else {
                report(Severity::Debug,"EvtGen") <<"No longer include decay factor for " << EvtPDL::name(thisPart).c_str() <<endl;
              }
            }
          }
          if(!lineShape.empty()) {
            EvtPDL::changeLS(thisPart,lineShape);
            if ( verbose )
              report(Severity::Debug,"EvtGen") <<"Change lineshape to " << lineShape << " for " << EvtPDL::name(thisPart).c_str() <<endl;
          }
          if(blattWeisskopfD != -1) {
            EvtPDL::reSetBlatt(thisPart,blattWeisskopfD);
            if ( verbose )
              report(Severity::Debug,"EvtGen") <<"Redefined Blatt-Weisskopf factor "
                                     << EvtPDL::name(thisPart).c_str() << " to be " << blattWeisskopfD << endl;
          }
          if(blattWeisskopfB != -1) {
            EvtPDL::reSetBlattBirth(thisPart,blattWeisskopfB);
            if ( verbose )
              report(Severity::Debug,"EvtGen") <<"Redefined Blatt-Weisskopf birth factor "
                                     << EvtPDL::name(thisPart).c_str() << " to be " << blattWeisskopfB << endl;
          }
        } else if(parser.getTagTitle() == "lineShapePW") {
          std::string parent = parser.readAttribute("parent");
          std::string daug1 = parser.readAttribute("daug1");
          std::string daug2 = parser.readAttribute("daug2");
          int pw = parser.readAttributeInt("pw");

          checkParticle(parent);
          checkParticle(daug1);
          checkParticle(daug2);

          EvtId thisPart = EvtPDL::getId(parent);
          EvtId thisD1 = EvtPDL::getId(daug1);
          EvtId thisD2 = EvtPDL::getId(daug2);

          EvtPDL::setPWForDecay(thisPart,pw,thisD1,thisD2);
          EvtPDL::setPWForBirthL(thisD1,pw,thisPart,thisD2);
          EvtPDL::setPWForBirthL(thisD2,pw,thisPart,thisD1);
          if ( verbose )
            report(Severity::Debug,"EvtGen") <<"Redefined Partial wave for " << parent.c_str() << " to "
                                   << daug1.c_str() << " " << daug2.c_str() << " ("<<pw<<")"<<endl;

        } else if(parser.getTagTitle() == "decay") { //start of a particle
          brfrSum = 0.;
          decayParent = parser.readAttribute("name");
          checkParticle(decayParent);
          ipar=EvtPDL::getId(decayParent);

          if (_decaytable[ipar.getAlias()].getNMode()!=0) {
            report(Severity::Debug,"EvtGen") <<"Redefined decay of "
                                   <<decayParent.c_str()<<endl;
            _decaytable[ipar.getAlias()].removeDecay();
          }

        } else if(parser.getTagTitle() == "copyDecay") {
          std::string particle = parser.readAttribute("particle");
          std::string copy = parser.readAttribute("copy");

          EvtId newipar=EvtPDL::getId(particle);
          EvtId oldipar=EvtPDL::getId(copy);

          checkParticle(particle);
          checkParticle(copy);

          if (_decaytable[newipar.getAlias()].getNMode()!=0) {
            report(Severity::Debug,"EvtGen") <<"Redefining decay of "
                                   <<particle<<endl;
            _decaytable[newipar.getAlias()].removeDecay();
          }
          _decaytable[newipar.getAlias()] = _decaytable[oldipar.getAlias()];

        } else if(parser.getTagTitle() == "removeDecay") {
          decayParent = parser.readAttribute("particle");
          checkParticle(decayParent);
          ipar=EvtPDL::getId(decayParent);

          if (_decaytable[ipar.getAlias()].getNMode()==0) {
            report(Severity::Debug,"EvtGen") << "No decays to delete for "
                                   << decayParent.c_str() << endl;
          } else {
            report(Severity::Debug,"EvtGen") <<"Deleting selected decays of "
                                   <<decayParent.c_str()<<endl;
          }

        } else if(parser.getTagTitle() == "pythiaParam") {
          Command command;
          command["GENERATOR"] = parser.readAttribute("generator");
          command["MODULE"]    = parser.readAttribute("module");
          command["PARAM"]     = parser.readAttribute("param");
          command["VALUE"]     = parser.readAttribute("value");
          command["VERSION"]   = "PYTHIA8";
          extGenCommands->addCommand("PYTHIA", command);

        } else if(parser.getTagTitle() == "pythia6Param") {
          Command command;
          command["GENERATOR"] = parser.readAttribute("generator");
          command["MODULE"]    = parser.readAttribute("module");
          command["PARAM"]     = parser.readAttribute("param");
          command["VALUE"]     = parser.readAttribute("value");
          command["VERSION"]   = "PYTHIA6";
          extGenCommands->addCommand("PYTHIA", command);

        } else if(parser.getTagTitle() == "/data") { //end of data
          endReached = true;
          parser.close();
          break;
        } else if(parser.getTagTitle() == "Title" || parser.getTagTitle() == "Details"
               || parser.getTagTitle() == "Author" || parser.getTagTitle() == "Version"
          //the above tags are expected to be in the XML decay file but are not used by EvtGen
               || parser.getTagTitle() == "dalitzDecay" || parser.getTagTitle() == "copyDalitz") {
          //the above tags are only used by EvtGenModels/EvtDalitzTable
        } else {  report(Severity::Info,"EvtGen") << "Unknown tag "<<parser.getTagTitle()
                  <<" found in XML decay file near line "<<parser.getLineNumber()<<". Tag will be ignored."<<endl;
        }
      //TAGS FOUND UNDER DECAY
      } else if(parser.getParentTagTitle() == "decay") {
        if(parser.getTagTitle() == "channel") { //start of a channel
          int nDaughters = 0;
          EvtId daughter[MAX_DAUG];

          EvtDecayBase* temp_fcn_new;
          std::string temp_fcn_new_model;
          std::vector<std::string> temp_fcn_new_args;

          double brfr = parser.readAttributeDouble("br");
          std::string daugStr = parser.readAttribute("daughters");
          std::istringstream daugStream(daugStr);
          std::string model = parser.readAttribute("model");
          std::string paramStr = parser.readAttribute("params");
          std::istringstream paramStream(paramStr);
          bool decVerbose = parser.readAttributeBool("verbose");
          bool decPhotos = parser.readAttributeBool("photos");
          bool decSummary = parser.readAttributeBool("summary");

          std::string daugh;
          while(std::getline(daugStream, daugh, ' ')) {
            checkParticle(daugh);
            daughter[nDaughters++] = EvtPDL::getId(daugh);
          }

          int modelAlias = -1;
          for(size_t iAlias=0;iAlias<modelAliasList.size();iAlias++){
            if ( modelAliasList[iAlias].matchAlias(model) ) {
              modelAlias=iAlias;
              break;
            }
          }

          if ( modelAlias==-1 ) {
            if(!modelist.isModel(model)){
              report(Severity::Error,"EvtGen") <<
                "Expected to find a model name near line "<<parser.getLineNumber()<<","<<
                "found:"<<model.c_str()<<endl;
              report(Severity::Error,"EvtGen") << "Will terminate execution!"<<endl;
              ::abort();
            }
          } else {
            model=modelAliasList[modelAlias].getName();
          }

          temp_fcn_new_model = model;
          temp_fcn_new = modelist.getFcn(model);

          if(decPhotos) temp_fcn_new->setPHOTOS();
          if(decVerbose) temp_fcn_new->setVerbose();
          if(decSummary) temp_fcn_new->setSummary();

          int ierr;
          if(modelAlias == -1) {
            std::string param;
            if(paramStr == "") {
              int i(0);
              std::string paramName = temp_fcn_new->getParamName(0); 
              while(paramName != "") {
                param = parser.readAttribute(paramName,temp_fcn_new->getParamDefault(i));
                if(param == "") break; //params must be added in order so we can't just skip the missing ones
                temp_fcn_new_args.push_back(EvtSymTable::get(param,ierr));
                if (ierr) {
                  report(Severity::Error,"EvtGen")
                    <<"Reading arguments near line "<<parser.getLineNumber()<<" and found:"<<
                    param.c_str()<<endl;
                  report(Severity::Error,"EvtGen")
                    << "Will terminate execution!"<<endl;
                  ::abort();
                }
                ++i;
                paramName = temp_fcn_new->getParamName(i);
              }

            } else {//if the params are not set seperately
              while(std::getline(paramStream, param, ' ')) {
                temp_fcn_new_args.push_back(EvtSymTable::get(param,ierr));
                if (ierr) {
                  report(Severity::Error,"EvtGen")
                    <<"Reading arguments near line "<<parser.getLineNumber()<<" and found:"<<
                    param.c_str()<<endl;
                  report(Severity::Error,"EvtGen")
                    << "Will terminate execution!"<<endl;
                  ::abort();
                }
              }
            }
          } else {
            std::vector<std::string> copyMe=modelAliasList[modelAlias].getArgList();
            temp_fcn_new_args=copyMe;
          }

          brfrSum+=brfr;

          temp_fcn_new->saveDecayInfo(ipar,nDaughters,
                                      daughter,
                                      temp_fcn_new_args.size(),
                                      temp_fcn_new_args,
                                      temp_fcn_new_model,
                                      brfr);

          double massMin=0.0;

          for (int i=0;i<temp_fcn_new->nRealDaughters();i++){
            if ( EvtPDL::getMinMass(daughter[i])>0.0001 ){
              massMin+=EvtPDL::getMinMass(daughter[i]);
            } else {
              massMin+=EvtPDL::getMeanMass(daughter[i]);
            }
          }

          _decaytable[ipar.getAlias()].addMode(temp_fcn_new,brfrSum,massMin);

        } else if(parser.getTagTitle() == "/decay") { //end of a particle
          _decaytable[ipar.getAlias()].finalize();
        } else report(Severity::Info,"EvtGen") << "Unexpected tag "<<parser.getTagTitle()
                                     <<" found in XML decay file near line "<<parser.getLineNumber()<<". Tag will be ignored."<<endl;
      //TAGS FOUND UNDER REMOVEDECAY
      } else if(parser.getParentTagTitle() == "removeDecay") {
        if(parser.getTagTitle() == "channel") { //start of a channel
          int nDaughters = 0;
          EvtId daughter[MAX_DAUG];

          std::string daugStr = parser.readAttribute("daughters");
          std::istringstream daugStream(daugStr);

          std::string daugh;
          while(std::getline(daugStream, daugh, ' ')) {
            checkParticle(daugh);
            daughter[nDaughters++] = EvtPDL::getId(daugh);
          }

          EvtDecayBase* temp_fcn_new = modelist.getFcn("PHSP");
          std::vector<std::string> temp_fcn_new_args;
          std::string temp_fcn_new_model("PHSP");
          temp_fcn_new->saveDecayInfo(ipar, nDaughters,
                                      daughter,
                                      0,
                                      temp_fcn_new_args,
                                      temp_fcn_new_model,
                                      0.);
          _decaytable[ipar.getAlias()].removeMode(temp_fcn_new);
        } else if(parser.getTagTitle() != "/removeDecay") {
          report(Severity::Info,"EvtGen") << "Unexpected tag "<<parser.getTagTitle()
                                <<" found in XML decay file near line "<<parser.getLineNumber()<<". Tag will be ignored."<<endl;
        }
      }
  }//while lines in file

  if(!endReached) {
    report(Severity::Info,"EvtGen") << "Either the decay file ended prematurely or the file is badly formed.\n"
                          <<"Error occured near line"<<parser.getLineNumber()<<endl;
    ::abort();
  }

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
        report(Severity::Info,"EvtGen") << "Given allowed decays, resetting minMass " << EvtPDL::name(temp).c_str() << " "
                              << EvtPDL::getMinMass(temp) << " to " << minMass << endl;
      EvtPDL::reSetMassMin(temp,minMass);
    }
  }
}

bool EvtDecayTable::stringToBoolean(std::string valStr) {
  return (valStr == "true" || valStr == "1" || valStr == "on" || valStr == "yes");
}

void EvtDecayTable::checkParticle(std::string particle) {
  if (EvtPDL::getId(particle)==EvtId(-1,-1)) {
    report(Severity::Error,"EvtGen") <<"Unknown particle name:"<<particle.c_str()<<endl;
    report(Severity::Error,"EvtGen") <<"Will terminate execution!"<<endl;
    ::abort();
  }
}

EvtDecayBase* EvtDecayTable::findDecayModel(EvtId id, int modeInt) {

  int aliasInt = id.getAlias();

  EvtDecayBase* theModel = this->findDecayModel(aliasInt, modeInt);

  return theModel;

}

EvtDecayBase* EvtDecayTable::findDecayModel(int aliasInt, int modeInt) {

  EvtDecayBase* theModel(0);

  if (aliasInt >= 0 && aliasInt < (int) EvtPDL::entries()) {

    theModel = _decaytable[aliasInt].getDecayModel(modeInt);

  }

  return theModel;

}

bool EvtDecayTable::hasPythia(EvtId id) {

  bool hasPythia = this->hasPythia(id.getAlias());
  return hasPythia;

}

bool EvtDecayTable::hasPythia(int aliasInt) {

  bool hasPythia(false);
  if (aliasInt >= 0 && aliasInt < (int) EvtPDL::entries()) {

    hasPythia = _decaytable[aliasInt].isJetSet();

  }
  
  return hasPythia;

}

int EvtDecayTable::getNModes(EvtId id) {

  int nModes = this->getNModes(id.getAlias());
  return nModes;

}

int EvtDecayTable::getNModes(int aliasInt) {

  int nModes(0);

  if (aliasInt >= 0 && aliasInt < (int) EvtPDL::entries()) {

    nModes = _decaytable[aliasInt].getNMode();
  }

  return nModes;

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
   
std::vector<std::string> EvtDecayTable::splitString(std::string& theString, 
						    std::string& splitter) {

  // Code from STLplus
  std::vector<std::string> result;

  if (!theString.empty() && !splitter.empty()) {

    for (std::string::size_type offset = 0;;) {

      std::string::size_type found = theString.find(splitter, offset);

      if (found != std::string::npos) {
	std::string tmpString = theString.substr(offset, found-offset);
        if (tmpString.size() > 0) {result.push_back(tmpString);}
        offset = found + splitter.size();
      } else {
	std::string tmpString = theString.substr(offset, theString.size()-offset);
        if (tmpString.size() > 0) {result.push_back(tmpString);}
        break;
      }
    }
  }

  return result;
}

     
