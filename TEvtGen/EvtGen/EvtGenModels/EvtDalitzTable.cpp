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
// Module: EvtGen/EvtGenericDalitz.hh
//
// Description: Model to describe a generic dalitz decay
//
// Modification history:
//
//    DCC     16 December, 2011         Module created
//
//------------------------------------------------------------------------

#include "EvtGenModels/EvtDalitzTable.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtParserXml.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtSpinType.hh"
#include "EvtGenBase/EvtDalitzPlot.hh"
#include "EvtGenBase/EvtCyclic3.hh"

#include <stdlib.h>
#include <sstream>

using std::endl;
using std::fstream;
using std::ifstream;

EvtDalitzTable::EvtDalitzTable() {
  _dalitztable.clear();
  _readFiles.clear();
}

EvtDalitzTable::~EvtDalitzTable() {
  _dalitztable.clear();
  _readFiles.clear();
}

EvtDalitzTable* EvtDalitzTable::getInstance(const std::string dec_name, bool verbose) { 

  static EvtDalitzTable* theDalitzTable = 0;

  if(theDalitzTable == 0) {
    theDalitzTable = new EvtDalitzTable();
  }

  if(!theDalitzTable->fileHasBeenRead(dec_name)) {
    theDalitzTable->readXMLDecayFile(dec_name,verbose);
  }

  return theDalitzTable;

}

bool EvtDalitzTable::fileHasBeenRead(const std::string dec_name) {
  std::vector<std::string>::iterator i = _readFiles.begin();
  for( ; i!=_readFiles.end(); i++) {
    if((*i).compare(dec_name) == 0) {
      return true;
    }
  }
  return false;
}

void EvtDalitzTable::readXMLDecayFile(const std::string dec_name, bool verbose){

  if (verbose) {
    report(Severity::Info,"EvtGen")<<"EvtDalitzTable: Reading in xml parameter file "<<dec_name<<endl;
  }

  _readFiles.push_back(dec_name);

  EvtDalitzDecayInfo* dalitzDecay = 0;
  double probMax = 0;
  EvtId ipar;
  std::string decayParent = "";
  std::string daugStr = "";
  EvtId daughter[3];

  EvtDalitzPlot dp;
  EvtComplex cAmp;
  std::vector< std::pair<EvtCyclic3::Pair,EvtCyclic3::Pair> > angAndResPairs;
  std::string shape("");
  EvtSpinType::spintype spinType(EvtSpinType::SCALAR);
  double mass(0.), width(0.), FFp(0.), FFr(0.);
  std::vector<EvtFlatteParam> flatteParams;
  //Nonres parameters
  double alpha(0.);
  //LASS parameters
  double aLass(0.), rLass(0.), BLass(0.), phiBLass(0.), RLass(0.), phiRLass(0.), cutoffLass(-1.);

  EvtParserXml parser;
  parser.open(dec_name);

  bool endReached = false;

  while(parser.readNextTag()) {
    //TAGS FOUND UNDER DATA
    if(parser.getParentTagTitle() == "data") {
      if(parser.getTagTitle() == "dalitzDecay") {
        int nDaughters = 0;

        decayParent = parser.readAttribute("particle");
        daugStr = parser.readAttribute("daughters");
        probMax = parser.readAttributeDouble("probMax",-1);

        checkParticle(decayParent);
        ipar=EvtPDL::getId(decayParent);

        std::istringstream daugStream(daugStr);

        std::string daugh;
        while(std::getline(daugStream, daugh, ' ')) {
          checkParticle(daugh);
          daughter[nDaughters++] = EvtPDL::getId(daugh);
        }

        if(nDaughters!=3) {
          report(Severity::Error,"EvtGen") <<
                "Expected to find three daughters for dalitzDecay of "<<decayParent<<" near line "
                <<parser.getLineNumber()<<", "<<"found "<<nDaughters<<endl;
              report(Severity::Error,"EvtGen") << "Will terminate execution!"<<endl;
              ::abort();
        }

        double m_d1 = EvtPDL::getMass(daughter[0]), m_d2 = EvtPDL::getMass(daughter[1]), m_d3 = EvtPDL::getMass(daughter[2]), M = EvtPDL::getMass(ipar);

        dp = EvtDalitzPlot( m_d1, m_d2, m_d3, M );

        dalitzDecay = new EvtDalitzDecayInfo(daughter[0],daughter[1],daughter[2]);

      } else if(parser.getTagTitle() == "copyDalitz") {
        int nDaughters = 0;
        EvtId daughter[3];
        int nCopyDaughters = 0;
        EvtId copyDaughter[3];

        decayParent = parser.readAttribute("particle");
        daugStr = parser.readAttribute("daughters");

        std::string copyParent = parser.readAttribute("copy");
        std::string copyDaugStr = parser.readAttribute("copyDaughters");

        checkParticle(decayParent);
        ipar=EvtPDL::getId(decayParent);

        checkParticle(copyParent);
        EvtId copypar=EvtPDL::getId(copyParent);

        std::istringstream daugStream(daugStr);
        std::istringstream copyDaugStream(copyDaugStr);

        std::string daugh;
        while(std::getline(daugStream, daugh, ' ')) {
          checkParticle(daugh);
          daughter[nDaughters++] = EvtPDL::getId(daugh);
        }
        while(std::getline(copyDaugStream, daugh, ' ')) {
          checkParticle(daugh);
          copyDaughter[nCopyDaughters++] = EvtPDL::getId(daugh);
        }

        if(nDaughters!=3 || nCopyDaughters!=3) {
          report(Severity::Error,"EvtGen") <<
                "Expected to find three daughters for copyDecay of "<<decayParent<<
                " from "<<copyParent<<" near line "<<parser.getLineNumber()<<
                ", "<<"found "<<nDaughters<<" and "<<nCopyDaughters<<endl;
          report(Severity::Error,"EvtGen") << "Will terminate execution!"<<endl;
              ::abort();
        }

        copyDecay(ipar, daughter, copypar, copyDaughter);

      } else if(parser.getTagTitle() == "/data") { //end of data
        endReached = true;
        parser.close();
        break;
      }
    //TAGS FOUND UNDER DALITZDECAY
    } else if(parser.getParentTagTitle() == "dalitzDecay") {
      if(parser.getTagTitle() == "resonance") {

        flatteParams.clear();

        //Amplitude
        EvtComplex ampFactor(parser.readAttributeDouble("ampFactorReal",1.),
                             parser.readAttributeDouble("ampFactorImag",0.));
        double mag = parser.readAttributeDouble("mag",-999.);
        double phase = parser.readAttributeDouble("phase",-999.);
        double real = parser.readAttributeDouble("real",-999.);
        double imag = parser.readAttributeDouble("imag",-999.);

        if((real!=-999. || imag!=-999.) && mag==-999. && phase==-999.) {
          if(real==-999.) { real = 0; }
          if(imag==-999.) { imag = 0; }
          mag = sqrt(real*real + imag*imag);
          phase = atan2(imag,real) * EvtConst::radToDegrees;
        }
        if( mag==-999. ) {
          mag = 1.;
        }
        if( phase==-999. ) {
          phase = 0.;
        }

        cAmp = ampFactor*EvtComplex(cos(phase*1.0/EvtConst::radToDegrees),sin(phase*1.0/EvtConst::radToDegrees))*mag;

        //Resonance particle properties
        mass = 0.;
        width = 0.;
        spinType = EvtSpinType::SCALAR;

        std::string particle = parser.readAttribute("particle");
        if(particle != "") {
          EvtId resId = EvtPDL::getId(particle);
          if(resId == EvtId(-1,-1)) {
            report(Severity::Error,"EvtGen") <<"Unknown particle name:"<<particle.c_str()<<endl;
            report(Severity::Error,"EvtGen") <<"Will terminate execution!"<<endl;
            ::abort();
          } else {
            mass = EvtPDL::getMeanMass(resId);
            width = EvtPDL::getWidth(resId);
            spinType = EvtPDL::getSpinType(resId);
          }
        }

        width = parser.readAttributeDouble("width",width);
        mass = parser.readAttributeDouble("mass",mass);
        switch(parser.readAttributeInt("spin",-1)) {
        case -1://not set here
          break;
        case 0:
          spinType = EvtSpinType::SCALAR;
          break;
        case 1:
          spinType = EvtSpinType::VECTOR;
          break;
        case 2:
          spinType = EvtSpinType::TENSOR;
          break;
        default:
          report(Severity::Error,"EvtGen") << "Unsupported spin near line "<<parser.getLineNumber()<<" of XML file."<<endl;
          ::abort();
        }

        //Shape and form factors
        shape = parser.readAttribute("shape");
        FFp = parser.readAttributeDouble("BlattWeisskopfFactorParent",0.0);
        FFr = parser.readAttributeDouble("BlattWeisskopfFactorResonance",1.5);

        //Shape specific attributes
        if(shape=="NonRes_Exp") {
          alpha = parser.readAttributeDouble("alpha",0.0);
        }
        if(shape=="LASS") {
          aLass = parser.readAttributeDouble("a",2.07);
          rLass = parser.readAttributeDouble("r",3.32);
          BLass = parser.readAttributeDouble("B",1.0);
          phiBLass = parser.readAttributeDouble("phiB",0.0);
          RLass = parser.readAttributeDouble("R",1.0);
          phiRLass = parser.readAttributeDouble("phiR",0.0);
          cutoffLass = parser.readAttributeDouble("cutoff",-1.0);
        }

        //Daughter pairs for resonance
        angAndResPairs.clear();

        std::string resDaugStr = parser.readAttribute("resDaughters");
        if(resDaugStr != "") {
          std::istringstream resDaugStream(resDaugStr);
          std::string resDaug;
          int nResDaug(0);
          EvtId resDaughter[2];
          while(std::getline(resDaugStream, resDaug, ' ')) {
            checkParticle(resDaug);
            resDaughter[nResDaug++] = EvtPDL::getId(resDaug);
          }
          if(nResDaug != 2) {
            report(Severity::Error,"EvtGen") << "Resonance must have exactly 2 daughters near line "<<parser.getLineNumber()<<" of XML file."<<endl;
            ::abort();
          }
          int nRes = getDaughterPairs(resDaughter, daughter, angAndResPairs);
          if(nRes==0) {
            report(Severity::Error,"EvtGen") << "Resonance daughters must match decay daughters near line "<<parser.getLineNumber()<<" of XML file."<<endl;
            ::abort();
          }
          if(parser.readAttributeBool("normalise",true)) cAmp /= sqrt(nRes);
        }
        if(angAndResPairs.empty()) {
          switch(parser.readAttributeInt("daughterPair")) {
          case 1:
            angAndResPairs.push_back(std::make_pair(EvtCyclic3::BC,EvtCyclic3::AB));
            break;
          case 2:
            angAndResPairs.push_back(std::make_pair(EvtCyclic3::CA,EvtCyclic3::BC));
            break;
          case 3:
            angAndResPairs.push_back(std::make_pair(EvtCyclic3::AB,EvtCyclic3::CA));
            break;
          default:
            if(shape == "NonRes") { //We don't expect a pair for non-resonant terms but add dummy values for convenience
              angAndResPairs.push_back(std::make_pair(EvtCyclic3::AB,EvtCyclic3::AB));
              break;
            }
            report(Severity::Error,"EvtGen") << "Daughter pair must be 1, 2 or 3 near line "<<parser.getLineNumber()<<" of XML file."<<endl;
            ::abort();
          }
        }

        if(parser.isTagInline()) {
          std::vector< std::pair<EvtCyclic3::Pair,EvtCyclic3::Pair> >::iterator it = angAndResPairs.begin();
          for( ; it != angAndResPairs.end(); it++) {
            std::pair<EvtCyclic3::Pair,EvtCyclic3::Pair> pairs = *it;
            EvtDalitzReso resonance = getResonance(shape, dp, pairs.first, pairs.second, spinType, mass, width, FFp, FFr, alpha, aLass, rLass, BLass, phiBLass, RLass, phiRLass, cutoffLass);
            dalitzDecay->addResonance(cAmp,resonance);
          }
        }
      } else if(parser.getTagTitle() == "/dalitzDecay") {
        if(probMax < 0) {
          report(Severity::Info,"EvtGen") << "probMax is not defined for " << decayParent << " -> " << daugStr << endl;
          report(Severity::Info,"EvtGen") << "Will now estimate probMax. This may take a while. Once probMax is calculated, update the XML file to skip this step in future." << endl;
          probMax = calcProbMax(dp,dalitzDecay);
        }
        dalitzDecay->setProbMax(probMax);
        addDecay(ipar, *dalitzDecay);
        delete dalitzDecay;
        dalitzDecay = 0;
      } else if(verbose) {
        report(Severity::Info,"EvtGen") << "Unexpected tag "<<parser.getTagTitle()
                              <<" found in XML decay file near line "
                              <<parser.getLineNumber()<<". Tag will be ignored."<<endl;
      }
    //TAGS FOUND UNDER RESONANCE
    } else if(parser.getParentTagTitle() == "resonance"){
      if(parser.getTagTitle() == "flatteParam") {
        EvtFlatteParam param(parser.readAttributeDouble("mass1"),
                             parser.readAttributeDouble("mass2"),
                             parser.readAttributeDouble("g"));
        flatteParams.push_back(param);
      } else if(parser.getTagTitle() == "/resonance") {
        std::vector< std::pair<EvtCyclic3::Pair,EvtCyclic3::Pair> >::iterator it = angAndResPairs.begin();
        for( ; it != angAndResPairs.end(); it++) {
          std::pair<EvtCyclic3::Pair,EvtCyclic3::Pair> pairs = *it;
          EvtDalitzReso resonance = getResonance(shape, dp, pairs.first, pairs.second, spinType, mass, width, FFp, FFr, alpha, aLass, rLass, BLass, phiBLass, RLass, phiRLass, cutoffLass);

          std::vector<EvtFlatteParam>::iterator flatteIt = flatteParams.begin();
          for( ; flatteIt != flatteParams.end(); flatteIt++) {
            resonance.addFlatteParam((*flatteIt));
          }

          dalitzDecay->addResonance(cAmp,resonance);
        }
      }
    }
  }

  if(!endReached) {
    report(Severity::Error,"EvtGen") << "Either the decay file ended prematurely or the file is badly formed.\n"
                          <<"Error occured near line"<<parser.getLineNumber()<<endl;
    ::abort();
  }
}

void EvtDalitzTable::checkParticle(std::string particle) {
  if (EvtPDL::getId(particle)==EvtId(-1,-1)) {
    report(Severity::Error,"EvtGen") <<"Unknown particle name:"<<particle.c_str()<<endl;
    report(Severity::Error,"EvtGen") <<"Will terminate execution!"<<endl;
    ::abort();
  }
}

void EvtDalitzTable::addDecay(EvtId parent, const EvtDalitzDecayInfo& dec) {
  if(_dalitztable.find(parent)!=_dalitztable.end()) {
    _dalitztable[parent].push_back(dec);
  } else {
    _dalitztable[parent].push_back(dec);
  }
}

void EvtDalitzTable::copyDecay(EvtId parent, EvtId* daughters,
                               EvtId copy, EvtId* copyd) {
  EvtDalitzDecayInfo decay(daughters[0],daughters[1],daughters[2]);
  std::vector<EvtDalitzDecayInfo> copyTable = getDalitzTable(copy);
  std::vector<EvtDalitzDecayInfo>::iterator i = copyTable.begin();
  for( ; i != copyTable.end(); i++) {
    EvtId daughter1 = (*i).daughter1();
    EvtId daughter2 = (*i).daughter2();
    EvtId daughter3 = (*i).daughter3();

    if((copyd[0] == daughter1 && copyd[1] == daughter2 && copyd[2] == daughter3) ||
       (copyd[0] == daughter1 && copyd[1] == daughter3 && copyd[2] == daughter2) ||
       (copyd[0] == daughter2 && copyd[1] == daughter1 && copyd[2] == daughter3) ||
       (copyd[0] == daughter2 && copyd[1] == daughter3 && copyd[2] == daughter1) ||
       (copyd[0] == daughter3 && copyd[1] == daughter1 && copyd[2] == daughter2) ||
       (copyd[0] == daughter3 && copyd[1] == daughter2 && copyd[2] == daughter1)) {
      decay.setProbMax((*i).getProbMax());
      std::vector<std::pair<EvtComplex, EvtDalitzReso> >::const_iterator j = (*i).getResonances().begin();
      for( ; j != (*i).getResonances().end(); j++) {
        decay.addResonance((*j));
      }
      addDecay(parent,decay);
      return;
    }
  }
  //if we get here then there was no match
  report(Severity::Error,"EvtGen") << "Did not find dalitz decays for particle:"
         <<copy<<"\n";
}

std::vector<EvtDalitzDecayInfo> EvtDalitzTable::getDalitzTable(const EvtId& parent) {
  std::vector<EvtDalitzDecayInfo> table;
  if ( _dalitztable.find(parent)!=_dalitztable.end() ) {
    table=_dalitztable[parent];
  }

  if (table.empty()){
    report(Severity::Error,"EvtGen") << "Did not find dalitz decays for particle:"
         <<parent<<"\n";
  }

  return table;
}


EvtDalitzReso EvtDalitzTable::getResonance(std::string shape, EvtDalitzPlot dp, EvtCyclic3::Pair angPair, EvtCyclic3::Pair resPair,
                                           EvtSpinType::spintype spinType, double mass, double width, double FFp, double FFr, double alpha,
                                           double aLass, double rLass, double BLass, double phiBLass, double RLass, double phiRLass, double cutoffLass) {
  if( shape=="RBW" || shape=="RBW_CLEO") {
    return EvtDalitzReso( dp, angPair, resPair, spinType, mass, width, EvtDalitzReso::RBW_CLEO, FFp, FFr );
  } else if( shape=="RBW_CLEO_ZEMACH" ) {
    return EvtDalitzReso( dp, angPair, resPair, spinType, mass, width, EvtDalitzReso::RBW_CLEO_ZEMACH, FFp, FFr );
  } else if( shape=="GS" || shape=="GS_CLEO" ) {
    return EvtDalitzReso( dp, angPair, resPair, spinType, mass, width, EvtDalitzReso::GS_CLEO, FFp, FFr );
  } else if( shape=="GS_CLEO_ZEMACH" ) {
    return EvtDalitzReso( dp, angPair, resPair, spinType, mass, width, EvtDalitzReso::GS_CLEO_ZEMACH, FFp, FFr );
  } else if( shape=="GAUSS" || shape=="GAUSS_CLEO" ) {
    return EvtDalitzReso( dp, angPair, resPair, spinType, mass, width, EvtDalitzReso::GAUSS_CLEO, FFp, FFr );
  } else if( shape=="GAUSS_CLEO_ZEMACH" ) {
    return EvtDalitzReso( dp, angPair, resPair, spinType, mass, width, EvtDalitzReso::GAUSS_CLEO_ZEMACH, FFp, FFr );
  } else if( shape=="Flatte" ) {
    return EvtDalitzReso( dp, resPair, mass );
  } else if( shape=="LASS" ) {
    return EvtDalitzReso( dp, resPair, mass, width, aLass, rLass, BLass, phiBLass, RLass, phiRLass, cutoffLass, true );
  } else if( shape=="NonRes" ) {
    return EvtDalitzReso( );
  } else if( shape=="NonRes_Linear" ) {
    return EvtDalitzReso( dp, resPair, EvtDalitzReso::NON_RES_LIN );
  } else if( shape=="NonRes_Exp" ) {
    return EvtDalitzReso( dp, resPair, EvtDalitzReso::NON_RES_EXP, alpha );
  } else { //NBW
    if( shape!="NBW") report(Severity::Warning,"EvtGen")<<"EvtDalitzTable: shape "<<shape<<" is unknown. Defaulting to NBW."<<endl;
    return EvtDalitzReso( dp, angPair, resPair, spinType, mass, width, EvtDalitzReso::NBW, FFp, FFr );
  }
}

int EvtDalitzTable::getDaughterPairs(EvtId* resDaughter, EvtId* daughter, std::vector< std::pair<EvtCyclic3::Pair,EvtCyclic3::Pair> >& angAndResPairs) {
  int n(0);
  if(resDaughter[0]==daughter[0] && resDaughter[1]==daughter[1]) {
    angAndResPairs.push_back(std::make_pair(EvtCyclic3::BC,EvtCyclic3::AB)); n++;
  } else if(resDaughter[0]==daughter[1] && resDaughter[1]==daughter[0]) {
    angAndResPairs.push_back(std::make_pair(EvtCyclic3::CA,EvtCyclic3::AB)); n++;
  }
  
  if(resDaughter[0]==daughter[1] && resDaughter[1]==daughter[2]) {
    angAndResPairs.push_back(std::make_pair(EvtCyclic3::CA,EvtCyclic3::BC)); n++;
  } else if(resDaughter[0]==daughter[2] && resDaughter[1]==daughter[1]) {
    angAndResPairs.push_back(std::make_pair(EvtCyclic3::AB,EvtCyclic3::BC)); n++;
  }

  if(resDaughter[0]==daughter[2] && resDaughter[1]==daughter[0]) {
    angAndResPairs.push_back(std::make_pair(EvtCyclic3::AB,EvtCyclic3::CA)); n++;
  } else if(resDaughter[0]==daughter[0] && resDaughter[1]==daughter[2]) {
    angAndResPairs.push_back(std::make_pair(EvtCyclic3::BC,EvtCyclic3::CA)); n++;
  }

  return n;
}

double EvtDalitzTable::calcProbMax(EvtDalitzPlot dp, EvtDalitzDecayInfo* model) {

  double factor = 1.2; //factor to increase our final answer by
  int nStep(1000);      //number of steps - total points will be 3*nStep*nStep

  double maxProb(0);
  double min(0), max(0), step(0), min2(0), max2(0), step2(0);

  //first do AB, BC
  min = dp.qAbsMin(EvtCyclic3::AB);
  max = dp.qAbsMax(EvtCyclic3::AB);
  step = (max-min)/nStep;
  for(int i=0; i<nStep; ++i) {
    double qAB = min + i*step;
    min2 = dp.qMin(EvtCyclic3::BC,EvtCyclic3::AB,qAB);
    max2 = dp.qMax(EvtCyclic3::BC,EvtCyclic3::AB,qAB);
    step2 = (max2-min2)/nStep;
    for(int j=0; j<nStep; ++j) {
      double qBC = min2+ j*step2;
      EvtDalitzCoord coord(EvtCyclic3::AB,qAB,EvtCyclic3::BC,qBC);
      EvtDalitzPoint point(dp,coord);
      double prob = calcProb(point,model);
      if(prob > maxProb) maxProb = prob;
    }
  }

  //next do BC, CA
  min = dp.qAbsMin(EvtCyclic3::BC);
  max = dp.qAbsMax(EvtCyclic3::BC);
  step = (max-min)/nStep;
  for(int i=0; i<nStep; ++i) {
    double qBC = min + i*step;
    min2 = dp.qMin(EvtCyclic3::CA,EvtCyclic3::BC,qBC);
    max2 = dp.qMax(EvtCyclic3::CA,EvtCyclic3::BC,qBC);
    step2 = (max2-min2)/nStep;
    for(int j=0; j<nStep; ++j) {
      double qCA = min2+ j*step2;
      EvtDalitzCoord coord(EvtCyclic3::BC,qBC,EvtCyclic3::CA,qCA);
      EvtDalitzPoint point(dp,coord);
      double prob = calcProb(point,model);
      if(prob > maxProb) maxProb = prob;
    }
  }

  //finally do CA, AB
  min = dp.qAbsMin(EvtCyclic3::CA);
  max = dp.qAbsMax(EvtCyclic3::CA);
  step = (max-min)/nStep;
  for(int i=0; i<nStep; ++i) {
    double qCA = min + i*step;
    min2 = dp.qMin(EvtCyclic3::AB,EvtCyclic3::CA,qCA);
    max2 = dp.qMax(EvtCyclic3::AB,EvtCyclic3::CA,qCA);
    step2 = (max2-min2)/nStep;
    for(int j=0; j<nStep; ++j) {
      double qAB = min2+ j*step2;
      EvtDalitzCoord coord(EvtCyclic3::CA,qCA,EvtCyclic3::AB,qAB);
      EvtDalitzPoint point(dp,coord);
      double prob = calcProb(point,model);
      if(prob > maxProb) maxProb = prob;
    }
  }
  report(Severity::Info,"EvtGen") << "Largest probability found was " << maxProb << endl;
  report(Severity::Info,"EvtGen") << "Setting probMax to " << factor*maxProb << endl;
  return factor*maxProb;
}

double EvtDalitzTable::calcProb(EvtDalitzPoint point, EvtDalitzDecayInfo* model) {

  std::vector<std::pair<EvtComplex,EvtDalitzReso> > resonances = model->getResonances();

  EvtComplex amp(0,0);
  std::vector<std::pair<EvtComplex,EvtDalitzReso> >::iterator i = resonances.begin();
  for( ; i!= resonances.end(); i++) {
    std::pair<EvtComplex,EvtDalitzReso> res = (*i);
    amp += res.first * res.second.evaluate( point );
  }
  return abs2(amp);
}
