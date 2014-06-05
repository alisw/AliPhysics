//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package. If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2011      University of Warwick, UK
//
// Module: EvtExternalGenFactory
//
// Description: A factory type method to create engines for external physics
// generators like Pythia.
//
// Modification history:
//
//    John Back       April 2011            Module created
//
//------------------------------------------------------------------------------
//

#ifndef EVTEXTERNALGENFACTORY_HH
#define EVTEXTERNALGENFACTORY_HH

#include "EvtGenModels/EvtAbsExternalGen.hh"

#include <map>

class EvtExternalGenFactory {

public:
  
  enum genId {PythiaGenId = 0, PhotosGenId, TauolaGenId};

  static EvtExternalGenFactory* getInstance();

  EvtAbsExternalGen* getGenerator(int genId = 0);

  void initialiseAllGenerators();

  void definePythiaGenerator(std::string xmlDir, bool convertPhysCodes, bool useEvtGenRandom = true);
  void definePhotosGenerator(std::string photonType = "gamma", bool useEvtGenRandom = true);
  void defineTauolaGenerator(bool useEvtGenRandom = true);

  //methods to add configuration commands to the pythia generators
  //void addPythiaCommand( std::string generator, std::string module, std::string param, std::string value);
  //void addPythia6Command(std::string generator, std::string module, std::string param, std::string value);

protected:

  EvtExternalGenFactory();
  ~EvtExternalGenFactory();

  typedef std::map<int, EvtAbsExternalGen*> ExtGenMap;
  typedef std::map<int, std::map<std::string, std::vector<std::string> > > ExtGenCommandMap;

private:

  EvtExternalGenFactory(const EvtExternalGenFactory&) {};

  ExtGenMap _extGenMap;
  ExtGenCommandMap _extGenCommandMap;

};

#endif
