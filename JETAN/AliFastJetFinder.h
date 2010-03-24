#ifndef ALIFASTJETFINDER_H
#define ALIFASTJETFINDER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */



//---------------------------------------------------------------------
// FastJet v2.3.4 finder algorithm interface
//
// Author: Rafael.Diaz.Valdes@cern.ch
//  
//---------------------------------------------------------------------

//FastJet classes 
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/JetDefinition.hh"
// get info on how fastjet was configured
#include "fastjet/config.h"
#ifdef ENABLE_PLUGIN_SISCONE
#include "fastjet/SISConePlugin.hh"
#endif


#include<sstream>  // needed for internal io
#include <vector> 
#include <cmath> 

#include "AliJetFinder.h"
#include "AliFastJetHeaderV1.h"

using namespace std;
class AliFastJetInput;
class AliJetBkg;

class AliFastJetFinder : public AliJetFinder
{
 public:

  AliFastJetFinder();
  ~AliFastJetFinder();

  virtual void    FindJets(); 
  void    RunTest(const char* datafile); // a simple test
  virtual void    WriteJHeaderToFile() const;
  Float_t EtaToTheta(Float_t arg);
  void    InitTask(TChain* tree);
  virtual Bool_t ProcessEvent();
  virtual Bool_t  ProcessEvent2();
  
      
  protected:
  AliFastJetFinder(const AliFastJetFinder& rfj);
  AliFastJetFinder& operator = (const AliFastJetFinder& rsfj);
  AliFastJetInput*                fInputFJ;     //! input particles array
  AliJetBkg*                      fJetBkg;      //! pointer to bkg class
  ClassDef(AliFastJetFinder,2)
};

#endif
