#ifndef ALISISCONEJETFINDER_H
#define ALISISCONEJETFINDER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */



//---------------------------------------------------------------------
// SISCone (FastJet v2.3.4) finder algorithm interface
//
// Author: swensy.jangal@ires.in2p3.fr
//  
//---------------------------------------------------------------------

// FastJet classes 
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"

// Get info on how fastjet was configured
#include "fastjet/config.h"
#ifdef ENABLE_PLUGIN_SISCONE
#include "fastjet/SISConePlugin.hh"
#endif

#include<sstream>  // needed for internal io
#include <vector> 
#include <cmath> 

#include "AliFastJetHeaderV1.h"
#include "AliJetFinder.h"

using namespace std;

class AliSISConeJetFinder : public AliJetFinder
{
 public:

  AliSISConeJetFinder();
  ~AliSISConeJetFinder();

  void    FindJets(); 

  // others
 
  void    WriteJHeaderToFile() const;
  Float_t EtaToTheta(Float_t arg);
  void    InitTask(TChain* tree);

  protected:
  AliSISConeJetFinder(const AliSISConeJetFinder& rfj);
  AliSISConeJetFinder& operator = (const AliSISConeJetFinder& rsfj);

  ClassDef(AliSISConeJetFinder,2)
};

#endif
