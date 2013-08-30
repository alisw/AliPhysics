#ifndef ALIFASTJETFINDER_H
#define ALIFASTJETFINDER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//---------------------------------------------------------------------
// FastJet v2.3.4 finder algorithm interface
//
// Author: Rafael.Diaz.Valdes@cern.ch
//
// ** 2011 magali.estienne@subatech.in2p3.fr &  alexandre.shabetai@cern.ch
// new implementation of background subtraction
// allowing to subtract bkg using a different algo than the one used for signal jets  
//---------------------------------------------------------------------

// FastJet classes 
#ifndef __CINT__
# include "fastjet/PseudoJet.hh"
# include "fastjet/ClusterSequenceArea.hh"
# include "fastjet/AreaDefinition.hh"
# include "fastjet/JetDefinition.hh"
#else
namespace fastjet {
  class PseudoJet;
  class ClusterSequenceArea;
  class AreaDefinition;
  class JetDefinition;
}
#endif

#include "AliJetFinder.h"

class AliFastJetInput;
class AliFastJetBkg;

using namespace std;

class AliFastJetFinder : public AliJetFinder
{
 public:

  AliFastJetFinder();
  ~AliFastJetFinder();

  virtual void      FindJets(); 
  void              RunTest(const char* datafile); // a simple test
  virtual void      WriteJHeaderToFile() const;
  virtual Bool_t    ProcessEvent();
      
  protected:
  AliFastJetFinder(const AliFastJetFinder& rfj);
  AliFastJetFinder& operator = (const AliFastJetFinder& rsfj);
  AliFastJetInput*  fInputFJ;  //! input particles array
  AliFastJetBkg*    fJetBkg;   //! pointer to bkg class

  ClassDef(AliFastJetFinder,3) //  Fastjet analysis class

};

#endif
