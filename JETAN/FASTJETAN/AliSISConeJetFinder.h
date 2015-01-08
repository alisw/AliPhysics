#ifndef ALISISCONEJETFINDER_H
#define ALISISCONEJETFINDER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//---------------------------------------------------------------------
// SISCone (FastJet v2.3.4) finder algorithm interface
//
// Author: swensy.jangal@ires.in2p3.fr
//  
// ** 2011 magali.estienne@subatech.in2p3.fr &  alexandre.shabetai@cern.ch
// Modified accordingly to reader/finder splitting and new handling of neutral information (via FastJetInput)
//---------------------------------------------------------------------

#include "AliJetFinder.h"

class AliFastJetHeaderV1;
class AliFastJetInput;
class AliFastJetBkg;
using namespace std;

class AliSISConeJetFinder : public AliJetFinder
{
 public:
  AliSISConeJetFinder();
  ~AliSISConeJetFinder();

  void    FindJets(); 

  // others
  Bool_t  ProcessEvent(); 
  void    WriteJHeaderToFile() const;

  protected:
  AliSISConeJetFinder(const AliSISConeJetFinder& rfj);
  AliSISConeJetFinder& operator = (const AliSISConeJetFinder& rsfj);

  AliFastJetInput*  fInputFJ;     //! input particles array
  AliFastJetBkg*    fJetBkg;      //! pointer to bkg class

  ClassDef(AliSISConeJetFinder,3) // SISCONE analysis class

};

#endif
