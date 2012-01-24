#ifndef ALIFASTJETINPUT_H
#define ALIFASTJETINPUT_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//---------------------------------------------------------------------
// Class for input particles
// manages the search for jets 
// Authors: Elena Bruna elena.bruna@yale.edu
//
// ** 2011 magali.estienne@subatech.in2p3.fr &  alexandre.shabetai@cern.ch             
// Modified accordingly to reader/finder splitting and new handling of neutral information
//---------------------------------------------------------------------

#include <TObject.h>

// FastJet classes 
#include "fastjet/PseudoJet.hh"

#include <vector> 

class AliJetCalTrkEvent;
class AliJetHeader;

class AliFastJetInput : public TObject
{
 public:
  AliFastJetInput();
  AliFastJetInput(const AliFastJetInput &input);
  AliFastJetInput& operator=(const AliFastJetInput& source);
  virtual                    ~AliFastJetInput() {;}
  void                       SetHeader(AliJetHeader *header)            {fHeader=header;}
  void                       SetCalTrkEvent(AliJetCalTrkEvent *caltrk)  {fCalTrkEvent=caltrk;}
  void                       FillInput();
  vector<fastjet::PseudoJet> GetInputParticles()   const                {return fInputParticles;}
  vector<fastjet::PseudoJet> GetInputParticlesCh() const                {return fInputParticlesCh;}
  static Double_t            Thermalspectrum(const Double_t *x, const Double_t *par);

 private:
  AliJetHeader *fHeader;                        //! header 
  AliJetCalTrkEvent *fCalTrkEvent;              //! caltrkevent
   
  vector<fastjet::PseudoJet> fInputParticles;   //! input particles for FastJet
  vector<fastjet::PseudoJet> fInputParticlesCh; //! input charged particles for FastJet

  ClassDef(AliFastJetInput, 2)                  //  fills input particles for FASTJET based analysis
    
};
 
#endif
