#ifndef ALIFASTJETINPUT_H
#define ALIFASTJETINPUT_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
  //---------------------------------------------------------------------
// Class for input particles
// manages the search for jets 
// Authors: Elena Bruna elena.bruna@yale.edu
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


#include<sstream>  // needed for internal io
#include <vector> 
#include <cmath> 
#endif
class AliJetFinder;
class AliJetHeader;
class AliJetReader;


class AliFastJetInput : public TObject
{
 public:
    AliFastJetInput();
    AliFastJetInput(const AliFastJetInput &input);
    AliFastJetInput& operator=(const AliFastJetInput& source);
    virtual ~AliFastJetInput() {;}
    void SetHeader(AliJetHeader *header)  {fHeader=header;}
    void SetReader(AliJetReader *reader)  {fReader=reader;}
    void FillInput();
    vector<fastjet::PseudoJet> GetInputParticles()   const {return fInputParticles;}
    vector<fastjet::PseudoJet> GetInputParticlesCh() const {return fInputParticlesCh;}
    Float_t  EtaToTheta(Float_t arg);
    static Double_t Thermalspectrum(const Double_t *x, const Double_t *par);

 private:
   AliJetReader *fReader;  //! reader 
   AliJetHeader *fHeader;  //! header 
   
    vector<fastjet::PseudoJet> fInputParticles;     //! input particles for FastJet
    vector<fastjet::PseudoJet> fInputParticlesCh;   //! input charged particles for FastJet

  ClassDef(AliFastJetInput, 1); // Analysis task for standard jet analysis
};
 
#endif
