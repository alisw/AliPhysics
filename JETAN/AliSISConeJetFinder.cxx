/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
 

//---------------------------------------------------------------------
// SISCone (FastJet v2.3.4) finder algorithm interface
//
// Author: swensy.jangal@ires.in2p3.fr 
//  
//---------------------------------------------------------------------

#include <Riostream.h>
#include <TArrayF.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TRandom.h>

#include "AliHeader.h"
#include "AliJet.h"
#include "AliJetKineReader.h"
#include "AliJetReader.h"
#include "AliJetReaderHeader.h"
#include "AliSISConeJetFinder.h"
#include "AliSISConeJetHeader.h"

#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"

// Get info on how fastjet was configured
#include "fastjet/config.h"

#ifdef ENABLE_PLUGIN_SISCONE
#include "fastjet/SISConePlugin.hh"
#endif

#include<sstream>  // Needed for internal io
#include<vector> 
#include <cmath> 

using namespace std;

ClassImp(AliSISConeJetFinder)

//____________________________________________________________________________

AliSISConeJetFinder::AliSISConeJetFinder():
AliJetFinder()
{
  // Constructor
}

//____________________________________________________________________________

AliSISConeJetFinder::~AliSISConeJetFinder()
{
  // destructor
}

//______________________________________________________________________________
void AliSISConeJetFinder::FindJets()
{

  Bool_t debug = kFALSE;
  
  // Pick up siscone header
  AliSISConeJetHeader *header = (AliSISConeJetHeader*)fHeader;

  // Check if we are reading AOD jets
  TRefArray *refs = 0;
  Bool_t fromAod = !strcmp(fReader->ClassName(),"AliJetAODReader");
  if (fromAod) { refs = fReader->GetReferences(); }
  

//******************************** SISCONE PLUGIN CONFIGURATION
// Here we look for SISCone parameters in the header and we define our plugin.  

  Double_t coneRadius                = header->GetConeRadius();                 // cone radius
  Double_t overlapThreshold          = header->GetOverlapThreshold();           // overlap parameter
  Int_t    nPassMax                  = header->GetNPassMax();                   // maximum number of passes
  Double_t ptProtoJetMin             = header->GetPtProtojetMin();              // pT min of protojets
  Double_t caching                   = header->GetCaching();                    // do we record found cones for this set of data?

  if (header->GetSplitMergeScale() == 0) fastjet::SISConePlugin::SplitMergeScale splitMergeScale = fastjet::SISConePlugin::SM_pttilde; // There's only one split merge scale

  fastjet::JetDefinition::Plugin * plugin;
  plugin = new fastjet::SISConePlugin(coneRadius, overlapThreshold, nPassMax, ptProtoJetMin, caching);

//********************************  READING OF INPUT PARTICLES
// Here we look for px, py pz and energy of each particle that we gather in a PseudoJet object, and we put all these PseudoJet in a vector of PseudoJets : input_particles. 

  TClonesArray *lvArray = fReader->GetMomentumArray();
  Int_t nIn = lvArray->GetEntries();

  // We check if lvArray is ok
  if(lvArray == 0)
  {
    cout << "Could not get the momentum array" << endl;
    return;
  }

  if(nIn == 0)// nIn = Number of particles in the event
  {
    if (debug) cout << "entries = 0 ; Event empty !!!" << endl ;
    return;
  }

  Int_t nJets = 0;        // Number of jets in this event
  fJets->SetNinput(nIn) ; // fJets = AliJet number of input objects
  Float_t px,py,pz,en;
  vector<fastjet::PseudoJet> inputParticles; 

  // Load input vectors
  for(Int_t i = 0; i < nIn; i++)
  { 
    TLorentzVector *lv = (TLorentzVector*) lvArray->At(i);
    px = lv->Px();
    py = lv->Py();
    pz = lv->Pz();
    en = lv->Energy();
      
    fastjet::PseudoJet inputPart(px,py,pz,en); 
    inputPart.set_user_index(i);
    inputParticles.push_back(inputPart); 

  }   
  
//******************************** JETS FINDING 

  fastjet::ClusterSequence clustSeq(inputParticles, plugin);

//***************************** JETS EXTRACTION AND CORRECTION

  // Here we extract inclusive jets with pt > ptmin, sorted by pt 
  Double_t ptMin = header->GetMinJetPt(); 
  vector<fastjet::PseudoJet> inclusiveJets = clustSeq.inclusive_jets(ptMin);
  vector<fastjet::PseudoJet> jets = sorted_by_pt(inclusiveJets);

  for (Int_t j = 0; j < jets.size(); j++)
  {
    cout<<"********************************** Reconstructed jet(s) (non corrected)"<<endl;
    cout<<"Jet number "<<j+1<<" : "<<"Rapidity : "<<jets[j].rap()<<" Phi : "<<jets[j].phi()<<" pT : "<<jets[j].perp()<<endl;
    cout<<"px = "<<jets[j].px()<<endl;
    cout<<"py = "<<jets[j].py()<<endl;
    cout<<"pz = "<<jets[j].pz()<<endl;
    cout<<"e = "<<jets[j].E()<<endl;
    cout<<"******************"<<endl;
    cout<<"******************"<<endl;
    cout<<"******************"<<endl;

    // Go to write AOD info
    AliAODJet aodjet (jets[j].px(), jets[j].py(), jets[j].pz(), jets[j].E());
    if(debug) aodjet.Print("");
    AddJet(aodjet);
  }
}
 
//____________________________________________________________________________

void AliSISConeJetFinder::WriteJHeaderToFile()
{
  fHeader->Write();
}
