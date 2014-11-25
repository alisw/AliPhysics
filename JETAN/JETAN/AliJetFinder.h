#ifndef ALIJETFINDER_H
#define ALIJETFINDER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
//---------------------------------------------------------------------
// Jet finder base class
// manages the search for jets 
// Authors: jgcn@mda.cinvestav.mx
//          andreas.morsch@cern.ch
//          magali.estienne@subatech.in2p3.fr
//          alexandre.shabetai@cern.ch
//---------------------------------------------------------------------

#include "AliJetCalTrk.h"
#include "AliAODJet.h"
#include "AliJetHeader.h"
#include "AliAODJetEventBackground.h"
class AliAODEvent;

class AliJetFinder : public TObject 
{
 public:
  AliJetFinder();
  virtual ~AliJetFinder();

  // Getters
  virtual AliJetCalTrkEvent* GetCalTrkEvent() const {return fCalTrkEvent;}
  virtual AliJetHeader *GetJetHeader() const {return fHeader;}
  virtual AliAODJetEventBackground* GetEventBackground() const {return fAODEvBkg;}
  // Setters
  virtual void              SetCalTrkEvent(AliJetCalTrkEvent& event) {fCalTrkEvent = &event;}
  virtual void              SetJetHeader(AliJetHeader* h) {fHeader=h;}
  virtual void              SetEventBackground(AliAODJetEventBackground* bkg) {fAODEvBkg = bkg;}
  // Others
  virtual void              AddJet(AliAODJet jet);
  virtual void              WriteHeaderToFile();
  virtual void		    WriteHeader();
  // the following have to be implemented for each specific finder
  virtual void              Init() {}
  virtual void              Reset() {fNAODjets = 0;}
  virtual void              FindJets() {}
  virtual void              ComputeBkgs() {}
  virtual void              CreateOutputObjects(TList * const /*histos*/) {} // Used by CDF for histo storage

  // some methods to allow steering from the outside
  virtual Bool_t            ProcessEvent();
  virtual void              ConnectAOD(const AliAODEvent* aod);
  virtual void              ConnectAODNonStd(AliAODEvent* aod,const char* bname);
  virtual void              AddHistosToList(TList */*list*/) {}

 protected:
  AliJetHeader*             fHeader;         //  pointer to header
  TClonesArray*             fAODjets;        //! reconstructed jets
  Int_t                     fNAODjets;       //  number of reconstructed jets
  AliAODJetEventBackground* fAODEvBkg;       //! bkg object to be store
  Int_t                     fDebug;          //  debug option, set through the header
  AliJetCalTrkEvent*        fCalTrkEvent;    //  pointer to AliJetCalTrkEvent object

 private:
  AliJetFinder(const AliJetFinder& rJetFinder); // not implemented
  AliJetFinder& operator = (const AliJetFinder& rhsf); // not implemented
 
  ClassDef(AliJetFinder,3)                   //  base class for any jet finder

};

#endif
