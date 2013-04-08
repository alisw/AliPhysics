//-*- Mode: C++ -*-
// $Id$

//* This file is property of and copyright by the ALICE Project        * 
//* ALICE Experiment at CERN, All rights reserved.                     *
//* See cxx source for full Copyright notice                           *

/// @file   AliDxHFEToolsMC.h
/// @author Hege Erdal, Matthias Richter
/// @date   2012-07-19
/// @brief  Common Tools for MC particle selection
///

#ifndef ALIDXHFETOOLSMC_H
#define ALIDXHFETOOLSMC_H

#include "TObject.h"
#include <vector>

class AliVEvent;
class AliVParticle;
class TH1;

using std::vector;

/**
 * @class AliDxHFEToolsMC
 * Common Tools for MC particle selection.
 */
class AliDxHFEToolsMC {
  public:
  /// constructor
  AliDxHFEToolsMC(const char* options="");
  /// destructor
  virtual ~AliDxHFEToolsMC();

  // different pdgs
  enum{
    kPDGnone=0,
    kPDGd=1,
    kPDGu=2,
    kPDGs=3,
    kPDGc=4,
    kPDGb=5,
    kPDGelectron=11,
    kPDGmuon=13,
    kPDGgluon=21,
    kPDGgamma=22,
    kPDGpi0=111,
    kPDGpion=211,
    kPDGeta=221,
    kPDGkaon=321,
    kPDGD0=421,
    kPDGJpsi=443,
    kPDGproton=2212
  };

  enum {
    kPDGLabelPositron,
    kPDGLabelElectron,
    kPDGLabelMuPlus,
    kPDGLabelMuMinus,
    kPDGLabelPiPlus,
    kPDGLabelPiMinus,
    kPDGLabelKPlus,
    kPDGLabelKMinus,
    kPDGLabelProton,
    kPDGLabelAntiproton,
    kPDGLabelOthers,
    kNofPDGLabels
  };

  enum {
    kPDGMotherLabelD,
    kPDGMotherLabelU,
    kPDGMotherLabelS,
    kPDGMotherLabelC,
    kPDGMotherLabelB,
    kPDGMotherLabelGluon,
    kPDGMotherLabelGamma,
    kPDGMotherLabelPi0,
    kPDGMotherLabelEta,
    kPDGMotherLabelProton,
    kPDGMotherLabelOthers,
    kNofPDGMotherLabels
  };

  enum {
    kMCFirst = 0,
    kMCLast
  };

  enum {
    kOriginNone=0,
    kOriginDown,
    kOriginUp,
    kOriginStrange,
    kOriginCharm,
    kOriginBeauty,
    kOriginGluon,
    kOriginGluonCharm,
    kOriginGluonBeauty,
    kNrOrginMother
  };

  enum {
    kGetOriginMother=0,
    kGetFirstMother
  };

  /// initialize according to options
  int Init(const char* /*option*/);

  /// init MC info from event object
  int InitMCParticles(const AliVEvent* pEvent);

  /// Returning MC Array
  TObjArray* GetMCArray() const {return fMCParticles;} 
  
  /// check state
  bool IsInitialized() const {return fMCParticles!=NULL;}

  /// clear internal memory
  virtual void Clear(const char* option="");

  /// flag indicating to do the selection on MC first
  bool MCFirst() const {return fSequence==kMCFirst;}
  /// flag indicating to do the selection on MC last
  bool MCLast() const {return fSequence==kMCLast;}

  /// Return the result on check on initial quark
  int GetOriginMother() const {return fOriginMother;}

  /// check if pdg should be rejected
  /// always false if pdg list is not initialized
  bool RejectByPDG(AliVParticle* p, bool doStatistics=true, int* pdgParticleResult=NULL);
  bool RejectByPDG(AliVParticle* p, int* pdgParticleResult) {
    return RejectByPDG(p, true, pdgParticleResult);
  }

  /// check if pdg should be rejected by mother
  /// always false if mother pdg list is not initialized
  bool RejectByMotherPDG(AliVParticle* p, bool doStatistics=true);

  /// Finds pdg of first or origin mother and returns value
  int FindMotherPDG(AliVParticle* p, bool bReturnFirstMother=false);

  /// step through tree and find the original mother particle
  /// TODO: this can possibly be const, however, a member is set inside
  /// check whether this is really necessary
  int FindPdgOriginMother(AliVParticle* p,bool bReturnFirstMother=false);

  // Compare pdg to quark and gluon
  void CheckOriginMother(int pdg);
  // Tests if particle have been marked as HF quark 
  Bool_t TestIfHFquark(int origin);

  //Tests of pdg corresponds to HF meson
  Bool_t TestMotherHFMeson(int pdg);

  // Setting MC label from outside
  void SetMClabel(int mclab){fMClabel=mclab;}

  /// TODO: want to have this function to be private again, currently
  /// used with an external vector
  /// check if pdg should be rejected, particle is not rejected
  /// if it is in the list, returns always false if list is empty
  bool RejectByPDG(int pdg, const vector<int> &list) const;

  int GetNrMCParticles() const {return fNrMCParticles;}
  int CheckMCParticle(AliVParticle* p);

  /// mapping of pdg code to enum
  int MapPDGLabel(int pdg) const;
  /// mapping of pdg code to enum
  int MapPDGMotherLabel(int pdg) const;

 protected:

 private:
  /// copy contructor prohibited
  AliDxHFEToolsMC(const AliDxHFEToolsMC&);
  /// assignment operator prohibited
  AliDxHFEToolsMC& operator=(const AliDxHFEToolsMC&);

  /// create control histogram
  TH1* CreateControlHistogram(const char* name,
			      const char* title,
			      int nBins,
			      const char** binLabels) const;

  static const char* fgkPDGBinLabels[];
  static const char* fgkPDGMotherBinLabels[];
  static const char* fgkStatisticsBinLabels[];

  int fSequence;           //  sequence of checks
  TObjArray* fMCParticles; //! pointer to external array of MC particles
  vector<int> fPDGs;       //  PDGs to be selected
  vector<int> fMotherPDGs; //  mother PDGs to be selected
  TH1* fHistPDG;           //  control histogram pdg of selected particle
  TH1* fHistPDGMother;     //  control histogram pdg of selected particle
  int fOriginMother;       //  Holds the origin motherquark (process)
  int fMClabel;            //  MClabel passed from outside (default =-1)
  int fNrMCParticles;      //  number of MC particles 
  bool fUseKine;           //  For looping over stack directly

  ClassDef(AliDxHFEToolsMC, 3);
};
#endif
