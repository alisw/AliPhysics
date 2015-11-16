#ifndef ALIGENEPOSREADER_H
#define ALIGENEPOSREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Realisation of the AliGenReader interface to be used with AliGenExFile.
// NextEvent() loops over events
// and NextParticle() loops over particles.
// This implementation reads EPOS v3.111 output format (from ROOT trees)
// Author: Igor Lakomov <Igor.Lakomov@cern.ch>

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <map>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TParticle.h>
#include <TString.h>
#include <TDatabasePDG.h>

#include "AliGenReader.h"
#include "AliGenEventHeader.h"
//#include "AliRun.h"
#include "AliStack.h" 
#include "AliLog.h"

class AliGenEposReader : public AliGenReader
{
 public:
  AliGenEposReader();

  AliGenEposReader(const AliGenEposReader &reader);
  virtual ~AliGenEposReader();
  AliGenEventHeader* GetGenEventHeader() const {return fGenEventHeader;};
  virtual void Init();
  virtual Int_t NextEvent();
  virtual TParticle* NextParticle();
  virtual void RewindEvent();
  virtual void ChangeFile(const Text_t *fNewFileName);
  AliGenEposReader & operator=(const AliGenEposReader &rhs);

 protected:
  Int_t fNcurrent;
  Int_t fNparticle; 
  Int_t fCurrentEvent;
  Int_t fCurrentParticle;
  TTree *fTreeNtuple;
  TTree *fTreeHeader;
  TFile *fFile;
  AliGenEventHeader* fGenEventHeader;   // AliGenEventHeader
  std::multimap<Int_t,Int_t> fMothersMap;

  // Parameters in EPOS Tree:
  //--------------teposevent---------------
  Int_t np;       //number of particles
  Float_t bim;    //impact parameter
  std::vector<Float_t> zus; //different meaning depending on ptl type:
     		            // partons: presently unused
                  	    // hadrons:  decay information :
	                    // -999 : hadron is decay product from decay 
        	            //        in cascade part (mother unknown)
                  	    //   -1 : hadron is decay product, mother not stored
	                    //   >0 : hadron is decay product, mother index = zus
        	            //   so, zus is the number n of the particle, not the id
           	            //   to get the id you need id(n)
       		            //   -2 : no mother
                  	    // phi 331, K+ 130
  std::vector<Float_t> px;  //px
  std::vector<Float_t> py;  //py   particle four momentum
  std::vector<Float_t> pz;  //pz
  std::vector<Float_t> e;   //energy pf particle
  std::vector<Float_t> x;   //x component of formation point
  std::vector<Float_t> y;   //y component of formation point
  std::vector<Float_t> z;   //z component of formation point
  std::vector<Float_t> t;   //formation time
  std::vector<Int_t> id;    //particle id
  std::vector<Int_t> ist;   //particle status (hadron last generation(0) or not(1))
  std::vector<Int_t> ity;   //type of particle origin (20-29 from soft strings, 30-39 from hard strings, 40-59 from remnants, 60 from fluid)
  std::vector<Int_t> ior;   //index of father (resonance decay products)
  std::vector<Int_t> jor;   //index of mother (mothers are needed for exemple for strings: the partons between ior and jor constitute the string)

//  --------------teposhead---------------
  Int_t fIversn;  //EPOS version number
  Int_t fLaproj;  //atomic number projectile
  Int_t fMaproj;  //mass number projectile
  Int_t fLatarg;  //atomic number target
  Int_t fMatarg;  //mass number target
  Float_t fEngy;  //energy in the CMS in GeV
  Int_t fNfull;   //number of full events
  Int_t fNfreeze; //number of freeze outs per full event

 private:
  Int_t EposToPdg(Int_t code);
  void Copy(TObject&) const;
  ClassDef(AliGenEposReader,3)
};
#endif
