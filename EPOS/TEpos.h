/*
 *###################################################################
 *#        EPOS 1.67     K. WERNER, T. PIEROG, S. PORTEBOEUF.       #
 *#                      Contact: werner@subatech.in2p3.fr          #
 *###################################################################
 *
 * TEpos.h
 * 
 * Wraper class for interfacing EPOS model, derived from ROOT's TGenerator.
 * It generates temporary input file for the model, providing user with
 * ability to add his/hers own lines to the input.
 * Output is read directly from common blocks.
 *
 *      Author: Piotr Ostrowski, postrow@if.pw.edu.pl
 */


#ifndef TEPOS_H
#define TEPOS_H

#ifndef ROOT_TGenerator
#include "TGenerator.h"
#endif

#include <vector>
class TObjArray;

class TEpos : public TGenerator {


public:

   TEpos();

   virtual            ~TEpos();

   virtual void        Initialize();

   virtual void        GenerateEvent();

   virtual Int_t       ImportParticles(TClonesArray *particles, Option_t *option="");
   virtual TObjArray*  ImportParticles(Option_t *option="");

   //Geters for event header

   Float_t GetBimevt(); //   bimevt ........ absolute value of impact parameter
   Float_t GetPhievt(); //   phievt ........ angle of impact parameter
   Int_t GetKolevt();   //   kolevt ........ number of collisions
   Int_t GetKoievt();   //   koievt ........ number of inelastic collisions
   Float_t GetPmxevt(); //   pmxevt ........ reference momentum
   Float_t GetEgyevt(); //   egyevt ........ pp cm energy (hadron) or string energy (lepton)
   Int_t GetNpjevt();   //   npjevt ........ number of primary projectile participants
   Int_t GetNtgevt();   //   ntgevt ........ number of primary target participants
   Int_t GetNpnevt();   //   npnevt ........ number of primary projectile neutron spectators
   Int_t GetNppevt();   //   nppevt ........ number of primary projectile proton spectators
   Int_t GetNtnevt();   //   ntnevt ........ number of primary target neutron spectators
   Int_t GetNtpevt();   //   ntpevt ........ number of primary target proton spectators
   Int_t GetJpnevt();   //   jpnevt ........ number of absolute projectile neutron spectators
   Int_t GetJppevt();   //   jppevt ........ number of absolute projectile proton spectators
   Int_t GetJtnevt();   //   jtnevt ........ number of absolute target neutron spectators
   Int_t GetJtpevt();   //   jtpevt ........ number of absolute target proton spectators
   Float_t GetXbjevt(); //   xbjevt ........ bjorken x for dis
   Float_t GetQsqevt(); //   qsqevt ........ q**2 for dis
   Int_t GetNglevt();   //   nglevt ........ number of collisions acc to  Glauber
   Float_t GetZppevt(); //   zppevt ........ average Z-parton-proj
   Float_t GetZptevt(); //   zptevt ........ average Z-parton-targ

   //Parameters for the generator:
   void SetLaproj(Int_t laproj) { fLaproj = laproj; }
   void SetMaproj(Int_t maproj) { fMaproj = maproj; }
   void SetLatarg(Int_t latarg) { fLatarg = latarg; }
   void SetMatarg(Int_t matarg) { fMatarg = matarg; }

   void SetBminim(Float_t bminim) { fBminim = bminim; }
   void SetBmaxim(Float_t bmaxim) { fBmaxim = bmaxim; }

   void SetPhimin(Float_t phimin) { fPhimin = phimin; }
   void SetPhimax(Float_t phimax) { fPhimax = phimax; }

   void SetEcms(Float_t ecms) { fEcms = ecms; }
   void SetSplitting(Bool_t splitting) { fSplitting = splitting; }

   void AddNoDecay(Int_t nodecay);
   void AddExtraInputLine(const char *);

   Int_t GetLaproj() const { return fLaproj; }
   Int_t GetMaproj() const { return fMaproj; }
   Int_t GetLatarg() const { return fLatarg; }
   Int_t GetMatarg() const { return fMatarg; }

   Double_t GetBminim() const { return fBminim; }
   Double_t GetBmaxim() const { return fBmaxim; }

   Double_t GetPhimin() const { return fPhimin; }
   Double_t GetPhimax() const { return fPhimax; }

   Double_t GetEcms() const { return fEcms; }
   Bool_t GetSplitting() const { return fSplitting; }

protected:
   virtual void GenerateInputFile();
   virtual const char * GetInputFileName() const { return "/tmp/epos.input"; }

   Int_t fLaproj;	//atomic number of projectile
   Int_t fMaproj;	//mass number of projectile
   Int_t fLatarg;	//atomic number of target
   Int_t fMatarg;	//mass number of target
   Float_t fBminim;     //lower limit for impact parameter
   Float_t fBmaxim;     //upper limit for impact parameter
   Float_t fPhimin;     //lower limit for reaction plane angle
   Float_t fPhimax;     //upper limit for reaction plane angle
   Float_t fEcms;       //energy in CMS 
   Bool_t fSplitting;   //consideration of parton splitting and fussion

   std::vector<Int_t> fNoDecays;		// list of user supplied ISAJET codes of particles that will not be decayed by the model
   std::vector<const char * > fExtraInputLines; // list of user supplied lines to be added to EPOS input


private:
   ClassDef(TEpos,1)  //Interface to EPOS Event Generator
};

#endif

