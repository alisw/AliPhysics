#ifndef ALITPYTHIA8_H
#define ALITPYTHIA8_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// TPythia8                                                                   //
//                                                                            //
// TPythia is an interface class to C++ version of Pythia 8.1                 //
// event generators, written by T.Sjostrand.                                  //
////////////////////////////////////////////////////////////////////////////////
/*
*------------------------------------------------------------------------------------*
 |                                                                                    |
 |  *------------------------------------------------------------------------------*  |
 |  |                                                                              |  |
 |  |                                                                              |  |
 |  |   PPP   Y   Y  TTTTT  H   H  III    A      Welcome to the Lund Monte Carlo!  |  |
 |  |   P  P   Y Y     T    H   H   I    A A     This is PYTHIA version 8.100      |  |
 |  |   PPP     Y      T    HHHHH   I   AAAAA    Last date of change: 20 Oct 2007  |  |
 |  |   P       Y      T    H   H   I   A   A                                      |  |
 |  |   P       Y      T    H   H  III  A   A    Now is 27 Oct 2007 at 18:26:53    |  |
 |  |                                                                              |  |
 |  |   Main author: Torbjorn Sjostrand; CERN/PH, CH-1211 Geneva, Switzerland,     |  |
 |  |     and Department of Theoretical Physics, Lund University, Lund, Sweden;    |  |
 |  |     phone: + 41 - 22 - 767 82 27; e-mail: torbjorn@thep.lu.se                |  |
 |  |   Author: Stephen Mrenna; Computing Division, Simulations Group,             |  |
 |  |     Fermi National Accelerator Laboratory, MS 234, Batavia, IL 60510, USA;   |  |
 |  |     phone: + 1 - 630 - 840 - 2556; e-mail: mrenna@fnal.gov                   |  |
 |  |   Author: Peter Skands; CERN/PH, CH-1211 Geneva, Switzerland,                |  |
 |  |     and Theoretical Physics Department,                                      |  |
 |  |     Fermi National Accelerator Laboratory, MS 106, Batavia, IL 60510, USA;   |  |
 |  |     phone: + 41 - 22 - 767 24 59; e-mail: skands@fnal.gov                    |  |
 |  |                                                                              |  |
 |  |   The main program reference is the 'Brief Introduction to PYTHIA 8.1',      |  |
 |  |   T. Sjostrand, S. Mrenna and P. Skands, arXiv:0710.3820                     |  |
 |  |                                                                              |  |
 |  |   The main physics reference is the 'PYTHIA 6.4 Physics and Manual',         |  |
 |  |   T. Sjostrand, S. Mrenna and P. Skands, JHEP05 (2006) 026 [hep-ph/0603175]. |  |
 |  |                                                                              |  |
 |  |   An archive of program versions and documentation is found on the web:      |  |
 |  |   http://www.thep.lu.se/~torbjorn/Pythia.html                                |  |
 |  |                                                                              |  |
 |  |   This program is released under the GNU General Public Licence version 2.   |  |
 |  |   Please respect the MCnet Guidelines for Event Generator Authors and Users. |  |
 |  |                                                                              |  |
 |  |   Disclaimer: this program comes without any guarantees.                     |  |
 |  |   Beware of errors and use common sense when interpreting results.           |  |
 |  |                                                                              |  |
 |  |   Copyright (C) 2007 Torbjorn Sjostrand                                      |  |
 |  |                                                                              |  |
 |  |                                                                              |  |
 |  *------------------------------------------------------------------------------*  |
 |                                                                                    |
 *------------------------------------------------------------------------------------*
*/

#include "TGenerator.h"
#include "Pythia8/Pythia.h"

using std::endl;
using std::cout; 

class TClonesArray;
class TObjArray;
class TH1;

class AliTPythia8 : public TGenerator
{
public:
   AliTPythia8();
   AliTPythia8(const char *xmlDir);
   virtual ~AliTPythia8();
   static AliTPythia8     *Instance();
   static void             SetXmldocPath(char* path) {fgXmldocPath = path;}
   Pythia8::Pythia        *Pythia8() {return fPythia;}

   // Interface
   virtual void            GenerateEvent();
   virtual Int_t           ImportParticles(TClonesArray *particles, Option_t *option="");
   virtual TObjArray      *ImportParticles(Option_t *option="");

   // Others
   void                    ReadString(const char* string) const;
   void                    ReadConfigFile(const char* string) const;
   Bool_t                  Initialize(Int_t idAin, Int_t idBin, Double_t ecms);
   void                    PrintStatistics() const;
   void                    EventListing() const;
   Int_t                   GetN() const;

   void                    SetPythiaSeed(UInt_t seed);
   void                    ActivateEventSuperposition(TH1 *histNcoll);
   Pythia8::Pythia*        GetPythiaPartonLevel();
    
    
protected:
   AliTPythia8(const AliTPythia8&);
   AliTPythia8 operator=(const AliTPythia8&);
   static  AliTPythia8     *fgInstance;             //! singleton instance
   static  char            *fgXmldocPath;           //! path to xmldoc
   Pythia8::Pythia         *fPythia;                //! The pythia8 instance
   Pythia8::Pythia         *fPythiaPartonLevel;     //! The pythia8 instance for the parton level object (used for superposition of events)
   Int_t                    fNumberOfParticles;     //! Number of particles
   TH1                     *fSuperPositionNcoll;    //! Ncoll distribution for superposition mode
   Int_t                    fLastNMPI;              //! last number of MPI
   Int_t                    fLastNSuperposition;    //! last randomized superpositions
  
   void                    AddEvent(Pythia8::Pythia* target, Pythia8::Pythia* source);
   Bool_t                  CheckEvent(Pythia8::Event);
    
private:
   void                    AddParticlesToPdgDataBase() const;
   void                    UpdateParticleProperties() const;
   ClassDef(AliTPythia8, 3)   // Interface class of Pythia8
};

#endif

