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

/*
$Log$
Revision 1.2  2002/10/14 14:57:40  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.1.2.8  2002/10/08 16:33:17  iglez2
LSOUIT is set to true before the second call to flukam.

Revision 1.1.2.7  2002/10/08 09:30:37  iglez2
Solved stupid missing ;

Revision 1.1.2.6  2002/10/07 13:40:22  iglez2
First implementations of the PDG <--> Fluka Id conversion routines

Revision 1.1.2.5  2002/09/26 16:26:03  iglez2
Added verbosity
Call to gAlice->Generator()->Generate()

Revision 1.1.2.4  2002/09/26 13:22:23  iglez2
Naive implementation of ProcessRun and ProcessEvent
Opening/Closing of input file (fInputFileName) with FORTRAN unit 5 before/after the first call to flukam inside Init()

Revision 1.1.2.3  2002/09/20 15:35:51  iglez2
Modification of LFDRTR. Value is passed to FLUKA !!!

Revision 1.1.2.2  2002/09/18 14:34:44  iglez2
Revised version with all pure virtual methods implemented

Revision 1.1.2.1  2002/07/24 08:49:41  alibrary
Adding TFluka to VirtualMC

Revision 1.1  2002/07/05 13:10:07  morsch
First commit of Fluka interface.

*/

#include <Riostream.h>

#include "TFluka.h"
#include "TCallf77.h"      //For the fortran calls
#include "Fdblprc.h"       //(DBLPRC) fluka common
#include "Fiounit.h"       //(IOUNIT) fluka common
#include "Fepisor.h"       //(EPISOR) fluka common
#include "TVirtualMC.h"

// Fluka methods that may be needed.
#ifndef WIN32 
# define flukam  flukam_
# define fluka_openinp fluka_openinp_
# define fluka_closeinp fluka_closeinp_
#else 
# define flukam  FLUKAM
# define fluka_openinp FLUKA_OPENINP
# define fluka_closeinp FLUKA_CLOSEINP
#endif

extern "C" 
{
  //
  // Prototypes for FLUKA functions
  //
  void type_of_call flukam(const int&);
  void type_of_call fluka_openinp(const int&, DEFCHARA);
  void type_of_call fluka_closeinp(const int&);
}

//
// Class implementation for ROOT
//
ClassImp(TFluka)

//
// TFluka methods.
//____________________________________________________________________________ 
TFluka::TFluka()
  :TVirtualMC(),
   fVerbosityLevel(0),
   fInputFileName("")
{ 
  //
  // Default constructor
  //
} 
 
//____________________________________________________________________________ 
TFluka::TFluka(const char *title, Int_t verbosity)
  :TVirtualMC("TFluka",title),
   fVerbosityLevel(verbosity),
   fInputFileName("")
{
  if (fVerbosityLevel >=3)
    cout << "==> TFluka::TFluka(" << title << ") constructor called." << endl;

  if (fVerbosityLevel >=3)
    cout << "<== TFluka::TFluka(" << title << ") constructor called." << endl;
}

//____________________________________________________________________________ 
void TFluka::Init() {
  if (fVerbosityLevel >=3)
    cout << "==> TFluka::Init() called." << endl;

  if (fVerbosityLevel >=2)
    cout << "\t* Changing lfdrtr = (" << (GLOBAL.lfdrtr?'T':'F')
	 << ") in fluka..." << endl;
  GLOBAL.lfdrtr = true;

  if (fVerbosityLevel >=2)
    cout << "\t* Opening file " << fInputFileName << endl;
  const char* fname = fInputFileName;
  fluka_openinp(lunin, PASSCHARA(fname));

  if (fVerbosityLevel >=2)
    cout << "\t* Calling flukam..." << endl;
  flukam(0);

  if (fVerbosityLevel >=2)
    cout << "\t* Closing file " << fInputFileName << endl;
  fluka_closeinp(lunin);

  if (fVerbosityLevel >=3)
    cout << "<== TFluka::Init() called." << endl;
}

//____________________________________________________________________________ 
void TFluka::ProcessEvent() {
  if (fVerbosityLevel >=3)
    cout << "==> TFluka::ProcessEvent() called." << endl;

  if (fVerbosityLevel >=3)
    cout << "<== TFluka::ProcessEvent() called." << endl;
}

//____________________________________________________________________________ 
void TFluka::ProcessRun(Int_t nevent) {
  if (fVerbosityLevel >=3)
    cout << "==> TFluka::ProcessRun(" << nevent << ") called." 
	 << endl;

  if (fVerbosityLevel >=2) {
    cout << "\t* GLOBAL.fdrtr = " << (GLOBAL.lfdrtr?'T':'F') << endl;
    cout << "\t* Calling flukam again..." << endl;
  }
  fApplication->GeneratePrimaries();
  EPISOR.lsouit = true;
  flukam(0);

  if (fVerbosityLevel >=3)
    cout << "<== TFluka::ProcessRun(" << nevent << ") called." 
	 << endl;
}

//_____________________________________________________________________________
Int_t TFluka::IdFromPDG(Int_t pdg) const 
{
  //
  // Return Geant3 code from PDG and pseudo ENDF code
  //
  for(Int_t i=0;i<fNPDGCodes;++i)
    if(pdg==fPDGCode[i])
      return i;
  return -99;
}

//_____________________________________________________________________________
Int_t TFluka::PDGFromId(Int_t id) const 
{
  //
  // Return PDG code and pseudo ENDF code from Geant3 code
  //
  if(id>0 && id<fNPDGCodes) 
    return fPDGCode[id];
  else 
    return -1;
}

//_____________________________________________________________________________
void TFluka::DefineParticles() 
{
  // Load standard numbers for GEANT particles and PDG conversion
  fPDGCode[fNPDGCodes++]=  -99; //  0 = Psudoparticle (Ray)
  fPDGCode[fNPDGCodes++]= 2212; //  1 = Proton
  fPDGCode[fNPDGCodes++]=-2212; //  2 = Anti Proton
  fPDGCode[fNPDGCodes++]=   11; //  3 = Electron
  fPDGCode[fNPDGCodes++]=  -11; //  4 = Positron
  fPDGCode[fNPDGCodes++]=   12; //  5 = Electron Neutrino
  fPDGCode[fNPDGCodes++]=  -12; //  6 = Electron Antineutrino
  fPDGCode[fNPDGCodes++]=   22; //  7 = Photon
  fPDGCode[fNPDGCodes++]= 2112; //  8 = Neutron
  fPDGCode[fNPDGCodes++]=-2112; //  9 = Anti Neutron
  fPDGCode[fNPDGCodes++]=  -13; // 10 = Mu+
  fPDGCode[fNPDGCodes++]=   13; // 11 = Mu-
  fPDGCode[fNPDGCodes++]=  130; // 12 = Kaon 0 long
  fPDGCode[fNPDGCodes++]=  211; // 13 = Pi+
  fPDGCode[fNPDGCodes++]= -211; // 14 = Pi-
  fPDGCode[fNPDGCodes++]=  321; // 15 = Kaon+
  fPDGCode[fNPDGCodes++]= -321; // 16 = Kaon-
  fPDGCode[fNPDGCodes++]= 3122; // 17 = Lambda
  fPDGCode[fNPDGCodes++]=-3122; // 18 = Anti Lambda
  fPDGCode[fNPDGCodes++]=  310; // 19 = Kaon 0 short
  fPDGCode[fNPDGCodes++]= 3112; // 20 = Sigma -
  fPDGCode[fNPDGCodes++]= 3222; // 21 = Sigma +
  fPDGCode[fNPDGCodes++]= 3212; // 22 = Sigma 0
  fPDGCode[fNPDGCodes++]=  111; // 23 = Pi0
  fPDGCode[fNPDGCodes++]=  311; // 24 = Kaon 0
  fPDGCode[fNPDGCodes++]= -311; // 25 = Antikaon 0
  fPDGCode[fNPDGCodes++]=  -99; // 26 = --Reserved
  fPDGCode[fNPDGCodes++]=   14; // 27 = Muon neutrino
  fPDGCode[fNPDGCodes++]=  -14; // 28 = Muon antineutrino
  fPDGCode[fNPDGCodes++]=  -99; // 29 = --Reserved
  fPDGCode[fNPDGCodes++]=  -99; // 30 = --Reserved
  fPDGCode[fNPDGCodes++]=-3222; // 31 = Antisigma -
  fPDGCode[fNPDGCodes++]=-3212; // 32 = Antisigma 0
  fPDGCode[fNPDGCodes++]=-3112; // 33 = Antisigma +
  fPDGCode[fNPDGCodes++]= 3322; // 34 = Xi 0
  fPDGCode[fNPDGCodes++]=-3322; // 35 = AntiXi 0
  fPDGCode[fNPDGCodes++]= 3312; // 36 = Xi -
  fPDGCode[fNPDGCodes++]=-3312; // 37 = Xi +
  fPDGCode[fNPDGCodes++]= 3334; // 38 = Omega -
  fPDGCode[fNPDGCodes++]=-3334; // 39 = Antiomega
  fPDGCode[fNPDGCodes++]=  -99; // 40 = --Reserved
  fPDGCode[fNPDGCodes++]=  -15; // 41 = Tau+
  fPDGCode[fNPDGCodes++]=   15; // 42 = Tau-
  fPDGCode[fNPDGCodes++]=   16; // 43 = Tau neutrino
  fPDGCode[fNPDGCodes++]=  -16; // 44 = Tau antineutrino
  fPDGCode[fNPDGCodes++]=  411; // 45 = D+
  fPDGCode[fNPDGCodes++]= -411; // 46 = D-
  fPDGCode[fNPDGCodes++]=  421; // 47 = D0
  fPDGCode[fNPDGCodes++]= -421; // 48 = AntiD 0
  fPDGCode[fNPDGCodes++]=  431; // 49 = D_s +
  fPDGCode[fNPDGCodes++]= -431; // 50 = D_s -
  fPDGCode[fNPDGCodes++]= 4122; // 51 = Lambda_c +
  fPDGCode[fNPDGCodes++]= 4232; // 52 = Xi_c +
  fPDGCode[fNPDGCodes++]= 4112; // 53 = Xi_c -
  fPDGCode[fNPDGCodes++]= 4322; // 54 = Xi'_c +
  fPDGCode[fNPDGCodes++]= 4312; // 55 = Xi'_c 0
  fPDGCode[fNPDGCodes++]= 4332; // 56 = Omega_c 0
  fPDGCode[fNPDGCodes++]=-4122; // 57 = Antilambda_c -
  fPDGCode[fNPDGCodes++]=-4232; // 58 = Antixsi_c -
  fPDGCode[fNPDGCodes++]=-4112; // 59 = Antixsi_c 0
  fPDGCode[fNPDGCodes++]=-4322; // 60 = AntiXi'_c -
  fPDGCode[fNPDGCodes++]=-4312; // 61 = AntiXi'_c 0
  fPDGCode[fNPDGCodes++]=-4332; // 62 = AntiOmega_c 0
  fPDGCode[fNPDGCodes++]=  -99; // 63 = --Reserved
  fPDGCode[fNPDGCodes++]=  -99; // 64 = --Reserved
}
