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
//____________________________________________________________________
//
// Utility class for analysing ESD data. 
// This class does sharing and background correction 
//
#include <AliESDFMD.h>
#include "AliFMDAnaESD.h"
#include "AliFMDAnaRing.h"
#include <TMath.h>
#include <TBrowser.h>
// clude <AliLog.h>

//____________________________________________________________________
AliFMDAnaESD::AliFMDAnaESD()
{
  // Constructor 
  // Parameters:
  //     lower cut. 
  //
  AddLoad(kESD);
  for (size_t i = 0; i < 5; i++) fRing[i] = 0;
}
//____________________________________________________________________
Bool_t
AliFMDAnaESD::Init() 
{ 
  // Called at beginning of run 
  // Parameters:
  // Return
  //     @c false on error 
  //
  for (int i = 0; i < 5; i++) if (fRing[i]) fRing[i]->Init();
  return AliFMDInput::Init();
}
//____________________________________________________________________
Bool_t
AliFMDAnaESD::Begin(Int_t ev) 
{
  // Begining of event
  // Parameters:
  //     ev Event number
  // Return
  //     @c false on error 
  //
  fNEvents++;
  for (int i = 0; i < 5; i++) if (fRing[i]) fRing[i]->Begin();
  Bool_t ret = AliFMDInput::Begin(ev);
  return ret;
}
//____________________________________________________________________
Bool_t
AliFMDAnaESD::ProcessESDs()
{
  // Loop over all ESD data, and call ProcessESD for each entry.
  // Parameters:
  // Return
  //      @c false on error  
  //

  // Process event summary data
  if (!fESD) return kFALSE;
  for (UShort_t det = 1; det <= 3; det++) {
    Char_t rings[] = { 'I', (det == 1 ? '\0' : 'O'), '\0' };
    for (Char_t* rng = rings; *rng != '\0'; rng++) {
      UShort_t   nsec = (*rng == 'I' ?  20 :  40);
      UShort_t   nstr = (*rng == 'I' ? 512 : 256);
      Int_t      ridx = FindRing(det, *rng);
      if (ridx < 0 || !fRing[ridx]) continue;
      for (UShort_t sec = 0; sec < nsec; sec++) {
	for (UShort_t str = 0; str < nstr; str++) {
	  Float_t eta    = fESD->Eta(det,*rng,sec,str);
	  Float_t mult   = fESD->Multiplicity(det,*rng,sec,str);
	  Float_t multp1 = (str == nstr-1 ? 0 : 
			    fESD->Multiplicity(det,*rng,sec,str+1));
	  if (!fESD->IsAngleCorrected()) {
	    double c=TMath::Abs(TMath::Cos(2.*TMath::ATan(TMath::Exp(-eta))));
	    mult    *= c;
	    multp1  *= c;
	  }

	  // Check signal isn't pedestal 
	  if (mult < 0.001) continue;

	  // Phi 
	  Double_t dPhi = 2*TMath::Pi() / fRing[ridx]->NSeq();
	  Double_t phi; 
	  if (det == 3) phi =  TMath::Pi() - (sec + .5) * dPhi;
	  else   	phi =  (sec + .5) * dPhi;
	  if (phi < 0)  phi += 2 * TMath::Pi();

	  if (fRing[ridx]->ProcessESD(phi, eta, mult, multp1)) str++;
	  // "mult" possibly updated 
	  Fill(phi, eta, mult);
	}
      }
    }
  }
  return kTRUE;
}
//____________________________________________________________________
Bool_t
AliFMDAnaESD::End() 
{
  for (int i = 0; i < 5; i++) if (fRing[i]) fRing[i]->End();
  return AliFMDInput::End();
}
//____________________________________________________________________
Bool_t
AliFMDAnaESD::Finish() 
{
  // Called at the end of run 
  // Parameters:
  // Return
  //     @c false in case of errors 
  //
  for (int i = 0; i < 5; i++) if (fRing[i]) fRing[i]->Finish();
  return AliFMDInput::Finish();
}
//__________________________________________________________________
void 
AliFMDAnaESD::AddRing(AliFMDAnaRing* ring)
{
  // Add a ring 
  // Parameters:
  //     ring Ring object 
  //
  fRing[FindRing(ring->Detector(),ring->Ring())] = ring;
}
//____________________________________________________________________
Int_t
AliFMDAnaESD::FindRing(UShort_t det, Char_t ring) const
{ 
  // Find a ring index
  // Parameters:
  //     det  Detector number 
  //     ring Ring id
  // Return
  //     Index of ring object 
  //
  Int_t idx = -1;
  switch (det) { 
  case 1: idx = 0; break;
  case 2: idx = 1 + (ring == 'O' || ring == 'o' ? 1 : 0); break;
  case 3: idx = 3 + (ring == 'O' || ring == 'o' ? 1 : 0); break;
  }
  return idx;
}

//____________________________________________________________________
void
AliFMDAnaESD::Browse(TBrowser* b)
{
  for (size_t i = 0; i < 5; i++)
    if (fRing[i]) b->Add(fRing[i], fRing[i]->Name());
}

//____________________________________________________________________
//
// EOF
//
