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
//
// Benchmarking class for V0 finder and PID. 
// Relies on MC information
// For more see source file
//

#ifndef ALIHFEV0PIDMC_H
#define ALIHFEV0PIDMC_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

class TList;

class AliMCEvent;

class AliHFEcollection;

class AliHFEV0pidMC : public TObject {

 public:
  AliHFEV0pidMC();
  ~AliHFEV0pidMC();
  
  void    Init();
  Bool_t  Process(TObjArray * const array, Int_t type);  

  void     SetMCEvent(AliMCEvent * const mc) { fMC = mc; };

  TList* GetListOfQAhistograms();

 private:
  AliHFEV0pidMC(const AliHFEV0pidMC &);
  AliHFEV0pidMC &operator=(const AliHFEV0pidMC &);
  Int_t PDGtoPIDdaughter(Int_t pdg) const;    // convert the PDG code to local PID
  Int_t PDGtoPIDmother(Int_t pdg) const;      // convert the PDG code to local PID
  Int_t PDGtoAliPID(Int_t pdg) const;         // convert PDG to AliPID

  AliMCEvent*         fMC;      // MC event
  AliHFEcollection*   fColl;    // Histogram collection

  UInt_t              fDestBits;    // status bits for the destructor

   ClassDef(AliHFEV0pidMC, 1)   // QA class for V0 PID
};
#endif
