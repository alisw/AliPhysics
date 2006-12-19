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
 
//------------------------------------------------------------------------
// Jet reader base class
// manages the reading of input for jet algorithms
// Authors: jgcn@mda.cinvestav.mx
//          magali.estienne@IReS.in2p3.fr
//------------------------------------------------------------------------

// root
#include <TClonesArray.h>
//AliRoot
#include "AliJetReader.h"
#include "AliJetReaderHeader.h"
#include "AliJetUnitArray.h"
#include "AliJetHadronCorrectionv1.h"

ClassImp(AliJetReader)

////////////////////////////////////////////////////////////////////////

AliJetReader::AliJetReader():
  fMomentumArray(new TClonesArray("TLorentzVector",2000)),
  fArrayMC(0),
  fFillUnitArray(new TTask("fillUnitArray","Fill unit array jet finder")),
  fReaderHeader(0),
  fSignalFlag(0),
  fCutFlag(0),
  fUnitArray(new AliJetUnitArray[60000]),     
  fUnitArrayNoCuts(new AliJetUnitArray[60000]),
  fArrayInitialised(0)
{
  // Default constructor
  fSignalFlag = TArrayI();
  fCutFlag    = TArrayI();
}

////////////////////////////////////////////////////////////////////////

AliJetReader::~AliJetReader()
{
  // Destructor
  if (fMomentumArray) {
      fMomentumArray->Delete();
      delete fMomentumArray;
  }
  
  if (fUnitArray) {
      fUnitArray->Delete();
      delete fUnitArray;
  }
  
  if (fUnitArrayNoCuts) {
    fUnitArrayNoCuts->Delete();
    delete fUnitArrayNoCuts;
  }

  if (fFillUnitArray) {
    fFillUnitArray->Delete();
    delete fFillUnitArray;
  }
  delete fArrayMC;
  
}


////////////////////////////////////////////////////////////////////////

void AliJetReader::ClearArray()

{
  if (fMomentumArray)  fMomentumArray->Clear();
  if (fFillUnitArray)  fFillUnitArray->Clear();
}
