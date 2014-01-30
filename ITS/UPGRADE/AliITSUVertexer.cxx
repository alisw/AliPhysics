/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
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
#include "AliITSUVertexer.h"
#include "AliESDVertex.h"
#include "AliLog.h"
#include <TString.h>
#include <TTree.h>

/////////////////////////////////////////////////////////////////
// THIS IS A DUMMY CLASS
////////////////////////////////////////////////////////////////

ClassImp(AliITSUVertexer)



//______________________________________________________________________
AliITSUVertexer::AliITSUVertexer():AliVertexer()
{
}


//______________________________________________________________________
AliITSUVertexer::~AliITSUVertexer() {
  // Destructor
}

//______________________________________________________________________
AliESDVertex* AliITSUVertexer::FindVertexForCurrentEvent(TTree *)
{
  Double_t startPos[3]={GetNominalPos()[0],GetNominalPos()[1],GetNominalPos()[2]};
  Double_t startCov[6]={GetNominalCov()[0],GetNominalCov()[1],GetNominalCov()[2],
			GetNominalCov()[3],GetNominalCov()[4],GetNominalCov()[5]};
  AliESDVertex* vtx = new AliESDVertex(startPos,startCov,99999.,-2);
  AliInfo("Creating dummy vertex from the mean vertex");
  vtx->Print();
  return vtx;
}  

