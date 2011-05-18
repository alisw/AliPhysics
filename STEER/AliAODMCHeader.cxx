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

//-------------------------------------------------------------------------
//                      Implementation of   Class AliAODMCHeader
//   Header data
//   for the ESD   
//   Origin: Christian Klein-Boesing, CERN, Christian.Klein-Boesing@cern.ch 
//-------------------------------------------------------------------------

#include "AliAODMCHeader.h"


ClassImp(AliAODMCHeader)

// Without a trailing dot root does not support
// direct drawing of some variables if the name is not unique on top label
// bnrach e.g. fEventType is found here and in AliAODHeader....
TString AliAODMCHeader::fgkStdBranchName("mcHeader");

//______________________________________________________________________________

AliAODMCHeader::AliAODMCHeader() :
  AliVHeader()
  ,fGenerator("")
  ,fImpactPar(0)
  ,fPtHard(0)
  ,fXsection(0)
  ,fTrials(0)
  ,fEventType(0)
  ,fReactionPlaneAngle(0)  
{
  // default constructor
  fVertex[0] = fVertex[1] = fVertex[2] = 0;  
  SetName(fgkStdBranchName.Data());
}


AliAODMCHeader::~AliAODMCHeader() 
{
  // destructor
}


AliAODMCHeader::AliAODMCHeader(const AliAODMCHeader &header) :
  AliVHeader(header)
  ,fGenerator(header.fGenerator)
  ,fImpactPar(header.fImpactPar)
  ,fPtHard(header.fPtHard)
  ,fXsection(0)
  ,fTrials(0)
  ,fEventType(header.fEventType)
  ,fReactionPlaneAngle(header.fReactionPlaneAngle)  
{
  // copy constructor
  for(int i = 0;i<3;++i)fVertex[i] = header.fVertex[i];
  SetName(header.fName);
  SetTitle(header.fTitle);
}

AliAODMCHeader& AliAODMCHeader::operator=(const AliAODMCHeader &header)
{ 
  // assigment operator
  if(this!=&header) {
    AliVHeader::operator=(header);
    fGenerator = header.fGenerator;
    for(int i = 0;i<3;++i)fVertex[i] = header.fVertex[i];
    fImpactPar = header.fImpactPar;
    fPtHard = header.fPtHard;
    fXsection = header.fXsection;
    fTrials = header.fTrials;
    fEventType = header.fEventType;
    fReactionPlaneAngle = header.fReactionPlaneAngle;
  } 
  return *this;
}

void AliAODMCHeader::Copy(TObject &obj) const {
  
  // this overwrites the virtual TOBject::Copy()
  // to allow run time copying without casting
  // in AliESDEvent

  if(this==&obj)return;
  AliAODMCHeader *robj = dynamic_cast<AliAODMCHeader*>(&obj);
  if(!robj)return; // not an AliAODMCHeader
  *robj = *this;

}



//______________________________________________________________________________
void AliAODMCHeader::Reset()
{
  // reset all data members
  fGenerator = "";
  fImpactPar = 0;
  fEventType = 0;
  fPtHard = 0;
  fXsection = 0;
  fTrials = 0;
  fVertex[0] = fVertex[1] = fVertex[2] = 0;  
  fReactionPlaneAngle = 0;
}

//______________________________________________________________________________
void AliAODMCHeader::Print(const Option_t *) const
{
  // Print some data members
  Printf("MC EventHeader Generator: %s # EventType %d  Vtx = (%3.3f,%3.3f,%3.3f) ptHard = %3.3f GeV Impact parameter %3.3f  \n",
	 GetGeneratorName(),
	 GetEventType(),
	 GetVtxX(),GetVtxY(),GetVtxZ(),GetPtHard(),
	 GetImpactParameter());
}

