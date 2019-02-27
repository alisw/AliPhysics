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

#include "TList.h"
#include "AliAODMCHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenCocktailEventHeader.h"




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
  ,fBGEventReused(0)
  ,fReactionPlaneAngle(0)  
  ,fHeaders(0)
{
  /// default constructor

  fVertex[0] = fVertex[1] = fVertex[2] = 0;  
  SetName(fgkStdBranchName.Data());
}


AliAODMCHeader::~AliAODMCHeader() 
{

  Reset();
  delete fHeaders;
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
  ,fBGEventReused(header.fBGEventReused)
  ,fReactionPlaneAngle(header.fReactionPlaneAngle)  
  ,fHeaders(0)
{
  /// copy constructor

  for(int i = 0;i<3;++i)fVertex[i] = header.fVertex[i];
  SetName(header.fName);
  SetTitle(header.fTitle);
}

AliAODMCHeader& AliAODMCHeader::operator=(const AliAODMCHeader &header)
{ 
  /// assigment operator

  if(this!=&header) {
    Reset();
    AliVHeader::operator=(header);
    fGenerator = header.fGenerator;
    for(int i = 0;i<3;++i)fVertex[i] = header.fVertex[i];
    fImpactPar = header.fImpactPar;
    fPtHard = header.fPtHard;
    fXsection = header.fXsection;
    fTrials = header.fTrials;
    fEventType = header.fEventType;
    fBGEventReused = header.fBGEventReused;
    fReactionPlaneAngle = header.fReactionPlaneAngle;

    if(header.fHeaders){
      for(int i = 0;i < header.fHeaders->GetEntries();i++){
	AddCocktailHeader(dynamic_cast<AliGenEventHeader*>(header.fHeaders->At(i)));
      }
    }
  } 
  return *this;
}

void AliAODMCHeader::AddCocktailHeader(const AliGenEventHeader* header)
{
/// Add a header to the list

  if(!header)return;
  if (!fHeaders){ 
    fHeaders = new TList();
    fHeaders->SetOwner(kTRUE);
  }
  fHeaders->Add(header->Clone());
}

void AliAODMCHeader::Copy(TObject &obj) const {
  
  /// this overwrites the virtual TOBject::Copy()
  /// to allow run time copying without casting
  /// in AliESDEvent

  if(this==&obj)return;
  AliAODMCHeader *robj = dynamic_cast<AliAODMCHeader*>(&obj);
  if(!robj)return; // not an AliAODMCHeader
  *robj = *this;

}



//______________________________________________________________________________
void AliAODMCHeader::Reset()
{
  /// reset all data members

  fGenerator = "";
  fImpactPar = 0;
  fEventType = 0;
  fBGEventReused = 0;
  fPtHard = 0;
  fXsection = 0;
  fTrials = 0;
  fVertex[0] = fVertex[1] = fVertex[2] = 0;  
  fReactionPlaneAngle = 0;
  if(fHeaders)fHeaders->Delete();
}

//______________________________________________________________________________
void AliAODMCHeader::Print(const Option_t *) const
{
  /// Print some data members

  Printf("MC EventHeader Generators: %s # EventType %d  Vtx = (%3.3f,%3.3f,%3.3f) ptHard = %3.3f GeV Impact parameter %3.3f  \n",
	 GetGeneratorName(),
	 GetEventType(),
	 GetVtxX(),GetVtxY(),GetVtxZ(),GetPtHard(),
	 GetImpactParameter());
  if(fHeaders){
    fHeaders->Print();
    for(int i = 0;i<fHeaders->GetEntries();++i){
      TObject *obj = fHeaders->At(i);
      if(obj){
	Printf(">> %d: %s %s",i,obj->GetName(),obj->GetTitle());
      }
    }
  }
}

AliGenEventHeader* AliAODMCHeader::GetCocktailHeader(Int_t i){
  if(i<0)return 0;
  return (AliGenEventHeader*)(fHeaders->At(i));
}

void  AliAODMCHeader::AddCocktailHeaders(AliGenEventHeader* header){
  AliGenCocktailEventHeader *cHeader = dynamic_cast<AliGenCocktailEventHeader*>(header);
  if(cHeader){
      TList *genHeaders = cHeader->GetHeaders();
      AliGenEventHeader* gH = 0;
      for (int i=0; i<genHeaders->GetEntries(); i++) {
	gH = (AliGenEventHeader*)genHeaders->At(i);
	if(gH){
	  AddGeneratorName(gH->GetName());
	  AddCocktailHeader(dynamic_cast<AliGenEventHeader*>(genHeaders->At(i)));	
	}
      }
  }
  else{
    // no cocktail header just addd the global header
    AddCocktailHeader(header);	
  }
}

void   AliAODMCHeader::AddGeneratorName(const char* c){
  if(fGenerator.Length()==0)fGenerator += c;
  else {
    fGenerator += " ";
    fGenerator += c;
  }
}
