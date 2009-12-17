#include <TNamed.h>
#include <TMath.h>

#include "AliVHeader.h"
#include "AliVVertex.h"
#include "AliMuonsHFHeader.h"

ClassImp(AliMuonsHFHeader)

//_____________________________________________________________________________
AliMuonsHFHeader::AliMuonsHFHeader() :
TNamed(),
fTriggerMask(0),
fNContributors(0),
fMultMuon(0),
fMultDimuon(0),
fCentrality(0.),
fUnrecoVertex(kFALSE)
{
  //
  // default constructor
  //
  for (Int_t i=0; i<3; i++) fPosition[i]=0.;
}

//_____________________________________________________________________________
AliMuonsHFHeader::~AliMuonsHFHeader()
{
  //
  // default destructor
  //
}

//_____________________________________________________________________________
void AliMuonsHFHeader::SetHeader(AliVHeader *header)
{
  fTriggerMask = header->GetTriggerMask();
  return;
}

//_____________________________________________________________________________
void AliMuonsHFHeader::SetVertex(AliVVertex *vertex)
{
  vertex->GetXYZ(fPosition);
  fNContributors = vertex->GetNContributors();
  fUnrecoVertex = (TMath::Abs(fPosition[0])<1e-6 && TMath::Abs(fPosition[1])<1e-6 &&
                   TMath::Abs(fPosition[2])<1e-6);
  return;
}
