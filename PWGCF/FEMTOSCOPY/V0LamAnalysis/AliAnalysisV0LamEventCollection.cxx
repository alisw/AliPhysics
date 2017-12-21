///
/// \file V0LamAnalysis/AliAnalysisV0LamEventCollection.cxx
///


#include "AliAnalysisV0LamEventCollection.h"

using namespace std;

//_____________________________________________________________________________
AliAnalysisV0LamEventCollection::~AliAnalysisV0LamEventCollection()
{
  if(fEvt) {
    delete [] fEvt;
    fEvt = NULL;
  }
}
//_____________________________________________________________________________
AliAnalysisV0LamEventCollection::AliAnalysisV0LamEventCollection():
  fEvt(NULL),
  fFifoSize(0)
{
}

//_____________________________________________________________________________
AliAnalysisV0LamEventCollection::AliAnalysisV0LamEventCollection(short eventBuffSize, int maxV0Mult)
{
  SetBuffSize(eventBuffSize);
  fEvt = new AliAnalysisV0LamEvent[fFifoSize];  //allocate pointer array of AliAnalysisV0LamEvents
  for(int ii = 0; ii < fFifoSize; ii++){ //Initialize particle table pointers to NULL
    (fEvt + ii)->fMaxV0Mult = maxV0Mult;
    (fEvt + ii)->fReconstructedV0 = NULL;
    (fEvt + ii)->fNumberCandidateV0 = 0;
    (fEvt + ii)->fReconstructedV0 = new AliReconstructedV0[maxV0Mult];
    for(int j=0; j < 3; j++){
      (fEvt + ii)->fPrimaryVertex(j) = 0.;
    }
  }
}


//_____________________________________________________________________________
AliAnalysisV0LamEventCollection::AliAnalysisV0LamEventCollection(const AliAnalysisV0LamEventCollection &eventCollection)
{
  SetBuffSize(eventCollection.GetFifoSize());
  fEvt = new AliAnalysisV0LamEvent[fFifoSize];
  for(Int_t i = 0; i < fFifoSize; i++) {
    fEvt[i] = AliAnalysisV0LamEvent(eventCollection.fEvt[i]);
  }
}



//_____________________________________________________________________________
AliAnalysisV0LamEventCollection& AliAnalysisV0LamEventCollection::operator=(const AliAnalysisV0LamEventCollection& eventCollection)
{
  if (&eventCollection == this) return *this;
  SetBuffSize(eventCollection.GetFifoSize());
  if(fEvt) {
    delete [] fEvt;
  }
  fEvt = new AliAnalysisV0LamEvent[fFifoSize];
  for(Int_t i = 0; i < fFifoSize; i++) {
    fEvt[i] = AliAnalysisV0LamEvent(eventCollection.fEvt[i]);
  }
  return *this;
}


//_____________________________________________________________________________
void AliAnalysisV0LamEventCollection::FifoShift()
{
  //Shift elements in FIFO by one and clear last element in FIFO
  for(unsigned short i=fFifoSize-1 ; i > 0; i--){
    for(int j=0; j<(fEvt + i-1)->fNumberCandidateV0; j++){
      (fEvt + i)->fReconstructedV0[j] = (fEvt + i-1)->fReconstructedV0[j];
    }
    (fEvt + i)->fNumberCandidateV0 = (fEvt + i-1)->fNumberCandidateV0;
    for(int j=0; j<3; j++){
      (fEvt + i)->fPrimaryVertex(j) = (fEvt + i-1)->fPrimaryVertex(j);
    }
  }
  (fEvt)->fNumberCandidateV0=0;
  for(int j=0; j<3; j++){
    (fEvt)->fPrimaryVertex(j) = 0.;
  }
}

//_____________________________________________________________________________
AliReconstructedV0::AliReconstructedV0()
{
}

//_____________________________________________________________________________
AliReconstructedV0::~AliReconstructedV0()
{
}

//_____________________________________________________________________________
AliAnalysisV0LamEvent::AliAnalysisV0LamEvent():
  fMaxV0Mult(1),
fNumberCandidateV0(0),
fPrimaryVertex(0,0,0),
fReconstructedV0(NULL)
{

}

//_____________________________________________________________________________
AliAnalysisV0LamEvent::AliAnalysisV0LamEvent(const AliAnalysisV0LamEvent &event):
  fMaxV0Mult(event.fMaxV0Mult),
  fNumberCandidateV0(event.fNumberCandidateV0),
  fPrimaryVertex(event.fPrimaryVertex)
{
  fReconstructedV0 = new AliReconstructedV0[fMaxV0Mult];
  for(Int_t i = 0; i < fNumberCandidateV0; i++) {
    fReconstructedV0[i] = AliReconstructedV0(event.fReconstructedV0[i]);
  }
}

//_____________________________________________________________________________
AliAnalysisV0LamEvent& AliAnalysisV0LamEvent::operator=(const AliAnalysisV0LamEvent &event)
{
  if (&event == this) return *this;
  fMaxV0Mult = event.fMaxV0Mult;
  fNumberCandidateV0 = event.fNumberCandidateV0;
  fPrimaryVertex = event.fPrimaryVertex;
  if(fReconstructedV0) {
    delete [] fReconstructedV0;
  }
  fReconstructedV0 = new AliReconstructedV0[fMaxV0Mult];
  for(Int_t i = 0; i < fNumberCandidateV0; i++) {
    fReconstructedV0[i] = AliReconstructedV0(event.fReconstructedV0[i]);
  }
  return *this;
}

//_____________________________________________________________________________
AliAnalysisV0LamEvent::~AliAnalysisV0LamEvent()
{
  if(fReconstructedV0) {
    delete [] fReconstructedV0;
    fReconstructedV0 = NULL;
  }
}
