#include "AliAnalysisV0LamEventCollection.h"
//_____________________________________________________________________________
AliAnalysisV0LamEventCollection::~AliAnalysisV0LamEventCollection(){
  for(int i = 0; i < fifo; i++){
    if((fEvt + i)->fReconstructedV0 != NULL){
      delete [] (fEvt + i)->fReconstructedV0;
    }	
  }
  delete [] fEvt;
}
//_____________________________________________________________________________
AliAnalysisV0LamEventCollection::AliAnalysisV0LamEventCollection(){}
//_____________________________________________________________________________
AliAnalysisV0LamEventCollection::AliAnalysisV0LamEventCollection(short eventBuffSize, int maxV0Mult){
  SetBuffSize(eventBuffSize);
  fEvt = new AliAnalysisV0LamEvent[fifo];  //allocate pointer array of AliAnalysisV0LamEvents
  for(int ii = 0; ii < fifo; ii++){ //Initialize particle table pointers to NULL
    (fEvt + ii)->fReconstructedV0 = NULL;
    (fEvt + ii)->fNumberCandidateV0 = 0;
    (fEvt + ii)->fReconstructedV0 = new AliReconstructedV0[maxV0Mult];
    for(int j=0; j < 3; j++){
      (fEvt + ii)->fPrimaryVertex[j] = 0.;
    }
  }
}
//_____________________________________________________________________________
void AliAnalysisV0LamEventCollection::FifoShift(){ //Shift elements in FIFO by one and clear last element in FIFO 
  for(unsigned short i=fifo-1 ; i > 0; i--){
    for(int j=0; j<(fEvt + i-1)->fNumberCandidateV0; j++){
      (fEvt + i)->fReconstructedV0[j] = (fEvt + i-1)->fReconstructedV0[j];
    }
    (fEvt + i)->fNumberCandidateV0 = (fEvt + i-1)->fNumberCandidateV0;
    for(int j=0; j<3; j++){
      (fEvt + i)->fPrimaryVertex[j] = (fEvt + i-1)->fPrimaryVertex[j];
    }
  }
  (fEvt)->fNumberCandidateV0=0;
  for(int j=0; j<3; j++){
    (fEvt)->fPrimaryVertex[j] = 0.;
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
