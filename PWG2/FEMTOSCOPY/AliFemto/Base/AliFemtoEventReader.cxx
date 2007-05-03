////////////////////////////////////////////////////////////////////////////////
/// AliFemtoEventReader - the pure virtual base class for the event reader   ///
/// All event readers must inherit from this one                             ///
////////////////////////////////////////////////////////////////////////////////
#include "Infrastructure/AliFemtoEvent.h"
#include "Base/AliFemtoEventCut.h"
#include "Base/AliFemtoTrackCut.h"
#include "Base/AliFemtoV0Cut.h"
#include "Base/AliFemtoXiCut.h"
#include "Base/AliFemtoKinkCut.h"
#include "Base/AliFemtoEventReader.h"

#ifdef __ROOT__
ClassImp(AliFemtoEventReader)
#endif

AliFemtoEventReader::AliFemtoEventReader(const AliFemtoEventReader& aReader):
  fEventCut(0),  
  fTrackCut(0),    
  fV0Cut(0),       
  fXiCut(0),       
  fKinkCut(0),    
  fReaderStatus(0),  
  fDebug(0)
{
  fEventCut = aReader.fEventCut;
  fTrackCut = aReader.fTrackCut;
  fV0Cut    = aReader.fV0Cut;
  fXiCut    = aReader.fXiCut;
  fKinkCut  = aReader.fKinkCut;
  fReaderStatus = aReader.fReaderStatus;
  fDebug = aReader.fDebug;
}

AliFemtoEventReader& AliFemtoEventReader::operator=(const AliFemtoEventReader& aReader)
{
  if (this == &aReader) 
    return *this;

  fEventCut = aReader.fEventCut;
  fTrackCut = aReader.fTrackCut;
  fV0Cut    = aReader.fV0Cut;
  fXiCut    = aReader.fXiCut;
  fKinkCut  = aReader.fKinkCut;
  fReaderStatus = aReader.fReaderStatus;
  fDebug = aReader.fDebug;

  return *this;
}


AliFemtoString AliFemtoEventReader::Report(){
  // Create a simple report from the workings of the reader
  AliFemtoString temp = "\n This is the base class AliFemtoEventReader reporting";
  temp += "\n---> EventCuts in Reader: ";
  if (fEventCut) {
    temp += fEventCut->Report();
  }
  else {
    temp += "NONE";
  }
  temp += "\n---> TrackCuts in Reader: ";
  if (fTrackCut) {
    temp += fTrackCut->Report();
  }
  else {
    temp += "NONE";
  }
  temp += "\n---> V0Cuts in Reader: ";
  if (fV0Cut) {
    temp += fV0Cut->Report();
  }
  else {
    temp += "NONE";
  }
  temp += "\n---> XiCuts in Reader: ";
  if (fXiCut) {
    temp += fXiCut->Report();
  }
  else {
    temp += "NONE";
  }
  temp += "\n---> KinkCuts in Reader: ";
  if (fKinkCut) {
    temp += fKinkCut->Report();
  }
  else {
    temp += "NONE";
  }
  temp += "\n";
  return temp;
}
//______________________________________
void AliFemtoEventReader::SetEventCut(AliFemtoEventCut* ecut){fEventCut=ecut;}
//______________________________________
void AliFemtoEventReader::SetTrackCut(AliFemtoTrackCut* pcut){cout << pcut << endl; fTrackCut=pcut;}
//______________________________________
void AliFemtoEventReader::SetV0Cut(AliFemtoV0Cut* pcut){fV0Cut=pcut;}
//______________________________________
void AliFemtoEventReader::SetXiCut(AliFemtoXiCut* pcut){fXiCut=pcut;}
//______________________________________
void AliFemtoEventReader::SetKinkCut(AliFemtoKinkCut* pcut){fKinkCut=pcut;}
//______________________________________
AliFemtoEventCut* AliFemtoEventReader::EventCut(){return fEventCut;}
//______________________________________
AliFemtoTrackCut* AliFemtoEventReader::TrackCut(){return fTrackCut;}
//______________________________________
AliFemtoV0Cut*    AliFemtoEventReader::V0Cut(){return fV0Cut;} 
//______________________________________
AliFemtoXiCut*    AliFemtoEventReader::XiCut(){return fXiCut;} 
//______________________________________
AliFemtoKinkCut*    AliFemtoEventReader::KinkCut(){return fKinkCut;}


