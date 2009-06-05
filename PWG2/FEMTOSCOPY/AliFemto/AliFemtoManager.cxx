///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoManager: main class managing femtoscopic analysis             //
// The Manager is the top-level object that coordinates activities       //
// and performs event, particle, and pair loops, and checks the          //
// various Cuts of the Analyses in its AnalysisCollection                //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "AliFemtoManager.h"
//#include "AliFemtoParticleCollection.h"
//#include "AliFemtoTrackCut.h"
//#include "AliFemtoV0Cut.h"
#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoManager)
#endif



//____________________________
AliFemtoManager::AliFemtoManager():
  fAnalysisCollection(0),
  fEventReader(0),
  fEventWriterCollection(0)
{
  // default constructor
  fAnalysisCollection = new AliFemtoAnalysisCollection;
  fEventWriterCollection = new AliFemtoEventWriterCollection;
  fEventReader = 0;
}
//____________________________
AliFemtoManager::AliFemtoManager(const AliFemtoManager& aManager) :
  fAnalysisCollection(0),
  fEventReader(0),
  fEventWriterCollection(0)
{
  // copy constructor
  fEventReader = aManager.fEventReader;
  AliFemtoSimpleAnalysisIterator tAnalysisIter;
  fAnalysisCollection = new AliFemtoAnalysisCollection;
  for (tAnalysisIter=aManager.fAnalysisCollection->begin();tAnalysisIter!=aManager.fAnalysisCollection->end();tAnalysisIter++){
    fAnalysisCollection->push_back(*tAnalysisIter);
  }
  AliFemtoEventWriterIterator tEventWriterIter;
  fEventWriterCollection = new AliFemtoEventWriterCollection;
  for (tEventWriterIter=aManager.fEventWriterCollection->begin();tEventWriterIter!=aManager.fEventWriterCollection->end();tEventWriterIter++){
    fEventWriterCollection->push_back(*tEventWriterIter);
  }
}

//____________________________
AliFemtoManager::~AliFemtoManager(){
  // destructor
  delete fEventReader;
  // now delete each Analysis in the Collection, and then the Collection itself
  AliFemtoSimpleAnalysisIterator tAnalysisIter;
  for (tAnalysisIter=fAnalysisCollection->begin();tAnalysisIter!=fAnalysisCollection->end();tAnalysisIter++){
    delete *tAnalysisIter;
    *tAnalysisIter = 0;
  }
  delete fAnalysisCollection;
  // now delete each EventWriter in the Collection, and then the Collection itself
  AliFemtoEventWriterIterator tEventWriterIter;
  for (tEventWriterIter=fEventWriterCollection->begin();tEventWriterIter!=fEventWriterCollection->end();tEventWriterIter++){
    delete *tEventWriterIter;
    *tEventWriterIter = 0;
  }
  delete fEventWriterCollection;
}
//____________________________
AliFemtoManager& AliFemtoManager::operator=(const AliFemtoManager& aManager)
{
  // assignment operator
  if (this == &aManager)
    return *this;

  fEventReader = aManager.fEventReader;
  AliFemtoSimpleAnalysisIterator tAnalysisIter;
  if (fAnalysisCollection) {
    for (tAnalysisIter=fAnalysisCollection->begin();tAnalysisIter!=fAnalysisCollection->end();tAnalysisIter++){
      delete *tAnalysisIter;
      *tAnalysisIter = 0;
    }
    delete fAnalysisCollection;
  }
  // now delete each EventWriter in the Collection, and then the Collection itself
  AliFemtoEventWriterIterator tEventWriterIter;
  if (fEventWriterCollection) {
    for (tEventWriterIter=fEventWriterCollection->begin();tEventWriterIter!=fEventWriterCollection->end();tEventWriterIter++){
      delete *tEventWriterIter;
      *tEventWriterIter = 0;
    }
    delete fEventWriterCollection;
  }

  fAnalysisCollection = new AliFemtoAnalysisCollection;
  for (tAnalysisIter=aManager.fAnalysisCollection->begin();tAnalysisIter!=aManager.fAnalysisCollection->end();tAnalysisIter++){
    fAnalysisCollection->push_back(*tAnalysisIter);
  }

  fEventWriterCollection = new AliFemtoEventWriterCollection;
  for (tEventWriterIter=aManager.fEventWriterCollection->begin();tEventWriterIter!=aManager.fEventWriterCollection->end();tEventWriterIter++){
    fEventWriterCollection->push_back(*tEventWriterIter);
  }
  return *this;
}

//____________________________
int AliFemtoManager::Init(){
  // Execute initialization procedures
  AliFemtoString readerMessage;
  readerMessage += "*** *** *** *** *** *** *** *** *** *** *** *** \n";
  // EventReader
  if (fEventReader) {
    if (fEventReader->Init("r",readerMessage)){
      cout << " AliFemtoManager::Init() - Reader initialization failed " << endl;
      return (1);
    }
    readerMessage += fEventReader->Report();
  }
  // EventWriters
  AliFemtoEventWriterIterator tEventWriterIter;
  for (tEventWriterIter=fEventWriterCollection->begin();tEventWriterIter!=fEventWriterCollection->end();tEventWriterIter++){
    //cout << "*EventWriterIter " << *EventWriterIter << endl;
    // The message (AliFemtoString) passed into Init will be at the file header.
    // for that reason take the readerReport, add my own report and pass as message 
    AliFemtoString writerMessage = readerMessage;
    writerMessage += "*** *** *** *** *** *** *** *** *** *** *** *** \n";
    writerMessage += (*tEventWriterIter)->Report();
    if (*tEventWriterIter) {
      if ( (*tEventWriterIter)->Init("w",writerMessage)){ // yes, the message from the reader is passed into the writer
	cout << " AliFemtoManager::Init() - Writer initialization failed " << endl;
	return (1);
      }
    }
  }
  
  
  return (0);
}
//____________________________
void AliFemtoManager::Finish(){
  // Initialize finish procedures
  // EventReader
  if (fEventReader) fEventReader->Finish();
  // EventWriters
  AliFemtoEventWriterIterator tEventWriterIter;
  AliFemtoEventWriter* currentEventWriter;
  for (tEventWriterIter=fEventWriterCollection->begin();tEventWriterIter!=fEventWriterCollection->end();tEventWriterIter++){
    currentEventWriter = *tEventWriterIter;
    currentEventWriter->Finish();
  }
  // Analyses
  AliFemtoSimpleAnalysisIterator tAnalysisIter;
  AliFemtoAnalysis* currentAnalysis;
  for (tAnalysisIter=fAnalysisCollection->begin();tAnalysisIter!=fAnalysisCollection->end();tAnalysisIter++){
    currentAnalysis = *tAnalysisIter;
    currentAnalysis->Finish();
  }
}
//____________________________
AliFemtoString AliFemtoManager::Report(){
  // Construct a report from all the classes
  string stemp;
  char ctemp[100];
  // EventReader
  stemp = fEventReader->Report();
  // EventWriters
  sprintf(ctemp,"\nAliFemtoManager Reporting %u EventWriters\n",(unsigned int) fEventWriterCollection->size());
  stemp += ctemp;
  AliFemtoEventWriterIterator tEventWriterIter;
  AliFemtoEventWriter* currentEventWriter;
  for (tEventWriterIter=fEventWriterCollection->begin();tEventWriterIter!=fEventWriterCollection->end();tEventWriterIter++){
    //    cout << "AliFemtoManager - asking for EventWriter Report" << endl;
    currentEventWriter = *tEventWriterIter;
    stemp+=currentEventWriter->Report();
  }
  // Analyses
  sprintf(ctemp,"\nAliFemtoManager Reporting %u Analyses\n",(unsigned int) fAnalysisCollection->size());
  stemp += ctemp;
  AliFemtoSimpleAnalysisIterator tAnalysisIter;
  AliFemtoAnalysis* currentAnalysis;
  for (tAnalysisIter=fAnalysisCollection->begin();tAnalysisIter!=fAnalysisCollection->end();tAnalysisIter++){
    //    cout << "AliFemtoManager - asking for Analysis Report" << endl;
    currentAnalysis = *tAnalysisIter;
    stemp+=currentAnalysis->Report();
  }

  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
AliFemtoAnalysis* AliFemtoManager::Analysis( int n ){  // return pointer to n-th analysis
  // return analysis number n
  if ( n<0 || n > (int) fAnalysisCollection->size() )
    return NULL;
  AliFemtoSimpleAnalysisIterator iter = fAnalysisCollection->begin();
  for (int i=0; i<n ;i++){
    iter++;
  }
  return *iter;
}
//____________________________
AliFemtoEventWriter* AliFemtoManager::EventWriter( int n ){  // return pointer to n-th analysis
  // return event writern number n
  if ( n<0 || n > (int) fEventWriterCollection->size() )
    return NULL;
  AliFemtoEventWriterIterator iter = fEventWriterCollection->begin();
  for (int i=0; i<n ;i++){
    iter++;
  }
  return *iter;
}
 //____________________________
int AliFemtoManager::ProcessEvent(){
  // process a single event by reading it and passing it to each
  // analysis and event writer
  //  cout << "AliFemtoManager::ProcessEvent" << endl;
  // NOTE - this ReturnHbtEvent makes a *new* AliFemtoEvent - delete it when done!
  AliFemtoEvent* currentHbtEvent = fEventReader->ReturnHbtEvent();
  //  cout << "Event reader has returned control to manager" << endl;
  
  // if no HbtEvent is returned, then we abort processing.
  // the question is now: do we try again next time (i.e. there may be an HbtEvent next time)
  // or are we at EOF or something?  If Reader says Status=0, then that means try again later.
  // so, we just return the Reader's Status.
  if (!currentHbtEvent){
    cout << "AliFemtoManager::ProcessEvent() - Reader::ReturnHbtEvent() has returned null pointer\n";
    return fEventReader->Status();
  }
  
  // loop over all the EventWriters
  AliFemtoEventWriterIterator tEventWriterIter;
  for (tEventWriterIter=fEventWriterCollection->begin();tEventWriterIter!=fEventWriterCollection->end();tEventWriterIter++){
#ifdef STHBRDEBUG
    cout << " *tEventWriterIter " <<  *tEventWriterIter << endl;
#endif
    (*tEventWriterIter)->WriteHbtEvent(currentHbtEvent);
  } 

  // loop over all the Analysis
  AliFemtoSimpleAnalysisIterator tAnalysisIter;
  for (tAnalysisIter=fAnalysisCollection->begin();tAnalysisIter!=fAnalysisCollection->end();tAnalysisIter++){
    (*tAnalysisIter)->ProcessEvent(currentHbtEvent);
  } 

  if (currentHbtEvent) delete currentHbtEvent;
#ifdef STHBRDEBUG
  cout << "AliFemtoManager::ProcessEvent() - return to caller ... " << endl;
#endif
  return 0;    // 0 = "good return"
}       // ProcessEvent
