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
  /// \cond CLASSIMP
  ClassImp(AliFemtoManager);
  /// \endcond
#endif



//____________________________
AliFemtoManager::AliFemtoManager():
  fAnalysisCollection(nullptr),
  fEventReader(nullptr),
  fEventWriterCollection(nullptr)
{
  // default constructor
  fAnalysisCollection = new AliFemtoAnalysisCollection;
  fEventWriterCollection = new AliFemtoEventWriterCollection;
}
//____________________________
AliFemtoManager::AliFemtoManager(const AliFemtoManager& aManager):
  fAnalysisCollection(new AliFemtoAnalysisCollection),
  fEventReader(aManager.fEventReader),
  fEventWriterCollection(new AliFemtoEventWriterCollection)
{
  // copy constructor
  for (auto *analysis : *aManager.fAnalysisCollection) {
    fAnalysisCollection->push_back(analysis);
  }
  for (auto *writer : *aManager.fEventWriterCollection) {
    fEventWriterCollection->push_back(writer);
  }
}

//____________________________
AliFemtoManager::~AliFemtoManager()
{
  // destructor
  delete fEventReader;
  // now delete each Analysis in the Collection, and then the Collection itself
  for (auto *analysis : *fAnalysisCollection) {
    delete analysis;
  }
  delete fAnalysisCollection;
  // now delete each EventWriter in the Collection, and then the Collection itself
  for (auto *writer : *fEventWriterCollection) {
    delete writer;
  }
  delete fEventWriterCollection;
}
//____________________________
AliFemtoManager& AliFemtoManager::operator=(const AliFemtoManager& aManager)
{
  // assignment operator
  if (this == &aManager) {
    return *this;
  }

  fEventReader = aManager.fEventReader;


  for (auto *analysis : *fAnalysisCollection) {
    delete analysis;
  }
  fAnalysisCollection->clear();
  for (auto *analysis : *aManager.fAnalysisCollection) {
    fAnalysisCollection->push_back(analysis);
  }

  // now delete each EventWriter in the Collection, and then the Collection itself
  for (auto *writer : *fEventWriterCollection) {
    delete writer;
  }
  fEventWriterCollection->clear();
  for (auto *writer : *aManager.fEventWriterCollection) {
    fEventWriterCollection->push_back(writer);
  }

  return *this;
}

//____________________________
int AliFemtoManager::Init()
{
  // Execute initialization procedures
  AliFemtoString readerMessage;
  readerMessage += "*** *** *** *** *** *** *** *** *** *** *** *** \n";
  // EventReader
  if (fEventReader) {
    if (fEventReader->Init("r",readerMessage)){
      cout << " AliFemtoManager::Init() - Reader initialization failed " << endl;
      return 1;
    }
    readerMessage += fEventReader->Report();
  }
  // EventWriters
  for (auto *event_writer : *fEventWriterCollection) {
    //cout << "*EventWriterIter " << *EventWriterIter << endl;
    // The message (AliFemtoString) passed into Init will be at the file header.
    // for that reason take the readerReport, add my own report and pass as message
    AliFemtoString writerMessage = readerMessage;
    writerMessage += "*** *** *** *** *** *** *** *** *** *** *** *** \n";
    writerMessage += event_writer->Report();
    if (event_writer) {
      if (event_writer->Init("w",writerMessage)) { // yes, the message from the reader is passed into the writer
        cout << " AliFemtoManager::Init() - Writer initialization failed\n";
        return 1;
      }
    }
  }
  return 0;
}
//____________________________
void AliFemtoManager::Finish()
{
  // Initialize finish procedures
  // EventReader
  if (fEventReader) {
    fEventReader->Finish();
  }

  // EventWriters
  for (auto *event_writer : *fEventWriterCollection) {
    event_writer->Finish();
  }
  // Analyses
  for (AliFemtoAnalysis *analysis : *fAnalysisCollection) {
    analysis->Finish();
  }
}
//____________________________
AliFemtoString AliFemtoManager::Report()
{
  // Construct a report from all the classes
  AliFemtoString report;

  // EventReader
  report = fEventReader->Report();

  // EventWriters
  report += Form("\nAliFemtoManager Reporting %lu EventWriters\n", fEventWriterCollection->size());
  for (auto *writer : *fEventWriterCollection) {
    report += writer->Report();
  }

  // Analyses
  report += Form("\nAliFemtoManager Reporting %lu Analyses\n", fAnalysisCollection->size());
  for (auto *analysis : *fAnalysisCollection) {
    report += analysis->Report();
  }

  return report;
}
//____________________________
AliFemtoAnalysis* AliFemtoManager::Analysis( int n )
{  // return pointer to n-th analysis
  if ( n < 0 || n > (int) fAnalysisCollection->size() ) {
    return nullptr;
  }

  AliFemtoSimpleAnalysisIterator iter = fAnalysisCollection->begin();
  std::advance(iter, n);
  return *iter;
}
//____________________________
AliFemtoEventWriter* AliFemtoManager::EventWriter( int n )
{ /// return pointer to n-th event writer
  if ( n < 0 || n > (int) fEventWriterCollection->size() ) {
    return nullptr;
  }

  AliFemtoEventWriterIterator iter = fEventWriterCollection->begin();
  std::advance(iter, n);
  return *iter;
}
 //____________________________
int AliFemtoManager::ProcessEvent()
{
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
    return fEventReader->Status();
  }

  // loop over all the EventWriters
  for (auto *writer : *fEventWriterCollection) {
    writer->WriteHbtEvent(currentHbtEvent);
  }

  // loop over all the Analysis
  for (auto *analysis : *fAnalysisCollection) {
    analysis->ProcessEvent(currentHbtEvent);
  }

  if (currentHbtEvent) {
    delete currentHbtEvent;
    currentHbtEvent = NULL;
  }

  return 0;    // 0 = "good return"
}       // ProcessEvent
