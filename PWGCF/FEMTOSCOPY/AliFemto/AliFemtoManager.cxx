///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoManager: main class managing femtoscopic analysis             //
// The Manager is the top-level object that coordinates activities       //
// and performs event, particle, and pair loops, and checks the          //
// various Cuts of the Analyses in its AnalysisCollection                //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "AliFemtoManager.h"
#include <TString.h>

//____________________________
AliFemtoManager::AliFemtoManager():
  fAnalysisCollection(),
  fEventReader(nullptr),
  fEventWriterCollection()
{
  // default constructor
}
//____________________________
AliFemtoManager::AliFemtoManager(const AliFemtoManager& aManager):
  fAnalysisCollection(aManager.fAnalysisCollection),
  fEventReader(aManager.fEventReader),
  fEventWriterCollection(aManager.fEventWriterCollection)
{
  // copy constructor
}
//____________________________
AliFemtoManager::~AliFemtoManager()
{
  // destructor
  delete fEventReader;

  for (auto analysis : fAnalysisCollection) {
    delete analysis;
  }

  for (auto writer : fEventWriterCollection) {
    delete writer;
  }
}
//____________________________
AliFemtoManager& AliFemtoManager::operator=(const AliFemtoManager& aManager)
{
  // assignment operator
  if (this == &aManager) {
    return *this;
  }

  fEventReader = aManager.fEventReader;

  // delete old analyses, then copy the collection
  for (auto analysis : fAnalysisCollection) {
    delete analysis;
  }
  fAnalysisCollection = aManager.fAnalysisCollection;

  // delete old event_writers, then copy the collection
  for (auto event_writer : fEventWriterCollection) {
    delete event_writer;
  }
  fEventWriterCollection = aManager.fEventWriterCollection;

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
    if (fEventReader->Init("r", readerMessage)) {
      cout << " AliFemtoManager::Init() - Reader initialization failed " << endl;
      return 1;
    }
    readerMessage += fEventReader->Report();
  }
  // EventWriters
  for (auto event_writer : fEventWriterCollection) {
    // The message (AliFemtoString) passed into Init will be at the file header.
    // for that reason take the readerReport, add my own report and pass as message
    AliFemtoString writerMessage = readerMessage;
    writerMessage += "*** *** *** *** *** *** *** *** *** *** *** *** \n";
    writerMessage += event_writer->Report();
    if (event_writer->Init("w", writerMessage)) { // yes, the message from the reader is passed into the writer
      cout << " AliFemtoManager::Init() - Writer initialization failed " << endl;
      return 1;
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
  for (auto event_writer : fEventWriterCollection) {
    event_writer->Finish();
  }

  // Analyses
  for (auto analysis : fAnalysisCollection) {
    analysis->Finish();
  }
}
//____________________________
AliFemtoString AliFemtoManager::Report()
{
  // Construct a report from all the classes
  TString report;

  // EventReader
  report = fEventReader->Report();

  // EventWriters
  report += TString::Format("\nAliFemtoManager Reporting %lu EventWriters\n", fEventWriterCollection.size());
  for (auto event_writer : fEventWriterCollection) {
    report += event_writer->Report();
  }

  // Analyses
  report += TString::Format("\nAliFemtoManager Reporting %lu Analyses\n", fAnalysisCollection.size());
  for (auto analysis : fAnalysisCollection) {
    report += analysis->Report();
  }

  return AliFemtoString(report);
}
//____________________________
AliFemtoAnalysis* AliFemtoManager::Analysis(int n)
{  // return pointer to n-th analysis
  if (n < 0 || static_cast<size_t>(n) > fAnalysisCollection.size()) {
    return nullptr;
  }

  AliFemtoSimpleAnalysisIterator iter = fAnalysisCollection.begin();
  std::advance(iter, n);
  return *iter;
}
//____________________________
AliFemtoEventWriter* AliFemtoManager::EventWriter(int n)
{ /// return pointer to n-th event writer
  if (n < 0 || static_cast<size_t>(n) > fEventWriterCollection.size()) {
    return nullptr;
  }
  AliFemtoEventWriterIterator iter = std::begin(fEventWriterCollection);
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
  if (currentHbtEvent == nullptr) {
    return fEventReader->Status();
  }

  // loop over all the EventWriters
  for (auto event_writer : fEventWriterCollection) {
    event_writer->WriteHbtEvent(currentHbtEvent);
  }

  // loop over all the Analysis
  for (auto analysis : fAnalysisCollection) {
    analysis->ProcessEvent(currentHbtEvent);
  }

  delete currentHbtEvent;

  return 0;    // 0 = "good return"
}
