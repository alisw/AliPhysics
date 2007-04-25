/***************************************************************************
 *
 * $Id$
 *
 * Author: Mike Lisa, Ohio State, lisa@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *   The Manager is the top-level object that coordinates activities
 *   and performs event, particle, and pair loops, and checks the
 *   various Cuts of the Analyses in its AnalysisCollection
 *
 ***************************************************************************
 *
 * $Log$
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.20  2001/06/21 19:15:46  laue
 * Modified fiels:
 *   CTH.h : new constructor added
 *   AliFemtoEvent, AliFemtoKink, AliFemtoTrack : constructors from the persistent
 *                                   (TTree) classes added
 *   AliFemtoLikeSignAnalysis : minor changes, for debugging
 *   AliFemtoTypes: split into different files
 * Added files: for the new TTree muDst's
 *   StExceptions.cxx StExceptions.h AliFemtoEnumeration.h
 *   AliFemtoHelix.h AliFemtoHisto.h AliFemtoString.h AliFemtoTFile.h
 *   AliFemtoTTreeEvent.cxx AliFemtoTTreeEvent.h AliFemtoTTreeKink.cxx
 *   AliFemtoTTreeKink.h AliFemtoTTreeTrack.cxx AliFemtoTTreeTrack.h
 *   AliFemtoTTreeV0.cxx AliFemtoTTreeV0.h AliFemtoVector.h
 *
 * Revision 1.19  2000/05/08 15:45:50  laue
 * Memory leak fixed. Current hbt event was not deleted
 *
 * Revision 1.18  2000/03/17 17:23:05  laue
 * Roberts new three particle correlations implemented.
 *
 * Revision 1.16  2000/02/26 19:04:52  laue
 * Some unnecessary includes removed.
 * StThreeVectorD replace by AliFemtoThreeVector.
 * AliFemtoCoulomb modified to compile without Root (ClassDef embraced into
 *   #ifdef __ROOT__  ..... #endif)
 * AliFemtoParticle now returns references (FourMomentum(),Helix(),
 *   DecayVertexPosiion())
 *
 * Revision 1.15  2000/02/18 21:32:24  laue
 * franksTrackCut changed. If mCharge is set to '0' there will be no cut
 * on charge. This is important for front-loaded cuts.
 *
 * copy constructor implemented for AliFemtoEvent, AliFemtoTrack and AliFemtoV0.
 *
 * franks1HistoD.cxx franks1HistoD.h franks2HistoD.cxx franks2HistoD.h
 * removed. We can now (CC5 on Solaris) use the versions (no D)
 *
 * Revision 1.14  2000/02/13 17:17:12  laue
 * Calls to the EventBegin() and EventEnd() functions implemented
 * The actual analysis is moved from AliFemtoManager to AliFemtoAnalysis
 *
 * Revision 1.13  2000/01/25 17:35:17  laue
 * I. In order to run the stand alone version of the AliFemtoMaker the following
 * changes have been done:
 * a) all ClassDefs and ClassImps have been put into #ifdef __ROOT__ statements
 * b) unnecessary includes of StMaker.h have been removed
 * c) the subdirectory AliFemtoMaker/doc/Make has been created including everything
 * needed for the stand alone version
 *
 * II. To reduce the amount of compiler warning
 * a) some variables have been type casted
 * b) some destructors have been declared as virtual
 *
 * Revision 1.12  1999/10/15 01:57:29  lisa
 * Important enhancement of AliFemtoMaker - implement Franks CutMonitors
 * ----------------------------------------------------------
 * This means 3 new files in Infrastructure area (CutMonitor),
 * several specific CutMonitor classes in the Cut area
 * and a new base class in the Base area (AliFemtoCutMonitor).
 * This means also changing all Cut Base class header files from .h to .h
 * so we have access to CutMonitor methods from Cint command line.
 * This last means
 * 1) files which include these header files are slightly modified
 * 2) a side benefit: the TrackCuts and V0Cuts no longer need
 * a SetMass() implementation in each Cut class, which was stupid.
 * Also:
 * -----
 * Include Franks AliFemtoAssociationReader
 * ** None of these changes should affect any user **
 *
 * Revision 1.11  1999/10/04 15:38:57  lisa
 * include Franks new accessor methods AliFemtoAnalysis::CorrFctn and AliFemtoManager::Analysis as well as McEnt example macro
 *
 * Revision 1.10  1999/09/17 22:38:02  lisa
 * first full integration of V0s into AliFemto framework
 *
 * Revision 1.9  1999/09/08 04:15:52  lisa
 * persistent microDST implementation tweaked to please fickle solaris details
 *
 * Revision 1.8  1999/09/05 02:58:11  lisa
 * add ASCII microDST reader/writer AND franksParticle cuts
 *
 * Revision 1.7  1999/09/04 04:41:01  lisa
 * AliFemtoEvent IO   --and--  AliFemtoEventWriter (microDST) method added to framework
 *
 * Revision 1.6  1999/09/03 22:39:15  lisa
 * Readers now MUST have Report() methods and MAY have WriteHbtEvent() methods
 *
 * Revision 1.5  1999/07/27 10:47:04  lisa
 * now works in dev on linux and solaris - mistake in deleting picoEvents fixed
 *
 * Revision 1.4  1999/07/26 16:21:26  lisa
 * always convert string to char when output - needed on solaris
 *
 * Revision 1.3  1999/07/22 18:49:10  lisa
 * Implement idea of Fabrice to not create and delete AliFemtoPair all the time
 *
 * Revision 1.2  1999/07/06 22:33:22  lisa
 * Adjusted all to work in pro and new - dev itself is broken
 *
 * Revision 1.1.1.1  1999/06/29 16:02:57  lisa
 * Installation of AliFemtoMaker
 *
 **************************************************************************/

#include "Infrastructure/AliFemtoManager.h"
//#include "Infrastructure/AliFemtoParticleCollection.h"
//#include "Base/AliFemtoTrackCut.h"
//#include "Base/AliFemtoV0Cut.h"
#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoManager)
#endif



//____________________________
AliFemtoManager::AliFemtoManager(){
  fAnalysisCollection = new AliFemtoAnalysisCollection;
  fEventWriterCollection = new AliFemtoEventWriterCollection;
  fEventReader = 0;
}
//____________________________
AliFemtoManager::~AliFemtoManager(){
  delete fEventReader;
  // now delete each Analysis in the Collection, and then the Collection itself
  AliFemtoAnalysisIterator AnalysisIter;
  for (AnalysisIter=fAnalysisCollection->begin();AnalysisIter!=fAnalysisCollection->end();AnalysisIter++){
    delete *AnalysisIter;
    *AnalysisIter = 0;
  }
  delete fAnalysisCollection;
  // now delete each EventWriter in the Collection, and then the Collection itself
  AliFemtoEventWriterIterator EventWriterIter;
  for (EventWriterIter=fEventWriterCollection->begin();EventWriterIter!=fEventWriterCollection->end();EventWriterIter++){
    delete *EventWriterIter;
    *EventWriterIter = 0;
  }
  delete fEventWriterCollection;
}
//____________________________
int AliFemtoManager::Init(){
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
  AliFemtoEventWriterIterator EventWriterIter;
  for (EventWriterIter=fEventWriterCollection->begin();EventWriterIter!=fEventWriterCollection->end();EventWriterIter++){
    //cout << "*EventWriterIter " << *EventWriterIter << endl;
    // The message (AliFemtoString) passed into Init will be at the file header.
    // for that reason take the readerReport, add my own report and pass as message 
    AliFemtoString writerMessage = readerMessage;
    writerMessage += "*** *** *** *** *** *** *** *** *** *** *** *** \n";
    writerMessage += (*EventWriterIter)->Report();
    if (*EventWriterIter) {
      if ( (*EventWriterIter)->Init("w",writerMessage)){ // yes, the message from the reader is passed into the writer
	cout << " AliFemtoManager::Init() - Writer initialization failed " << endl;
	return (1);
      }
    }
  }
  
  
  return (0);
}
//____________________________
void AliFemtoManager::Finish(){
  // EventReader
  if (fEventReader) fEventReader->Finish();
  // EventWriters
  AliFemtoEventWriterIterator EventWriterIter;
  AliFemtoEventWriter* currentEventWriter;
  for (EventWriterIter=fEventWriterCollection->begin();EventWriterIter!=fEventWriterCollection->end();EventWriterIter++){
    currentEventWriter = *EventWriterIter;
    currentEventWriter->Finish();
  }
  // Analyses
  AliFemtoAnalysisIterator AnalysisIter;
  AliFemtoBaseAnalysis* currentAnalysis;
  for (AnalysisIter=fAnalysisCollection->begin();AnalysisIter!=fAnalysisCollection->end();AnalysisIter++){
    currentAnalysis = *AnalysisIter;
    currentAnalysis->Finish();
  }
}
//____________________________
AliFemtoString AliFemtoManager::Report(){
  string stemp;
  char ctemp[100];
  // EventReader
  stemp = fEventReader->Report();
  // EventWriters
  sprintf(ctemp,"\nAliFemtoManager Reporting %u EventWriters\n",fEventWriterCollection->size());
  stemp += ctemp;
  AliFemtoEventWriterIterator EventWriterIter;
  AliFemtoEventWriter* currentEventWriter;
  for (EventWriterIter=fEventWriterCollection->begin();EventWriterIter!=fEventWriterCollection->end();EventWriterIter++){
    cout << "AliFemtoManager - asking for EventWriter Report" << endl;
    currentEventWriter = *EventWriterIter;
    stemp+=currentEventWriter->Report();
  }
  // Analyses
  sprintf(ctemp,"\nAliFemtoManager Reporting %u Analyses\n",fAnalysisCollection->size());
  stemp += ctemp;
  AliFemtoAnalysisIterator AnalysisIter;
  AliFemtoBaseAnalysis* currentAnalysis;
  for (AnalysisIter=fAnalysisCollection->begin();AnalysisIter!=fAnalysisCollection->end();AnalysisIter++){
    cout << "AliFemtoManager - asking for Analysis Report" << endl;
    currentAnalysis = *AnalysisIter;
    stemp+=currentAnalysis->Report();
  }

  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
AliFemtoBaseAnalysis* AliFemtoManager::Analysis( int n ){  // return pointer to n-th analysis
  if ( n<0 || n > (int) fAnalysisCollection->size() )
    return NULL;
  AliFemtoAnalysisIterator iter = fAnalysisCollection->begin();
  for (int i=0; i<n ;i++){
    iter++;
  }
  return *iter;
}
//____________________________
AliFemtoEventWriter* AliFemtoManager::EventWriter( int n ){  // return pointer to n-th analysis
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
  cout << "AliFemtoManager::ProcessEvent" << endl;
  // NOTE - this ReturnHbtEvent makes a *new* AliFemtoEvent - delete it when done!
  AliFemtoEvent* currentHbtEvent = fEventReader->ReturnHbtEvent();
  cout << "Event reader has returned control to manager" << endl;
  
  // if no HbtEvent is returned, then we abort processing.
  // the question is now: do we try again next time (i.e. there may be an HbtEvent next time)
  // or are we at EOF or something?  If Reader says Status=0, then that means try again later.
  // so, we just return the Reader's Status.
  if (!currentHbtEvent){
    cout << "AliFemtoManager::ProcessEvent() - Reader::ReturnHbtEvent() has returned null pointer\n";
    return fEventReader->Status();
  }
  
  // loop over all the EventWriters
  AliFemtoEventWriterIterator EventWriterIter;
  for (EventWriterIter=fEventWriterCollection->begin();EventWriterIter!=fEventWriterCollection->end();EventWriterIter++){
#ifdef STHBRDEBUG
    cout << " *EventWriterIter " <<  *EventWriterIter << endl;
#endif
    (*EventWriterIter)->WriteHbtEvent(currentHbtEvent);
  } 

  // loop over all the Analysis
  AliFemtoAnalysisIterator AnalysisIter;
  for (AnalysisIter=fAnalysisCollection->begin();AnalysisIter!=fAnalysisCollection->end();AnalysisIter++){
    (*AnalysisIter)->ProcessEvent(currentHbtEvent);
  } 

  if (currentHbtEvent) delete currentHbtEvent;
#ifdef STHBRDEBUG
  cout << "AliFemtoManager::ProcessEvent() - return to caller ... " << endl;
#endif
  return 0;    // 0 = "good return"
}       // ProcessEvent
