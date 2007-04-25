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
 * Revision 1.10  2000/03/17 17:23:05  laue
 * Roberts new three particle correlations implemented.
 *
 * Revision 1.9  2000/02/18 21:32:24  laue
 * franksTrackCut changed. If mCharge is set to '0' there will be no cut
 * on charge. This is important for front-loaded cuts.
 *
 * copy constructor implemented for AliFemtoEvent, AliFemtoTrack and AliFemtoV0.
 *
 * franks1HistoD.cxx franks1HistoD.h franks2HistoD.cxx franks2HistoD.h
 * removed. We can now (CC5 on Solaris) use the versions (no D)
 *
 * Revision 1.8  2000/01/25 17:35:17  laue
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
 * Revision 1.7  1999/10/04 15:38:58  lisa
 * include Franks new accessor methods AliFemtoAnalysis::CorrFctn and AliFemtoManager::Analysis as well as McEvent example macro
 *
 * Revision 1.6  1999/09/24 01:23:12  fisyak
 * Reduced Include Path
 *
 * Revision 1.5  1999/09/08 04:15:52  lisa
 * persistent microDST implementation tweaked to please fickle solaris details
 *
 * Revision 1.4  1999/09/05 02:58:12  lisa
 * add ASCII microDST reader/writer AND franksParticle cuts
 *
 * Revision 1.3  1999/09/04 04:41:02  lisa
 * AliFemtoEvent IO   --and--  AliFemtoEventWriter (microDST) method added to framework
 *
 * Revision 1.2  1999/07/06 22:33:22  lisa
 * Adjusted all to work in pro and new - dev itself is broken
 *
 * Revision 1.1.1.1  1999/06/29 16:02:57  lisa
 * Installation of AliFemtoMaker
 *
 **************************************************************************/

#ifndef AliFemtoManager_hh
#define AliFemtoManager_hh


#include "Infrastructure/AliFemtoTypes.h"
#include "Infrastructure/AliFemtoAnalysisCollection.h"
#include "Infrastructure/AliFemtoEventWriterCollection.h"
#include "Infrastructure/AliFemtoEvent.h"
#include "Base/AliFemtoBaseAnalysis.h"
#include "Base/AliFemtoEventReader.h"
#include "Base/AliFemtoEventWriter.h"

class AliFemtoManager{

private:
  AliFemtoAnalysisCollection* fAnalysisCollection;
  AliFemtoEventReader*        fEventReader;
  AliFemtoEventWriterCollection* fEventWriterCollection;

public:
  AliFemtoManager();
  virtual ~AliFemtoManager();

  // Gets and Sets...
  AliFemtoAnalysisCollection* AnalysisCollection();
  AliFemtoBaseAnalysis* Analysis(int n);  // Access to Analysis within Collection
  void AddAnalysis(AliFemtoBaseAnalysis*);

  AliFemtoEventWriterCollection* EventWriterCollection();
  AliFemtoEventWriter* EventWriter(int n);// Access to EventWriter within Collection
  void SetEventWriter(AliFemtoEventWriter*);  // just for historic reasons
  void AddEventWriter(AliFemtoEventWriter*);

  AliFemtoEventReader* EventReader();
  void SetEventReader(AliFemtoEventReader*);


  int Init();
  int ProcessEvent();   // a "0" return value means success - otherwise quit
  void Finish();

  AliFemtoString Report(); //!
#ifdef __ROOT__
  ClassDef(AliFemtoManager, 0)
#endif
};

inline AliFemtoAnalysisCollection* AliFemtoManager::AnalysisCollection(){return fAnalysisCollection;}
inline void AliFemtoManager::AddAnalysis(AliFemtoBaseAnalysis* anal){fAnalysisCollection->push_back(anal);}

inline AliFemtoEventWriterCollection* AliFemtoManager::EventWriterCollection(){return fEventWriterCollection;}
inline void AliFemtoManager::AddEventWriter(AliFemtoEventWriter* writer){fEventWriterCollection->push_back(writer);}
inline void AliFemtoManager::SetEventWriter(AliFemtoEventWriter* writer){fEventWriterCollection->push_back(writer);}

inline AliFemtoEventReader* AliFemtoManager::EventReader(){return fEventReader;}
inline void AliFemtoManager::SetEventReader(AliFemtoEventReader* reader){fEventReader = reader;}


#endif

