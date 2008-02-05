///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoManager: main class managing femtoscopic analysis             //
// The Manager is the top-level object that coordinates activities       //
// and performs event, particle, and pair loops, and checks the          //
// various Cuts of the Analyses in its AnalysisCollection                //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOMANAGER_H
#define ALIFEMTOMANAGER_H


#include "AliFemtoTypes.h"
#include "AliFemtoAnalysisCollection.h"
#include "AliFemtoEventWriterCollection.h"
#include "AliFemtoEvent.h"
#include "AliFemtoAnalysis.h"
#include "AliFemtoEventReader.h"
#include "AliFemtoEventWriter.h"

class AliFemtoManager{

private:
  AliFemtoAnalysisCollection* fAnalysisCollection;       // Collection of analyzes
  AliFemtoEventReader*        fEventReader;              // Event reader
  AliFemtoEventWriterCollection* fEventWriterCollection; // Event writer collection

public:
  AliFemtoManager();
  AliFemtoManager(const AliFemtoManager& aManager);
  virtual ~AliFemtoManager();

  AliFemtoManager& operator=(const AliFemtoManager& aManager);

  // Gets and Sets...
  AliFemtoAnalysisCollection* AnalysisCollection();
  AliFemtoAnalysis* Analysis(int n);  // Access to Analysis within Collection
  void AddAnalysis(AliFemtoAnalysis* a);

  AliFemtoEventWriterCollection* EventWriterCollection();
  AliFemtoEventWriter* EventWriter(int n);// Access to EventWriter within Collection
  void SetEventWriter(AliFemtoEventWriter* w);  // just for historic reasons
  void AddEventWriter(AliFemtoEventWriter* w);

  AliFemtoEventReader* EventReader();
  void SetEventReader(AliFemtoEventReader* r);

  int Init();
  int ProcessEvent();   // a "0" return value means success - otherwise quit
  void Finish();

  AliFemtoString Report(); //!
#ifdef __ROOT__
  ClassDef(AliFemtoManager, 0)
#endif
};

inline AliFemtoAnalysisCollection* AliFemtoManager::AnalysisCollection(){return fAnalysisCollection;}
inline void AliFemtoManager::AddAnalysis(AliFemtoAnalysis* anal){fAnalysisCollection->push_back(anal);}

inline AliFemtoEventWriterCollection* AliFemtoManager::EventWriterCollection(){return fEventWriterCollection;}
inline void AliFemtoManager::AddEventWriter(AliFemtoEventWriter* writer){fEventWriterCollection->push_back(writer);}
inline void AliFemtoManager::SetEventWriter(AliFemtoEventWriter* writer){fEventWriterCollection->push_back(writer);}

inline AliFemtoEventReader* AliFemtoManager::EventReader(){return fEventReader;}
inline void AliFemtoManager::SetEventReader(AliFemtoEventReader* reader){fEventReader = reader;}


#endif

