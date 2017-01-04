///
/// \file  AliFemtoManager.h
///

#pragma once

#ifndef ALIFEMTOMANAGER_H
#define ALIFEMTOMANAGER_H

#include "AliFemtoTypes.h"
#include "AliFemtoAnalysisCollection.h"
#include "AliFemtoEventWriterCollection.h"
#include "AliFemtoAnalysis.h"
#include "AliFemtoEventReader.h"
#include "AliFemtoEventWriter.h"


/// \class AliFemtoManager
/// \brief Main class for managing femtoscopic analyses
///
/// The Manager is the top-level object containing an EventReader
/// (the input), and collections of EventWriters and AliFemtoAnalyses
/// (the outputs).
///
/// A manager object is owned by an AliAnalysistaskFemto object which
/// calls the `ProcessEvent()` method, which reads an AliFemtoEvent
/// from the input files, and forwards it to the `ProcessEvent` method
/// in each output analysis, which is responsible for carrying out the
/// actual cuts & computation.
///
/// AliFemtoManager objects "own" the EventReader, Analyses, and
/// EventWriters added to them, and is responsible for deleting them
/// upon its own destruction.
///
/// AliFemtoManager objects are not copyable, as the AliFemtoAnalysis
/// objects they contain have no means of copying/cloning.
/// Denying copyability by making the copy constructor and assignment
/// operator private prevents potential dangling pointer (segfault)
/// errors.
///
class AliFemtoManager {

private:
  AliFemtoAnalysisCollection fAnalysisCollection;        ///< Collection of analyzes
  AliFemtoEventReader*        fEventReader;              ///< Event reader
  AliFemtoEventWriterCollection fEventWriterCollection;  ///< Event writer collection

  /// private to avoid memory bugs
  AliFemtoManager(const AliFemtoManager& aManager);

  /// private to avoid memory bugs
  AliFemtoManager& operator=(const AliFemtoManager& aManager);

public:
  AliFemtoManager();
  virtual ~AliFemtoManager();

  // Gets and Sets...
  AliFemtoAnalysisCollection* AnalysisCollection();
  AliFemtoAnalysis* Analysis(int n);            ///< Access to Analysis within Collection
  void AddAnalysis(AliFemtoAnalysis* a);

  AliFemtoEventWriterCollection* EventWriterCollection();
  AliFemtoEventWriter* EventWriter(int n);      ///< Access to EventWriter within Collection
  void SetEventWriter(AliFemtoEventWriter* w);  ///< just for historic reasons
  void AddEventWriter(AliFemtoEventWriter* w);

  AliFemtoEventReader* EventReader();
  void SetEventReader(AliFemtoEventReader* r);

  /// Calls `Init()` on all owned EventWriters
  ///
  /// Returns 0 for success, 1 for failure.
  ///
  int Init();

  /// Reads event from EventReader and passes to each EventWriters'
  /// `WriteHbtEvent` and each Analyses' `ProcessEvent` methods.
  ///
  /// Returns 0 if successful, otherwise the 'Status' of the
  /// EventReader after failing to return an AliFemtoEvent.
  ///
  int ProcessEvent();

  /// Calls `Finish()` on the EventReader, EventWriters, and the Analyses.
  void Finish();

  /// Returns combined report of all readers, writers, and analyses.
  AliFemtoString Report();
};

inline AliFemtoAnalysisCollection* AliFemtoManager::AnalysisCollection(){return &fAnalysisCollection;}
inline void AliFemtoManager::AddAnalysis(AliFemtoAnalysis* anal){fAnalysisCollection.push_back(anal);}

inline AliFemtoEventWriterCollection* AliFemtoManager::EventWriterCollection(){return &fEventWriterCollection;}
inline void AliFemtoManager::AddEventWriter(AliFemtoEventWriter* writer){fEventWriterCollection.push_back(writer);}
inline void AliFemtoManager::SetEventWriter(AliFemtoEventWriter* writer){fEventWriterCollection.push_back(writer);}

inline AliFemtoEventReader* AliFemtoManager::EventReader(){return fEventReader;}
inline void AliFemtoManager::SetEventReader(AliFemtoEventReader* reader){fEventReader = reader;}

#endif
