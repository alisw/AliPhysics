///
/// \file AliFemto/AliFemtoEventReaderESD.h
///

#ifndef ALIFEMTOEVENTREADERESD_H
#define ALIFEMTOEVENTREADERESD_H

#include "AliFemtoEventReader.h"
#include "AliFemtoEnumeration.h"

#include <TFile.h>
#include <TChain.h>

#include <string>

class AliESDEvent;
class AliFemtoEvent;


/// \class AliFemtoEventReaderESD
/// \brief Event reader class for the Alice ESD
///
/// Reads in ESD information and converts it into internal AliFemtoEvent
/// Reads in AliESDfriend to create shared hit/quality information
///
/// ### History
///
/// * Revision 1.2.2.2  2007/10/04 13:10:52  akisiel
///   Add Kink index storageAliFemtoEventReaderESD.cxx AliFemtoTrack.cxx AliFemtoTrack.h
///
/// * Revision 1.1.2.1  2007/09/30 11:38:59  akisiel
///   Adapt the readers to the new AliESDEvent structure
///
/// * Revision 1.1  2007/05/16 10:22:11  akisiel
///   Making the directory structure of AliFemto flat. All files go into one common directory
///
/// * Revision 1.4  2007/05/03 09:45:20  akisiel
///   Fixing Effective C++ warnings
///
/// * Revision 1.3  2007/04/27 07:25:16  akisiel
///   Make revisions needed for compilation from the main AliRoot tree
///
/// * Revision 1.1.1.1  2007/04/25 15:38:41  panos
///   Importing the HBT code dir
///
/// \author Marek Chojnacki mchojnacki@knf.pw.edu.pl
/// \author Adam Kisiel kisiel@mps.ohio-state.edu
///
class AliFemtoEventReaderESD : public AliFemtoEventReader
{
public:
  AliFemtoEventReaderESD();
  AliFemtoEventReaderESD(const AliFemtoEventReaderESD &aReader);
  virtual ~AliFemtoEventReaderESD();

  AliFemtoEventReaderESD &operator=(const AliFemtoEventReaderESD &aReader);

  virtual AliFemtoEvent* ReturnHbtEvent();
  AliFemtoString Report();
  //void SetFileName(const char* fileName);

  /// Load rootfiles into TChain
  ///
  /// Opens specified file containing list of root filenames.
  /// Those containing "esdTree" TTree objects are added to fTree
  /// TChain member.
  ///
  void SetInputFile(const char *inputFile);

  /// Sets whether to use constrained momentum
  void SetConstrained(const bool constrained);
  bool GetConstrained() const;

  /// Sets whether to create an AliFemtoModelHiddenInfo from TPCInnerParam
  /// Note - pion data is hardcoded into PDG-code & mass info
  void SetReadTPCInner(const bool readinner);
  bool GetReadTPCInner() const;

protected:

private:
  /// name of input file with ESD filenames
  std::string fInputFile;

  /// name of current ESD file
  std::string fFileName;

  /// flag to set which momentum from ESD file will be use
  bool fConstrained;

  /// flag to set if one wants to read TPC-only momentum
  /// instead of the global one
  bool fReadInner;

  /// number of Events in ESD file
  int fNumberofEvent;

  /// number of current event
  int fCurEvent;

  /// ESD tree
  TChain *fTree;

  /// ESD file
  TFile *fEsdFile;

  /// ESD event
  AliESDEvent *fEvent;

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoEventReaderESD, 11);
  /// \endcond
#endif
};

#endif
