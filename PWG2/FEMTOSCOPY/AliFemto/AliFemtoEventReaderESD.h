////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoEventReaderESD - the reader class for the Alice ESD              ///
/// Reads in ESD information and converts it into internal AliFemtoEvent     ///
/// Reads in AliESDfriend to create shared hit/quality information           ///
/// Authors: Marek Chojnacki mchojnacki@knf.pw.edu.pl                        ///
///          Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////

/*
 *$Id$
 *$Log$
 *Revision 1.1.2.1  2007/09/30 11:38:59  akisiel
 *Adapt the readers to the new AliESDEvent structure
 *
 *Revision 1.1  2007/05/16 10:22:11  akisiel
 *Making the directory structure of AliFemto flat. All files go into one common directory
 *
 *Revision 1.4  2007/05/03 09:45:20  akisiel
 *Fixing Effective C++ warnings
 *
 *Revision 1.3  2007/04/27 07:25:16  akisiel
 *Make revisions needed for compilation from the main AliRoot tree
 *
 *Revision 1.1.1.1  2007/04/25 15:38:41  panos
 *Importing the HBT code dir
 *
 */
  

#ifndef ALIFEMTOEVENTREADERESD_H
#define ALIFEMTOEVENTREADERESD_H
#include "AliFemtoEventReader.h"
#include "AliFemtoEnumeration.h"

#include <string>
#include <vector>
#include "TTree.h"
#include "TChain.h"
#include "AliESDEvent.h"
#include <list>

class AliFemtoEvent;

class AliFemtoEventReaderESD : public AliFemtoEventReader 
{
 public:
  AliFemtoEventReaderESD();
  AliFemtoEventReaderESD(const AliFemtoEventReaderESD &aReader);
  ~AliFemtoEventReaderESD();

  AliFemtoEventReaderESD& operator=(const AliFemtoEventReaderESD& aReader);

  virtual AliFemtoEvent* ReturnHbtEvent();
  AliFemtoString Report();
  //void SetFileName(const char* fileName);
  void SetInputFile(const char* inputFile);
  void SetConstrained(const bool constrained);
  bool GetConstrained() const;
  void SetReadTPCInner(const bool readinner);
  bool GetReadTPCInner() const;

 protected:

 private:
  string         fInputFile;        // name of input file with ESD filenames
  string         fFileName;         // name of current ESD file
  bool           fConstrained;      // flag to set which momentum from ESD file will be use
  bool           fReadInner;        // flag to set if one wants to read TPC-only momentum
                                    // instead of the global one
  int            fNumberofEvent;    // number of Events in ESD file
  int            fCurEvent;         // number of current event
  TChain*        fTree;             // ESD tree
  TFile*         fEsdFile;          // ESD file 
  AliESDEvent*   fEvent;            // ESD event
		
#ifdef __ROOT__
  ClassDef(AliFemtoEventReaderESD, 11)
#endif

    };
  
#endif


