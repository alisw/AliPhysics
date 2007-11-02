////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoEventReaderESDKine - the reader class for the Alice ESD            //
// Reads in ESD information and converts it into internal AliFemtoEvent       //
// Reads in AliESDfriend to create shared hit/quality information             //
// Reads in Kine information and stores it in the hidden info                 //
// Authors: Marek Chojnacki mchojnacki@knf.pw.edu.pl                          //
//          Adam Kisiel kisiel@mps.ohio-state.edu                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

/*
 *$Id$
 *$Log$
 *Revision 1.1  2007/05/25 12:42:54  akisiel
 *Adding a reader for the Kine information
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
  

#ifndef ALIFEMTOEVENTREADERESDKINE_H
#define ALIFEMTOEVENTREADERESDKINE_H
#include "AliFemtoEventReader.h"
#include "AliFemtoEnumeration.h"

#include <string>
#include <vector>
#include "TChain.h"
#include "AliESDEvent.h"
#include <list>
#include "AliRunLoader.h"
#include "AliFemtoModelHiddenInfo.h"

class AliFemtoEvent;

class AliFemtoEventReaderESDKine : public AliFemtoEventReader 
{
 public:
  AliFemtoEventReaderESDKine();
  AliFemtoEventReaderESDKine(const AliFemtoEventReaderESDKine &aReader);
  ~AliFemtoEventReaderESDKine();

  AliFemtoEventReaderESDKine& operator=(const AliFemtoEventReaderESDKine& aReader);

  AliFemtoEvent* ReturnHbtEvent();
  AliFemtoString Report() const;
  //void SetFileName(const char* fileName);
  void SetInputFile(const char* inputFile);
  void SetConstrained(const bool constrained);
  bool GetConstrained() const;

 protected:

 private:
  TString        fInputFile;        // name of input file with ESD filenames
  TString        fFileName;         // name of current ESD file
  bool           fConstrained;      // flag to set which momentum from ESD file will be use
  int            fNumberofEvent;    // number of Events in ESD file
  int            fCurEvent;         // number of current event
  int            fCurRLEvent;       // Current simulated event
  TChain*        fTree;             // ESD tree
  AliESDEvent*   fEvent;            // ESD event
  AliRunLoader*  fRunLoader;        // Run loader for kine reading 
		
#ifdef __ROOT__
  ClassDef(AliFemtoEventReaderESDKine, 1)
#endif

    };
  
#endif


