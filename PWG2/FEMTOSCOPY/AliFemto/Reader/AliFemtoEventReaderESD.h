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
 *Revision 1.3  2007/04/27 07:25:16  akisiel
 *Make revisions needed for compilation from the main AliRoot tree
 *
 *Revision 1.1.1.1  2007/04/25 15:38:41  panos
 *Importing the HBT code dir
 *
 */
  

#ifndef AliFemtoEventReaderESD_hh
#define AliFemtoEventReaderESD_hh
#include "Base/AliFemtoEventReader.h"
#include "Infrastructure/AliFemtoEnumeration.h"

#include <string>
#include <vector>
#include "TTree.h"
#include "AliESD.h"
#include "AliESDfriend.h"
#include <list>

class AliFemtoEvent;

class AliFemtoEventReaderESD : public AliFemtoEventReader 
{
 public:
  AliFemtoEventReaderESD();
  AliFemtoEventReaderESD(const AliFemtoEventReaderESD &aReader);
  ~AliFemtoEventReaderESD();

  AliFemtoEventReaderESD& operator=(const AliFemtoEventReaderESD& aReader);

  AliFemtoEvent* ReturnHbtEvent();
  AliFemtoString Report();
  //void SetFileName(const char* fileName);
  void SetInputFile(const char* inputFile);
  void SetConstrained(const bool constrained);
  bool GetConstrained() const;

 protected:

 private:
  bool           GetNextFile();     // setting next file to read 

  string         fInputFile;        // name of input file with ESD filenames
  string         fFileName;         // name of current ESD file
  bool           fConstrained;      // flag to set which momentum from ESD file will be use
  int            fNumberofEvent;    // number of Events in ESD file
  int            fCurEvent;         // number of current event
  unsigned int   fCurFile;          // number of current file
  vector<string> fListOfFiles;      // list of ESD files 		
  TTree*         fTree;             // ESD tree
  AliESD*        fEvent;            // ESD event
  TFile*         fEsdFile;          // ESD file 
  AliESDfriend*  fEventFriend;      // ESD friend informaion

  list<Int_t>  **fSharedList;       //! Table (one list per padrow) of clusters which are shared
  list<Int_t>  **fClusterPerPadrow; //! Table (one list per padrow) of clusters in each padrow
		
#ifdef __ROOT__
  ClassDef(AliFemtoEventReaderESD, 10)
#endif

    };
  
#endif


