/*
 *$Id$
 *$Log$
 *Revision 1.1.1.1  2007/04/25 15:38:41  panos
 *Importing the HBT code dir
 *
 *Revision 1.3  2007/03/13 15:30:03  mchojnacki
 *adding reader for simulated data
 *
 *Revision 1.2  2007/03/07 13:36:17  mchojnacki
 *Add some comments
 *
 *Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 *First version on CVS
 *
 */
  

#ifndef AliFemtoEventReaderESD_hh
#define AliFemtoEventReaderESD_hh
//Reader for ESD files for StHbt version 10 with hidden info part
//made by Marek Chojnacki mchojnacki@knf.pw.edu.pl
// Version 11 Added AliESDfriend reading
// Adam Kisiel kisiel@mps.ohio-state.edu
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
  ~AliFemtoEventReaderESD();
  AliFemtoEvent* ReturnHbtEvent();
  AliFemtoString Report();
  //void SetFileName(const char* fileName);
  void SetInputFile(const char* inputFile);
  void SetConstrained(const bool constrained);
  bool GetConstrained() const;

 protected:

 private:
  bool           GetNextFile();//setting next file to read 

  string         fInputFile; //name of input file
  string         fFileName; //name of current ESD file
  bool           fConstrained; //flag to set which momentum from ESD file will be use
  int            fNumberofEvent;//number of Events in ESD file
  int            fCurEvent; //number of current event
  unsigned int   fCurFile; //number of current file
  vector<string> fListOfFiles;//list of ESD files 		
  TTree*         fTree;//ESD tree
  AliESD*        fEvent;//ESD event
  TFile*         fEsdFile;//ESD file 
  AliESDfriend*  fEventFriend;

  list<Int_t>  **fSharedList;
  list<Int_t>  **fClusterPerPadrow;
		
#ifdef __ROOT__
  ClassDef(AliFemtoEventReaderESD, 10)
#endif

    };
  
#endif


