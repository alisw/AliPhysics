////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoEventReaderAOD - the reader class for the Alice AOD                //
// Reads in AOD information and converts it into internal AliFemtoEvent       //
// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOEVENTREADERAOD_H
#define ALIFEMTOEVENTREADERAOD_H
#include "AliFemtoEventReader.h"
#include "AliFemtoEnumeration.h"

#include <string>
#include <vector>
#include "TTree.h"
#include "TChain.h"
#include "TBits.h"
#include "AliAODEvent.h"
#include <list>
#include "AliPWG2AODTrack.h"

class AliFemtoEvent;
class AliFemtoTrack;

class AliFemtoEventReaderAOD : public AliFemtoEventReader 
{
 public:
  AliFemtoEventReaderAOD();
  AliFemtoEventReaderAOD(const AliFemtoEventReaderAOD &aReader);
  virtual ~AliFemtoEventReaderAOD();

  AliFemtoEventReaderAOD& operator=(const AliFemtoEventReaderAOD& aReader);

  virtual AliFemtoEvent* ReturnHbtEvent();
  AliFemtoString Report();
  void SetInputFile(const char* inputFile);
  void SetFilterBit(UInt_t ibit);

 protected:
  virtual void CopyAODtoFemtoEvent(AliFemtoEvent *tEvent);
  virtual void CopyAODtoFemtoTrack(const AliAODTrack *tAodTrack, 
				   AliFemtoTrack *tFemtoTrack, 
				   AliPWG2AODTrack *tPWG2AODTrack);

  int            fNumberofEvent;    // number of Events in AOD file
  int            fCurEvent;         // number of current event
  AliAODEvent*   fEvent;            // AOD event
  TBits          fAllTrue;          // Bit set with all true bits
  TBits          fAllFalse;         // Bit set with all false bits
  UInt_t         fFilterBit;        // Bitmap bit for AOD filters
  TClonesArray*  fPWG2AODTracks;    // Link to PWG2 specific AOD information (if it exists)

 private:

  string         fInputFile;        // name of input file with AOD filenames
  string         fFileName;         // name of current AOD file
  TChain*        fTree;             // AOD tree
  TFile*         fAodFile;          // AOD file 

#ifdef __ROOT__
  ClassDef(AliFemtoEventReaderAOD, 11)
#endif

    };
  
#endif


