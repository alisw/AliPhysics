#ifndef ALIFMDQADATAMAKERREC_H
#define ALIFMDQADATAMAKERREC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
#include "AliQADataMakerRec.h"
#include "TClonesArray.h"
class TH1F; 
class TH1I; 
class TList; 


//_____________________________________________________________________
// This class implements the AliQADataMakerRec for the FMD. Some
// functions are not implemented yet. 
// Author : Hans Hjersing Dalsgaard, hans.dalsgaard@cern.ch
//_____________________________________________________________________

class AliFMDQADataMakerRec: public AliQADataMakerRec 
{
public:
  AliFMDQADataMakerRec();
  AliFMDQADataMakerRec(const AliFMDQADataMakerRec& qadm);
  AliFMDQADataMakerRec& operator = (const AliFMDQADataMakerRec& qadm) ;
  virtual ~AliFMDQADataMakerRec();
private:
  virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX_t, TObjArray ** list);
  virtual void   InitESDs(); 
  virtual void   InitDigits(); 
  virtual void   InitRecPoints(); 
  virtual void   InitRaws(); 
  virtual void   MakeESDs(AliESDEvent * esd);
  virtual void   MakeDigits(TClonesArray * digits); 
  virtual void   MakeDigits(TTree * digitTree); 
  virtual void   MakeRecPoints(TTree * recpoTree); 
  virtual void   MakeRaws(AliRawReader* rawReader); 
  virtual void   StartOfDetectorCycle(); 
  Int_t GetHalfringIndex(UShort_t det, Char_t ring, UShort_t board, UShort_t monitor = 0);
  ClassDef(AliFMDQADataMakerRec,0)  // description 
  TClonesArray fDigitsArray;
  TClonesArray fRecPointsArray;

};

#endif // AliFMDQADataMakerRec_H
//____________________________________________________________________
//
// Local Variables: 
//  mode: c++
// End:
//

