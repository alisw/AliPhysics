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
  /** 
   * Constructor
   */
  AliFMDQADataMakerRec();
  /** 
   * Copy constructor 
   * 
   * @param qadm What to copy from
   */
  AliFMDQADataMakerRec(const AliFMDQADataMakerRec& qadm);
  /** 
   * Assignment operator 
   * 
   * @param qadm What to assign from 
   * 
   * @return Reference to this
   */
  AliFMDQADataMakerRec& operator = (const AliFMDQADataMakerRec& qadm) ;
  /** 
   * Destrcutor 
   */
  virtual ~AliFMDQADataMakerRec();
private:
  /** 
   * Called at end of monitor cycle 
   * 
   * @param TASKINDEX_t Task
   * @param list        Output list
   */
  virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX_t, TObjArray ** list);
  /** 
   * Intialize for ESD
   */
  virtual void   InitESDs(); 
  /** 
   * Intialize for Digits
   */
  virtual void   InitDigits(); 
  /** 
   * Intialize for RecPoints
   */
  virtual void   InitRecPoints(); 
  /** 
   * Initialise for raw 
   */
  virtual void   InitRaws(); 
  /** 
   * Analyse ESD event
   * 
   * @param esd ESD event
   */
  virtual void   MakeESDs(AliESDEvent * esd);
  /** 
   * Analyse digits 
   */
  virtual void   MakeDigits(); 
  /** 
   * Analyse digits
   * 
   * @param digitTree Tree of digits
   */
  virtual void   MakeDigits(TTree * digitTree); 
  /** 
   * Analyse rec points
   * 
   * @param recpoTree Tree of RecPoints
   */
  virtual void   MakeRecPoints(TTree * recpoTree); 
  /** 
   * Analyse raw 
   * 
   * @param rawReader Raw reader
   */
  virtual void   MakeRaws(AliRawReader* rawReader); 
  /** 
   * Called at start of a cycle 
   * 
   */
  virtual void   StartOfDetectorCycle(); 
  /** 
   * Get the half-ring index
   * 
   * @param det      Detector
   * @param ring     Ring
   * @param board    Board number
   * @param monitor  Monitor 
   * 
   * @return Half ring index
   */
  Int_t GetHalfringIndex(UShort_t det, Char_t ring, 
			 UShort_t board, UShort_t monitor = 0) const;
  ClassDef(AliFMDQADataMakerRec,0)  // description 
  TClonesArray fRecPointsArray; // Rec points
};

#endif // AliFMDQADataMakerRec_H
//____________________________________________________________________
//
// Local Variables: 
//  mode: c++
// End:
//

