#ifndef ALIQUALASSDATAMAKER_H
#define ALIQUALASSDATAMAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

/*
  Base Class:
  Produces the data needed to calculate the quality assurance. 
  All data must be mergeable objects.
  Y. Schutz CERN July 2007
*/


// --- ROOT system ---
#include <TNamed.h>  
class TFile;  
class TDirectory;
class TObject; 
class TTree; 
class AliESDEvent;

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQualAss.h"

class AliQualAssDataMaker: public TNamed {
  
public:
  
  AliQualAssDataMaker(const char * name="", const char * title="") ;          // ctor
  AliQualAssDataMaker(const AliQualAssDataMaker& qadm) ;   
  AliQualAssDataMaker& operator = (const AliQualAssDataMaker& qadm) ;
  virtual ~AliQualAssDataMaker() {;} // dtor
  
  virtual void        Exec(AliQualAss::TASKINDEX, TObject * data) ;
  void                Finish(AliQualAss::TASKINDEX task) const ; 
  static const char * GetDetectorDirName() { return fDetectorDirName.Data() ; }
  void                Init(AliQualAss::TASKINDEX) ;
 
protected: 

  virtual void   InitDigits()        {AliInfo("To ne implemented by detectors");}
  virtual void   InitESDs()          {AliInfo("To ne implemented by detectors");}
  virtual void   InitHits()          {AliInfo("To ne implemented by detectors");}
  virtual void   InitRecParticles()  {AliInfo("To ne implemented by detectors");}
  virtual void   InitRecPoints()     {AliInfo("To ne implemented by detectors");}
  virtual void   InitSDigits()       {AliInfo("To ne implemented by detectors");}
  virtual void   InitTrackSegments() {AliInfo("To ne implemented by detectors");}
  virtual void   MakeESDs(AliESDEvent * )          {AliInfo("To ne implemented by detectors");} 
  virtual void   MakeHits(TObject * )              {AliInfo("To ne implemented by detectors");} 
  virtual void   MakeDigits(TObject * )            {AliInfo("To ne implemented by detectors");} 
  //  virtual void   MakeRecParticles(TClonesArray * ) {AliInfo("To ne implemented by detectors");} 
  virtual void   MakeRecPoints(TTree * )           {AliInfo("To ne implemented by detectors");} 
  virtual void   MakeSDigits(TObject * )           {AliInfo("To ne implemented by detectors");} 
  //virtual void   MakeTrackSegments(TTree * )       {AliInfo("To ne implemented by detectors");} 

  TFile *       fOutput ;      //! output root file
  TDirectory *  fDetectorDir ; //! directory for the given detector in the file
  static TString fDetectorDirName ; //! detector directory name in the quality assurance data file
  ClassDef(AliQualAssDataMaker,1)  // description 

};

#endif // AliQualAssDataMaker_H
