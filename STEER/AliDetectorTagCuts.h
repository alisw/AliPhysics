#ifndef ALIDETECTORTAGCUTS_H
#define ALIDETECTORTAGCUTS_H
/*  See cxx source for full Copyright notice */


/* $Id$ */

//-------------------------------------------------------------------------
//                       Class AliDetectorTagCuts
//              This is the class for the cuts in run tags
//
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------

#include <TObject.h>
#include <TString.h>

class AliDetectorTag;

//___________________________________________________________________________
class AliDetectorTagCuts : public TObject {
 public:
  AliDetectorTagCuts();
  ~AliDetectorTagCuts();
   
 //____________________________________________________//
  void SetListOfDetectors(const TString& detectors) {fDetectors = detectors; fDetectorsFlag = kTRUE;}
 
  Bool_t IsAccepted(AliDetectorTag *lhcTag) const;

  //____________________________________________________//
 private:
  Bool_t  IsSelected(TString detName, TString& detectors) const;

  TString fDetectors; //detectors active
  Bool_t  fDetectorsFlag; //cut used or not
  
  ClassDef(AliDetectorTagCuts, 1)
};

#endif
