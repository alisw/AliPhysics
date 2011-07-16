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
#include <AliDAQ.h>

class AliDetectorTag;

//___________________________________________________________________________
class AliDetectorTagCuts : public TObject {
 public:
  AliDetectorTagCuts();
  ~AliDetectorTagCuts();
   
 //____________________________________________________//
  void SetListOfDetectors(const TString& detectors) {fDetectorsDAQ = AliDAQ::DetectorPattern(detectors); fDetectorsReco = AliDAQ::DetectorPattern(detectors); fDetectorsFlag = kTRUE;}
  void SetListOfDetectorsDAQ(const TString& detectors) {fDetectorsDAQ = AliDAQ::DetectorPattern(detectors); fDetectorsFlag = kTRUE;}
  void SetListOfDetectorsReco(const TString& detectors) {fDetectorsReco = AliDAQ::DetectorPattern(detectors); fDetectorsFlag = kTRUE;}
  void SetDetectorValidityValue(TString det, UShort_t val);
 
  Bool_t IsAccepted(AliDetectorTag *detTag) const;

  //____________________________________________________//
 private:
  //  Bool_t  IsSelected(TString detName, TString& detectors) const;

  //  TString fDetectors; //detectors active
  UInt_t fDetectorsReco;  //selected detector pattern for Reco
  UInt_t fDetectorsDAQ;   //selected detector pattern for DAQ
  Bool_t fDetectorsFlag; //cut used or not
  UShort_t   fDetectorValidityMatch[AliDAQ::kHLTId]; // Detector validity to match
  Bool_t     fDetectorValidityFlag[AliDAQ::kHLTId];  // Flag if validity match is to be used
  
  ClassDef(AliDetectorTagCuts, 3)
};

#endif
