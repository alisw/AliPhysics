#ifndef ALIMUONFACTORYV2_H
#define ALIMUONFACTORYV2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup sim
/// \class AliMUONFactoryV2
/// \brief Factory for muon chambers, segmentations and response
///
////////////////////////////////////////////////////////////
///  Factory for muon chambers, segmentations and response 
///  The number 2 is refering to new segmentation
///
////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TNamed.h>

class AliMUON;
class AliMUONResponseV0;

class AliMUONFactoryV2 : public  TNamed {

  public:
    AliMUONFactoryV2(const char* name);
    AliMUONFactoryV2();
    virtual ~AliMUONFactoryV2();
    
    void Build(AliMUON* where, const char* what);
    void BuildStation(AliMUON* where, Int_t stationNumber);

  protected:
    AliMUONFactoryV2(const AliMUONFactoryV2& rhs);
    AliMUONFactoryV2& operator=(const AliMUONFactoryV2& rhs);

  private:
    Bool_t IsGeometryDefined(Int_t ichamber);
    void BuildCommon();
    void BuildStation1();
    void BuildStation2();
    void BuildStation3();
    void BuildStation4();
    void BuildStation5();
    void BuildStation6();

    // data members	
    AliMUON*           fMUON;           // MUON detector 
    AliMUONResponseV0* fResponse0;      // default response 
    TObjArray*         fDESegmentations;// DE segmentations

  ClassDef(AliMUONFactoryV2,0)  // MUON Factory for Chambers and Segmentation
};
#endif















