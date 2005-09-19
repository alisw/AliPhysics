#ifndef ALIMUONFACTORYV3_H
#define ALIMUONFACTORYV3_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////
//  Factory for muon chambers, segmentations and response 
////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TNamed.h>

class AliMUON;
class AliMUONResponseV0;

class AliMUONFactoryV3 : public  TNamed {

  public:
    AliMUONFactoryV3(const char* name);
    AliMUONFactoryV3();
    virtual ~AliMUONFactoryV3();
    
    void Build(AliMUON* where, const char* what);
    void BuildStation(AliMUON* where, Int_t stationNumber);

  protected:
    AliMUONFactoryV3(const AliMUONFactoryV3& rhs);
    AliMUONFactoryV3& operator=(const AliMUONFactoryV3& rhs);

  private:
    Bool_t IsGeometryDefined(Int_t ichamber);
    void BuildCommon();
    void BuildStation1();
    void BuildStation2();
    void BuildStation3();
    void BuildStation4();
    void BuildStation5();
    void BuildStation6();
    void BuildChamber345(Int_t firstDetElemId, Int_t lastDetElemId);

    // data members	
    AliMUON*           fMUON;           // MUON detector 
    AliMUONResponseV0* fResponse0;      // default response 
    TObjArray*         fDESegmentations;// DE segmentations

  ClassDef(AliMUONFactoryV3,0)  // MUON Factory for Chambers and Segmentation
};

#endif















