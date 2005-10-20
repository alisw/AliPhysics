#ifndef ALIMUONFactoryV4_H
#define ALIMUONFactoryV4_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup sim
/// \class AliMUONFactoryV4
/// \brief Factory for muon chambers, segmentations and response
///
/// New class to create AliMUONTriggerSegmentationV2 objects for St 6 & 7  
/// (the rest as in V3).

////////////////////////////////////////////////////////////
//  Factory for muon chambers, segmentations and response 
////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TNamed.h>

class AliMUON;
class AliMUONResponseV0;

class AliMUONFactoryV4 : public  TNamed {

  public:
    AliMUONFactoryV4(const char* name);
    AliMUONFactoryV4();
    virtual ~AliMUONFactoryV4();
    
    void Build(AliMUON* where, const char* what);
    void BuildStation(AliMUON* where, Int_t stationNumber);

  protected:
    AliMUONFactoryV4(const AliMUONFactoryV4& rhs);
    AliMUONFactoryV4& operator=(const AliMUONFactoryV4& rhs);

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
    void BuildChamberTrigger(Int_t firstDetElemId, Int_t lastDetElemId);
    
    // data members	
    AliMUON*           fMUON;           // MUON detector 
    AliMUONResponseV0* fResponse0;      // default response 
    TObjArray*         fDESegmentations;// DE segmentations

  ClassDef(AliMUONFactoryV4,0)  // MUON Factory for Chambers and Segmentation
};

#endif















