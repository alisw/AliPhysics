#ifndef ALIMUONFACTORY_H
#define ALIMUONFACTORY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

////////////////////////////////////////////////////////////
//  Factory for muon chambers, segmentations and response //
////////////////////////////////////////////////////////////

#include <TObject.h>

class AliMUON;
class AliMUONResponseV0;

class AliMUONFactory : public  TObject {

  public:
    AliMUONFactory();
    virtual ~AliMUONFactory();
    
    void Build(AliMUON* where, const char* what);
    void BuildStation(AliMUON* where, Int_t stationNumber);

  protected:
    AliMUONFactory(const AliMUONFactory& rhs);
    AliMUONFactory& operator=(const AliMUONFactory& rhs);

  private:
    void BuildCommon();
    void BuildStation1();
    void BuildStation2();
    void BuildStation3();
    void BuildStation4();
    void BuildStation5();
    void BuildStation6();

    // data members	
    AliMUON*           fMUON;      // MUON detector 
    AliMUONResponseV0* fResponse0; // default response 

  ClassDef(AliMUONFactory,0)  // MUON Factory for Chambers and Segmentation
};
#endif















