/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup sim
/// \class AliMUONResponseFactory
/// \brief Factory for muon response
///
/// Factory for building response, separated from AliMUONFactoryV4

#ifndef ALI_MUON_RESPONSE_FACTORY_H
#define ALI_MUON_RESPONSE_FACTORY_H

#include <TObject.h>
#include <TNamed.h>

class AliMUON;
class AliMUONSegmentation;
class AliMUONResponseV0;

class AliMUONResponseFactory : public  TNamed {

  public:
    AliMUONResponseFactory(const char* name);
    AliMUONResponseFactory();
    virtual ~AliMUONResponseFactory();
    
    void Build(AliMUON* where);
    void BuildStation(AliMUON* where, Int_t stationNumber);

  private:
    AliMUONResponseFactory(const AliMUONResponseFactory& rhs);
    AliMUONResponseFactory& operator=(const AliMUONResponseFactory& rhs);

    void BuildCommon();
    void BuildStation1();
    void BuildStation2();
    void BuildStation3();
    void BuildStation4();
    void BuildStation5();
    void BuildStation6();
    
    // data members	
    AliMUON*             fMUON;           ///< MUON detector 
    AliMUONResponseV0*   fResponse0;      ///< default response 

  ClassDef(AliMUONResponseFactory,0)  // MUON Factory for Chambers and Segmentation
};

#endif //ALI_MUON_RESPONSE_FACTORY_H















