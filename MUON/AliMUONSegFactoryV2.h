/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup sim
/// \class AliMUONSegFactoryV2
/// \brief Factory for muon segmentations
///
/// New class separated from AliMUONFactoryV2 in order to get
/// building of segmentations independent from AliMUON and AliMUONChamber
/// objects

#ifndef ALI_MUON_SEG_FACTORY_V2_H
#define ALI_MUON_SEG_FACTORY_V2_H

#include <TNamed.h>

class AliMUONSegmentation;
class AliMUONGeometryTransformer;

class AliMUONSegFactoryV2 : public  TNamed {

  public:
    AliMUONSegFactoryV2(const char* name);
    AliMUONSegFactoryV2();
    virtual ~AliMUONSegFactoryV2();
    
    // Build methods
    void Build(const AliMUONGeometryTransformer* geometry);
    void BuildStation(const AliMUONGeometryTransformer*, Int_t stationNumber);
    
    // Access method
    AliMUONSegmentation* GetSegmentation() const;

  protected:
    AliMUONSegFactoryV2(const AliMUONSegFactoryV2& rhs);
    AliMUONSegFactoryV2& operator=(const AliMUONSegFactoryV2& rhs);

  private:
    Bool_t IsGeometryDefined(Int_t ichamber);

    void BuildStation1();
    void BuildStation2();
    void BuildStation3();
    void BuildStation4();
    void BuildStation5();
    void BuildStation6();
    
    // data members	
    AliMUONSegmentation*  fSegmentation;   // Segmentation container 
    const AliMUONGeometryTransformer* fkGeomTransformer; // Geometry parametrisation

  ClassDef(AliMUONSegFactoryV2,0)  // MUON Factory for Chambers and Segmentation
};

// inline functions

inline AliMUONSegmentation* AliMUONSegFactoryV2::GetSegmentation() const
{ return fSegmentation; }


#endif //ALI_MUON_SEG_FACTORY_V3_H















