/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup sim
/// \class AliMUONSegFactoryV3
/// \brief Factory for muon segmentations
///
/// New class separated from AliMUONFactoryV3 in order to get
/// building of segmentations independent from AliMUON and AliMUONChamber
/// objects

#ifndef ALI_MUON_SEG_FACTORY_V3_H
#define ALI_MUON_SEG_FACTORY_V3_H

#include <TNamed.h>

class AliMpExMap;

class AliMUONSegmentation;
class AliMUONGeometryTransformer;

class AliMUONSegFactoryV3 : public  TNamed {

  public:
    AliMUONSegFactoryV3(const char* name);
    AliMUONSegFactoryV3();
    virtual ~AliMUONSegFactoryV3();
    
    // Build methods
    void Build(const AliMUONGeometryTransformer* geometry);
    void BuildStation(const AliMUONGeometryTransformer*, Int_t stationNumber);
    
    // Access method
    AliMUONSegmentation* GetSegmentation() const;

  protected:
    AliMUONSegFactoryV3(const AliMUONSegFactoryV3& rhs);
    AliMUONSegFactoryV3& operator=(const AliMUONSegFactoryV3& rhs);

  private:
    Bool_t IsGeometryDefined(Int_t ichamber);
    Bool_t ReadDENames(const TString& fileName, AliMpExMap& map);

    void BuildChamber345(Int_t firstDetElemId, Int_t lastDetElemId,
                         const AliMpExMap& deNamesMap);

    void BuildStation1();
    void BuildStation2();
    void BuildStation3(const AliMpExMap& deNamesMap);
    void BuildStation4(const AliMpExMap& deNamesMap);
    void BuildStation5(const AliMpExMap& deNamesMap);
    void BuildStation6();
    
    // data members	
    AliMUONSegmentation*  fSegmentation;   // Segmentation container 
    const AliMUONGeometryTransformer* fkGeomTransformer; // Geometry parametrisation

  ClassDef(AliMUONSegFactoryV3,0)  // MUON Factory for Chambers and Segmentation
};

// inline functions

inline AliMUONSegmentation* AliMUONSegFactoryV3::GetSegmentation() const
{ return fSegmentation; }


#endif //ALI_MUON_SEG_FACTORY_V3_H















