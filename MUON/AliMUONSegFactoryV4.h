/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup sim
/// \class AliMUONSegFactoryV4
/// \brief Factory for muon segmentations
///
/// New class separated from AliMUONFactoryV4 in order to get
/// building of segmentations independent from AliMUON and AliMUONChamber
/// objects

#ifndef ALI_MUON_SEG_FACTORY_V4_H
#define ALI_MUON_SEG_FACTORY_V4_H

#include <TNamed.h>

class AliMpExMap;

class AliMUONSegmentation;
class AliMUONGeometryTransformer;

class AliMUONSegFactoryV4 : public  TNamed {

  public:
    AliMUONSegFactoryV4(const char* name);
    AliMUONSegFactoryV4();
    virtual ~AliMUONSegFactoryV4();
    
    // Build methods
    void Build(const AliMUONGeometryTransformer* geometry);
    void BuildStation(const AliMUONGeometryTransformer*, Int_t stationNumber);
    
    // Access method
    AliMUONSegmentation* GetSegmentation() const;

  protected:
    AliMUONSegFactoryV4(const AliMUONSegFactoryV4& rhs);
    AliMUONSegFactoryV4& operator=(const AliMUONSegFactoryV4& rhs);

  private:
    Bool_t IsGeometryDefined(Int_t ichamber);
    Bool_t ReadDENames(const TString& fileName, AliMpExMap& map);

    void BuildChamber345(Int_t firstDetElemId, Int_t lastDetElemId,
                         const AliMpExMap& deNamesMap);
    void BuildChamberTrigger(Int_t firstDetElemId, Int_t lastDetElemId,
                         const AliMpExMap& deNamesMap);

    void BuildStation1();
    void BuildStation2();
    void BuildStation3(const AliMpExMap& deNamesMap);
    void BuildStation4(const AliMpExMap& deNamesMap);
    void BuildStation5(const AliMpExMap& deNamesMap);
    void BuildStation6(const AliMpExMap& deNamesMap);
    
    // data members	
    AliMUONSegmentation* fSegmentation;   // Segmentation container 
    const AliMUONGeometryTransformer* fkGeomTransformer; // Geometry parametrisation

  ClassDef(AliMUONSegFactoryV4,0)  // MUON Factory for Chambers and Segmentation
};

// inline functions

inline AliMUONSegmentation* AliMUONSegFactoryV4::GetSegmentation() const
{ return fSegmentation; }

#endif //ALI_MUON_SEG_FACTORY_V4_H















