/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup base
/// \class AliMUONSegFactory
/// \brief New factory for building segmentations at all levels
///
/// The factory is associated with the AliMUONGeometryTransformer
/// object, used in geometry (module) segmentations for performing
/// trasformation between the global reference frame and the local DE one.
/// This transformer object can be passed by pointer or can be created 
/// by a factory and filled from the transformations data file.               
/// The transformer need not to be set if factory is used only
/// to create mapping segmentation:                                           \n
///
/// Construction:
/// - AliMUONSegFactory  factory(kTransformer);
/// - AliMUONSegFactory  factory("volpaths.dat", "transform.dat");
/// - AliMUONSegFactory  factory(0);                                          \n
///
/// All created objects are registered in the AliMUONSegmentation
/// object, which can be accessed via GetSegmentation() method.
/// A repetetive call to the same Create.. method does not create
/// a new object but returns the existing one.                                \n                
/// 
/// Factory does not delete the created segmentation objects.
/// They have to be deleted in the client code via the AliMUONSegmentation
/// container:                                                                \n
/// delete factory.GetSegmentation();
///
/// \author Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_SEG_FACTORY_H
#define ALI_MUON_SEG_FACTORY_H

#include "AliMpSegFactory.h"
#include "AliMpStringObjMap.h"

#include <TObject.h>

class AliMpVSegmentation;
class AliMUONVGeometryDESegmentation;
class AliMUONGeometrySegmentation;
class AliMUONSegmentation;
class AliMUONGeometryTransformer;

class AliMUONSegFactory : public  TObject {

  public:
    AliMUONSegFactory(const AliMUONGeometryTransformer* geometry);
    AliMUONSegFactory(const TString& volPathsFileName,
                      const TString& transformsFileName);
    AliMUONSegFactory();
    virtual ~AliMUONSegFactory();
    
    //
    // Build methods
    //
    
    AliMpVSegmentation*              
      CreateMpSegmentation(Int_t detElemId, Int_t cath);
              // Create mapping segmentation only 

    AliMUONSegmentation*  
      CreateSegmentation(const TString& option = "default"); 
              // Create segmentations on all levels and return their container.
    
    //
    // Get method
    //
    AliMUONSegmentation* GetSegmentation() const;
              // Returned segmentation contains all the lower level segmentations
	      // created with the factory

  protected:
    AliMUONSegFactory(const AliMUONSegFactory& rhs);
    AliMUONSegFactory& operator=(const AliMUONSegFactory& rhs);

  private:
    AliMUONVGeometryDESegmentation*  
      CreateDESegmentation(Int_t detElemId, Int_t cath);
              // Create DE segmentation, operating in local reference frame
    
    void
      CreateModuleSegmentations(Int_t chamberId, Int_t cath); 
              // Create module segmentation(s) for a given chamber, operating 
	      // in global reference frame

    // methods
    Bool_t IsGeometryDefined(Int_t ichamber);
    AliMUONSegmentation* Segmentation();
    
    // data members	
    AliMpSegFactory       fMpSegFactory;   ///< Mapping segmentation factory
    AliMpStringObjMap     fDESegmentations;///< Map of DE segmentations to DE names
    AliMUONSegmentation*  fSegmentation;   ///< Segmentation container 
    const AliMUONGeometryTransformer* fkTransformer; ///< Geometry transformer

  ClassDef(AliMUONSegFactory,0)  // MUON Factory for Chambers and Segmentation
};

/// Return segmentation
inline AliMUONSegmentation* AliMUONSegFactory::GetSegmentation() const
{ return fSegmentation; }

#endif //ALI_MUON_SEG_FACTORY_H















