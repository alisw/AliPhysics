/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup geometry
/// \class AliMUONSegmentation
/// \brief Container class for modules segmentations
///
/// It provides access to segmentations on all levels:
/// - mapping segmentation
/// - DE segmentation (operating in local DE reference frame)
/// - module segmentation (operating in global reference frame)
///
/// As some detection elements are sharing the same objects
/// (AliMpVSegmentation, AliMUONVGeometryDESegmentation),
/// all segmentations objects have to be always deleted
/// altogether via deleting this container object. 
/// 
/// \author Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_SEGMENTATION_H
#define ALI_MUON_SEGMENTATION_H

#include <TObject.h>
#include <TGeoMatrix.h>

class TObjArray;

class AliMpVSegmentation;

class AliMUONGeometrySegmentation;
class AliMUONVGeometryDESegmentation;

class AliMUONSegmentation : public TObject
{
  public:
    AliMUONSegmentation(Int_t nofModules);
    AliMUONSegmentation();
    virtual  ~AliMUONSegmentation();
    
    // methods
    void  AddMpSegmentation(AliMpVSegmentation* segmentation);
    void  AddDESegmentation(AliMUONVGeometryDESegmentation* segmentation);

    void  AddModuleSegmentation(Int_t moduleId, Int_t cathod,
                            AliMUONGeometrySegmentation* segmentation);
    void  Init();
            // This function should not be needed;
	    // the segmentations should be built in a valid state
	    // To be revised			    

    //
    // get methods
    //
    
    // Geometry segmentations
    //

    AliMUONGeometrySegmentation* GetModuleSegmentationByDEId(
                     Int_t detElemId, Int_t cathod, Bool_t warn = true) const;

    // DE segmentations
    //
    const AliMUONVGeometryDESegmentation* GetDESegmentation(
                     Int_t detElemId, Int_t cathod, Bool_t warn = true) const;

    /** Mapping segmentations access by cathode number.
      cathod can be 0 or 1. Note that there's no trivial relationship
      between the cathod number and whether the corresponding plane
      is a Bending or NonBending one.
      **/
    const AliMpVSegmentation* GetMpSegmentation(
                     Int_t detElemId, Int_t cathod, Bool_t warn = true) const;
		         
    // DE properties
    //
    Bool_t   HasDE(Int_t detElemId, Int_t cathod = 0) const;
    TString  GetDEName(Int_t detElemId, Int_t cathod = 0) const;

  protected:
    AliMUONSegmentation(const AliMUONSegmentation& right);
    AliMUONSegmentation&  operator = (const AliMUONSegmentation& right);
     
  private:
    AliMUONGeometrySegmentation* GetModuleSegmentation(
                     Int_t moduleId, Int_t cathod, Bool_t warn = true) const;

    // data members
    TObjArray*  fMpSegmentations;        ///< array of mapping segmentations
    TObjArray*  fDESegmentations;        ///< array of DE segmentations
    TObjArray*  fModuleSegmentations[2]; ///< \brief array of module segmentations
                                         /// for two cathods         
  ClassDef(AliMUONSegmentation,2)  // Container class for module segmentations
};

#endif //ALI_MUON_SEGMENTATION_H







