/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONGeometryDEIndexing
/// \brief Conversion between the detection element Id and the array index
///
/// The class that provides static methods for conversions
/// between the module & detection element Id
/// and the index in the array.
/// Used in storing DE transformations and segmentations.
/// The detection elements numbering:
///    DetElemId = chamberId*100 + detElemNum
///                where  chamberId  = 1, 2, ..., 14
///                       detElemNum = 0, 1, ...
///
/// Author: Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_GEOMETRY_DE_INDEXING_H
#define ALI_MUON_GEOMETRY_DE_INDEXING_H

class AliMUONGeometryDEIndexing 
{
  public:
    // static methods
    static  Int_t GetModuleId(Int_t detElemId);
    static  Int_t GetDEIndex(Int_t detElemId);
    static  Int_t GetDEId(Int_t moduleId, Int_t detElemIndex);

  private:
    AliMUONGeometryDEIndexing() {}
    ~AliMUONGeometryDEIndexing() {}
    
    // data members
    static  const Int_t fgkSeparator; 
};

#endif //ALI_MUON_GEOMETRY_DE_INDEXING_H
