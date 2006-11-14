/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup shuttle
/// \class AliMUONGMSSubprocessor
/// \brief The shuttle subprocessor for GMS data
///
/// \author Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_GMS_SUBPROCESSOR_H
#define ALI_MUON_GMS_SUBPROCESSOR_H

#ifndef ALIMUONVSUBPROCESSOR_H
#  include "AliMUONVSubprocessor.h"
#endif

#ifndef ALI_MUON_GEOMETRY_TRANSFORMER_H
  #include "AliMUONGeometryTransformer.h"
#endif

class AliMUONPreprocessor; 

class AliMUONGMSSubprocessor : public AliMUONVSubprocessor
{
  public:
    AliMUONGMSSubprocessor(AliMUONPreprocessor* master);
    virtual ~AliMUONGMSSubprocessor();

    // methods
    virtual UInt_t Process(TMap* /*dcsAliasMap*/);

  private:
    // static data members
    static const Int_t    fgkSystem;
    static const TString  fgkDataId;
    static const TString  fgkMatrixArrayName;
  
    // data members
    AliMUONGeometryTransformer fTransformer;///< Geometry transformer (used to get vo

    ClassDef(AliMUONGMSSubprocessor, 1); // Shuttle sub-processor for GMS
};

#endif
