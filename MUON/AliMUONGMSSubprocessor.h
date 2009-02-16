/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup shuttle
/// \class AliMUONGMSSubprocessor
/// \brief The shuttle subprocessor for GMS data
///
/// \author Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUONGMS_SUBPROCESSOR_H
#define ALI_MUONGMS_SUBPROCESSOR_H

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
    virtual void   Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
    virtual UInt_t Process(TMap* /*dcsAliasMap*/);

  private:
    /// Not implemented
    AliMUONGMSSubprocessor(const AliMUONGMSSubprocessor&);
    /// Not implemented
    AliMUONGMSSubprocessor& operator=(const AliMUONGMSSubprocessor&);

    UInt_t ProcessFile(const TString& filename);

    // static data members
    static const Int_t    fgkSystem;           ///< The data system
    static const TString  fgkDataId;           ///< The data Id
    static const TString  fgkMatrixArrayName;  ///< The fixed matrix array name
  
    // data members
    AliMUONGeometryTransformer* fTransformer; ///< Geometry transformer

    ClassDef(AliMUONGMSSubprocessor, 1) // Shuttle sub-processor for GMS
};

#endif //ALI_MUONGMS_SUBPROCESSOR_H
