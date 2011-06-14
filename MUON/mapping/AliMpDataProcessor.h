/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup basic
/// \class AliMpDataProcessor
/// \brief Class for converting ASCII data files in the map of string
///
/// \author Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MP_DATA_PROCESSOR_H
#define ALI_MP_DATA_PROCESSOR_H

#include <TObject.h>
#include <TString.h>

#include <fstream>

class AliMpDataMap;

class AliMpDataProcessor : public TObject
{
  public:
    AliMpDataProcessor();
    virtual ~AliMpDataProcessor();
    
    // methods
    AliMpDataMap* CreateDataMap(const TString& dataDir = "data" );
    Bool_t        GenerateData(AliMpDataMap* dataMap,
                               const TString& outputDataDir = "data_new" );
    Bool_t        GenerateCode(AliMpDataMap* dataMap);
    
  private:
    // static methods
    static const TString& GetHeaderFileName();
    static const TString& GetImplFileName();  
 
    // methods  
    void ProcessDirectory(const TString& path, AliMpDataMap* map);
    void ProcessFile(const TString& path, AliMpDataMap* map );
    void GenerateFunction(const TString& path, const TString& data);
    void GenerateFileCode(const TString& path);
    void GenerateFill();

    // data
    Int_t     fCounter;    ///< data files counter
    ofstream  fHeaderFile; ///< header file
    ofstream  fImplFile;   ///< implementation file
 
  ClassDef(AliMpDataProcessor,0)  // Helper class for sorted integer array
};

#endif //ALI_MP_DATA_PROCESSOR_H

