#ifndef ALIOCDBTOOLKIT_H
#define ALIOCDBTOOLKIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */



#include <TObject.h>
class TStopwatch;
class TTree;
class TMemStatManager;
using std::fstream;

class AliOCDBtoolkit : public TObject {
public:  
  static void MakeDiffExampleUseCase();    // example usage
  static void DumpOCDBAsTxt(const TString fInput,const TString fType, const TString outfile);
  static void MakeSnapshotFromTxt(const TString fInput,const TString outfile, Bool_t singleKeys); 
  //
  static Bool_t   ParseInfoFromOcdbString(TString ocdbString, TString &ocdbPath, Int_t &run0, Int_t &run1, Int_t &version, Int_t &subVersion);   
  static Bool_t   ParseInfoFromOcdbString(TString ocdbString, AliCDBId &cdbId);
  static TList  * ConvertListStringToCDBId(const TList */*cdbList0*/);   
  //
  // Load OCDB entries 
  //
  static void SetStorage(const TMap *cdbMap0);   
  static void LoadOCDBFromMap(const TMap *cdbMap, const TList *cdbList);
  static void LoadOCDBFromLog(const char *logName, Int_t verbose);
  static void LoadOCDBFromESD(const char *fname="AliESDs.root");
  static void LoadOCDBFromList(const char */*ocdbList*/){;} // to be implemented  

  //
  // Dump object functionality
  //  
  static void DumpOCDB(const TMap *cdbMap0, const TList *cdbList0, const TString outfile);
  static void DumpObjectRecursive(TObject *obj);
  static void DumpObjectRecursive(TObject *obj, TString prefix, Int_t &counterRec);
  static void DumpOCDBFile(const char *finput , const char *foutput, Bool_t dumpMetaData, TString  printOption);
  static void MakeDiff(const TMap *cdbMap0, const TList *cdbList0, const TMap *cdbMap1, const TList *cdbList1, Int_t verbose);
  //
  // addopt OCDB entry
  //
  static Bool_t AddoptOCDBEntry( const char *finput, const char *output,  Int_t ustartRun, Int_t uendRun);
private:
  AliOCDBtoolkit(const AliOCDBtoolkit& source);
  AliOCDBtoolkit& operator= (const AliOCDBtoolkit& rec);

  ClassDef(AliOCDBtoolkit,0)
};

#endif
