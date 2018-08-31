#ifndef ALIFILEMERGER_H
#define ALIFILEMERGER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \class AliFileMerger
/// \brief AliFileMerger - Utilities for file merging
///
/// Utilities for file merging.
/// Additional functionality on top of the standard TFileMerger:
///
/// 1. Possibility to Set the reject/accept list.
///   1.a)  Only entries selected in accept list are merged. By default all entries are selected
/// use AddAccept 0 to specify your desired entry
///   1.b)  Entries selected in reject list are not merged. By default the reject list is empty.
///
/// 2. syswatch.log is created diring mergin procedure.
///   Memeory consumption - for reading and for merging can be monitored
///
///  RS: Changed merger to respect the structure of files being merged (directories, collections...)
///      Additional option: SetNoTrees (default false) to not merge any tree
///      The code mostly taken from root's hadd.cxx
///
/// Usage:
/// Libraries for all classes to be merged should be loaded before using the class
///  gSystem->Load("libANALYSIS");
///  gSystem->Load("libANALYSIScalib");
///  gSystem->Load("libTPCcalib"); 
///  TH1::AddDirectory(0);
///
/// Example usage starting from the input data list in text file:
///
///  AliFileMerger merger;
///  merger.AddReject("esdFriend");
///  merger.IterTXT("calib.list","CalibObjects.root",kFALSE);
/// \author marian.ivanov@cern.ch


class TObjString;

#include "TNamed.h"

class AliFileMerger : public TNamed
{

 public:
  AliFileMerger();
  AliFileMerger(const char* name);
  virtual ~AliFileMerger();
  void Merge(TFile* fileIn, TObjArray * array);

  void IterTXT( const char * fileList,  const char* outputFileName,Bool_t dontOverwrite=kFALSE);
  void IterAlien(const char* outputDir, const char* outputFileName = "CalibObjects.root" , const char* pattern = "AliESDfriends_v1.root", Bool_t dontOverwrite=kFALSE);
  void IterList(const TList* namesList, const char* outputFileName, Bool_t dontOverwrite=kFALSE);
  //
  void StoreResults(TObjArray * array, const char* outputFileName);
  void StoreSeparateResults(TObjArray * array, const char* outputFileName);
  //
  Bool_t IsAccepted(TString name);
  Bool_t IsRejected(TString name);
  void AddReject(const char *reject);
  void AddAccept(const char *accept);
  void SetNoTrees(Bool_t v=kTRUE)        {fNoTrees = v;}
  Bool_t IsNoTrees()               const {return fNoTrees;}
  void SetMaxFilesOpen(Int_t n)          {fMaxFilesOpen = n<3 ? 3 : n;}
  Int_t GetMaxFilesOpen()          const {return fMaxFilesOpen;}
  Bool_t      GetCheckTitle()      const {return fCheckTitle;}
  void        SetCheckTitle(Bool_t v=kTRUE) {fCheckTitle = v;}
  //
protected:
  int AddFile(TList* sourcelist, std::string entry);
  int MergeRootfile( TDirectory *target, TList *sourceNames, Bool_t nameFiltering=kTRUE);
  int OpenNextChunks(const TList* namesList, TList* filesList, Int_t from, Int_t to);
  void CheckTitle(TObject* tgt, TObject* src);
protected:
  TObjArray * fRejectMask;  ///< mask of the objects to be rejected
  TObjArray * fAcceptMask;    ///< mask of the objects to be accepted
  Int_t       fMaxFilesOpen;  ///< max number of open files
  Bool_t      fNoTrees;       ///< do we merge trees
  Bool_t      fCheckTitle;    ///< if source obj. title is empty, override it by eventual valid title
private:
  AliFileMerger(const AliFileMerger&);
  AliFileMerger& operator=(const AliFileMerger& other);
  
  ClassDef(AliFileMerger, 2); // File merger utilities for AliRoot
};

#endif
