#ifndef ALIESDTAGCREATOR_H
#define ALIESDTAGCREATOR_H
/*  See cxx source for full Copyright notice */


/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliESDTagCreator
//   This is the AliESDTagCreator class for the tag creation (post process)
//
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                        AliESDTagCreator                              //
//                                                                      //
//           Implementation of the tag creation mechanism.              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


//ROOT
#include <TSystem.h>
class TChain;
class TList;

#include <AliTagCreator.h>

class TFile;
class TGridResult;


//___________________________________________________________________________
class AliESDTagCreator : public AliTagCreator {

 public:
  AliESDTagCreator();
  ~AliESDTagCreator(); 

  void CreateESDTags(Int_t fFirstEvent, Int_t fLastEvent);

  Bool_t ReadGridCollection(TGridResult *result);
  Bool_t ReadLocalCollection(const char *localpath);
  Bool_t ReadCAFCollection(const char *filename);
  
 protected:  
  void CreateTag(TChain *chain, const char* type);
  void CreateTag(TFile* file, const char *guid, const char *md5, const char *turl, Long64_t size, Int_t Counter);
  void CreateTag(TFile* file, const char *filepath, Int_t Counter);

 private:
  TChain *fChain; //chain of esd files
  TList  *fGUIDList; //TList of guid TObjString
  TList  *fMD5List; //TList of md5 TObjString
  TList  *fTURLList; //TList of turl TObjString
  MemInfo_t *meminfo; //mem info

  ClassDef(AliESDTagCreator,0)  
};

#endif

