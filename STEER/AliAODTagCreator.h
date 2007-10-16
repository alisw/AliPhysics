#ifndef ALIAODTAGCREATOR_H
#define ALIAODTAGCREATOR_H
/*  See cxx source for full Copyright notice */


/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliAODTagCreator
//   This is the AliAODTagCreator class for the tag creation (post process)
//
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                        AliAODTagCreator                              //
//                                                                      //
//           Implementation of the tag creation mechanism.              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


//ROOT
//#include <TObject.h>

#include <AliTagCreator.h>

class TChain;
class TFile;
class TGridResult;


//___________________________________________________________________________
class AliAODTagCreator : public AliTagCreator {

 public:
  AliAODTagCreator();
  ~AliAODTagCreator(); 

  void CreateAODTags(Int_t fFirstEvent, Int_t fLastEvent);

  Bool_t ReadGridCollection(TGridResult *result);
  Bool_t ReadLocalCollection(const char *localpath);
  Bool_t ReadCAFCollection(const char *filename);
  
 protected:  
  void CreateTag(TChain *chain, const char* type);
  void CreateTag(TFile* file, const char *guid, const char *md5, const char *turl, Long64_t size, Int_t Counter);
  void CreateTag(TFile* file, const char *filepath, Int_t Counter);
  
 private:
  TChain *fChain; //chain of esd files

  ClassDef(AliAODTagCreator,0)  
};

#endif

