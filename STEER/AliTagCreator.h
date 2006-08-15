#ifndef ALITAGCREATOR_H
#define ALITAGCREATOR_H
/*  See cxx source for full Copyright notice */


/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliTagCreator
//   This is the AliTagCreator class for the tag creation (post process)
//
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                        AliTagCreator                                 //
//                                                                      //
//           Implementation of the tag creation mechanism.              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


//ROOT
#include <TObject.h>

class TFile;
class TGridResult;


//___________________________________________________________________________
class AliTagCreator : public TObject {

 public:
  AliTagCreator();
  ~AliTagCreator(); 

  //____________________________________________________//
  Bool_t MergeTags();
  Bool_t MergeTags(TGridResult *result);

  void SetSE(const char *se){fSE = se;}
  void SetStorage(Int_t storage);
  void SetGridPath(const char *gridpath){fgridpath = gridpath;}

  Bool_t ReadGridCollection(TGridResult *result);
  Bool_t ReadLocalCollection(const char *localpath);
  Bool_t ReadCAFCollection(const char *filename);

  //____________________________________________________//
 protected:
  TString fSE;   //the defined storage element
  TString fgridpath;   //the alien location of the tag files
  Int_t fStorage;  //0:local - 1:grid
  
  void CreateTag(TFile* file, const char *guid, const char *md5, const char *turl, Long64_t size, Int_t Counter);
  void CreateTag(TFile* file, const char *filepath, Int_t Counter);
 
  ClassDef(AliTagCreator,0)  
};

#endif

