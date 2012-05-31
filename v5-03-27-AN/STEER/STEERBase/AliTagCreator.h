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

class TGridResult;
class TChain;

//___________________________________________________________________________
class AliTagCreator : public TObject {

 public:
  AliTagCreator();
  ~AliTagCreator(); 

  //____________________________________________________//
  Bool_t MergeTags(const char* type);
  Bool_t MergeTags(const char* type, TGridResult *result);
  Bool_t MergeTags(const char* type, const char *inflist);

  Bool_t MergeTagsForRun(const char* type);
  Bool_t MergeTagsForRun(const char* type, TGridResult *result);
  Bool_t MergeTagsForRun(const char* type, const char *inflist);

  void SetSE(const char *se){fSE = se;}
  void SetStorage(Int_t storage);
  void SetGridPath(const char *gridpath){fgridpath = gridpath;}

  //____________________________________________________//
 protected:

  Bool_t MergeToSingleRunTag(TChain *chain, const char *filename);

  TString fSE;   //the defined storage element
  TString fgridpath;   //the alien location of the tag files
  Int_t fStorage;  //0:local - 1:grid
   
  ClassDef(AliTagCreator,0)  
};

#endif

