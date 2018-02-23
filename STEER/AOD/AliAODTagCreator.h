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
class TTree;
class TFile;
class TGridResult;

class AliAODEvent;
class AliRunTag;
class AliEventTag;
class AliFileTag;

//___________________________________________________________________________
class AliAODTagCreator : public AliTagCreator {

 public:
  AliAODTagCreator();
  ~AliAODTagCreator(); 

  void CreateAODTags(Int_t fFirstEvent, Int_t fLastEvent, TList *grpList);

  Bool_t ReadGridCollection(TGridResult *result);
  Bool_t ReadLocalCollection(const char *localpath, const char* pattern = "AliAOD.root");
  Bool_t ReadCAFCollection(const char *filename);
  void FillEventTag(AliAODEvent* aod, AliEventTag* evtTag);  
 protected:  
  void CreateTag(TChain *chain, const char* type);
  void CreateTags(const char* type = "");

  
 private:
  AliAODTagCreator(const AliAODTagCreator& creator);             
  AliAODTagCreator& operator=(const AliAODTagCreator& creator);    
 private:
  TChain       *fChain;     //! Chain of esd files
  AliAODEvent  *fAODEvent;  //! AOD Event 
  TTree        *fTreeT;     //! Tag Tree
  AliRunTag    *fRunTag;    //! Run tag
  TChain       *fTreeTEsd;  //! ESD tag Tree
  AliRunTag    *fRunTagEsd; //! ESD run tag
  
  ClassDef(AliAODTagCreator,0)  
};

#endif

