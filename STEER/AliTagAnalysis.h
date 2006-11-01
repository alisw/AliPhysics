#ifndef ALITAGANALYSIS_H
#define ALITAGANALYSIS_H
/*  See cxx source for full Copyright notice */


/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliTagAnalysis
//   This is the AliTagAnalysis class for the tag analysis
//
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                        AliTagAnalysis                                //
//                                                                      //
//           Implementation of the tag analysis mechanism.              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


//ROOT
#include <TObject.h>

class AliEventTag;
class TChain;
class AliEventTagCuts;
class AliRunTagCuts;
class TGridResult;
class TTreeFormula;

class AliTagAnalysis : public TObject {
 public:
  AliTagAnalysis();
  ~AliTagAnalysis(); 

  void ChainLocalTags(const char *dirname);
  void ChainGridTags(TGridResult *result);
  TChain *QueryTags(const char *fRunCut, const char *fEventCut);
  TChain *QueryTags(AliRunTagCuts *RunTagCuts, AliEventTagCuts *EvTagCuts);
  Bool_t CreateXMLCollection(const char* name, AliRunTagCuts *RunTagCuts, AliEventTagCuts *EvTagCuts);
  Bool_t CreateXMLCollection(const char* name, const char *fRunCut, const char *fEventCut);
   
 protected:
  TGridResult *ftagresult; //the results from the tag grid query     
  TString fTagDirName; //the location of the locally stored tags
    
  static TChain *fgChain; //tag chain 
  TChain *fChain; //tag chain 

 private:
  AliTagAnalysis(const AliTagAnalysis & source);
  AliTagAnalysis & operator=(const AliTagAnalysis & source);
       
  ClassDef(AliTagAnalysis,0)  
};

#endif

