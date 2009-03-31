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
#include <TString.h>
#include <TEntryList.h>

class AliEventTag;
class TChain;
class AliEventTagCuts;
class AliDetectorTagCuts;
class AliLHCTagCuts;
class AliRunTagCuts;
class TGridResult;
class TTreeFormula;

//____________________________________________________//
class AliTagAnalysis : public TObject {
 public:
  AliTagAnalysis();
  AliTagAnalysis(const char* type);
  ~AliTagAnalysis(); 
  
  void SetType(const char* type) {fAnalysisType = type;}
  const char* GetType() {return fAnalysisType.Data();}
  Bool_t AddTagsFile(const char *alienUrl);
  void ChainLocalTags(const char *dirname);
  TChain *ChainGridTags(TGridResult *result);
  
  TChain *QueryTags(AliRunTagCuts *runTagCuts, 
		    AliLHCTagCuts *lhcTagCuts, 
		    AliDetectorTagCuts *detTagCuts, 
		    AliEventTagCuts *evTagCuts);
  TChain *QueryTags(const char *fRunCut, 
		    const char *fLHCCut, 
		    const char *fDetectorCut, 
		    const char *fEventCut);  

  Bool_t CreateXMLCollection(const char* name, 
			     AliRunTagCuts *runTagCuts, 
			     AliLHCTagCuts *lhcTagCuts, 
			     AliDetectorTagCuts *detTagCuts, 
			     AliEventTagCuts *evTagCuts);
  Bool_t CreateXMLCollection(const char* name, 
			     const char *fRunCut, 
			     const char *fLHCCut, 
			     const char *fDetectorCut, 
			     const char *fEventCut);

  Bool_t CreateAsciiCollection(const char* name, 
			       AliRunTagCuts *runTagCuts, 
			       AliLHCTagCuts *lhcTagCuts, 
			       AliDetectorTagCuts *detTagCuts, 
			       AliEventTagCuts *evTagCuts);
  Bool_t CreateAsciiCollection(const char* name, 
			       const char *fRunCut, 
			       const char *fLHCCut, 
			       const char *fDetectorCut, 
			       const char *fEventCut);

  TChain *GetInputChain(const char* system, const char *wn);
  TChain *GetChainFromCollection(const char* collectionname, 
				 const char* treename);
  
  TEntryList *GetGlobalList() {return fGlobalList;}
  //____________________________________________________//
 protected:
  TGridResult *ftagresult; //the results from the tag grid query     
  TString fTagDirName; //the location of the locally stored tags
  
  TChain *fChain; //tag chain 
  
  TString fAnalysisType; //define the type of analysis (esd or aod)

  TEntryList *fGlobalList; //global TEntryList
  
  //____________________________________________________//
 private:
  AliTagAnalysis(const AliTagAnalysis & source);
  AliTagAnalysis & operator=(const AliTagAnalysis & source);
       
  ClassDef(AliTagAnalysis,0)  
};

#endif

