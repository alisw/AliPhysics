#ifndef ALIDIELECTRONHISTOS_H
#define ALIDIELECTRONHISTOS_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////////////////
//                                                                                       //
// Generic Histogram container with support for groups and filling of groups by passing  //
// a vector of data                                                                      //
//                                                                                       //
// Authors:                                                                              //
//   Jens Wiechula <Jens.Wiechula@cern.ch>                                               //
//   Julian Book   <Julian.Book@cern.ch>                                                 //
//                                                                                       //
///////////////////////////////////////////////////////////////////////////////////////////

#include <Rtypes.h>

#include <TNamed.h>
// #include <TCollection.h>
#include <THashList.h>
#include <TVectorDfwd.h>
#include <THnBase.h>
#include <TBits.h>

class TH1;
class TString;
class TList;
// class TVectorT<double>;

class AliDielectronHistos : public TNamed {
public:

  AliDielectronHistos();
  AliDielectronHistos(const char* name, const char* title);
  virtual ~AliDielectronHistos();

  enum {kNoAutoFill=1000000000, kNoProfile=999, kNoWeights=998};

  void UserProfile(const char* histClass,const char *name, const char* title,
		   UInt_t valTypeP,
		   Int_t nbinsX, Double_t xmin, Double_t xmax,
		   UInt_t valTypeX, Bool_t logBinX=kFALSE, TString option="",
		   UInt_t valTypeW=kNoWeights);

  void UserProfile(const char* histClass,const char *name, const char* title,
		   UInt_t valTypeP,
		   Int_t nbinsX, Double_t xmin, Double_t xmax,
		   Int_t nbinsY, Double_t ymin, Double_t ymax,
		   UInt_t valTypeX, UInt_t valTypeY,
		   Bool_t logBinX=kFALSE, Bool_t logBinY=kFALSE, TString option="",
		   UInt_t valTypeW=kNoWeights);
  
  void UserProfile(const char* histClass,const char *name, const char* title,
		   UInt_t valTypeP,
		   Int_t nbinsX, Double_t xmin, Double_t xmax,
		   Int_t nbinsY, Double_t ymin, Double_t ymax,
		   Int_t nbinsZ, Double_t zmin, Double_t zmax,
		   UInt_t valTypeX, UInt_t valTypeY, UInt_t valTypeZ,
		   Bool_t logBinX=kFALSE, Bool_t logBinY=kFALSE, Bool_t logBinZ=kFALSE, TString option="",
		   UInt_t valTypeW=kNoWeights);
  
  void UserProfile(const char* histClass,const char *name, const char* title,
		   UInt_t valTypeP,
		   const char* binning, UInt_t valTypeX, TString option="",
		   UInt_t valTypeW=kNoWeights);
  
  void UserProfile(const char* histClass,const char *name, const char* title,
		   UInt_t valTypeP,
		   const TVectorD * const binsX, UInt_t valTypeX, TString option="",
		   UInt_t valTypeW=kNoWeights);
  
  void UserProfile(const char* histClass,const char *name, const char* title,
		   UInt_t valTypeP,
		   const TVectorD * const binsX, const TVectorD * const binsY,
		   UInt_t valTypeX, UInt_t valTypeY, TString option="",
		   UInt_t valTypeW=kNoWeights);
  
  void UserProfile(const char* histClass,const char *name, const char* title,
		   UInt_t valTypeP,
		   const TVectorD * const binsX, const TVectorD * const binsY, const TVectorD * const binsZ,
		   UInt_t valTypeX, UInt_t valTypeY, UInt_t valTypeZ, TString option="",
		   UInt_t valTypeW=kNoWeights);

  void UserHistogram(const char* histClass, Int_t ndim, Int_t *bins, Double_t *mins, Double_t *maxs, UInt_t *vars, UInt_t valTypeW=kNoWeights);
  void UserSparse(   const char* histClass, Int_t ndim, Int_t *bins, Double_t *mins, Double_t *maxs, UInt_t *vars, UInt_t valTypeW=kNoWeights);
  void UserHistogram(const char* histClass, Int_t ndim, TObjArray *limits, UInt_t *vars, UInt_t valTypeW=kNoWeights);
  void UserSparse(   const char* histClass, Int_t ndim, TObjArray *limits, UInt_t *vars, UInt_t valTypeW=kNoWeights);

  void UserSparse(   const char* histClass, const char *name, const char* title, Int_t ndim, Int_t *bins, Double_t *mins, Double_t *maxs, UInt_t *vars, UInt_t valTypeW=kNoWeights);
  void UserSparse(   const char* histClass, const char* name, const char* title, Int_t ndim, TObjArray *limits, UInt_t *vars, UInt_t valTypeW=kNoWeights);

  void UserHistogram(const char* histClass,const char *name, const char* title,
                     Int_t nbinsX, Double_t xmin, Double_t xmax, UInt_t valTypeX, Bool_t logBinX=kFALSE,
		     UInt_t valTypeW=kNoWeights)
  { UserProfile(histClass,name,title,kNoProfile,nbinsX,xmin,xmax,valTypeX,logBinX,"",valTypeW); }

  void UserHistogram(const char* histClass,const char *name, const char* title,
                     Int_t nbinsX, Double_t xmin, Double_t xmax, Int_t nbinsY, Double_t ymin, Double_t ymax,
                     UInt_t valTypeX, UInt_t valTypeY, Bool_t logBinX=kFALSE, Bool_t logBinY=kFALSE,
		     UInt_t valTypeW=kNoWeights)
  { UserProfile(histClass,name,title,kNoProfile,nbinsX,xmin,xmax,nbinsY,ymin,ymax,valTypeX,valTypeY,logBinX,logBinY,"",valTypeW); }

  void UserHistogram(const char* histClass,const char *name, const char* title,
                     Int_t nbinsX, Double_t xmin, Double_t xmax, Int_t nbinsY, Double_t ymin, Double_t ymax,
                     Int_t nbinsZ, Double_t zmin, Double_t zmax, UInt_t valTypeX, UInt_t valTypeY, UInt_t valTypeZ,
                     Bool_t logBinX=kFALSE, Bool_t logBinY=kFALSE, Bool_t logBinZ=kFALSE,
		     UInt_t valTypeW=kNoWeights)
  { UserProfile(histClass,name,title,kNoProfile,nbinsX,xmin,xmax,nbinsY,ymin,ymax,nbinsZ,zmin,zmax,valTypeX,valTypeY,valTypeZ,logBinX,logBinY,logBinZ,"",valTypeW); }

  void UserHistogram(const char* histClass,const char *name, const char* title,
                     const char* binning, UInt_t valTypeX, UInt_t valTypeW=kNoWeights)
  { UserProfile(histClass,name,title,kNoProfile,binning,valTypeX,"",valTypeW); }
  //  void UserHistogram(const char* histClass,const char *name, const char* title,
  //                     const TVectorD * const binsX, UInt_t valTypeX)
  //  { UserProfile(histClass,name,title,kNoProfile,binsX,valTypeX); }

  void UserHistogram(const char* histClass,const char *name, const char* title,
		     const TVectorD * const binsX, UInt_t valTypeX, UInt_t valTypeW=kNoWeights)
  { UserProfile(histClass,name,title,kNoProfile,binsX,valTypeX,"",valTypeW); }

  void UserHistogram(const char* histClass,const char *name, const char* title,
                     const TVectorD * const binsX, const TVectorD * const binsY, UInt_t valTypeX, UInt_t valTypeY,
		     UInt_t valTypeW=kNoWeights)
  { UserProfile(histClass,name,title,kNoProfile,binsX,binsY,valTypeX,valTypeY,"",valTypeW); }

  void UserHistogram(const char* histClass,const char *name, const char* title,
                     const TVectorD * const binsX, const TVectorD * const binsY, const TVectorD * const binsZ,
                     UInt_t valTypeX, UInt_t valTypeY, UInt_t valTypeZ,
		     UInt_t valTypeW=kNoWeights)
  { UserProfile(histClass,name,title,kNoProfile,binsX,binsY,binsZ,valTypeX,valTypeY,valTypeZ,"",valTypeW); }

  void UserHistogram(const char* histClass, TObject* hist, UInt_t valTypes=kNoAutoFill);

  // some functions needed to garantuee backward compatibility
  void UserHistogram(const char* histClass,const char *name, const char* title,
                     Int_t nbinsX, Double_t xmin, Double_t xmax)
  { UserProfile(histClass,name,title,kNoProfile,nbinsX,xmin,xmax,kNoAutoFill); }
  void UserHistogram(const char* histClass,const char *name, const char* title,
                     Int_t nbinsX, Double_t xmin, Double_t xmax,
                     Int_t nbinsY, Double_t ymin, Double_t ymax)
  { UserProfile(histClass,name,title,kNoProfile,nbinsX,xmin,xmax,nbinsY,ymin,ymax,kNoAutoFill,kNoAutoFill); }


  void Fill(const char* histClass, const char* name, Double_t xval);
  void Fill(const char* histClass, const char* name, Double_t xval, Double_t yval);
  void Fill(const char* histClass, const char* name, Double_t xval, Double_t yval, Double_t zval);
  
//   void FillClass(const char* histClass, const TVectorD &vals);
  void FillClass(const char* histClass, Int_t nValues, const Double_t *values);
  
  TObject* GetHist(const char* histClass, const char* name) const;
  TH1* GetHistogram(const char* histClass, const char* name) const;
  TObject* GetHist(const char* cutClass, const char* histClass, const char* name) const;
  TH1* GetHistogram(const char* cutClass, const char* histClass, const char* name) const;

  void SetHistogramList(THashList &list, Bool_t setOwner=kTRUE);
  void ResetHistogramList(){fHistoList.Clear();}
  const THashList* GetHistogramList() const {return &fHistoList;}

  void SetList(TList * const list) { fList=list; }
  TList *GetList() const { return fList; }
	TBits *GetUsedVars() const { return fUsedVars; }

  void AddClass(const char* histClass);

  void DumpToFile(const char* file="histos.root");
  void ReadFromFile(const char* file="histos.root");
  
  virtual void Print(const Option_t* option = "") const;
  virtual void Draw(const Option_t* option = "");
  virtual void DrawSame(const char* histName, const Option_t *opt="leg can");

  void SetReservedWords(const char* words);

//   virtual void       Add(TObject *obj) {};
//   virtual void       Clear(Option_t *option="") {};
//   virtual void       Delete(Option_t *option="") {};
//   virtual TObject  **GetObjectRef(const TObject *obj) const { return 0; }
//   virtual TIterator *MakeIterator(Bool_t dir = kIterForward) const ;
//   virtual TObject   *Remove(TObject *obj) { return 0; }

  Bool_t SetCutClass(const char* cutClass);
  static void StoreVariables(TObject *obj, UInt_t valType[20]);
  static void StoreVariables(TH1 *obj, UInt_t valType[20]);
  static void StoreVariables(THnBase *obj, UInt_t valType[20]);
  static void AdaptNameTitle(TH1 *hist, const char* histClass);
  static Int_t GetPrecision(Double_t value);
  static void FillValues(TObject *obj, const Double_t *values);
  static void FillValues(TH1 *obj, const Double_t *values);
  static void FillValues(THnBase *obj, const Double_t *values);

private:

  void FillVarArray(TObject *obj, UInt_t *valType);

  THashList fHistoList;             //-> list of histograms
  TList    *fList;                  //! List of list of histograms
	TBits     *fUsedVars;            // list of used variables

  TString *fReservedWords;          //! list of reserved words
  void UserHistogramReservedWords(const char* histClass, const TObject *hist, UInt_t valTypes);
  void FillClass(THashTable *classTable, Int_t nValues, Double_t *values);
  
  void PrintPDF(Option_t* opt);
  void PrintStructure() const;

  Bool_t IsHistogramOk(const char* classTable, const char* name);
  

  AliDielectronHistos(const AliDielectronHistos &hist);
  AliDielectronHistos& operator = (const AliDielectronHistos &hist);

  ClassDef(AliDielectronHistos,4)
};

#endif

