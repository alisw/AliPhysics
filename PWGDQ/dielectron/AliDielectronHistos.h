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
//                                                                                       //
///////////////////////////////////////////////////////////////////////////////////////////

#include <Rtypes.h>

#include <TNamed.h>
// #include <TCollection.h>
#include <THashList.h>
#include <TVectorDfwd.h>
#include <THn.h>
#include <THnSparse.h>

//class THn;
class TH1;
class TString;
class TList;
// class TVectorT<double>;

class AliDielectronHistos : public TNamed {
public:

  AliDielectronHistos();
  AliDielectronHistos(const char* name, const char* title);
  virtual ~AliDielectronHistos();

  
  void UserProfile(const char* histClass,const char *name, const char* title,
		   UInt_t valTypeP,
		   Int_t nbinsX, Double_t xmin, Double_t xmax,
		   UInt_t valTypeX=kNoAutoFill, Bool_t logBinX=kFALSE, TString option="");

  void UserProfile(const char* histClass,const char *name, const char* title,
		   UInt_t valTypeP,
		   Int_t nbinsX, Double_t xmin, Double_t xmax,
		   Int_t nbinsY, Double_t ymin, Double_t ymax,
		   UInt_t valTypeX=kNoAutoFill, UInt_t valTypeY=0,
		   Bool_t logBinX=kFALSE, Bool_t logBinY=kFALSE, TString option="");
  
  void UserProfile(const char* histClass,const char *name, const char* title,
		   UInt_t valTypeP,
		   Int_t nbinsX, Double_t xmin, Double_t xmax,
		   Int_t nbinsY, Double_t ymin, Double_t ymax,
		   Int_t nbinsZ, Double_t zmin, Double_t zmax,
		   UInt_t valTypeX=kNoAutoFill, UInt_t valTypeY=0, UInt_t valTypeZ=0,
		   Bool_t logBinX=kFALSE, Bool_t logBinY=kFALSE, Bool_t logBinZ=kFALSE, TString option="");
  
  void UserProfile(const char* histClass,const char *name, const char* title,
		   UInt_t valTypeP,
		   const char* binning, UInt_t valTypeX=kNoAutoFill, TString option="");
  
  void UserProfile(const char* histClass,const char *name, const char* title,
		   UInt_t valTypeP,
		   const TVectorD * const binsX, UInt_t valTypeX=kNoAutoFill, TString option="");
  
  void UserProfile(const char* histClass,const char *name, const char* title,
		   UInt_t valTypeP,
		   const TVectorD * const binsX, const TVectorD * const binsY,
		   UInt_t valTypeX=kNoAutoFill, UInt_t valTypeY=0, TString option="");
  
  void UserProfile(const char* histClass,const char *name, const char* title,
		   UInt_t valTypeP,
		   const TVectorD * const binsX, const TVectorD * const binsY, const TVectorD * const binsZ,
		   UInt_t valTypeX=kNoAutoFill, UInt_t valTypeY=0, UInt_t valTypeZ=0, TString option="");

  void UserHistogram(const char* histClass, Int_t ndim, Int_t *bins, Double_t *mins, Double_t *maxs, UInt_t *vars);
  void UserSparse(const char* histClass, Int_t ndim, Int_t *bins, Double_t *mins, Double_t *maxs, UInt_t *vars);

  void UserHistogram(const char* histClass,const char *name, const char* title,
                     Int_t nbinsX, Double_t xmin, Double_t xmax, UInt_t valTypeX=kNoAutoFill, Bool_t logBinX=kFALSE)
  { UserProfile(histClass,name,title,999,nbinsX,xmin,xmax,valTypeX,logBinX); }

  void UserHistogram(const char* histClass,const char *name, const char* title,
                     Int_t nbinsX, Double_t xmin, Double_t xmax, Int_t nbinsY, Double_t ymin, Double_t ymax,
                     UInt_t valTypeX=kNoAutoFill, UInt_t valTypeY=0, Bool_t logBinX=kFALSE, Bool_t logBinY=kFALSE)
  { UserProfile(histClass,name,title,999,nbinsX,xmin,xmax,nbinsY,ymin,ymax,valTypeX,valTypeY,logBinX,logBinY); }

  void UserHistogram(const char* histClass,const char *name, const char* title,
                     Int_t nbinsX, Double_t xmin, Double_t xmax, Int_t nbinsY, Double_t ymin, Double_t ymax,
                     Int_t nbinsZ, Double_t zmin, Double_t zmax, UInt_t valTypeX=kNoAutoFill, UInt_t valTypeY=0, UInt_t valTypeZ=0,
                     Bool_t logBinX=kFALSE, Bool_t logBinY=kFALSE, Bool_t logBinZ=kFALSE)
  { UserProfile(histClass,name,title,999,nbinsX,xmin,xmax,nbinsY,ymin,ymax,nbinsZ,zmin,zmax,valTypeX,valTypeY,valTypeZ,logBinX,logBinY,logBinZ); }

  void UserHistogram(const char* histClass,const char *name, const char* title,
                     const char* binning, UInt_t valTypeX=kNoAutoFill)
  { UserProfile(histClass,name,title,999,binning,valTypeX); }

  void UserHistogram(const char* histClass,const char *name, const char* title,
                     const TVectorD * const binsX, UInt_t valTypeX=kNoAutoFill)
  { UserProfile(histClass,name,title,999,binsX,valTypeX); }

  void UserHistogram(const char* histClass,const char *name, const char* title,
                     const TVectorD * const binsX, const TVectorD * const binsY, UInt_t valTypeX=kNoAutoFill, UInt_t valTypeY=0)
  { UserProfile(histClass,name,title,999,binsX,binsY,valTypeX,valTypeY); }

  void UserHistogram(const char* histClass,const char *name, const char* title,
                     const TVectorD * const binsX, const TVectorD * const binsY, const TVectorD * const binsZ,
                     UInt_t valTypeX=kNoAutoFill, UInt_t valTypeY=0, UInt_t valTypeZ=0)
  { UserProfile(histClass,name,title,999,binsX,binsY,binsZ,valTypeX,valTypeY,valTypeZ); }
  
  void UserHistogram(const char* histClass, TObject* hist, UInt_t valTypes=kNoAutoFill);

  void Fill(const char* histClass, const char* name, Double_t xval);
  void Fill(const char* histClass, const char* name, Double_t xval, Double_t yval);
  void Fill(const char* histClass, const char* name, Double_t xval, Double_t yval, Double_t zval);
  
//   void FillClass(const char* histClass, const TVectorD &vals);
  void FillClass(const char* histClass, Int_t nValues, const Double_t *values);
  
  void StoreVariables(TObject *obj, UInt_t valType[20]);
  void StoreVariables(TH1 *obj, UInt_t valType[20]);
  void StoreVariables(THn *obj, UInt_t valType[20]);
  void StoreVariables(THnSparse *obj, UInt_t valType[20]);
  void AdaptNameTitle(TH1 *hist, const char* histClass);
  void FillVarArray(TObject *obj, UInt_t *valType);
  void FillValues(TObject *obj, const Double_t *values);
  void FillValues(TH1 *obj, const Double_t *values);
  void FillValues(THn *obj, const Double_t *values);
  void FillValues(THnSparse *obj, const Double_t *values);

  TObject* GetHist(const char* histClass, const char* name) const;
  TH1* GetHistogram(const char* histClass, const char* name) const;
  TObject* GetHist(const char* cutClass, const char* histClass, const char* name) const;
  TH1* GetHistogram(const char* cutClass, const char* histClass, const char* name) const;
  
  void SetHistogramList(THashList &list, Bool_t setOwner=kTRUE);
  void ResetHistogramList(){fHistoList.Clear();}
  const THashList* GetHistogramList() const {return &fHistoList;}

  void SetList(TList * const list) { fList=list; }
  TList *GetList() const { return fList; }

  void AddClass(const char* histClass);

  void DumpToFile(const char* file="histos.root");
  void ReadFromFile(const char* file="histos.root");
  
  virtual void Print(const Option_t* option = "") const;
  virtual void Draw(const Option_t* option = "");
  virtual void DrawSame(const char* histName, const Option_t *opt="leg can");

  void SetReservedWords(const char* words);
  //  void StoreVarForProfile(TObject *obj, UInt_t valType);
//   virtual void       Add(TObject *obj) {};
//   virtual void       Clear(Option_t *option="") {};
//   virtual void       Delete(Option_t *option="") {};
//   virtual TObject  **GetObjectRef(const TObject *obj) const { return 0; }
//   virtual TIterator *MakeIterator(Bool_t dir = kIterForward) const ;
//   virtual TObject   *Remove(TObject *obj) { return 0; }

  Bool_t SetCutClass(const char* cutClass);
  
private:
  THashList fHistoList;             //-> list of histograms
  TList    *fList;                  //! List of list of histograms

  TString *fReservedWords;          //! list of reserved words
  void UserHistogramReservedWords(const char* histClass, const TObject *hist, UInt_t valTypes);
  void FillClass(THashTable *classTable, Int_t nValues, Double_t *values);
  
  void PrintPDF(Option_t* opt);
  void PrintStructure() const;

  Bool_t IsHistogramOk(const char* classTable, const char* name);
  
  enum {kNoAutoFill=1000000000};

  AliDielectronHistos(const AliDielectronHistos &hist);
  AliDielectronHistos& operator = (const AliDielectronHistos &hist);

  ClassDef(AliDielectronHistos,2)
};

#endif

