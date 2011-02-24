#ifndef ALIDIELECTRONHISTOS_H
#define ALIDIELECTRONHISTOS_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

///////////////////////////////////////////////////////////////////////////////////////////
//                                                                                       //
// Generic Histogram container with support for groups and filling of groups by passing  //
// a vector of data                                                                      //
//                                                                                       //
// Authors:                                                                              //
//   Jens Wiechula <Jens.Wiechula@cern.ch>                                               //
//                                                                                       //
///////////////////////////////////////////////////////////////////////////////////////////


#include <TNamed.h>
// #include <TCollection.h>
#include <THashList.h>
#include <TVectorDfwd.h>

class TH1;
class TString;
class TList;
// class TVectorT<double>;

class AliDielectronHistos : public TNamed {
public:
  AliDielectronHistos();
  AliDielectronHistos(const char* name, const char* title);
  virtual ~AliDielectronHistos();

  
  void UserHistogram(const char* histClass,const char *name, const char* title,
                     Int_t nbinsX, Double_t xmin, Double_t xmax,
                     UInt_t valTypeX=kNoAutoFill, Bool_t logBinX=kFALSE);
  void UserHistogram(const char* histClass,const char *name, const char* title,
                     Int_t nbinsX, Double_t xmin, Double_t xmax,
                     Int_t nbinsY, Double_t ymin, Double_t ymax,
                     UInt_t valTypeX=kNoAutoFill, UInt_t valTypeY=0,
                     Bool_t logBinX=kFALSE, Bool_t logBinY=kFALSE);
  void UserHistogram(const char* histClass,const char *name, const char* title,
                     Int_t nbinsX, Double_t xmin, Double_t xmax,
                     Int_t nbinsY, Double_t ymin, Double_t ymax,
                     Int_t nbinsZ, Double_t zmin, Double_t zmax,
                     UInt_t valTypeX=kNoAutoFill, UInt_t valTypeY=0, UInt_t valTypeZ=0,
                     Bool_t logBinX=kFALSE, Bool_t logBinY=kFALSE, Bool_t logBinZ=kFALSE);
  
  void UserHistogram(const char* histClass,const char *name, const char* title,
                     const char* binning,
                     UInt_t valTypeX=kNoAutoFill);

  void UserHistogram(const char* histClass,const char *name, const char* title,
                     const TVectorD * const binsX,
                     UInt_t valTypeX=kNoAutoFill);
  void UserHistogram(const char* histClass,const char *name, const char* title,
                     const TVectorD * const binsX, const TVectorD * const binsY,
                     UInt_t valTypeX=kNoAutoFill, UInt_t valTypeY=0);
  void UserHistogram(const char* histClass,const char *name, const char* title,
                     const TVectorD * const binsX, const TVectorD * const binsY, const TVectorD * const binsZ,
                     UInt_t valTypeX=kNoAutoFill, UInt_t valTypeY=0, UInt_t valTypeZ=0);
  
  void UserHistogram(const char* histClass, TH1* hist, UInt_t valTypes=kNoAutoFill);


  void Fill(const char* histClass, const char* name, Double_t xval);
  void Fill(const char* histClass, const char* name, Double_t xval, Double_t yval);
  void Fill(const char* histClass, const char* name, Double_t xval, Double_t yval, Double_t zval);
  
//   void FillClass(const char* histClass, const TVectorD &vals);
  void FillClass(const char* histClass, Int_t nValues, const Double_t *values);

  
  TH1* GetHistogram(const char* histClass, const char* name) const;
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
  void UserHistogramReservedWords(const char* histClass, const TH1 *hist, UInt_t valTypes);
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

