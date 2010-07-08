// This class is used to print a table which can be pasted on a latex
// document

#ifndef ALILATEXTABLE_H
#define ALILATEXTABLE_H

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TObject.h"
#include "TString.h"
class TObjArray;

#endif


using namespace std;

class AliLatexTable : public TObject {

public:
  AliLatexTable() : fNcol(1), fFormat("c"), fRows(0), fCols(0), fNcolReady(0) {;}
  AliLatexTable(Int_t ncol, TString format);
  ~AliLatexTable();

  // first you set the value of each column
  void SetNextCol(Int_t val); 
  void SetNextCol(Int_t val, Int_t err); 

  void SetNextCol(Double_t val, Int_t scientificNotation = -1); // if different from -1 gives significant digits
  void SetNextCol(Double_t val, Double_t err, Int_t scientificNotation = -1); 
  void SetNextCol(Double_t val, Double_t err, Double_t errSyst, Int_t scientificNotation = -1); 

  void SetNextCol(TString val); 
//   // allows to use printf syntax
//   void SetNextColPrintf(const char *va_(fmt), ...);

  // Then you add the row (it's up to user make sure all columns are
  // there)
  void InsertRow();


  // insert a row without building it up with the methods. May be
  // usefull for header or multcol rows
  void InsertCustomRow(TString row);


  void InsertHline();

  void PrintTable(Option_t * opt = "");
  //  void ParseExponent(TString &expo);
  void GetMantissaAndExpBase10(Double_t num, Double_t &man, Double_t &exp) ;

  // used to print columns with correct width
  Int_t * GetColWidths();
  void StripLatex(TString &row, TString format) ;
  

private:

  Int_t fNcol;     // number of columns
  TString fFormat; // latex format (es "c|ccc")
  
  TObjArray * fRows; // rows 
  TObjArray * fCols; // columns
  
  Int_t fNcolReady; // number of cols ready to be insert


  AliLatexTable(const AliLatexTable&);
  AliLatexTable& operator=(const AliLatexTable&);

  ClassDef(AliLatexTable, 1)


};

#endif

