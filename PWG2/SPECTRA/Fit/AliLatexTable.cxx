// ------------------------------------------------------------------
//
//                           AliLatexTable
//
// This class is used to produce tables using latex sintax. It can
// than print the table as Latex source code (to be pasted on a paper)
// or ASCII.
//
// The basic idea is that you add columns one after the other with the
// SetNextCol method and than you insert this row. The SetNextCol
// method comes in different flavours to alow the insertion of
// different type of values (numbers with or without errors,
// strings...).
//
// TODO:
// 1. Make the class drawable
// 2. Implement vertical lines in ascii print
// 3. Print output in HTML format
//
// Author: Michele Floris, CERN
// ------------------------------------------------------------------

#include "AliLatexTable.h"
#include "TString.h"
#include "TRegexp.h"
#include "TObjString.h"
#include "TObjArray.h"
#include <stdarg.h>
#include "snprintf.h"
#include "Varargs.h"
#include "TPRegexp.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include "AliLog.h"

using namespace std;

ClassImp(AliLatexTable)

AliLatexTable::AliLatexTable() : fNcol(0), fFormat(""), fRows(0), fCols(0), fNcolReady(0){
  // default constructor
  fNcol = 1;
  fFormat = "c";
  fNcolReady = 0;
  fRows = new TObjArray;
  fCols = new TObjArray;

  fRows->SetOwner();
  fCols->SetOwner();

}

AliLatexTable::AliLatexTable(Int_t ncol, TString format) : fNcol(0), fFormat(""), fRows(0), fCols(0), fNcolReady(0){
  // constructor, specify number of cols 
  fNcol = ncol;
  fFormat = format;
  fNcolReady = 0;
  fRows = new TObjArray;
  fCols = new TObjArray;

  fRows->SetOwner();
  fCols->SetOwner();
}


AliLatexTable::~AliLatexTable() {

  // dtor
  if (fRows) delete fRows;
  if (fCols) delete fCols;

}

void AliLatexTable::SetNextCol(Int_t val){
  // Set next column in current row - integer
  char col[200];
  snprintf(col, 200, " %d ", val);
  SetNextCol(col);

} 

void AliLatexTable::SetNextCol(Int_t val, Int_t err){

  // Set next column in current row - int +- int
  char col[200];
  snprintf(col, 200, " $%d \\pm %d$ ", val, err);
  SetNextCol(col);

} 

void AliLatexTable::SetNextCol(Double_t val, Int_t scientificNotation, Bool_t rounding){

  // Set next column in current row - double, specify resolution

  char col[200];
  if(rounding) {
    if(scientificNotation >= 0) {
      char format[100];

      //    cout << format << endl;    
      Double_t mantissa, exp;
      GetMantissaAndExpBase10(val, mantissa, exp);
    
      if (exp == 0) {
	snprintf(format, 100," $%%%d.%df $", scientificNotation, scientificNotation);
	snprintf(col, 200, format, mantissa);      
      }
      else  {
	snprintf(format, 100," $%%%d.%df \\cdot 10^{%%0.0f} $", scientificNotation, scientificNotation);
	snprintf(col, 200,format, mantissa, exp);
      }


    

      SetNextCol(col);

    } else {
      char format[100];
      snprintf(format, 100," $%%%d.%df $", -scientificNotation,-scientificNotation);    
      snprintf(col, 200, format , TMath::Nint(val*TMath::Power(10,-scientificNotation))/TMath::Power(10,-scientificNotation));
      SetNextCol(col);
    } 
  }else {
    snprintf(col,200, " %f ", val);
    SetNextCol(col);
  }

} 
void AliLatexTable::SetNextCol(Double_t val, Double_t err, Int_t scientificNotation, Bool_t rounding){

  // Set next column in current row - double +- double, specify resolution

  // scientific notation is used to determine number of
  // digits in error 

  //if it is > 0 exp notation is used 


  char col[200];
  if(rounding) {
    if(scientificNotation >=0 ) {
      
      Double_t mantissa, exp;
      GetMantissaAndExpBase10(val, mantissa, exp);
      
      Double_t mantissa_err, exp_err;
      GetMantissaAndExpBase10(err, mantissa_err, exp_err);

      Int_t nSigDigits =TMath::Nint(exp - exp_err);
      if(scientificNotation != 0) nSigDigits =  nSigDigits + scientificNotation - 1;


      char format[100];
      if (exp == 0) {
	snprintf(format, 100," $%%%d.%df \\pm %%%d.%df $", nSigDigits, nSigDigits, nSigDigits, nSigDigits);
	snprintf(col, 200,format, mantissa, mantissa_err/TMath::Power(10,exp - exp_err));
      }
      else  {
	snprintf(format, 100, " $%%%d.%df \\pm %%%d.%df \\cdot 10^{%%0.0f}$", nSigDigits, nSigDigits, nSigDigits, nSigDigits);
	snprintf(col, 200, format, mantissa, mantissa_err/TMath::Power(10,exp - exp_err), exp);
      }


      //cout << format << endl;
    
      SetNextCol(col);

    } else  {
      char format[100];
      snprintf(format, 100, " $%%%d.%df \\pm %%%d.%df $", -scientificNotation,-scientificNotation,-scientificNotation,-scientificNotation);    
      snprintf(col, 200, format , TMath::Nint(val*TMath::Power(10,-scientificNotation))/TMath::Power(10,-scientificNotation), TMath::Nint(err*TMath::Power(10,-scientificNotation))/TMath::Power(10,-scientificNotation));
      SetNextCol(col);
    }
  }
  else {
    snprintf(col, 200, " $%f \\pm %f$ ", val, err);
    SetNextCol(col);
  }
} 

void AliLatexTable::SetNextCol(Double_t val, Double_t err, Double_t errSyst, Int_t scientificNotation, Bool_t rounding){

  // Set next column in current row - double +- double +- double, specify resolution

  // if scientific notation is != -1 is used to determine number of
  // digits in error 

  //if it is > 0 exp notation is used 

  //if it is < -1 number is truncated to the number of digits after
  //point dictated by scientificNotation,

  char col[200];
  if (rounding) {
    if(scientificNotation >=0 ) {

      Double_t mantissa, exp;
      GetMantissaAndExpBase10(val, mantissa, exp);

      Double_t mantissa_err, exp_err;
      GetMantissaAndExpBase10(err, mantissa_err, exp_err);

      Double_t mantissa_errsyst, exp_errsyst;
      GetMantissaAndExpBase10(errSyst, mantissa_errsyst, exp_errsyst);

      Int_t nSigDigits =TMath::Nint(exp - exp_err);
      if(scientificNotation != 0) nSigDigits =  nSigDigits + scientificNotation - 1;


      char format[100];
      if (exp == 0) {
	snprintf(format, 100, " $%%%d.%df \\pm %%%d.%df \\pm %%%d.%df $", 
		nSigDigits, nSigDigits, nSigDigits, nSigDigits, nSigDigits, nSigDigits);
	snprintf(col, 200, format, mantissa, 
		mantissa_err/TMath::Power(10,exp - exp_err), 
		mantissa_errsyst/TMath::Power(10,exp - exp_errsyst));
      }
      else  {
	snprintf(format, 100, " $%%%d.%df \\pm %%%d.%df  \\pm %%%d.%df \\cdot 10^{%%0.0f}$", 
		nSigDigits, nSigDigits, nSigDigits, nSigDigits, nSigDigits, nSigDigits);
	snprintf(col, 200, format, mantissa, 
		mantissa_err/TMath::Power(10,exp - exp_err), 
		mantissa_errsyst/TMath::Power(10,exp - exp_errsyst),
		exp);
      }


      //cout << format << endl;
    
      SetNextCol(col);

    } else  {
      char format[100];
      snprintf(format, 100, " $%%%d.%df \\pm %%%d.%df \\pm %%%d.%df $", -scientificNotation,-scientificNotation,-scientificNotation,-scientificNotation, -scientificNotation,-scientificNotation);    
      snprintf(col, 200, format ,TMath::Nint(val*TMath::Power(10,-scientificNotation))/TMath::Power(10,-scientificNotation), TMath::Nint(err*TMath::Power(10,-scientificNotation))/TMath::Power(10,-scientificNotation), TMath::Nint(errSyst*TMath::Power(10,-scientificNotation))/TMath::Power(10,-scientificNotation));
      SetNextCol(col);
    } 
  }
  else {
    snprintf(col, 200, " $%f \\pm %f  \\pm %f$ ", val, err, errSyst);
    SetNextCol(col);
  }
} 

// void AliLatexTable::SetNextColPrintf(const char *va_(fmt), ...){


//   const Int_t buf_size = 1024;  // hope it's long enough
//   char colBuffer[buf_size];

//   va_list ap;
//   va_start(ap, va_(fmt));
//   Int_t n = vsnprintf(colBuffer, buf_size, fmt, ap);

//   if (n == -1 || n >= buf_size) {
//     Error("SetNextCol", "vsnprintf error!");
//   }

//   va_end(ap);
     
// //   for(int iobj=0; iobj<num; iobj++){   //Loop until all objects are added
// //     TH1 * obj = (TH1 *) va_arg(arguments, void *);
// //     returnValue = returnValue && AddObject(obj,(char*)fOptMany.Data(), EColor(iobj+1));
// //   }
//   SetNextCol(colBuffer);
// }


void AliLatexTable::SetNextCol(TString val){

  // Set next column in current row - tstring

  fCols->Add(new TObjString(val));
  fNcolReady++;

} 

void AliLatexTable::InsertCustomRow(TString row){
  // insert a full row from string. Make sure you have all the columns
  fRows->Add(new TObjString(row));
}



void AliLatexTable::InsertRow(){

  // Insert the row, based on all the individual columns

  if ( fNcolReady != fNcol) {
    Warning("InsertRow", "Wrong number of cols: %d (!= %d)", fNcolReady, fNcol); 
  }

  TString row = "";
  for(Int_t icol = 0; icol < fNcol; icol++){
    row = row + ((TObjString*) fCols->At(icol))->String();
    if(icol != (fNcol-1)) row = row + " & ";
  }
  row+="\\\\";

  fRows->Add(new TObjString(row));

  fNcolReady = 0;
  fCols->Clear();

}

void AliLatexTable::InsertHline(){

  // insert an horizontal line
  fRows->Add(new TObjString("\\hline"));
  

}

Int_t * AliLatexTable::GetColWidths() {

  // Compute the width of columns, for nice ascii printing

  Int_t * col_widths= new Int_t[fNcol];
  Int_t nrow = fRows->GetEntriesFast();

  for(Int_t icol = 0; icol < fNcol; icol++){
    col_widths[icol] = 0;
  }
  

  for(Int_t irow = 0; irow < nrow; irow++){
    TString row(((TObjString*) fRows->At(irow))->String());
    if(row.Contains("\\hline"))continue;
    
    StripLatex(row, "ASCII");

    TObjArray * cols = row.Tokenize("&");
    if(cols->GetEntries() != fNcol) {
      Error("GetColWidths", "Wrong number of cols in row %s: %d - %d", row.Data(), cols->GetEntries(), fNcol);      
    }
    for(Int_t icol = 0; icol < fNcol; icol++){
      Int_t w = ((TObjString *) cols->At(icol))->String().Length();
      if (w>col_widths[icol]) col_widths[icol] = w;
    }
    
  }

  return col_widths;
}

void AliLatexTable::PrintTable(Option_t * opt){

  // Print the table on screen. You can specify the format with opt.
  // Currently supported:
  // "" -> LaTeX
  // "ASCII" -> plain ASCII
  // "HTML"  -> HTML, to be improved
  // "CSV"   -> skips hline, usefult for importing in excell 
  // "TWIKI" -> skips hline, usefult for importing in TWIKI

  if(TString(opt) == "ASCII" || TString(opt)=="HTML" ||  TString(opt)=="CSV" ||  TString(opt)=="TWIKI") {
    
    Int_t nrow = fRows->GetEntriesFast();

    Int_t * colWidths = GetColWidths();

    Int_t total_lenght = 0;
    for(Int_t icol = 0; icol < fNcol; icol++) total_lenght = total_lenght + colWidths[icol] + 2 ;
    

    for(Int_t irow = 0; irow < nrow; irow++){
      TString row = ((TObjString*) fRows->At(irow))->String();
      if (row.Contains("\\hline")){	
	if (TString(opt)!="CSV" && TString(opt)!="TWIKI") {
	  for(Int_t il = 0; il < total_lenght; il++) printf("-");
	  printf("\n");	  
	}
	continue;
      }
      StripLatex(row, opt);
      TObjArray * cols = row.Tokenize("&");
      if (TString(opt)=="TWIKI") printf(" | ");
      for(Int_t icol = 0; icol < fNcol; icol++){
	TString strTmp = ((TObjString *) cols->At(icol))->String();
	if(TString(opt)=="TWIKI" || TString(opt)=="HTML"){
	  strTmp.ReplaceAll("AMPER","&");
	}
	const char * colstr = strTmp.Data();
	char format [200];
	if (TString(opt)!="CSV") {
	  snprintf(format, 200, "%%%ds", colWidths[icol] + 2);	
	} else {
	  snprintf(format, 200, "%%s");	
	}
	printf(format, colstr);
	if (TString(opt)=="CSV") printf(", ");
	if (TString(opt)=="TWIKI") printf(" | ");

      }
      printf ("\n");
    }
    
    delete [] colWidths;
    return;
  }
  

  cout << "\\begin{tabular}{"<<fFormat<<"}"<<endl;

  Int_t nrow = fRows->GetEntriesFast();

  for(Int_t irow = 0; irow < nrow; irow++){
    cout << ((TObjString*) fRows->At(irow))->String() << endl;
  }
  
  cout << "\\end{tabular}" << endl;


}

// void AliLatexTable::ParseExponent(TString &expo){
// //    TString parseExponent = col;
//   TRegexp re = "e[+-][0-9][0-9]";
//   //  cout << col << endl;
  
//   if(expo.Contains(re)){
//     Int_t index = expo.Index(re);
//     TString replacement = "\\cdot 10^{";
//     replacement = replacement + expo(index+1, 3) +"}";
//     //    cout << replacement <<" --- "<< endl;
//     replacement.ReplaceAll("+","");
// //     cout << "B: " << expo << endl;      
// //     cout << "RE " << expo(re) << endl;
    
//     expo.ReplaceAll(expo(re), replacement);
//     //    cout << "A: " << expo << endl;      
//     // recursion to parse all exponents
//     if (expo.Contains(re)) ParseExponent(expo);
//   }
//   else Warning("", "Error parsing exponent");
// }

void AliLatexTable::GetMantissaAndExpBase10(Double_t num, Double_t &man, Double_t &exp) {

  // Helper used to get mantissa and exponent

  exp = TMath::Floor(TMath::Log10(TMath::Abs(num)));  
  man = num / TMath::Power(10, exp);

//   cout << "" << endl;
//   cout << num << " = " << man << " x10^{"<<exp<<"} "   << endl;
  

}

void AliLatexTable::StripLatex(TString &text, TString format) {

  // Strip latex away for ascii and html printing. Replaces latex
  // command with corresponding text/tags

  text.ReplaceAll("\\cdot", "x");
  text.ReplaceAll("$", "");
  if (format == "ASCII") {
    text.ReplaceAll("\\right>", ">");
    text.ReplaceAll("\\left<", "<");
    text.ReplaceAll("\\rangle", ">");
    text.ReplaceAll("\\langle", "<");
    text.ReplaceAll("\\pm", "+-");
  } else if (format == "HTML" || format == "TWIKI") {
    // the & is used to tokenize... Have to cheat here
    text.ReplaceAll("\\right>", "AMPERrang;");
    text.ReplaceAll("\\left<",  "AMPERlang;");
    text.ReplaceAll("\\rangle", "AMPERrang;");
    text.ReplaceAll("\\langle", "AMPERlang;");
    text.ReplaceAll("\\pm",     "AMPERplusmn;");
  } 
  if(text.Contains("multicolumn")) {
    //    cout << "col " << text.Data() << endl;
    // in case of multicol span, replace first column with content and
    // add n empty cols
    TObjArray * array = TPRegexp("multicolumn\\{([^}]*)\\}\\{[^}]*\\}\\{([^]]*)\\}").MatchS(text); // Added double \\ because gcc 4 triggers hard error on unknown escape sequence. Hope it still works...
    const TString content = ((TObjString *)array->At(2))->GetString();
    Int_t nspan   = ((TObjString *)array->At(1))->GetString().Atoi();
    text = content;
    // cout << "ns " << nspan <<  ((TObjString *)array->At(1))->GetString() << endl;
    // cout << "t " << text << endl;
    
    for(Int_t ispan = 1; ispan < nspan; ispan++){
      text+=" & ";
    }
    //    cout << "t " << text << endl;
    
  
  }

  text.ReplaceAll("\\mathrm", "");

  text.ReplaceAll("\\", "");
  text.Strip(TString::EStripType(1), ' ');

}

void AliLatexTable::LoadTeXFromFileAndPrintASCII(const char * filename) {

  // opens a file containing only a latex table and prints it on screen as ASCII

  ifstream file (filename);
  if (!file.is_open()) {
    AliError(Form("Cannot open file %s", filename));
  }
  TString line;
  while (line.ReadLine(file)) {
    if (line.Contains("begin") && line.Contains("tabular")) {
      // We need to get and parse the format
      //      TPRegexp re("\\begin\\{tabular\\}\\{([^\\}]*)\\}");
      TPRegexp re(".*begin{tabular}{(.*)}");
      TObjArray * arr = re.MatchS(line);
      if (arr->GetLast() > 0){ 
	//	cout << "Size: " << arr->GetSize() << " " << arr->GetLast() << endl;
	
	TString subStr = ((TObjString *)arr->At(1))->GetString();
	subStr.ReplaceAll("|","");
	subStr.ReplaceAll(" ","");
	subStr.ReplaceAll("\t","");
	//	subStr.ReplaceAll(" ","");	
	//	cout << subStr.Data() << " " << subStr.Length()<< endl;
	fNcol = subStr.Length();
	delete arr;
      }
    }

    // Skip stuff we don't parse
    if (line.Contains("begin")) continue;
    if (line.Contains("end")) continue;
    if (line.Contains("tabular")) continue;
    // add line
    InsertCustomRow(line.Data());
    
  }
  PrintTable("ASCII");
}
