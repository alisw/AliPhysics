#ifndef UtilFiles_h
#define UtilFiles_h

#include "TString.h"

class TFile;
class TKey;
class TH1F;
class TH1D;
class TH2F;
class TH2I;
class TH3I;

//_________________________________________________________________________________________________
TFile * GetFile(TString fname, TString opt = "READ");

//_________________________________________________________________________________________________
void GetListFromFile(TFile *fin, TString name, TList *& lin);

//_________________________________________________________________________________________________
void GetListFromDirectory(TDirectory *dir, TString name, TList *& lin);

//_________________________________________________________________________________________________
void GetDirectoryFromFile(TFile *fin, TString name, TDirectory *& dir);

//_________________________________________________________________________________________________
void GetHistogram(TFile *fin, const TString hname, TH1F *& histo);///Gets the histogram from the file

//_________________________________________________________________________________________________
void GetHistogram(TFile *fin, const TString hname, TH1D *& histo);///Gets the histogram from the file

//_________________________________________________________________________________________________
void GetHistogram(TFile *fin, const TString hname, TH2I *& histo);///Gets the histogram from the file

//_________________________________________________________________________________________________
void GetHistogram(TFile *fin, const TString hname, TH2F *& histo);///Gets the histogram from the file

//_________________________________________________________________________________________________
void GetHistogram(TList *lin, const TString hname, TH1F *& histo);///Gets the histogram from the list

//_________________________________________________________________________________________________
void GetHistogram(TList *lin, const TString hname, TH1D *& histo);///Gets the histogram from the list

//_________________________________________________________________________________________________
void GetHistogram(TList *lin, const TString hname, TH2F *& histo);///Gets the histogram from the list

//_________________________________________________________________________________________________
void GetHistogram(TList *lin, const TString hname, TH2I *& histo);///Gets the histogram from the list

//_________________________________________________________________________________________________
void GetHistogram(TList *lin, const TString hname, TH3I *& histo);///Gets the histogram from the list

//_________________________________________________________________________________________________
void GetHistogram(TDirectory *dir, const TString hname, TH1F *& histo);///Gets the histogram from the directory

//_________________________________________________________________________________________________
TList *ReduceList(TList *lin, const TString criteria);///Macro to produce multiple lists from one -> Useful for writing to file

//_________________________________________________________________________________________________
TList *FormListFromFile(TFile *fin, const TString criteria = "", const TString checklists = "");
  
#endif
