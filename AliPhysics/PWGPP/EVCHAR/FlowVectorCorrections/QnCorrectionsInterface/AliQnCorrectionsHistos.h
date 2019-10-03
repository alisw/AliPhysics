/*
***********************************************************
    Histogram manager for event plane framework QA
    Contact: Jaap Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
    Histogram manager inspired from PWGDQ/dielectron/AliDielectronHistos by J.Wiechula
    and based on work of Ionut-Cristian Arsene
***********************************************************
*/

#ifndef ALIQNCORRECTIONS_HISTOS_H
#define ALIQNCORRECTIONS_HISTOS_H

#include <iostream>
#include <fstream>

#include <TObjArray.h>
#include <TFile.h>
#include <TDirectory.h>
#include <THashList.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include <THn.h>
#include <TIterator.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TAxis.h>
#include <TVector.h>
#include <TTimeStamp.h>


class TH1;
class TString;
class TList;


class AliQnCorrectionsHistos : public TNamed{
public:

    
  AliQnCorrectionsHistos();
  AliQnCorrectionsHistos(const char* name, const char* title);
  virtual ~AliQnCorrectionsHistos();
  void FillHistClass(const Char_t* className, Float_t* values);// Bool_t* usedvalues);
  void AddHistClass(const Char_t* histClass);
  void AddHistogram(const Char_t* histClass,
		    const Char_t* name, const Char_t* title, Bool_t isProfile,
        Int_t nXbins, Double_t xmin, Double_t xmax, Int_t varX,
		    Int_t nYbins=0, Double_t ymin=0, Double_t ymax=0, Int_t varY=kNothing,
		    Int_t nZbins=0, Double_t zmin=0, Double_t zmax=0, Int_t varZ=kNothing,
                    const Char_t* xLabels="", const Char_t* yLabels="", const Char_t* zLabels="",
                    Int_t varT=kNothing, Int_t varW=kNothing);
  void AddHistogram(const Char_t* histClass,
		    const Char_t* name, const Char_t* title, Bool_t isProfile,
                    Int_t nXbins, Double_t* xbins, Int_t varX,
		    Int_t nYbins=0, Double_t* ybins=0x0, Int_t varY=kNothing,
		    Int_t nZbins=0, Double_t* zbins=0x0, Int_t varZ=kNothing,
                    const Char_t* xLabels="", const Char_t* yLabels="", const Char_t* zLabels="",
                    Int_t varT=kNothing, Int_t varW=kNothing);
  void AddHistogram(const Char_t* histClass,
                    const Char_t* name, const Char_t* title,
                    Int_t nDimensions, Int_t* vars,
                    Int_t* nBins, Double_t* xmin, Double_t* xmax,
                    TString* axLabels=0x0,
                    Int_t varW=kNothing, Int_t axisA=-1, Double_t * newbinsA=0x0,
		    Int_t axisB=-1, Double_t * newbinsB=0x0, Int_t axisC=-1, Double_t * newbinsC=0x0);
  void AddHistogram(const Char_t* histClass,
                    const Char_t* name, const Char_t* title,
                    Int_t nDimensions, Int_t* vars,
                    TArrayD* binLimits,
                    TString* axLabels=0x0,
                    Int_t varW=kNothing);
  void MakeAxisLabels(TAxis* ax, const Char_t* labels);
  void WriteOutput(TFile* save);

  //void SetHistList(TObjArray* histList) const {fHistLists=histList;}
  TObjArray* HistList()  {return fHistLists;}

  static AliQnCorrectionsHistos* Instance() {if(!fgEventPlaneHistos) fgEventPlaneHistos = new AliQnCorrectionsHistos(); return fgEventPlaneHistos;}
  
  static AliQnCorrectionsHistos* fgEventPlaneHistos;

private:
    
  static const Int_t kNothing=-1;
  AliQnCorrectionsHistos(const AliQnCorrectionsHistos &hist);
  AliQnCorrectionsHistos& operator = (const AliQnCorrectionsHistos &hist);


  Int_t fBinsAllocated;

  TObjArray* fHistLists;   // main histogram list for the current running process
  TFile* fHistFile;      // pointer to a TFile opened for reading


  ClassDef(AliQnCorrectionsHistos,1)

};  


#endif
