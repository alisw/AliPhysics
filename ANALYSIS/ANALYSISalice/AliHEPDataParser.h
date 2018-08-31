
#ifndef ALIHEPDATAPARSER_H
#define ALIHEPDATAPARSER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \class AliHEPDataParser
/// \brief Implementation of Class AliHEPDataParser
///
/// This class is used to save the content of hisograms/graphs in the
/// HEP data format and viceversa. The HEP data format is not strictly
/// defined and there are many variants, the class only support a few
/// of them. More will be added as needed.  The input can be a set of
/// 2 TH1, TGraphAsymmErrors or TGraphErrors (one for the stat and one
/// for the syst error). If the second one is a null pointer, only the
/// stat error is printed. The class can also import hepdata ascii
/// file (very preliminary)
///
/// \author Michele Floris, CERN

#include "TObject.h"
#include "TString.h"

class TH1;
class TGraph;

class AliHEPDataParser : public TObject{


public:
  AliHEPDataParser();
  AliHEPDataParser(TH1 * hStat, TH1 * hSyst);
  AliHEPDataParser(TGraph * grStat, TGraph * grSyst);
  AliHEPDataParser(const char * hepfileName);

  ~AliHEPDataParser();
  
  TH1 * GetHistoStat() { return fHistStat;}
  TH1 * GetHistoSyst() { return fHistSyst;}
  TGraph * GetGraphStat() { return fGraphStat;}
  TGraph * GetGraphSyst() { return fGraphSyst;}
  void SaveHEPDataFile(const char * hepfileName, Bool_t trueUseGraphFalesUseHisto = 0);

  void SetName(const char * name) { fValueName = name;}
  void SetXaxisName (TString val) {fXaxisName = val ;}
  void SetTitle(TString val) {fTitle = val;}
  void SetReaction(TString val) {fReaction = val;}
  void SetEnergy(TString val) {fEnergy = val;}
  void SetRapidityRange(TString val) {fRapidityRange = val;}
  void SetPrecision(Int_t   val) {fPrecision = val;}


protected:

  Double_t RoundToSignificantFigures(double num, int n) ;
  TString GetFixWidthCol(Double_t number, Int_t width) ;


  TH1 * fHistStat; ///< statistical errors (hist)
  TH1 * fHistSyst; ///< systematic errors (hist)
  TGraph * fGraphStat; ///< statistical errors (hist)
  TGraph * fGraphSyst; ///< systematic errors (hist)
  TObjArray * fHEPDataFileLines;///< TClones array of TObjString
  TString fValueName; ///< title for the y axis on the ascii file
  TString fXaxisName; ///< title for the y axis
  TString fTitle; ///< title for the HEP DATA file
  TString fReaction; ///< Raction ,e.g. RE : Pb + Pb --> pbar + X
  TString fEnergy; ///< Raction ,e.g. sqrts : 2760 GeV
  TString fRapidityRange; ///< Rapidity ABS(YRAP) : 0.5'
  Int_t   fPrecision; ///< number of significant figures for rounding

  ClassDef(AliHEPDataParser, 2);
    
private:


  AliHEPDataParser(const AliHEPDataParser&);
  AliHEPDataParser& operator=(const AliHEPDataParser&);
};

#endif
