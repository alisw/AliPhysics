
#ifndef ALIHEPDATAPARSER_H
#define ALIHEPDATAPARSER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                      Implementation of   Class AliHEPDataParser
//
//  This class is used to save the content of hisograms/graphs in the
//  HEP data format and viceversa
//  Author: Michele Floris, CERN
//-------------------------------------------------------------------------

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

protected:

  TH1 * fHistStat; // statistical errors (hist)
  TH1 * fHistSyst; // systematic errors (hist) 
  TGraph * fGraphStat; // statistical errors (hist)
  TGraph * fGraphSyst; // systematic errors (hist)  
  TObjArray * fHEPDataFileLines;// TClones array of TObjString
  TString fValueName; // title for the y axis on the ascii file


  ClassDef(AliHEPDataParser, 1);
    
private:


  AliHEPDataParser(const AliHEPDataParser&);
  AliHEPDataParser& operator=(const AliHEPDataParser&);
};

#endif
