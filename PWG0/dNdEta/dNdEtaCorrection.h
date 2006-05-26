/* $Id$ */

#ifndef DNDETACORRECTION_H
#define DNDETACORRECTION_H


// ------------------------------------------------------
//
// Class to handle corrections for dN/dEta measurements
//
// ------------------------------------------------------
//
// TODO:
// - add documentation
// - add status: generate or use maps
// - add functionality to set the bin sizes
// - add histograms with errors (for error visualization)
// 

#ifndef ROOT_TObject
#include "TObject.h"
#endif
#ifndef ROOT_TFile
#include "TFile.h"
#endif
#ifndef ROOT_TH2
#include "TH2.h"
#endif


class dNdEtaCorrection : public TObject
{
protected:
  
  TString  fName; 
  
  TH2F*    hEtaVsVtx_meas;
  TH2F*    hEtaVsVtx_gene;

  TH2F*    hEtaVsVtx_corr; 
  TH2F*    hEtaVsVtx_ratio;

public:
  dNdEtaCorrection(Char_t* name="dndeta_correction");

  TH2F* GetGeneratedHistogram() { return hEtaVsVtx_gene; }
  TH2F* GetMeasuredHistogram() { return hEtaVsVtx_meas; }

  void SetGeneratedHistogram(TH2F* aGeneratedHistogram) { hEtaVsVtx_gene = aGeneratedHistogram; }
  void SetMeasuredHistogram(TH2F* aMeasuredHistogram) { hEtaVsVtx_meas = aMeasuredHistogram; }

  void FillMeas(Float_t vtx, Float_t eta) {hEtaVsVtx_meas->Fill(vtx, eta);}
  void FillGene(Float_t vtx, Float_t eta) {hEtaVsVtx_gene->Fill(vtx, eta);}

  void Finish();

  void    SaveHistograms();
  Bool_t  LoadHistograms(Char_t* fileName, Char_t* dir = "dndeta_correction");
  Bool_t  LoadCorrection(Char_t* fileName, Char_t* dir = "dndeta_correction") 
    {return LoadHistograms(fileName, dir);}
  
  void DrawHistograms();
  
  void    RemoveEdges(Float_t cut=2, Int_t nBinsVtx=0, Int_t nBinsEta=0);
  
  Float_t GetCorrection(Float_t vtx, Float_t eta) 
    {return hEtaVsVtx_corr->GetBinContent(hEtaVsVtx_corr->FindBin(vtx,eta));}
  
  ClassDef(dNdEtaCorrection,0)
};

#endif
