#ifndef CORRECTIONMATRIX2D_H
#define CORRECTIONMATRIX2D_H

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
#ifndef ROOT_TCanvas
#include "TCanvas.h"
#endif
#ifndef ROOT_TError
#include "TError.h"
#endif

class CorrectionMatrix2D : public TObject 
{
protected:
  
  TString  fName; 
  
  TH2F*    hmeas;
  TH2F*    hgene;

  TH2F*    hcorr; 
  TH2F*    hratio;

public:
  CorrectionMatrix2D(Char_t* name="CorrectionMatrix");

  TH2F* GetGeneratedHistogram() { return hgene; }
  TH2F* GetMeasuredHistogram()  { return hmeas; }

  void SetGeneratedHistogram(TH2F* agene) { hgene = agene; }
  void SetMeasuredHistogram(TH2F* ameas) { hmeas = ameas; }


  void FillMeas(Float_t x, Float_t y) {hmeas->Fill(x,y);}
  void FillGene(Float_t x, Float_t y) {hgene->Fill(x,y);}

  void Finish();
  
  //Char_t *name , Char_t *title , 
               
  void SetHist(Char_t* title ,Int_t nBinX=10, Float_t Xmin=0., Float_t Xmax=10.,
	       Int_t nBinY=10, Float_t Ymin=0., Float_t Ymax=10.);
  void SetHist(Char_t* title ,Int_t nBinX, Float_t *X, Int_t nBinY, Float_t *Y);
   
  void SetHistTitle(Char_t* titleX=" ", Char_t* titleY=" ");
  
  void SaveHistograms();
  void DrawHistograms();  

  Bool_t  LoadHistograms(Char_t* fileName, Char_t* dir);
  Bool_t  LoadCorrection(Char_t* fileName, Char_t* dir) {return LoadHistograms(fileName, dir);}
  
  void    RemoveEdges(Float_t cut=2, Int_t nBinsVtx=0, Int_t nBinsEta=0);
  
  Float_t GetCorrection(Float_t x, Float_t y) {return hcorr->GetBinContent(hcorr->FindBin(x,y));}
  
  ClassDef(CorrectionMatrix2D,0)
};

#endif

