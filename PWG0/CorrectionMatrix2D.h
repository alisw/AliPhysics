#ifndef CORRECTIONMATRIX2D_H
#define CORRECTIONMATRIX2D_H

// ------------------------------------------------------
//
// Class to handle 2d - corrections 
//
// ------------------------------------------------------

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif
#ifndef ROOT_TFile
#include "TFile.h"
#endif
#ifndef ROOT_TH2
#include "TH2.h"
#endif
#ifndef ROOT_TError
#include "TError.h"
#endif
#ifndef ROOT_TCanvas
#include "TCanvas.h"
#endif

class CorrectionMatrix2D : public TNamed 
{
public:
  CorrectionMatrix2D(Char_t* name="CorrectionMatrix", Char_t* title="",
		     Int_t nBinX=10, Float_t Xmin=0., Float_t Xmax=10.,
		     Int_t nBinY=10, Float_t Ymin=0., Float_t Ymax=10.);
  
  CorrectionMatrix2D(Char_t* name, Char_t* title,
		     Int_t nBinX, Float_t *X, Int_t nBinY, Float_t *Y);

  virtual ~CorrectionMatrix2D(); 

  TH2F* GetGeneratedHistogram() { return fhGene; }
  TH2F* GetMeasuredHistogram()  { return fhMeas; }

  void SetGeneratedHistogram(TH2F* agene) { fhGene = agene; }
  void SetMeasuredHistogram(TH2F* ameas)  { fhMeas  = ameas; }

  void FillMeas(Float_t ax, Float_t ay) {fhMeas->Fill(ax,ay);}
  void FillGene(Float_t ax, Float_t ay) {fhGene->Fill(ax,ay);}

  void Finish();  
                  
  void SetAxisTitles(Char_t* titleX="", Char_t* titleY="");
  
  void SaveHistograms();
  void DrawHistograms();  

  Bool_t  LoadHistograms(Char_t* fileName, Char_t* dir);
  Bool_t  LoadCorrection(Char_t* fileName, Char_t* dir) {return LoadHistograms(fileName, dir);}
  
  void    RemoveEdges(Float_t cut=2, Int_t nBinsVtx=0, Int_t nBinsEta=0);
  
  Float_t GetCorrection(Float_t ax, Float_t ay) {return fhCorr->GetBinContent(fhCorr->FindBin(ax,ay));}
  
protected:
  
  TH2F*    fhMeas;
  TH2F*    fhGene;

  TH2F*    fhCorr; 
  TH2F*    fhRatio;


  ClassDef(CorrectionMatrix2D,0)
};

#endif

