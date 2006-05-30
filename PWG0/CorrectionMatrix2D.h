// ------------------------------------------------------
//
// Class to handle 2d-corrections. 
//
// ------------------------------------------------------
//
// TODO:
//
// - change the finish method (should not be called finish)
// - add option in draw method
// 

/* $Id$ */

#ifndef CORRECTIONMATRIX2D_H
#define CORRECTIONMATRIX2D_H


#include <TNamed.h>
#include <TH2.h>

class TFile;
class TCanvas;
class AliLog;

class CorrectionMatrix2D : public TNamed 
{
public:
  CorrectionMatrix2D(const CorrectionMatrix2D& c);
  CorrectionMatrix2D(Char_t* name, Char_t* title,
		     Int_t nBinX=10, Float_t Xmin=0., Float_t Xmax=10.,
		     Int_t nBinY=10, Float_t Ymin=0., Float_t Ymax=10.);
  
  CorrectionMatrix2D(Char_t* name, Char_t* title,
		     Int_t nBinX, Float_t *X, Int_t nBinY, Float_t *Y);

  virtual ~CorrectionMatrix2D(); 

  CorrectionMatrix2D& operator=(const CorrectionMatrix2D& corrMatrix);

  virtual void Copy(TObject& c) const;

  TH2F* GetGeneratedHistogram() { return fhGene; }
  TH2F* GetMeasuredHistogram()  { return fhMeas; }

  TH1F* Get1DCorrection(Char_t* opt="x");

  void SetGeneratedHistogram(TH2F* agene) { fhGene = agene; }
  void SetMeasuredHistogram(TH2F* ameas)  { fhMeas = ameas; }

  void FillMeas(Float_t ax, Float_t ay) {fhMeas->Fill(ax,ay);}
  void FillGene(Float_t ax, Float_t ay) {fhGene->Fill(ax,ay);}
  
  void Divide();  
  
  virtual Long64_t Merge(TCollection* list);
                
  void SetAxisTitles(Char_t* titleX="", Char_t* titleY="");
  
  void SaveHistograms();
  void DrawHistograms();  

  Bool_t  LoadHistograms(Char_t* fileName, Char_t* dir = ".");
  
  void    RemoveEdges(Float_t cut=2, Int_t nBinsX=0, Int_t nBinsY=0);
  
  Float_t GetCorrection(Float_t ax, Float_t ay) {return fhCorr->GetBinContent(fhCorr->FindBin(ax,ay));}
  
protected:
  
  TH2F*    fhMeas;  // histogram of measured particles (or tracks)
  TH2F*    fhGene;  // histogram of generated particles

  TH2F*    fhCorr;  // correction histogram (ratio generated/measured)

  ClassDef(CorrectionMatrix2D,1)
};

#endif

