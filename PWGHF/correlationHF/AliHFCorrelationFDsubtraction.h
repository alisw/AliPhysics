#ifndef ALIHFCORRELATIONFDSUBTRACTION_H
#define ALIHFCORRELATIONFDSUBTRACTION_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

/////////////////////////////////////////////////////////////
// Class for subtracting D-hadron correlations from secondary D from B meson decay from D-hadron correlations of prompt+secondary D mesons
//       n.b. requires the evaluation of the fraction of prompt D mesons using the same receipes used for D meson cross-section and RAA in D2H pag
//       (-> fprompt obtained from FONLL predictions of prompt and secondary D, prompt and secondary D meson efficiencies and a range of hypotheses for
//       RAA(DfromB)/RAA(promptD) (if needed)
//
// Author: A. Rossi, andrea.rossi@cern.ch
/////////////////////////////////////////////////////////////


#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>

class AliHFCorrelationFDsubtraction : public TNamed {
  
 public:

  AliHFCorrelationFDsubtraction();
  AliHFCorrelationFDsubtraction(const char* name, const char* title);
  ~AliHFCorrelationFDsubtraction();

  void SetUncorrectedHistogram(TH1D *h){fhDataUncorrected=(TH1D*)h->Clone();}
  Bool_t AddTemplateHisto(TH1D *h);
  void SetDptRange(Double_t ptmin,Double_t ptmax){fptmin=ptmin;fptmax=ptmax;}
  Bool_t Init();
  void SetNTemplates(Int_t ntempl){fnTemplates=ntempl;}
  void SetFpromptGraphFc(TGraphAsymmErrors *gr){fgrConservativeFc=(TGraphAsymmErrors*)gr->Clone();}
  void SetFpromptGraphNb(TGraphAsymmErrors *gr){fgrConservativeNb=(TGraphAsymmErrors*)gr->Clone();}
  void SetMethod(Int_t method){fmethod=method;}
  void SetSystOption(Int_t systOpt){fSystUseRMSfromFlatDistr=systOpt;}
  Int_t GetSystOption(){return fSystUseRMSfromFlatDistr;}
  void SubtractFeedDown(TH1D *hFDTempl);

  TGraphAsymmErrors* GetUncGraphFromHistos(TH1D *hRef,TH1D *hMin,TH1D *hMax);
  void CalculateEnvelope();
  TGraphAsymmErrors* GetGraphEnvelope(){return fgrEnvelope;}
  TGraphAsymmErrors* GetGraphEnvelopeRatio(){return fgrEnvelopeRatio;}
  TH1D* GetHistoEnvelopeRatioMin(){return fhEnvelopeMinRatio;}
  TH1D* GetHistoEnvelopeRatioMax(){return fhEnvelopeMaxRatio;}
  TH1D* GetHistoRelSystUncMin();
  TH1D* GetHistoRelSystUncMax();
  TH1D* GetHistoEnvelopeMin(){return fhEnvelopeMin;}
  TH1D* GetHistoEnvelopeMax(){return fhEnvelopeMax;}
  TH1D* GetCentralSubtractedPlot(){return fhDataCorrected[0];}
  TH1D* GetTemplate(Int_t i){       
    if (i>=fLastTemplAdded){Printf("Get Template: Error"); return 0;}
    else return fhTemplates[i];
  }
  TH1D* ReflectHisto(TH1D *h,Double_t scale=1.);
  TH1D* DuplicateHistoTo2piRange(TH1D *h,Double_t scale=0.5);

 private:

  Double_t fptmin;                                      //  min pt (D meson pt range)
  Double_t fptmax;                                      //  max pt (D meson pt range)
  TH1D* fhDataUncorrected;                              // Input correlation histogram
  Int_t fmethod;                                        // method: 1= fc, 2= Nb, 3= both
  TGraphAsymmErrors* fgrConservativeFc;                 // fc graph, fc method
  TGraphAsymmErrors* fgrConservativeNb;                 // fc graph, Nb method
  Int_t fnTemplates;                                   // maximum number of templates that will be included
  Int_t fLastTemplAdded;                                // counter of templates included
  TH1D **fhTemplates;                                  //  Array with template histo, size fnTemplates
  TH1D **fhDataCorrected;                               // Array with corrected histo, size 3*fnTemplates
  TH1D **fhRatioSameTemplate;                            // Array with ratio of histo based on the same template
  TH1D **fhRatioAllToCentralTempl;                                // array with ratio of all histos with respect to the central hypo
  TCanvas **fcCompare;                                            // array of canvases with comparison (template by template)
  TCanvas **fcRatio;                                          // array of canvases with ratio (templ by templ)
  Int_t fCountTempl;                                               // internal counter
  TCanvas *fcAllRatio;                                         // canva with all ratios with respect to central hypo
  TH1D* fhEnvelopeMax;                                         // envelope with max variation
  TH1D* fhEnvelopeMin;                                          // envelope with min variation
  TH1D* fhEnvelopeMaxRatio;                                     // envelope with max relative variation
  TH1D* fhEnvelopeMinRatio;                                      // envelope with min relative variation
  TGraphAsymmErrors* fgrEnvelope;                                // graph with envelope
  TGraphAsymmErrors* fgrEnvelopeRatio;                           // graph with envelope of ratios
  Int_t fSystUseRMSfromFlatDistr;                      // different option to extract the systematic uncertainty. See .cxx, GetHistoRelSystUncMin method for a description
  ClassDef(AliHFCorrelationFDsubtraction,2);
  
};
#endif
