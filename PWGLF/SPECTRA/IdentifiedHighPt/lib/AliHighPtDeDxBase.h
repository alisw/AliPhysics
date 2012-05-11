#ifndef ALIHIGHPTDEDXBASE_H
#define ALIHIGHPTDEDXBASE_H

#include <TNamed.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TF1.h>
#include <TCanvas.h>

class AliHighPtDeDxBase : public TNamed {
 public:
  AliHighPtDeDxBase(); // default constructor  
  AliHighPtDeDxBase(const char* name, const char* title); // named constructor    
  virtual ~AliHighPtDeDxBase(); // default destructor
  
  void MakeNice1dHisto(TH1* hist, TVirtualPad* c1);
  void MakeNice2dHisto(TH2* hist, TVirtualPad* c1, Bool_t colz=kFALSE);

  TCanvas* DrawPhiCutHistograms();
  TCanvas* FindCanvas(const Char_t* canvasName, 
		      Int_t xwidth, Int_t ywidth);
  TCanvas* DrawNice(TH1* hist, const Char_t* canvasName,
		    Int_t xwidth, Int_t ywidth, const Char_t* option);
  
  Double_t GetEtaLow()  { return fEtaLow; }    
  Double_t GetEtaHigh() { return fEtaHigh; }    
  Bool_t   IsMc()       { return fIsMc; } 

  virtual void SetIsMc        (Bool_t  value) { fIsMc = value; } 

  virtual void SetUseRunCut   (Bool_t  value) { fUseRunCut = value; } 
  virtual void SetRun	      (Int_t   value) { fRun = value; }	      
  virtual void SetUseEtaCut(Bool_t  value);
  virtual void SetUseEtaCutAbs(Bool_t  value);
  virtual void SetEtaLow      (Double_t value) { fEtaLow = value; }    
  virtual void SetEtaHigh     (Double_t value) { fEtaHigh = value; }   
  virtual void SetUseFilterCut(Bool_t  value) { fUseFilterCut = value; }
  virtual void SetFilter      (Int_t   value) { fFilter = value; }    
  virtual void SetUsePhiCut   (Bool_t  value) { fUsePhiCut = value; } 
  virtual void SetPhiCutLow   (TF1*    value) { fPhiCutLow = value; } 
  virtual void SetPhiCutHigh  (TF1*    value) { fPhiCutHigh = value; }
  
  virtual void SetEventVtxStatus(Int_t value)   { fEventVtxStatus = value; }
  virtual void SetEventVtxStatusMc(Int_t value) { fEventVtxStatusMc = value; }
  virtual void SetEventRun(Int_t value)         { fEventRun = value; }
  virtual void SetEventMag(Double_t value)      { fEventMag = value; }
  virtual void SetEventTrigger(Int_t value)     { fEventTrigger = value; }
  virtual void SetTrackCharge(Int_t value)      { fTrackCharge = value; }
  virtual void SetTrackEta(Double_t value)      { fTrackEta = value; }
  virtual void SetTrackP(Double_t value)        { fTrackP = value; }  
  virtual void SetTrackPt(Double_t value)       { fTrackPt = value; } 
  virtual void SetTrackFilter(Int_t value)      { fTrackFilter = value; }
  virtual void SetTrackPhi(Double_t value)      { fTrackPhi = value; }
  virtual void SetTrackDeDx(Double_t value)     { fTrackDeDx = value; }
  virtual void SetTrackNcl(Int_t value)         { fTrackNcl = value; }
  virtual void SetTrackBeta(Double_t value)     { fTrackBeta = value; }
  
  virtual void SetTrackPidMc(Int_t value)       { fTrackPidMc = value; }
  virtual void SetTrackPrimaryMc(Int_t value)   { fTrackPrimaryMc = value; }
  
  TF1* GetStandardPhiCutLow();
  TF1* GetStandardPhiCutHigh();

  TH1D* GetHistNevents()    { return hNevents; };
  TH1D* GetHistVtxStatus()  { return hVtxStatus; };
  TH1D* GetHistPt()         { return hPt; };
  TH1D* GetHistEta()        { return hEta; };
  TProfile* GetHistMeanPt() { return hMeanPt; };
  TH3F* GetHistNclVsPhiVsPtBefore() { return hNclVsPhiVsPtBefore; }
  TH3F* GetHistNclVsPhiVsPtAfter()  { return hNclVsPhiVsPtAfter; }
  
  virtual void   FillEventInfo();
  virtual void   FillTrackInfo(Float_t weight=1);
  virtual Bool_t EventAccepted();
  virtual Bool_t TrackAccepted();
  virtual void Init(Int_t nPtBins, Double_t* ptBins);

  void Print(Option_t* option="") const;
  
 protected: // so we can use the variables in derived classes
  // cut ranges
  Bool_t   fIsMc;                 // kTRUE is data is Mc (e.g. pid info is available)
  Bool_t   fUseRunCut;            // kTRUE to only sture data from fRun
  Int_t    fRun;
  Bool_t   fUseEtaCut;            // kTRUE to require fEtaLow < eta < fEtaHigh
  Bool_t   fUseEtaCutAbs;         // kTRUE to require fEtaLow < eta < fEtaHigh
  Double_t fEtaLow;
  Double_t fEtaHigh;
  Bool_t   fUseFilterCut;         // kTRUE to require that fFilter (bit) was set
  Int_t    fFilter;
  Bool_t   fUsePhiCut;            // kTRUE to use pT > 2.0 edge cut
  TF1*     fPhiCutLow;            // Function for phi cut low
  TF1*     fPhiCutHigh;           // Function for phi cut low

  // Actual values for the event and track
  Int_t    fEventVtxStatus;      //! event vtx status
  Int_t    fEventVtxStatusMc;    //! event vtx status based on Mc vtx
  Int_t    fEventRun;            //! event run
  Double_t fEventMag;            //! event magnetic field
  Double_t fEventTrigger;        //! 1 is triggered / 0 if not (Mc only)

  Int_t    fTrackCharge;         //! charge (+1 or -1)
  Double_t fTrackEta;            //! eta
  Double_t fTrackP;              //! p
  Double_t fTrackPt;             //! pt
  Int_t    fTrackFilter;         //! filter
  Double_t fTrackPhi;            //! phi
  Double_t fTrackPhiCutVariable; //! derived phi variable used for cutting
  Double_t fTrackDeDx;           //! dE/dx for track
  Int_t    fTrackNcl;            //! N cl used for dE/dx for track
  Double_t fTrackBeta;           //! Beta

  Int_t    fTrackPidMc;          //! Mc pid information
  Int_t    fTrackPrimaryMc;      //! Mc information if track is primary

 private:
  // Debug histograms
  TH1D*     hVtxStatus;
  TH1D*     hNevents;
  TH1D*     hPhi;
  TH1D*     hEta;
  TH1D*     hPt;
  TProfile* hMeanPt;
  TH3F*     hNclVsPhiVsPtBefore;
  TH3F*     hNclVsPhiVsPtAfter;
  ClassDef(AliHighPtDeDxBase, 4)  // AliHighPtDeDxBase information
    };

#endif
	
