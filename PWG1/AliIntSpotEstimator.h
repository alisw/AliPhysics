#ifndef ALIINTSPOTESTIMATOR_H
#define ALIINTSPOTESTIMATOR_H

#include <TNamed.h>
#include <TObjArray.h>
#include <TMath.h>
#include <TH2F.h>
#include <TAxis.h>

class TH1;
class TCanvas;
class AliESDEvent;
class AliESDtrack;
class AliESDVertex;
class AliVertexerTracks;

class AliIntSpotEstimator : public TNamed {
  //
 public:
  //
  AliIntSpotEstimator(Bool_t initDef=kFALSE);
  AliIntSpotEstimator(const char* name, Double_t outcut=1e-4,
		      Int_t nPhiBins=12,Int_t nestb=500,
		      Double_t estmin=-2e-2,Double_t estmax=6e-2,
		      Int_t ntrBins=10,Int_t ntMn=2,Int_t ntMx=32,
		      Int_t nPBins=14,Double_t pmn=0.2,Double_t pmx=3.);
  ~AliIntSpotEstimator();
  AliIntSpotEstimator &operator += (const AliIntSpotEstimator &src);  
  //
  void          InitEstimators(Int_t nPhiBins=12,Int_t nestb=500,
			       Double_t estmin=-2e-2,Double_t estmax=6e-2,
			       Int_t ntrBins=10,Int_t ntMn=2,Int_t ntMx=32,
			       Int_t nPBins=14,Double_t pmn=0.2,Double_t pmx=3.);
  //
  Bool_t        ProcessEvent(const AliESDEvent* esd, const AliESDVertex* vtx=0);
  Bool_t        ProcessEvent(const TObjArray* tracks);
  //
  Double_t      GetIPCenter(Int_t id)                 const  {return fIPCenterStat>0 ? fIPCenter[id]/fIPCenterStat:0;}
  Double_t      GetIPSigma(Int_t phibin=0)            const;
  Double_t      GetVtxSigma(int ntr)                  const;
  Double_t      GetDCASigma(double p)                 const;
  //
  Int_t         GetEventsAccepted()                   const  {return fIPCenterStat;} 
  Int_t         GetEventsProcessed()                  const  {return fEvProc;} 
  //
  void          SetOutlierCut(Double_t v=1e-4)               {fOutlierCut = v;}
  //
  TH2F*         GetHistoIP()                          const  {return fEstimIP;}
  TH2F*         GetHistoVtx()                         const  {return fEstimVtx;}
  TH2F*         GetHistoTrc()                         const  {return fEstimTrc;}
  AliVertexerTracks* GetVertexer()                    const  {return fVertexer;}
  //
  Int_t         GetNPhiBins()    const {return !IsValid() ? 0:fEstimIP->GetXaxis()->GetNbins();}
  Int_t         GetNTrackBins()  const {return !IsValid() ? 0:fEstimVtx->GetXaxis()->GetNbins();}
  Int_t         GetMinTracks()   const {return !IsValid() ? 0:TMath::Nint(fEstimVtx->GetXaxis()->GetXmin());}
  Int_t         GetMaxTracks()   const {return !IsValid() ? 0:TMath::Nint(fEstimVtx->GetXaxis()->GetXmax());}
  Int_t         GetNPBins()      const {return !IsValid() ? 0:fEstimTrc->GetXaxis()->GetNbins();}
  Double_t      GetTrackMinP()   const {return !IsValid() ? 0:1./fEstimTrc->GetXaxis()->GetXmax();}
  Double_t      GetTrackMaxP()   const {return !IsValid() ? 0:1./fEstimTrc->GetXaxis()->GetXmin();}
  //
  TCanvas*      CreateReport(const char* outname = 0);
  virtual void  Print(Option_t *opt="")               const;
  virtual void  Clear(Option_t *opt="");
  virtual Long64_t Merge(TCollection *coll);
  static Double_t CalcMean(TH1* histo, Double_t ctfact);
  //
 protected:
  Bool_t        IsValid()       const {return fEstimIP!=0;}
  Bool_t        IsZero(Double_t v,Double_t thresh=1e-15) const {return TMath::Abs(v)<thresh;}
  Bool_t        ProcessTracks();
  Bool_t        ProcessIPCenter(const AliESDVertex* vtx);
  Bool_t        ProcessEstimators(const AliESDEvent* esd);
  void          UpdateEstimators(double rvD, double rtD, double nTracks, double pTrack, double phiTrack);
  //
 private :
  AliIntSpotEstimator(const AliIntSpotEstimator &src);
  AliIntSpotEstimator &operator=(const AliIntSpotEstimator &src);  
  //
 protected:
  //
  Int_t         fEvProc;                         // number of events processed
  Int_t         fIPCenterStat;                   // number of events used for fIPCenter
  Double_t      fOutlierCut;                     // cut on outliers
  Double_t      fIPCenter[3];                    // mean IP position XYZ
  //
  TH2F         *fEstimIP;                        // distribution of IP estimator
  TH2F         *fEstimVtx;                       // distribution of VTRes estimator
  TH2F         *fEstimTrc;                       // distribution of DCA res estimator
  //
  AliVertexerTracks* fVertexer;                  //! vertex fitter
  TObjArray    *fTracks;                         //! storage for processed tracks
  //
 ClassDef(AliIntSpotEstimator,1)
};


#endif
