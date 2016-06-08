// -*- C++ -*-

#ifndef _ALI_LUMIMNOUS_REGION_FIT_H_
#define _ALI_LUMIMNOUS_REGION_FIT_H_

#include <ostream>

#include <TObject.h>
#include <TFile.h>
#include <TTree.h>
#include <TList.h>
#include <TString.h>
#include <TCut.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#include <TVectorD.h>
#include <TMatrixD.h>

class AliLuminousRegionFit : public TObject {
public:
  AliLuminousRegionFit(Int_t   fillNumber,
		       Int_t   minNumberOfTracks,
		       TString vtxFileName,
		       TString sepFileName)
    : TObject()
    , fFillNumber(fillNumber)
    , fMinNumberOfTracks(minNumberOfTracks)
    , f(TFile::Open(vtxFileName))
    , fL(NULL)
    , fTE(NULL)
    , fTSep(new TTree)
    , fN(0)
    , fListSave(NULL)
  {
    if (NULL == f)   AliFatal("NULL == f");
    f->cd("Vertex_Performance");

    fL  = dynamic_cast<TList*>(gDirectory->Get("cOutputVtxESD"));
    if (NULL == fL)  AliFatal("NULL == fL");

    fTE = dynamic_cast<TTree*>(fL->FindObject("fTreeBeamSpot"));
    if (NULL == fTE) AliFatal("NULL == fTE");

    if (!SetupTreeSep(sepFileName))
      AliFatal("SetupTreeSep failed");
  }

  virtual ~AliLuminousRegionFit() {
    if (fTSep) {
      fTSep->ResetBranchAddresses();
      delete fTSep; 
      fTSep = NULL;
    }
    if (fL) {
      delete fL;
      fL = NULL;
    }
    if (f) f->Close();
    if (fListSave) {
      delete fListSave;
      fListSave = NULL;
    }
  }

  Bool_t DoFit(TString  scanName,
	       Double_t tMin, Double_t tMax,
	       Int_t    scanType,
	       Double_t offset,
	       Int_t    bcSel=-1);  // <0: no selection on BCID

  Double_t MinuitFunction(const Double_t *par);
    
protected:
  Bool_t SetupTreeSep(TString sepFileName) {
    if (!fTSep->ReadFile(sepFileName, "timeStart/D:timeEnd:sep"))
      return kFALSE;
    return kTRUE;
  }
  void ComputeMoments(TTree *t,
		      const TCut& sel, const TVectorD &mu, const TMatrixDSym &cov,
		      TVectorD &x, TMatrixDSym &cx);

  TGraph*       ConfGraph   (TGraph       *g, const char *name, Int_t color=kBlack, Int_t marker=kFullDotLarge, Int_t lineWidth=1);
  TGraphErrors* ConfGraphErr(TGraphErrors *g, const char *name, Int_t color=kRed,   Int_t marker=kFullDiamond);

  static TTree* RemoveOutliersRobust(TTree *t, TVectorD &mu, TMatrixDSym &cov);
  static TTree* RemoveOutliersOld   (TTree *t, TCut &sel, TVectorD &mu, TMatrixDSym &cov);
  static TTree* SkimTTree(TTree *t, const TArrayI *idxOutliers);

  static Bool_t CutOutliers(TTree *t, Bool_t doCut, Double_t sigmaThreshold,
			    TCut &sel, TVectorD &mu, TMatrixDSym &cov);

private:
  AliLuminousRegionFit(const AliLuminousRegionFit&);
  AliLuminousRegionFit operator=(const AliLuminousRegionFit&);

  Int_t   fFillNumber;        //
  Int_t   fMinNumberOfTracks; //

  TFile *f;            //
  TList *fL;           //
  TTree *fTE;          // TTree containing unconstrained vertex data
  TTree *fTSep;        // TTree containing beam separation vs. time

  Int_t     fN;        //!
  Double_t *fX[3];     //!
  Double_t *fCov[3];   //!

  TList    *fListSave; //! list of TGraph{,Errors}s

  ClassDef(AliLuminousRegionFit, 1);
};

#endif // _ALI_LUMIMNOUS_REGION_FIT_H_
