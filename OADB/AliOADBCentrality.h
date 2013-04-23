#ifndef AliOADBCentrality_H
#define AliOADBCentrality_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     OADB class for run dependent centrality scaling factors and 
//     data for centrality determination
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------

#include <TNamed.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>


class AliOADBCentrality : public TNamed {

 public :
  AliOADBCentrality();
  AliOADBCentrality(char* name);
  virtual ~AliOADBCentrality();
  Float_t V0MScaleFactor()   const  {return fV0MScaleFactor;}
  Float_t SPDScaleFactor()   const  {return fSPDScaleFactor;}
  Float_t TPCScaleFactor()   const  {return fTPCScaleFactor;}
  Float_t V0MScaleFactorMC() const  {return fV0MScaleFactorMC;}

  Float_t V0MSPDOutlierPar0()      const  {return fV0MSPDOutlierPar0      ;}
  Float_t V0MSPDOutlierPar1()      const  {return fV0MSPDOutlierPar1      ;}
  Float_t V0MTPCOutlierPar0()      const  {return fV0MTPCOutlierPar0      ;}
  Float_t V0MTPCOutlierPar1()      const  {return fV0MTPCOutlierPar1      ;}

  Float_t V0MSPDSigmaOutlierPar0() const  {return fV0MSPDSigmaOutlierPar0 ;}
  Float_t V0MSPDSigmaOutlierPar1() const  {return fV0MSPDSigmaOutlierPar1 ;}
  Float_t V0MSPDSigmaOutlierPar2() const  {return fV0MSPDSigmaOutlierPar2 ;}
  Float_t V0MTPCSigmaOutlierPar0() const  {return fV0MTPCSigmaOutlierPar0 ;}
  Float_t V0MTPCSigmaOutlierPar1() const  {return fV0MTPCSigmaOutlierPar1 ;}
  Float_t V0MTPCSigmaOutlierPar2() const  {return fV0MTPCSigmaOutlierPar2 ;}

  Float_t V0MZDCOutlierPar0()      const  {return fV0MZDCOutlierPar0      ;}
  Float_t V0MZDCOutlierPar1()      const  {return fV0MZDCOutlierPar1      ;}
  Float_t V0MZDCEcalOutlierPar0()  const  {return fV0MZDCEcalOutlierPar0  ;}
  Float_t V0MZDCEcalOutlierPar1()  const  {return fV0MZDCEcalOutlierPar1  ;}

  Float_t ZVCut()       const {return fZVCut;}
  Float_t OutliersCut() const {return fOutliersCut;}
  Bool_t UseScaling()   const {return fUseScaling;}
  Bool_t UseCleaning()  const {return fUseCleaning;}

  TH1F*   V0hist()         const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultV0M_percentile")));}
  TH1F*   V0Ahist()        const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultV0A_percentile")));}
  TH1F*   V0Chist()        const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultV0C_percentile")));}
  TH1F*   V0Eqhist()       const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultV0MEq_percentile")));}
  TH1F*   V0AEqhist()      const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultV0AEq_percentile")));}
  TH1F*   V0CEqhist()      const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultV0CEq_percentile")));}
  TH1F*   TPChist()        const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultTRK_percentile")));}
  TH1F*   CNDhist()        const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultCND_percentile")));}
  TH1F*   SPDhist()        const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultCL1_percentile")));}
  TH1F*   FMDhist()        const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultFMD_percentile")));}
  TH1F*   ZNAhist()        const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultZNA_percentile")));}
  TH1F*   ZNChist()        const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultZNC_percentile")));}
  TH1F*   ZPAhist()        const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultZPA_percentile")));}
  TH1F*   ZPChist()        const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultZPC_percentile")));}
  TH2F*   ZEMvsZDChist()   const  {return ((TH2F*) (Hists2D()->FindObject("fHOutMultZEMvsZDC")));}

  TH1F*   NPAhist()        const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultNPA_percentile")));}
  TH1F*   NPAhistDPM()     const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultNPADPM_percentile")));}

  TH1F*   V0histtrue()     const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultV0Mtrue_percentile")));}
  TH1F*   V0Ahisttrue()    const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultV0Atrue_percentile")));}
  TH1F*   V0Chisttrue()    const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultV0Ctrue_percentile")));}
  TH1F*   V0Eqhisttrue()   const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultV0MEqtrue_percentile")));}
  TH1F*   V0AEqhisttrue()  const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultV0AEqtrue_percentile")));}
  TH1F*   V0CEqhisttrue()  const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultV0CEqtrue_percentile")));}
  TH1F*   TPChisttrue()    const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultTRKtrue_percentile")));}
  TH1F*   CNDhisttrue()    const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultCNDtrue_percentile")));}
  TH1F*   SPDhisttrue()    const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultCL1true_percentile")));}
  TH1F*   FMDhisttrue()    const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultFMDtrue_percentile")));}
  TH1F*   ZNAhisttrue()    const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultZNAtrue_percentile")));}
  TH1F*   ZNChisttrue()    const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultZNCtrue_percentile")));}
  TH1F*   ZPAhisttrue()    const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultZPAtrue_percentile")));}
  TH1F*   ZPChisttrue()    const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultZPCtrue_percentile")));}
 
  TH1F*   V0histtrueDPM()     const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultV0MtrueDPM_percentile")));}
  TH1F*   V0AhisttrueDPM()    const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultV0AtrueDPM_percentile")));}
  TH1F*   V0ChisttrueDPM()    const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultV0CtrueDPM_percentile")));}
  TH1F*   V0EqhisttrueDPM()   const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultV0MEqtrueDPM_percentile")));}
  TH1F*   V0AEqhisttrueDPM()  const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultV0AEqtrueDPM_percentile")));}
  TH1F*   V0CEqhisttrueDPM()  const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultV0CEqtrueDPM_percentile")));}
  TH1F*   TPChisttrueDPM()    const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultTRKtrueDPM_percentile")));}
  TH1F*   CNDhisttrueDPM()    const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultCNDtrueDPM_percentile")));}
  TH1F*   SPDhisttrueDPM()    const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultCL1trueDPM_percentile")));}
  TH1F*   FMDhisttrueDPM()    const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultFMDtrueDPM_percentile")));}
  TH1F*   ZNAhisttrueDPM()    const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultZNAtrueDPM_percentile")));}
  TH1F*   ZNChisttrueDPM()    const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultZNCtrueDPM_percentile")));}
  TH1F*   ZPAhisttrueDPM()    const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultZPAtrueDPM_percentile")));}
  TH1F*   ZPChisttrueDPM()    const  {return ((TH1F*) (Hists1D()->FindObject("fHOutMultZPCtrueDPM_percentile")));}
  TList*  Hists1D()        const  {return f1DHistos;}
  TList*  Hists2D()        const  {return f2DHistos;}

  void    SetScaleFactors(Float_t v0m, Float_t spd, Float_t tpc, Float_t v0mMC)
  {fV0MScaleFactor = v0m; fSPDScaleFactor = spd; fTPCScaleFactor = tpc; fV0MScaleFactorMC = v0mMC;}
  
  void    SetOutlierV0MSPDFactors(Float_t a1, Float_t a2, Float_t a3, Float_t a4, Float_t a5)
  {fV0MSPDOutlierPar0=a1;fV0MSPDOutlierPar1=a2;fV0MSPDSigmaOutlierPar0=a3;fV0MSPDSigmaOutlierPar1=a4;fV0MSPDSigmaOutlierPar2=a5;}
  
  void    SetOutlierV0MTPCFactors(Float_t a1, Float_t a2, Float_t a3, Float_t a4, Float_t a5)
  {fV0MTPCOutlierPar0=a1;fV0MTPCOutlierPar1=a2;fV0MTPCSigmaOutlierPar0=a3;fV0MTPCSigmaOutlierPar1=a4;fV0MTPCSigmaOutlierPar2=a5;}

  void    SetOutlierV0MZDCFactors(Float_t a1, Float_t a2)
  {fV0MZDCOutlierPar0=a1;fV0MZDCOutlierPar1=a2;}
  
  void    SetOutlierV0MZDCEcalFactors(Float_t a1, Float_t a2)
  {fV0MZDCEcalOutlierPar0=a1;fV0MZDCEcalOutlierPar1=a2;}
  
  void    SetHistReferences(TList* l1, TList* l2)
  {f1DHistos = l1; f2DHistos = l2;}

  void SetZVCut(Float_t z)
  {fZVCut=z;}

  void SetOutliersCut(Float_t o)
  {fOutliersCut=o;}

  void SetUseScaling(Bool_t x)
  {fUseScaling=x;}

  void SetUseCleaning(Bool_t x)
  {fUseCleaning=x;}


 private:
  AliOADBCentrality(const AliOADBCentrality& cont); 
  AliOADBCentrality& operator=(const AliOADBCentrality& cont);

 private:
  Float_t fV0MScaleFactor;     // V0  scale factor
  Float_t fSPDScaleFactor;     // SPD scale factor
  Float_t fTPCScaleFactor;     // TPC scale factor
  Float_t fV0MScaleFactorMC;   // V0  scale factor for MC

  Float_t fV0MSPDOutlierPar0;  // V0-SPD outlier parameterisation Par0
  Float_t fV0MSPDOutlierPar1;  // Par1
  Float_t fV0MTPCOutlierPar0;  // Par2
  Float_t fV0MTPCOutlierPar1;  // Par3

  Float_t fV0MSPDSigmaOutlierPar0; // V0-SPD Sigma outlier parameterisation Par0
  Float_t fV0MSPDSigmaOutlierPar1; // Par1
  Float_t fV0MSPDSigmaOutlierPar2; // Par2
  Float_t fV0MTPCSigmaOutlierPar0; // Par3
  Float_t fV0MTPCSigmaOutlierPar1; // Par4
  Float_t fV0MTPCSigmaOutlierPar2; // Par5

  Float_t fV0MZDCOutlierPar0;     // V0-ZDC outlier parameterisation Par0
  Float_t fV0MZDCOutlierPar1;     // Par1
  Float_t fV0MZDCEcalOutlierPar0; // Par2
  Float_t fV0MZDCEcalOutlierPar1; // Par3

  Float_t fZVCut;                 // zV-cut
  Float_t fOutliersCut;           // outlier cuts
  Bool_t fUseScaling;             // Flag for scaling
  Bool_t fUseCleaning;            // Flag for cleaning

  TList*    f1DHistos; // Reference to list of 1D Centrality histos 
  TList*    f2DHistos; // Reference to list of 2D Centrality histos
  ClassDef(AliOADBCentrality, 3);
};

#endif
