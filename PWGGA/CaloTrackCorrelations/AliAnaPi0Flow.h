#ifndef ALIANAPI0FLOW_H
#define ALIANAPI0FLOW_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliAnaPi0Flow
/// \ingroup CaloTrackCorrelationsAnalysis 
/// \brief Flow analysis
///
/// Class for the analysis of pi0 and eta flow.
///
/// \author Qyie Shou, <qiye.shou@cern.ch> Wuhan 
///_________________________________________________________________________

// Root
class TList;
class TH3F;
class TH2F;
class THnSparse;
class TObjString;
class AliAODEvent;
class AliESDEvent;
class AliAODPWG4Particle;
class AliOADBContainer;
class AliEMCALGeometry;
class AliESDEvent;
class AliESDtrack;
class AliESDCaloCells;
class AliAODEvent;
class AliAODCaloCells;
class AliVCluster;
class AliCentrality;
class AliEPFlattener;

// Analysis
#include "AliAnaCaloTrackCorrBaseClass.h"

class AliAnaPi0Flow : public AliAnaCaloTrackCorrBaseClass 
{
  
 public:   
  
  AliAnaPi0Flow();
  virtual ~AliAnaPi0Flow();
  
  TObjString             *GetAnalysisCuts();  
  TList                  *GetCreateOutputObjects(); 
  void                    Print(const Option_t * opt) const;
  void                    MakeAnalysisFillHistograms();
  void                    InitParameters();
  
  void                    IsPHOSCali(Bool_t e) {isPhosCali = e;}
  void                    IsCentFlat(Bool_t e) {isCentFlat = e;}

 private:
  
   Bool_t                 IsCentAccepted();
   Int_t                  GetInternalRunNum(Int_t runnumber);
   void                   GetVZEROEventPlane(Bool_t isFlatten);
   Double_t               ApplyFlatteningTPC(Double_t phi, Double_t c);
   Double_t               ApplyFlatteningV0A(Double_t phi, Double_t c);
   Double_t               ApplyFlatteningV0C(Double_t phi, Double_t c);

   Bool_t                 isPhosCali;
   Bool_t                 isCentFlat;

   AliVEvent             *fInputEvent;
   AliEventplane         *fEventPlane;
   Double_t               fCentrality;
   Int_t                  fRunNumber;
   Int_t                  fInternalRunNum;
   AliOADBContainer      *fFlatContainer;

   AliEPFlattener        *fTPCFlat;
   AliEPFlattener        *fV0AFlat;
   AliEPFlattener        *fV0CFlat;
   Double_t               fEPTPC;
   Double_t               fEPTPCResolution;
   Double_t               fEPV0;
   Double_t               fEPV0A;
   Double_t               fEPV0C;
   Double_t               fEPV0AR;
   Double_t               fEPV0CR;
   Double_t               fEPV0R;
   Double_t               fEPV0AR4;
   Double_t               fEPV0AR5;
   Double_t               fEPV0AR6;
   Double_t               fEPV0AR7;
   Double_t               fEPV0CR0;
   Double_t               fEPV0CR1;
   Double_t               fEPV0CR2;
   Double_t               fEPV0CR3;

   //
   // hists
   //
   TH1D                  *fHistStatCentrality;
   TH1D                  *fHistStatCentralityCorrected;
   TH1I                  *fHistStatRunNum;
   TH2F                  *fHistEPTPC;
   TH2F                  *fHistEPTPCResolution;
   TH2F                  *fHistEPV0;
   TH2F                  *fHistEPV0A;
   TH2F                  *fHistEPV0C;
   TH2F                  *fHistEPV0AR;
   TH2F                  *fHistEPV0CR;
   TH2F                  *fHistEPV0R;
   TH2F                  *fHistEPV0AR4;
   TH2F                  *fHistEPV0AR7;
   TH2F                  *fHistEPV0CR0;
   TH2F                  *fHistEPV0CR3;
   TH2F                  *fHistEPTPCFlatten;
   TH2F                  *fHistEPV0AFlatten;
   TH2F                  *fHistEPV0CFlatten;
   TH2F                  *fHistEPDiffV0A_V0CR0;
   TH2F                  *fHistEPDiffV0A_V0CR3;
   TH2F                  *fHistEPDiffV0CR0_V0CR3;
   TH2F                  *fHistEPDiffV0C_V0AR4;
   TH2F                  *fHistEPDiffV0C_V0AR7;
   TH2F                  *fHistEPDiffV0AR4_V0AR7;
   TH2F                  *fHistEPDiffV0AR_V0CR;
   TH2F                  *fHistClusterEtaPhi;
   TH2F                  *fHistClusterEN;
   TH2F                  *fHistClusterEtN;
   TH2F                  *fHistClusterEM02;
   TH2F                  *fHistClusterEtM02;

   //
   // effective data
   //
   THnSparse             *fDataV0;
   THnSparse             *fDataV0A;
   THnSparse             *fDataV0C;
   THnSparse             *fDataTPC;
  
  /// Copy constructor not implemented.
  AliAnaPi0Flow(              const AliAnaPi0Flow & api0) ;
   
  /// Assignment operator not implemented.
  AliAnaPi0Flow & operator = (const AliAnaPi0Flow & api0) ;
  
  /// \cond CLASSIMP
  ClassDef(AliAnaPi0Flow,31) ;
  /// \endcond
  
} ;

#endif //ALIANAPI0FLOW_H



