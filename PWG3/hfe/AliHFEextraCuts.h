/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
//
// Extra Cuts from the ALICE HFE Group
// Container for cuts which are currently not implemented by
// the Correction Framework
//
#ifndef ALIHFEEXTRACUTS_H
#define ALIHFEEXTRACUTS_H

// #ifndef ALICFCUTBASE_H
#include "AliCFCutBase.h"
// #endif

class TList;

class AliVEvent;
class AliVParticle;
class AliVTrack;

class AliHFEextraCuts : public AliCFCutBase{
  public:
    typedef enum{
      kFirst = 0,
      kSecond = 1,
      kBoth = 2,
      kNone = 3,
      kAny = 4,
      kExclusiveSecond = 5
    } ITSPixel_t;
    typedef enum{
      kFound = 0,
      kFoundIter1 = 1,
      kCrossedRows = 2
    } ETPCclusterDef_t;
    typedef enum{
      kFoundOverFindable = 0,
      kFoundOverFindableIter1 = 1,
      kFoundOverCR = 2,
      kCROverFindable = 3
    } ETPCclrDef_t;
    AliHFEextraCuts(const Char_t *name, const Char_t *title);
    AliHFEextraCuts(const AliHFEextraCuts &c);
    AliHFEextraCuts &operator=(const AliHFEextraCuts &c);
    virtual ~AliHFEextraCuts();
    
    virtual Bool_t IsSelected(TObject *o);
    virtual Bool_t IsSelected(TList *) { return kTRUE; };
    virtual void SetRecEventInfo(const TObject *event);

    inline void SetClusterRatioTPC(Double_t ratio, ETPCclrDef_t def);
    inline void SetRequireITSpixel(ITSPixel_t pixel);
    inline void SetMinImpactParamR(Double_t impactParam);
    inline void SetMaxImpactParamR(Double_t impactParam);
    inline void SetMinImpactParamZ(Double_t impactParam);
    inline void SetMaxImpactParamZ(Double_t impactParam);
    inline void SetMinHFEImpactParamR(Float_t ipcutParam[4], Bool_t issigmacut);
    inline void SetMinTrackletsTRD(Int_t minTracklets);
    inline void SetMinNClustersTPC(Int_t minclusters, ETPCclusterDef_t def);
    void SetTOFPID(Bool_t tofPid) { tofPid ? SETBIT(fRequirements, kTOFPID) : CLRBIT(fRequirements, kTOFPID); }
    void SetTOFMISMATCH(Bool_t tofMismatch) { tofMismatch ? SETBIT(fRequirements, kTOFmismatch) : CLRBIT(fRequirements, kTOFmismatch); }
    void SetTPCPIDCleanUp(Bool_t tpcPIDCleanUp) { tpcPIDCleanUp ? SETBIT(fRequirements, kTPCPIDCleanUp) : CLRBIT(fRequirements, kTPCPIDCleanUp); }
    void SetMaxImpactParameterRpar(Bool_t maxImpactParameterRpar) { maxImpactParameterRpar ? SETBIT(fRequirements, kMaxImpactParameterRpar) : CLRBIT(fRequirements, kMaxImpactParameterRpar); }
    void SetFractionOfTPCSharedClusters(Double_t fractionShared) { fFractionTPCShared= fractionShared; SETBIT(fRequirements, kTPCfractionShared); }

    void SetCheckITSstatus(Bool_t check) { fCheck = check; };
    void SetDebugLevel(Int_t level) { fDebugLevel = level; };

    Bool_t GetCheckITSstatus() const { return fCheck; };
    Int_t GetDebugLevel() const { return fDebugLevel; };
    void GetHFEImpactParameters(AliVTrack *track, Double_t &dcaxy, Double_t &dcansigmaxy); // temporary moved from protected to publich for IP QA 
    
  protected:
    virtual void AddQAHistograms(TList *qaList);
    Bool_t CheckRecCuts(AliVTrack *track);
    Bool_t CheckMCCuts(AliVParticle * /*track*/) const;
    Bool_t CheckITSstatus(Int_t itsStatus) const;
    void FillQAhistosRec(AliVTrack *track, UInt_t when);
//     void FillQAhistosMC(AliMCParticle *track, UInt_t when);
    void FillCutCorrelation(ULong64_t survivedCut);
    void PrintBitMap(Int_t bitmap);
    
    // Getter Functions for ESD/AOD compatible mode
    Int_t GetITSstatus(AliVTrack *track, Int_t layer);
    UInt_t GetTPCncls(AliVTrack *track);
    UInt_t GetTPCnclusdEdx(AliVTrack *track);
    Bool_t GetTPCCountSharedMapBitsAboveThreshold(AliVTrack *track);
    Double_t GetTPCclusterRatio(AliVTrack *track);
    void GetImpactParameters(AliVTrack *track, Float_t &radial, Float_t &z);
    //void GetHFEImpactParameters(AliVTrack *track, Double_t &dcaxy, Double_t &dcansigmaxy);
    void GetHFEImpactParameterCuts(AliVTrack *track, Double_t &hfeimpactRcut, Double_t &hfeimpactnsigmaRcut);
    void GetMaxImpactParameterCutR(AliVTrack *track, Double_t &maximpactRcut);
    Float_t GetTPCsharedClustersRatio(AliVTrack *track);

  private:
    typedef enum{
      kMinImpactParamR = 0,
      kMaxImpactParamR = 1,
      kMinImpactParamZ = 2,
      kMaxImpactParamZ = 3,
      kClusterRatioTPC = 4,
      kMinTrackletsTRD = 5,
      kPixelITS = 6,
      kMinHFEImpactParamR = 7,
      kMinHFEImpactParamNsigmaR = 8,
      kMinNClustersTPC = 9,
      kTPCfractionShared = 10,
      kTOFPID = 11,
      kTOFmismatch = 12,
      kTPCPIDCleanUp = 13,
      kEMCALmatch = 14,
      kMaxImpactParameterRpar = 15,
      kNcuts = 16
    } Cut_t;
    enum{
      //
      // Common Constants
      //
      kBeforeCuts =0,
      kAfterCuts = 1
    };
    static const Int_t fgkNQAhistos;    // Number of QA histos
    AliVEvent *fEvent;                //! working event
    ULong64_t fCutCorrelation;		    // Cut Correlation
    ULong64_t fRequirements;		      // Cut Requirements
    Float_t fImpactParamCut[4];		    // Impact Parmameter Cut
    Float_t fIPcutParam[4];		        // Parmameter of impact parameter cut parametrization
    UInt_t fMinNClustersTPC;          // Minimum TPC clusters cut
    Float_t fClusterRatioTPC;		      // Ratio of findable vs. found clusters in TPC
    UChar_t fMinTrackletsTRD;		      // Min. Number of Tracklets inside TRD
    UChar_t fPixelITS;                // Cut on ITS Pixels
    UChar_t fTPCclusterDef;           // TPC cluster definition Bitmap
    UChar_t fTPCclusterRatioDef;      // TPC cluster ratio definition Bitmap
    Double_t  fFractionTPCShared;     // Cut on fraction of shared clusters

    Bool_t  fCheck;                     // check
    TList *fQAlist;			//! Directory for QA histograms
    Int_t   fDebugLevel;                // Debug Level
  
  ClassDef(AliHFEextraCuts, 1)      // Additional cuts implemented by the ALICE HFE group
};

//__________________________________________________________
void AliHFEextraCuts::SetClusterRatioTPC(Double_t ratio, ETPCclrDef_t def) {
  SETBIT(fRequirements, kClusterRatioTPC);
  SETBIT(fTPCclusterRatioDef, def);
  fClusterRatioTPC = ratio; 
}

//__________________________________________________________
void AliHFEextraCuts::SetRequireITSpixel(ITSPixel_t pixel) {
  SETBIT(fRequirements, kPixelITS);
  fPixelITS = pixel; 
}

//__________________________________________________________
void AliHFEextraCuts::SetMinImpactParamR(Double_t impactParam){
  SETBIT(fRequirements, kMinImpactParamR);
  fImpactParamCut[0] = impactParam;
}

//__________________________________________________________
void AliHFEextraCuts::SetMaxImpactParamR(Double_t impactParam){
  SETBIT(fRequirements, kMaxImpactParamR);
  fImpactParamCut[2] = impactParam;
}

//__________________________________________________________
void AliHFEextraCuts::SetMinImpactParamZ(Double_t impactParam){
  SETBIT(fRequirements, kMinImpactParamZ);
  fImpactParamCut[1] = impactParam;
}

//__________________________________________________________
void AliHFEextraCuts::SetMaxImpactParamZ(Double_t impactParam){
  SETBIT(fRequirements, kMaxImpactParamZ);
  fImpactParamCut[3] = impactParam;
}

//__________________________________________________________
void AliHFEextraCuts::SetMinHFEImpactParamR(Float_t ipcutParam[4], Bool_t isSigmacut){
  if(isSigmacut) SETBIT(fRequirements, kMinHFEImpactParamNsigmaR);
  else SETBIT(fRequirements, kMinHFEImpactParamR);
  fIPcutParam[0]=ipcutParam[0];
  fIPcutParam[1]=ipcutParam[1];
  fIPcutParam[2]=ipcutParam[2];
  fIPcutParam[3]=ipcutParam[3];
}

//__________________________________________________________
void AliHFEextraCuts::SetMinTrackletsTRD(Int_t minTracklets){
  SETBIT(fRequirements, kMinTrackletsTRD);
  fMinTrackletsTRD = minTracklets;
}

//__________________________________________________________
void AliHFEextraCuts::SetMinNClustersTPC(Int_t minClusters, ETPCclusterDef_t tpcdef){
  SETBIT(fRequirements, kMinNClustersTPC);
  SETBIT(fTPCclusterDef, tpcdef);
  fMinNClustersTPC = minClusters;
}
#endif
