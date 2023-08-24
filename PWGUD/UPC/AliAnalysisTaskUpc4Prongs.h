/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskUpc4Prongs_H
#define AliAnalysisTaskUpc4Prongs_H

class TClonesArray;
class TFile;
class TTree;
class TList;
class TH1;
class TList;
class AliPIDResponse;
class AliESDEvent;
class AliESDTrack;
class TBits;

#include "AliAnalysisTaskSE.h"
#include <vector>

class AliAnalysisTaskUpc4Prongs : public AliAnalysisTaskSE
{
public:
    AliAnalysisTaskUpc4Prongs(); // = delete;
    AliAnalysisTaskUpc4Prongs(const char* name);
    virtual ~AliAnalysisTaskUpc4Prongs();

    virtual void Init();
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t* option);
    virtual void Terminate(Option_t*) {};

    void SetTrigger(TString _fTriggerName) { fTriggerName = _fTriggerName; }

private:
    Bool_t Is0STPfired(Int_t*, Int_t*);

    TString fTriggerName;

    // tree
    TTree* fRhoTree;
    // tree variables and branches
    Int_t RunNum;
    UShort_t BunchCrossNumber;
    UInt_t OrbitNumber;
    UInt_t PeriodNumber;
    Float_t Mass;
    Float_t Pt;
    Short_t Q;
    Float_t Rapidity;
    Int_t V0Adecision;
    Int_t V0Cdecision;
    Int_t ADAdecision;
    Int_t ADCdecision;
    Bool_t UBAfired;
    Bool_t UBCfired;
    Bool_t VBAfired;
    Bool_t VBCfired;
    Float_t ZNAenergy;
    Float_t ZNCenergy;
    Float_t ZPAenergy;
    Float_t ZPCenergy;
    Int_t VtxContrib;
    Float_t VtxChi2, VtxNDF;
    Int_t SpdVtxContrib;
    Int_t nTracklets;
    Int_t nTracks;
    Float_t Phi;
    Float_t Vertex[3];
    Float_t SpdVertex[3];
    Float_t ZDCAtime[4];
    Float_t ZDCCtime[4];

    std::vector<Float_t> T_NumberOfSigmaITSPion;
    std::vector<Float_t> T_NumberOfSigmaITSElectron;
    std::vector<Float_t> T_NumberOfSigmaTPCPion;
    std::vector<Float_t> T_NumberOfSigmaTPCElectron;
    std::vector<Int_t>   T_TPCsignal;
    std::vector<Int_t>   T_TPCNCls;
    std::vector<Int_t>   T_ITSNCls;
    std::vector<Float_t> T_P;
    std::vector<Float_t> T_Eta;
    std::vector<Float_t> T_Phi;
    std::vector<Float_t> T_Px;
    std::vector<Float_t> T_Py;
    std::vector<Float_t> T_Pz;
    std::vector<Float_t> T_Dca0;
    std::vector<Float_t> T_Dca1;
    std::vector<Short_t> T_Q;
    std::vector<Bool_t>  T_TPCRefit;
    std::vector<Bool_t>  T_ITSRefit;
    std::vector<Bool_t>  T_HasPointOnITSLayer0;
    std::vector<Bool_t>  T_HasPointOnITSLayer1;
    std::vector<Int_t>   T_ITSModuleInner;
    std::vector<Int_t>   T_ITSModuleOuter;
    std::vector<Float_t> T_Lets_Theta;
    std::vector<Float_t> T_Lets_Phi;
    std::vector<Int_t>   T_ITSSensorNum;

    AliPIDResponse* fPIDResponse;

    AliAnalysisTaskUpc4Prongs(
        const AliAnalysisTaskUpc4Prongs&); // not implemented
    AliAnalysisTaskUpc4Prongs&
        operator=(const AliAnalysisTaskUpc4Prongs&); // not implemented

    ClassDef(AliAnalysisTaskUpc4Prongs, 2);
};

#endif
