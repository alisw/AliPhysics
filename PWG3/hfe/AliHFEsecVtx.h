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
//  Secondary vertexing construction Class
//  Construct secondary vertex from Beauty hadron with electron and
//  hadrons, then apply selection criteria
//

#ifndef ALIHFESECVTX_H
#define ALIHFESECVTX_H

#ifndef ROOT_TObject
//#include <TObject.h>
#endif

#ifndef ROOT_THnSparse
#include <THnSparse.h>
#endif

class TH1F;
class TH2F;
class TString;
class AliESDEvent;
class AliAODEvent;
class AliVTrack;
class AliESDtrack;
class AliAODTrack;
class AliMCEvent;
class AliHFEtrackFilter;
class AliHFEpairs;
class AliHFEsecVtxs;
class AliKFParticle;
class AliHFEmcQA;

//________________________________________________________________
class AliHFEsecVtx : public TObject {

  public: 
    AliHFEsecVtx();
    AliHFEsecVtx(const AliHFEsecVtx &p); // copy constructor
    AliHFEsecVtx &operator=(const AliHFEsecVtx &); // assignment operator
    virtual ~AliHFEsecVtx();

    void CreateHistograms(TList * const qaList);

    Bool_t Process(AliVTrack *track);

    Bool_t HasMCData() const { return TestBit(kHasMCData); };
    Bool_t IsAODanalysis() const { return TestBit(kAODanalysis); };
    Bool_t IsESDanalysis() const { return !TestBit(kAODanalysis); };

		void SetHasMCData(Bool_t hasMCdata = kTRUE) { SetBit(kHasMCData,hasMCdata); };
    void SetAODAnalysis() { SetBit(kAODanalysis, kTRUE); };
    void SetESDAnalysis() { SetBit(kAODanalysis, kFALSE); };
    void SetEvent(AliESDEvent* const ESD){fESD1=ESD;};    // set ESD pointer
    void SetEventAOD(AliAODEvent* const AOD){fAOD1=AOD;}; // set ESD pointer
    void SetMCEvent(AliMCEvent* const mcEvent){fMCEvent=mcEvent;};  // set stack pointer
    void SetMCArray(TClonesArray* const mcarry){fMCArray=mcarry;} // set mcarray pointer
    void SetUseMCPID(Bool_t usemcpid){fUseMCPID=usemcpid;};
    void SetMCQA(AliHFEmcQA * const mcqa){fMCQA=mcqa;};    // set mcqa pointer


    Int_t GetMCPID(const AliESDtrack *track); // return MC pid
		Int_t GetMCPDG(const AliVTrack *track);   // return MC pid
    Int_t GetPairOriginESD(AliESDtrack* track1, AliESDtrack* track2); // return pair origin as a pdg code
    Int_t GetPairOriginAOD(AliAODTrack* track1, AliAODTrack* track2); // return pair origin as a pdg code
    Int_t GetPairCode(const AliVTrack* const track1, const AliVTrack* const track2); // return corresponding pair code to pdg code
    Int_t GetElectronSource(Int_t mclabel); // return origin of the electron
    Int_t GetPDG(AliVTrack *track);     // return pdg 
    void GetESDPID(const AliESDtrack *track, Int_t &recpid, Double_t &recprob); //return esd pid likelihood
    void GetPrimaryCondition();
    void RecalcPrimvtx(Int_t nkftrk, const Int_t * const, const AliKFParticle * const); //recalculate primary vertex

    TClonesArray *HFEpairs();
    TClonesArray *HFEsecvtxs();

    void AddHFEpairToArray(const AliHFEpairs* const pair);
    void AddHFEsecvtxToArray(const AliHFEsecVtxs* const secvtx);

    void InitHFEpairs();
    void InitHFEsecvtxs();

    void DeleteHFEpairs();
    void DeleteHFEsecvtxs();

    void PairAnalysis(AliVTrack* ESDtrack1, AliVTrack* ESDtrack2, Int_t index2); // do e-h analysis
    void RunSECVTX(AliVTrack *track); // run secondary vertexing algorithm

    void MakeContainer(); // make containers
    void MakeHistos(Int_t step); // make histograms for different steps
    void FillHistos(Int_t step, const AliESDtrack *track); // fill histograms for different steps

  protected:
    void Init();
    void FindSECVTXCandid(AliVTrack *track);
    void CalcSECVTXProperty(AliVTrack* track1, AliVTrack* track2, AliVTrack* track3); // calculated distinctive variables
    void CalcSECVTXProperty(AliVTrack* track1, AliVTrack* track2, AliVTrack* track3, AliVTrack* track4); // calculated distinctive variables

    void Fill4TrkSECVTX(AliVTrack* track, Int_t ipair, Int_t jpair, Int_t kpair);
    void Fill3TrkSECVTX(AliVTrack* track, Int_t ipair, Int_t jpair);
    void Fill2TrkSECVTX(AliVTrack* track, const AliHFEpairs *pair);

  private:
    enum{
      kHasMCData = BIT(15),     // bitset for mc data usage
      kAODanalysis = BIT(16)    // bitset for aod analysis
    };
    enum {kAll, kDirectCharm, kDirectBeauty, kBeautyCharm, kGamma, kPi0, kElse, kBeautyGamma, kBeautyPi0, kBeautyElse, kMisID}; // electron origin 
    enum {kCharm=4, kBeauty=5}; // quark flavor

    AliHFEtrackFilter *fFilter; // filter Tracks to combine the signal track with
    AliESDEvent* fESD1; // ESD pointer             
    AliAODEvent* fAOD1; // AOD pointer             
    AliMCEvent* fMCEvent;   // MCEvent pointer              

    AliHFEmcQA* fMCQA;  // mcqa pointer

    Bool_t fUseMCPID;   // if use MC pid 

    TString fkSourceLabel[10]; // electron source label

    Int_t fNparents;           // number of heavy hadrons to be considered
    Int_t fParentSelect[2][7]; // heavy hadron species
		Double_t fPtRng[7];        // pt ranges to consider pt dependant dca cut
		Double_t fDcaCut[6];       // pt dependant dca cut

    Int_t fNoOfHFEpairs;       // number of e-h pairs  
    Int_t fNoOfHFEsecvtxs;     // number of secondary vertexes
    Int_t fArethereSecVtx;     // checker

    TClonesArray *fHFEpairs;   //! Array of pair 
    TClonesArray *fHFEsecvtxs; //! Array of secondary vertexes 
    TClonesArray *fMCArray;    //! mc array pointer

		Double_t fPVx;          // primary vertex copy x 
		Double_t fPVy;          // primary vertex copy y
		Double_t fPVx2;         // recalculated primary vertex x 
		Double_t fPVy2;         // recalculated primary vertex y
    Double_t fCosPhi;       // cos of opening angle of two pair vertex
    Double_t fSignedLxy;    // signed Lxy of secondary vertex
    Double_t fSignedLxy2;   // signed Lxy of secondary vertex based on recalculated primary vertex
    Double_t fKFchi2;       // chi2 of secondary vertex
    Double_t fInvmass;      // invariant mass of secondary vertex
    Double_t fInvmassSigma; // invariant mass sigma of secondary vertex
    Double_t fKFip;         // impact parameter of secondary vertex track
    Double_t fKFip2;        // impact parameter of secondary vertex track based on recalculated primary vertex

    Int_t fNsectrk2prim;    // # of secvtx tracks contributing to primvtx calculation

    Double_t fVtxchi2Tightcut; // pair vertex chi2 cut
    Double_t fVtxchi2Loosecut; // secvtx vertex chi2 cut

    THnSparseF *fPairQA;    // qa histos for pair analysis 
    THnSparseF *fSecvtxQA;  // qa histos for secvtx
    TList *fSecVtxList;     // list for secondary vertexing outputs

  ClassDef(AliHFEsecVtx,0);
};

#endif
