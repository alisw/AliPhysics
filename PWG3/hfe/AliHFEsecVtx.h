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
class AliStack;
class AliHFEtrackFilter;
class AliHFEpairs;
class AliHFEsecVtxs;

//________________________________________________________________
class AliHFEsecVtx : public TObject {

  public: 
    AliHFEsecVtx();
    AliHFEsecVtx(const AliHFEsecVtx &p); // copy constructor
    AliHFEsecVtx &operator=(const AliHFEsecVtx &); // assignment operator
    virtual ~AliHFEsecVtx();

    void CreateHistograms(TList *qaList);

    void Process(AliVTrack *track);

    Bool_t HasMCData() const { return TestBit(kHasMCData); };
    Bool_t IsAODanalysis() const { return TestBit(kAODanalysis); };
    Bool_t IsESDanalysis() const { return !TestBit(kAODanalysis); };

		void SetHasMCData(Bool_t hasMCdata = kTRUE) { SetBit(kHasMCData,hasMCdata); };
    void SetAODAnalysis() { SetBit(kAODanalysis, kTRUE); };
    void SetESDAnalysis() { SetBit(kAODanalysis, kFALSE); };
    void SetEvent(AliESDEvent* const ESD){fESD1=ESD;};    // set ESD pointer
    void SetEventAOD(AliAODEvent* const AOD){fAOD1=AOD;}; // set ESD pointer
    void SetStack(AliStack* const stack){fStack=stack;};  // set stack pointer
    void SetMCArray(TClonesArray* const mcarry){fMCArray=mcarry;} // set mcarray pointer
    void SetUseMCPID(Bool_t usemcpid){fUseMCPID=usemcpid;};


    Int_t GetMCPID(AliESDtrack *track); // return MC pid
		Int_t GetMCPDG(AliVTrack *track);   // return MC pid
    Int_t GetPairOriginESD(AliESDtrack* track1, AliESDtrack* track2); // return pair origin as a pdg code
    Int_t GetPairOriginAOD(AliAODTrack* track1, AliAODTrack* track2); // return pair origin as a pdg code
    Int_t GetPairCode(const AliVTrack* const track1, const AliVTrack* const track2); // return corresponding pair code to pdg code
    Int_t GetElectronSource(Int_t mclabel); // return origin of the electron
    Int_t GetPDG(AliVTrack *track);     // return pdg 
		void GetESDPID(AliESDtrack *track, Int_t &recpid, Double_t &recprob); //return esd pid likelihood
    void GetPrimaryCondition();

    TClonesArray *HFEpairs();
    TClonesArray *HFEsecvtxs();

    void AddHFEpairToArray(const AliHFEpairs* const pair);
    void AddHFEsecvtxToArray(const AliHFEsecVtxs* const secvtx);

    void InitHFEpairs();
    void InitHFEsecvtxs();

    void DeleteHFEpairs();
    void DeleteHFEsecvtxs();

    Bool_t SingleTrackCut(AliESDtrack* track1) const; // single track cut
    void PairAnalysis(AliVTrack* ESDtrack1, AliVTrack* ESDtrack2, Int_t index2); // do e-h analysis
    void RunSECVTX(AliVTrack *track); // run secondary vertexing algorithm

    void MakeContainer(); // make containers

  protected:
    void Init();
    void FindSECVTXCandid(AliVTrack *track);
    void CalcSECVTXProperty(AliVTrack* track1, AliVTrack* track2, AliVTrack* track3); // calculated distinctive variables


  private:
    enum{
      kHasMCData = BIT(15),     // bitset for mc data usage
      kAODanalysis = BIT(16)    // bitset for aod analysis
    };
    enum {kAll, kDirectCharm, kDirectBeauty, kBeautyCharm, kGamma, kPi0, kElse, kBeautyGamma, kBeautyPi0, kBeautyElse}; // electron origin 
    enum {kCharm=4, kBeauty=5}; // quark flavor

    AliHFEtrackFilter *fFilter; // filter Tracks to combine the signal track with
    AliESDEvent* fESD1; // ESD pointer             
    AliAODEvent* fAOD1; // AOD pointer             
    AliStack* fStack;   // stack pointer              

    Bool_t fUseMCPID;   // if use MC pid 

    TString fkSourceLabel[10]; // electron source label

    Int_t fNparents;           // number of heavy hadrons to be considered
    Int_t fParentSelect[2][7]; // heavy hadron species
		Double_t fPtRng[7];        // pt ranges to consider pt dependant dca cut
		Double_t fDcaCut[6];       // pt dependant dca cut

    Int_t fNoOfHFEpairs;       // number of e-h pairs  
    Int_t fNoOfHFEsecvtxs;     // number of secondary vertexes

    TClonesArray *fHFEpairs;   //! Array of pair 
    TClonesArray *fHFEsecvtxs; //! Array of secondary vertexes 
    TClonesArray *fMCArray;    //! mc array pointer

		Double_t fPVx;          // primary vertex copy x 
		Double_t fPVy;          // primary vertex copy y
    Double_t fCosPhi;       // cos of opening angle of two pair vertex
    Double_t fSignedLxy;    // signed Lxy of secondary vertex
    Double_t fKFchi2;       // chi2 of secondary vertex
    Double_t fInvmass;      // invariant mass of secondary vertex
    Double_t fInvmassSigma; // invariant mass sigma of secondary vertex
    Double_t fKFip;         // impact parameter of secondary vertex track

    THnSparseF *fPairQA;    // qa histos for pair analysis 
    THnSparseF *fSecvtxQA;  // qa histos for secvtx
    TList *fSecVtxList;     // list for secondary vertexing outputs

  ClassDef(AliHFEsecVtx,0);
};

#endif
