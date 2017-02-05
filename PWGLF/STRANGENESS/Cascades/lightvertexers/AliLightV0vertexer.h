#ifndef AliLightV0vertexer_H
#define AliLightV0vertexer_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//------------------------------------------------------------------
//                    V0 Vertexer Class
//            reads tracks writes out V0 vertices
//   Origin: Iouri Belikov, IReS, Strasbourg, Jouri.Belikov@cern.ch
//------------------------------------------------------------------

#include "TObject.h"

class TTree;
class AliESDEvent;

//_____________________________________________________________________________
class AliLightV0vertexer : public TObject {
public:
    AliLightV0vertexer();
    void SetCuts(const Double_t cuts[7]);
    static void SetDefaultCuts(const Double_t cuts[7]);
    
    Int_t Tracks2V0vertices(AliESDEvent *event);
    
    void GetCuts(Double_t cuts[7]) const;
    static void GetDefaultCuts(Double_t cuts[7]);
    
    static void SetDefaultMaxEta(Double_t lMaxEta);
    static void SetDefaultMinClusters(Double_t lMaxEta);
    void SetMaxEta(Double_t lMaxEta);
    void SetMinClusters(Double_t lMaxEta);
    
    //Experimental implementation of V0 refit functionality 
    void SetDoRefit( Bool_t lDoRefit ) { fkDoRefit = lDoRefit; }
    
private:
    static
    Double_t fgChi2max;      // maximal allowed chi2
    static
    Double_t fgDNmin;        // min allowed impact parameter for the 1st daughter
    static
    Double_t fgDPmin;        // min allowed impact parameter for the 2nd daughter
    static
    Double_t fgDCAmax;       // maximal allowed DCA between the daughter tracks
    static
    Double_t fgCPAmin;       // minimal allowed cosine of V0's pointing angle
    static
    Double_t fgRmin, fgRmax; // max & min radii of the fiducial volume
    
    static Double_t fgMaxEta;       // maximum eta value for track pre-selection
    static Double_t fgMinClusters;  // minimum single-track clusters value (>=)
    
    Double_t fChi2max;      // maximal allowed chi2
    Double_t fDNmin;        // min allowed impact parameter for the 1st daughter
    Double_t fDPmin;        // min allowed impact parameter for the 2nd daughter
    Double_t fDCAmax;       // maximal allowed DCA between the daughter tracks
    Double_t fCPAmin;       // minimal allowed cosine of V0's pointing angle
    Double_t fRmin, fRmax;  // max & min radii of the fiducial volume
    
    Double_t fMaxEta;       // maximum eta value for track pre-selection
    Double_t fMinClusters;  // minimum single-track clusters value (>=)
    
    Bool_t fkDoRefit; //improve precision with a V0 refit (+ calculate chi2)
    
    ClassDef(AliLightV0vertexer,3)  // V0 verterxer
};

inline AliLightV0vertexer::AliLightV0vertexer() :
TObject(),
fChi2max(fgChi2max),
fDNmin(fgDNmin),
fDPmin(fgDPmin),
fDCAmax(fgDCAmax),
fCPAmin(fgCPAmin),
fRmin(fgRmin),
fRmax(fgRmax),
fMaxEta(fgMaxEta),
fMinClusters(fgMinClusters),
fkDoRefit(kTRUE)
{
}

inline void AliLightV0vertexer::SetCuts(const Double_t cuts[7]) {
    fChi2max=cuts[0];
    fDNmin=cuts[1];   fDPmin=cuts[2];
    fDCAmax=cuts[3];  fCPAmin=cuts[4];
    fRmin=cuts[5];    fRmax=cuts[6];
}

inline void AliLightV0vertexer::SetDefaultCuts(const Double_t cuts[7]) {
    fgChi2max=cuts[0];
    fgDNmin=cuts[1];   fgDPmin=cuts[2];
    fgDCAmax=cuts[3];  fgCPAmin=cuts[4];
    fgRmin=cuts[5];    fgRmax=cuts[6];
}

inline void AliLightV0vertexer::GetCuts(Double_t cuts[7]) const {
    cuts[0]=fChi2max;
    cuts[1]=fDNmin;   cuts[2]=fDPmin;
    cuts[3]=fDCAmax;  cuts[4]=fCPAmin;
    cuts[5]=fRmin;    cuts[6]=fRmax;
}

inline void AliLightV0vertexer::GetDefaultCuts(Double_t cuts[7]) {
    cuts[0]=fgChi2max;
    cuts[1]=fgDNmin;   cuts[2]=fgDPmin;
    cuts[3]=fgDCAmax;  cuts[4]=fgCPAmin;
    cuts[5]=fgRmin;    cuts[6]=fgRmax; 
}

inline void AliLightV0vertexer::SetDefaultMaxEta(Double_t lMaxEta) {
    fgMaxEta = lMaxEta;
}
inline void AliLightV0vertexer::SetDefaultMinClusters(Double_t lMinClusters) {
    fgMinClusters = lMinClusters;
}
inline void AliLightV0vertexer::SetMaxEta(Double_t lMaxEta) {
    fMaxEta = lMaxEta;
}
inline void AliLightV0vertexer::SetMinClusters(Double_t lMinClusters) {
    fMinClusters = lMinClusters;
}

#endif


