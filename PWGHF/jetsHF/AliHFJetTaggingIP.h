#ifndef ALIHFJETTAGGINGIP_H
#define ALIHFJETTAGGINGIP_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
#include <utility>
#include <vector>

class AliVEvent;
class AliVVertex;
class AliVTrack;
class AliAODTrack;
class AliEmcalJet;
class AliParticleContainer;
class TF1;

class AliHFJetTaggingIP
{
public:
    enum { S_ITSNCLS, S_MINNTRACKS, S_PTTRACK, S_DECAYLENGTH, S_TRANSVERSEIP, S_TRACKCHI2, S_MAXDCAJETTRACK,S_MAXDCAZ };

    AliHFJetTaggingIP();
    virtual ~AliHFJetTaggingIP();
    void SetImpactParameterThreshold(Double_t threshold);
    void SetImpactParameterThreshold(TF1* threshld_fct);
    void SetEvent(AliVEvent* bEvent)
    {
	fEvent = bEvent;
    };
    void SetParticleContainer(AliParticleContainer* particles)
    {
	fParticles = particles;
    };
    void InitTrackSelectionParams(Double_t* params = 0x0);
    Bool_t GetJetDiscriminator(AliEmcalJet* jet, Double_t* discriminator, Bool_t* check_discr);
    Bool_t
    GetJetDiscriminatorQualityClass(Int_t qtyclass, AliEmcalJet* jet, Double_t* discriminator, Bool_t* check_discr);
    Bool_t DoTagging(AliEmcalJet* jetrec, Double_t& n0, Double_t& n1, Double_t& n2);
    void SetAnalysisTypeAOD(Bool_t IsAnaTypeAOD = kTRUE){fAnaTypeAOD = IsAnaTypeAOD;}
    void SetUseSignDefinitionAtlas(Bool_t val = kTRUE){fUseSignAtlas = val;};
    void SetUseSignificance(Bool_t val = kTRUE){fUseSignificance = val;};
    void SetUse3DIP(Bool_t val = kTRUE){fUse3DsIP = val;};
    virtual AliVTrack* GetCurrentTrack(int i = 0);
    virtual Double_t GetCurrentTrackDCAz(int i = 0);

private:
    Bool_t GetImpactParameter(AliVTrack* bTrack, Double_t* bSign, Double_t* bIp2d, Double_t* zDCA);
    Bool_t PassedCuts(AliVTrack* bTrack, Double_t ip,Double_t ipz);
    Bool_t IsV0DaughterRadius(AliVTrack* track, Double_t& Radius);
    Bool_t IsInQualityClass(AliAODTrack* track, int qclass, double ip, double ipz);
    Double_t GetDecayLength(AliVTrack* bTrack);
    Double_t CalculateTrackProbability(AliVTrack* bTrack);
    Double_t GetSignAtlasDefinition(Double_t* xDCA, Double_t* pDCA, Double_t* xVtx, Double_t* pJet);
    static bool mysort(const std::pair<Int_t, Double_t>& i, const std::pair<Int_t, Double_t>& j);
    void SetCurrentDCAz();
    Bool_t fUseThresholdFuction;
    Bool_t fAnaTypeAOD;
    Bool_t fUseSignAtlas;
    Bool_t fUseSignificance;
    Bool_t fUse3DsIP;
    Int_t fCurrentTrack[3];        //[3];
    Double_t fCurrentTrackDCAz[3]; //[3];

    Double_t fSelectionCuts[8]; //[7]
    Double_t fCurrentDCA;
    Double_t fThreshold;
    std::vector<Double_t> fDiscriminators;              //
    std::vector<std::pair<Int_t, Double_t> > fTrackDCA; //
    std::vector<unsigned int> fJetIndices; // //Indices to get matched MC jet
	
    TF1* fThresholdFuction;           //!
    AliVEvent* fEvent;                //!
    AliVVertex* fVertex;              //!
    AliVVertex* fVertexRecalculated;  //!
    AliEmcalJet* fJet;                //!
    AliParticleContainer* fParticles; //! Particle container containing AliVTracks

    ClassDef(AliHFJetTaggingIP, 1);
};
#endif
