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
  enum {S_ITSNCLS,S_MINNTRACKS,S_PTTRACK,S_DECAYLENGTH,S_TRANSVERSEIP,S_TRACKCHI2,S_MAXDCAJETTRACK};

  AliHFJetTaggingIP();
  ~AliHFJetTaggingIP();
  void SetImpactParameterThreshold(Double_t threshold);
  void SetImpactParameterThreshold(TF1 * threshld_fct);
  void SetEvent(AliVEvent * bEvent){fEvent = bEvent;};
  void SetParticleContainer(AliParticleContainer * particles){fParticles = particles;};
  void InitTrackSelectionParams(Double_t * params=0x0);
  Bool_t GetJetDiscriminator(AliEmcalJet * jet,Double_t *discriminator,Bool_t * check_discr);
  Bool_t GetJetDiscriminatorQualityClass(Int_t qtyclass, AliEmcalJet * jet,Double_t *discriminator,Bool_t * check_discr);
  Bool_t DoTagging(AliEmcalJet * jetrec,Double_t  &n0,Double_t  &n1,Double_t  &n2);
  void SetAnalysisTypeAOD(Bool_t IsAnaTypeAOD = kTRUE){fAnaTypeAOD=IsAnaTypeAOD;}

  
private:
  Bool_t GetImpactParameter	(AliVTrack * bTrack, Double_t *bSign, Double_t *bIp2d);
  Bool_t PassedCuts		(AliVTrack * bTrack, Double_t ip);
  Bool_t IsV0DaughterRadius(AliVTrack * track,Double_t &Radius);
  Bool_t IsInQualityClass(AliAODTrack * track, int qclass);
  Double_t GetDecayLength 	(AliVTrack * bTrack);
  Double_t CalculateTrackProbability(AliVTrack * bTrack);
  static bool mysort( const std::pair<Int_t, Double_t>& i, const  std::pair<Int_t, Double_t>& j );

  
  
  Bool_t        fUseThresholdFuction;
  Bool_t        fAnaTypeAOD;
  Double_t      fSelectionCuts[7]; //[7]
  Double_t      fCurrentDCA ;
  Double_t      fThreshold;
  std::vector<Double_t>        fDiscriminators;//
  std::vector<unsigned int>    fJetIndices;// //Indices to get matched MC jet
  
  TF1 *          fThresholdFuction;//!  
  AliVEvent     *fEvent;//! 
  AliVVertex	*fVertex;//!
  AliVVertex	*fVertexRecalculated;//!
  AliEmcalJet	*fJet;//!

  AliParticleContainer * fParticles ; //! Particle container containing AliVTracks
  
  ClassDef(AliHFJetTaggingIP,1);
};
#endif
