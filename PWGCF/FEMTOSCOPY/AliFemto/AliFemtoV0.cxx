///
/// \file AliFemtoV0.cxx
///

#include "AliFemtoV0.h"
#include "phys_constants.h"

// -----------------------------------------------------------------------
AliFemtoV0::AliFemtoV0():
  fDecayLengthV0(0), fDecayVertexV0(0), fPrimaryVertex(0),
  fDcaV0Daughters(0), fDcaV0ToPrimVertex(0),
  fDcaPosToPrimVertex(0), fDcaNegToPrimVertex(0),
  fMomPos(0), fMomNeg(0),
  fTpcHitsPos(0), fTpcHitsNeg(0), fOnFlyStatusV0(0),
  fChi2V0(0),  fClV0(0),  fChi2Pos(0),  fClPos(0),  fChi2Neg(0),  fClNeg(0), fCosPointingAngle(0),
  fDedxPos(0),  fErrDedxPos(0),  fLenDedxPos(0),
  fDedxNeg(0),  fErrDedxNeg(0),  fLenDedxNeg(0),
  fNufDedxPos(0), fNufDedxNeg(0),
  fHelixPos(), fHelixNeg(),
  fMomV0(0), fEtaV0(0), fPhiV0(0), fYV0(0),
  fAlphaV0(0),  fPtArmV0(0),
  fELambda(0),  fEK0Short(0),
  fEPosProton(0),  fEPosPion(0),
  fENegProton(0),  fENegPion(0),
  fMassLambda(0),  fMassAntiLambda(0),
  fMassK0Short(0),  fRapLambda(0),
  fRapK0Short(0),  fCTauLambda(0),
  fCTauK0Short(0),  fPtV0(0),  fPtotV0(0),
  fPtPos(0),  fPtotPos(0),
  fPtNeg(0),  fPtotNeg(0),
  fEtaPos(0), fEtaNeg(0), fTPCNclsPos(0), fTPCNclsNeg(0), fClustersPos(0), fClustersNeg(0), fSharingPos(0), fSharingNeg(0), fNdofPos(0), fNdofNeg(0), fStatusPos(0), fStatusNeg(0),
  fPosNSigmaTPCK(0), fPosNSigmaTPCPi(0), fPosNSigmaTPCP(0), fNegNSigmaTPCK(0), fNegNSigmaTPCPi(0), fNegNSigmaTPCP(0),
  fPosNSigmaTOFK(0), fPosNSigmaTOFPi(0), fPosNSigmaTOFP(0), fNegNSigmaTOFK(0), fNegNSigmaTOFPi(0), fNegNSigmaTOFP(0),
  fKeyNeg(0),   fKeyPos(0),
  fNominalTpcEntrancePointPos(0,0,0),fNominalTpcExitPointPos(0,0,0),
  fNominalTpcEntrancePointNeg(0,0,0),fNominalTpcExitPointNeg(0,0,0),
  fNominalTpcPointPosShifted(0,0,0),fNominalTpcPointNegShifted(0,0,0),
  fTPCMomentumPos(0), fTPCMomentumNeg(0),
  fTOFProtonTimePos(0), fTOFPionTimePos(0), fTOFKaonTimePos(0),
  fTOFProtonTimeNeg(0), fTOFPionTimeNeg(0), fTOFKaonTimeNeg(0),
  fImpactDprimPos(-999), fImpactDweakPos(-999), fImpactDmatPos(-999), fImpactDprimNeg(-999), fImpactDweakNeg(-999), fImpactDmatNeg(-999),
  fCorrLam(0.0),
  fCorrLamMinus(0.0),
  fRadiusV0(0.0),
  fMultiplicity(0),
  fZvtx(0),
  fHiddenInfo(NULL)
{
  // Default empty constructor
  fTrackTopologyMapPos[0] = 0;
  fTrackTopologyMapPos[1] = 0;
  fTrackTopologyMapNeg[0] = 0;
  fTrackTopologyMapNeg[1] = 0;

  for(int i=0;i<9;i++)
    {
      fNominalTpcPointsPos[i].SetX(0);
      fNominalTpcPointsPos[i].SetY(0);
      fNominalTpcPointsPos[i].SetZ(0);
      fNominalTpcPointsNeg[i].SetX(0);
      fNominalTpcPointsNeg[i].SetY(0);
      fNominalTpcPointsNeg[i].SetZ(0);

    }
}
// -----------------------------------------------------------------------
AliFemtoV0::AliFemtoV0(const AliFemtoV0& v) :
  fDecayLengthV0(v.fDecayLengthV0), fDecayVertexV0(v.fDecayVertexV0), fPrimaryVertex(v.fPrimaryVertex),
  fDcaV0Daughters(v.fDcaV0Daughters), fDcaV0ToPrimVertex(v.fDcaV0ToPrimVertex),
  fDcaPosToPrimVertex(v.fDcaPosToPrimVertex), fDcaNegToPrimVertex(v.fDcaNegToPrimVertex),
  fMomPos(v.fMomPos), fMomNeg(v.fMomNeg),
  fTpcHitsPos(v.fTpcHitsPos), fTpcHitsNeg(v.fTpcHitsNeg), fOnFlyStatusV0(v.fOnFlyStatusV0),
  fChi2V0(v.fChi2V0),  fClV0(v.fClV0),  fChi2Pos(v.fChi2Pos),  fClPos(v.fClPos),  fChi2Neg(v.fChi2Neg),  fClNeg(v.fClNeg),
  fCosPointingAngle(v.fCosPointingAngle),
  fDedxPos(v.fDedxPos),  fErrDedxPos(v.fErrDedxPos),  fLenDedxPos(v.fLenDedxPos),
  fDedxNeg(v.fDedxNeg),  fErrDedxNeg(v.fErrDedxNeg),  fLenDedxNeg(v.fLenDedxNeg),
  fNufDedxPos(v.fNufDedxPos), fNufDedxNeg(v.fNufDedxNeg),
  fHelixPos(v.fHelixPos), fHelixNeg(v.fHelixNeg),
  fMomV0(0), fEtaV0(v.fEtaV0), fPhiV0(v.fPhiV0), fYV0(v.fYV0),
  fAlphaV0(0),  fPtArmV0(0),
  fELambda(0),  fEK0Short(0),
  fEPosProton(0),  fEPosPion(0),
  fENegProton(0),  fENegPion(0),
  fMassLambda(0),  fMassAntiLambda(0),
  fMassK0Short(0),  fRapLambda(0),
  fRapK0Short(0),  fCTauLambda(0),
  fCTauK0Short(0),  fPtV0(0),  fPtotV0(0),
  fPtPos(0),  fPtotPos(0),
  fPtNeg(0),  fPtotNeg(0),
  fEtaPos(v.fEtaPos), fEtaNeg(v.fEtaNeg), fTPCNclsPos(v.fTPCNclsPos), fTPCNclsNeg(v.fTPCNclsNeg),
  fClustersPos(v.fClustersPos), fClustersNeg(v.fClustersNeg),
  fSharingPos(v.fSharingPos), fSharingNeg(v.fSharingNeg),
  fNdofPos(v.fNdofPos), fNdofNeg(v.fNdofNeg),
  fStatusPos(v.fStatusPos), fStatusNeg(v.fStatusNeg),
  fPosNSigmaTPCK(v.fPosNSigmaTPCK), fPosNSigmaTPCPi(v.fPosNSigmaTPCPi),
  fPosNSigmaTPCP(v.fPosNSigmaTPCP), fNegNSigmaTPCK(v.fNegNSigmaTPCK),
  fNegNSigmaTPCPi(v.fNegNSigmaTPCPi), fNegNSigmaTPCP(v.fNegNSigmaTPCP),
  fPosNSigmaTOFK(v.fPosNSigmaTOFK), fPosNSigmaTOFPi(v.fPosNSigmaTOFPi),
  fPosNSigmaTOFP(v.fPosNSigmaTOFP), fNegNSigmaTOFK(v.fNegNSigmaTOFK),
  fNegNSigmaTOFPi(v.fNegNSigmaTOFPi), fNegNSigmaTOFP(v.fNegNSigmaTOFP),
  fKeyNeg(v.fKeyNeg),   fKeyPos(v.fKeyPos),
  fNominalTpcEntrancePointPos(v.fNominalTpcEntrancePointPos),
  fNominalTpcExitPointPos(v.fNominalTpcExitPointPos),
  fNominalTpcEntrancePointNeg(v.fNominalTpcEntrancePointNeg),
  fNominalTpcExitPointNeg(v.fNominalTpcExitPointNeg),
  fNominalTpcPointPosShifted(v.fNominalTpcPointPosShifted),
  fNominalTpcPointNegShifted(v.fNominalTpcPointNegShifted),
  fTPCMomentumPos(v.fTPCMomentumPos), fTPCMomentumNeg(v.fTPCMomentumNeg),
  fTOFProtonTimePos(v.fTOFProtonTimePos), fTOFPionTimePos(v.fTOFPionTimePos), fTOFKaonTimePos(v.fTOFKaonTimePos),
  fTOFProtonTimeNeg(v.fTOFProtonTimeNeg), fTOFPionTimeNeg(v.fTOFPionTimeNeg), fTOFKaonTimeNeg(v.fTOFKaonTimeNeg),
  fImpactDprimPos(v.fImpactDprimPos), fImpactDweakPos(v.fImpactDweakPos),
  fImpactDmatPos(v.fImpactDmatPos), fImpactDprimNeg(v.fImpactDprimNeg),
  fImpactDweakNeg(v.fImpactDweakNeg), fImpactDmatNeg(v.fImpactDmatNeg),
  fCorrLam(v.fCorrLam), fCorrLamMinus(v.fCorrLamMinus), fRadiusV0(v.fRadiusV0),
  fMultiplicity(v.fMultiplicity),
  fZvtx(v.fZvtx),
  fHiddenInfo( v.fHiddenInfo ? v.fHiddenInfo->Clone() : NULL)  /***/
{
  // copy constructor
  fTrackTopologyMapPos[0] = v.fTrackTopologyMapPos[0];
  fTrackTopologyMapPos[1] = v.fTrackTopologyMapPos[1];
  fTrackTopologyMapNeg[0] = v.fTrackTopologyMapNeg[0];
  fTrackTopologyMapNeg[1] = v.fTrackTopologyMapNeg[1];

  for (int i = 0; i < 9; i++) {
    fNominalTpcPointsPos[i] = v.fNominalTpcPointsPos[i];
    fNominalTpcPointsNeg[i] = v.fNominalTpcPointsNeg[i];
  }

  UpdateV0();
}
AliFemtoV0& AliFemtoV0::operator=(const AliFemtoV0& aV0)
{
  // assignment operator
  if (this == &aV0)
    return *this;
  fDecayLengthV0 = aV0.fDecayLengthV0;
  fDecayVertexV0 = aV0.fDecayVertexV0;
  fDcaV0Daughters = aV0.fDcaV0Daughters;
  fDcaV0ToPrimVertex = aV0.fDcaV0ToPrimVertex;
  fDcaPosToPrimVertex = aV0.fDcaPosToPrimVertex;
  fDcaNegToPrimVertex = aV0.fDcaNegToPrimVertex;
  fMomPos = aV0.fMomPos;
  fMomNeg = aV0.fMomNeg;

  fTrackTopologyMapPos[0] = aV0.fTrackTopologyMapPos[0];
  fTrackTopologyMapPos[1] = aV0.fTrackTopologyMapPos[1];
  fTrackTopologyMapNeg[0] = aV0.fTrackTopologyMapNeg[0];
  fTrackTopologyMapNeg[1] = aV0.fTrackTopologyMapNeg[1];

  fKeyPos = aV0.fKeyPos;
  fKeyNeg = aV0.fKeyNeg;

  fEtaPos = aV0.fEtaPos;
  fEtaNeg = aV0.fEtaNeg;
  fTPCNclsPos = aV0.fTPCNclsPos;
  fTPCNclsNeg = aV0.fTPCNclsNeg;
  fClustersPos = aV0.fClustersPos;
  fClustersNeg = aV0.fClustersNeg;
  fSharingPos = aV0.fSharingPos;
  fSharingNeg = aV0.fSharingNeg;
  fNdofPos = aV0.fNdofPos;
  fNdofNeg = aV0.fNdofNeg;
  fStatusPos = aV0.fStatusPos;
  fStatusNeg = aV0.fStatusNeg;
  fOnFlyStatusV0 = aV0.fOnFlyStatusV0;

  fPosNSigmaTPCK =  aV0.fPosNSigmaTPCK;
  fPosNSigmaTPCPi = aV0.fPosNSigmaTPCPi;
  fPosNSigmaTPCP = aV0.fPosNSigmaTPCP;
  fNegNSigmaTPCK = aV0.fNegNSigmaTPCK;
  fNegNSigmaTPCPi = aV0.fNegNSigmaTPCPi;
  fNegNSigmaTPCP = aV0.fNegNSigmaTPCP;
  fPosNSigmaTOFK = aV0.fPosNSigmaTOFK;
  fPosNSigmaTOFPi = aV0.fPosNSigmaTOFPi;
  fPosNSigmaTOFP = aV0.fPosNSigmaTOFP;
  fNegNSigmaTOFK = aV0.fNegNSigmaTOFK;
  fNegNSigmaTOFPi = aV0.fNegNSigmaTOFPi;
  fNegNSigmaTOFP = aV0.fNegNSigmaTOFP;

  fEtaV0 = aV0.fEtaV0;
  fPhiV0 = aV0.fPhiV0;
  fYV0 = aV0.fYV0;
  fCosPointingAngle = aV0.fCosPointingAngle;

  fTpcHitsPos = aV0.fTpcHitsPos;
  fTpcHitsNeg = aV0.fTpcHitsNeg;

  fChi2V0 = aV0.fChi2V0;
  fClV0 = aV0.fClV0;
  fChi2Pos = aV0.fChi2Pos;
  fClPos = aV0.fClPos;
  fChi2Neg = aV0.fChi2Neg;
  fClNeg = aV0.fClNeg;
  fDedxPos = aV0.fDedxPos;
  fErrDedxPos = aV0.fErrDedxPos;//Gael 04Fev2002
  fLenDedxPos = aV0.fLenDedxPos;//Gael 04Fev2002
  fDedxNeg = aV0.fDedxNeg;
  fErrDedxNeg = aV0.fErrDedxNeg;//Gael 04Fev2002
  fLenDedxNeg = aV0.fLenDedxNeg;//Gael 04Fev2002

  fNufDedxPos = aV0.fNufDedxPos;
  fNufDedxNeg = aV0.fNufDedxNeg;

  fHelixPos = aV0.fHelixPos;// Gael 12 Sept
  fHelixNeg = aV0.fHelixNeg;// Gael 12 Sept

  fNominalTpcEntrancePointPos = aV0.fNominalTpcEntrancePointPos;
  fNominalTpcExitPointPos = aV0.fNominalTpcExitPointPos;
  fNominalTpcEntrancePointPos = aV0.fNominalTpcEntrancePointPos;
  fNominalTpcExitPointPos = aV0.fNominalTpcExitPointPos;
  fNominalTpcEntrancePointNeg = aV0.fNominalTpcEntrancePointNeg;
  fNominalTpcExitPointNeg = aV0.fNominalTpcExitPointNeg;
  fNominalTpcPointPosShifted = aV0.fNominalTpcPointPosShifted;
  fNominalTpcPointNegShifted = aV0.fNominalTpcPointNegShifted;

  fTPCMomentumPos = aV0.fTPCMomentumPos;
  fTPCMomentumNeg = aV0.fTPCMomentumNeg;

  fTOFProtonTimePos=aV0.fTOFProtonTimePos; fTOFPionTimePos=aV0.fTOFPionTimePos; fTOFKaonTimePos=aV0.fTOFKaonTimePos;
  fTOFProtonTimeNeg=aV0.fTOFProtonTimeNeg; fTOFPionTimeNeg=aV0.fTOFPionTimeNeg; fTOFKaonTimeNeg=aV0.fTOFKaonTimeNeg;

  for (int i = 0; i < 9; i++) {
    fNominalTpcPointsPos[i] = aV0.fNominalTpcPointsPos[i];
    fNominalTpcPointsNeg[i] = aV0.fNominalTpcPointsNeg[i];
  }

  fImpactDprimPos = aV0.fImpactDprimPos;
  fImpactDweakPos = aV0.fImpactDweakPos;
  fImpactDmatPos = aV0.fImpactDmatPos;
  fImpactDprimNeg = aV0.fImpactDprimNeg;
  fImpactDweakNeg = aV0.fImpactDweakNeg;
  fImpactDmatNeg = aV0.fImpactDmatNeg;

  fCorrLam =  aV0.fCorrLam;
  fCorrLamMinus =  aV0.fCorrLamMinus;

  fRadiusV0 = aV0.fRadiusV0;

  fMultiplicity = aV0.fMultiplicity;
  fZvtx = aV0.fZvtx;
  
  if (fHiddenInfo) delete fHiddenInfo;
  fHiddenInfo = aV0.fHiddenInfo? aV0.fHiddenInfo->Clone() : NULL;// GR 11 DEC 02
  UpdateV0();

  return *this;
}

// -----------------------------------------------------------------------
void AliFemtoV0::UpdateV0()
{
  // Calc. derived members of the v0 class

  fMomV0      = fMomPos + fMomNeg;
  fPtV0       = fMomV0.Perp();
  fPtotV0     = fMomV0.Mag();
  fPtPos      = fMomPos.Perp();
  fPtotPos    = fMomPos.Mag();
  fPtNeg      = fMomNeg.Perp();
  fPtotNeg    = fMomNeg.Mag();
  fELambda    = ::sqrt(fPtotV0*fPtotV0+kMLAMBDA*kMLAMBDA);
  fEK0Short   = ::sqrt(fPtotV0*fPtotV0+kMKAON0SHORT*kMKAON0SHORT);
  fEPosProton = ::sqrt(fPtotPos*fPtotPos+kMPROTON*kMPROTON);
  fENegProton = ::sqrt(fPtotNeg*fPtotNeg+kMPROTON*kMPROTON);
  fEPosPion   = ::sqrt(fPtotPos*fPtotPos+kMPIONPLUS*kMPIONPLUS);
  fENegPion   = ::sqrt(fPtotNeg*fPtotNeg+kMPIONMINUS*kMPIONMINUS);

  Float_t tMomNegAlongV0 = fMomNeg*fMomV0 / ::sqrt(::pow(fPtotV0,2)),
          tMomPosAlongV0 = fMomPos*fMomV0 / ::sqrt(::pow(fPtotV0,2));

  if (tMomPosAlongV0 + tMomNegAlongV0 != 0) {
    fAlphaV0 = (tMomPosAlongV0-tMomNegAlongV0)/(tMomPosAlongV0+tMomNegAlongV0);
  }
   //printf("%1.15f\n",fPtotPos);
   //printf("%1.15f\n",tMomPosAlongV0);

  if (fPtotPos < tMomPosAlongV0) {
    fPtArmV0 = 0.0;
  } else {
    fPtArmV0 = ::sqrt(fPtotPos*fPtotPos - tMomPosAlongV0*tMomPosAlongV0);
  }

  fMassLambda = ::sqrt(::pow(fEPosProton+fENegPion,2)-::pow(fPtotV0,2));
  fMassAntiLambda = ::sqrt(::pow(fENegProton+fEPosPion,2)-::pow(fPtotV0,2));
  fMassK0Short = ::sqrt(::pow(fENegPion+fEPosPion,2)-::pow(fPtotV0,2));
  fRapLambda = 0.5*::log( (fELambda+fMomV0.z()) / (fELambda-fMomV0.z()) );

  fCTauLambda = kMLAMBDA*(fDecayLengthV0) / ::sqrt( ::pow((double)fMomV0.Mag(),2.) );

  fRapK0Short = 0.5*::log( (fEK0Short+fMomV0.z()) / (fEK0Short-fMomV0.z()) );
  fCTauK0Short = kMKAON0SHORT*(fDecayLengthV0) / ::sqrt( ::pow((double)fMomV0.Mag(),2.) );
}
// -----------------------------------------------------------------------
#ifndef __NO_STAR_DEPENDENCE_ALLOWED__
#ifdef __ROOT__
#include "StStrangeMuDstMaker/StV0MuDst.h"
AliFemtoV0::AliFemtoV0( StV0MuDst& v){ // from strangess micro dst structure
  fDecayLengthV0 = v.decayLengthV0();
  fDecayVertexV0 = AliFemtoThreeVector( v.decayVertexV0X(), v.decayVertexV0Y(), v.decayVertexV0Z() );
  fDcaV0Daughters = v.dcaV0Daughters();
  fDcaV0ToPrimVertex = v.dcaV0ToPrimVertex();
  fDcaPosToPrimVertex = v.dcaPosToPrimVertex();
  fDcaNegToPrimVertex = v.dcaNegToPrimVertex();
  fMomPos = AliFemtoThreeVector( v.momPosX(), v.momPosY(), v.momPosZ() );
  fMomNeg = AliFemtoThreeVector( v.momNegX(), v.momNegY(), v.momNegZ() );
#ifdef STHBTDEBUG
  cout << " hist pos ";
  cout << v.topologyMapPos().numberOfHits(kTpcId);
  cout << " hist neg ";
  cout << v.topologyMapNeg().numberOfHits(kTpcId) << endl;
#endif
  fTpcHitsPos = ( v.topologyMapPos().numberOfHits(kTpcId) );
  fTpcHitsNeg = ( v.topologyMapNeg().numberOfHits(kTpcId) );
  fTrackTopologyMapPos[0] = ( v.topologyMapPos().data(0) );
  fTrackTopologyMapPos[1] = ( v.topologyMapPos().data(1) );
  fTrackTopologyMapNeg[0] = ( v.topologyMapNeg().data(0) );
  fTrackTopologyMapNeg[1] = ( v.topologyMapNeg().data(1) );
  fKeyPos = v.keyPos();
  fKeyNeg = v.keyNeg();
  fChi2V0 = v.chi2V0();
  fClV0 = v.clV0();
  fChi2Pos = v.chi2Pos();
  fClPos = v.clPos();
  fChi2Neg = v.chi2Neg();
  fClNeg = v.clNeg();
  fDedxPos = v.dedxPos();
  fErrDedxPos = v.errDedxPos();//Gael 04Fev2002
  fLenDedxPos = v.lenDedxPos();//Gael 04Fev2002
  fDedxNeg = v.dedxNeg();
  fErrDedxNeg = v.errDedxNeg();//Gael 04Fev2002
  fLenDedxNeg = v.lenDedxNeg();//Gael 04Fev2002
  fNufDedxPos = v.nufDedxPos();
  fNufDedxNeg = v.nufDedxNeg();
  fHiddenInfo =  0;//GR 11 DEC 02

#ifdef STHBTDEBUG
  cout << " keyPos " << v.keyPos() << endl;
  cout << " keyNeg " << v.keyNeg() << endl;
#endif
  fMomV0 = AliFemtoThreeVector( v.momV0X(), v.momV0Y(), v.momV0Z() );
#ifdef STHBTDEBUG
  cout << " alpha  ";
  cout << v.alphaV0();
  cout << " ptArm  ";
  cout << v.ptArmV0() << endl;
#endif
  fAlphaV0 = v.alphaV0();
  fPtArmV0 = v.ptArmV0();
  fELambda = v.eLambda();
  fEK0Short = v.eK0Short();
  fEPosProton = v.ePosProton();
  fEPosPion = v.ePosPion();
  fENegPion = v.eNegPion();
  fENegProton = v.eNegProton();
  fMassLambda = v.massLambda();
  fMassAntiLambda = v.massAntiLambda();
  fMassK0Short = v.massK0Short();
  fRapLambda = v.rapLambda();
  fRapK0Short = v.rapK0Short();
  fCTauLambda = v.cTauLambda();
  fCTauK0Short = v.cTauK0Short();
  fPtV0 = v.ptV0();
  fPtotV0 = v.ptotV0();
  fPtPos = v.ptPos();
  fPtotPos = v.ptotPos();
  fDedxPos = v.dedxPos();
  fPtNeg = v.ptNeg();
  fPtotNeg = v.ptotNeg();
  fDedxNeg = v.dedxNeg();
}
#endif // __ROOT__
#endif  //  __NO_STAR_DEPENDENCE_ALLOWED__



void AliFemtoV0::SetHelixPos(const AliFmPhysicalHelixD& h){fHelixPos = h;}// Gael 12 Sept 02
const AliFmPhysicalHelixD& AliFemtoV0::HelixPos() const {return fHelixPos;}// Gael 12 Sept 02
void AliFemtoV0::SetHelixNeg(const AliFmPhysicalHelixD& h){fHelixNeg = h;}// Gael 12 Sept 02
const AliFmPhysicalHelixD& AliFemtoV0::HelixNeg() const {return fHelixNeg;}// Gael 12 Sept 02

void AliFemtoV0::SetCorrectionLambdas(const double& x){fCorrLam=x;}
float AliFemtoV0::CorrectionLambda() const {return fCorrLam;}
void AliFemtoV0::SetCorrectionLambdasMinus(const double& x){fCorrLamMinus=x;}
float AliFemtoV0::CorrectionLambdaMinus() const {return fCorrLamMinus;}

void AliFemtoV0::SetHiddenInfo(AliFemtoHiddenInfo* aHiddenInfo) {fHiddenInfo=aHiddenInfo;}
bool AliFemtoV0::ValidHiddenInfo() const { if (fHiddenInfo) return true; else return false; }
AliFemtoHiddenInfo* AliFemtoV0::GetHiddenInfo() const {return fHiddenInfo;}

AliFemtoThreeVector AliFemtoV0::NominalTpcPointPos(int i) const
{
  if (i < 0)
    return fNominalTpcPointsPos[0];
  if (i > 8)
    return fNominalTpcPointsPos[8];
  return fNominalTpcPointsPos[i];
}

AliFemtoThreeVector AliFemtoV0::NominalTpcPointNeg(int i) const
{
  if (i < 0)
    return fNominalTpcPointsNeg[0];
  if (i > 8)
    return fNominalTpcPointsNeg[8];
  return fNominalTpcPointsNeg[i];
}

//______________________

void AliFemtoV0::SetMultiplicity(int mult)
{
  fMultiplicity=mult;
}
//______________________
void AliFemtoV0::SetZvtx(double vtx)
{
  fZvtx=vtx;
}
//______________________
