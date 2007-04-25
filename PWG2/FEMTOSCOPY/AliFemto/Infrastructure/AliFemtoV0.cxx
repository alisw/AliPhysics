#include "AliFemtoV0.h"
#include "phys_constants.h"

// -----------------------------------------------------------------------
AliFemtoV0::AliFemtoV0(const AliFemtoV0& v){ // copy constructor
  fDecayLengthV0 = v.fDecayLengthV0;
  fDecayVertexV0 = v.fDecayVertexV0;
  fDcaV0Daughters = v.fDcaV0Daughters;
  fDcaV0ToPrimVertex = v.fDcaV0ToPrimVertex;
  fDcaPosToPrimVertex = v.fDcaPosToPrimVertex;
  fDcaNegToPrimVertex = v.fDcaNegToPrimVertex;
  fMomPos = v.fMomPos;
  fMomNeg = v.fMomNeg;

  fTrackTopologyMapPos[0] = v.fTrackTopologyMapPos[0];
  fTrackTopologyMapPos[1] = v.fTrackTopologyMapPos[1];
  fTrackTopologyMapNeg[0] = v.fTrackTopologyMapNeg[0];
  fTrackTopologyMapNeg[1] = v.fTrackTopologyMapNeg[1];
   
  fKeyPos = v.fKeyPos;
  fKeyNeg = v.fKeyNeg;
     
  fTpcHitsPos = v.fTpcHitsPos;
  fTpcHitsNeg = v.fTpcHitsNeg;

  fChi2V0 = v.fChi2V0;
  fClV0 = v.fClV0;
  fChi2Pos = v.fChi2Pos;
  fClPos = v.fClPos;
  fChi2Neg = v.fChi2Neg;
  fClNeg = v.fClNeg;
  fDedxPos = v.fDedxPos;
  fErrDedxPos = v.fErrDedxPos;//Gael 04Fev2002
  fLenDedxPos = v.fLenDedxPos;//Gael 04Fev2002
  fDedxNeg = v.fDedxNeg;
  fErrDedxNeg = v.fErrDedxNeg;//Gael 04Fev2002
  fLenDedxNeg = v.fLenDedxNeg;//Gael 04Fev2002

  fNufDedxPos = v.fNufDedxPos;
  fNufDedxNeg = v.fNufDedxNeg;

  fHelixPos = v.fHelixPos;// Gael 12 Sept
  fHelixNeg = v.fHelixNeg;// Gael 12 Sept
  fHiddenInfo = v.fHiddenInfo? v.fHiddenInfo->clone() : 0;// GR 11 DEC 02
  UpdateV0();
}
// -----------------------------------------------------------------------
void AliFemtoV0::UpdateV0(){
  //Calc. derived memebers of the v0 class
  float MomNegAlongV0, MomPosAlongV0;

   fMomV0  = fMomPos + fMomNeg;
   fPtV0   = fMomV0.perp();
   fPtotV0 = fMomV0.mag();
   fPtPos  = fMomPos.perp();
   fPtotPos= fMomPos.mag();
   fPtNeg  = fMomNeg.perp();
   fPtotNeg= fMomNeg.mag();
   fELambda= ::sqrt(fPtotV0*fPtotV0+M_LAMBDA*M_LAMBDA);
   fEK0Short= ::sqrt(fPtotV0*fPtotV0+M_KAON_0_SHORT*M_KAON_0_SHORT);
   fEPosProton = ::sqrt(fPtotPos*fPtotPos+M_PROTON*M_PROTON);
   fENegProton = ::sqrt(fPtotNeg*fPtotNeg+M_PROTON*M_PROTON);
   fEPosPion = ::sqrt(fPtotPos*fPtotPos+M_PION_PLUS*M_PION_PLUS);
   fENegPion = ::sqrt(fPtotNeg*fPtotNeg+M_PION_MINUS*M_PION_MINUS);
  
   MomNegAlongV0 =  fMomNeg*fMomV0 / ::sqrt(::pow(fPtotV0,2));
   MomPosAlongV0 =  fMomPos*fMomV0 / ::sqrt(::pow(fPtotV0,2));

   fAlphaV0 = (MomPosAlongV0-MomNegAlongV0)/(MomPosAlongV0+MomNegAlongV0);
   fPtArmV0 =  ::sqrt(fPtotPos*fPtotPos - MomPosAlongV0*MomPosAlongV0);
   fMassLambda = ::sqrt(::pow(fEPosProton+fENegPion,2)-::pow(fPtotV0,2));
   fMassAntiLambda = ::sqrt(::pow(fENegProton+fEPosPion,2)-::pow(fPtotV0,2));
   fMassK0Short = ::sqrt(::pow(fENegPion+fEPosPion,2)-::pow(fPtotV0,2));

   fRapLambda = 0.5*::log( (fELambda+fMomV0.z()) / (fELambda-fMomV0.z()) );
   fCTauLambda = M_LAMBDA*(fDecayLengthV0) / ::sqrt( ::pow((double)fMomV0.mag(),2.) );
   
   fRapK0Short = 0.5*::log( (fEK0Short+fMomV0.z()) / (fEK0Short-fMomV0.z()) );
   fCTauK0Short = M_KAON_0_SHORT*(fDecayLengthV0) / ::sqrt( ::pow((double)fMomV0.mag(),2.) );

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

void AliFemtoV0::SetHiddenInfo(AliFemtoHiddenInfo* aHiddenInfo) {fHiddenInfo=aHiddenInfo;}
bool AliFemtoV0::ValidHiddenInfo() const { if (fHiddenInfo) return true; else return false; }
AliFemtoHiddenInfo* AliFemtoV0::getHiddenInfo() const {return fHiddenInfo;}

