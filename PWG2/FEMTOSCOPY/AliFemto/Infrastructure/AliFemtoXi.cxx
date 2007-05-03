#include "Infrastructure/AliFemtoXi.h"
#include "phys_constants.h"

// -----------------------------------------------------------------------
AliFemtoXi::AliFemtoXi():
  fCharge(0), fDecayLengthXi(0),
  fDecayVertexXi(0),
  fDcaXiDaughters(0), fDcaXiToPrimVertex(0), fDcaBachelorToPrimVertex(0),
  fMomBachelor(0), fKeyBachelor(0),
  fTpcHitsBac(0), fChi2Xi(0), fClXi(0), fChi2Bachelor(0), fClBachelor(0),
  fDedxBachelor(0), fNufDedxBachelor(0), fMomXi(0),
  fAlphaXi(0), fPtArmXi(0),
  fEXi(0), fEOmega(0), fEBacPion(0), fEBacKaon(0),
  fMassXi(0), fMassOmega(0), fRapXi(0), fRapOmega(0),
  fCTauXi(0), fCTauOmega(0),
  fPtXi(0), fPtotXi(0), fPtBac(0), fPtotBac(0),
  fKeyBac(0)
{/* no-op */}
// -----------------------------------------------------------------------
void AliFemtoXi::UpdateXi(){
  //Calc. derived members of the xi class
  float MomV0AlongXi, MomBacAlongXi;

   fMomXi  = momV0() + momBac(); 
   fPtXi   = fMomXi.perp();
   fPtotXi = fMomXi.mag();
   fPtBac  = momBac().perp();
   fPtotBac= momBac().mag();
   fEXi= ::sqrt(fPtotXi*fPtotXi+M_XI_MINUS*M_XI_MINUS);
   fEOmega= ::sqrt(fPtotXi*fPtotXi+M_OMEGA_MINUS*M_OMEGA_MINUS);
   fEBacPion = ::sqrt(ptotBac()*ptotBac()+M_PION_MINUS*M_PION_MINUS);
   fEBacKaon = ::sqrt(ptotBac()*ptotBac()+M_KAON_MINUS*M_KAON_MINUS);

   MomV0AlongXi  =  momV0()*fMomXi / ::sqrt(::pow(fPtotXi,2));
   MomBacAlongXi =  momBac()*fMomXi / ::sqrt(::pow(fPtotXi,2));

   fAlphaXi = (MomBacAlongXi-MomV0AlongXi)/(MomBacAlongXi+MomV0AlongXi);
   fPtArmXi =  ::sqrt(ptotBac()*ptotBac() - MomBacAlongXi*MomBacAlongXi);
   fMassXi = ::sqrt(::pow(eBacPion()+eLambda(),2)-::pow(fPtotXi,2));
   fMassOmega = ::sqrt(::pow(eBacKaon()+eLambda(),2)-::pow(fPtotXi,2));

   fRapXi = 0.5*::log( (eXi()+fMomXi.z()) / (eXi()-fMomXi.z()) );
   fCTauXi = M_XI_MINUS*(fDecayLengthXi) / ::sqrt( ::pow((double)fMomXi.mag(),2.) );
   
   fRapOmega = 0.5*::log( (eOmega()+fMomXi.z()) / (eOmega()-fMomXi.z()) );// eO,
   fCTauOmega = M_OMEGA_MINUS*(fDecayLengthXi) / ::sqrt( ::pow((double)fMomXi.mag(),2.) );
}
// -----------------------------------------------------------------------
#ifdef __ROOT__
#ifndef __NO_STAR_DEPENDENCE_ALLOWED__
#include "StStrangeMuDstMaker/StXiMuDst.h"
AliFemtoXi::AliFemtoXi( StXiMuDst& xiFromMuDst)  : AliFemtoV0(xiFromMuDst) { // from strangess micro dst structure
  UpdateV0(); // the v0 stuff


  fCharge = xiFromMuDst.charge();
  fDecayLengthXi = xiFromMuDst.decayLengthXi(); // 12/07/2001 Gael
  fDecayVertexXi.setX(xiFromMuDst.decayVertexXiX());
  fDecayVertexXi.setY(xiFromMuDst.decayVertexXiY());
  fDecayVertexXi.setZ(xiFromMuDst.decayVertexXiZ());
  fDcaXiDaughters = xiFromMuDst.dcaXiDaughters();
  fDcaBachelorToPrimVertex = xiFromMuDst.dcaBachelorToPrimVertex();
  fDcaXiToPrimVertex = xiFromMuDst.dcaXiToPrimVertex();
  fMomBachelor.setX(xiFromMuDst.momBachelorX());
  fMomBachelor.setY(xiFromMuDst.momBachelorY());
  fMomBachelor.setZ(xiFromMuDst.momBachelorZ());
  
  fKeyBachelor = xiFromMuDst.keyBachelor();
  fTopologyMapBachelor[0] = xiFromMuDst.topologyMapBachelor().data(1);
  fTopologyMapBachelor[1] = xiFromMuDst.topologyMapBachelor().data(2);
  fTpcHitsBac = xiFromMuDst.topologyMapBachelor().numberOfHits(kTpcId); // 12/07/2001 Gael

  fChi2Xi = xiFromMuDst.chi2Xi();//nulle
  fClXi = xiFromMuDst.clXi();//nulle
  fChi2Bachelor = xiFromMuDst.chi2Bachelor();
  fClBachelor = xiFromMuDst.clBachelor();
  
  fDedxBachelor = xiFromMuDst.dedxBachelor();
  fNufDedxBachelor = xiFromMuDst.nufDedxBachelor();

  UpdateXi(); // the xi stuff
  
}

#endif // __NO_STAR_DEPENDENCE_ALLOWED__
#endif // __ROOT__

