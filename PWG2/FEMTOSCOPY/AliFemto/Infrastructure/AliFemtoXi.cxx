#include "Infrastructure/AliFemtoXi.h"
#include "phys_constants.h"

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

