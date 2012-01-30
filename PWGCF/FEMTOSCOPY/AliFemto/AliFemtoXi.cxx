///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoXi: special type of particle desling with the specifics       //
// of the Xi type of particle                                            //
// It stores the information both about the Xi itself and about it's     //
// daughters, as well as the bachelor particle from first decay vertex   //
// so that the caut betwen the daughter characteristics and the bachelor //
// is possible.                                                          //
//                                                                       //
///////////////////////////////////////////////////////////////////////////
#include "AliFemtoXi.h"
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
{
  fTopologyMapBachelor[0] = 0;
  fTopologyMapBachelor[1] = 0;
}
// -----------------------------------------------------------------------
void AliFemtoXi::UpdateXi(){
  //Calc. derived members of the xi class
  float tMomV0AlongXi, tMomBacAlongXi;

   fMomXi  = MomV0() + MomBac(); 
   fPtXi   = fMomXi.Perp();
   fPtotXi = fMomXi.Mag();
   fPtBac  = MomBac().Perp();
   fPtotBac= MomBac().Mag();
   fEXi= ::sqrt(fPtotXi*fPtotXi+kMXIMINUS*kMXIMINUS);
   fEOmega= ::sqrt(fPtotXi*fPtotXi+kMOMEGAMINUS*kMOMEGAMINUS);
   fEBacPion = ::sqrt(PtotBac()*PtotBac()+kMPIONMINUS*kMPIONMINUS);
   fEBacKaon = ::sqrt(PtotBac()*PtotBac()+kMKAONMINUS*kMKAONMINUS);

   tMomV0AlongXi  =  MomV0()*fMomXi / ::sqrt(::pow(fPtotXi,2));
   tMomBacAlongXi =  MomBac()*fMomXi / ::sqrt(::pow(fPtotXi,2));

   fAlphaXi = (tMomBacAlongXi-tMomV0AlongXi)/(tMomBacAlongXi+tMomV0AlongXi);
   fPtArmXi =  ::sqrt(PtotBac()*PtotBac() - tMomBacAlongXi*tMomBacAlongXi);
   fMassXi = ::sqrt(::pow(EBacPion()+ELambda(),2)-::pow(fPtotXi,2));
   fMassOmega = ::sqrt(::pow(EBacKaon()+ELambda(),2)-::pow(fPtotXi,2));

   fRapXi = 0.5*::log( (EXi()+fMomXi.z()) / (EXi()-fMomXi.z()) );
   fCTauXi = kMXIMINUS*(fDecayLengthXi) / ::sqrt( ::pow((double)fMomXi.Mag(),2.) );
   
   fRapOmega = 0.5*::log( (EOmega()+fMomXi.z()) / (EOmega()-fMomXi.z()) );// eO,
   fCTauOmega = kMOMEGAMINUS*(fDecayLengthXi) / ::sqrt( ::pow((double)fMomXi.Mag(),2.) );
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

