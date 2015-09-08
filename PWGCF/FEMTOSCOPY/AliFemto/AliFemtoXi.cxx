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
  fKeyBac(0), fCosPointingAngleXi(0), fPhiXi(0),
  fEtaXi(0), fTPCNclsBac(0), fNdofBac(0), fStatusBac(0),
  fEtaBac(0), fIdBac(0), fBacNSigmaTPCK(-999),fBacNSigmaTPCPi(-999),
  fBacNSigmaTPCP(-999), fBacNSigmaTOFK(-999), fBacNSigmaTOFPi(-999),
  fBacNSigmaTOFP(-999), fTPCMomentumBac(0), fTOFProtonTimeBac(0), fTOFPionTimeBac(0),
  fTOFKaonTimeBac(0)
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


