#include "AliJetShape.h"

#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVector3.h"
#include "TVector2.h"
#include "AliFJWrapper.h"
using namespace std;

#ifdef FASTJET_VERSION

//________________________________________________________________________
Double32_t AliJetShapeGRNum::result(const fastjet::PseudoJet &jet) const {
  if (!jet.has_constituents())
    return 0; //AliFatal("Angular structure can only be applied on jets for which the constituents are known.");

  Double_t A = 0.;
  std::vector<fastjet::PseudoJet> constits = jet.constituents();
  for(UInt_t ic = 0; ic < constits.size(); ++ic) {
    /* Int_t uid = constits[ic].user_index(); */
    /* if (uid == -1) //skip ghost particle */
    /*   continue; */
    for(UInt_t jc = ic+1; jc < constits.size(); ++jc) {
      /* Int_t uid = constits[jc].user_index(); */
      /* if (uid == -1) //skip ghost particle */
      /* 	continue; */
      Double_t dphi = constits[ic].phi()-constits[jc].phi();
      if(dphi<-1.*TMath::Pi()) dphi+=TMath::TwoPi();
      if(dphi>TMath::Pi()) dphi-=TMath::TwoPi();
      Double_t dr2 = (constits[ic].eta()-constits[jc].eta())*(constits[ic].eta()-constits[jc].eta()) + dphi*dphi;
      if(dr2>0.) {
	Double_t dr = TMath::Sqrt(dr2);
	Double_t x = fR-dr;
	//noisy function
	Double_t noise = TMath::Exp(-x*x/(2*fDRStep*fDRStep))/(TMath::Sqrt(2.*TMath::Pi())*fDRStep);
	A += constits[ic].perp()*constits[jc].perp()*dr2*noise;
      }
    }
  }
  return A;
}

Double32_t AliJetShapeGRDen::result(const fastjet::PseudoJet &jet) const {
  if (!jet.has_constituents())
    return 0; //AliFatal("Angular structure can only be applied on jets for which the constituents are known.");

  Double_t A = 0.;
  std::vector<fastjet::PseudoJet> constits = jet.constituents();
  for(UInt_t ic = 0; ic < constits.size(); ++ic) {
    /* Int_t uid = constits[ic].user_index(); */
    /* if (uid == -1) //skip ghost particle */
    /*   continue; */
    for(UInt_t jc = ic+1; jc < constits.size(); ++jc) {
      /* Int_t uid = constits[jc].user_index(); */
      /* if (uid == -1) //skip ghost particle */
      /* 	continue; */
      Double_t dphi = constits[ic].phi()-constits[jc].phi();
      if(dphi<-1.*TMath::Pi()) dphi+=TMath::TwoPi();
      if(dphi>TMath::Pi()) dphi-=TMath::TwoPi();
      Double_t dr2 = (constits[ic].eta()-constits[jc].eta())*(constits[ic].eta()-constits[jc].eta()) + dphi*dphi;
      if(dr2>0.) {
	Double_t dr = TMath::Sqrt(dr2);
	Double_t x = fR-dr;
	//error function
	Double_t erf = 0.5*(1.+TMath::Erf(x/(TMath::Sqrt(2.)*fDRStep)));
	A += constits[ic].perp()*constits[jc].perp()*dr2*erf;
      }
    }
  }
  return A;
}

Double32_t AliJetShapeAngularity::result(const fastjet::PseudoJet &jet) const {
  if (!jet.has_constituents())
    return 0; 
  Double_t den=0.;
  Double_t num = 0.;
  std::vector<fastjet::PseudoJet> constits = jet.constituents();
  for(UInt_t ic = 0; ic < constits.size(); ++ic) {
    Double_t dphi = constits[ic].phi()-jet.phi();
    if(dphi<-1.*TMath::Pi()) dphi+=TMath::TwoPi();
    if(dphi>TMath::Pi()) dphi-=TMath::TwoPi();
    Double_t dr2 = (constits[ic].eta()-jet.eta())*(constits[ic].eta()-jet.eta()) + dphi*dphi;
    Double_t dr = TMath::Sqrt(dr2);
    num=num+constits[ic].perp()*dr;
    den=den+constits[ic].perp();
  }
  return num/den;
}

Double32_t AliJetShapepTD::result(const fastjet::PseudoJet &jet) const {
  if (!jet.has_constituents())
    return 0; 
  Double_t den=0;
  Double_t num = 0.;
  std::vector<fastjet::PseudoJet> constits = jet.constituents();
  for(UInt_t ic = 0; ic < constits.size(); ++ic) {
    num=num+constits[ic].perp()*constits[ic].perp();
    den=den+constits[ic].perp();
  }
  return TMath::Sqrt(num)/den;
}

Double32_t AliJetShapeCircularity::result(const fastjet::PseudoJet &jet) const {
  if (!jet.has_constituents())
    return 0;
  Double_t mxx    = 0.;
  Double_t myy    = 0.;
  Double_t mxy    = 0.;
  int  nc     = 0;
  Double_t sump2  = 0.;
  Double_t pxjet=jet.px();
  Double_t pyjet=jet.py();
  Double_t pzjet=jet.pz();
         
  //2 general normalized vectors perpendicular to the jet
  TVector3  ppJ1(pxjet, pyjet, pzjet);
  TVector3  ppJ3(- pxjet* pzjet, - pyjet * pzjet, pxjet * pxjet + pyjet * pyjet);
  ppJ3.SetMag(1.);
  TVector3  ppJ2(-pyjet, pxjet, 0);
  ppJ2.SetMag(1.);
    
  std::vector<fastjet::PseudoJet> constits = jet.constituents();
  for(UInt_t ic = 0; ic < constits.size(); ++ic) {
    TVector3 pp(constits[ic].px(), constits[ic].py(), constits[ic].pz());
    //local frame
    TVector3 pLong = pp.Dot(ppJ1) / ppJ1.Mag2() * ppJ1;
    TVector3 pPerp = pp - pLong;
    //projection onto the two perpendicular vectors defined above
    Float_t ppjX = pPerp.Dot(ppJ2);
    Float_t ppjY = pPerp.Dot(ppJ3);
    Float_t ppjT = TMath::Sqrt(ppjX * ppjX + ppjY * ppjY);
    if(ppjT<=0) return 0;
    mxx += (ppjX * ppjX / ppjT);
    myy += (ppjY * ppjY / ppjT);
    mxy += (ppjX * ppjY / ppjT);
    nc++;
    sump2 += ppjT;
  }
  if(nc<2) return 0;
  if(sump2==0) return 0;
  // Sphericity Matrix
  Double_t ele[4] = {mxx / sump2, mxy / sump2, mxy / sump2, myy / sump2};
  TMatrixDSym m0(2,ele);
      
  // Find eigenvectors
  TMatrixDSymEigen m(m0);
  TVectorD eval(2);
  TMatrixD evecm = m.GetEigenVectors();
  eval  = m.GetEigenValues();
  // Largest eigenvector
  int jev = 0;
  if (eval[0] < eval[1]) jev = 1;
  TVectorD evec0(2);
  // Principle axis
  evec0 = TMatrixDColumn(evecm, jev);
  Double_t compx=evec0[0];
  Double_t compy=evec0[1];
  TVector2 evec(compx, compy);
  Double_t circ=0;
  if(jev==1) circ=2*eval[0];
  if(jev==0) circ=2*eval[1];
    
  return circ;
}

Double32_t AliJetShapeSigma2::result(const fastjet::PseudoJet &jet) const {
  if (!jet.has_constituents())
    return 0;
  Double_t mxx    = 0.;
  Double_t myy    = 0.;
  Double_t mxy    = 0.;
  int  nc     = 0;
  Double_t sump2  = 0.;

  std::vector<fastjet::PseudoJet> constits = jet.constituents();
  for(UInt_t ic = 0; ic < constits.size(); ++ic) {
    Double_t ppt=constits[ic].perp();
    Double_t dphi = constits[ic].phi()-jet.phi();
    if(dphi<-1.*TMath::Pi()) dphi+=TMath::TwoPi();
    if(dphi>TMath::Pi()) dphi-=TMath::TwoPi();
    Double_t deta = constits[ic].eta()-jet.eta();
    mxx += ppt*ppt*deta*deta;
    myy += ppt*ppt*dphi*dphi;
    mxy -= ppt*ppt*deta*TMath::Abs(dphi);
    nc++;
    sump2 += ppt*ppt;
  }
  if(nc<2) return 0;
  if(sump2==0) return 0;
  // Sphericity Matrix
  Double_t ele[4] = {mxx , mxy, mxy, myy };
  TMatrixDSym m0(2,ele);
      
  // Find eigenvectors
  TMatrixDSymEigen m(m0);
  TVectorD eval(2);
  TMatrixD evecm = m.GetEigenVectors();
  eval  = m.GetEigenValues();
  // Largest eigenvector
  int jev = 0;
  if (eval[0] < eval[1]) jev = 1;
  TVectorD evec0(2);
  // Principle axis
  evec0 = TMatrixDColumn(evecm, jev);
  Double_t compx=evec0[0];
  Double_t compy=evec0[1];
  TVector2 evec(compx, compy);
  Double_t sigma2=0;
  if(jev==1) sigma2=TMath::Sqrt(TMath::Abs(eval[0])/sump2);
  if(jev==0) sigma2=TMath::Sqrt(TMath::Abs(eval[1])/sump2);
    
  return sigma2;
}


Double32_t AliJetShape1subjettiness_kt::result(const fastjet::PseudoJet &jet) const {
  if (!jet.has_constituents()) 
    return 0;
  AliFJWrapper *fFastjetWrapper = new AliFJWrapper("FJWrapper", "FJWrapper");
  Double32_t Result = fFastjetWrapper->AliFJWrapper::NSubjettinessDerivativeSub(1,0,0.2,1.0,0.4,jet,0);
  fFastjetWrapper->Clear();
  delete fFastjetWrapper;
  return Result;
}


Double32_t AliJetShape2subjettiness_kt::result(const fastjet::PseudoJet &jet) const {
  if (!jet.has_constituents()) 
    return 0;
  AliFJWrapper *fFastjetWrapper = new AliFJWrapper("FJWrapper", "FJWrapper");
  Double32_t Result = fFastjetWrapper->AliFJWrapper::NSubjettinessDerivativeSub(2,0,0.2,1.0,0.4,jet,0);
  fFastjetWrapper->Clear();
  delete fFastjetWrapper;
  return Result;
}

Double32_t AliJetShape3subjettiness_kt::result(const fastjet::PseudoJet &jet) const {
  if (!jet.has_constituents()) 
    return 0;
  AliFJWrapper *fFastjetWrapper = new AliFJWrapper("FJWrapper", "FJWrapper");
  Double32_t Result = fFastjetWrapper->AliFJWrapper::NSubjettinessDerivativeSub(3,0,0.2,1.0,0.4,jet,0);
  fFastjetWrapper->Clear();
  delete fFastjetWrapper;
  return Result;
}

Double32_t AliJetShapeOpeningAngle_kt::result(const fastjet::PseudoJet &jet) const {
  if (!jet.has_constituents()) 
    return 0;
  AliFJWrapper *fFastjetWrapper = new AliFJWrapper("FJWrapper", "FJWrapper");
  Double32_t Result = fFastjetWrapper->AliFJWrapper::NSubjettinessDerivativeSub(2,0,0.2,1.0,0.4,jet,1);
  fFastjetWrapper->Clear();
  delete fFastjetWrapper;
  return Result;
}

#endif

