/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// Implementation for correlation analysis

#include "AliJCorrelations.h"
#include "AliJDataManager.h"

#include "AliJHistos.h"
#include "AliJBaseTrack.h"
#include "AliJCard.h"
#include "AliJHistos.h"
#include "AliJRunTable.h"



//==========================================================================
// don't count the trigger here ! You'll miss the not associated triggers
// but you  will miss them in the main code too, unles you loop over those
// fevent with no cgl (set the masMix to 1)
//==========================================================================
// blah

AliJCorrelations::AliJCorrelations( AliJCard *cardIn, AliJHistos *histosIn) :
  AliJCorrelationInterface(),
  fcard(cardIn),
  fhistos(histosIn),
  fnReal(0),
  fnMix(0),
  fsumTriggerAndAssoc(0),
  fsamplingMethod(0),    // flat by default
  fIsHeavyIon(0),
  fawayPhiGap(0),
  fmaxEtaRange(0),
  fRSignalBin(0),
  frandom(0x0),
  fptt(0),
  fpta(0),
  fTrackPairEfficiency(0),
  fIsIsolatedTrigger(false),
  fpttBin(0),
  fptaBin(0),
  fPhiTrigger(0),
  fPhiAssoc(0),
  fEtaTrigger(0),
  fEtaAssoc(0),
  fDeltaPhi(0),
  fDeltaPhiPiPi(0),
  fDeltaEta(0),
  fXlong(0),
  fNearSide(true),
  fNearSide3D(true),
  fEtaGapBin(0),
  fPhiGapBinNear(0),
  fPhiGapBinAway(0),
  fRGapBinNear(0),
  fRGapBinAway(0),
  fCentralityBin(0),
  fXlongBin(0),
  fIsLikeSign(false),
  fGeometricAcceptanceCorrection(1),
  fGeometricAcceptanceCorrection3D(1)
{
  // constructor
  
  fsumTriggerAndAssoc = int( fcard->Get("sumTriggerAndAssoc") );
  fmaxEtaRange    = fcard->Get("EtaRange");
  fRSignalBin = int(fcard->Get("EtaGapSignalBin"));
  
  fDPhiUERegion[0]    = fcard->Get("DPhiUERegion",0);
  fDPhiUERegion[1]    = fcard->Get("DPhiUERegion",1);
  cout << fmaxEtaRange <<" fDPhiUERegion[0]="<< fDPhiUERegion[0] <<" fDPhiUERegion[1]="<< fDPhiUERegion[1] <<endl;
  fIsHeavyIon = AliJRunTable::GetInstance().IsHeavyIon();
  fAcceptanceCorrection = new AliJAcceptanceCorrection(cardIn);

  // -----------------------------------------------------------------------------------------------
  // HARD CODED NUMBERS - VIOLATIONS - BREAKS the code when only on bin used in card.input!!!
  // Who did this? Further more it us unused variable !!!!
  // Jan
  //fawayPhiGap = fcard->GetBinBorder(kEtaGapType, 2);
  // -----------------------------------------------------------------------------------------------
  
  for(int iRGap=0; iRGap < fcard->GetNoOfBins(kRGapType); iRGap++)
    fRGap[iRGap]   = fcard->GetBinBorder( kRGapType, iRGap);
  
  //dPhiRange = fhistos->GetDphiRange();
  frandom = new TRandom3(); //FK// frandom generator for jt flow UE
  frandom->SetSeed(0); //FK//
  
}

AliJCorrelations::AliJCorrelations() :
  AliJCorrelationInterface(),
  fcard(0x0),
  fhistos(0x0),
  fAcceptanceCorrection(0x0),
  fnReal(0),
  fnMix(0),
  fsumTriggerAndAssoc(0),
  fsamplingMethod(0),    // flat by default
  fIsHeavyIon(0),
  fawayPhiGap(0),
  fmaxEtaRange(0),
  fRSignalBin(0),
  frandom(0x0),
  fptt(0),
  fpta(0),
  fTrackPairEfficiency(0),
  fIsIsolatedTrigger(false),
  fpttBin(0),
  fptaBin(0),
  fPhiTrigger(0),
  fPhiAssoc(0),
  fEtaTrigger(0),
  fEtaAssoc(0),
  fDeltaPhi(0),
  fDeltaPhiPiPi(0),
  fDeltaEta(0),
  fXlong(0),
  fNearSide(true),
  fNearSide3D(true),
  fEtaGapBin(0),
  fPhiGapBinNear(0),
  fPhiGapBinAway(0),
  fRGapBinNear(0),
  fRGapBinAway(0),
  fCentralityBin(0),
  fXlongBin(0),
  fIsLikeSign(false),
  fGeometricAcceptanceCorrection(1),
  fGeometricAcceptanceCorrection3D(1)
{
  // default constructor
}

AliJCorrelations::AliJCorrelations(const AliJCorrelations& in) :
  fcard(in.fcard),
  fhistos(in.fhistos),
  fAcceptanceCorrection(in.fAcceptanceCorrection),
  fnReal(in.fnReal),
  fnMix(in.fnMix),
  fsumTriggerAndAssoc(in.fsumTriggerAndAssoc),
  fsamplingMethod(in.fsamplingMethod),
  fIsHeavyIon(in.fIsHeavyIon),
  fawayPhiGap(in.fawayPhiGap),
  fmaxEtaRange(in.fmaxEtaRange),
  fRSignalBin(in.fRSignalBin),
  frandom(in.frandom),
  fptt(in.fptt),
  fpta(in.fpta),
  fTrackPairEfficiency(in.fTrackPairEfficiency),
  fIsIsolatedTrigger(in.fIsIsolatedTrigger),
  fpttBin(in.fpttBin),
  fptaBin(in.fptaBin),
  fPhiTrigger(in.fPhiTrigger),
  fPhiAssoc(in.fPhiAssoc),
  fEtaTrigger(in.fEtaTrigger),
  fEtaAssoc(in.fEtaAssoc),
  fDeltaPhi(in.fDeltaPhi),
  fDeltaPhiPiPi(in.fDeltaPhiPiPi),
  fDeltaEta(in.fDeltaEta),
  fXlong(in.fXlong),
  fNearSide(in.fNearSide),
  fNearSide3D(in.fNearSide3D),
  fEtaGapBin(in.fEtaGapBin),
  fPhiGapBinNear(in.fPhiGapBinNear),
  fPhiGapBinAway(in.fPhiGapBinAway),
  fRGapBinNear(in.fRGapBinNear),
  fRGapBinAway(in.fRGapBinAway),
  fCentralityBin(in.fCentralityBin),
  fXlongBin(in.fXlongBin),
  fIsLikeSign(in.fIsLikeSign),
  fGeometricAcceptanceCorrection(in.fGeometricAcceptanceCorrection),
  fGeometricAcceptanceCorrection3D(in.fGeometricAcceptanceCorrection3D)
{
  // The pointers to card and histos are just copied. I think this is safe, since they are not created by
  // AliJCorrelations and thus should not disappear if the AliJCorrelation managing them is destroyed.
  
  fDPhiUERegion[0] = in.fDPhiUERegion[0];
  fDPhiUERegion[1] = in.fDPhiUERegion[1];
  
  for(int iRGap=0; iRGap < fcard->GetNoOfBins(kRGapType); iRGap++){
    fRGap[iRGap] = in.fRGap[iRGap];
  }
  
  frandom = new TRandom3(); // frandom generator for jt flow UE
  frandom->SetSeed(0);
}

AliJCorrelations& AliJCorrelations::operator=(const AliJCorrelations& in){
  // Assingment operator
  
  if (&in==this) return *this;
  
  fptt = in.fptt;
  fpta = in.fpta;
  fTrackPairEfficiency = in.fTrackPairEfficiency;
  fIsIsolatedTrigger = in.fIsIsolatedTrigger;
  fpttBin = in.fpttBin;
  fptaBin = in.fptaBin;
  fPhiTrigger = in.fPhiTrigger;
  fPhiAssoc = in.fPhiAssoc;
  fEtaTrigger = in.fEtaTrigger;
  fEtaAssoc = in.fEtaAssoc;
  fDeltaPhi = in.fDeltaPhi;
  fDeltaPhiPiPi = in.fDeltaPhiPiPi;
  fDeltaEta = in.fDeltaEta;
  fXlong = in.fXlong;
  fNearSide = in.fNearSide;
  fNearSide3D = in.fNearSide3D;
  fEtaGapBin = in.fEtaGapBin;
  fPhiGapBinNear = in.fPhiGapBinNear;
  fPhiGapBinAway = in.fPhiGapBinAway;
  fRGapBinNear = in.fRGapBinNear;
  fRGapBinAway = in.fRGapBinAway;
  fCentralityBin = in.fCentralityBin;
  fXlongBin = in.fXlongBin;
  fIsLikeSign = in.fIsLikeSign;
  fGeometricAcceptanceCorrection = in.fGeometricAcceptanceCorrection;
  fGeometricAcceptanceCorrection3D = in.fGeometricAcceptanceCorrection3D;
  fnReal = in.fnReal;
  fnMix = in.fnMix;
  fsumTriggerAndAssoc = in.fsumTriggerAndAssoc;
  fsamplingMethod = in.fsamplingMethod;
  fIsHeavyIon = in.fIsHeavyIon;
  fawayPhiGap = in.fawayPhiGap;
  fmaxEtaRange = in.fmaxEtaRange;
  fRSignalBin = in.fRSignalBin;
  
  // The pointers to card and histos are just copied. I think this is safe, since they are not created by
  // AliJCorrelations and thus should not disappear if the AliJCorrelation managing them is destroyed.
  fcard = in.fcard;
  fhistos = in.fhistos;
  fAcceptanceCorrection = in.fAcceptanceCorrection;
  
  fDPhiUERegion[0] = in.fDPhiUERegion[0];
  fDPhiUERegion[1] = in.fDPhiUERegion[1];
  
  for(int iRGap=0; iRGap < fcard->GetNoOfBins(kRGapType); iRGap++){
    fRGap[iRGap] = in.fRGap[iRGap];
  }
  
  frandom = new TRandom3(); // frandom generator for jt flow UE
  frandom->SetSeed(0);
  
  return *this;
  // copy constructor
}


void AliJCorrelations::FillHisto(corrFillType cFTyp, fillType fTyp, int cBin, int zBin, AliJBaseTrack *ftk1, AliJBaseTrack *ftk2){
  // histo filler
  if( cFTyp == kAzimuthFill )
    FillAzimuthHistos( fTyp, cBin, zBin, ftk1, ftk2);
  
}

//=============================================================================================
void AliJCorrelations::FillAzimuthHistos(fillType fTyp, int CentBin, int ZBin, AliJBaseTrack *ftk1, AliJBaseTrack *ftk2)
//=============================================================================================
{
  // histo filler
  bool twoTracks = false;
  if(ftk1->GetParticleType()==kJHadron && ftk2->GetParticleType()==kJHadron) twoTracks =true;
  
  //double-counting check
  if(fTyp == kReal && twoTracks && ftk1->GetID()==ftk2->GetID()) return;
  
  // Check the signs of the paired particles
  fIsLikeSign = false;
  if(ftk1->GetCharge() > 0 && ftk2->GetCharge() > 0) fIsLikeSign = true;
  if(ftk1->GetCharge() < 0 && ftk2->GetCharge() < 0) fIsLikeSign = true;
  
  //----------------------------------------------------------------
  fptt = ftk1->Pt();
  fpta = ftk2->Pt();
  
  fTrackPairEfficiency = 1./( ftk1->GetTrackEff() * ftk2->GetTrackEff() );
  
  fIsIsolatedTrigger =  ftk1->GetIsIsolated()>0  ? true : false; //FK// trigger particle is isolated hadron
  
  TLorentzVector vTrigger = ftk1->GetLorentzVector(), vAssoc = ftk2->GetLorentzVector();    // Lorentz vectors for trigger and associated particles
  
  //phit= ftk1->Phi();                                 //for RP
  fpttBin       = ftk1->GetTriggBin();
  fptaBin       = ftk2->GetAssocBin();
  fPhiTrigger   = ftk1->Phi();
  fPhiAssoc     = ftk2->Phi();
  fDeltaPhi     = DeltaPhi(fPhiTrigger, fPhiAssoc);  //radians
  fDeltaPhiPiPi = atan2(sin(fPhiTrigger-fPhiAssoc), cos(fPhiTrigger-fPhiAssoc));
  fDeltaEta     = ftk1->Eta() - ftk2->Eta();
  fEtaTrigger   = ftk1->Eta();
  fEtaAssoc     = ftk2->Eta();
  //double dEtaFar   = ftk1->Eta() + ftk2->Eta();
  
  fNearSide     = cos(fPhiTrigger-fPhiAssoc) > 0 ? true : false;  // Traditional near side definition using deltaPhi
  fNearSide3D   = vTrigger.Vect().Dot(vAssoc.Vect()) > 0 ? true : false; // Near side definition using half ball around the trigger

  fEtaGapBin = fcard->GetBin( kEtaGapType, fabs(fDeltaEta));
  fPhiGapBinNear = fcard->GetBin( kEtaGapType, fabs(fDeltaPhiPiPi) );
  fPhiGapBinAway = fcard->GetBin( kEtaGapType, fabs(fDeltaPhi-kJPi) ); //here the angle must be 0-2pi and not (-pi,pi)
  fRGapBinNear   = fcard->GetBin( kRGapType, fabs(ftk1->DeltaR(*ftk2)));
  fRGapBinAway   = fcard->GetBin( kRGapType, sqrt(pow(fDeltaPhi-kJPi,2)+fDeltaEta*fDeltaEta) );
  fCentralityBin = CentBin;
  
  
  fXlong = vTrigger.Vect().Dot(vAssoc.Vect())/pow(vTrigger.P(),2);
  fXlongBin = fcard->GetBin(kXeType, TMath::Abs(fXlong));
  
  //if( rGapBin != fRGapBinNear ) cout<<"dR vs fRGapBinNear = "<<rGapBin<<"\t"<<fRGapBinNear<<endl;
  //if( rGapBin != fRGapBinAway ) cout<<"dR vs fRGapBinAway = "<<rGapBin<<"\t"<<fRGapBinAway<<endl;
  //----------------------------------------------------------------
  
  //acceptance correction  triangle  or mixed fevent
  //  fGeometricAcceptanceCorrection = 1;
  fGeometricAcceptanceCorrection = fAcceptanceCorrection->GetAcceptanceCorrectionTraditional(fsamplingMethod, fDeltaEta, fDeltaPhiPiPi, fCentralityBin, fpttBin);
  fGeometricAcceptanceCorrection3D = fAcceptanceCorrection->GetAcceptanceCorrection3DNearSide(fsamplingMethod, fDeltaEta, fDeltaPhiPiPi, fCentralityBin, fpttBin);
  
  if(fpttBin<0 || fptaBin<0 || fEtaGapBin<0 ){
    cout<<"Error in FillAzimuthHistos: some pT or eta out of bin. pttBin="<<fpttBin<<" pTaBin="<<fptaBin <<" etaGapBin="<< fEtaGapBin << endl;
    ftk1->Print();
    ftk2->Print();
    exit(-1);
    //return;
  }
  
  if(fDeltaPhi==0) cout <<" fdphi=0; fptt="<<  fptt<<"   fpta="<<fpta<<"  TID="<<ftk1->GetID()<<"  AID="<<ftk2->GetID() <<" tphi="<< fPhiTrigger <<" aphi="<< fPhiAssoc << endl;
  
  // ===================================================================
  // =====================  Fill Histograms  ===========================
  // ===================================================================
  
  bool fill2DBackgroundQualityControlHistograms = false;  // Choose whether to fill the DeltaPhi DeltaEta histograms for jT background
  
  // If quality control level is high enough, fill 2D quality control histograms
  if(fcard->Get("QualityControlLevel")>1) fill2DBackgroundQualityControlHistograms = true;
  
  //if(fhistos->fhCosThetaStar.Dimension()>0) FillPairPtAndCosThetaStarHistograms(fTyp, ftk1, ftk2);  // Fill the pair pT and cos(theta*) histograms TODO: Does not work! Needs debugging
  if(fhistos->fhxEF.Dimension()>0) FillXeHistograms(fTyp);  // Fill the xE and xLong histograms
  FillDeltaEtaHistograms(fTyp, ZBin);  // Fill all the delta eta histograms
  FillDeltaPhiHistograms(fTyp);  // Fill the azimuthal correlation functions
  FillDeltaEtaDeltaPhiHistograms(fTyp, ZBin);  // Fill the 2D correlation functions
  FillPtaHistograms(fTyp);  // Fill various pTa histograms
  if(fhistos->fhDRNearPtMoon.Dimension()>0) FillIAAAndMoonHistograms(fTyp, ZBin);  // Fill the I_AA and moon histograms
  
}


void AliJCorrelations::FillPairPtAndCosThetaStarHistograms(fillType fTyp, AliJBaseTrack *ftk1, AliJBaseTrack *ftk2)
{
  // This method fills the pair pT and Cos(thata*) histograms
  
  TVector3 v3trigg = ftk1->GetLorentzVector().Vect();
  TVector3 v3assoc = ftk2->GetLorentzVector().Vect();
  double pairMass = sqrt(2*v3trigg.Mag()*v3assoc.Mag()*(1-cos(v3trigg.Angle(v3assoc))));
  int imass = fcard->GetBin(kMassType, pairMass);
  TVector3 v3pairPt = v3trigg + v3assoc;
  double pairPt = v3pairPt.Perp();
  double wPairPt = pairPt>0 ? 1/pairPt : 0;
  
  //================================
  // Pair pT and cos(Theta^star) ===  // TODO: Needs debugging, the code below does not work!
  //================================
  
  if(!fNearSide) {
    fhistos->fhPairPt[fTyp][fpttBin][fptaBin]->Fill( pairPt, fTrackPairEfficiency);
    double pairDphi = v3pairPt.DeltaPhi(v3trigg);
    if(pairDphi<-kJPi/3.0) pairDphi += kJTwoPi;
    fhistos->fhPairPtDphi[fTyp][fpttBin][fptaBin]->Fill( pairDphi/kJPi );
  }
  
  if ( fTyp == kReal ) {
    //cout<<" fptt="<< fptt <<" fpta="<< fpta <<" fdphi="<< fDeltaPhi <<" pairPt="<< pairPt << endl;
    
    int ipairPt = fcard->IsLessThanUpperPairPtCut(pairPt);
    //cout<<"ppt="<<pairPt<<" mass="<<pairMass<<" ip="<<ipairPt<<" im="<<imass<<endl;
    if(imass>=0) {
      fhistos->fhPairPtMass[imass]->Fill( pairPt, fTrackPairEfficiency*wPairPt );
      fhistos->fhPairDPhi[imass]->Fill( fDeltaPhi/kJPi, fTrackPairEfficiency );
      fhistos->fhPairDpT[imass]->Fill( fptt-fpta, fTrackPairEfficiency );
    }
    if(imass>=0 && ipairPt>=0 ){
      double cmsRap = (ftk1->Eta()+ftk2->Eta())/2.;
      double cosThetaStar = fabs(cos(2*atan(exp(-ftk1->Eta()+cmsRap))));
      fhistos->fhCosThetaStar[fTyp][ipairPt][imass]->Fill( cosThetaStar, fTrackPairEfficiency );
      
      //frandom boost, Let's try just trigger phi as a cms boost.
      cosThetaStar = fabs(cos(2*atan(exp(-ftk1->Eta()+fPhiTrigger))));
      fhistos->fhCosThetaStar[kMixed][ipairPt][imass]->Fill( cosThetaStar, fTrackPairEfficiency );
      
      fhistos->fhInvMass[ipairPt]->Fill(pairMass, fTrackPairEfficiency);
      fhistos->fhCMSrap[ipairPt][imass]->Fill( cosThetaStar, cmsRap*fTrackPairEfficiency );
      fhistos->fpCMSrap->Fill( cosThetaStar, fabs(cmsRap)*fTrackPairEfficiency );
    }
  }
  
  //---------------------------------
  //if(pairMass>10 && fcard->IsLessThanUpperPairPtCut(pairPt) ){
  //    cout<<pairMass<<" e1="<<ftk1->Eta()-cmsRap<<" e2="<<ftk2->Eta()-cmsRap<<" cms="<<cmsRap;
  //    cout<<" ppt="<<pairPt<<" tphi="<<fPhiTrigger<<" aphi="<<fPhiAssoc<<" t1="<<fptt<<" t2="<<fpta<<endl;
  //    double thetaStar = fabs(2*atan(exp(-ftk1->Eta()-cmsRap)));
  //    cout<<thetaStar<<" "<<fabs(cos(thetaStar))<<endl;
  //}
  
  //---------------------------------
}

void AliJCorrelations::FillXeHistograms(fillType fTyp)
{
  // This method fills the xE and xLong histograms
  
  double xe = -fpta*cos(fPhiTrigger-fPhiAssoc)/fptt;
  
  if( fTyp == kReal ) {
    fhistos->fhxEPtBin[0][fpttBin][fptaBin]->Fill(fXlong, fGeometricAcceptanceCorrection * fTrackPairEfficiency);
    if( fNearSide ) {
      fhistos->fhxEPtBin[1][fpttBin][fptaBin]->Fill(fXlong, fGeometricAcceptanceCorrection * fTrackPairEfficiency);
    } else {
      fhistos->fhxEPtBin[2][fpttBin][fptaBin]->Fill(fXlong, fGeometricAcceptanceCorrection * fTrackPairEfficiency);
    }
  }
  
  if(fNearSide) {
    fhistos->fhxEN [fTyp][fpttBin]->Fill(-xe, fGeometricAcceptanceCorrection * fTrackPairEfficiency);
  } else {
    fhistos->fhxEF [fTyp][fpttBin]->Fill(xe, fGeometricAcceptanceCorrection * fTrackPairEfficiency);
    if(fIsIsolatedTrigger) fhistos->fhxEFIsolTrigg[fTyp][fpttBin]->Fill(xe, fGeometricAcceptanceCorrection * fTrackPairEfficiency);
  }
}

void AliJCorrelations::FillDeltaEtaHistograms(fillType fTyp, int ZBin)
{
  // This method fills the DeltaEta histograms
  
  if( fNearSide ){ //one could check the phiGapBin, but in the pi/2 <1.6 and thus phiGap is always>-1
    if( fTyp == 0 ) {
      fhistos->fhDEtaNear[fCentralityBin][ZBin][fPhiGapBinNear][fpttBin][fptaBin]->Fill( fDeltaEta , fGeometricAcceptanceCorrection * fTrackPairEfficiency );
    } else {
      fhistos->fhDEtaNearM[fCentralityBin][ZBin][fPhiGapBinNear][fpttBin][fptaBin]->Fill( fDeltaEta , fGeometricAcceptanceCorrection * fTrackPairEfficiency );
      fhistos->fhDetaNearMixAcceptance[fCentralityBin][fpttBin][fptaBin]->Fill( fDeltaEta, fTrackPairEfficiency);
    }
  } else {
    if(fPhiGapBinAway<=3) fhistos->fhDEtaFar[fTyp][fCentralityBin][fpttBin]->Fill( fDeltaEta, fGeometricAcceptanceCorrection * fTrackPairEfficiency );
  }
  
  // Different near side definition for xlong bins
  if( fNearSide3D ){
    if( fTyp == 0 ) {
      if(fPhiGapBinNear>=0 && fXlongBin >= 0) fhistos->fhDEtaNearXEbin[fCentralityBin][ZBin][fPhiGapBinNear][fpttBin][fXlongBin]->Fill( fDeltaEta , fGeometricAcceptanceCorrection3D * fTrackPairEfficiency );
    } else {
      if(fPhiGapBinNear>=0 && fXlongBin >= 0){
        fhistos->fhDEtaNearMXEbin[fCentralityBin][ZBin][fPhiGapBinNear][fpttBin][fXlongBin]->Fill( fDeltaEta , fGeometricAcceptanceCorrection3D * fTrackPairEfficiency );
        fhistos->fhDeta3DNearMixAcceptance[fCentralityBin][fpttBin][fXlongBin]->Fill( fDeltaEta, fTrackPairEfficiency);
      }
    }
  }
}

void AliJCorrelations::FillDeltaPhiHistograms(fillType fTyp)
{
  // This method fills the DeltaPhi histograms
  
  // When hists are filled for thresholds they are not properly normalized and need to be subtracted
  // This induced improper errors - subtraction of not-independent entries
  
  fhistos->fhDphiAssoc[fTyp][fCentralityBin][fEtaGapBin][fpttBin][fptaBin]->Fill( fDeltaPhi/kJPi , fGeometricAcceptanceCorrection * fTrackPairEfficiency);
  if(fXlongBin>=0 && fNearSide3D) fhistos->fhDphiAssocXEbin[fTyp][fCentralityBin][fEtaGapBin][fpttBin][fXlongBin]->Fill( fDeltaPhi/kJPi , fGeometricAcceptanceCorrection3D * fTrackPairEfficiency);
  
  if(fIsIsolatedTrigger) fhistos->fhDphiAssocIsolTrigg[fTyp][fCentralityBin][fpttBin][fptaBin]->Fill( fDeltaPhi/kJPi , fGeometricAcceptanceCorrection * fTrackPairEfficiency); //FK//
}

void AliJCorrelations::FillDeltaEtaDeltaPhiHistograms(fillType fTyp, int zBin)
{
  // This method fills the two dimensional DeltaEta,DeltaPhi histograms
  // No acceptance correction here, since we want to see tha structure caused by acceptance effects
  
  // Fill the histogram in pTa bins
  if(fNearSide){
    fhistos->fhDphiDetaPta[fTyp][fCentralityBin][zBin][fpttBin][fptaBin]->Fill(fDeltaEta, fDeltaPhiPiPi, fTrackPairEfficiency);
  }
  
  // Fill the histogram in xlong bins
  if(fNearSide3D && fXlongBin >= 0){
    fhistos->fhDphiDetaXlong[fTyp][fCentralityBin][zBin][fpttBin][fXlongBin]->Fill(fDeltaEta, fDeltaPhiPiPi, fTrackPairEfficiency);
  }
  
}

void AliJCorrelations::FillPtaHistograms(fillType fTyp)
{
  // This method fills various pta histograms
  
  if ( fTyp == kReal ) {
    //must be here, not in main, to avoid counting triggers
    fhistos->fhAssocPtBin[fCentralityBin][fpttBin][fptaBin]->Fill(fpta ); //I think It should not be weighted by Eff
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++
    // in order to get mean pTa in the jet peak one has
    // to fill fhMeanPtAssoc in |DeltaEta|<0.4
    // +++++++++++++++++++++++++++++++++++++++++++++++++
    if(fEtaGapBin>=0 && fEtaGapBin<2){
      fhistos->fhMeanPtAssoc[fCentralityBin][fpttBin][fptaBin]->Fill( fDeltaPhi/kJPi , fpta );
      fhistos->fhMeanZtAssoc[fCentralityBin][fpttBin][fptaBin]->Fill( fDeltaPhi/kJPi , fpta/fptt);
    }
    
    //UE distribution
    if(fabs(fDeltaPhiPiPi/kJPi)>fDPhiUERegion[0] && fabs(fDeltaPhiPiPi/kJPi)<fDPhiUERegion[1]){
      for(int iEtaGap=0; iEtaGap<=fEtaGapBin; iEtaGap++)  //FK// UE Pta spectrum for different eta gaps
        fhistos->fhPtAssocUE[fCentralityBin][iEtaGap][fpttBin]->Fill(fpta, fTrackPairEfficiency);
      if(fIsIsolatedTrigger){ //FK// trigger is isolated hadron
        fhistos->fhPtAssocUEIsolTrigg[fpttBin]->Fill(fpta, fTrackPairEfficiency); //FK//
      }
    }
    if(fabs(fDeltaPhi/kJPi)<0.15) fhistos->fhPtAssocN[fpttBin]->Fill(fpta, fTrackPairEfficiency);
    if(fabs(fDeltaPhi/kJPi-1)<0.15) fhistos->fhPtAssocF[fpttBin]->Fill(fpta, fTrackPairEfficiency);
    
    fnReal++;
  } else { // only mix
    fnMix++;
  }
}

void AliJCorrelations::FillIAAAndMoonHistograms(fillType fTyp, int ZBin)
{
  // This method fills the I_AA and moon histograms
  
  if(fhistos->Is2DHistosEnabled()) fhistos->fhDphiAssoc2DIAA[fTyp][fCentralityBin][ZBin][fpttBin][fptaBin]->Fill( fDeltaEta, fDeltaPhi/kJPi, fTrackPairEfficiency);
  
  if(fRGapBinNear>=0){
    if(fRGapBinNear <= fRSignalBin) fhistos->fhDRNearPt[fTyp][fCentralityBin][ZBin][fRGapBinNear][fpttBin]->Fill( fpta, fGeometricAcceptanceCorrection * fTrackPairEfficiency );
    // - moon -
    if(fRGapBinNear>0){
      for( int irs = 0; irs <= fRSignalBin;irs++ ){
        if( fRGapBinNear < irs ) continue;
        if( fPhiGapBinNear > irs ) continue;
        for ( int ir1=0;ir1<= fRGapBinNear;ir1++ ){
          if( irs > ir1 ) continue;
          if( ir1 > fRGapBinNear ) continue;
          double rGapS= fRGap[irs+1];
          double rGap1= fRGap[ir1+1];
          if( rGap1 < TMath::Abs(fDeltaPhiPiPi) ) continue;
          if( rGapS < TMath::Abs(fDeltaPhiPiPi) ) continue;
          double dEtaMin = sqrt(rGap1*rGap1-fDeltaPhiPiPi*fDeltaPhiPiPi);
          double dEtaMax = dEtaMin + sqrt(rGapS*rGapS-fDeltaPhiPiPi*fDeltaPhiPiPi);
          double eta = TMath::Abs( fDeltaEta );
          if( eta > dEtaMin && eta < dEtaMax ){
            //if( eta > dEtaMin && eta < dEtaMax && fDeltaEta <0 ){
            // xxx
            // fhistos->hDRNearPtMoon[fTyp][fCentralityBin][ZBin][ir1][irs][fpttBin]->Fill( fpta, fGeometricAcceptanceCorrection * fTrackPairEfficiency );
            
            if( fTyp == 0 )
              fhistos->fhDRNearPtMoon[fCentralityBin][ZBin][ir1][irs][fpttBin]->Fill( fpta, fGeometricAcceptanceCorrection * fTrackPairEfficiency );
            else
              fhistos->fhDRNearPtMoonM[fCentralityBin][ZBin][ir1][irs][fpttBin]->Fill( fpta, fGeometricAcceptanceCorrection * fTrackPairEfficiency );
            if(fTyp == kReal && fhistos->Is2DHistosEnabled())       fhistos->fhDphiAssoc2D[ir1][irs]->Fill( fDeltaEta, fDeltaPhi/kJPi, fGeometricAcceptanceCorrection * fTrackPairEfficiency );
          }
        }
      }
    }
  }
  
  if(fRGapBinAway>=0){
    if(fRGapBinAway <= fRSignalBin) fhistos->fhDRFarPt[fTyp][fCentralityBin][ZBin][fRGapBinAway][fpttBin]->Fill( fpta, fGeometricAcceptanceCorrection * fTrackPairEfficiency );
    // - moon -
    if(fRGapBinAway>0){
      for( int irs = 0; irs <= fRSignalBin;irs++ ){
        if( fRGapBinAway < irs ) continue;
        if( fPhiGapBinAway > irs ) continue;
        for ( int ir1=0;ir1<= fRGapBinAway;ir1++ ){
          if( irs > ir1 ) continue;
          if( ir1 > fRGapBinAway ) continue;
          double rGapS= fRGap[irs+1];
          double rGap1= fRGap[ir1+1];
          if( rGap1 < TMath::Abs(fDeltaPhi-kJPi) ) continue;
          if( rGapS < TMath::Abs(fDeltaPhi-kJPi) ) continue;
          double dEtaMin = sqrt(rGap1*rGap1- (fDeltaPhi-kJPi)*(fDeltaPhi-kJPi) );
          double dEtaMax = dEtaMin + sqrt(rGapS*rGapS-(fDeltaPhi-kJPi)*(fDeltaPhi-kJPi) );
          double eta = TMath::Abs( fDeltaEta );
          if( eta > dEtaMin && eta < dEtaMax ){
            //if( eta > dEtaMin && eta < dEtaMax && fDeltaEta <0 ){
            // xxx
            //                        fhistos->hDRFarPtMoon[fTyp][fCentralityBin][ZBin][ir1][irs][fpttBin]->Fill( fpta, fGeometricAcceptanceCorrection * fTrackPairEfficiency );
            if( fTyp == 0 )
              fhistos->fhDRFarPtMoon[fCentralityBin][ZBin][ir1][irs][fpttBin]->Fill( fpta, fGeometricAcceptanceCorrection * fTrackPairEfficiency );
            else
              fhistos->fhDRFarPtMoonM[fCentralityBin][ZBin][ir1][irs][fpttBin]->Fill( fpta, fGeometricAcceptanceCorrection * fTrackPairEfficiency );
            
            if(fTyp == kReal)       fhistos->fhDphiAssoc2D[ir1][irs]->Fill( fDeltaEta, fDeltaPhi/kJPi, fGeometricAcceptanceCorrection * fTrackPairEfficiency );
          }
        }
      }
    }
  }
}


double AliJCorrelations::GetGeoAccCorrFlat(double deltaEta){
  //FK// calculate acceptance correction on pseudorapidity triangle
  double absDEta = fabs(deltaEta);
  
  double denominator = 1 - absDEta/(2*fmaxEtaRange);
  //double denominator = 1 - (absDEta - ftriggFiducCut)/(2*fmaxEtaRange-2*ftriggFiducCut);//When Fid>0 max_Deta je mensi nez 2*EtaMax
  //double denominator = 1 - (absDEta - ftriggFiducCut)/(2*fmaxEtaRange-ftriggFiducCut);
  if(denominator > 1e-6)
    return 1.0/denominator;
  else
    return 0;
}

double AliJCorrelations::GetGeoAccCorrIncl(double deltaEta, int assocBin, int assocType){
  // histo filler
  
  if(assocBin < 0) return GetGeoAccCorrFlat(deltaEta);
  
  // Define the acceptance histogram
  TH1D *acceptanceHistogram;
  
  // Choose different acceptance histogram for xlong bins and other bins
  if(assocType == 1){
    acceptanceHistogram = fhistos->fhDEta3DNearMixFromFile[fCentralityBin][fpttBin][assocBin];
  } else {
    acceptanceHistogram = fhistos->fhDEtaNearMixFromFile[fCentralityBin][fpttBin][assocBin];
  }
  
  if(acceptanceHistogram->GetEntries()<1000) return GetGeoAccCorrFlat(deltaEta);
  
  int bin =  acceptanceHistogram->FindBin(deltaEta);
  double denominator  =  acceptanceHistogram->GetBinContent(bin);
  if(denominator > 1e-6)
    return 1.0/denominator;
  else
    return 0;
  
}

double AliJCorrelations::DeltaPhi(double phi1, double phi2) {
  // dphi
  double res =  atan2(sin(phi1-phi2), cos(phi1-phi2));
  //double res = phi1 - phi2;
  //return res>-kJPi/3.0 ? res : kJTwoPi+res ;
  return res>-kJPi*9./20. ? res : kJTwoPi+res ;
}





















