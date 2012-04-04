/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliMuonTrackCuts.h"

#include "TMath.h"
#include "TList.h"
#include "TArrayD.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TFile.h"
#include "TParameter.h"
#include "TKey.h"
#include "TVector3.h"
#include "TSystem.h"

#include "AliLog.h"
#include "AliVParticle.h"
#include "AliESDMuonTrack.h"
#include "AliAODTrack.h"
#include "AliAnalysisManager.h"

/// \cond CLASSIMP
ClassImp(AliMuonTrackCuts) // Class implementation in ROOT context
/// \endcond


//________________________________________________________________________
AliMuonTrackCuts::AliMuonTrackCuts() :
  AliAnalysisCuts(),
  fIsESD(kFALSE),
  fIsMC(kFALSE),
  fUseCustomParam(kFALSE),
  fSharpPtCut(kFALSE),
  fParameters(TArrayD(kNParameters))
{
  /// Default ctor.
  fParameters.Reset();
}

//________________________________________________________________________
AliMuonTrackCuts::AliMuonTrackCuts(const char* name, const char* title ) :
AliAnalysisCuts(name, title),
fIsESD(kFALSE),
fIsMC(kFALSE),
fUseCustomParam(kFALSE),
fSharpPtCut(kFALSE),
fParameters(TArrayD(kNParameters))
{
  /// Constructor
  fParameters.Reset();
  SetDefaultFilterMask();
}


//________________________________________________________________________
AliMuonTrackCuts::AliMuonTrackCuts(const char* name, const char* title, Bool_t isESD ) :
  AliAnalysisCuts(name, title),
  fIsESD(isESD),
  fIsMC(kFALSE),
  fUseCustomParam(kFALSE),
  fSharpPtCut(kFALSE),
  fParameters(TArrayD(kNParameters))
{
  /// Obsolete Constructor
  AliWarning(Form("\n\n *****  This constructor is obsolete and will be removed!\nPlease use AliMuonTrackCuts(name, title) instead!\n"));
  fParameters.Reset();
  SetDefaultFilterMask();
}


//________________________________________________________________________
AliMuonTrackCuts::AliMuonTrackCuts(const AliMuonTrackCuts& obj) :
  AliAnalysisCuts(obj),
  fIsESD(obj.fIsESD),
  fIsMC(obj.fIsMC),
  fUseCustomParam(obj.fUseCustomParam),
  fSharpPtCut(obj.fSharpPtCut),
  fParameters(obj.fParameters)
{
  /// Copy constructor
}


//________________________________________________________________________
AliMuonTrackCuts& AliMuonTrackCuts::operator=(const AliMuonTrackCuts& obj)
{
  /// Assignment operator
  if ( this != &obj ) { 
    AliAnalysisCuts::operator=(obj);
    fIsESD = obj.fIsESD;
    fIsMC = obj.fIsMC;
    fUseCustomParam = obj.fUseCustomParam;
    fSharpPtCut = obj.fSharpPtCut;
    fParameters = obj.fParameters;
  }
  return *this;
}


//________________________________________________________________________
AliMuonTrackCuts::~AliMuonTrackCuts()
{
  /// Destructor
}

//________________________________________________________________________
Bool_t AliMuonTrackCuts::IsESDTrack( const AliVParticle* track ) const
{
  /// Check if track is from ESD or AOD
  return ( track->IsA() != AliAODTrack::Class() );
}

//________________________________________________________________________
Bool_t AliMuonTrackCuts::RunMatchesRange( Int_t runNumber, const Char_t* objName ) const
{
  /// Check if the object contains the run
  TString sname(objName);
  TObjArray* array = sname.Tokenize("_");
  array->SetOwner();
  Int_t runRange[2] = { -1, -1 };
  if ( array->GetEntries() >= 3 ) {
    for ( Int_t irun=0; irun<2; ++irun ) {
      TString currRun = array->At(irun+1)->GetName();
      if ( currRun.IsDigit() ) runRange[irun] = currRun.Atoi();
    }
  }
  delete array;
  return ( runNumber >= runRange[0] && runNumber <= runRange[1]);
}

//________________________________________________________________________
void AliMuonTrackCuts::SetUseCustomParam( Bool_t useCustomParam, Int_t runNumber  )
{
  /// Flag to select custom parameters
  /// It first searches the default parameters in OADB
  /// then disables the access to the OADB
  /// and allows to manually modify parameters

  if ( ! fUseCustomParam && useCustomParam ) SetRun(runNumber);
  fUseCustomParam = useCustomParam;
}

//________________________________________________________________________
Bool_t AliMuonTrackCuts::SetRun( Int_t runNumber )
{
  /// Get parameters from OADB for runNumber
  
  if ( fUseCustomParam ) return kFALSE;
  return StreamParameters(runNumber, -1);
}


//________________________________________________________________________
Bool_t AliMuonTrackCuts::StreamParameters( Int_t runNumber,  Int_t runMax )
{
  if ( runMax > 0 ) { // Stream to OADB
    if ( ! fUseCustomParam ) {
      AliError("Users are not allowed to update OADB. Use SetUseCustomParam() instead");
      return kFALSE;
    }
  }

  TString filename = Form("%s/PWG/MUON/MuonTrackCuts.root",AliAnalysisManager::GetOADBPath());
  if ( fIsMC ) filename.ReplaceAll(".root", "_MC.root");

  TString parNames[kNParameters];
  parNames[kMeanDcaX]       = "MeanDcaX";
  parNames[kMeanDcaY]       = "MeanDcaY";
  parNames[kMeanDcaZ]       = "MeanDcaZ";
  parNames[kMeanPCorr23]    = "MeanPCorr23";
  parNames[kMeanPCorr310]   = "MeanPCorr310";
  parNames[kSigmaPdca23]    = "SigmaPdca23";
  parNames[kSigmaPdca310]   = "SigmaPdca310";
  parNames[kNSigmaPdcaCut]  = "NSigmaPdcaCut";
  parNames[kChi2NormCut]    = "Chi2NormCut";
  parNames[kRelPResolution] = "RelPResolution";
  parNames[kSharpPtApt]     = "SharpPtApt";
  parNames[kSharpPtLpt]     = "SharpPtLpt";
  parNames[kSharpPtHpt]     = "SharpPtHpt";

  TObjArray* paramList = 0x0;

  if ( runMax < 0 ) { // Get from OADB
    TFile* file = TFile::Open(filename.Data(), "READ");
    if ( ! file ) {
      AliError(Form("OADB file %s not found!", filename.Data()));
      return kFALSE;
    }

    TList* listOfKeys = file->GetListOfKeys();
    TIter next(listOfKeys);
    TObject* key = 0x0;
    Bool_t foundMatch = kFALSE;
    TObject* defaultObj = 0x0;
    while ( ( key = next() ) ) {
      TString objName = key->GetName();
      objName.ToUpper();
      if ( RunMatchesRange(runNumber, objName.Data()) ) {
        paramList = static_cast<TObjArray*>(file->Get(key->GetName()));
        foundMatch = kTRUE;
        break;
      }
      if ( objName.Contains("DEFAULT") ) defaultObj = file->Get(key->GetName());
    }

    if ( ! foundMatch ) {
      AliWarning("Run number not found in OADB: using default");
      if ( defaultObj ) paramList = static_cast<TObjArray*>(defaultObj);
      else {
        file->Close();
        AliError("Default parameters not found in OADB!");
        return kFALSE;
      }
    }

    AliInfo(Form("Required run %i. Param. set: %s", runNumber, paramList->GetName()));

    for ( Int_t ipar=0; ipar<kNParameters; ++ipar ) {
      TParameter<Double_t>* param = static_cast<TParameter<Double_t>*>(paramList->FindObject(parNames[ipar].Data()));
      if ( ! param ) {
        AliWarning(Form("Parameter %s not set", parNames[ipar].Data()));
        continue;
      }
      fParameters[ipar] = param->GetVal();
    }

    file->Close();
  }
  else {
    if ( ! paramList ) {
      paramList = new TObjArray(kNParameters);
      paramList->SetOwner();
    }
    for ( Int_t ipar=0; ipar<kNParameters; ++ipar ) {
      TParameter<Double_t>* param= new TParameter<Double_t>(parNames[ipar].Data(), fParameters[ipar]);
      paramList->AddAt(param, ipar);
    }

    TString paramListName = "MuonCuts_";
    paramListName += ( runNumber < 0 ) ? "Default" : Form("%i_%i",runNumber, runMax);
    AliInfo(Form("Adding %s to file %s", paramListName.Data(), filename.Data()));
    paramList->SetName(paramListName.Data());
    TFile* file = TFile::Open(filename.Data(), "UPDATE");
    paramList->Write(paramListName.Data(), TObject::kSingleKey);
    file->Close();
    delete paramList;
  }
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliMuonTrackCuts::SetParameter(Int_t iparam, Float_t value)
{
  /// Set parameter
  if ( fUseCustomParam ) {
    fParameters.SetAt(value, iparam);
    return kTRUE;
  }

  AliWarning("Parameters automatically taken from OADB. If you want to use with custom parameters, use SetUseCustomParam()");
  return kFALSE;
}


//________________________________________________________________________
Bool_t AliMuonTrackCuts::IsSelected( TObject* obj )
{
  /// Track is selected
  UInt_t filterMask = GetFilterMask();
  UInt_t selectionMask = GetSelectionMask(obj);
  
  return ( ( selectionMask & filterMask ) == filterMask );
}


//________________________________________________________________________
UInt_t AliMuonTrackCuts::GetSelectionMask( const TObject* obj )
{
  /// Get selection mask
  
  const AliVParticle* track = static_cast<const AliVParticle*> ( obj );
  
  Bool_t isESD = IsESDTrack(track);

  UInt_t selectionMask = 0;

  Bool_t isMuon = ( isESD ) ? ((AliESDMuonTrack*)track)->ContainTrackerData() :  ((AliAODTrack*)track)->IsMuonTrack();

  if ( ! isMuon ) return selectionMask;

  Double_t eta = track->Eta();
  if ( eta > -4. && eta < -2.5 ) selectionMask |= kMuEta;

  Double_t rAbsEnd =  ( isESD ) ? ((AliESDMuonTrack*)track)->GetRAtAbsorberEnd() : ((AliAODTrack*)track)->GetRAtAbsorberEnd();
  Double_t thetaAbsEndDeg = TMath::ATan( rAbsEnd / 505. ) * TMath::RadToDeg();

  if ( thetaAbsEndDeg > 2. && thetaAbsEndDeg < 10. ) selectionMask |= kMuThetaAbs;

  Int_t matchTrig = ( isESD ) ? ((AliESDMuonTrack*)track)->GetMatchTrigger() : ((AliAODTrack*)track)->GetMatchTrigger();
  Int_t cutLevel[3] = {kMuMatchApt, kMuMatchLpt, kMuMatchHpt};
  Double_t pt = track->Pt();
  for ( Int_t ilevel=0; ilevel<3; ilevel++ ) {
    if ( matchTrig < ilevel+1 ) break;
    if ( fSharpPtCut && pt < GetSharpPtCut(ilevel) ) break;
    selectionMask |= cutLevel[ilevel];
  }

  Double_t chi2norm = ( isESD ) ? ((AliESDMuonTrack*)track)->GetNormalizedChi2() : ((AliAODTrack*)track)->Chi2perNDF();
  if ( chi2norm < GetChi2NormCut() ) selectionMask |= kMuTrackChiSquare;

  TVector3 dcaAtVz = GetCorrectedDCA(track);
  Double_t pTotMean = GetAverageMomentum(track);

  Double_t pDca = pTotMean * dcaAtVz.Mag();
    
  Double_t pTot = track->P();
  
  // Momentum resolution only
  // The cut depends on the momentum resolution. In particular:
  // Sigma_pDCA = Sqrt( Sigma_pDCA_measured^2 + Sigma_p/p * pDCA)
  // The relative momentum distribution Sigma_p/p is estimated from the track Delta_p,
  // and from the error on sagitta Delta_s:
  // Delta_p/p = p*Delta_s/(1+p*Delta_s)
  // A cut at N sigmas requres a cut in N Delta_s, so
  // Sigma_p/p(N) = p*N*Delta_s/(1+p*N*Delta_s)
  //Double_t pResolutionEffect = pDca * pTot * GetRelPResolution() / ( 1. + GetNSigmaPdca() * GetRelPResolution()*pTot );
  
  //Double_t pResolutionEffect = 0.4 * pTot;  // Values used in 2010 data
  //Double_t pResolutionEffect = 0.32 * pTot; // Values in 2011
  //Double_t sigmaPdca = GetSigmaPdca(rAbsEnd);
  //Double_t sigmaPdcaWithRes = TMath::Sqrt( sigmaPdca*sigmaPdca + pResolutionEffect*pResolutionEffect );
  
  // Momentum resolution and slope resolution 
  // Due to the momentum resolution, the measured momentum is biased
  // Since we want to keep as much signal as possible, we want to avoid
  // that a measured pxDCA is rejected since the momentum is overestimated
  // p_true = p_meas - Delta_p
  // p_true = p_meas - N*Delta_s*p_meas / (1+n*Delta_s*p_meas)
  // Hence:
  // p_true x DCA < N * Sigma_pDCA_meas =>
  // p_meas x DCA < N * Sigma_pDCA_meas / ( 1 - N*Delta_s*p_meas / (1+n*Delta_s*p_meas))
  // Finally the cut value has to be summed in quadrature with the error on DCA,
  // which is given by the slope resolution
  // p_meas x DCA < N * Sqrt( ( Sigma_pDCA_meas / ( 1 - N*Delta_s*p_meas / (1+n*Delta_s*p_meas)) )^2 + (distance * sigma_slope * p_meas )^2)
  Double_t nrp = GetNSigmaPdca() * GetRelPResolution() * pTot;
  Double_t pResolutionEffect = GetSigmaPdca(rAbsEnd) / ( 1. - nrp / ( 1. + nrp ) );
  Double_t slopeResolutionEffect = 535. * GetSlopeResolution() * pTot;
  
  Double_t sigmaPdcaWithRes = TMath::Sqrt( pResolutionEffect*pResolutionEffect + slopeResolutionEffect*slopeResolutionEffect );
  
  if ( pDca < GetNSigmaPdca() * sigmaPdcaWithRes ) selectionMask |= kMuPdca;
  
  AliDebug(1, Form("Selection mask 0x%x\n", selectionMask));

  return selectionMask;
}


//________________________________________________________________________
Bool_t AliMuonTrackCuts::IsSelected( TList* /* list */)
{
  /// Not implemented
  AliError("Function not implemented: Use IsSelected(TObject*)");
  return kFALSE;
}


//________________________________________________________________________
Int_t AliMuonTrackCuts::GetThetaAbsBin ( Double_t rAtAbsEnd ) const
{
  /// Get theta abs bin
  Double_t thetaAbsEndDeg = TMath::ATan( rAtAbsEnd / 505. ) * TMath::RadToDeg();
  Int_t thetaAbsBin = ( thetaAbsEndDeg < 3. ) ? kThetaAbs23 : kThetaAbs310;
  return thetaAbsBin;
}


//________________________________________________________________________
TVector3 AliMuonTrackCuts::GetCorrectedDCA ( const AliVParticle* track ) const
{
  /// Get corrected DCA

  Bool_t isESD = IsESDTrack(track);
  
  Double_t vtxPos[3];
  if ( isESD ) {
    vtxPos[0] = ((AliESDMuonTrack*)track)->GetNonBendingCoor();
    vtxPos[1] = ((AliESDMuonTrack*)track)->GetBendingCoor();
    vtxPos[2] = ((AliESDMuonTrack*)track)->GetZ();
  }
  else ((AliAODTrack*)track)->GetXYZ(vtxPos);
  
  TVector3 vertex(vtxPos);
  
  TVector3 dcaTrack(0.,0., vtxPos[2]);
  dcaTrack.SetX( ( isESD ) ? ((AliESDMuonTrack*)track)->GetNonBendingCoorAtDCA() : ((AliAODTrack*)track)->XAtDCA() );
  dcaTrack.SetY( ( isESD ) ? ((AliESDMuonTrack*)track)->GetBendingCoorAtDCA() : ((AliAODTrack*)track)->YAtDCA() );
  
  TVector3 dcaAtVz = dcaTrack - vertex - GetMeanDCA();

  return dcaAtVz;
}

//________________________________________________________________________
Double_t AliMuonTrackCuts::GetAverageMomentum ( const AliVParticle* track ) const
{
  /// Get average momentum before and after the absorber

  Bool_t isESD = IsESDTrack(track);
  
  Double_t pTotMean = 0.;
  Double_t pTot = track->P();
  //if ( isESD ) pTotMean = 0.5 * ( pTot + ((AliESDMuonTrack*)track)->PUncorrected() );
  if ( isESD ) pTotMean = ((AliESDMuonTrack*)track)->PUncorrected(); // Increased stability if using uncorrected value
  else {
    pTotMean = pTot - GetMeanPCorr(((AliAODTrack*)track)->GetRAtAbsorberEnd());
  }

  return pTotMean;
}


//________________________________________________________________________
void AliMuonTrackCuts::SetMeanDCA ( Double_t xAtDca, Double_t yAtDca, Double_t zAtDca )
{
    /// Set mean DCA from track
  SetParameter(kMeanDcaX, xAtDca);
  SetParameter(kMeanDcaY, yAtDca);
  SetParameter(kMeanDcaZ, zAtDca);
}

//________________________________________________________________________
TVector3 AliMuonTrackCuts::GetMeanDCA () const
{ 
    /// Get mean DCA from track
  return TVector3(fParameters[kMeanDcaX], fParameters[kMeanDcaY], fParameters[kMeanDcaZ]);
}

//________________________________________________________________________
void AliMuonTrackCuts::SetMeanPCorr ( Double_t pCorrThetaAbs23, Double_t pCorrThetaAbs310 )
{
  /// Set mean p correction
  SetParameter(kMeanPCorr23, pCorrThetaAbs23);
  SetParameter(kMeanPCorr310, pCorrThetaAbs310);
}

//________________________________________________________________________
Double_t AliMuonTrackCuts::GetMeanPCorr ( Double_t rAtAbsEnd ) const
{
  /// Get mean p correction
  return fParameters[kMeanPCorr23+GetThetaAbsBin(rAtAbsEnd)];
}

//________________________________________________________________________
void AliMuonTrackCuts::SetSigmaPdca ( Double_t sigmaThetaAbs23, Double_t sigmaThetaAbs310 )
{ 
  /// Set sigma pdca
  SetParameter(kSigmaPdca23, sigmaThetaAbs23);
  SetParameter(kSigmaPdca310, sigmaThetaAbs310);
}

//________________________________________________________________________
Double_t AliMuonTrackCuts::GetSigmaPdca ( Double_t rAtAbsEnd ) const
{ 
  /// Get mean pdca
  return fParameters[kSigmaPdca23+GetThetaAbsBin(rAtAbsEnd)];
}

//________________________________________________________________________
void AliMuonTrackCuts::SetNSigmaPdca ( Double_t nSigmas )
{ 
  /// Set N sigma pdca cut
  SetParameter(kNSigmaPdcaCut, nSigmas);
}

//________________________________________________________________________
Double_t AliMuonTrackCuts::GetNSigmaPdca () const
{
  /// Get N sigma pdca cut
  return fParameters[kNSigmaPdcaCut];
}

//________________________________________________________________________
void AliMuonTrackCuts::SetChi2NormCut ( Double_t chi2normCut )
{
  /// Set cut on normalized chi2 of tracks
  SetParameter(kChi2NormCut,chi2normCut);
}

//________________________________________________________________________
Double_t AliMuonTrackCuts::GetChi2NormCut () const
{
  /// Get cut on normalized chi2 of tracks
  return fParameters[kChi2NormCut];
}

//________________________________________________________________________
void AliMuonTrackCuts::SetRelPResolution ( Double_t relPResolution )
{
  /// Set relative momentum resolution
  SetParameter(kRelPResolution,relPResolution);
}

//________________________________________________________________________
Double_t AliMuonTrackCuts::GetRelPResolution () const
{
  /// Get relative momentum resolution
  return fParameters[kRelPResolution];
}


//________________________________________________________________________
void AliMuonTrackCuts::SetSlopeResolution ( Double_t slopeResolution )
{
  /// Set slope resolution
  SetParameter(kSlopeResolution,slopeResolution);
}

//________________________________________________________________________
Double_t AliMuonTrackCuts::GetSlopeResolution () const
{
  /// Get slope resolution
  return fParameters[kSlopeResolution];
}

//________________________________________________________________________
void AliMuonTrackCuts::SetSharpPtCut ( Double_t valueApt, Double_t valueLpt, Double_t valueHpt  )
{
  /// Set sharp tracker cut matching the trigger level

  SetParameter(kSharpPtApt, valueApt);
  SetParameter(kSharpPtLpt, valueLpt);
  SetParameter(kSharpPtHpt, valueHpt);
}

//________________________________________________________________________
Double_t AliMuonTrackCuts::GetSharpPtCut ( Int_t trigPtCut, Bool_t warn ) const
{
  /// Get sharp tracker cut matching the trigger level
  /// trigPtCut can be 0 (Apt), 1 (Lpt) or 2 (Hpt)
  if ( trigPtCut < 0 || trigPtCut > 2 ) {
    if ( warn ) AliError("Allowed values for trigPtCut are 0 (Apt), 1 (Lpt), 2 (Hpt)");
    return 0.;
  }
  Int_t ipar = kSharpPtApt + trigPtCut;
  return fParameters[ipar];
}

//________________________________________________________________________
void AliMuonTrackCuts::SetDefaultFilterMask ()
{
  /// Standard cuts for single muon
  SetFilterMask ( kMuEta | kMuThetaAbs | kMuPdca | kMuMatchApt );
}  

//________________________________________________________________________
void AliMuonTrackCuts::Print(Option_t* option) const
{
  //
  /// Print info
  //
  TString sopt(option);
  sopt.ToLower();
  if ( sopt.IsNull() || sopt.Contains("*") || sopt.Contains("all") ) sopt = "mask param";
  UInt_t filterMask = GetFilterMask();
  Int_t cutLevel[3] = {kMuMatchApt, kMuMatchLpt, kMuMatchHpt};
  TString cutLevelName[3] = {"Apt", "Lpt", "Hpt"};
  if ( sopt.Contains("mask") ) {
    printf(" *** Muon track filter mask: *** \n");
    printf("  0x%x\n", filterMask);
    if ( filterMask & kMuEta ) printf("  -4 < eta < -2.5\n");
    if ( filterMask & kMuThetaAbs ) printf("  2 < theta_abs < 10 deg\n");
    if ( filterMask & kMuPdca ) printf("  pxDCA cut\n");
    for ( Int_t ilevel=0; ilevel<3; ilevel++ ) {
      if ( filterMask & cutLevel[ilevel] ) {
        printf("  match %s", cutLevelName[ilevel].Data());
        if ( fSharpPtCut ) printf(" && sharp pt from tracker");
        printf("\n");
      }
    }
    if ( filterMask & kMuTrackChiSquare ) printf("  Chi2 cut on track\n");
    printf(" ******************** \n");
  }
  if ( sopt.Contains("param") ) {
    printf(" *** Muon track parameter summary: ***\n");
    printf("  Mean vertex DCA: (%g, %g, %g)\n", fParameters[kMeanDcaX], fParameters[kMeanDcaY], fParameters[kMeanDcaZ]);
    printf("  Mean p correction (GeV/c): theta2-3 = %g  theta3-10 = %g\n", fParameters[kMeanPCorr23], fParameters[kMeanPCorr310]);
    printf("  Sigma p x DCA (cm x GeV/c): theta2-3 = %g  theta3-10 = %g\n", fParameters[kSigmaPdca23], fParameters[kSigmaPdca310]);
    printf("  Cut p x DCA in units of sigma: %g\n", fParameters[kNSigmaPdcaCut]);
    printf("  Cut on track chi2/NDF: %g\n", fParameters[kChi2NormCut]);
    printf("  Momentum resolution: %g\n", fParameters[kRelPResolution]);
    printf("  Slope resolution: %g\n", fParameters[kSlopeResolution]);
    printf("  Sharp pt cut: %g (Apt)  %g (Lpt)  %g (Hpt)\n", fParameters[kSharpPtApt], fParameters[kSharpPtLpt], fParameters[kSharpPtHpt]);
    printf(" ********************************\n");
  }
}
