#include "AliEmcalJetUtilityGenSubtractor.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliEmcalJetTask.h"

ClassImp(AliEmcalJetUtilityGenSubtractor)

//______________________________________________________________________________
AliEmcalJetUtilityGenSubtractor::AliEmcalJetUtilityGenSubtractor() :
AliEmcalJetUtility(),
  fDoGenericSubtractionJetMass(kFALSE),
  fDoGenericSubtractionGR(kFALSE),
  fDoGenericSubtractionExtraJetShapes(kFALSE),
  fUseExternalBkg(kFALSE),
  fRhoName(""),
  fRhomName(""),
  fRho(0),
  fRhom(0),
  fRMax(0.4),
  fDRStep(0.04),
  fPtMinGR(40.),
  fRhoParam(0),
  fRhomParam(0)
{
  // Dummy constructor.

}

//______________________________________________________________________________
AliEmcalJetUtilityGenSubtractor::AliEmcalJetUtilityGenSubtractor(const char* name) :
  AliEmcalJetUtility(name),
  fDoGenericSubtractionJetMass(kFALSE),
  fDoGenericSubtractionGR(kFALSE),
  fDoGenericSubtractionExtraJetShapes(kFALSE),
  fUseExternalBkg(kFALSE),
  fRhoName(""),
  fRhomName(""),
  fRho(0),
  fRhom(0),
  fRMax(0.4),
  fDRStep(0.04),
  fPtMinGR(40.),
  fRhoParam(0),
  fRhomParam(0)
{
  // Default constructor.
}

//______________________________________________________________________________
AliEmcalJetUtilityGenSubtractor::AliEmcalJetUtilityGenSubtractor(const AliEmcalJetUtilityGenSubtractor &other) :
  AliEmcalJetUtility(other),
  fDoGenericSubtractionJetMass(other.fDoGenericSubtractionJetMass),
  fDoGenericSubtractionGR(other.fDoGenericSubtractionGR),
  fDoGenericSubtractionExtraJetShapes(other.fDoGenericSubtractionExtraJetShapes),
  fUseExternalBkg(other.fUseExternalBkg),
  fRhoName(other.fRhoName),
  fRhomName(other.fRhomName),
  fRho(other.fRho),
  fRhom(other.fRhom),
  fRMax(other.fRMax),
  fDRStep(other.fDRStep),
  fPtMinGR(other.fPtMinGR),
  fRhoParam(other.fRhoParam),
  fRhomParam(other.fRhomParam)
{
  // Copy constructor.
}

//______________________________________________________________________________
AliEmcalJetUtilityGenSubtractor& AliEmcalJetUtilityGenSubtractor::operator=(const AliEmcalJetUtilityGenSubtractor &other)
{
  // Assignment.

  if (&other == this) return *this;
  AliEmcalJetUtility::operator=(other);
  fDoGenericSubtractionJetMass = other.fDoGenericSubtractionJetMass;
  fDoGenericSubtractionGR = other.fDoGenericSubtractionGR;
  fDoGenericSubtractionExtraJetShapes = other.fDoGenericSubtractionExtraJetShapes;
  fUseExternalBkg = other.fUseExternalBkg;
  fRhoName = other.fRhoName;
  fRhomName = other.fRhomName;
  fRho = other.fRho;
  fRhom = other.fRhom;
  fRMax = other.fRMax;
  fDRStep = other.fDRStep;
  fPtMinGR = other.fPtMinGR;
  fRhoParam = other.fRhoParam;
  fRhomParam = other.fRhomParam;
  return *this;
}

//______________________________________________________________________________
void AliEmcalJetUtilityGenSubtractor::Init()
{
  // Initialize the utility.

  if (!fRhoName.IsNull() && !fRhoParam) { // get rho from the event
    if(!fJetTask) return;
    fRhoParam = dynamic_cast<AliRhoParameter*>(fJetTask->InputEvent()->FindListObject(fRhoName));
    if (!fRhoParam) {
      AliError(Form("%s: Could not retrieve rho %s!", GetName(), fRhoName.Data()));
      return;
    }
  }
  
  if (!fRhomName.IsNull() && !fRhomParam) { // get rhom from the event
    fRhomParam = dynamic_cast<AliRhoParameter*>(fJetTask->InputEvent()->FindListObject(fRhomName));
    if (!fRhomParam) {
      AliError(Form("%s: Could not retrieve rho_m %s!", GetName(), fRhomName.Data()));
      return;
    }
  }

  fInit = kTRUE;
}

//______________________________________________________________________________
void AliEmcalJetUtilityGenSubtractor::Prepare(AliFJWrapper& fjw)
{
  // Prepare the utility.

  if (!fInit) return;
  
  if (fRhoParam) fRho = fRhoParam->GetVal();
  if (fRhomParam) fRhom = fRhomParam->GetVal();

  //run generic subtractor
  if (fDoGenericSubtractionJetMass) {
    fjw.SetUseExternalBkg(fUseExternalBkg,fRho,fRhom);
    fjw.DoGenericSubtractionJetMass();
  }
 
 if (fDoGenericSubtractionExtraJetShapes) {
   fjw.SetUseExternalBkg(fUseExternalBkg,fRho,fRhom);
   fjw.DoGenericSubtractionJetAngularity();
   fjw.DoGenericSubtractionJetpTD();
   fjw.DoGenericSubtractionJetCircularity();
   fjw.DoGenericSubtractionJetSigma2();
   fjw.DoGenericSubtractionJetConstituent();
   fjw.DoGenericSubtractionJetLeSub();
 }
}

//______________________________________________________________________________
void AliEmcalJetUtilityGenSubtractor::ProcessJet(AliEmcalJet* jet, Int_t ij, AliFJWrapper& fjw)
{
  // Proceess each jet.

  if (!fInit) return;

#ifdef FASTJET_VERSION

  if (fDoGenericSubtractionJetMass) {
    std::vector<fastjet::contrib::GenericSubtractorInfo> jetMassInfo = fjw.GetGenSubtractorInfoJetMass();
    Int_t n = (Int_t)jetMassInfo.size();
    if(n > ij && n > 0) {
      jet->SetFirstDerivative(jetMassInfo[ij].first_derivative());
      jet->SetSecondDerivative(jetMassInfo[ij].second_derivative());
      jet->SetFirstOrderSubtracted(jetMassInfo[ij].first_order_subtracted());
      jet->SetSecondOrderSubtracted(jetMassInfo[ij].second_order_subtracted());
    }
  }

  //here do generic subtraction for angular structure function
  Double_t ptcorr = jet->Pt()-fjw.GetJetArea(ij)*fRho;
  if (fDoGenericSubtractionGR && ptcorr>fPtMinGR) {
    fjw.SetUseExternalBkg(fUseExternalBkg, fRho, fRhom);
    fRMax = fJetTask->GetRadius()+0.2;
    fjw.SetRMaxAndStep(fRMax, fDRStep);
    fjw.DoGenericSubtractionGR(ij);
    std::vector<double> num = fjw.GetGRNumerator();
    std::vector<double> den = fjw.GetGRDenominator();
    std::vector<double> nums = fjw.GetGRNumeratorSub();
    std::vector<double> dens = fjw.GetGRDenominatorSub();
    //pass this to AliEmcalJet
    jet->SetGRNumSize(num.size());
    jet->SetGRDenSize(den.size());
    jet->SetGRNumSubSize(nums.size());
    jet->SetGRDenSubSize(dens.size());
    Int_t nsize = (Int_t)num.size();
    for (Int_t g = 0; g < nsize; ++g) {
      jet->AddGRNumAt(num[g],g);
      jet->AddGRNumSubAt(nums[g],g);
    }
    Int_t dsize = (Int_t)den.size();
    for (Int_t g = 0; g < dsize; ++g) {
      jet->AddGRDenAt(den[g], g);
      jet->AddGRDenSubAt(dens[g], g);
    }
  }

  if (fDoGenericSubtractionExtraJetShapes) {
    std::vector<fastjet::contrib::GenericSubtractorInfo> jetAngularityInfo = fjw.GetGenSubtractorInfoJetAngularity();
    Int_t na = (Int_t)jetAngularityInfo.size();
    if(na > ij && na > 0) {
      jet->SetFirstDerivativeAngularity(jetAngularityInfo[ij].first_derivative());
      jet->SetSecondDerivativeAngularity(jetAngularityInfo[ij].second_derivative());
      jet->SetFirstOrderSubtractedAngularity(jetAngularityInfo[ij].first_order_subtracted());
      jet->SetSecondOrderSubtractedAngularity(jetAngularityInfo[ij].second_order_subtracted());
    }

    std::vector<fastjet::contrib::GenericSubtractorInfo> jetpTDInfo = fjw.GetGenSubtractorInfoJetpTD();
    Int_t np = (Int_t)jetpTDInfo.size();
    if(np > ij && np > 0) {
      jet->SetFirstDerivativepTD(jetpTDInfo[ij].first_derivative());
      jet->SetSecondDerivativepTD(jetpTDInfo[ij].second_derivative());
      jet->SetFirstOrderSubtractedpTD(jetpTDInfo[ij].first_order_subtracted());
      jet->SetSecondOrderSubtractedpTD(jetpTDInfo[ij].second_order_subtracted());
    }

    std::vector<fastjet::contrib::GenericSubtractorInfo> jetCircularityInfo = fjw.GetGenSubtractorInfoJetCircularity();
    Int_t nc = (Int_t)jetCircularityInfo.size();
    if(nc > ij && nc > 0) {
      jet->SetFirstDerivativeCircularity(jetCircularityInfo[ij].first_derivative());
      jet->SetSecondDerivativeCircularity(jetCircularityInfo[ij].second_derivative());
      jet->SetFirstOrderSubtractedCircularity(jetCircularityInfo[ij].first_order_subtracted());
      jet->SetSecondOrderSubtractedCircularity(jetCircularityInfo[ij].second_order_subtracted());
    }

    std::vector<fastjet::contrib::GenericSubtractorInfo> jetSigma2Info = fjw.GetGenSubtractorInfoJetSigma2();
    Int_t ns = (Int_t)jetSigma2Info.size();
    if (ns > ij && ns > 0) {
      jet->SetFirstDerivativeSigma2(jetSigma2Info[ij].first_derivative());
      jet->SetSecondDerivativeSigma2(jetSigma2Info[ij].second_derivative());
      jet->SetFirstOrderSubtractedSigma2(jetSigma2Info[ij].first_order_subtracted());
      jet->SetSecondOrderSubtractedSigma2(jetSigma2Info[ij].second_order_subtracted());
    }


    std::vector<fastjet::contrib::GenericSubtractorInfo> jetConstituentInfo = fjw.GetGenSubtractorInfoJetConstituent();
    Int_t nco = (Int_t)jetConstituentInfo.size();
    if(nco > ij && nco > 0) {
      jet->SetFirstDerivativeConstituent(jetConstituentInfo[ij].first_derivative());
      jet->SetSecondDerivativeConstituent(jetConstituentInfo[ij].second_derivative());
      jet->SetFirstOrderSubtractedConstituent(jetConstituentInfo[ij].first_order_subtracted());
      jet->SetSecondOrderSubtractedConstituent(jetConstituentInfo[ij].second_order_subtracted());
    }

    std::vector<fastjet::contrib::GenericSubtractorInfo> jetLeSubInfo = fjw.GetGenSubtractorInfoJetLeSub();
    Int_t nlsub = (Int_t)jetLeSubInfo.size();
    if(nlsub > ij && nlsub > 0) {
      jet->SetFirstDerivativeLeSub(jetLeSubInfo[ij].first_derivative());
      jet->SetSecondDerivativeLeSub(jetLeSubInfo[ij].second_derivative());
      jet->SetFirstOrderSubtractedLeSub(jetLeSubInfo[ij].first_order_subtracted());
      jet->SetSecondOrderSubtractedLeSub(jetLeSubInfo[ij].second_order_subtracted());
    }
  }

#endif
}

//______________________________________________________________________________
void AliEmcalJetUtilityGenSubtractor::Terminate(AliFJWrapper& /*fjw*/)
{
  // Run termination of the utility (after each event).
}
