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
  fDoGenericSubtractionNsubjettiness(kFALSE),  
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
  fDoGenericSubtractionNsubjettiness(kFALSE), 
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
  fDoGenericSubtractionNsubjettiness(other.fDoGenericSubtractionNsubjettiness), 
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
  fDoGenericSubtractionNsubjettiness = other.fDoGenericSubtractionNsubjettiness;
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
 
 if  (fDoGenericSubtractionNsubjettiness) {
   fjw.SetUseExternalBkg(fUseExternalBkg,fRho,fRhom);
   fjw.DoGenericSubtractionJet1subjettiness_kt();
   fjw.DoGenericSubtractionJet2subjettiness_kt();
   fjw.DoGenericSubtractionJet3subjettiness_kt();
   fjw.DoGenericSubtractionJetOpeningAngle_kt();
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
      jet->GetShapeProperties()->SetFirstDerivative(jetMassInfo[ij].first_derivative());
      jet->GetShapeProperties()->SetSecondDerivative(jetMassInfo[ij].second_derivative());
      jet->GetShapeProperties()->SetFirstOrderSubtracted(jetMassInfo[ij].first_order_subtracted());
      jet->GetShapeProperties()->SetSecondOrderSubtracted(jetMassInfo[ij].second_order_subtracted());
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
    jet->GetShapeProperties()->SetGRNumSize(num.size());
    jet->GetShapeProperties()->SetGRDenSize(den.size());
    jet->GetShapeProperties()->SetGRNumSubSize(nums.size());
    jet->GetShapeProperties()->SetGRDenSubSize(dens.size());
    Int_t nsize = (Int_t)num.size();
    for (Int_t g = 0; g < nsize; ++g) {
      jet->GetShapeProperties()->AddGRNumAt(num[g],g);
      jet->GetShapeProperties()->AddGRNumSubAt(nums[g],g);
    }
    Int_t dsize = (Int_t)den.size();
    for (Int_t g = 0; g < dsize; ++g) {
      jet->GetShapeProperties()->AddGRDenAt(den[g], g);
      jet->GetShapeProperties()->AddGRDenSubAt(dens[g], g);
    }
  }

  if (fDoGenericSubtractionExtraJetShapes) {
    std::vector<fastjet::contrib::GenericSubtractorInfo> jetAngularityInfo = fjw.GetGenSubtractorInfoJetAngularity();
    Int_t na = (Int_t)jetAngularityInfo.size();
    if(na > ij && na > 0) {
      jet->GetShapeProperties()->SetFirstDerivativeAngularity(jetAngularityInfo[ij].first_derivative());
      jet->GetShapeProperties()->SetSecondDerivativeAngularity(jetAngularityInfo[ij].second_derivative());
      jet->GetShapeProperties()->SetFirstOrderSubtractedAngularity(jetAngularityInfo[ij].first_order_subtracted());
      jet->GetShapeProperties()->SetSecondOrderSubtractedAngularity(jetAngularityInfo[ij].second_order_subtracted());
    }

    std::vector<fastjet::contrib::GenericSubtractorInfo> jetpTDInfo = fjw.GetGenSubtractorInfoJetpTD();
    Int_t np = (Int_t)jetpTDInfo.size();
    if(np > ij && np > 0) {
      jet->GetShapeProperties()->SetFirstDerivativepTD(jetpTDInfo[ij].first_derivative());
      jet->GetShapeProperties()->SetSecondDerivativepTD(jetpTDInfo[ij].second_derivative());
      jet->GetShapeProperties()->SetFirstOrderSubtractedpTD(jetpTDInfo[ij].first_order_subtracted());
      jet->GetShapeProperties()->SetSecondOrderSubtractedpTD(jetpTDInfo[ij].second_order_subtracted());
    }

    std::vector<fastjet::contrib::GenericSubtractorInfo> jetCircularityInfo = fjw.GetGenSubtractorInfoJetCircularity();
    Int_t nc = (Int_t)jetCircularityInfo.size();
    if(nc > ij && nc > 0) {
      jet->GetShapeProperties()->SetFirstDerivativeCircularity(jetCircularityInfo[ij].first_derivative());
      jet->GetShapeProperties()->SetSecondDerivativeCircularity(jetCircularityInfo[ij].second_derivative());
      jet->GetShapeProperties()->SetFirstOrderSubtractedCircularity(jetCircularityInfo[ij].first_order_subtracted());
      jet->GetShapeProperties()->SetSecondOrderSubtractedCircularity(jetCircularityInfo[ij].second_order_subtracted());
    }

    std::vector<fastjet::contrib::GenericSubtractorInfo> jetSigma2Info = fjw.GetGenSubtractorInfoJetSigma2();
    Int_t ns = (Int_t)jetSigma2Info.size();
    if (ns > ij && ns > 0) {
      jet->GetShapeProperties()->SetFirstDerivativeSigma2(jetSigma2Info[ij].first_derivative());
      jet->GetShapeProperties()->SetSecondDerivativeSigma2(jetSigma2Info[ij].second_derivative());
      jet->GetShapeProperties()->SetFirstOrderSubtractedSigma2(jetSigma2Info[ij].first_order_subtracted());
      jet->GetShapeProperties()->SetSecondOrderSubtractedSigma2(jetSigma2Info[ij].second_order_subtracted());
    }


    std::vector<fastjet::contrib::GenericSubtractorInfo> jetConstituentInfo = fjw.GetGenSubtractorInfoJetConstituent();
    Int_t nco = (Int_t)jetConstituentInfo.size();
    if(nco > ij && nco > 0) {
      jet->GetShapeProperties()->SetFirstDerivativeConstituent(jetConstituentInfo[ij].first_derivative());
      jet->GetShapeProperties()->SetSecondDerivativeConstituent(jetConstituentInfo[ij].second_derivative());
      jet->GetShapeProperties()->SetFirstOrderSubtractedConstituent(jetConstituentInfo[ij].first_order_subtracted());
      jet->GetShapeProperties()->SetSecondOrderSubtractedConstituent(jetConstituentInfo[ij].second_order_subtracted());
    }
    
    std::vector<fastjet::contrib::GenericSubtractorInfo> jetLeSubInfo = fjw.GetGenSubtractorInfoJetLeSub();
    Int_t nlsub = (Int_t)jetLeSubInfo.size();
    if(nlsub > ij && nlsub > 0) {
      jet->GetShapeProperties()->SetFirstDerivativeLeSub(jetLeSubInfo[ij].first_derivative());
      jet->GetShapeProperties()->SetSecondDerivativeLeSub(jetLeSubInfo[ij].second_derivative());
      jet->GetShapeProperties()->SetFirstOrderSubtractedLeSub(jetLeSubInfo[ij].first_order_subtracted());
      jet->GetShapeProperties()->SetSecondOrderSubtractedLeSub(jetLeSubInfo[ij].second_order_subtracted());
    }
  }

  if (fDoGenericSubtractionNsubjettiness) {
    std::vector<fastjet::contrib::GenericSubtractorInfo> jet1subjettinessktInfo = fjw.GetGenSubtractorInfoJet1subjettiness_kt();
    Int_t n1subjettiness_kt = (Int_t)jet1subjettinessktInfo.size();
    if(n1subjettiness_kt > ij && n1subjettiness_kt > 0) {
      jet->GetShapeProperties()->SetFirstDerivative1subjettiness_kt(jet1subjettinessktInfo[ij].first_derivative());
      jet->GetShapeProperties()->SetSecondDerivative1subjettiness_kt(jet1subjettinessktInfo[ij].second_derivative());
      jet->GetShapeProperties()->SetFirstOrderSubtracted1subjettiness_kt(jet1subjettinessktInfo[ij].first_order_subtracted());
      jet->GetShapeProperties()->SetSecondOrderSubtracted1subjettiness_kt(jet1subjettinessktInfo[ij].second_order_subtracted());
    }
          
    std::vector<fastjet::contrib::GenericSubtractorInfo> jet2subjettinessktInfo = fjw.GetGenSubtractorInfoJet2subjettiness_kt();
    Int_t n2subjettiness_kt = (Int_t)jet2subjettinessktInfo.size();
    if(n2subjettiness_kt > ij && n2subjettiness_kt > 0) {
      jet->GetShapeProperties()->SetFirstDerivative2subjettiness_kt(jet2subjettinessktInfo[ij].first_derivative());
      jet->GetShapeProperties()->SetSecondDerivative2subjettiness_kt(jet2subjettinessktInfo[ij].second_derivative());
      jet->GetShapeProperties()->SetFirstOrderSubtracted2subjettiness_kt(jet2subjettinessktInfo[ij].first_order_subtracted());
      jet->GetShapeProperties()->SetSecondOrderSubtracted2subjettiness_kt(jet2subjettinessktInfo[ij].second_order_subtracted());
    }

    std::vector<fastjet::contrib::GenericSubtractorInfo> jet3subjettinessktInfo = fjw.GetGenSubtractorInfoJet3subjettiness_kt();
    Int_t n3subjettiness_kt = (Int_t)jet3subjettinessktInfo.size();
    if(n3subjettiness_kt > ij && n3subjettiness_kt > 0) {
      jet->GetShapeProperties()->SetFirstDerivative3subjettiness_kt(jet3subjettinessktInfo[ij].first_derivative());
      jet->GetShapeProperties()->SetSecondDerivative3subjettiness_kt(jet3subjettinessktInfo[ij].second_derivative());
      jet->GetShapeProperties()->SetFirstOrderSubtracted3subjettiness_kt(jet3subjettinessktInfo[ij].first_order_subtracted());
      jet->GetShapeProperties()->SetSecondOrderSubtracted3subjettiness_kt(jet3subjettinessktInfo[ij].second_order_subtracted());
    }

    std::vector<fastjet::contrib::GenericSubtractorInfo> jetOpeningAnglektInfo = fjw.GetGenSubtractorInfoJetOpeningAngle_kt();
    Int_t nOpeningAngle_kt = (Int_t)jetOpeningAnglektInfo.size();
    if(nOpeningAngle_kt > ij && nOpeningAngle_kt > 0) {
      jet->GetShapeProperties()->SetFirstDerivativeOpeningAngle_kt(jetOpeningAnglektInfo[ij].first_derivative());
      jet->GetShapeProperties()->SetSecondDerivativeOpeningAngle_kt(jetOpeningAnglektInfo[ij].second_derivative());
      jet->GetShapeProperties()->SetFirstOrderSubtractedOpeningAngle_kt(jetOpeningAnglektInfo[ij].first_order_subtracted());
      jet->GetShapeProperties()->SetSecondOrderSubtractedOpeningAngle_kt(jetOpeningAnglektInfo[ij].second_order_subtracted());
    }
  }

#endif
}

//______________________________________________________________________________
void AliEmcalJetUtilityGenSubtractor::Terminate(AliFJWrapper& /*fjw*/)
{
  // Run termination of the utility (after each event).
}
