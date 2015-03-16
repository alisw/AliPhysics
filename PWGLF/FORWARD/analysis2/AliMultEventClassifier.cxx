#include <TH1.h>
#include <TList.h>
#include "AliMultEventClassifier.h"
#include "AliAODMultEventClass.h"
#include "AliPPVsMultUtils.h"
#include <AliESDEvent.h>
#include <TArrayD.h>
#include <TH2.h>
#include <AliCentrality.h>
#include <AliESDtrackCuts.h>

//____________________________________________________________________
const char*
AliMultEventClassifier::GetCentName(UShort_t which) const
{
  Bool_t isEq = (which & AliAODMultEventClass::kEq);
  if      (which & AliAODMultEventClass::kV0M) return isEq ? "V0MEq" : "V0M";
  else if (which & AliAODMultEventClass::kV0A) return isEq ? "V0AEq" : "V0A";
  else if (which & AliAODMultEventClass::kV0C) return isEq ? "V0CEq" : "V0C";
  return 0;
}

//____________________________________________________________________
TH2*
AliMultEventClassifier::GetCorr(UShort_t which) const
{
  Bool_t isEq = (which & AliAODMultEventClass::kEq);
  if      (which & AliAODMultEventClass::kV0M) 
    return isEq ? fCorrV0MEq : fCorrV0M;
  else if (which & AliAODMultEventClass::kV0A) 
    return isEq ? fCorrV0AEq : fCorrV0A;
  else if (which & AliAODMultEventClass::kV0C) 
    return isEq ? fCorrV0CEq : fCorrV0C;
  return 0;
}  
//____________________________________________________________________
TH2*
AliMultEventClassifier::GetVs(UShort_t which) const
{
  Bool_t isEq = (which & AliAODMultEventClass::kEq);
  if      (which & AliAODMultEventClass::kV0M) 
    return isEq ? fMultV0MEq : fMultV0M;
  else if (which & AliAODMultEventClass::kV0A) 
    return isEq ? fMultV0AEq : fMultV0A;
  else if (which & AliAODMultEventClass::kV0C) 
    return isEq ? fMultV0CEq : fMultV0C;
  return 0;
}  
//____________________________________________________________________
void
AliMultEventClassifier::GetCentrality(AliESDEvent* esd,
				      AliAODMultEventClass* data,
				      Int_t mult,
				      UShort_t which)
{
  const char*    meth = GetCentName(which);
  Float_t        util = fUtil->GetMultiplicityPercentile(esd,meth);
  Float_t        sel  = -1;
  AliCentrality* cObj = esd->GetCentrality();
  if (cObj)      sel  = cObj->GetCentralityPercentile(meth);  
  TH2*           corr = GetCorr(which);
  TH2*           vs   = GetVs(which);

  if (corr) corr->Fill(util, sel);
  if (vs)   vs  ->Fill(mult, util);

  if (data) return;
  data->SetCentrality(which, true, util);
  data->SetCentrality(which, false, sel);
  
}
//____________________________________________________________________
TH2* 
AliMultEventClassifier::MakeCorr(UShort_t which)
{
  const char*  meth = GetCentName(which);
  TH2*         corr = new TH2D(Form("corr%s", meth),
			       Form("Correlation of %s estimators", meth),
			       101, -1.5, 100.5, 101, -1.5, 100.5);
  corr->SetDirectory(0);
  corr->SetXTitle("Centrality from AliPPVsMultUtils [%]");
  corr->SetYTitle("Centrality from AliCentralitySelector [%]");

  Bool_t isEq = (which & AliAODMultEventClass::kEq);
  if      (which & AliAODMultEventClass::kV0M)
    if (isEq) fCorrV0MEq = corr; else fCorrV0M = corr;
  else if (which & AliAODMultEventClass::kV0A) 
    if (isEq) fCorrV0AEq = corr; else fCorrV0A = corr;
  else if (which & AliAODMultEventClass::kV0C) 
    if (isEq) fCorrV0CEq = corr; else fCorrV0C = corr;

  fList->Add(corr);
  return corr;
}
  
//____________________________________________________________________
TH2* 
AliMultEventClassifier::MakeVs(UShort_t which, const TArrayD& bins)
{
  const char*  meth = GetCentName(which);
  TH2*         vs   = new TH2D(Form("vs%s", meth),
			       Form("Refernce multiplicity vs %s estimator",
				    meth),
			       bins.GetSize()-1,bins.GetArray(),
			       101, -1.5, 100.5);
  vs->SetDirectory(0);
  vs->SetXTitle("Reference multiplicity");
  vs->SetYTitle("Centrality from AliPPVsMultUtils [%]");

  Bool_t isEq = (which & AliAODMultEventClass::kEq);
  if      (which & AliAODMultEventClass::kV0M)
    if (isEq) fMultV0MEq = vs; else fMultV0M = vs;
  else if (which & AliAODMultEventClass::kV0A) 
    if (isEq) fMultV0AEq = vs; else fMultV0A = vs;
  else if (which & AliAODMultEventClass::kV0C) 
    if (isEq) fMultV0CEq = vs; else fMultV0C = vs;
  

  fList->Add(vs);
  return vs;
}
  
//____________________________________________________________________
void
AliMultEventClassifier::CreateOutputObjects(TList* l)
{
  fList = new TList;
  fList->SetOwner(true);
  if (l) l->Add(fList);

  const Int_t* arr   = AliAODMultEventClass::GetBins();
  Int_t*       tmp   = const_cast<Int_t*>(arr);
  Int_t        n     = 0;
  while ((*tmp >= 0)) { n++; tmp++; }

  TArrayD bins(n+2);
  bins[0] = -1;
  Int_t i = 1;
  tmp     = const_cast<Int_t*>(arr);
  while ((*tmp >= 0)) {
    bins[i] = *tmp;
    tmp++;
    i++;
  }
  fMax    = bins[i-1];
  bins[i] = fMax + 10;

  MakeCorr(AliAODMultEventClass::kV0M);
  MakeCorr(AliAODMultEventClass::kV0A);
  MakeCorr(AliAODMultEventClass::kV0C);
  MakeCorr(AliAODMultEventClass::kV0M|AliAODMultEventClass::kEq);
  MakeCorr(AliAODMultEventClass::kV0A|AliAODMultEventClass::kEq);
  MakeCorr(AliAODMultEventClass::kV0C|AliAODMultEventClass::kEq);
  MakeVs(AliAODMultEventClass::kV0M, bins);
  MakeVs(AliAODMultEventClass::kV0A, bins);
  MakeVs(AliAODMultEventClass::kV0C, bins);
  MakeVs(AliAODMultEventClass::kV0M|AliAODMultEventClass::kEq, bins);
  MakeVs(AliAODMultEventClass::kV0A|AliAODMultEventClass::kEq, bins);
  MakeVs(AliAODMultEventClass::kV0C|AliAODMultEventClass::kEq, bins);

  fUtil = new AliPPVsMultUtils;
}
//____________________________________________________________________
void
AliMultEventClassifier::Process(AliESDEvent* esd,
				AliAODMultEventClass* data)
{
  if (data) data->Clear();
  Int_t mult =
    AliESDtrackCuts
    ::GetReferenceMultiplicity(esd,
			       AliESDtrackCuts::kTrackletsITSTPC, 0.8 ) ;
  if (data) data->SetMult(mult);
  Int_t fill = (mult >= fMax ? fMax+1 : mult);

  GetCentrality(esd, data, fill, AliAODMultEventClass::kV0M);
  GetCentrality(esd, data, fill, AliAODMultEventClass::kV0A);
  GetCentrality(esd, data, fill, AliAODMultEventClass::kV0C);
  GetCentrality(esd, data, fill, AliAODMultEventClass::kV0M|
		AliAODMultEventClass::kEq);
  GetCentrality(esd, data, fill, AliAODMultEventClass::kV0A|
		AliAODMultEventClass::kEq);
  GetCentrality(esd, data, fill, AliAODMultEventClass::kV0C|
		AliAODMultEventClass::kEq);
}

//____________________________________________________________________
//
// EOF
//

  
  
