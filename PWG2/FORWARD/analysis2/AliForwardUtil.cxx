#include "AliForwardUtil.h"
#include <TH2D.h>
#include <TMath.h>

//====================================================================
AliForwardUtil::Histos::~Histos()
{
  if (fFMD1i) delete fFMD1i;
  if (fFMD2i) delete fFMD2i;
  if (fFMD2o) delete fFMD2o;
  if (fFMD3i) delete fFMD3i;
  if (fFMD3o) delete fFMD3o;
}

//____________________________________________________________________
TH2D*
AliForwardUtil::Histos::Make(UShort_t d, Char_t r, 
			     const TAxis& etaAxis) const
{
  Int_t ns = (r == 'I' || r == 'i') ? 20 : 40;
  TH2D* hist = new TH2D(Form("FMD%d%c_cache", d, r), 
			Form("FMD%d%c cache", d, r),
			etaAxis.GetNbins(), etaAxis.GetXmin(), 
			etaAxis.GetXmax(), ns, 0, 2*TMath::Pi());
  hist->SetXTitle("#eta");
  hist->SetYTitle("#phi [radians]");
  hist->SetZTitle("d^{2}N_{ch}/d#etad#phi");
  hist->Sumw2();
  hist->SetDirectory(0);

  return hist;
}
//____________________________________________________________________
void
AliForwardUtil::Histos::Init(const TAxis& etaAxis)
{
  fFMD1i = Make(1, 'I', etaAxis);
  fFMD2i = Make(2, 'I', etaAxis);
  fFMD2o = Make(2, 'O', etaAxis);
  fFMD3i = Make(3, 'I', etaAxis);
  fFMD3o = Make(3, 'O', etaAxis);
}
//____________________________________________________________________
void
AliForwardUtil::Histos::Clear(Option_t* option)
{
  fFMD1i->Reset(option);
  fFMD2i->Reset(option);
  fFMD2o->Reset(option);
  fFMD3i->Reset(option);
  fFMD3o->Reset(option);
}

//____________________________________________________________________
TH2D*
AliForwardUtil::Histos::Get(UShort_t d, Char_t r) const
{
  switch (d) { 
  case 1: return fFMD1i;
  case 2: return (r == 'I' || r == 'i' ? fFMD2i : fFMD2o);
  case 3: return (r == 'I' || r == 'i' ? fFMD3i : fFMD3o);
  }
  return 0;
}
#if 0
//____________________________________________________________________
TH2D*
AliForwardUtil::Histos::Get(UShort_t d, Char_t r)
{
  switch (d) { 
  case 1: return fFMD1i;
  case 2: return (r == 'I' || r == 'i' ? fFMD2i : fFMD2o);
  case 3: return (r == 'I' || r == 'i' ? fFMD3i : fFMD3o);
  }
  default: return 0;
}
#endif

//
// EOF
//
