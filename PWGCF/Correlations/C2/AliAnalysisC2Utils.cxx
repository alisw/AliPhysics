#include <iostream>

#include "TAxis.h"
#include "THn.h"

#include "AliAODTrack.h"
#include "AliVEvent.h"
#include "AliAnalysisC2Utils.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisC2Utils)

/// Modulo for float numbers
///
/// \param x nominator
/// \param y denominator
///
/// \return Rest of the rounded down division
Double_t AliAnalysisC2Utils::Mod(Double_t x, Double_t y) {
  if (0 == y)
    return x;
  return x - y * floor(x/y);
}

/// Wrap angle around 0 and 2pi
Double_t AliAnalysisC2Utils::Wrap02pi(Double_t angle) {
  const Double_t two_pi = 6.283185307179586;
  Double_t lower_edge = 0;
  Double_t interval = two_pi;
  if (lower_edge <= angle && angle < two_pi) {
    return angle;
  }
  return Mod(angle - lower_edge, interval) + lower_edge;
}

/// Set the range of each axis such that it excludes over and under flow
Int_t AliAnalysisC2Utils::ExcludeOverUnderFlow(THn *h) {
  for (Int_t iaxis = 0; iaxis < h->GetNdimensions(); iaxis++) {
    TAxis* axis = h->GetAxis(iaxis);
    axis->SetRange(1, axis->GetNbins());
    axis->SetBit(TAxis::kAxisRange);
  }
  return 0;
}

/// Copy the properties of one histogram's axis to the axis of another histogram
void AliAnalysisC2Utils::CopyAxisFromHist(Int_t from_axis, THn *from_hist, Int_t to_axis, THn *to_hist){
  /*
        Copy the properties of one axis from one hist to another
  */
  TAxis *ax_from = from_hist->GetAxis(from_axis);
  TAxis *ax_to = to_hist->GetAxis(to_axis);
  Double_t edges[ax_from->GetNbins() + 1]; 
  for (Int_t iedge = 0; iedge < ax_from->GetNbins(); iedge++){
    // Beware of the ROOT trap! Bin numbering starts at 1!
    // highest value of iedge in the loop is GetNbins() - 1
    edges[iedge] = ax_from->GetBinLowEdge(iedge + 1);
  }
  edges[ax_from->GetNbins()] = ax_from->GetBinUpEdge(ax_from->GetNbins());
  ax_to->Set(ax_from->GetNbins(), edges);
  ax_to->SetTitle(ax_from->GetTitle());
}

/// Transform the given array to the "sum and difference" coordinates
void AliAnalysisC2Utils::TransformPoints(Double_t *points_in, Double_t *points_out){
  // transform the given coordinates expecting array of dim 6
  points_out[0] = points_in[0] - points_in[1];
  points_out[1] = (points_in[0] + points_in[1]) / 2;
  points_out[2] = points_in[2] - points_in[3];
  points_out[3] = (points_in[2] + points_in[3]) / 2;
  points_out[4] = points_in[4];
  points_out[5] = points_in[5];
}


/// Create a histogram with corresponding settings as the given one, but in the "sum and differences coordinates"
THnS* AliAnalysisC2Utils::CreateTransformedHist(THn *h){
  Int_t eta_bar_nbins = h->GetAxis(0)->GetNbins() * 2 + 1;
  Double_t eta_bar_binwidth = h->GetAxis(0)->GetBinWidth(1) / 2.;
  Double_t eta_bar_min = ((h->GetAxis(0)->GetBinLowEdge(1) + h->GetAxis(1)->GetBinLowEdge(1)) / 2 - eta_bar_binwidth / 2.);
  Double_t eta_bar_max = ((h->GetAxis(0)->GetBinUpEdge(h->GetAxis(0)->GetNbins()) +
			  h->GetAxis(1)->GetBinUpEdge(h->GetAxis(1)->GetNbins())) / 2 + eta_bar_binwidth / 2.);

  Int_t delta_eta_nbins = h->GetAxis(0)->GetNbins() * 2 - 1;
  Double_t delta_eta_binwidth = h->GetAxis(0)->GetBinWidth(1);
  Double_t delta_eta_min = (h->GetAxis(0)->GetBinLowEdge(1) * 2 + delta_eta_binwidth / 2.);
  Double_t delta_eta_max = (h->GetAxis(0)->GetBinUpEdge(h->GetAxis(0)->GetNbins()) * 2 - delta_eta_binwidth / 2.);

  Int_t phi_bar_nbins = h->GetAxis(0)->GetNbins();
  Double_t phi_bar_min = (h->GetAxis(2)->GetBinLowEdge(1) + h->GetAxis(3)->GetBinLowEdge(1)) / 2;
  Double_t phi_bar_max = (h->GetAxis(2)->GetBinUpEdge(h->GetAxis(2)->GetNbins()) +
			 h->GetAxis(3)->GetBinUpEdge(h->GetAxis(3)->GetNbins())) / 2;
  if ((h->GetAxis(2)->GetNbins() != h->GetAxis(3)->GetNbins()) ||
      (h->GetAxis(2)->GetNbins() % 2 != 0) ||
      (h->GetAxis(2)->GetNbins() % 4 == 0))  // must be divisible by 2 but not by 4!
    {
      cout << "Number of phi bins is incompatible for shift by pi/2!" << endl;
    }
  Double_t delta_phi_binwidth = h->GetAxis(2)->GetBinWidth(1);
  Int_t delta_phi_nbins = h->GetAxis(2)->GetNbins();
  // we want to shift by pi/2 which is 1/4th of the full interval.
  Double_t delta_phi_offset = h->GetAxis(2)->GetNbins() / 4. * delta_phi_binwidth;
  Double_t delta_phi_min = h->GetAxis(2)->GetBinLowEdge(1) - delta_phi_offset;
  Double_t delta_phi_max = h->GetAxis(2)->GetBinUpEdge(h->GetAxis(2)->GetNbins()) - delta_phi_offset;

  const Int_t nbins[] = {delta_eta_nbins,
			 eta_bar_nbins,
			 delta_phi_nbins,
			 phi_bar_nbins,
			 h->GetAxis(4)->GetNbins(),         // mult
			 h->GetAxis(5)->GetNbins()};       // zvtx
  const Double_t xmin[] = {delta_eta_min,
			  eta_bar_min,
			  delta_phi_min,
			  phi_bar_min,
			  0,         // mult
			  0};       // zvtx

  const Double_t xmax[] = {delta_eta_max,
			  eta_bar_max,
			  delta_phi_max,
			  phi_bar_max,
			  1,
			  1};

  THnS *hist_out = new THnS(TString(h->GetName()) + TString("_trans"),
			    TString(h->GetTitle())
			    + TString(";#Delta#eta;#bar{#eta};#Delta#phi;#bar{#phi};mult;z_{vtx};"),
			    h->GetNdimensions(),
			    nbins,
			    xmin,
			    xmax);
  // Copy the axis for mult and zvtx
  AliAnalysisC2Utils::CopyAxisFromHist(4, h, 4, hist_out);  // mult
  AliAnalysisC2Utils::CopyAxisFromHist(5, h, 5, hist_out);  // zvtx
  return hist_out;
}

void AliAnalysisC2Utils::GetDCA(Double_t& DCAtang, Double_t& DCAlong, AliAODTrack* AODt){
  DCAtang = -999;
  DCAlong=-999;
  if(AODt->TestBit(AliAODTrack::kIsDCA)) {
    DCAtang = AODt->DCA();
    DCAlong = AODt->ZAtDCA();
  }
  else if(AODt->GetEvent()->GetPrimaryVertex()) {
    Double_t dca[2];
    Double_t dcacov[3];
    Double_t maxDistance = 1000000;
    // Warning: This might change some internal parameters of the track!
    // Might be safer to clone the track before Propagation!
    Bool_t propagatedSuccessfully = AODt->PropagateToDCA(AODt->GetEvent()->GetPrimaryVertex(),
							 AODt->GetEvent()->GetMagneticField(),
							 maxDistance, dca, dcacov);
    if(propagatedSuccessfully){
      DCAtang=dca[0];
      DCAlong = dca[1];
    }
  }
}


Int_t AliAnalysisC2Utils::ComputePtPairBin(Int_t pt1Bin, Int_t pt2Bin){
  if (pt1Bin <= 0 || pt2Bin <= 0) {
    return -1;
  }
  Int_t idx = 0;
  for (Int_t i = 1; i <= pt2Bin; i++){
    idx += i;
  }
  idx += (pt1Bin - pt2Bin);
  idx -= 1;
  return idx;
}

Bool_t AliAnalysisC2Utils::EventFitsTrigger(UInt_t trig)
{
  UInt_t maskIsSelected =
    dynamic_cast<AliInputEventHandler*>
    (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())
    ->IsEventSelected();
  if (maskIsSelected & trig) {
    return true;
  } else {
    return false;
  }
}
