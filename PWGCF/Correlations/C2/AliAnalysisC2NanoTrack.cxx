#include "TObject.h"

#include "AliAnalysisC2NanoTrack.h"

AliAnalysisC2NanoTrack::AliAnalysisC2NanoTrack()
: TObject(),
  eta(0),
  phi(0),
  pt(0),
  weight(0)
{
}

AliAnalysisC2NanoTrack::AliAnalysisC2NanoTrack(Double_t eta_in, Double_t phi_in,
					       Double_t pt_in, Double_t weight)
  : TObject(),
    eta(eta_in),
    phi(phi_in),
    pt(pt_in),
    weight(weight)
{

}
