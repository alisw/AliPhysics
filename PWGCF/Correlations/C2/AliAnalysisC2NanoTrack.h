/// \class AliAnalysisC2NanoTrack
/// \brief A small class to hold information of tracks. This makes it easier to add or remove
/// tracks for debugging purposes.
///
/// \author Christian Bourjau <cbourjau@cern.ch>, University of Copenhagen, Denmark
/// \date Jul 18, 2016

#ifndef AliAnalysisC2NanoTrack_cxx
#define AliAnalysisC2NanoTrack_cxx

#include "TObject.h"

class AliAnalysisC2NanoTrack : public TObject {
 public:
  Double_t eta;
  Double_t phi;
  Double_t pt;
  // Default constructor
  AliAnalysisC2NanoTrack();
  AliAnalysisC2NanoTrack(Double_t eta_in, Double_t phi_in, Double_t pt_in);
  virtual ~AliAnalysisC2NanoTrack() {};
 private:
  ClassDef(AliAnalysisC2NanoTrack,1);
};

#endif
