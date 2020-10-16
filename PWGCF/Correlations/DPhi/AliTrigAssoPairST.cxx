#include "AliTrigAssoPairST.h"

ClassImp(AliTrigAssoPairST)

//_____________________________________________________________________________
AliTrigAssoPairST::AliTrigAssoPairST(Short_t charge, Float_t eta, Float_t phi, Float_t pt, Float_t pt2,
                       Int_t ID, Int_t ID1, Int_t ID2, Short_t candidate,
                       Double_t multiplicity, Double_t deta_pairs, Double_t dphi_pairs) : 
fCharge(charge),
fEta(eta),
fPhi(phi),
fpT(pt),
fpT_Asso(pt2),
fID(ID),
fID1(ID1),
fID2(ID2),
fCandidate(candidate),
fMultiplicity(multiplicity),
fdeta_pairs(deta_pairs),
fdphi_pairs(dphi_pairs) 
{
//
// AliTrigAssoPairST::AliTrigAssoPairST(Short_t charge, Float_t eta, Float_t phi, Float_t pt, Float_t pt2, Int_t ID, Int_t ID1, Int_t ID2, Short_t candidate, Double_t multiplicity, Double_t deta_pairs, Double_t dphi_pairs) :
//
}

AliTrigAssoPairST::~AliTrigAssoPairST() 
{
//
// AliTrigAssoPairST::~AliTrigAssoPairST()
//
}

