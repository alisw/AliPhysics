
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliFast Track class                                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "AliFTrack.h"
#include "AliFast.h"

ClassImp(AliFTrack)

//_____________________________________________________________________________
AliFTrack::AliFTrack(Int_t code, Double_t charge,
                     Double_t pT, Double_t eta, Double_t phi,
                     Double_t v11, Double_t v22,Double_t v33,
                     Double_t v12, Double_t v13, Double_t v23, Int_t iFlag)
{
   fIdTrack    = code;
   fChTrack    = charge;
   fPT         = pT;
   fEta        = eta;
   fPhi        = phi;
   fV11        = v11;
   fV22        = v22;
   fV33        = v33;
   fV13        = v13;
   fV12        = v12;
   fV23        = v23;
   fIFlag      = iFlag;
}


//_____________________________________________________________________________
void AliFTrack::Draw(Option_t *)
{

}

//_____________________________________________________________________________
void AliFTrack::Paint(Option_t *)
{

}

