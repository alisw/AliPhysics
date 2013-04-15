//
// Author: Francesca Bellini (fbellini@cern.ch)
//
//

#include <Riostream.h>

#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliRsnCutSetDaughterParticle.h"

class AliRsnCutPIDNSigma;
class AliRsnPIDRange;
class AliRsnCutPhi;

ClassImp(AliRsnCutSetDaughterParticle)

//__________________________________________________________________________________________________
AliRsnCutSetDaughterParticle::AliRsnCutSetDaughterParticle() :
   AliRsnCutSet("AliRsnCutSetDaughterParticle", AliRsnTarget::kDaughter),
   fPID(AliPID::kPion),
   fAppliedCutSetID(AliRsnCutSetDaughterParticle::kNDaughterCuts),
   fNsigmaTPC(1E20),
   fNsigmaTOF(1E20),
   fCutQuality(0x0),
   fAODTrkCutFilterBit(0)
{
   //
   // Default constructor
}

//__________________________________________________________________________________________________
AliRsnCutSetDaughterParticle::AliRsnCutSetDaughterParticle(const char *name, AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutSetID, AliPID::EParticleType pid, Float_t nSigmaFast = -1.0, Int_t AODfilterBit = 0) :
   AliRsnCutSet(name, AliRsnTarget::kDaughter),
   fPID(pid),
   fAppliedCutSetID(cutSetID),
   fNsigmaTPC(1E20),
   fNsigmaTOF(1E20),
   fCutQuality(new AliRsnCutTrackQuality("CutQuality")),
   fAODTrkCutFilterBit(AODfilterBit)
{
   //
   // Constructor
   //
   if ( (nSigmaFast<=0) &&
        ((cutSetID == AliRsnCutSetDaughterParticle::kFastTPCpidNsigma) || (cutSetID == AliRsnCutSetDaughterParticle::kFastTOFpidNsigma) || (cutSetID == AliRsnCutSetDaughterParticle::kTOFMatchTPCpidNsigma)) ) {
      AliError("Requested fast n-sigma PID with invalid value for n. Setting n = 1E20");
   } else {
      if (cutSetID == AliRsnCutSetDaughterParticle::kFastTPCpidNsigma) {
         fNsigmaTPC = nSigmaFast;
      }
      if (cutSetID == AliRsnCutSetDaughterParticle::kTOFMatchTPCpidNsigma) {
	fNsigmaTPC = nSigmaFast;
      }
      if ( (cutSetID == AliRsnCutSetDaughterParticle::kFastTOFpidNsigma) ||
           (cutSetID == AliRsnCutSetDaughterParticle::kTOFpidKstarPbPb2010) ) {
         fNsigmaTOF = nSigmaFast;
      }
   }

   Init();
}

//__________________________________________________________________________________________________
AliRsnCutSetDaughterParticle::AliRsnCutSetDaughterParticle(const AliRsnCutSetDaughterParticle &copy) :
   AliRsnCutSet(copy),
   fPID(copy.fPID),
   fAppliedCutSetID(copy.fAppliedCutSetID),
   fNsigmaTPC(copy.fNsigmaTPC),
   fNsigmaTOF(copy.fNsigmaTOF),
   fCutQuality(copy.fCutQuality),
   fAODTrkCutFilterBit(copy.fAODTrkCutFilterBit)
{
   //
   // copy constructor
}

//__________________________________________________________________________________________________
AliRsnCutSetDaughterParticle &AliRsnCutSetDaughterParticle::operator=(const AliRsnCutSetDaughterParticle &copy)
{
   //
   // operator =
   //
   AliRsnCutSet::operator=(copy);
   if (this == &copy)
      return *this;
   fPID=copy.fPID;
   fAppliedCutSetID=copy.fAppliedCutSetID;
   fNsigmaTPC=copy.fNsigmaTPC;
   fNsigmaTOF=copy.fNsigmaTOF;
   fAODTrkCutFilterBit=copy.fAODTrkCutFilterBit;
   fCutQuality=copy.fCutQuality;
   return (*this);
}

//----------------------------------------------------------------------------
AliRsnCutSetDaughterParticle::~AliRsnCutSetDaughterParticle()
{
   //
   //destructor
   //
   if (fCutQuality)
      delete fCutQuality;
}
//----------------------------------------------------------------------------
void AliRsnCutSetDaughterParticle::Init()
{
   //
   // init cut sets by setting variable params
   //
   fCutQuality->SetDefaults2010();
   fCutQuality->SetPtRange(0.15, 20.0);
   fCutQuality->SetEtaRange(-0.8, 0.8);
   fCutQuality->SetDCARPtFormula("0.0182+0.0350/pt^1.01");
   fCutQuality->SetDCAZmax(2.0);
   fCutQuality->SetSPDminNClusters(1);
   fCutQuality->SetITSminNClusters(0);
   fCutQuality->SetITSmaxChi2(36);
   fCutQuality->SetTPCminNClusters(70);
   fCutQuality->SetTPCmaxChi2(4.0);
   fCutQuality->SetRejectKinkDaughters();
   fCutQuality->SetAODTestFilterBit(fAODTrkCutFilterBit);
   //fCutQuality->SetITSmaxChi2(36);
   //fCutQuality->SetMaxChi2TPCConstrainedGlobal(36);

   AliRsnCutTOFMatch  *iCutTOFMatch     = new AliRsnCutTOFMatch("CutTOFMatch");
   AliRsnCutPIDNSigma *iCutTPCNSigma    = new AliRsnCutPIDNSigma("CutTPCNSigma", fPID, AliRsnCutPIDNSigma::kTPC);//, AliRsnCutPIDNSigma::kTPCinnerP );
   AliRsnCutPIDNSigma *iCutTPCTOFNSigma = new AliRsnCutPIDNSigma("CutTPCTOFNSigma", fPID, AliRsnCutPIDNSigma::kTPC);//, AliRsnCutPIDNSigma::kTPCinnerP );
   AliRsnCutPIDNSigma *iCutTOFNSigma    = new AliRsnCutPIDNSigma("CutTOFNSigma", fPID, AliRsnCutPIDNSigma::kTOF);//, AliRsnCutPIDNSigma::kP );
   AliRsnCutPhi  *iCutPhiTRD2010        = new AliRsnCutPhi("CutPhiTRD2010","InTRD");
   AliRsnCutPhi  *iCutPhiNoTRD2010      = new AliRsnCutPhi("CutPhiNoTRD2010","OutTRD");
   
   switch (fAppliedCutSetID)
   {
      case AliRsnCutSetDaughterParticle::kNoCuts :
         AliInfo("No cuts applied to daughter particle");
         break;

      case AliRsnCutSetDaughterParticle::kQualityStd2010 :
         AddCut(fCutQuality);
         SetCutScheme(fCutQuality->GetName());
         break;

      case AliRsnCutSetDaughterParticle::kQualityStd2011:
         //fCutQuality->SetAODTestFilterBit(10);     //1024
         AddCut(fCutQuality);
         SetCutScheme(fCutQuality->GetName());
         break;

      case AliRsnCutSetDaughterParticle::kTOFMatch :
         AddCut(fCutQuality);
         AddCut(iCutTOFMatch);
         SetCutScheme( Form("%s&(%s)",fCutQuality->GetName(), iCutTOFMatch->GetName()) );
         break;

      case    AliRsnCutSetDaughterParticle::kFastTPCpidNsigma :
         if (fNsigmaTPC <= 0.0) {
            AliWarning(Form("Invalid number of sigmas required for %s. Setting default nSigma = 10",iCutTPCNSigma->GetName()));
            SetNsigmaForFastTPCpid(10.0);
         }
         iCutTPCNSigma->SinglePIDRange(fNsigmaTPC);
         AddCut(fCutQuality);
         AddCut(iCutTPCNSigma);
         SetCutScheme( Form("%s&%s",fCutQuality->GetName(), iCutTPCNSigma->GetName()) );
         break;

      case    AliRsnCutSetDaughterParticle::kFastTOFpidNsigma :
         if (fNsigmaTOF <= 0.0) {
            AliWarning(Form("Invalid number of sigmas required for %s. Setting default nSigma = 10",iCutTOFNSigma->GetName()));
            SetNsigmaForFastTOFpid(10.0);
         }
         iCutTOFNSigma->SinglePIDRange(fNsigmaTOF);
         //iCutTOFNSigma->AddPIDRange(3.0, 0.35, 1E20);
         AddCut(fCutQuality);
         AddCut(iCutTOFNSigma);
         SetCutScheme( Form("%s&%s",fCutQuality->GetName(), iCutTOFNSigma->GetName()) );
         break;

      case    AliRsnCutSetDaughterParticle::kTPCTOFpidKstarPP2010 :
         /* Set TPC  PID (if no TOF)*/
         // all   below   350         MeV: 5sigma
         // all   between 350 and 500 MeV: 3sigma
         // pios above   500         MeV: 2sigma
         // kaons between 500 and 700 MeV: 2sigma
         // kaons above   700         MeV: rejected
         iCutTPCNSigma->AddPIDRange(5.0, 0.0, 0.35);
         iCutTPCNSigma->AddPIDRange(3.0, 0.35, 0.5);
         if (fPID==AliPID::kPion)
            iCutTPCNSigma->AddPIDRange(2.0, 0.5, 1E20);
         if (fPID==AliPID::kKaon)
            iCutTPCNSigma->AddPIDRange(2.0, 0.5, 0.7);

         AddCut(fCutQuality);
         AddCut(iCutTOFMatch);
         AddCut(iCutTPCNSigma);

         /* set TPC+TOF PID*/
         iCutTPCTOFNSigma->SinglePIDRange(5.0);
         iCutTOFNSigma->AddPIDRange(3.0, 0.0, 1.5);
         iCutTOFNSigma->AddPIDRange(2.0, 1.5, 1E20);

         AddCut(iCutTPCTOFNSigma);
         AddCut(iCutTOFNSigma);

         // scheme:
         // quality & [ (TOF & TPCTOF) || (!TOFmatch & TPConly) ]
         SetCutScheme( Form("%s&((%s&%s)|((!%s)&%s))",fCutQuality->GetName(), iCutTPCTOFNSigma->GetName(), iCutTOFNSigma->GetName(), iCutTOFMatch->GetName(), iCutTPCNSigma->GetName()) ) ;
         break;

      case AliRsnCutSetDaughterParticle::kTOFpidKstarPbPb2010:
         if (fNsigmaTOF <= 0.0) {
            AliWarning(Form("Invalid number of sigmas required for %s. Setting default nSigma = 10",iCutTOFNSigma->GetName()));
            SetNsigmaForFastTOFpid(10.0);
         }
         iCutTOFNSigma->SinglePIDRange(fNsigmaTOF);
         //iCutTOFNSigma->AddPIDRange(3.0, 0.35, 1E20);
         iCutTPCTOFNSigma->SinglePIDRange(5.0); //5-sigma veto on tpc signal

         AddCut(fCutQuality);
         AddCut(iCutTOFNSigma);
         AddCut(iCutTPCTOFNSigma);
         SetCutScheme( Form("%s&%s&%s",fCutQuality->GetName(), iCutTOFNSigma->GetName(), iCutTPCTOFNSigma->GetName()) );
         break;

      case AliRsnCutSetDaughterParticle::kTOFTPCmismatchKstarPbPb2010:
         iCutTPCTOFNSigma->SinglePIDRange(5.0); //5-sigma veto on tpc signal
         AddCut(fCutQuality);
         AddCut(iCutTOFMatch);
         AddCut(iCutTPCTOFNSigma);
         SetCutScheme( Form("%s&%s&%s",fCutQuality->GetName(), iCutTOFMatch->GetName(), iCutTPCTOFNSigma->GetName()) );
         break;

   case AliRsnCutSetDaughterParticle::kTOFMatchTRD2010 :
         AddCut(fCutQuality);
         AddCut(iCutTOFMatch);
	 AddCut(iCutPhiTRD2010);
	 SetCutScheme( Form("%s&(%s)&(%s)",fCutQuality->GetName(), iCutTOFMatch->GetName(), iCutPhiTRD2010->GetName()) );
         break;
	 
   case AliRsnCutSetDaughterParticle::kTOFMatchNoTRD2010 :
         AddCut(fCutQuality);
         AddCut(iCutTOFMatch);
	 AddCut(iCutPhiNoTRD2010);
	 SetCutScheme( Form("%s&(%s)&(%s)",fCutQuality->GetName(), iCutTOFMatch->GetName(), iCutPhiNoTRD2010->GetName()) );
         break;

      case AliRsnCutSetDaughterParticle::kTOFpidKstarPbPbTRD2010:
         if (fNsigmaTOF <= 0.0) {
            AliWarning(Form("Invalid number of sigmas required for %s. Setting default nSigma = 10",iCutTOFNSigma->GetName()));
            SetNsigmaForFastTOFpid(10.0);
         }
         iCutTOFNSigma->SinglePIDRange(fNsigmaTOF);
         //iCutTOFNSigma->AddPIDRange(3.0, 0.35, 1E20);
         iCutTPCTOFNSigma->SinglePIDRange(5.0); //5-sigma veto on tpc signal

         AddCut(fCutQuality);
         AddCut(iCutTOFNSigma);
         AddCut(iCutTPCTOFNSigma);
	 AddCut(iCutPhiTRD2010);
         SetCutScheme( Form("%s&%s&%s%s",fCutQuality->GetName(), iCutTOFNSigma->GetName(), iCutTPCTOFNSigma->GetName(),iCutPhiTRD2010->GetName()) );
         break;

 case AliRsnCutSetDaughterParticle::kTOFpidKstarPbPbNoTRD2010:
         if (fNsigmaTOF <= 0.0) {
            AliWarning(Form("Invalid number of sigmas required for %s. Setting default nSigma = 10",iCutTOFNSigma->GetName()));
            SetNsigmaForFastTOFpid(10.0);
         }
         iCutTOFNSigma->SinglePIDRange(fNsigmaTOF);
         //iCutTOFNSigma->AddPIDRange(3.0, 0.35, 1E20);
         iCutTPCTOFNSigma->SinglePIDRange(5.0); //5-sigma veto on tpc signal

         AddCut(fCutQuality);
         AddCut(iCutTOFNSigma);
         AddCut(iCutTPCTOFNSigma);
	 AddCut(iCutPhiNoTRD2010);
         SetCutScheme( Form("%s&%s&%s%s",fCutQuality->GetName(), iCutTOFNSigma->GetName(), iCutTPCTOFNSigma->GetName(),iCutPhiNoTRD2010->GetName()) );
         break;

   case AliRsnCutSetDaughterParticle::kTOFMatchTPCpidNsigma :
         if (fNsigmaTPC <= 0.0) {
            AliWarning(Form("Invalid number of sigmas required for %s. Setting default nSigma = 10",iCutTPCNSigma->GetName()));
            SetNsigmaForFastTPCpid(10.0);
         }
         iCutTPCNSigma->SinglePIDRange(fNsigmaTPC);
         AddCut(fCutQuality);
         AddCut(iCutTPCNSigma);
	 AddCut(iCutTOFMatch);	
         SetCutScheme( Form("%s&%s&%s",fCutQuality->GetName(), iCutTPCNSigma->GetName(), iCutTOFMatch->GetName()) );
         break;
    default :
         break;
   }

}



