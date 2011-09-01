/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

////////////////////////////////////////////////////////////////////////////////
//
//  This class contains all code which is used to compute any of the values
//  which can be of interest within a resonance analysis. Besides the obvious
//  invariant mass, it allows to compute other utility values on all possible
//  targets, in order to allow a wide spectrum of binning and checks.
//  When needed, this object can also define a binning in the variable which
//  it is required to compute, which is used for initializing axes of output
//  histograms (see AliRsnFunction).
//  The value computation requires this object to be passed the object whose
//  informations will be used. This object can be of any allowed input type
//  (track, pair, event), then this class must inherit from AliRsnTarget.
//  Then, when value computation is attempted, a check on target type is done
//  and computation is successful only if expected target matches that of the
//  passed object.
//  In some cases, the value computation can require a support external object,
//  which must then be passed to this class. It can be of any type inheriting
//  from TObject.
//
//  authors: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
//           M. Vala (martin.vala@cern.ch)
//
////////////////////////////////////////////////////////////////////////////////

#include "Riostream.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "AliLog.h"
#include "AliPID.h"
#include "AliPIDResponse.h"

#include "AliRsnDaughter.h"
#include "AliRsnEvent.h"

#include "AliRsnMiniMonitor.h"

ClassImp(AliRsnMiniMonitor)

//__________________________________________________________________________________________________
AliRsnMiniMonitor::AliRsnMiniMonitor() :
   TNamed(), 
   fType(kTypes), 
   fCutID(-1), 
   fCharge(0),
   fListID(-1), 
   fList(0x0)
{
//
// Dummy constructor
//
}

//__________________________________________________________________________________________________
AliRsnMiniMonitor::AliRsnMiniMonitor(const char *name, EType type, Int_t cutID) :
   TNamed(name, ""), 
   fType(type), 
   fCutID(cutID), 
   fCharge(0),
   fListID(-1), 
   fList(0x0)
{
//
// Default constructor
//
}
   
//__________________________________________________________________________________________________
AliRsnMiniMonitor::AliRsnMiniMonitor(const AliRsnMiniMonitor& copy) :
   TNamed(copy), 
   fType(copy.fType), 
   fCutID(copy.fCutID), 
   fCharge(copy.fCharge),
   fListID(copy.fListID), 
   fList(copy.fList)
{
//
// Copy constructor
//
}
   
//__________________________________________________________________________________________________
AliRsnMiniMonitor& AliRsnMiniMonitor::operator=(const AliRsnMiniMonitor& copy) 
{
//
// Assignment operator
//
   
   TNamed::operator=(copy); 
   fType = copy.fType; 
   fCutID = copy.fCutID; 
   fCharge = copy.fCharge;
   fListID = copy.fListID; 
   fList = copy.fList; 
   
   return (*this); 
}

//__________________________________________________________________________________________________
Bool_t AliRsnMiniMonitor::Init(const char *name, TList *list)
{
//
// Initialize this output histogram and put into the passed list
//

   // name
   TString sname(name);
   sname += '_';
   sname += GetName();
   
   // check list
   fList = list;
   if (!list) {
      AliError("No list!");
      return kFALSE;
   }
   
   // reset histogram
   TH1 *histogram = 0x0;

   switch (fType) {
      case kTrackPt:
         sname += "_TrackPt_cut";
         sname += fCutID;
         histogram = new TH1F(sname.Data(), "", 100, 0.0, 10.0);
         break;
      case kdEdxTPCvsP: 
         sname += "_TPCsignal";
         histogram = new TH2F(sname.Data(), "", 500, 0.0, 5.0, 1000, 0.0, 1000.0);
         break;
      case ktimeTOFvsPPion:
         sname += "_TOFsignalPi";
         histogram = new TH2F(sname.Data(), "", 500, 0.0, 5.0, 1000, 0.0, 200000.0);
      case ktimeTOFvsPKaon:
         sname += "_TOFsignalK";
         histogram = new TH2F(sname.Data(), "", 500, 0.0, 5.0, 1000, 0.0, 200000.0);
         break;
      case ktimeTOFvsPProton:
         sname += "_TOFsignalP";
         histogram = new TH2F(sname.Data(), "", 500, 0.0, 5.0, 1000, 0.0, 200000.0);
         break;
      default:
         AliError("Wrong enum type");
         return kFALSE;
   }
   
   // add to list
   if (histogram && fList) {
      histogram->Sumw2();
      fList->Add(histogram);
      fListID = fList->IndexOf(histogram);
      AliInfo(Form("Histogram '%s' added to list in slot #%d", histogram->GetName(), fListID));
      return kTRUE;
   }
   
   return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliRsnMiniMonitor::Fill(AliRsnDaughter *track, AliRsnEvent *event)
{
//
// Fill the histogram
//

   // retrieve object from list
   if (!fList) {
      AliError("List pointer is NULL");
      return kFALSE;
   }
   TObject *obj = fList->At(fListID);
   if (!obj) {
      AliError("List object is NULL");
      return kFALSE;
   }

   Double_t valueX, valueY;
   AliVTrack *vtrack = track->Ref2Vtrack();
   
   AliPIDResponse *pid = event->GetPIDResponse();

   switch (fType) {
      case kTrackPt: 
         if (!vtrack) {
            AliWarning("Required vtrack for this value");
            return kFALSE;
         }
         if (fCharge == '+' && vtrack->Charge() <= 0) return kFALSE;
         if (fCharge == '-' && vtrack->Charge() >= 0) return kFALSE;
         if (fCharge == '0' && vtrack->Charge() != 0) return kFALSE;
         valueX = vtrack->Pt();
         ((TH1F*)obj)->Fill(valueX);
         return kTRUE;
      case kdEdxTPCvsP: 
         if (!vtrack) {
            AliWarning("Required vtrack for this value");
            return kFALSE;
         }
         valueX = vtrack->GetTPCmomentum();
         valueY = vtrack->GetTPCsignal();
         ((TH2F*)obj)->Fill(valueX, valueY);
         return kTRUE;
      case ktimeTOFvsPPion:
         if (!vtrack) {
            AliWarning("Required vtrack for this value");
            return kFALSE;
         }
         valueX = vtrack->P();
         //valueY = vtrack->GetTOFsignal();
         valueY = 1E20;
         if (pid) valueY = pid->NumberOfSigmasTOF(vtrack, AliPID::kPion);
         ((TH2F*)obj)->Fill(valueX, valueY);
         return kTRUE;
      case ktimeTOFvsPKaon:
         if (!vtrack) {
            AliWarning("Required vtrack for this value");
            return kFALSE;
         }
         valueX = vtrack->P();
         //valueY = vtrack->GetTOFsignal();
         valueY = 1E20;
         if (pid) valueY = pid->NumberOfSigmasTOF(vtrack, AliPID::kKaon);
         ((TH2F*)obj)->Fill(valueX, valueY);
         return kTRUE;
      case ktimeTOFvsPProton:
         if (!vtrack) {
            AliWarning("Required vtrack for this value");
            return kFALSE;
         }
         valueX = vtrack->P();
         //valueY = vtrack->GetTOFsignal();
         valueY = 1E20;
         if (pid) valueY = pid->NumberOfSigmasTOF(vtrack, AliPID::kProton);
         ((TH2F*)obj)->Fill(valueX, valueY);
         return kTRUE;
      default:
         AliError("Invalid value type");
         return kFALSE;
   }
}
