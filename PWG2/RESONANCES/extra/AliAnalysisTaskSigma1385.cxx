/**************************************************************************
 * Authors : Massimo Venaruzzo (massimo.venaruzzo@ts.infn.it)             *
 *      Enrico Fragiacomo (enrico.fragiacomo@ts.infn.it)             *
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

//-----------------------------------------------------------------
//                 AliAnalysisTaskSigma1385 class
//-----------------------------------------------------------------

class TTree;
class TParticle;
class TVector3;

#include "AliAnalysisManager.h"
#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <AliStack.h>

class AliESDVertex;
class AliESDv0;
class AliAODv0;

#include <iostream>

#include "TList.h"
#include "TH1.h"
#include "TNtuple.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TChain.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliCascadeVertexer.h"
#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliAnalysisTaskSigma1385.h"
#include "AliESDtrackCuts.h"
#include "AliESDpid.h"
#include "AliTOFT0maker.h"
#include "AliTOFcalib.h"
#include "AliCDBManager.h"
#include "AliESDtrackCuts.h"

ClassImp(AliAnalysisTaskSigma1385)

//________________________________________________________________________
AliAnalysisTaskSigma1385::AliAnalysisTaskSigma1385()
   : AliAnalysisTaskSE(),
     
     //-------------------------------- For PID
     
     fisMC(0),
     fIsMC(fisMC),
     fCheckITS(kTRUE),
     fCheckTPC(kTRUE),
     fCheckTOF(kTRUE),
     fUseGlobal(kTRUE),
     fUseITSSA(kTRUE),
     fMaxITSband(3.0),
     fTPCpLimit(0.35),
     fMinTPCband(3.0),
     fMaxTPCband(5.0),
     fESDpid(0x0),
     fTOFmaker(0x0),
     fTOFcalib(0x0),
     fTOFcalibrateESD(!fisMC),
     fTOFcorrectTExp(kTRUE),
     fTOFuseT0(kTRUE),
     fTOFtuneMC(fisMC),
     fTOFresolution(100.0),
     fMinTOF(-2.5),
     fMaxTOF(3.5),
     fLastRun(-1),
     fAnalysisType("ESD"), fCollidingSystems(0), fDataType("REAL"), fListHistCascade(0), 
     fHistEventMultiplicity(0), fHistEventMultiplicityRAVS(0), 
     fNtuple1(0), fNtuple2(0), fNtuple3(0), fNtuple4(0)
 


     //--------------------------------

{
   // Dummy Constructor
   Int_t i;
   for (i = 0; i < 5; i++) fOkTrack[i] = kFALSE;
}

//________________________________________________________________________
AliAnalysisTaskSigma1385::AliAnalysisTaskSigma1385(const char *name)
   : AliAnalysisTaskSE(name),
       

     //-------------------------------- For PID

     fisMC(0),
     fIsMC(fisMC),
     fCheckITS(kTRUE),
     fCheckTPC(kTRUE),
     fCheckTOF(kTRUE),
     fUseGlobal(kTRUE),
     fUseITSSA(kTRUE),
     fMaxITSband(3.0),
     fTPCpLimit(0.35),
     fMinTPCband(3.0),
     fMaxTPCband(5.0),
     fESDpid(0x0),
     fTOFmaker(0x0),
     fTOFcalib(0x0),
     fTOFcalibrateESD(!fisMC),
     fTOFcorrectTExp(kTRUE),
     fTOFuseT0(kTRUE),
     fTOFtuneMC(fisMC),
     fTOFresolution(100.0),
     fMinTOF(-2.5),
     fMaxTOF(3.5),
     fLastRun(-1),
     fAnalysisType("ESD"), fCollidingSystems(0), fDataType("REAL"), fListHistCascade(0), 
     fHistEventMultiplicity(0), fHistEventMultiplicityRAVS(0), 
     fNtuple1(0), fNtuple2(0), fNtuple3(0), fNtuple4(0)

     //--------------------------------
{

   // Output slot #0 writes into a TList container (Cascade)
   DefineOutput(1, TList::Class());
   
   Int_t i;
   for (i = 0; i < 5; i++) fOkTrack[i] = kFALSE;
}

//________________________________________________________________________
void AliAnalysisTaskSigma1385::UserCreateOutputObjects()
{
   fListHistCascade = new TList();

   if (! fHistEventMultiplicity) {
      fHistEventMultiplicity   = new TH1F("fHistEventMultiplicity" , "Nb of Events" , 4, -1.0, 3.0);
      fListHistCascade->Add(fHistEventMultiplicity);
   }

   if (! fHistEventMultiplicityRAVS) {
      fHistEventMultiplicityRAVS   = new TH1F("fHistEventMultiplicityRAVS" , "Nb of Events Rejected After Vertex selection" , 4, -1.0, 3.0);
      fListHistCascade->Add(fHistEventMultiplicityRAVS);
   }

   if (! fNtuple1) {
      fNtuple1 = new TNtuple("fNtuple1", "Ntuple1", "TrkNmb");
      fNtuple1->SetDirectory(0);
      fListHistCascade->Add(fNtuple1);
   }

   if (! fNtuple2) {
      fNtuple2 = new TNtuple("fNtuple2", "Ntuple2", "s:dcal:lCosPoinAn:lDaugDCA:lambdap:lambdapt:lambdamass");
      fNtuple2->SetDirectory(0);
      fListHistCascade->Add(fNtuple2);
   }

   if (! fNtuple3) {
      fNtuple3 = new TNtuple("fNtuple3", "Ntuple3", "c:dcapi:ppi:ptpi:bachphi:bachtheta:okPiTPC:okPiTOF");
      fNtuple3->SetDirectory(0);
      fListHistCascade->Add(fNtuple3);
   }

   if (! fNtuple4) {
      fNtuple4 = new TNtuple("fNtuple4", "Ntuple4", "dca:mc:phi:theta:eta:y:pt:p:opang:invmass");
      fListHistCascade->Add(fNtuple4);
   }


}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskSigma1385::UserExec(Option_t *)
{

   // Main loop
   // Called for each event


   Info("AliAnalysisTaskSigma1385", "Starting UserExec");

   AliMCEventHandler* eventHandler;
   AliMCEvent* mcEvent = 0;

   if (fDataType == "SIM") {

      eventHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
      if (!eventHandler) {
         Printf("ERROR: Could not retrieve MC event handler");
         return;
      }

      mcEvent = eventHandler->MCEvent();
      if (!mcEvent) {
         Printf("ERROR: Could not retrieve MC event");
         return;
      }

   }

   AliStack* stack = 0;
   if (fDataType == "SIM") {stack = mcEvent->Stack(); fIsMC = 1; fisMC = 1;}

   AliESDEvent *lESDevent = 0x0;
   AliAODEvent *lAODevent = 0x0;


   // Connect to the InputEvent
   Int_t ncascades = -1;
   if (fAnalysisType == "ESD") {
      lESDevent = dynamic_cast<AliESDEvent*>(InputEvent());
      if (lESDevent) ncascades = lESDevent->GetNumberOfCascades();
   } else if (fAnalysisType == "AOD") {
      lAODevent = dynamic_cast<AliAODEvent*>(InputEvent());
      if (lAODevent) ncascades = lAODevent->GetNumberOfCascades();
   } else {
      Printf("ERROR: neither lESDevent nor lAODevent are available \n");
      return;
   }


   //------------------------------for PID

   SetCheckITS(kTRUE);
   SetCheckTPC(kTRUE);
   SetCheckTOF(kTRUE);

// ----> set TPC range for PID and calibration
   SetTPCrange(5.0, 3.0);
   SetTPCpLimit(0.35);

   // ----> set ITS range for PID
   SetITSband(4.0);

   // ----> set TPC calibration
   if (fDataType == "SIM") SetTPCpar(2.15898 / 50.0, 1.75295E1, 3.40030E-9, 1.96178, 3.91720);
   else       SetTPCpar(1.41543 / 50.0, 2.63394E1, 5.0411E-11, 2.12543, 4.88663);

   // ----> set the TOF calibration depending on type of input (sim/data)
   SetTOFcorrectTExp(kTRUE);
   SetTOFuseT0(kTRUE);
   SetTOFresolution(100.0);
   if (fDataType == "SIM") {
      SetTOFcalibrateESD(kFALSE);
      SetTOFtuneMC(kTRUE);
   } else {

      SetTOFcalibrateESD(kTRUE);
      SetTOFtuneMC(kFALSE);
   }


   if (!fESDpid) {
      fESDpid = new AliESDpid;
      fESDpid->GetTPCResponse().SetBetheBlochParameters(fTPCpar[0], fTPCpar[1], fTPCpar[2], fTPCpar[3], fTPCpar[4]);
   }

   // initialize DB to current run
   Int_t run = lESDevent->GetRunNumber();
   if (run != fLastRun) {
      //cout << "Run = " << run << " -- LAST = " << fLastRun << endl;
      fLastRun = run;

      // setup TOF maker & calibration
      if (!fTOFcalib) fTOFcalib = new AliTOFcalib;
      fTOFcalib->SetCorrectTExp(fTOFcorrectTExp);
      if (!fTOFmaker) fTOFmaker = new AliTOFT0maker(fESDpid, fTOFcalib);
      fTOFmaker->SetTimeResolution(fTOFresolution);

      AliCDBManager *cdb = AliCDBManager::Instance();
      cdb->ClearCache(); // suggestion by Annalisa
      cdb->Clear();      // suggestion by Annalisa
      cdb->SetDefaultStorage("raw://");
      cdb->SetRun(run);
      fTOFcalib->SetCorrectTExp(fTOFcorrectTExp);
      fTOFcalib->Init();
   }

   // if required, calibrate the TOF t0 maker with current event
   if (fTOFcalibrateESD) fTOFcalib->CalibrateESD(lESDevent);
   if (fTOFtuneMC) fTOFmaker->TuneForMC(lESDevent);
   if (fTOFuseT0) {
      fTOFmaker->ComputeT0TOF(lESDevent);
      fTOFmaker->ApplyT0TOF(lESDevent);
      fESDpid->MakePID(lESDevent, kFALSE, 0.);
   }


   //--------------------------------------------------------


   fHistEventMultiplicity->Fill(1);

   //Some Quantities to characterize the event

   Double_t b = lESDevent->GetMagneticField();
   Int_t trackNumber = lESDevent->GetNumberOfTracks();


   //---------------------------
   // new part from AliCascadeVertexer (use vertexer for reconstructing Sigma(1385)
   //
   //const AliESDVertex *vtxT3D=lESDevent->GetPrimaryVertex();
   const AliESDVertex *vtxT3D = lESDevent->GetPrimaryVertexTracks();
   if (vtxT3D->GetNContributors() < 1) {
      // SPD vertex
      vtxT3D = lESDevent->GetPrimaryVertexSPD();
      if (vtxT3D->GetNContributors() < 1) {


         fHistEventMultiplicityRAVS->Fill(1);
         return;
      }


   }

   Double_t xPrimaryVertex = vtxT3D->GetXv();
   Double_t yPrimaryVertex = vtxT3D->GetYv();
   Double_t zPrimaryVertex = vtxT3D->GetZv();

   if (zPrimaryVertex > 10 || zPrimaryVertex < -10)  return;

   //-------------------------------------------------------
   // New Part about tracks global selection criteria
   //-------------------------------------------------------

   AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");

   esdTrackCuts->SetAcceptKinkDaughters(0); // 0 = kFalse
   esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
   esdTrackCuts->SetMaxChi2PerClusterTPC(4);
   esdTrackCuts->SetMinNClustersTPC(70);


   //-------------------------------------------------------
   // loops over V0s
   //stores relevant V0s in an array
   Int_t nV0 = (Int_t)lESDevent->GetNumberOfV0s();

   TObjArray vtcs(nV0);
   for (Int_t i = 0; i < nV0; i++) {
      AliESDv0 *v = lESDevent->GetV0(i);
      if (v->GetOnFlyStatus()) continue; // if kTRUE, then this V0 is recontructed

      vtcs.AddLast(v);
   }
   nV0 = vtcs.GetEntriesFast();

   //-------------------------------------------------------
   // loops over bachelor tracks
   // stores relevant tracks in another array

   Int_t nentr = (Int_t)lESDevent->GetNumberOfTracks();
   TArrayI trk(nentr); Int_t ntr = 0;

   for (Int_t i = 0; i < nentr; i++) {
      AliESDtrack *esdtr = lESDevent->GetTrack(i);

      if (!(esdtr->GetStatus() & AliESDtrack::kITSrefit)) continue;
      if (!(esdtr->GetStatus() & AliESDtrack::kTPCrefit)) continue;


      if (!(esdTrackCuts->AcceptTrack(esdtr))) continue;

      trk[ntr++] = i;
   }



   //-----------------------------------------------------
   // nested loops over V0s and bachelors
   Float_t massLambda = 1.11568;
   Int_t ncasc = 0;
   AliCascadeVertexer *cascvert = new AliCascadeVertexer();
   for (Int_t i = 0; i < nV0; i++) { //loop on V0s

      // Lambda
      AliESDv0 *v = (AliESDv0*)vtcs.UncheckedAt(i);
      Int_t strness = 0;
      Double_t lambdaMass = 0;
      Float_t lambdaP = 0;
      Float_t lambdaPt = 0;
      Float_t lambdaDCA = 0; // DCA between Lambda and Primary Vertex
      Double_t v0cospointangle = 0;
      Float_t v0daughtersDCA = 0;


      // Lambda quality cuts
      UInt_t lIdxPosXi    = (UInt_t) TMath::Abs(v->GetPindex());
      UInt_t lIdxNegXi    = (UInt_t) TMath::Abs(v->GetNindex());

      AliESDtrack *pTrackXi     = lESDevent->GetTrack(lIdxPosXi);
      AliESDtrack *nTrackXi     = lESDevent->GetTrack(lIdxNegXi);

      // Filter like-sign V0
      if ( (TMath::Abs(pTrackXi->GetSign()) - TMath::Abs(nTrackXi->GetSign()) ) < 0.1) continue;

      // WARNING: the following selections cannot be done for AOD yet...


      v->ChangeMassHypothesis(kLambda0); // the v0 must be Lambda


      if ((TMath::Abs(v->GetEffMass() - massLambda)) < 0.2) {

         if (!(pTrackXi->GetStatus() & AliESDtrack::kTPCrefit))  continue;
         if (!(nTrackXi->GetStatus() & AliESDtrack::kTPCrefit))  continue;


         if ((pTrackXi->GetTPCNcls()) < 70) continue;
         if ((nTrackXi->GetTPCNcls()) < 70) continue;

         strness = 1; lambdaMass = v->GetEffMass(); lambdaP = v->P(); lambdaPt = v->Pt(); lambdaDCA = v->GetD(xPrimaryVertex, yPrimaryVertex, zPrimaryVertex); v0cospointangle = v->GetV0CosineOfPointingAngle(); v0daughtersDCA = v->GetDcaV0Daughters();



      }


      v->ChangeMassHypothesis(kLambda0Bar); // the v0 must be Anti Lambda



      if ((TMath::Abs(v->GetEffMass() - massLambda)) < 0.2) {

         if (!(pTrackXi->GetStatus() & AliESDtrack::kTPCrefit))  continue;
         if (!(nTrackXi->GetStatus() & AliESDtrack::kTPCrefit))  continue;


         if ((pTrackXi->GetTPCNcls()) < 70) continue;
         if ((nTrackXi->GetTPCNcls()) < 70) continue;

         Int_t temp = strness + 1; strness = -1 * temp; lambdaMass = v->GetEffMass(); lambdaP = v->P(); lambdaPt = v->Pt(); lambdaDCA = v->GetD(xPrimaryVertex, yPrimaryVertex, zPrimaryVertex); v0cospointangle = v->GetV0CosineOfPointingAngle(); v0daughtersDCA = v->GetDcaV0Daughters();


      }

      if (strness == 0) continue;


      for (Int_t j = 0; j < ntr; j++) { //loop on tracks

         // Pion bachelor
         Int_t bidx = trk[j];
         Int_t bachTPCcls = 0;


         if ((bidx == v->GetIndex(0)) || (bidx == v->GetIndex(1))) continue; // bach and V0's daughter's must be different!

         AliESDtrack *btrk = lESDevent->GetTrack(bidx);

         if (!(btrk->GetStatus() & AliESDtrack::kTPCrefit)) continue;

         if ((bachTPCcls = btrk->GetTPCNcls()) < 70) continue;

         //Bool_t *IsOkTrack = IsSelected(btrk); //for Alberto's PID

         //Bool_t *okTrack = new Bool_t[5];

         //for (Int_t k = 0; k < 5; k++) okTrack[k] = IsOkTrack[k];
         
         IsSelected(btrk);


         Int_t bachCharge = btrk->Charge();
         Float_t pionDCA = TMath::Abs(btrk->GetD(xPrimaryVertex, yPrimaryVertex, b)); // DCA between Bachelor and Primary Vertex

         Double_t bachphi = btrk->Phi();
         Double_t bachtheta = btrk->Theta();
         Double_t bachmass = 0.13957;
         Double_t bachp = btrk->GetP();
         Double_t bachpt = btrk->Pt();


         // Distance between Lambda and Pion bachelor
         AliESDv0 v0(*v), *pv0 = &v0;
         AliExternalTrackParam bt(*btrk), *pbt = &bt;
         Double_t dca = cascvert->PropagateToDCA(pv0, pbt, b); // distance bwn V0 and bach (in cm)
         if (dca > 10.0) continue;  // Note: was 0.1! Enlarged->further filter in second pass analysis


         AliESDcascade cascade(*pv0, *pbt, bidx); //constucts a sigma1385 candidate
         AliESDcascade *xi = &cascade;


         UInt_t lBachIdx   = (UInt_t) TMath::Abs(xi->GetBindex());
         AliESDtrack *bachTrackXi   = lESDevent->GetTrack(lBachIdx);



         Short_t mcTrue=0;
         if (fDataType == "SIM") {

            Int_t pLabel    = TMath::Abs(pTrackXi->GetLabel());
            Int_t nLabel    = TMath::Abs(nTrackXi->GetLabel());
            Int_t bachLabel = TMath::Abs(bachTrackXi->GetLabel());

            Int_t motherpLabel = stack->Particle(pLabel)->GetFirstMother();
            Int_t mothernLabel = stack->Particle(nLabel)->GetFirstMother();
            Int_t motherbachLabel = stack->Particle(bachLabel)->GetFirstMother();


            //cout<< "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
            //cout<< "plabel " << pLabel << " nlabel " << nLabel << " mother p " << motherpLabel << " mother n " << mothernLabel << " bachlabel " << bachLabel << " mother bach " << motherbachLabel << endl;
            //cout<< "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;

            if (motherpLabel > -1 && mothernLabel > -1) {
               TParticle *plambda = stack->Particle(motherpLabel);
               Int_t grandmother = plambda->GetFirstMother();

               motherpLabel = TMath::Abs(motherpLabel);
               mothernLabel = TMath::Abs(mothernLabel);
               //motherbachLabel = TMath::Abs(motherbachLabel);

               //cout<< "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
               //cout<< "plabel " << pLabel << " nlabel " << nLabel << " mother p " << motherpLabel << " mother n " << mothernLabel << " mother lambda " << grandmother << " bachlabel " << bachLabel << " mother bach " << motherbachLabel << endl;
               //cout<< "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;

               if (motherbachLabel > -1) {


                  //cout<< "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
                  //cout<< "mother lambda " << grandmother << " mother bach " << motherbachLabel << " PDG Code "<< stack->Particle(grandmother)->GetPdgCode() << endl;
                  //cout<< "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;

                  if ((motherpLabel == mothernLabel) && (grandmother == motherbachLabel) &&
                      ((TMath::Abs(stack->Particle(grandmother)->GetPdgCode()) == 3114) ||
                       (TMath::Abs(stack->Particle(grandmother)->GetPdgCode()) == 3224))) mcTrue = 1;
		       
		  if( (motherpLabel == mothernLabel) && (grandmother == motherbachLabel) &&
			 (TMath::Abs(stack->Particle(grandmother)->GetPdgCode())==3312) ) mcTrue = 2;

               }

            }

         }

         Double_t lBachMomX       = 0., lBachMomY  = 0., lBachMomZ   = 0.;
         Double_t lPMom[3] = {0., 0., 0.,};
         Double_t lNMom[3] = {0., 0., 0.,};
         Float_t lLambdaMomX       = 0., lLambdaMomY  = 0., lLambdaMomZ   = 0.;
         xi->GetBPxPyPz(lBachMomX,  lBachMomY,  lBachMomZ);
         xi->GetPPxPyPz(lPMom[0],  lPMom[1],  lPMom[2]);
         xi->GetNPxPyPz(lNMom[0],  lNMom[1],  lNMom[2]);
         lLambdaMomX = lPMom[0] + lNMom[0];
         lLambdaMomY = lPMom[1] + lNMom[1];
         lLambdaMomZ = lPMom[2] + lNMom[2];

         Float_t  lRapXi    = xi->RapXi();
         Float_t  lEta      = xi->Eta();
         Float_t  lTheta    = xi->Theta();
         Float_t  lPhi      = xi->Phi();
         Float_t  lPt       = xi->Pt();
         Float_t   lP     = xi->P();

         // Support variables for invariant mass calculation
         TLorentzVector lambda, pion, sigma;
         Double_t eLambda = TMath::Sqrt(lLambdaMomX * lLambdaMomX + lLambdaMomY * lLambdaMomY + lLambdaMomZ * lLambdaMomZ + lambdaMass * lambdaMass);
         Double_t ePion = TMath::Sqrt(lBachMomX * lBachMomX + lBachMomY * lBachMomY + lBachMomZ * lBachMomZ + bachmass * bachmass) ;

         lambda.SetPxPyPzE(lLambdaMomX, lLambdaMomY, lLambdaMomZ, eLambda);
         pion.SetPxPyPzE(lBachMomX, lBachMomY, lBachMomZ, ePion);

         sigma = lambda + pion;

         Double_t openingangle, invmass = 0;

         openingangle = TMath::ACos((lBachMomX * lLambdaMomX + lBachMomY * lLambdaMomY + lBachMomZ * lLambdaMomZ) / ((TMath::Sqrt(lBachMomX * lBachMomX + lBachMomY * lBachMomY + lBachMomZ * lBachMomZ)) * (TMath::Sqrt(lLambdaMomX * lLambdaMomX + lLambdaMomY * lLambdaMomY +
                                    lLambdaMomZ * lLambdaMomZ))));

         invmass = sigma.M();


         fNtuple1->Fill(trackNumber);
         fNtuple2->Fill((1.*strness), lambdaDCA, v0cospointangle, v0daughtersDCA , lambdaP, lambdaPt, lambdaMass);
         //fNtuple3->Fill((1.*bachCharge), pionDCA, bachp, bachpt, bachphi, bachtheta, okTrack[1], okTrack[2]);
         fNtuple3->Fill((1.*bachCharge), pionDCA, bachp, bachpt, bachphi, bachtheta, fOkTrack[1], fOkTrack[2]);
         fNtuple4->Fill(dca, (1.*mcTrue), lPhi, lTheta, lEta, lRapXi, lPt, lP, openingangle, invmass);

         ncasc++;

         //delete IsOkTrack;
         //delete okTrack;


      } // end loop tracks
   } // end loop V0s

   Info("AliAnalysisTaskSigma1385", "Number of reconstructed Sigma(1385): %d", ncasc);


   // Post output data.
   PostData(1, fListHistCascade);

}

//________________________________________________________________________

//Bool_t *AliAnalysisTaskSigma1385::IsSelected(AliESDtrack *track)
void AliAnalysisTaskSigma1385::IsSelected(AliESDtrack *track)
{
//
//
   //Bool_t okTrack[5];

   //for (Int_t i = 0; i < 5; i++) okTrack[i] = kFALSE;
   for (Int_t i = 0; i < 5; i++) fOkTrack[i] = kFALSE;


   AliITSPIDResponse itsrsp(fIsMC);


   ULong_t  status;
   Int_t    k, nITS;
   Double_t times[10], tpcNSigma, tpcMaxNSigma, itsSignal, itsNSigma, mom, tofTime, tofSigma, tofRef, tofRel;
   Bool_t   okTOF = kFALSE;
   Bool_t   okTPC = kFALSE;
   Bool_t   okITS = kFALSE;
   Bool_t   isTPC = kFALSE;
   Bool_t   isITSSA = kFALSE;
   Bool_t   isTOF = kFALSE;
   UChar_t  itsCluMap;

   // get commonly used variables
   status  = (ULong_t)track->GetStatus();
   mom     = track->P();
   isTPC   = ((status & AliESDtrack::kTPCin)  != 0);
   isITSSA = (!isTPC && (status & AliESDtrack::kITSrefit) != 0 && (status & AliESDtrack::kITSpureSA) == 0 && (status & AliESDtrack::kITSpid) != 0);
   isTOF   = (((status & AliESDtrack::kTOFout) != 0) && ((status & AliESDtrack::kTIME) != 0) /* && mom > TMath::Max(b1, b2)*/);


   // check if the track type matches what is required
   if (!isTPC && !isITSSA) {
      AliDebug(AliLog::kDebug + 2, "Track is not either a TPC track or a ITS standalone. Rejected");
      //return okTrack;
      return;

   } else if (isTPC && !fUseGlobal) {
      AliDebug(AliLog::kDebug + 2, "Global tracks not used. Rejected");
      //return okTrack;
      return;

   } else if (isITSSA && !fUseITSSA) {
      AliDebug(AliLog::kDebug + 2, "ITS standalone not used. Rejected");
      //return okTrack;
      return;
   }

   // does a preliminary check on TOF values, if necessary
   // then, if the reference time or TOF signal are meaningless
   // even if the 'isTOF' flag is true, switch it to false
   if (isTOF) {
      track->GetIntegratedTimes(times);
      tofTime  = (Double_t)track->GetTOFsignal();
      tofSigma = fTOFmaker->GetExpectedSigma(mom, times[AliPID::kPion], AliPID::ParticleMass(AliPID::kPion));
      tofRef   = times[AliPID::kPion];
      if (tofRef <= 0.0 && tofSigma <= 0.0) isTOF = kFALSE;
   }



   // check TPC dE/dx
   if (isTPC) { // this branch is entered by all global tracks
      // check TPC dE/dx:
      if (fCheckTPC) {
         tpcNSigma = TMath::Abs(fESDpid->NumberOfSigmasTPC(track, AliPID::kPion));
         if (track->GetInnerParam()->P() > fTPCpLimit) tpcMaxNSigma = fMinTPCband; else tpcMaxNSigma = fMaxTPCband;
         okTPC = (tpcNSigma <= tpcMaxNSigma);
         AliDebug(AliLog::kDebug + 2, Form("TPC nsigma = %f -- max = %f -- cut %s", tpcNSigma, tpcMaxNSigma, (okTPC ? "passed" : "failed")));
      } else {
         // if TPC is not checked, it is as if all tracks do pass the cut
         okTPC = kTRUE;
         AliDebug(AliLog::kDebug + 2, "TPC not checked, track accepted");
      }

      // check TOF (only if flags are OK)
      if (fCheckTOF) {
         if (isTOF) {
            // TOF can be checked only when track is matched there
            track->GetIntegratedTimes(times);
            tofTime  = (Double_t)track->GetTOFsignal();
            tofSigma = fTOFmaker->GetExpectedSigma(mom, times[AliPID::kPion], AliPID::ParticleMass(AliPID::kPion));
            tofRef   = times[AliPID::kPion];

            tofRel   = (tofTime - tofRef) / tofSigma;
            okTOF    = ((tofRel >= fMinTOF) && (tofRel <= fMaxTOF));
            AliDebug(AliLog::kDebug + 2, Form("TOF nsigma = %f -- range = %f %f -- cut %s", tofRel, fMinTOF, fMaxTOF, (okTOF ? "passed" : "failed")));
         } else {
            // if TOF is not matched, the answer depends on TPC:
            // - if TPC is required, track is checked only there and TOF simply ignored
            // - if TPC is not required, track is rejected when TOF does not match it, if TOF check is required
            if (fCheckTPC) okTOF = kTRUE; else okTOF = kFALSE;
         }
      } else {
         okTOF = kTRUE;
      }

   } else if (isITSSA) { // this branch is entered by all ITS standalone tracks
      // check dE/dx only if this is required, otherwise ITS standalone are just added but not checked for PID
      if (fCheckITS) {
         itsSignal = track->GetITSsignal();
         itsCluMap = track->GetITSClusterMap();
         nITS      = 0;
         for (k = 2; k < 6; k++) if (itsCluMap & (1 << k)) ++nITS;
         //if (nITS < 3) return kFALSE;
         //itsNSigma = itsrsp.GetNumberOfSigmas(mom, itsSignal, AliPID::kPion, nITS, kTRUE);
         //okITS = ((TMath::Abs(itsNSigma)) <= fMaxITSband);
         if (nITS < 3) 
            okITS = kFALSE;
         else {
            itsNSigma = itsrsp.GetNumberOfSigmas(mom, itsSignal, AliPID::kPion, nITS, kTRUE);
            okITS = ((TMath::Abs(itsNSigma)) <= fMaxITSband);
            AliDebug(AliLog::kDebug + 2, Form("ITS nsigma = %f -- max = %f -- cut %s", itsNSigma, fMaxITSband, (okITS ? "passed" : "failed")));
         }
      } else {
         okITS = kTRUE;
      }
   } else {
      // if we are here, the track is surely bad
      okITS = kFALSE;
      okTPC = kFALSE;
      okTOF = kFALSE;
   }


   //okTrack[0] = okITS;
   //okTrack[1] = okTPC;
   //okTrack[2] = okTOF;
   //okTrack[3] = isTPC;
   //okTrack[4] = isITSSA;
   
   fOkTrack[0] = okITS;
   fOkTrack[1] = okTPC;
   fOkTrack[2] = okTOF;
   fOkTrack[3] = isTPC;
   fOkTrack[4] = isITSSA;

   //cout<<"##########################################"<<endl;
   //cout<<"isITSSA "<<isITSSA<< " isTPC "<<isTPC<<endl;
   //cout<<"##########################################"<<endl;


   //return okTrack;
}











//________________________________________________________________________
void AliAnalysisTaskSigma1385::Terminate(Option_t *)
{
   // Draw result to the screen
   // Called once at the end of the query
}



