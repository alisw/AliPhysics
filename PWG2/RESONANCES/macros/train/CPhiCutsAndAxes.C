//
// Initialization of some example cuts of common use
//

#if !defined(__CINT__) || defined(__MAKECINT__)

   #include "AliRsnCutPrimaryVertex.h"
   #include "AliRsnCutTrackQuality.h"
   #include "AliRsnCutPIDTPC.h"
   #include "AliRsnCutPIDTOF.h"
   #include "AliRsnCutPIDITS.h"
   #include "AliRsnCutValue.h"
   
   #include "AliRsnValue.h"
   #include "AliRsnPairDef.h"
   
#endif

//__________________________________________________________________________________________________
//
// Track quality for ITS standalone:
// this cut is used to select tracks of good quality, irrespective of the PID.
// When adding status flags, the second argument tells if each considered flag
// must be active or not in the track status, since the ITS-SA tracks need that
// some of them are OFF (e.g.: kTPCin)
//
AliRsnCutTrackQuality* TrackQualityITS()
{
   AliRsnCutTrackQuality *cutQualityITS = new AliRsnCutTrackQuality("cutQualityITS");
   cutQualityITS->AddStatusFlag(AliESDtrack::kITSin    , kTRUE);
   cutQualityITS->AddStatusFlag(AliESDtrack::kTPCin    , kFALSE);
   cutQualityITS->AddStatusFlag(AliESDtrack::kITSrefit , kTRUE);
   cutQualityITS->AddStatusFlag(AliESDtrack::kTPCrefit , kFALSE);
   cutQualityITS->AddStatusFlag(AliESDtrack::kITSpureSA, kFALSE);
   cutQualityITS->AddStatusFlag(AliESDtrack::kITSpid   , kTRUE);
   cutQualityITS->SetPtRange(0.15, 1E+20);
   cutQualityITS->SetEtaRange(-0.8, 0.8);
   cutQualityITS->SetDCARPtFormula("0.0595+0.0182/pt^1.55");
   cutQualityITS->SetDCAZmax(2.0);
   cutQualityITS->SetSPDminNClusters(1);
   cutQualityITS->SetITSminNClusters(4);
   cutQualityITS->SetITSmaxChi2(2.0);
   cutQualityITS->SetTPCminNClusters(0);
   cutQualityITS->SetTPCmaxChi2(1E+10);
   cutQualityITS->SetRejectKinkDaughters();
   
   return cutQualityITS;
}
   
//__________________________________________________________________________________________________
//
// Track quality for TPC+ITS:
// works exactly like the one above, but has settings for selecting TPC+ITS tracks
// in this case, the flags required are all necessary, so here the procedure is simpler
//

AliRsnCutTrackQuality *TrackQualityTPC()
{
   AliRsnCutTrackQuality *cutQualityTPC = new AliRsnCutTrackQuality("cutQualityTPC");
   cutQualityTPC->AddStatusFlag(AliESDtrack::kTPCin   , kTRUE);
   cutQualityTPC->AddStatusFlag(AliESDtrack::kTPCrefit, kTRUE);
   cutQualityTPC->AddStatusFlag(AliESDtrack::kITSrefit, kTRUE);
   cutQualityTPC->SetPtRange(0.15, 1E+20);
   cutQualityTPC->SetEtaRange(-0.8, 0.8);
   cutQualityTPC->SetDCARPtFormula("0.0182+0.0350/pt^1.01");
   cutQualityTPC->SetDCAZmax(2.0);
   cutQualityTPC->SetSPDminNClusters(1);
   cutQualityTPC->SetITSminNClusters(0);
   cutQualityTPC->SetITSmaxChi2(1E+20);
   cutQualityTPC->SetTPCminNClusters(70);
   cutQualityTPC->SetTPCmaxChi2(4.0);
   cutQualityTPC->SetRejectKinkDaughters();
   
   return cutQualityTPC;
}
   
//__________________________________________________________________________________________________
//
// The ITS PID is done with a 3sigma band in the dE/dx vs. Bethe-Bloch comparison,
// in the whole momentum range, with a set of parameters which are different in
// data and MonteCarlo, so we need to know on what we are running.
//
AliRsnCutPIDITS *PIDITS
(Bool_t isMC, AliPID::EParticleType pid, Double_t pmin, Double_t pmax, Double_t nsigma, const char *name)
{
   AliRsnCutPIDITS *cutPIDITS = new AliRsnCutPIDITS(name, pid, -nsigma, nsigma);
   cutPIDITS->SetMC(isMC);
   cutPIDITS->SetMomentumRange(pmin, pmax);
   
   return cutPIDITS;
}
   
//__________________________________________________________________________________________________  
//
// For the TPC PID we define two zones: one below 350 MeV/c in momentum (at the TPC barrel)
// where we apply a 5sigma cut, and another above 350 MeV/c where the cut is tightened 
// to the range of 3sigma.
// Even in this case we have different parameters for data and MC, but here we have to set
// them manually (they are passed to the AliESDpid object which is in the cut)
//
AliRsnCutPIDTPC *PIDTPC
(Bool_t isMC, AliPID::EParticleType pid, Double_t pmin, Double_t pmax, Double_t nsigma, const char *name)
{   
   AliRsnCutPIDTPC *cutPIDTPC = new AliRsnCutPIDTPC(name, pid, -nsigma, nsigma);
   
   // BB parameterization depends on data sample (MC, data)
   // the momentum range is passed and tracks outside it are rejected
   Double_t bbPar[5];
   if (isMC) {
      bbPar[0] = 2.15898 / 50.0;
      bbPar[1] = 1.75295E1;
      bbPar[2] = 3.40030E-9;
      bbPar[3] = 1.96178;
      bbPar[4] = 3.91720;
   } else {
      bbPar[0] = 1.41543 / 50.0;
      bbPar[1] = 2.63394E1;
      bbPar[2] = 5.0411E-11;
      bbPar[3] = 2.12543;
      bbPar[4] = 4.88663;
   }
   cutPIDTPC->SetBBParam(bbPar);
   cutPIDTPC->SetMomentumRange(pmin, pmax);
   cutPIDTPC->SetRejectOutside(kTRUE);
   
   return cutPIDTPC;
}
   
//__________________________________________________________________________________________________
//
// The TOF PID is simply a 3sigma cut.
// Since it is not doing a comparison with a reference function,
// we don't need to know if we are on data or MC.
// It is important to choose if this cut must reject tracks not matched in TOF.
// Usually, if TPC pid is used, we can accept them, otherwise we must reject.
//
AliRsnCutPIDTOF *PIDTOF
(AliPID::EParticleType pid, Double_t nsigma, Bool_t rejectUnmatched, const char *name)
{
   AliRsnCutPIDTOF *cutPIDTOF = new AliRsnCutPIDTOF(name, pid, -nsigma, nsigma);
   cutPIDTOF->SetRejectUnmatched(rejectUnmatched);
   
   return cutPIDTOF;
}

//__________________________________________________________________________________________________ 
//
// A cut is applied on the rapidity of the pair.
// This is a vaolue which can be computed for filling histograms,
// and all such values can be cut using the generic AliRsnCutValue class.
// Since rapidity needs a mass hypothesis, we need to add a support object
// which is a reference pair definition (we can use any of those defined above)
//
AliRsnCutValue *RapidityRange
(AliRsnPairDef *def, Double_t min, Double_t max, const char *name = "cutY")
{   
   AliRsnCutValue *cutRapidity = new AliRsnCutValue(name, AliRsnValue::kPairY, min, max);
   cutRapidity->GetValueObj()->SetSupportObject(def);
   
   return cutRapidity;
}

//__________________________________________________________________________________________________  
//
// Define all axes
//
AliRsnValue *AxisIM(Double_t min = 0.9, Double_t max = 1.4, Double_t step = 0.001)
{
   AliRsnValue *axis = new AliRsnValue("IM", AliRsnValue::kPairInvMass, min, max, step);
   return axis;
}

AliRsnValue *AxisRes(Double_t min = -0.5, Double_t max = 0.5, Double_t step = 0.001)
{
   AliRsnValue *axis = new AliRsnValue("IM", AliRsnValue::kPairInvMassRes, min, max, step);
   return axis;
}

AliRsnValue *AxisPt(Double_t min = 0.0, Double_t max = 5.0, Double_t step = 0.1)
{
   AliRsnValue *axis = new AliRsnValue("PT", AliRsnValue::kPairPt, min, max, step);
   return axis;
}

AliRsnValue *AxisY(Double_t min = -1.5, Double_t max = 1.5, Double_t step = 0.1)
{
   AliRsnValue *axis = new AliRsnValue("Y", AliRsnValue::kPairPt, min, max, step);
   return axis;
}

AliRsnValue *AxisMultSPD()
{
   Double_t mult[] = { 0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9., 10., 11., 12., 13.,  14.,  15.,  16.,  17.,  18.,  19., 
                      20., 21., 22., 23., 24., 25., 30., 35., 40., 50., 60., 70., 80., 90., 100., 120., 140., 160., 180., 200., 500.};
   Int_t    nmult  = sizeof(mult) / sizeof(mult[0]);
   
   AliRsnValue *axis = new AliRsnValue("MSPD", AliRsnValue::kEventMultSPD, nmult, mult);
   return axis;
}

AliRsnValue *AxisMultMC()
{
   Double_t mult[] = { 0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9., 10., 11., 12., 13.,  14.,  15.,  16.,  17.,  18.,  19., 
                      20., 21., 22., 23., 24., 25., 30., 35., 40., 50., 60., 70., 80., 90., 100., 120., 140., 160., 180., 200., 500.};
   Int_t    nmult  = sizeof(mult) / sizeof(mult[0]);
   
   AliRsnValue *axis = new AliRsnValue("MSPD", AliRsnValue::kEventMultMC, nmult, mult);
   return axis;
}
