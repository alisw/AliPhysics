/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TBits.h>
#include <TArrayD.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TSystem.h>
#include <TProfile.h>
#include <TKey.h>
#include "TVariableBinning.h"
#include "TLinearBinning.h"
#include "THnSparse.h"

#include <AliAnalysisManager.h>
#include <AliVEventHandler.h>
#include <AliVCluster.h>
#include <AliVParticle.h>
#include <AliLog.h>
#include "AliParticleContainer.h"
#include "AliAnalysisHelperJetTasks.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliAODEvent.h"
#include "AliMCEventHandler.h"
#include "AliGenPythiaEventHeader.h"

#include "AliGenPythiaEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliMCEvent.h"

#include "AliAnalysisTaskIPCalib.h"

ClassImp(AliAnalysisTaskIPCalib);

//_______________________________________________________________________________________
AliAnalysisTaskIPCalib::AliAnalysisTaskIPCalib() : AliAnalysisTaskEmcalJet(), fHistManager(), fReadMC(kFALSE), fCorrectRes(kFALSE)
{
  //
  // def ctor
  //

  for (int i = 0; i < 5; i++)
  {
    fCorrectionFactors[i] = 0x0;
  }
}

//_______________________________________________________________________________________
AliAnalysisTaskIPCalib::AliAnalysisTaskIPCalib(const char *name) : AliAnalysisTaskEmcalJet(name, kTRUE), fHistManager(name), fReadMC(kFALSE), fCorrectRes(kFALSE)
{
  //
  SetMakeGeneralHistograms(kTRUE);

  for (int i = 0; i < 5; i++)
  {
    fCorrectionFactors[i] = 0x0;
  }
}

//_______________________________________________________________________________________
AliAnalysisTaskIPCalib::~AliAnalysisTaskIPCalib()
{
  //
}

//_______________________________________________________________________________________
void AliAnalysisTaskIPCalib::UserCreateOutputObjects()
{
  //
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  AllocateTrackHistograms();

  TIter next(fHistManager.GetListOfHistograms());
  TObject *obj = 0;
  while ((obj = next()))
  {
    TH1 *hh = dynamic_cast<TH1 *>(obj);
    if (hh)
    {
      hh->Sumw2();
    }
    THnSparse *hn = dynamic_cast<THnSparse *>(obj);
    if (hn)
    {
      hn->Sumw2();
    }
    fOutput->Add(obj);
  }

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

//_______________________________________________________________________________________
void AliAnalysisTaskIPCalib::AllocateTrackHistograms()
{
  TString histname;
  TString histtitle;
  TString groupname;

  AliParticleContainer *partCont = 0;
  TIter next(&fParticleCollArray);
  while ((partCont = static_cast<AliParticleContainer *>(next())))
  {
    groupname = partCont->GetName();
    // Protect against creating the histograms twice
    if (fHistManager.FindObject(groupname))
    {
      AliWarning(TString::Format("%s: Found groupname %s in hist manager. The track containers will be filled into the same histograms.",
                                 GetName(), groupname.Data()));
      continue;
    }

    fHistManager.CreateHistoGroup(groupname);

    Int_t ttype[2] = {4, 9};

    for (Int_t i = 0; i < 2; i++)
    {
      histname = TString::Format("%s/histTrackPt_%d", groupname.Data(), ttype[i]);
      histtitle = TString::Format("%s;#it{p}_{T,track} (GeV/#it{c});counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, 500, 0, 100, "s");

      histname = TString::Format("%s/histTrackPhi_%d", groupname.Data(), ttype[i]);
      histtitle = TString::Format("%s;#it{#phi}_{track};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, 200, 0, TMath::TwoPi(), "s");

      histname = TString::Format("%s/histTrackEta_%d", groupname.Data(), ttype[i]);
      histtitle = TString::Format("%s;#it{#eta}_{track};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, 200, -1, 1, "s");

      histname = TString::Format("%s/histSd0Vms_%d", groupname.Data(), ttype[i]);
      histtitle = TString::Format("%s;1/(#it{p}^{2}sin^{3}#theta);#sigma_{d_{0}}^{2} (#mum^{2});counts",
                                  histname.Data());
      fHistManager.CreateTH2(histname, histtitle, 100, 0, 100, 1000, 0, 1e7, "s");

      histname = TString::Format("%s/histSd0Pms_%d", groupname.Data(), ttype[i]);
      histtitle = TString::Format("%s;1/(#it{p}_{T}^{2}sin#theta);#sigma_{d_{0}}^{2} (#mum^{2});counts",
                                  histname.Data());
      fHistManager.CreateTH2(histname, histtitle, 100, 0, 100, 1000, 0, 1e7, "s");

      histname = TString::Format("%s/histSd0PmsLog_%d", groupname.Data(), ttype[i]);
      histtitle = TString::Format("%s;Log_{10}(1/(#it{p}_{T}^{2}sin#theta));Log_{10}(#sigma_{d_{0}}^{2}) (#mum^{2});counts",
                                  histname.Data());
      fHistManager.CreateTH2(histname, histtitle, 500, -4, 3, 500, 2, 8, "s");

      histname = TString::Format("%s/histSd0Pt_%d", groupname.Data(), ttype[i]);
      histtitle = TString::Format("%s;#it{p}_{T};#sigma_{d_{0}} (#mum);counts",
                                  histname.Data());
      fHistManager.CreateTH2(histname, histtitle, 500, 0, 50, 1000, 0, 500, "s");

      histname = TString::Format("%s/hist1D0XY_%d", groupname.Data(), ttype[i]);
      histtitle = TString::Format("%s;d_{0} (mm);counts",
                                  histname.Data());
      fHistManager.CreateTH1(histname, histtitle, 2000, -10, 10, "s");

      Int_t bins = 15;
      TArrayD binedges(bins);

      Double_t from = 0.5;
      Double_t to = 50;

      for (int j = 0; j < bins; j++)
      {
        binedges[j] = pow(10, log10(from) + (log10(to) - log10(from)) / double(bins) * double(j));
      }

      TLinearBinning *xbinning = new TLinearBinning(400, -40, 40);
      TVariableBinning *ybinning = new TVariableBinning(binedges);
      TLinearBinning nbinning(100, 0, 100);
      const TBinning *histsmbinning[3] = {xbinning, ybinning, &nbinning};

      for (Int_t j = 1; j <= 6; j++)
      {
        histname = TString::Format("%s/histSIPPscat_%d_%d", groupname.Data(), ttype[i], j);
        histtitle = TString::Format("%s;S_{d_{0}};p(sin#theta)^{3/2};counts",
                                    histname.Data());
        fHistManager.CreateTH2(histname, histtitle, *xbinning, *ybinning, "s");

        histname = TString::Format("%s/ScorrTypeSigPscaNctr_%d_%d", groupname.Data(), ttype[i], j);
        histtitle = TString::Format("%s;S_{d_{0}};p(sin#theta)^{3/2};N_{PV}",
                                    histname.Data());
        fHistManager.CreateTHnSparse(histname, histtitle, 3, histsmbinning, "s");
      }

      histname = TString::Format("%s/hist2D0XY_%d", groupname.Data(), ttype[i]);
      histtitle = TString::Format("%s;d_{0} (mm);counts",
                                  histname.Data());
      fHistManager.CreateTH1(histname, histtitle, 2000, -10, 10, "s");

      histname = TString::Format("%s/hist3D0XY_%d", groupname.Data(), ttype[i]);
      histtitle = TString::Format("%s;d_{0} (mm);counts",
                                  histname.Data());
      fHistManager.CreateTH1(histname, histtitle, 2000, -10, 10, "s");

      histname = TString::Format("%s/histD0Z_%d", groupname.Data(), ttype[i]);
      histtitle = TString::Format("%s;d_{0} (mm);counts",
                                  histname.Data());
      fHistManager.CreateTH1(histname, histtitle, 2000, -10, 10, "s");

      histname = TString::Format("%s/histTrackNTPCcls_%d", groupname.Data(), ttype[i]);
      histtitle = TString::Format("%s;N_{TPC,cls} (mm);counts",
                                  histname.Data());
      fHistManager.CreateTH1(histname, histtitle, 200, 0, 200, "s");

      histname = TString::Format("%s/histTrackChi2perNDF_%d", groupname.Data(), ttype[i]);
      histtitle = TString::Format("%s;#chi^{2}/NDF;counts",
                                  histname.Data());
      fHistManager.CreateTH1(histname, histtitle, 200, 0, 20, "s");

      histname = TString::Format("%s/histTrackNITScls_%d", groupname.Data(), ttype[i]);
      histtitle = TString::Format("%s;N_{ITS,cls};counts",
                                  histname.Data());
      fHistManager.CreateTH1(histname, histtitle, 10, 0, 10, "s");
    }

    histname = TString::Format("%s/histNTracks", groupname.Data());
    histtitle = TString::Format("%s;number of tracks;events", histname.Data());
    if (fForceBeamType != kpp)
    {
      fHistManager.CreateTH1(histname, histtitle, 500, 0, 5000, "s");
    }
    else
    {
      fHistManager.CreateTH1(histname, histtitle, 200, 0, 200, "s");
    }
  }

  histname = "fHistSumNTracks";
  histtitle = TString::Format("%s;Sum of n tracks;events", histname.Data());
  if (fForceBeamType != kpp)
  {
    fHistManager.CreateTH1(histname, histtitle, 500, 0, 5000, "s");
  }
  else
  {
    fHistManager.CreateTH1(histname, histtitle, 200, 0, 200, "s");
  }
}

//_______________________________________________________________________________________
Bool_t AliAnalysisTaskIPCalib::FillHistograms()
{
  DoTrackLoop();

  return kTRUE;
}

//_______________________________________________________________________________________
void AliAnalysisTaskIPCalib::DoTrackLoop()
{
  AliAODVertex *vtx = (AliAODVertex *)InputEvent()->GetPrimaryVertex();
  Double_t Bz = InputEvent()->GetMagneticField();

  AliClusterContainer *clusCont = GetClusterContainer(0);

  TString histname;
  TString groupname;

  UInt_t sumAcceptedTracks = 0;

  TClonesArray *mcArray = 0x0;
  if (fReadMC)
  {
    mcArray = dynamic_cast<TClonesArray *>(((AliAODEvent *)InputEvent())->FindListObject(AliAODMCParticle::StdBranchName()));
  }

  AliParticleContainer *partCont = 0;
  TIter next(&fParticleCollArray);
  while ((partCont = static_cast<AliParticleContainer *>(next())))
  {
    groupname = partCont->GetName();
    UInt_t count = 0;
    for (auto part : partCont->accepted())
    {
      if (!part)
        continue;

      count++;

      AliAODTrack *track = dynamic_cast<AliAODTrack *>(part);

      Double_t d0z0[2], covd0z0[3];
      track->PropagateToDCA(vtx, Bz, 100, d0z0, covd0z0);
      Double_t d0 = TMath::Abs(d0z0[0]);
      Double_t sd0 = TMath::Sqrt(covd0z0[0]);

      Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z for the vertex cut
      track->GetImpactParameters(dca, cov);

      //TBits filterMap;
      //UInt_t fm = track->GetFilterMap();
      //filterMap.Set(16, &fm);
      //Int_t fb = (filterMap.TestBitNumber(4)) ? 4 : 9;
      Int_t fb = track->TestFilterBit(4) ? 4 : 9;

      histname = TString::Format("%s/histTrackNTPCcls_%d", groupname.Data(), fb);
      fHistManager.FillTH1(histname, track->GetTPCNcls());

      histname = TString::Format("%s/histTrackChi2perNDF_%d", groupname.Data(), fb);
      fHistManager.FillTH1(histname, track->Chi2perNDF());

      histname = TString::Format("%s/hist3D0XY_%d", groupname.Data(), fb);
      fHistManager.FillTH1(histname, 10 * d0z0[0]); // mm

      histname = TString::Format("%s/histTrackNITScls_%d", groupname.Data(), fb);
      fHistManager.FillTH1(histname, track->GetITSNcls());

      if (
          track->Pt() < 0.5 ||
          track->GetTPCNcls() < 80 ||
          track->Chi2perNDF() > 5 ||
          track->GetITSNcls() <= 1 ||
          TMath::Abs(d0z0[0]) > 1 ||
          TMath::Abs(d0z0[1]) > 2)
        continue;

      histname = TString::Format("%s/histTrackPt_%d", groupname.Data(), fb);
      fHistManager.FillTH1(histname, part->Pt());

      histname = TString::Format("%s/histTrackPhi_%d", groupname.Data(), fb);
      fHistManager.FillTH1(histname, part->Phi());

      histname = TString::Format("%s/histTrackEta_%d", groupname.Data(), fb);
      fHistManager.FillTH1(histname, part->Eta());

      Double_t vms = 1. / (TMath::Power(track->P(), 2) * TMath::Power(TMath::Sin(track->Theta()), 3));

      Double_t psc = track->P() * TMath::Power(TMath::Sin(track->Theta()), 1.5);

      Double_t pms = 1. / (TMath::Power(track->Pt(), 2) * TMath::Sin(track->Theta()));

      Double_t sip = d0z0[0] / sd0;

      if (fCorrectRes)
        sip = CorrectIPs(sip, psc, track->GetITSNcls());

      histname = TString::Format("%s/histSIPPscat_%d_%d", groupname.Data(), fb, track->GetITSNcls());
      fHistManager.FillTH2(histname, sip, psc);

      double point[3] = {sip, psc, static_cast<double>(vtx->GetNContributors())};
      histname = TString::Format("%s/ScorrTypeSigPscaNctr_%d_%d", groupname.Data(), fb, track->GetITSNcls());
      fHistManager.FillTHnSparse(histname, point);

      histname = TString::Format("%s/hist1D0XY_%d", groupname.Data(), fb);
      fHistManager.FillTH1(histname, 10 * d0z0[0]); // mm

      histname = TString::Format("%s/hist2D0XY_%d", groupname.Data(), fb);
      fHistManager.FillTH1(histname, 10 * dca[0]); // mm

      histname = TString::Format("%s/histD0Z_%d", groupname.Data(), fb);
      fHistManager.FillTH1(histname, 10 * d0z0[1]); // mm

      histname = TString::Format("%s/histSd0Vms_%d", groupname.Data(), fb);
      fHistManager.FillTH2(histname, vms, 1e8 * covd0z0[0]); // um^2

      histname = TString::Format("%s/histSd0Pms_%d", groupname.Data(), fb);
      fHistManager.FillTH2(histname, pms, 1e8 * covd0z0[0]); // um^2

      histname = TString::Format("%s/histSd0PmsLog_%d", groupname.Data(), fb);
      fHistManager.FillTH2(histname, TMath::Log10(pms), TMath::Log10(1e8 * covd0z0[0])); // um^2

      histname = TString::Format("%s/histSd0Pt_%d", groupname.Data(), fb);
      fHistManager.FillTH2(histname, track->Pt(), 1e4 * sd0); // um
    }

    sumAcceptedTracks += count;

    histname = TString::Format("%s/histNTracks", groupname.Data());
    fHistManager.FillTH1(histname, count);
  }

  histname = "fHistSumNTracks";
  fHistManager.FillTH1(histname, sumAcceptedTracks);
}

//_______________________________________________________________________________________
Double_t AliAnalysisTaskIPCalib::CorrectIPs(Double_t sIP, Double_t pScat, Int_t nITS)
{
  Double_t AlphaPull = 1.;
  if (pScat > 36.7)
    pScat = 36.7;

  for (int i = 0; i < 3; i++)
  {
    if (!fCorrectionFactors[i])
    {
      std::cout << "Correction function doesn't exists, Returning the original value of IPs" << std::endl;
      return sIP;
    }
  }

  for (int i = 3; i < 5; i++)
  {
    if (nITS > 4 && !fCorrectionFactors[i])
    {
      std::cout << "Correction function doesn't exists for ITS hits > 4, Returning the original value of IPs" << std::endl;
      return sIP;
    }
  }

  AlphaPull = fCorrectionFactors[nITS - 2]->Eval(pScat);

  Double_t CorrIPs = sIP / AlphaPull;

  return CorrIPs;
}

//_______________________________________________________________________________________
void AliAnalysisTaskIPCalib::ExecOnce()
{
  AliAnalysisTaskEmcalJet::ExecOnce();
}

/**
 * Run analysis code here, if needed.
 * It will be executed before FillHistograms().
 * If this function return kFALSE, FillHistograms() will *not*
 * be executed for the current event
 * @return Always kTRUE
 */
//_______________________________________________________________________________________
Bool_t AliAnalysisTaskIPCalib::Run()
{
  return kTRUE;
}

/**
 * This function is called once at the end of the analysis.
 */
//_______________________________________________________________________________________
void AliAnalysisTaskIPCalib::Terminate(Option_t *)
{
  //
}

/**
 * This function adds the task to the analysis manager. Often, this function is called
 * by an AddTask C macro. However, by compiling the code, it ensures that we do not
 * have to deal with difficulties caused by CINT.
 */
//_______________________________________________________________________________________
AliAnalysisTaskIPCalib *AliAnalysisTaskIPCalib::AddTaskIPCalib(const char *ntracks,
                                                                   TString pathToCorrFunc,
                                                                   const char *suffix)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalJetSample", "No analysis manager to connect to.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler *handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskEmcalJetSample", "This task requires an input event handler");
    return 0;
  }

  enum EDataType_t
  {
    kUnknown,
    kESD,
    kAOD
  };

  EDataType_t dataType = kUnknown;

  if (handler->InheritsFrom("AliESDInputHandler"))
  {
    dataType = kESD;
  }
  else if (handler->InheritsFrom("AliAODInputHandler"))
  {
    dataType = kAOD;
  }

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString trackName(ntracks);

  if (trackName == "usedefault")
  {
    if (dataType == kESD)
    {
      trackName = "Tracks";
    }
    else if (dataType == kAOD)
    {
      trackName = "tracks";
    }
    else
    {
      trackName = "";
    }
  }

  TString name("AliAnalysisTaskIPCalib");
  if (!trackName.IsNull())
  {
    name += "_";
    name += trackName;
  }
  if (strcmp(suffix, "") != 0)
  {
    name += "_";
    name += suffix;
  }

  AliAnalysisTaskIPCalib *ipTask = new AliAnalysisTaskIPCalib(name);
  //   ipTask->SetVzRange(-10,10);

  if (trackName == "mcparticles")
  {
    ipTask->AddMCParticleContainer(trackName);
  }
  else if (trackName == "tracks" || trackName == "Tracks")
  {
    ipTask->AddTrackContainer(trackName);
  }
  else if (!trackName.IsNull())
  {
    ipTask->AddParticleContainer(trackName);
  }

  if (!pathToCorrFunc.IsNull())
  {
    TFile *file = TFile::Open(pathToCorrFunc.Data());
    for (int i = 0; i < 5; i++)
    {
      if ((TH1D *)file->Get(Form("fno_%i", i)))
      {
        TF1 *CorrectionFunction = (TF1 *)file->Get(Form("fno_%i", i));
        ipTask->SetCorrectionFunction(CorrectionFunction, i);
      }
    }
  }

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(ipTask);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  TString contname(name);
  contname += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
                                                            TList::Class(), AliAnalysisManager::kOutputContainer,
                                                            Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput(ipTask, 0, cinput1);
  mgr->ConnectOutput(ipTask, 1, coutput1);

  return ipTask;
}
