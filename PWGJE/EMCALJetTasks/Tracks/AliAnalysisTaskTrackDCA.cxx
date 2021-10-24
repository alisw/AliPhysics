/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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
#include <map>
#include <string>
#include <vector>

#include <TArrayD.h>
#include <THashList.h>
#include <THistManager.h>
#include <TList.h>
#include <TMath.h>

#include "AliAnalysisUtils.h"
#include "AliInputEventHandler.h"
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"

#include "AliAnalysisTaskTrackDCA.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskTrackDCA)

using namespace PWGJE::EMCALJetTasks;

/**
 * Default constructor
 */
AliAnalysisTaskTrackDCA::AliAnalysisTaskTrackDCA():
    AliAnalysisTaskSE(),
    fAnalysisUtils(NULL),
    fStandardCuts(NULL),
    fHistos(NULL)
{
}

/**
 * Named constructor - setting name of the task and defining output container
 * \param name Name of the task
 */
AliAnalysisTaskTrackDCA::AliAnalysisTaskTrackDCA(const char *name):
    AliAnalysisTaskSE(name),
    fAnalysisUtils(NULL),
    fStandardCuts(NULL),
    fHistos(NULL)
{
  DefineOutput(1, TList::Class());
}

/**
 * Destructor, cleaning up
 */
AliAnalysisTaskTrackDCA::~AliAnalysisTaskTrackDCA() {
  if(fAnalysisUtils) delete fAnalysisUtils;
  if(fStandardCuts) delete fStandardCuts;
  if(fHistos) delete fHistos;
}

/**
 * Creating output histograms:
 * - Event counter for each trigger class
 * - DCA vs. \f$ p_{t} \f$ for each trigger class (min. bias, EMCAL gamma, EMCAL jet)
 */
void AliAnalysisTaskTrackDCA::UserCreateOutputObjects(){
  fAnalysisUtils = new AliAnalysisUtils;

  std::string triggerclasses[] = {"MinBiasINT7", "MinBiasINT8", "EGAINT7", "EGAINT8", "EJEINT7", "EJEINT8"};
  TArrayD ptbinning, dcarbinning, dcazbinning;
  CreatePtBinning(ptbinning);
  CreateLinearBinning(dcarbinning, 200, -0.1, 0.1);
  CreateLinearBinning(dcazbinning, 500, -5, 5);
  fHistos = new THistManager("trackDCA");
  for(std::string *trgiter = triggerclasses; trgiter < triggerclasses + sizeof(triggerclasses)/sizeof(std::string); trgiter++){
    fHistos->CreateTH1(Form("hEvents%s", trgiter->c_str()), Form("Event counter for trigger class %s", trgiter->c_str()), 1, 0.5, 1.5);
    fHistos->CreateTH2(Form("hDCAr%s", trgiter->c_str()), Form("DCA distribution vs. p_{t} for trigger class %s", trgiter->c_str()), ptbinning, dcarbinning);
    fHistos->CreateTH2(Form("hDCAz%s", trgiter->c_str()), Form("DCA distribution vs. p_{t} for trigger class %s", trgiter->c_str()), ptbinning, dcazbinning);
  }

  PostData(1, fHistos->GetListOfHistograms());

}

/**
 * Event loop
 * - Apply event selection (triggers and vertex)
 * - Loop over tracks and apply track selection
 * - Fill DCA r and z distribution for selected tracks
 *
 * \param Not needed
 */
void AliAnalysisTaskTrackDCA::UserExec(Option_t *){
  AliESDEvent *esdev = dynamic_cast<AliESDEvent *>(fInputEvent);
  if(!esdev) return;

  // Check triggers - split them up in INT7 (V0) and INT8 (T0)
  std::vector<std::string> triggers;
  if(fInputHandler->IsEventSelected() & AliVEvent::kINT7) triggers.push_back("MinBiasINT7");
  if(fInputHandler->IsEventSelected() & AliVEvent::kINT8) triggers.push_back("MinBiasINT8");
  TString triggerstring = esdev->GetFiredTriggerClasses();
  if(triggerstring.Contains("EGA")){
    if(triggerstring.Contains("CEMC7")) triggers.push_back("EGAINT7");
    else if(triggerstring.Contains("CEMC8")) triggers.push_back("EGAINT8");
  } else if(triggerstring.Contains("EJE")){
    if(triggerstring.Contains("CEMC7")) triggers.push_back("EJEINT7");
    else if(triggerstring.Contains("CEMC8")) triggers.push_back("EJEINT8");
  }
  if(!triggers.size()) return;

  // Apply event selection via AliAnalysisUtils
  if(!fAnalysisUtils->IsVertexSelected2013pA(fInputEvent)) return;
  const AliESDVertex *primVertex = esdev->GetPrimaryVertex();
  if(TMath::Abs(primVertex->GetZ()) > 10.) return;

  for(std::vector<std::string>::iterator trgiter = triggers.begin(); trgiter != triggers.end(); ++trgiter)
    fHistos->FillTH1(Form("hEvents%s", trgiter->c_str()), 1.);
  /*
   * In case AliRoot will at some point support c++-11 (in the far future)
   * for(auto trgit : triggers) fHistos->FillTH1(Form("hEvents%s", trgit.c_str()), 1.);
   */

  // Loop over tracks, apply track selection, fill histograms
  const AliESDtrack *track(NULL);
  for(Int_t ipart = 0; ipart < esdev->GetNumberOfTracks(); ipart++){
    track = esdev->GetTrack(ipart);
    AliESDtrack copytrack(*track);
    if(!fStandardCuts->IsSelected(&copytrack)) continue;
    Float_t dcaR, dcaZ;
    Double_t pt = TMath::Abs(copytrack.Pt());
    copytrack.GetImpactParameters(dcaR, dcaZ);
    for(std::vector<std::string>::iterator trgit = triggers.begin(); trgit != triggers.end(); ++trgit){
      fHistos->FillTH1(Form("hDCAr%s", trgit->c_str()), pt, dcaR);
      fHistos->FillTH1(Form("hDCAz%s", trgit->c_str()), pt, dcaZ);
    }
    /*
     * In case AliRoot will at some point support c++-11 (in the far future)
     *
     * for(auto trgit : triggers){
     *   fHistos->FillTH1(Form("hDCAr%s", trgit.c_str()), dcaR);
     *   fHistos->FillTH1(Form("hDCAz%s", trgit.c_str()), dcaZ);
     * }
     */
  }

  PostData(1, fHistos->GetListOfHistograms());
}

/**
 * Create \f$ p_{t} \f$ binning used in the \f$ R_{AA} \f$ analysis:
 *
 * Definitions are:
 * - from 0.15 to 1 GeV/c: 0.05 GeV/c bins
 * - from 1 to 2 GeV/c: 0.1 GeV/c bins
 * - from 2 to 4 GeV/c: 0.2 GeV/c bins
 * - from 4 to 7 GeV/c: 0.5 GeV/c bins
 * - from 7 to 16 GeV/c: 1 GeV/c bins
 * - from 16 to 36 GeV/c: 2 GeV/c bins
 * - from 36 to 40 GeV/c: 4 GeV/c bins
 * - from 40 to 50 GeV/c: 5 GeV/c bins
 * - from 50 to 100 GeV/c: 10 GeV/c bins
 *
 * \param binning: Array where to store the results
 */
void AliAnalysisTaskTrackDCA::CreatePtBinning(TArrayD& binning) const {
  std::vector<double> mybinning;
  std::map<double,double> definitions;
  definitions.insert(std::pair<double, double>(1, 0.05));
  definitions.insert(std::pair<double, double>(2, 0.1));
  definitions.insert(std::pair<double, double>(4, 0.2));
  definitions.insert(std::pair<double, double>(7, 0.5));
  definitions.insert(std::pair<double, double>(16, 1));
  definitions.insert(std::pair<double, double>(36, 2));
  definitions.insert(std::pair<double, double>(40, 4));
  definitions.insert(std::pair<double, double>(50, 5));
  definitions.insert(std::pair<double, double>(100, 10));
  definitions.insert(std::pair<double, double>(200, 20));
  double currentval = 0.;
  mybinning.push_back(currentval);
  for(std::map<double,double>::iterator id = definitions.begin(); id != definitions.end(); ++id){
    double limit = id->first, binwidth = id->second;
    while(currentval < limit){
      currentval += binwidth;
      mybinning.push_back(currentval);
    }
  }
  binning.Set(mybinning.size());
  int ib = 0;
  for(std::vector<double>::iterator it = mybinning.begin(); it != mybinning.end(); ++it)
    binning[ib++] = *it;
}

/**
 * Create any kind of linear binning from given ranges and stores it in the binning array.
 *
 * \param binning output array
 * \param nbins Number of bins
 * \param min lower range
 * \param max upper range
 */
void AliAnalysisTaskTrackDCA::CreateLinearBinning(TArrayD& binning, int nbins, double min, double max) const {
  double binwidth = (max-min)/static_cast<double>(nbins);
  binning.Set(nbins+1);
  binning[0] = min;
  double currentlimit = min + binwidth;
  for(int ibin = 0; ibin < nbins; ibin++){
    binning[ibin+1] = currentlimit;
    currentlimit += binwidth;
  }
}
