#include "AliMUONBusPatchEvolution.h"

/// \class AliMUONBusPatchEvolution
///
/// This utility class can be used to perform two operations on  the output of the MCHBPEVO DA
/// (which is, on purpose, limited to the strict minimum) : extract alarm information
/// from it, and augment it with derived information.
///
/// The GetListOfFaultyBusPatches() can be used to get the list of bus patches with an occupancy above
/// a given threshold.
///
/// The Augment() method goes from the list of histograms (in the form of an AliMergeableCollection)
/// describing the time evolution of the bus patches hit count, this class
/// adds histograms to get the bus patches occupancy (=hit count / nevents / npads).
/// Then the occupancy is also computed at different levels of the detector : detection
/// element, DDL, chamber, station, in order to allow a global view of issues (if any)
///
/// \author Laurent Aphecetche, Subatech

///\cond CLASSIMP
ClassImp(AliMUONBusPatchEvolution)
///\endcond

#include "AliCDBManager.h"
#include "AliDAQ.h"
#include "AliLog.h"
#include "AliMergeableCollection.h"
#include "AliMpBusPatch.h"
#include "AliMpCDB.h"
#include "AliMpDCSNamer.h"
#include "AliMpDDL.h"
#include "AliMpDDLStore.h"
#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpDetElement.h"
#include "TH1F.h"
#include "TList.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#include <cassert>
#include <iostream>
#include <set>

  //_________________________________________________________________________________________________
  AliMUONBusPatchEvolution::AliMUONBusPatchEvolution(AliMergeableCollection& hc)
  : TObject(), fBPEVO(hc), fNofPads()
{
  ComputeNumberOfPads();
}

//_________________________________________________________________________________________________
AliMUONBusPatchEvolution::AliMUONBusPatchEvolution(AliMergeableCollection& hc,
                                                   const std::map<int, int>& nofPadsPerBusPatch)
  : TObject(), fBPEVO(hc), fNofPads(nofPadsPerBusPatch)
{
}

//_________________________________________________________________________________________________
void AliMUONBusPatchEvolution::Augment()
{
  /// From the list of histograms (in the form of an AliMergeableCollection)
  /// describing the time evolution of the bus patches hit count, this class
  /// adds histograms to get the bus patches occupancy (=hit count / nevents / npads).
  /// Then the occupancy is also computed at different levels of the detector : detection
  /// element, DDL, chamber, station, in order to allow a global view of issues (if any)
  ///

  if (FillNumberOfPads()) {
    std::vector<int> timeResolutions;

    GetTimeResolutions(timeResolutions);

    AliInfo(Form("Number of time resolutions found : %lu", timeResolutions.size()));

    ShrinkTimeAxis();

    for (std::vector<int>::size_type i = 0; i < timeResolutions.size(); ++i) {
      AliInfo(Form("TimeResolution: %ds", timeResolutions[i]));
      const int tr = timeResolutions[i];
      GroupByDE(tr);
      GroupByDDL(tr);
      GroupByChamber(tr);
      GroupByStation(tr);
    }
    Normalize();
  }
}

//_________________________________________________________________________________________________
void AliMUONBusPatchEvolution::ComputeNumberOfPads()
{
  /// Compute the number of pads in each bus patch
  /// from the mapping

  if (!fNofPads.empty())
    return;

  AliCDBManager* cdbm = AliCDBManager::Instance();

  if (!cdbm->IsDefaultStorageSet()) {
    Long_t id, size, flags, modtime;
    TDatime now;

    TString testPath;

    testPath.Form("/cvmfs/alice-ocdb.cern.ch/calibration/data/%d/OCDB", now.GetYear());

    if (!gSystem->GetPathInfo(testPath.Data(), &id, &size, &flags, &modtime)) {
      if (flags & 0x1) {
        cdbm->SetDefaultStorage(Form("local://%s", testPath.Data()));
      }
    }
    else {
      cdbm->SetDefaultStorage("raw://");
    }
  }

  cdbm->SetRun(0);

  AliMpCDB::LoadAll();

  TIter next(AliMpDDLStore::Instance()->CreateBusPatchIterator());
  AliMpBusPatch* bp;

  Int_t total(0);

  while ((bp = static_cast<AliMpBusPatch*>(next()))) {
    Int_t npads(0);

    Int_t detElemId = AliMpDDLStore::Instance()->GetDEfromBus(bp->GetId());

    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);

    for (Int_t i = 0; i < bp->GetNofManus(); ++i) {
      Int_t manuId = bp->GetManuId(i);

      npads += de->NofChannelsInManu(manuId);
    }

    fNofPads[bp->GetId()] = npads;

    total += npads;
  }

  assert(total == 1064008);
}

//_________________________________________________________________________________________________
TH1* AliMUONBusPatchEvolution::ExpandTimeAxis(const TH1& h, Int_t expansionTime, Int_t timeResolution)
{
  /// Extend (by increasing the number of bins) the time range of h by nbins
  /// if expansion < 0 expansion is from before the start of x-axis, otherwise from the end.

  if (timeResolution < 0) {
    timeResolution = GetTimeResolution(h);
  }

  if (timeResolution <= 0)
    return 0x0;

  const TAxis* timeAxis = h.GetXaxis();

  Int_t extraBins = TMath::Nint(TMath::Abs(expansionTime) / timeResolution);

  Int_t nbins = timeAxis->GetNbins() + extraBins;
  Double_t xmin = timeAxis->GetXmin();
  Double_t xmax = timeAxis->GetXmax();

  TTimeStamp origin;
  GetTimeOffset(h, origin);

  Int_t makeRoomBefore = 0;

  if (expansionTime < 0) {
    makeRoomBefore = 1;
    xmin += expansionTime;
  }
  else {
    xmax += expansionTime;
  }

  TH1* hnew = new TH1F(h.GetName(), h.GetTitle(), nbins, xmin, xmax);
  hnew->SetDirectory(0);
  hnew->GetXaxis()->SetTimeDisplay(1);
  hnew->GetXaxis()->SetTimeFormat(timeAxis->GetTimeFormat());

  for (Int_t i = 1; i <= timeAxis->GetNbins(); ++i) {
    hnew->SetBinContent(i + extraBins*makeRoomBefore, h.GetBinContent(i));
    hnew->SetBinError(i + extraBins*makeRoomBefore, h.GetBinError(i));
  }

  if (expansionTime < 0) {
    /// Must update the time offset
    TTimeStamp origin;
    GetTimeOffset(h, origin);
    UInt_t ot = origin.GetSec();
    ot += expansionTime;
    origin.SetSec(ot);
    hnew->GetXaxis()->SetTimeOffset(origin.GetSec(), "gmt");
    GetTimeOffset(*hnew, origin);
  }

  return hnew;
}

//_________________________________________________________________________________________________
Bool_t AliMUONBusPatchEvolution::FillNumberOfPads()
{
  /// Fill histograms with the number of pads per element
  /// (BP, DE, DDL, chamber, station)

  assert(fNofPads.size() == 888);

  if (fBPEVO.Histo("/BUSPATCH/NPADS/BP0001")) {
    // work already done
    AliWarning("Already done. Not re-doing it...");
    return kFALSE;
  }

  Int_t total(0);

  AliMpDCSNamer dcs("TRACKER");

  std::map<int, int>::const_iterator it;

  for (it = fNofPads.begin(); it != fNofPads.end(); ++it) {
    int busPatchId = it->first;
    int npads = it->second;

    TH1* h = new TH1F(Form("BP%04d", busPatchId), "number of pads", 1, 0, 1);

    Int_t detElemId = AliMpDDLStore::Instance()->GetDEfromBus(busPatchId);

    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);

    h->Fill(0.0, 1.0 * npads);

    TH1* hde = fBPEVO.Histo(Form("/DE/NPADS/DE%04d", detElemId));
    if (!hde) {
      hde = new TH1F(Form("DE%04d", detElemId), "number of pads", 1, 0, 1);
      fBPEVO.Adopt("/DE/NPADS/", hde);
    }
    hde->Fill(0.0, 1.0 * npads);

    Int_t ddlId = de->GetDdlId() + AliDAQ::DdlIDOffset("MUONTRK");

    TH1* hddl = fBPEVO.Histo(Form("/DDL/NPADS/DDL%d", ddlId));
    if (!hddl) {
      hddl = new TH1F(Form("DDL%d", ddlId), "number of pads", 1, 0, 1);
      fBPEVO.Adopt("/DDL/NPADS/", hddl);
    }
    hddl->Fill(0.0, 1.0 * npads);

    Int_t chamberId = 1 + AliMpDEManager::GetChamberId(detElemId);
    Int_t stationId = 1 + AliMpDEManager::GetChamberId(detElemId) / 2;

    TH1* hchamberSide(0x0);
    TH1* hstationSide(0x0);

    if (dcs.DCSAliasName(detElemId).Contains("Left")) {
      hchamberSide = fBPEVO.Histo(Form("/CHAMBER/NPADS/CHAMBER%dLEFT", chamberId));
      if (!hchamberSide) {
        hchamberSide = new TH1F(Form("CHAMBER%dLEFT", chamberId), "number of pads", 1, 0, 1);
        fBPEVO.Adopt("/CHAMBER/NPADS/", hchamberSide);
      }
      hstationSide = fBPEVO.Histo(Form("/STATION/NPADS/STATION%dLEFT", stationId));
      if (!hstationSide) {
        hstationSide = new TH1F(Form("STATION%dLEFT", stationId), "number of pads", 1, 0, 1);
        fBPEVO.Adopt("/STATION/NPADS/", hstationSide);
      }
    }
    else {
      if (dcs.DCSAliasName(detElemId).Contains("Right")) {
        hchamberSide = fBPEVO.Histo(Form("/CHAMBER/NPADS/CHAMBER%dRIGHT", chamberId));
        if (!hchamberSide) {
          hchamberSide = new TH1F(Form("CHAMBER%dRIGHT", chamberId), "number of pads", 1, 0, 1);
          fBPEVO.Adopt("/CHAMBER/NPADS/", hchamberSide);
        }
        hstationSide = fBPEVO.Histo(Form("/STATION/NPADS/STATION%dRIGHT", stationId));
        if (!hstationSide) {
          hstationSide = new TH1F(Form("STATION%dRIGHT", stationId), "number of pads", 1, 0, 1);
          fBPEVO.Adopt("/STATION/NPADS/", hstationSide);
        }
      }
    }

    hchamberSide->Fill(0.0, 1.0 * npads);
    hstationSide->Fill(0.0, 1.0 * npads);

    TH1* hchamber = fBPEVO.Histo(Form("/CHAMBER/NPADS/CHAMBER%d", chamberId));
    if (!hchamber) {
      hchamber = new TH1F(Form("CHAMBER%d", chamberId), "number of pads", 1, 0, 1);
      fBPEVO.Adopt("/CHAMBER/NPADS/", hchamber);
    }
    hchamber->Fill(0.0, 1.0 * npads);

    TH1* hstation = fBPEVO.Histo(Form("/STATION/NPADS/STATION%d", stationId));
    if (!hstation) {
      hstation = new TH1F(Form("STATION%d", stationId), "number of pads", 1, 0, 1);
      fBPEVO.Adopt("/STATION/NPADS/", hstation);
    }
    hstation->Fill(0.0, 1.0 * npads);

    total += npads;

    fBPEVO.Adopt("/BUSPATCH/NPADS/", h);
  }

  return kTRUE;
}

//_________________________________________________________________________________________________
Bool_t AliMUONBusPatchEvolution::GetTimeOffset(const TH1& h, TTimeStamp& origin)
{
  /// Assuming histogram h has a time axis, fills origin with its time offset and returns
  /// kTRUE. Otherwise returns kFALSE
  /// Assumes that the time offset is indicated in GMT

  const TAxis* x = h.GetXaxis();

  if (x->GetTimeDisplay()) {
    TString tf = x->GetTimeFormat();
    Int_t ix = tf.Index('F');
    Int_t year, month, day;
    Int_t hour, minute, second;
    sscanf(tf(ix + 1, tf.Length() - ix - 1).Data(), "%4d-%02d-%02d %02d:%02d:%02d", &year, &month, &day, &hour, &minute,
           &second);
    origin.Set(year, month, day, hour, minute, second, 0, kTRUE, 0);
    return kTRUE;
  }
  else {
    return kFALSE;
  }
}

//_________________________________________________________________________________________________
int AliMUONBusPatchEvolution::GetTimeResolution(const TH1& h)
{
  // we assume the title ends with "xxx s bins"
  // and we try to find xxx ...

  TString title = h.GetTitle();

  int timeResolution(-1);

  TObjArray* a = title.Tokenize(" ");
  TString last = static_cast<TObjString*>(a->Last())->String();
  TString m = static_cast<TObjString*>(a->At(a->GetLast() - 1))->String();

  if (last != "bins" || m != "s") {
    std::cerr << "Histogram title does not match the expected pattern (... xxx s bins) : " << std::endl
              << h.GetTitle() << std::endl;
  }
  else {
    timeResolution = static_cast<TObjString*>(a->At(a->GetLast() - 2))->String().Atoi();
  }
  delete a;

  return timeResolution;
}

//_________________________________________________________________________________________________
void AliMUONBusPatchEvolution::GetTimeResolutions(std::vector<int>& timeResolutions)
{
  TIter next(fBPEVO.CreateIterator());

  TH1F* h(0x0);

  while ((h = static_cast<TH1F*>(next()))) {
    if (TString(h->GetName()).BeginsWith("Nevents")) {
      int tr = GetTimeResolution(*h);
      if (tr > 0) {
        timeResolutions.push_back(tr);
      }
    }
  }
}

//______________________________________________________________________________
Bool_t AliMUONBusPatchEvolution::GetFaultyBusPatches(int timeResolution, int requiredEvents, float occupancyThreshold,
                                                     std::map<int, double>& faultyBusPatchOccupancies)
{
  /// Loop over all bus patches and compute their occupancy in the last nevents
  /// Return kFALSE if the number of required events was not reached
  ///

  // find how many bins should be considered, by finding in the Nevents histogram
  // how many of the latest bins are required to get an integral >= requiredEvents
  TH1* hnevents = fBPEVO.Histo(Form("Nevents%ds", timeResolution));

  int lastbin(0);

  for (int i = hnevents->GetXaxis()->GetNbins() - 1; i > 0 && lastbin == 0; --i) {
    if (hnevents->GetBinContent(i) > 1)
      lastbin = i;
  }

  float nevents(0);

  int firstbin = lastbin;

  for (firstbin = lastbin; firstbin > 0; --firstbin) {
    nevents += hnevents->GetBinContent(firstbin);
    if (nevents > requiredEvents)
      break;
  }

  if (nevents < requiredEvents) {
    return kFALSE;
  }

  double dnevents = hnevents->Integral(firstbin, lastbin);

  for (std::map<int, int>::const_iterator it = fNofPads.begin(); it != fNofPads.end(); ++it) {
    const int& buspatchId = it->first;
    const int& npads = it->second;

    TString bpName = Form("BP%04d", buspatchId);

    TH1* hbp = fBPEVO.Histo(Form("/BUSPATCH/HITS/%ds/%s", timeResolution, bpName.Data()));

    double occ = hbp->Integral(firstbin, lastbin) / dnevents / npads;

    if (occ > occupancyThreshold) {
      faultyBusPatchOccupancies[buspatchId] = occ;
    }
  }

  return kTRUE;
}

//_________________________________________________________________________________________________
void AliMUONBusPatchEvolution::GroupByStation(int timeResolution)
{
  /// From the histograms at chamber level, make the corresponding ones at station level

  int station(1);

  for (Int_t ich = 1; ich < 10; ich += 2) {
    TH1* h = fBPEVO.Histo(Form("/CHAMBER/HITS/%ds/CHAMBER%d", timeResolution, ich));
    TH1* h1 = fBPEVO.Histo(Form("/CHAMBER/HITS/%ds/CHAMBER%d", timeResolution, ich + 1));

    TH1* hstation = static_cast<TH1*>(h->Clone(Form("STATION%d", station)));

    hstation->Add(h1);
    fBPEVO.Adopt(Form("/STATION/HITS/%ds", timeResolution), hstation);

    h = fBPEVO.Histo(Form("/CHAMBER/HITS/%ds/CHAMBER%dLEFT", timeResolution, ich));
    h1 = fBPEVO.Histo(Form("/CHAMBER/HITS/%ds/CHAMBER%dLEFT", timeResolution, ich + 1));

    hstation = static_cast<TH1*>(h->Clone(Form("STATION%dLEFT", station)));

    hstation->Add(h1);
    fBPEVO.Adopt(Form("/STATION/HITS/%ds", timeResolution), hstation);

    h = fBPEVO.Histo(Form("/CHAMBER/HITS/%ds/CHAMBER%dRIGHT", timeResolution, ich));
    h1 = fBPEVO.Histo(Form("/CHAMBER/HITS/%ds/CHAMBER%dRIGHT", timeResolution, ich + 1));

    hstation = static_cast<TH1*>(h->Clone(Form("STATION%dRIGHT", station)));

    hstation->Add(h1);
    fBPEVO.Adopt(Form("/STATION/HITS/%ds", timeResolution), hstation);

    ++station;
  }
}

//_________________________________________________________________________________________________
void AliMUONBusPatchEvolution::GroupByChamber(int timeResolution)
{
  /// From the histograms at DE level, make the corresponding ones at chamber level
  /// with left and right separated as well.

  for (Int_t ich = 1; ich <= 10; ++ich) {
    AliMpDEIterator it;

    it.First(ich - 1);

    TH1* hchamberLeft(0x0);
    TH1* hchamberRight(0x0);
    TList listLeft;
    TList listRight;
    listLeft.SetOwner(kFALSE);
    listRight.SetOwner(kFALSE);

    AliMpDCSNamer dcs("TRACKER");

    while (!it.IsDone()) {
      Int_t detElemId = it.CurrentDEId();

      AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);

      TH1* h = fBPEVO.Histo(Form("/DE/HITS/%ds/DE%04d", timeResolution, detElemId));

      if (dcs.DCSAliasName(detElemId).Contains("Left")) {
        if (!hchamberLeft) {
          hchamberLeft = static_cast<TH1*>(h->Clone(Form("CHAMBER%dLEFT", ich)));
        }
        else {
          listLeft.Add(h);
        }
      }
      else {
        if (!hchamberRight) {
          hchamberRight = static_cast<TH1*>(h->Clone(Form("CHAMBER%dRIGHT", ich)));
        }
        else {
          listRight.Add(h);
        }
      }

      it.Next();
    }

    hchamberLeft->Merge(&listLeft);
    hchamberRight->Merge(&listRight);

    hchamberLeft->SetLineColor(4);
    hchamberRight->SetLineColor(2);

    fBPEVO.Adopt(Form("/CHAMBER/HITS/%ds", timeResolution), hchamberLeft);
    fBPEVO.Adopt(Form("/CHAMBER/HITS/%ds", timeResolution), hchamberRight);
    TH1* hchamber = static_cast<TH1*>(hchamberLeft->Clone(Form("CHAMBER%d", ich)));
    hchamber->Add(hchamberRight);
    hchamber->SetLineColor(1);
    fBPEVO.Adopt(Form("/CHAMBER/HITS/%ds", timeResolution), hchamber);
  }
}

//_________________________________________________________________________________________________
void AliMUONBusPatchEvolution::GroupByDDL(int timeResolution)
{
  /// From the histograms at DE level, make the corresponding ones at DDL level

  Int_t nddls = AliDAQ::NumberOfDdls("MUONTRK");
  Int_t offset = AliDAQ::DdlIDOffset("MUONTRK");

  for (Int_t i = 0; i < nddls; ++i) {
    Int_t ddlId = offset + i;

    AliMpDDL* ddl = AliMpDDLStore::Instance()->GetDDL(i);

    TH1* hddl(0x0);
    TList list;
    list.SetOwner(kFALSE);

    for (Int_t ide = 0; ide < ddl->GetNofDEs(); ++ide) {
      Int_t detElemId = ddl->GetDEId(ide);

      TH1* h = fBPEVO.Histo(Form("/DE/HITS/%ds/DE%04d", timeResolution, detElemId));

      if (!hddl) {
        hddl = static_cast<TH1*>(h->Clone(Form("DDL%d", ddlId)));
      }
      else {
        list.Add(h);
      }
    }

    hddl->Merge(&list);
    fBPEVO.Adopt(Form("/DDL/HITS/%ds", timeResolution), hddl);
  }
}

//_________________________________________________________________________________________________
void AliMUONBusPatchEvolution::GroupByDE(int timeResolution)
{
  /// From the histograms at Bus Patch level (which is the original DA output)
  /// make the corresponding histograms at DE level

  AliMpDEIterator it;

  it.First();

  while (!it.IsDone()) {
    Int_t detElemId = it.CurrentDEId();

    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);

    TList list;
    list.SetOwner(kFALSE);
    TH1* hde(0x0);

    if (de->GetStationType() != AliMp::kStationTrigger) {
      for (Int_t i = 0; i < de->GetNofBusPatches(); ++i) {
        Int_t busPatchId = de->GetBusPatchId(i);

        TH1* h = fBPEVO.Histo(Form("/BUSPATCH/HITS/%ds/BP%04d", timeResolution, busPatchId));

        h->SetLineColor(1);

        if (!hde) {
          hde = static_cast<TH1*>(h->Clone());
          hde->SetName(Form("DE%04d", detElemId));
        }
        else {
          list.Add(h);
        }
      }

      hde->Merge(&list);
      fBPEVO.Adopt(Form("/DE/HITS/%ds", timeResolution), hde);
    }

    it.Next();
  }
}

//_________________________________________________________________________________________________
void AliMUONBusPatchEvolution::Normalize()
{
  TObjArray* a = fBPEVO.SortAllIdentifiers();
  TIter nextId(a);
  TObjString* sid;

  while ((sid = static_cast<TObjString*>(nextId()))) {

    if (!sid->String().Contains("HITS"))
      continue;

    TObjArray* parts = sid->String().Tokenize("/");
    TString npadsId("/");

    npadsId += static_cast<TObjString*>(parts->At(0))->String();
    npadsId += "/NPADS";

    TString sduration = static_cast<TObjString*>(parts->Last())->String();

    delete parts;

    Int_t duration = sduration.Atoi();

    TList* list = fBPEVO.CreateListOfObjectNames(sid->String().Data());
    TIter nextObject(list);
    TObjString* sobject;
    while ((sobject = static_cast<TObjString*>(nextObject()))) {
      TString hname;

      hname.Form("%s%s", sid->String().Data(), sobject->String().Data());

      TH1* hnhits = fBPEVO.Histo(hname.Data());

      TString padname;

      padname.Form("%s/%s", npadsId.Data(), sobject->String().Data());

      TH1* hnpads = fBPEVO.Histo(padname.Data());

      Double_t npads = hnpads->GetBinContent(1);

      TH1* hocchz = static_cast<TH1*>(hnhits->Clone());

      hocchz->SetDirectory(0);

      hocchz->SetTitle("Occupancy (Hz)");

      hocchz->Scale(1.0 / npads / duration);

      TString occ = sid->String();

      occ.ReplaceAll("HITS", "OCCHZ");

      fBPEVO.Adopt(occ.Data(), hocchz);

      TH1* hocc = static_cast<TH1*>(hnhits->Clone());

      hocc->SetDirectory(0);

      hocc->SetTitle("Occupancy");

      occ.ReplaceAll("HZ", "");

      TH1* hnevents = fBPEVO.Histo(Form("Nevents%s", sduration.Data()));

      hocc->GetXaxis()->SetTimeDisplay(0); // to avoid bin inflation in Divide...

      hocc->Divide(hnhits, hnevents, 1.0 / npads, 1.0);

      hocc->GetXaxis()->SetTimeDisplay(1);

      fBPEVO.Adopt(occ.Data(), hocc);
    }

    delete list;
  }

  delete a;
}

//_____________________________________________________________________________
void AliMUONBusPatchEvolution::ShrinkTimeAxis()
{
  /// remove unused bins (after the last used bin) in all histograms that
  /// have a time axis. Note that we _assume_ that at this
  /// point the collection indeed contains only histograms that have time bins

  TIter next(fBPEVO.CreateIterator());

  TH1F* h(0x0);

  std::map<int, int> lastbins;
  std::map<int, int> firstbins;
  std::map<int, int> nbins;
  std::map<int, std::vector<TH1F*> > histos;

  std::set<int> timeResolutions;

  // get the upper bin per time resolution (only one in order to stick to a common bin
  // definition for all histograms

  while ((h = static_cast<TH1F*>(next()))) {

    if (!h->GetXaxis()->GetTimeDisplay())
      continue;

    int timeResolution = GetTimeResolution(*h);

    if (timeResolution <= 0)
      continue;

    timeResolutions.insert(timeResolution);

    if (h->GetEntries()) {
      lastbins[timeResolution] = TMath::Max(lastbins[timeResolution], h->FindLastBinAbove(0.0));

      if (firstbins.count(timeResolution)) {
        firstbins[timeResolution] = TMath::Min(firstbins[timeResolution], h->FindFirstBinAbove(0.0));
      }
      else {
        firstbins[timeResolution] = h->FindFirstBinAbove(0.0);
      }
    }

    int nb = h->GetNbinsX();

    if (nbins.count(timeResolution)) {
      assert(nbins[timeResolution] = nb);
    }
    else {
      nbins[timeResolution] = nb;
    }

    histos[timeResolution].push_back(h);
  }

  std::set<int, int>::const_iterator it;

  TH1::AddDirectory(kFALSE);

  for (it = timeResolutions.begin(); it != timeResolutions.end(); ++it) {

    int timeRes = *it;
    int firstbin = firstbins[timeRes];
    int lastbin = lastbins[timeRes];
    int nb = nbins[timeRes];
    const std::vector<TH1F*>& v = histos[timeRes];

    std::map<int, std::vector<TH1F*> >::const_iterator hit;
    for (std::vector<TH1F*>::size_type iv = 0; iv < v.size(); ++iv) {
      TH1F* hold = static_cast<TH1F*>(v[iv]);

      TH1F* hnew = new TH1F(hold->GetName(), hold->GetTitle(), lastbin - firstbin + 1, hold->GetBinLowEdge(firstbin),
                            hold->GetBinLowEdge(lastbin) + hold->GetBinWidth(lastbin));

      Double_t nentries = hold->GetEntries();

      for (int i = firstbin; i <= lastbin; ++i) {
        hnew->SetBinContent(i - firstbin + 1, hold->GetBinContent(i));
      }

      TString timeFormat = hold->GetXaxis()->GetTimeFormat();

      hnew->Copy(*hold);

      hold->SetEntries(nentries);

      hold->GetXaxis()->SetTimeFormat(timeFormat.Data());
      hold->GetXaxis()->SetTimeDisplay(1);

      delete hnew;
    }
  }
}

