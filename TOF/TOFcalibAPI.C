/*
 * calibapi.cxx
 * set of functions to read TOFcompactCalib data and
 * deal with TOF calibration in a centralized way
 */

#include <stdio.h>
#include <signal.h>

#define MAXHITS 100000
#define MAXPOINTS 100000

/* input data */
TFile *datafile = NULL;
TTree *datatree = NULL;
Int_t nevents = 0;
Int_t curev = 0;
Int_t runNb = 0;
UInt_t timestamp = 0;
Float_t vertex = 0.;
Float_t timezero = 0.;
Int_t nhits = 0;
Float_t momentum[MAXHITS];
Float_t length[MAXHITS];
Int_t index[MAXHITS];
Float_t time[MAXHITS];
Float_t tot[MAXHITS];
Float_t texp[MAXHITS];

/* calib histos */
enum EParam_t {
    kTRM,
    kFEA,
    kChannel,
    kNParams
};
const Char_t *paramName[kNParams] = {"trm", "fea", "channel"};
Int_t paramBins[kNParams] = {720, 6552, 157248};
Float_t paramMin[kNParams] = {0., 0., 0.};
Float_t paramMax[kNParams] = {720., 6552., 157248.};
Int_t deltatBins = 201;
Float_t deltatMin = -2440. - 12.2;
Float_t deltatMax = 2440. + 12.2;
Int_t totBins = 400;
Float_t totMin = 4.88;
Float_t totMax = 24.4;
UInt_t timeZeroSampling = 600; /* seconds */
UInt_t timeMin = 0;
UInt_t timeMax = 0;
UInt_t timeBins = 0;
TH2F *hFEAHitMap = NULL;
TH2F *hChannelHitMap = NULL;
TH2F *hTimeZeroFillHisto = NULL;
TH1F *hTimeZeroFill_mean = NULL;
TH1F *hTimeZeroFill_sigma = NULL;
TProfile *hTimePressureHisto = NULL;
TH2F *hCalibHisto[kNParams] = {NULL, NULL, NULL};
TH1F *hCalibParam_mean[kNParams] = {NULL, NULL, NULL};
TH1F *hCalibParam_sigma[kNParams] = {NULL, NULL, NULL};
TH1F *hChannelParam_mean = NULL;
TH1F *hChannelParam_sigma = NULL;

TH2F *hPerfHistoDeltaT = NULL;
TH2F *hPerfHistoBeta = NULL;

AliTOFcalibHisto *calibHisto = NULL;

/* other */
#define MAXRUNS 1000000
AliLHCClockPhase *lhcClockPhase[MAXRUNS];
AliDCSSensor *cavernPressure[MAXRUNS];

/* monitoring */
TStopwatch stopwatch;

//_____________________________________________________________

Bool_t init()
{
    
    AliCDBManager *cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("raw://");
    
    /* loop over events */
    printf("init: loop over events\n");
    for (curev = 0; curev < nevents; curev++) {
        /* get and check event */
        datatree->GetEvent(curev);
        if (cdb->GetRun() == runNb) continue;
        /* set run and read GRP */
        cdb->SetRun(runNb);
        /* get LHC clock-phase */
        AliCDBEntry *cdbe = (AliCDBEntry *)cdb->Get("GRP/Calib/LHCClockPhase");
        lhcClockPhase[runNb] = (AliLHCClockPhase *)cdbe->GetObject();
    }
    
    return kFALSE;
    
}

//_____________________________________________________________

Bool_t
runCalib(const Char_t *filename, const Char_t *paramfilename = NULL, Bool_t useTimeZeroTOF = kFALSE, Bool_t useLHCClockPhase = kFALSE, Bool_t runFix = kTRUE, Int_t nSteps = 2, Int_t evMax = kMaxInt)
{
    
    /* open data */
    if (openData(filename, evMax))
        return kTRUE;
    /* open calib params */
    if (paramfilename)
        if (openCalibParams(paramfilename))
            return kTRUE;
    
    /* init if needed */
    if (useLHCClockPhase)
        init();
    
    /* fill and fit time-zero fill */
    deltatMin = -24400. - 122.;
    deltatMax = +24400. + 122.;
    fillTimeZeroFillHisto(useTimeZeroTOF, useLHCClockPhase);
    fitTimeZeroFillHisto();
    deltatMin = -2440. - 12.2;
    deltatMax = +2440. + 12.2;
    fillTimeZeroFillHisto(useTimeZeroTOF, useLHCClockPhase);
    fitTimeZeroFillHisto();
    writeCalibHistos("calibhistos.step-timezero.root");
    writeCalibParams("calibparams.step-timezero.root");
    
    /* fill with extended range and fix what is far away */
    if (runFix) {
        /* fix up to 250 ns shift */
        deltatMin = -244000. - 1220.;
        deltatMax = +244000. + 1220.;
        fillCalibHistos(useTimeZeroTOF, useLHCClockPhase);
        writeCalibHistos("calibhistos.step-fix-250ns.root");
        fitCalibHistos(10., 10., 100, kTRUE);
        writeCalibParams("calibparams.step-fix-250ns.root");
        /* fix up to 25 ns shift */
        deltatMin = -24400. - 122.;
        deltatMax = +24400. + 122.;
        fillCalibHistos(useTimeZeroTOF, useLHCClockPhase);
        writeCalibHistos("calibhistos.step-fix-25ns.root");
        fitCalibHistos(10., 10., 100, kTRUE);
        writeCalibParams("calibparams.step-fix-25ns.root");
    }
    
    /* fill and fit with limited range */
    deltatMin = -2440. - 12.2;
    deltatMax = +2440. + 12.2;
    for (Int_t istep = 0; istep < nSteps; istep++) {
        fillCalibHistos(useTimeZeroTOF, useLHCClockPhase);
        writeCalibHistos(Form("calibhistos.step-%d.root", istep));
        fitCalibHistos(3., 2., 100, kFALSE);
        writeCalibParams(Form("calibparams.step-%d.root", istep));
    }
    
    /* fill and fit common shift */
    //  fillCalibHistos(kFALSE, useLHCClockPhase);
    //  writeCalibHistos("calibhistos.step-common.root");
    //  fitCommonShift();
    //  writeCalibParams("calibparams.step-common.root");
    
    /* final resuls */
    fillCalibHistos(useTimeZeroTOF, useLHCClockPhase);
    writeCalibHistos("calibhistos.root");
    writeCalibParams("calibparams.root");

    return kFALSE;
}

//_____________________________________________________________

Bool_t
checkCalib(const Char_t *filename, const Char_t *paramfilename = NULL, Bool_t useTimeZeroTOF = kFALSE, Bool_t useLHCClockPhase = kFALSE, Bool_t runExtended = kTRUE, Int_t evMax = kMaxInt)
{
    
    /* open data */
    if (openData(filename, evMax))
        return kTRUE;
    /* open calib params */
    if (paramfilename)
        if (openCalibParams(paramfilename))
            return kTRUE;
    
    /* fill and fit time-zero fill */
    deltatMin = -24400.;
    deltatMax = +24400.;
    fillTimeZeroFillHisto(useTimeZeroTOF, useLHCClockPhase);
    fitTimeZeroFillHisto();
    deltatMin = -2440.;
    deltatMax = +2440.;
    fillTimeZeroFillHisto(useTimeZeroTOF, useLHCClockPhase);
    fitTimeZeroFillHisto();
    
    /* fill with extended range */
    if (runExtended) {
        deltatMin = -244000.;
        deltatMax = +244000.;
        fillCalibHistos(useTimeZeroTOF, useLHCClockPhase);
        writeCalibHistos("checkcalib.extended-250ns.root");
        deltatMin = -24400.;
        deltatMax = +24400.;
        fillCalibHistos(useTimeZeroTOF, useLHCClockPhase);
        writeCalibHistos("checkcalib.extended-25ns.root");
    }
    
    /* fill with limited range */
    deltatMin = -2440.;
    deltatMax = +2440.;
    fillCalibHistos(useTimeZeroTOF, useLHCClockPhase);
    writeCalibHistos("checkcalib.root");
    writeCalibParams("checkcalib.calibparams.root");
    
    return kFALSE;
}

//_____________________________________________________________

Bool_t
checkPerf(const Char_t *filename, const Char_t *paramfilename = NULL, Bool_t useTimeZeroTOF = kFALSE, Bool_t useLHCClockPhase = kFALSE, Int_t evMax = kMaxInt)
{
    
    /* open data */
    if (openData(filename, evMax))
        return kTRUE;
    /* open calib params */
    if (paramfilename)
        if (openCalibParams(paramfilename))
            return kTRUE;
    
    fillTimeZeroFillHisto(useTimeZeroTOF, useLHCClockPhase);
    fitTimeZeroFillHisto();
    fillPerfHistos(useTimeZeroTOF, useLHCClockPhase);
    writePerfHistos("checkperf.root");
    
    return kFALSE;
}

//_____________________________________________________________

Bool_t
openData(const Char_t *filename, Int_t evMax = kMaxInt)
{
    /*
     * open data
     */
    
    /* open file */
    printf("openData: opening file: %s\n", filename);
    datafile = TFile::Open(filename);
    if (!datafile || !datafile->IsOpen()) {
        printf("openData: cannot open file %s\n", filename)
        return kTRUE;
    }
    /* get tree */
    datatree = (TTree *)datafile->Get("aodTree");
    if (!datatree) {
        printf("openData: cannot find \'aodTree\' tree in %s\n", filename);
        return kTRUE;
    }
    nevents = datatree->GetEntries();
    printf("openData: got \'aodTree\' tree: %d entries\n", nevents);
    if (nevents > evMax) {
        printf("openData: setting max event to %d as requested\n", evMax);
        nevents = evMax;
    }
    /* connect inputs */
    datatree->SetBranchAddress("run", &runNb);
    datatree->SetBranchAddress("timestamp", &timestamp);
    datatree->SetBranchAddress("vertex", &vertex);
    datatree->SetBranchAddress("timezero", &timezero);
    datatree->SetBranchAddress("nhits", &nhits);
    datatree->SetBranchAddress("momentum", &momentum);
    datatree->SetBranchAddress("length", &length);
    datatree->SetBranchAddress("index", &index);
    datatree->SetBranchAddress("time", &time);
    datatree->SetBranchAddress("tot", &tot);
    datatree->SetBranchAddress("texp", &texp);
    
    /* get first event, and retrieve the year */
    datatree->GetEvent(0);
    TTimeStamp ts(timestamp);
    UInt_t year;
    ts.GetDate(kTRUE, 0, &year);
    TTimeStamp firstTimestamp(year, 1, 1, 0, 0, 0, 0, kTRUE);
    TTimeStamp lastTimestamp(year + 1, 1, 1, 0, 0, 0, 0, kTRUE);
    timeMin = firstTimestamp.GetTimeSpec().tv_sec;
    timeMax = lastTimestamp.GetTimeSpec().tv_sec;
    timeBins = (timeMax - timeMin) / timeZeroSampling;
    printf("openData: calibration running on %d data\n", year);
    
    return kFALSE;
}

//_____________________________________________________________

Bool_t
acceptEvent(Bool_t useTimeZeroTOF = kFALSE)
{
    /*
     * acceptEvent
     */
    
    if (useTimeZeroTOF && timezero == 999999.) return kFALSE;
    return kTRUE;
}

//_____________________________________________________________

Float_t
getTimeZeroFill(UInt_t ts)
{
    /*
     * getTimeZeroFill
     */
    
    if (!hTimeZeroFill_mean) return 0.;
    Int_t tsbin = hTimeZeroFill_mean->FindBin(ts);
    return hTimeZeroFill_mean->GetBinContent(tsbin);
}

//_____________________________________________________________

Float_t
getLHCClockPhase(UInt_t ts)
{
    /*
     * getLHCClockPhase
     */
    
    return lhcClockPhase[runNb] ? 1.e3 * lhcClockPhase[runNb]->GetPhase(ts) : 0.;
}

//_____________________________________________________________

Float_t
getCalib(Int_t idx)
{
    /*
     * getCalib
     */
    
    if (!hChannelParam_mean) return 0.;
    return hChannelParam_mean->GetBinContent(idx + 1);
}

//_____________________________________________________________

Float_t
getDeltaT(Int_t ihit, Bool_t useTimeZeroTOF = kFALSE, Bool_t useLHCClockPhase = kFALSE)
{
    /*
     * getDeltaT
     */
    
    Float_t val = time[ihit] - getCalib(index[ihit]) - getTimeZeroFill(timestamp) - texp[ihit];
    if (useTimeZeroTOF) val -= timezero ;
    if (useLHCClockPhase) val += getLHCClockPhase(timestamp);
    return val;
}

//_____________________________________________________________

Float_t
getBeta(Int_t ihit, Bool_t useTimeZeroTOF = kFALSE, Bool_t useLHCClockPhase = kFALSE)
{
    /*
     * getBeta
     */
    
    Float_t tof = time[ihit] - getCalib(index[ihit]) - getTimeZeroFill(timestamp);
    if (useTimeZeroTOF) tof -= timezero ;
    if (useLHCClockPhase) tof += getLHCClockPhase(timestamp);
    Float_t beta = length[ihit] / tof / 2.99792458000000000e-2;
    return beta;
}

//_____________________________________________________________

Float_t
getMass(Int_t ihit, Bool_t useTimeZeroTOF = kFALSE, Bool_t useLHCClockPhase = kFALSE)
{
    /*
     * getMass
     */
    
    Float_t tof = time[ihit] - getCalib(index[ihit]) - getTimeZeroFill(timestamp);
    if (useTimeZeroTOF) tof -= timezero ;
    if (useLHCClockPhase) tof += getLHCClockPhase(timestamp);
    Float_t beta = length[ihit] / tof / 2.99792458000000000e-2;
    if (beta > 1.) beta *= -1.
        Float_t mass = momentum[ihit] * TMath::Sqrt(1. / (beta * beta) - 1.)
        return mass;
}

//_____________________________________________________________

Int_t
getIndex(Int_t param, Int_t idx)
{
    /*
     * getIndex
     */
    
    if (!calibHisto) {
        calibHisto = new AliTOFcalibHisto();
        calibHisto->LoadCalibHisto();
    }
    
    /* switch param */
    Int_t sector, ddl, trm, strip, padx, paramIndex;
    switch (param) {
        case kTRM:
            ddl = calibHisto->GetCalibMap(AliTOFcalibHisto::kDDL, idx);
            trm = calibHisto->GetCalibMap(AliTOFcalibHisto::kTRM, idx);
            paramIndex = trm + 10 * ddl;
            break;
        case kFEA:
            sector = calibHisto->GetCalibMap(AliTOFcalibHisto::kSector, idx);
            strip = calibHisto->GetCalibMap(AliTOFcalibHisto::kSectorStrip, idx);
            padx = calibHisto->GetCalibMap(AliTOFcalibHisto::kPadX, idx);
            paramIndex = padx / 12 + 4 * strip + 364 * sector;
            break;
        case kChannel:
            paramIndex = idx;
            break;
    }
    
    return paramIndex;
}

//_____________________________________________________________

Int_t
getHitMapXY(Int_t param, Int_t idx, Float_t *hitmap)
{
    /*
     * getHitMapXY
     */
    
    if (!calibHisto) {
        calibHisto = new AliTOFcalibHisto();
        calibHisto->LoadCalibHisto();
    }
    
    /* get from calib histo */
    Int_t sector, strip, padx, padz, fea;
    sector = calibHisto->GetCalibMap(AliTOFcalibHisto::kSector, idx);
    strip = calibHisto->GetCalibMap(AliTOFcalibHisto::kSectorStrip, idx);
    padx = calibHisto->GetCalibMap(AliTOFcalibHisto::kPadX, idx);
    padz = calibHisto->GetCalibMap(AliTOFcalibHisto::kPadZ, idx);
    fea = padx / 12;
    /* switch param */
    switch (param) {
        case kFEA:
            hitmap[0] = sector + ((Double_t)(3 - fea) + 0.5) / 4.;
            hitmap[1] = strip;
            break;
        case kChannel:
            hitmap[0] = sector + ((Double_t)(47 - padx) + 0.5) / 48.;
            hitmap[1] = strip + ((Double_t)(padz) + 0.5) / 2.;
            break;
        default:
            hitmap[0] = 0.;
            hitmap[1] = 0.;
            break;
    }
    
}

//_____________________________________________________________

void
fillTimeZeroFillHisto(Bool_t useTimeZeroTOF = kFALSE, Bool_t useLHCClockPhase = kFALSE)
{
    /*
     * fillTimeZeroFillHisto
     */
    
    /* create/reset histos */
    if (hTimeZeroFillHisto) delete hTimeZeroFillHisto;
    hTimeZeroFillHisto = new TH2F("hTimeZeroFillHisto", "", timeBins, timeMin, timeMax, deltatBins, deltatMin, deltatMax);
    if (hTimePressureHisto) delete hTimePressureHisto;
    hTimePressureHisto = new TProfile("hTimePressureHisto", "", timeBins, timeMin, timeMax);
    
    /* reset and start stopwatch */
    stopwatch.Reset();
    stopwatch.Start();
    
    /* loop over events */
    printf("fillTimeZeroFillHisto: loop over events\n");
    if (useTimeZeroTOF)
        printf("fillTimeZeroFillHisto: time-zero TOF requested\n");
    else
        printf("fillTimeZeroFillHisto: not using time-zero TOF\n");
    if (useLHCClockPhase)
        printf("fillTimeZeroFillHisto: BPTX clock-phase requested\n");
    else
        printf("fillTimeZeroFillHisto: not using BPTX clock-phase\n");
    Float_t hitmap[2], deltat;
    for (curev = 0; curev < nevents; curev++) {
        /* get and check event */
        datatree->GetEvent(curev);
        if (!acceptEvent(useTimeZeroTOF)) continue;
        /* loop over hits */
        for (Int_t ihit = 0; ihit < nhits; ihit++) {
            deltat = getDeltaT(ihit, useTimeZeroTOF, useLHCClockPhase);
            /* fill histos */
            hTimeZeroFillHisto->Fill(timestamp, deltat);
        } /* end of loop over hits */
    } /* end of loop over events */
    
    /* print monitor */
    monitor();
    
}

//_____________________________________________________________

void
fillCalibHistos(Bool_t useTimeZeroTOF = kFALSE, Bool_t useLHCClockPhase = kFALSE)
{
    /*
     * fillCalibHistos
     */
    
    /* create/reset histos */
    if (!hFEAHitMap)
        hFEAHitMap = new TH2F("hFEAHitMap", "", 72, 0., 18., 91, 0., 91.);
    else
        hFEAHitMap->Reset();
    
    if (!hChannelHitMap)
        hChannelHitMap = new TH2F("hChannelHitMap", "", 864, 0., 18., 91, 0., 91.);
    else
        hChannelHitMap->Reset();
    
    for (Int_t iparam = 0; iparam < kNParams; iparam++) {
        if (hCalibHisto[iparam]) delete hCalibHisto[iparam];
        hCalibHisto[iparam] = new TH2F(Form("hCalibHisto_%s", paramName[iparam]), "", paramBins[iparam], paramMin[iparam], paramMax[iparam], deltatBins, deltatMin, deltatMax);
    }
    
    /* reset and start stopwatch */
    stopwatch.Reset();
    stopwatch.Start();
    
    /* loop over events */
    printf("fillCalibHistos: loop over events\n");
    if (useTimeZeroTOF)
        printf("fillCalibHistos: time-zero TOF requested\n");
    else
        printf("fillCalibHistos: not using time-zero TOF\n");
    if (useLHCClockPhase)
        printf("fillCalibHistos: BPTX clock-phase requested\n");
    else
        printf("fillCalibHistos: not using BPTX clock-phase\n");
    Float_t hitmap[2], deltat;
    for (curev = 0; curev < nevents; curev++) {
        /* get and check event */
        datatree->GetEvent(curev);
        if (!acceptEvent(useTimeZeroTOF)) continue;
        /* loop over hits */
        for (Int_t ihit = 0; ihit < nhits; ihit++) {
            deltat = getDeltaT(ihit, useTimeZeroTOF, useLHCClockPhase);
            /* fill histos */
            getHitMapXY(kFEA, index[ihit], hitmap);
            hFEAHitMap->Fill(hitmap[0], hitmap[1]);
            getHitMapXY(kChannel, index[ihit], hitmap);
            hChannelHitMap->Fill(hitmap[0], hitmap[1]);
            for (Int_t iparam = 0; iparam < kNParams; iparam++) {
                hCalibHisto[iparam]->Fill(getIndex(iparam, index[ihit]), deltat);
            }
        } /* end of loop over hits */
    } /* end of loop over events */
    
    /* print monitor */
    monitor();
    
}

//_____________________________________________________________

void
fillPerfHistos(Bool_t useTimeZeroTOF = kFALSE, Bool_t useLHCClockPhase = kFALSE)
{
    /*
     * fillPerfHistos
     */
    
    /* create/reset histos */
    if (!hFEAHitMap)
        hFEAHitMap = new TH2F("hFEAHitMap", "", 72, 0., 18., 91, 0., 91.);
    else
        hFEAHitMap->Reset();
    
    if (!hChannelHitMap)
        hChannelHitMap = new TH2F("hChannelHitMap", "", 864, 0., 18., 91, 0., 91.);
    else
        hChannelHitMap->Reset();
    
    if (!hPerfHistoDeltaT)
        hPerfHistoDeltaT = new TH2F("hPerfHistoDeltaT", "", 250, 0., 5., 20000, -244000., 244000.);
    else
        hPerfHistoDeltaT->Reset();
    
    if (!hPerfHistoBeta)
        hPerfHistoBeta = new TH2F("hPerfHistoBeta", "", 1000, 0., 10., 2200, 0., 1.1);
    else
        hPerfHistoBeta->Reset();
    
    
    /* reset and start stopwatch */
    stopwatch.Reset();
    stopwatch.Start();
    
    /* loop over events */
    printf("fillPerfHistos: loop over events\n");
    if (useTimeZeroTOF)
        printf("fillPerfHistos: time-zero TOF requested\n");
    else
        printf("fillPerfHistos: not using time-zero TOF\n");
    if (useLHCClockPhase)
        printf("fillPerfHisto: BPTX clock-phase requested\n");
    else
        printf("fillPerfHisto: not using BPTX clock-phase\n");
    Float_t hitmap[2], deltat, beta, mass;
    for (curev = 0; curev < nevents; curev++) {
        /* get and check event */
        datatree->GetEvent(curev);
        if (!acceptEvent(useTimeZeroTOF)) continue;
        /* loop over hits */
        for (Int_t ihit = 0; ihit < nhits; ihit++) {
            deltat = getDeltaT(ihit, useTimeZeroTOF, useLHCClockPhase);
            beta = getBeta(ihit, useTimeZeroTOF, useLHCClockPhase);
            
            //      mass = getMass(ihit, useTimeZeroTOF, useLHCClockPhase);
            /* fill histos */
            getHitMapXY(kFEA, index[ihit], hitmap);
            hFEAHitMap->Fill(hitmap[0], hitmap[1]);
            getHitMapXY(kChannel, index[ihit], hitmap);
            hChannelHitMap->Fill(hitmap[0], hitmap[1]);
            hPerfHistoDeltaT->Fill(momentum[ihit], deltat);
            hPerfHistoBeta->Fill(momentum[ihit], beta);
            //      hPerfHistoMass->Fill(p, mass);
        } /* end of loop over hits */
    } /* end of loop over events */
    
    /* print monitor */
    monitor();
    
}

//_____________________________________________________________

Bool_t
writeCalibHistos(const Char_t *filename)
{
    /*
     * writeCalibHistos
     */
    
    /* open output file and write */
    TFile *fileout = TFile::Open(filename, "RECREATE");
    if (!fileout || !fileout->IsOpen()) {
        printf("writeCalibHistos: error while opening output file %s\n", filename);
        return kTRUE;
    }
    if (hFEAHitMap) hFEAHitMap->Write();
    if (hChannelHitMap) hChannelHitMap->Write();
    if (hTimeZeroFillHisto) hTimeZeroFillHisto->Write();
    if (hTimePressureHisto) hTimePressureHisto->Write();
    for (Int_t iparam = 0; iparam < kNParams; iparam++) {
        if (hCalibHisto[iparam]) hCalibHisto[iparam]->Write();
    }
    fileout->Close();
    printf("writeCalibHistos: output written on %s\n", filename);
    return kFALSE;
}

//_____________________________________________________________

Bool_t
writePerfHistos(const Char_t *filename)
{
    /*
     * writePerfHistos
     */
    
    /* open output file and write */
    TFile *fileout = TFile::Open(filename, "RECREATE");
    if (!fileout || !fileout->IsOpen()) {
        printf("writePerfHistos: error while opening output file %s\n", filename);
        return kTRUE;
    }
    if (hFEAHitMap) hFEAHitMap->Write();
    if (hChannelHitMap) hChannelHitMap->Write();
    if (hPerfHistoDeltaT) hPerfHistoDeltaT->Write();
    if (hPerfHistoBeta) hPerfHistoBeta->Write();
    fileout->Close();
    printf("writePerfHistos: output written on %s\n", filename);
    return kFALSE;
}

//_____________________________________________________________

Bool_t
writeCalibParams(const Char_t *filename)
{
    /*
     * writeCalibParams
     */
    
    /* open output file and write */
    TFile *fileout = TFile::Open(filename, "RECREATE");
    if (!fileout || !fileout->IsOpen()) {
        printf("writeCalibParams: error while opening output file %s\n", filename);
        return kTRUE;
    }
    if (hTimeZeroFill_mean) hTimeZeroFill_mean->Write();
    if (hTimeZeroFill_sigma) hTimeZeroFill_sigma->Write();
    for (Int_t iparam = 0; iparam < kNParams; iparam++) {
        if (hCalibParam_mean[iparam]) hCalibParam_mean[iparam]->Write();
        if (hCalibParam_sigma[iparam]) hCalibParam_sigma[iparam]->Write();
    }
    if (hChannelParam_mean) hChannelParam_mean->Write();
    if (hChannelParam_sigma) hChannelParam_sigma->Write();
    fileout->Close();
    printf("writeCalibParams: output written on %s\n", filename);
    return kFALSE;
}

//_____________________________________________________________

Bool_t
openCalibHistos(const Char_t *filename)
{
    /*
     * openCalibHistos
     */
    
    /* open  file and write */
    printf("openCalibHistos: opening file: %s\n", filename);
    TFile *filein = TFile::Open(filename);
    if (!filein || !filein->IsOpen()) {
        printf("openCalibHistos: cannot open file %s\n", filename);
        return kTRUE;
    }
    /* time-zero fill */
    if (hTimeZeroFillHisto) {
        printf("openCalibHistos: histo \'hTimeZeroFillHisto\' will be replaced\n");
        delete hTimeZeroFillHisto;
    }
    hTimeZeroFillHisto = (TH2F *)filein->Get("hTimeZeroFillHisto");
    if (!hTimeZeroFillHisto) {
        printf("openCalibHistos: cannot get \'hTimeZeroFillHisto\' from file %s\n", filename);
        return kTRUE;
    }
    printf("openCalibHistos: got \'hTimeZeroFillHisto\' from file\n");
    for (Int_t iparam = 0; iparam < kNParams; iparam++) {
        /* calib */
        if (hCalibHisto[iparam]) {
            printf("openCalibHistos: histo \'hCalibHisto_%s\' will be replaced\n", paramName[iparam]);
            delete hCalibHisto[iparam];
        }
        hCalibHisto[iparam] = (TH2F *)filein->Get(Form("hCalibHisto_%s", paramName[iparam]));
        if (!hCalibHisto[iparam]) {
            printf("openCalibHistos: cannot get \'hCalibHisto_%s\' from file %s\n", paramName[iparam], filename);
            return kTRUE;
        }
        printf("openCalibHistos: got \'hCalibHisto_%s\' from file\n", paramName[iparam]);
    }
    
    return kFALSE;
}

//_____________________________________________________________

Bool_t
openCalibParams(const Char_t *filename)
{
    /*
     * openCalibParams
     */
    
    /* open file */
    printf("openCalibParams: opening file: %s\n", filename);
    TFile *filein = TFile::Open(filename);
    if (!filein || !filein->IsOpen()) {
        printf("openCalibParams: cannot open file %s\n", filename);
        return kTRUE;
    }
    /* mean */
    if (hChannelParam_mean) {
        printf("openCalibParams: histo \'hChannelParam_mean\' will be replaced\n");
        delete hChannelParam_mean;
    }
    hChannelParam_mean = (TH1F *)filein->Get("hChannelParam_mean");
    if (!hChannelParam_mean) {
        printf("openCalibParams: cannot get \'hChannelParam_mean\' from file %s\n", filename);
        return kTRUE;
    }
    printf("openCalibParams: got \'hChannelParam_mean\' from file\n");
    /* sigma */
    if (hChannelParam_sigma) {
        printf("openCalibParams: histo \'hChannelParam_sigma\' will be replaced\n");
        delete hChannelParam_sigma;
    }
    hChannelParam_sigma = (TH1F *)filein->Get("hChannelParam_sigma");
    if (!hChannelParam_sigma) {
        printf("openCalibParams: cannot get \'hChannelParam_sigma\' from file %s\n", filename);
        return kTRUE;
    }
    printf("openCalibParams: got \'hChannelParam_sigma\' from file\n");
    
    return kFALSE;
}

//_____________________________________________________________

Bool_t
fitTimeZeroFillHisto(Float_t nSigmaMin = 2., Float_t nSigmaMax = 1., Int_t minIntegral = 1000, Bool_t useMaxBin = kFALSE)
{
    /*
     * fitTimeZeroFillHisto
     */
    
    //  if (!hChannelParam_mean || !hCalibHisto[kChannel]) return kTRUE;
    
    if (!hTimeZeroFill_mean) hTimeZeroFill_mean = new TH1F("hTimeZeroFill_mean", "", timeBins, timeMin, timeMax);
    if (!hTimeZeroFill_sigma) hTimeZeroFill_sigma = new TH1F("hTimeZeroFill_sigma", "", timeBins, timeMin, timeMax);
    
    if (useMaxBin)
        printf("fitTimeZeroFill: fitting time-zero fill (using max bin, intMin=%d)\n", minIntegral);
    else
        printf("fitTimeZeroFill: fitting time-zero fill (sigmaMin=%.1f, sigmaMax=%.1f, intMin=%d)\n", nSigmaMin, nSigmaMax, minIntegral);
    
    TF1 *fitFunc = (TF1 *)gROOT->GetFunction("gaus");
    for (Int_t ibin = 0; ibin < hTimeZeroFillHisto->GetNbinsX(); ibin++) {
        TH1D *hpy = hTimeZeroFillHisto->ProjectionY("hpy", ibin + 1, ibin +1);
        if (hpy->GetEntries() <= minIntegral) {
            delete hpy;
            continue;
        }
        if (fitPeak(fitFunc, hpy, nSigmaMin, nSigmaMax) != 0) {
            delete hpy;
            continue;
        }
        hTimeZeroFill_mean->AddBinContent(ibin + 1, fitFunc->GetParameter(1));
        hTimeZeroFill_mean->SetBinError(ibin + 1, fitFunc->GetParError(1));
        hTimeZeroFill_sigma->SetBinContent(ibin + 1, fitFunc->GetParameter(2));
        hTimeZeroFill_sigma->SetBinError(ibin + 1, fitFunc->GetParError(2));
        delete hpy;
    }
    
    return kFALSE;
    
}

//_____________________________________________________________

Bool_t
fitCommonShift(Float_t nSigmaMin = 2., Float_t nSigmaMax = 1., Int_t minIntegral = 1000, Bool_t useMaxBin = kFALSE)
{
    /*
     * fitCommonShift
     */
    
    if (!hChannelParam_mean || !hCalibHisto[kChannel]) return kTRUE;
    
    if (useMaxBin)
        printf("fitCommonShift: fitting common shift (using max bin, intMin=%d)\n", minIntegral);
    else
        printf("fitCommonShift: fitting common shift (sigmaMin=%.1f, sigmaMax=%.1f, intMin=%d)\n", nSigmaMin, nSigmaMax, minIntegral);
    
    TF1 *fitFunc = (TF1 *)gROOT->GetFunction("gaus");
    TH1D *hpy = hCalibHisto[kChannel]->ProjectionY("hpy");
    if (hpy->GetEntries() <= minIntegral) {
        delete hpy;
        continue;
    }
    if (useMaxBin) {
        for (Int_t idx = 0; idx < 157248; idx++)
            hChannelParam_mean->AddBinContent(idx + 1, hpy->GetBinCenter(hpy->GetMaximumBin()));
        delete hpy;
        return kFALSE;
    }
    if (fitPeak(fitFunc, hpy, nSigmaMin, nSigmaMax) != 0) {
        delete hpy;
        return kTRUE;
    }
    for (Int_t idx = 0; idx < 157248; idx++)
        hChannelParam_mean->AddBinContent(idx + 1, fitFunc->GetParameter(1));
    delete hpy;
    
    return kFALSE;
    
}

//_____________________________________________________________

Bool_t
fitCalibHistos(Float_t nSigmaMin, Float_t nSigmaMax, Int_t minIntegral, Bool_t useMaxBin)
{
    /*
     * fitCalibHistos
     */
    
    for (Int_t iparam = 0; iparam < kNParams; iparam++)
        fitCalibHisto(iparam, nSigmaMin, nSigmaMax, minIntegral, useMaxBin);
    updateChannelParams();
    
}

//_____________________________________________________________

Bool_t
updateChannelParams()
{
    
    if (!hChannelParam_mean) hChannelParam_mean = new TH1F("hChannelParam_mean", "", paramBins[kChannel], paramMin[kChannel], paramMax[kChannel]);
    if (!hChannelParam_sigma) hChannelParam_sigma = new TH1F("hChannelParam_sigma", "", paramBins[kChannel], paramMin[kChannel], paramMax[kChannel]);
    
    Double_t mean, meanerr, sigma, sigmaerr;
    for (Int_t idx = 0; idx < 157248; idx++) {
        for (Int_t iparam = kNParams - 1; iparam >= 0; iparam--) {
            if (!hCalibParam_mean[iparam]) continue;
            if (hCalibParam_mean[iparam]->GetBinError(getIndex(iparam, idx) + 1) == 0.) continue;
            mean = hCalibParam_mean[iparam]->GetBinContent(getIndex(iparam, idx) + 1);
            meanerr = hCalibParam_mean[iparam]->GetBinError(getIndex(iparam, idx) + 1);
            sigma = hCalibParam_sigma[iparam]->GetBinContent(getIndex(iparam, idx) + 1);
            sigmaerr = hCalibParam_sigma[iparam]->GetBinError(getIndex(iparam, idx) + 1);
            hChannelParam_mean->AddBinContent(idx + 1, mean);
            hChannelParam_mean->SetBinError(idx + 1, meanerr);
            hChannelParam_sigma->SetBinContent(idx + 1, sigma);
            hChannelParam_sigma->SetBinError(idx + 1, sigmaerr);
            break;
        }
    }
    
}

//_____________________________________________________________

Bool_t
fitCalibHisto(Int_t param, Float_t nSigmaMin, Float_t nSigmaMax, Int_t minIntegral, Bool_t useMaxBin)
{
    /*
     * fitCalibHisto
     */
    
    /* check data histo */
    if (!hCalibHisto[param]) {
        printf("fitCalibHisto: cannot get \'hCalibHisto_%s\'\n", paramName[param]);
        return kTRUE;
    }
    /* check mean histo */
    if (!hCalibParam_mean[param]) hCalibParam_mean[param] = new TH1F(Form("hCalibParam_%s_mean", paramName[param]), "", paramBins[param], paramMin[param], paramMax[param]);
    else hCalibParam_mean[param]->Reset();
    /* check sigma histo */
    if (!hCalibParam_sigma[param]) hCalibParam_sigma[param] = new TH1F(Form("hCalibParam_%s_sigma", paramName[param]), "", paramBins[param], paramMin[param], paramMax[param]);
    else hCalibParam_sigma[param]->Reset();
    /* fit */
    return fitCalibHisto(hCalibHisto[param], hCalibParam_mean[param], hCalibParam_sigma[param], nSigmaMin, nSigmaMax, minIntegral, useMaxBin);
}

//_____________________________________________________________

Bool_t
fitCalibHisto(TH2F *hdata, TH1F *hmean, TH1F *hsigma, Float_t nSigmaMin, Float_t nSigmaMax, Int_t minIntegral, Bool_t useMaxBin)
{
    /*
     * fitCalibHisto
     */
    
    if (!hdata || !hmean || !hsigma) return kTRUE;
    if (useMaxBin)
        printf("fitCalibHisto: fitting \'%s\' (using max bin, intMin=%d)\n", hdata->GetName(), minIntegral);
    else
        printf("fitCalibHisto: fitting \'%s\' (sigmaMin=%.1f, sigmaMax=%.1f, intMin=%d)\n", hdata->GetName(), nSigmaMin, nSigmaMax, minIntegral);
    TF1 *fitFunc = (TF1 *)gROOT->GetFunction("gaus");
    Int_t nDone = 0;
    for (Int_t ibin = 0; ibin < hdata->GetNbinsX(); ibin++) {
        TH1D *hpy = hdata->ProjectionY("hpy", ibin + 1, ibin + 1);
        if (hpy->GetEntries() <= minIntegral) {
            delete hpy;
            continue;
        }
        if (useMaxBin) {
            hmean->SetBinContent(ibin + 1, hpy->GetBinCenter(hpy->GetMaximumBin()));
            hmean->SetBinError(ibin + 1, hpy->GetBinWidth(hpy->GetMaximumBin()));
            hsigma->SetBinContent(ibin + 1, hpy->GetBinWidth(hpy->GetMaximumBin()));
            hsigma->SetBinError(ibin + 1, hpy->GetBinWidth(hpy->GetMaximumBin()));
            //      hmean->SetBinContent(ibin + 1, hpy->GetMean());
            //      hmean->SetBinError(ibin + 1, hpy->GetMeanError());
            //      hsigma->SetBinContent(ibin + 1, hpy->GetRMS());
            //      hsigma->SetBinError(ibin + 1, hpy->GetRMSError());
            nDone++;
            delete hpy;
            continue;
        }
        /* fit and check result */
        if (fitPeak(fitFunc, hpy, nSigmaMin, nSigmaMax) != 0) {
            delete hpy;
            continue;
        }
#if 0
        /* check mean value larger than error */
        if (TMath::Abs(fitFunc->GetParameter(1)) < fitFunc->GetParameter(2)) {
            delete hpy;
            continue;
        }
#endif
        hmean->SetBinContent(ibin + 1, fitFunc->GetParameter(1));
        hmean->SetBinError(ibin + 1, fitFunc->GetParError(1));
        hsigma->SetBinContent(ibin + 1, fitFunc->GetParameter(2));
        hsigma->SetBinError(ibin + 1, fitFunc->GetParError(2));
        nDone++;
        delete hpy;
    }
    printf("fitCalibHisto: %d done\n", nDone);
    return kFALSE;
}

//_____________________________________________________________

problematicChannels(const Char_t *filename, Float_t nsigma = 5., Int_t minHits = 0, Int_t startRun = 0, Int_t endRun = AliCDBRunRange::Infinity(), const Char_t *dbString = "local://$HOME/OCDB")
{
    
    TFile *fin = TFile::Open(filename);
    TH2 *hin = (TH2 *)fin->Get("hCalibHisto_channel");
    TH2 *hingood = (TH2 *)hin->Clone("hingood");
    TH2 *hinbad = (TH2 *)hin->Clone("hinbad");
    
    TH1D *hpy;
    
    TH1F *hEntries = new TH1F("hEntries", "", 10. * hin->GetEntries() / 157248, 0., 5. * hin->GetEntries() / 157248);
    for (Int_t ich = 0; ich < 157248; ich++) {
        hpy = hin->ProjectionY("hpy", ich + 1, ich + 1);
        hEntries->Fill(hpy->GetEntries());
    }
    
    new TCanvas("cEntries");
    hEntries->Draw();
    
    TH1D *hinpy_all = hin->ProjectionY("hinpy_all");
    TF1 *fitFunc = (TF1 *)gROOT->GetFunction("gaus");
    fitPeak(fitFunc, hinpy_all, 2., 1.);
    Double_t mean = fitFunc->GetParameter(1);
    Double_t sigma = fitFunc->GetParameter(2);
    Double_t intmin = mean - nsigma * sigma;
    Double_t intmax = mean + nsigma * sigma;
    Int_t binmin = hinpy_all->FindBin(intmin);
    Int_t binmax = hinpy_all->FindBin(intmax);
    
    Double_t integral_all = hinpy_all->Integral(binmin, binmax);
    Double_t entries_all = hinpy_all->GetEntries();
    Double_t fraction_all = integral_all / entries_all;
    
    printf("found %f hits in +- %d sigma\n", fraction_all, nsigma);
    
    TH1D *hFraction = new TH1D("hFraction", "", 2000, 0., 2.);
    
    TH1C *obj = new TH1C("hProblematic", "", 157248, 0., 157248.);
    
    TH2F *hFEAHitMap = new TH2F("hFEAHitMap", "", 72, 0., 18., 91, 0., 91.);
    
    Float_t hitmap[2];
    
    for (Int_t ich = 0; ich < 157248; ich++) {
        hpy = hin->ProjectionY("hpy", ich + 1, ich + 1);
        if (hpy->GetEntries() <= minHits) {
            delete hpy;
            continue;
        }
        Double_t integral_ch = hpy->Integral(binmin, binmax);
        Double_t entries_ch = hpy->GetEntries();
        Double_t fraction_ch = integral_ch / entries_ch / fraction_all;
        hFraction->Fill(fraction_ch);
        
    }
    
    TF1 *dgaus = new TF1("dgaus", "[0] * ((x < [1]) * TMath::Gaus(x, [1], [2]) + (x > [1]) * TMath::Gaus(x, [1], [3]))", 0., 2.);
    dgaus->SetParameter(0, 1.);
    dgaus->SetParameter(1, 1.);
    dgaus->SetParameter(2, 1.);
    dgaus->SetParameter(3, 1.);
    dgaus->SetParLimits(2, 0., 1.);
    dgaus->SetParLimits(3, 0., 1.);
  
    hFraction->Fit(dgaus, "0", "");
    hFraction->Fit(dgaus, "0", "", dgaus->GetParameter(1) - 3. * dgaus->GetParameter(2), dgaus->GetParameter(1) + 3. * dgaus->GetParameter(3));
    Double_t cutLimit = dgaus->GetParameter(1) - 3. * dgaus->GetParameter(2);
    Double_t cutLimit2 = dgaus->GetParameter(1) - 3. * dgaus->GetParameter(3);
    
    TLine *lcutfrac = new TLine(cutLimit, 0, cutLimit, hFraction->GetMaximum());
    TLine *lcutfrac2 = new TLine(cutLimit2, 0, cutLimit2, hFraction->GetMaximum());
    
    Int_t ntotch = 0, nprobch = 0;
    
    for (Int_t ich = 0; ich < 157248; ich++) {
        hpy = hin->ProjectionY("hpy", ich + 1, ich + 1);
        if (hpy->GetEntries() <= minHits) {
            delete hpy;
            continue;
        }
        Double_t integral_ch = hpy->Integral(binmin, binmax);
        Double_t entries_ch = hpy->GetEntries();
        Double_t fraction_ch = integral_ch / entries_ch / fraction_all;
        
        ntotch++;
        if (fraction_ch < cutLimit) {
            nprobch++;
            obj->SetBinContent(ich + 1, 0x1);
            getHitMapXY(kFEA, ich, hitmap);
            hFEAHitMap->Fill(hitmap[0], hitmap[1]);
            for (Int_t i = 0; i < hingood->GetNbinsY(); i++)
                hingood->SetBinContent(ich + 1, i + 1, 0);
        }
        else {
            for (Int_t i = 0; i < hinbad->GetNbinsY(); i++)
                hinbad->SetBinContent(ich + 1, i + 1, 0);
        }
    }
    
    printf("flagged %d problematic over %d checked channels\n", nprobch, ntotch);
    
    new TCanvas("cFraction");
    hFraction->Draw();
    dgaus->Draw("same");
    lcutfrac->Draw("same");
    lcutfrac2->Draw("same");
    gPad->Update();
    new TCanvas("cGood");
    hingood->Draw("colz");
    gPad->Update();
    new TCanvas("cBad");
    hinbad->Draw("colz");
    gPad->Update();
    new TCanvas;
    hin->Draw("colz");
    new TCanvas("hFEA");
    hFEAHitMap->Draw("colz");
    gPad->Update();
    
    TFile *fout = TFile::Open("problematicChannels.root", "RECREATE");
    hFraction->Write();
    hingood->Write();
    hinbad->Write();
    hFEAHitMap->Write("hBadFEAMap");
    fout->Close();
    
    /* create cdb info */
    AliCDBId id("TOF/Calib/Problematic", startRun, endRun);
    AliCDBMetaData *md = new AliCDBMetaData();
    md->SetResponsible("Roberto Preghenella");
    md->SetComment("Problematic");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    md->SetBeamPeriod(0);
    
    /* put object in cdb */
    AliCDBManager *cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage(dbString);
    cdb->GetDefaultStorage()->Put(obj, id, md);
    
}

//_____________________________________________________________

problematicFEAs(const Char_t *filename, Float_t nsigma = 3.)
{
    
    /* map FEA-channels */
    Float_t hitmap[6552][2];
    for (Int_t ich = 0; ich < 157248; ich++) {
        Int_t ifea = getIndex(kFEA, ich);
        getHitMapXY(kFEA, ich, hitmap[ifea]);
    }
    
    TFile *fin = TFile::Open(filename);
    TH2 *hin = (TH2 *)fin->Get("hCalibHisto_fea");
    TH2 *hingood = (TH2 *)hin->Clone("hingood");
    TH2 *hinbad = (TH2 *)hin->Clone("hinbad");
    
    TH1D *hinpy_all = hin->ProjectionY("hinpy_all");
    TF1 *fitFunc = (TF1 *)gROOT->GetFunction("gaus");
    fitPeak(fitFunc, hinpy_all, 2., 1.);
    Double_t mean = fitFunc->GetParameter(1);
    Double_t sigma = fitFunc->GetParameter(2);
    Double_t intmin = mean - nsigma * sigma;
    Double_t intmax = mean + nsigma * sigma;
    Int_t binmin = hinpy_all->FindBin(intmin);
    Int_t binmax = hinpy_all->FindBin(intmax);
    
    Double_t integral_all = hinpy_all->Integral(binmin, binmax);
    Double_t entries_all = hinpy_all->GetEntries();
    Double_t fraction_all = integral_all / entries_all;
    
    printf("found %f hits in +- 3 sigma\n", fraction_all);
    
    TH1D *hFraction = new TH1D("hFraction", "", 2000, 0., 2.);
    
    TH1C *obj = new TH1C("hProblematic", "", 157248, 0., 157248.);
    
    TH2F *hFEAHitMap = new TH2F("hFEAHitMap", "", 72, 0., 18., 91, 0., 91.);
    TH2F *hBadFEAHitMap = new TH2F("hBadFEAHitMap", "", 72, 0., 18., 91, 0., 91.);
    
    Float_t hitmap[2];
    
    TH1D *hpy;
    for (Int_t ich = 0; ich < 6552; ich++) {
        hpy = hin->ProjectionY("hpy", ich + 1, ich + 1);
        if (hpy->GetEntries() <= 0) {
            delete hpy;
            continue;
        }
        Double_t integral_ch = hpy->Integral(binmin, binmax);
        Double_t entries_ch = hpy->GetEntries();
        Double_t fraction_ch = integral_ch / entries_ch;// / fraction_all;
        hFraction->Fill(fraction_ch);
        
    }
    TF1 *gaus = (TF1 *)gROOT->GetFunction("gaus");
    hFraction->Fit(gaus, "WW");
    hFraction->Fit(gaus, "WW", "", gaus->GetParameter(1) - gaus->GetParameter(2), 2.);
    fraction_all = gaus->GetParameter(1);
    Float_t cutLimit = gaus->GetParameter(1) - 3. * gaus->GetParameter(2);
    cutLimit /= fraction_all;
    
    hFraction->Reset();
    for (Int_t ich = 0; ich < 6552; ich++) {
        hpy = hin->ProjectionY("hpy", ich + 1, ich + 1);
        if (hpy->GetEntries() <= 0) {
            delete hpy;
            continue;
        }
        Double_t integral_ch = hpy->Integral(binmin, binmax);
        Double_t entries_ch = hpy->GetEntries();
        Double_t fraction_ch = integral_ch / entries_ch / fraction_all;
        hFraction->Fill(fraction_ch);
        
        hFEAHitMap->Fill(hitmap[ich][0], hitmap[ich][1], fraction_ch);
        if (fraction_ch < cutLimit)
            hBadFEAHitMap->Fill(hitmap[ich][0], hitmap[ich][1]);
    }
    
    new TCanvas;
    hFraction->Draw();
    new TCanvas;
    hFEAHitMap->Draw("colz");
    new TCanvas;
    hBadFEAHitMap->Draw("colz");
    
    TFile *fout = TFile::Open("problematicFEAs.root", "RECREATE");
    hFraction->Write();
    hFEAHitMap->Write("hFEAFractionMap");
    hBadFEAHitMap->Write("hBadFEAMap");
    fout->Close();
    
    return;
    
    TF1 *gaus = (TF1 *)gROOT->GetFunction("gaus");
    hFraction->Fit(gaus, "WW");
    Double_t cutLimit = gaus->GetParameter(1) - 3. * gaus->GetParameter(2);
    
    for (Int_t ich = 0; ich < 157248; ich++) {
        hpy = hin->ProjectionY("hpy", ich + 1, ich + 1);
        if (hpy->GetEntries() <= 0) {
            delete hpy;
            continue;
        }
        Double_t integral_ch = hpy->Integral(binmin, binmax);
        Double_t entries_ch = hpy->GetEntries();
        Double_t fraction_ch = integral_ch / entries_ch / fraction_all;
        
        if (fraction_ch < cutLimit) {
            obj->SetBinContent(ich + 1, 0x1);
            getHitMapXY(kFEA, ich, hitmap);
            hFEAHitMap->Fill(hitmap[0], hitmap[1]);
            for (Int_t i = 0; i < hingood->GetNbinsY(); i++)
                hingood->SetBinContent(ich + 1, i + 1, 0);
        }
        else {
            for (Int_t i = 0; i < hinbad->GetNbinsY(); i++)
                hinbad->SetBinContent(ich + 1, i + 1, 0);
        }
    }
    
    hFraction->Draw();
    new TCanvas;
    hingood->Draw("colz");
    new TCanvas;
    hinbad->Draw("colz");
    new TCanvas;
    hin->Draw("colz");
    new TCanvas;
    hFEAHitMap->Draw("colz");
    
    TFile *fout = TFile::Open("problematicChannels.root", "RECREATE");
    hFraction->Write();
    hingood->Write();
    hinbad->Write();
    hFEAHitMap->Write("hBadFEAMap");
    fout->Close();
    
    /* create cdb info */
    AliCDBId id("TOF/Calib/Problematic", startRun, endRun);
    AliCDBMetaData *md = new AliCDBMetaData();
    md->SetResponsible("Roberto Preghenella");
    md->SetComment("Problematic");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    md->SetBeamPeriod(0);
    
    /* put object in cdb */
    AliCDBManager *cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage(dbString);
    cdb->GetDefaultStorage()->Put(obj, id, md);
    
}

//_____________________________________________________________

mergeProblematicMaps(const Char_t *mapfilename1, const Char_t *mapfilename2, Int_t startRun = 0, Int_t endRun = AliCDBRunRange::Infinity(), const Char_t *dbString = "local://$HOME/OCDB")
{
    TFile *file1 = TFile::Open(mapfilename1);
    TFile *file2 = TFile::Open(mapfilename2);
    
    AliCDBEntry *cdbe1 = (AliCDBEntry *)file1->Get("AliCDBEntry");
    AliCDBEntry *cdbe2 = (AliCDBEntry *)file2->Get("AliCDBEntry");
    
    TH1C *h1 = (TH1C *)cdbe1->GetObject();
    TH1C *h2 = (TH1C *)cdbe2->GetObject();
    TH1C *obj = new TH1C("hProblematic", "", 157248, 0., 157248.);
    
    for (Int_t ich = 0; ich < 157248; ich++) {
        if (h1->GetBinContent(ich + 1) || h2->GetBinContent(ich + 1))
	        obj->SetBinContent(ich + 1, 0x1);
    }
    
    /* create cdb info */
    AliCDBId id("TOF/Calib/Problematic", startRun, endRun);
    AliCDBMetaData *md = new AliCDBMetaData();
    md->SetResponsible("Roberto Preghenella");
    md->SetComment("Problematic");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    md->SetBeamPeriod(0);
    
    /* put object in cdb */
    AliCDBManager *cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage(dbString);
    cdb->GetDefaultStorage()->Put(obj, id, md);
}

//_____________________________________________________________

Int_t
fitPeak(TF1 *f, TH1D *h, Float_t nSigmaMin, Float_t nSigmaMax)
{
    /*
     * fitPeak
     */
    
    /* rebin if needed */
    Double_t ent = h->GetEntries();
    //  while (h->GetBinContent(h->GetMaximumBin()) < 0.1 * h->GetEntries())
    //    h->Rebin(2);
    Double_t fitCent = h->GetBinCenter(h->GetMaximumBin());
    Double_t fitMin = fitCent - nSigmaMin * 300.;
    Double_t fitMax = fitCent + nSigmaMax * 300.;
    if (fitMin < h->GetXaxis()->GetXmin()) fitMin = h->GetXaxis()->GetXmin();
    if (fitMax > h->GetXaxis()->GetXmax()) fitMax = h->GetXaxis()->GetXmax();
    f->SetParameter(1, fitCent);
    f->SetParameter(2, 150.);
    Int_t fitres = h->Fit(f, "WWq0", "", fitMin, fitMax);
    if (fitres != 0) return fitres;
    /* refit with better range */
    for (Int_t i = 0; i < 3; i++) {
        fitCent = f->GetParameter(1);
        fitMin = fitCent - nSigmaMin * f->GetParameter(2);
        fitMax = fitCent + nSigmaMax * f->GetParameter(2);
        if (fitMin < h->GetXaxis()->GetXmin()) fitMin = h->GetXaxis()->GetXmin();
        if (fitMax > h->GetXaxis()->GetXmax()) fitMax = h->GetXaxis()->GetXmax();
        fitres = h->Fit(f, "q0", "", fitMin, fitMax);
        if (fitres != 0) return fitres;
    }
    return fitres;
}

//_____________________________________________________________

void
monitor()
{
    /*
     * monitor
     */
    
    Float_t evfrac = (Float_t)curev / (Float_t)nevents;
    Float_t eltime = stopwatch.RealTime();
    stopwatch.Continue();
    Float_t tottime = evfrac > 0. ? eltime / evfrac : eltime;
    Float_t speed = eltime > 0. ? (Float_t)curev / eltime : 0.;
    printf("%4.2f \% | elapsed: %d s | total: %d s | remaining: %d s | speed: %d ev/s\n", 100. * evfrac, eltime, tottime, tottime - eltime, speed);
    
}

//_____________________________________________________________

Bool_t
updateParOfflineEntry(const Char_t *paramfilename, const Char_t *ocdbfilename, Int_t startrun, Int_t endrun = AliCDBRunRange::Infinity(), const Char_t *dbString = "alien://?folder=/alice/cern.ch/user/r/rpreghen/OCDB")
{
    
    TFile *paramfile = TFile::Open(paramfilename);
    if (!paramfile || !paramfile->IsOpen()) {
        return kTRUE;
    }
    TH1D *hshiftparams = (TH1D *)paramfile->Get("hChannelParam_mean");
    if (!hshiftparams) {
        return kTRUE;
    }
    TFile *ocdbfile = TFile::Open(ocdbfilename);
    if (!ocdbfile || !ocdbfile->IsOpen()) {
        return kTRUE;
    }
    AliCDBEntry *cdbe = (AliCDBEntry *)ocdbfile->Get("AliCDBEntry");
    if (!cdbe) {
        return kTRUE;
    }
    TObjArray *currentpararray = (TObjArray *)cdbe->GetObject();
    if (!currentpararray) {
        return kTRUE;
    }
    
    /* create calib structure */
    AliTOFcalib *tofcalib = new AliTOFcalib();
    tofcalib->CreateCalArrays();
    TObjArray *newpararray = (TObjArray *) tofcalib->GetTOFCalArrayOffline();
    AliTOFChannelOffline *currentparchannel, *newparchannel;
    Float_t currentslewpar, newslewpar;
    /* loop over channels */
    for (Int_t idx = 0; idx < 157248; idx++) {
        currentparchannel = (AliTOFChannelOffline *)currentpararray->At(idx);
        newparchannel = (AliTOFChannelOffline *)newpararray->At(idx);
        /* loop over slew params */
        for (Int_t ipar = 0; ipar < 6; ipar++) {
            currentslewpar = currentparchannel->GetSlewPar(ipar);
            if (ipar == 0) /* add shift to par[0] */
                newslewpar = currentslewpar + 1.e-3 * hshiftparams->GetBinContent(idx + 1);
            else /* leave other params unchanged */
                newslewpar = currentslewpar;
            newparchannel->SetSlewPar(ipar, newslewpar);
        }
    }
    
    AliCDBManager *man = AliCDBManager::Instance();
    man->SetDefaultStorage("local://$HOME/OCDB");
    man->SetSpecificStorage("TOF/Calib/ParOffline", dbString);
    tofcalib->WriteParOfflineOnCDB("TOF/Calib", "valid", startrun, endrun);
    
}

//_____________________________________________________________

Bool_t
updateRunParamsEntry(const Char_t *paramfilename, Int_t startrun, Int_t endrun = AliCDBRunRange::Infinity(), const Char_t *dbString = "alien://?folder=/alice/cern.ch/user/r/rpreghen/OCDB?user=rpreghen?se=ALICE::CERN::SE")
{
    
    TFile *paramfile = TFile::Open(paramfilename);
    if (!paramfile || !paramfile->IsOpen()) {
        return kTRUE;
    }
    TH1D *htimezerofill = (TH1D *)paramfile->Get("hTimeZeroFill_mean");
    if (!htimezerofill) {
        return kTRUE;
    }
    
    /* get t0-fill data */
    UInt_t timestamp[MAXPOINTS];
    Float_t t0[MAXPOINTS];
    Float_t t0spread[MAXPOINTS];
    Float_t tofreso[MAXPOINTS];
    Int_t npoints = 0;
    for (Int_t ibin = 0; ibin < htimezerofill->GetNbinsX(); ibin++) {
        if (htimezerofill->GetBinError(ibin + 1) == 0.) continue;
        timestamp[npoints] = (UInt_t)htimezerofill->GetBinCenter(ibin + 1);
        t0[npoints] = htimezerofill->GetBinContent(ibin + 1);
        t0spread[npoints] = -1.;
        tofreso[npoints] = -1.;
        npoints++;
    }
    
    /* setup dummy run */
    UInt_t run[1] = {-1};
    UInt_t runFirstPoint[1] = {0};
    UInt_t runLastPoint[1] = {npoints - 1};
    
    /* create runParams object */
    AliTOFRunParams obj(npoints, 1);
    obj.SetTimestamp(timestamp);
    obj.SetT0(t0);
    obj.SetTOFResolution(tofreso);
    obj.SetT0Spread(t0spread);
    obj.SetRunNb(run);
    obj.SetRunFirstPoint(runFirstPoint);
    obj.SetRunLastPoint(runLastPoint);
    
    /* install run params object in OCDB */
    AliCDBManager *cdb = AliCDBManager::Instance();
    AliCDBStorage *sto = cdb->GetStorage(dbString);
    if (!sto) {
        printf("cannot get storage %s", dbString);
        return kTRUE;
    }
    AliCDBId id("TOF/Calib/RunParams", startrun, endrun);
    AliCDBMetaData md;
    md.SetResponsible("Roberto Preghenella");
    md.SetComment("offline TOF run parameters (emergency entry, no distinction between runs)");
    md.SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    md.SetBeamPeriod(0);
    if (!sto->Put(&obj, id, &md)) {
        printf("error while putting object in storage %s", dbString);
        return kTRUE;
    }
    
    return kFALSE;
}

//_____________________________________________________________

inspectTRM(const Char_t *filename, Int_t trmIndex)
{
    
    TFile *fin = TFile::Open(filename);
    TFile *fout = TFile::Open(Form("inspectTRM_%d.%s", trmIndex, filename), "RECREATE");
    TH2 *hin = (TH2 *)fin->Get("hCalibHisto_channel");
    TH1 *hpy;
    
    for (Int_t i = 0; i < 157248; i++) {
        if (getIndex(kTRM, i) != trmIndex) continue;
        hpy = hin->ProjectionY(Form("hChannel_%d", i), i + 1, i + 1);
        fout->cd();
        hpy->Write();
    }
    
    fout->Close();
}
