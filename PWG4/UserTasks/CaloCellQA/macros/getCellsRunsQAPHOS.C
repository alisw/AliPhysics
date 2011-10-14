/** Welcome!
 *
 *  This macro is intended for the following tasks:
 *    1. find bad (switched off/dead/noisy/strange) cells;
 *    2. find out the status/quality of analysed data (missing RCUs, run quality);
 *    3. find the extent of problems related to bad cells: required for
 *       systematics estimation for a physical quantity related to a cluster.
 *
 *  For impatient, change the "infile" line below and execute
 *    root -l getCellsRunsQA.C
 *
 *  For curios, continue reading getCellsRunsQA() code: it is self-documented.
 *
 *  The macro options are tuned for a user (and pp runs), and in most cases no user
 *  intervension is necessary. Still, it is likely that you will have to edit
 *  nexc/excells[] in the parameters section below and run the macro several times.
 *  Consult with the output from this macro.
 *  In case of PbPb runs, a small modification is necessary:
 *   1) change ExcludeSmallRuns() line by putting a smaller number, e.g. 5-10k events;
 *   2) change FindDeadNoisyCellsPerRun() factor thresholds to a more narrow region, e.g. 0.07-2;
 *   3) probably, change fit region in FitPi0().
 *  Also, do not forget to adjust cluster cut for pi0s in AliCaloCellsQA. The value 2.5GeV
 *  is currently reasonable.
 *
 *  Generally, a QA expert uncomments all the functions which return (print to stdout)
 *  bad cell candidates and checks them by hand.
 *
 *  Detector-specific parts require to run this macro with aliroot instead of root,
 *  they are commented in getCellsRunsQA() by default.
 *
 *  This macro is written as a number of small functions and is designed both
 *  for EMCAL and for PHOS (and DCAL in future).
 *  Drawing options are chosen with respect to PPRStyle() drawing defaults.
 *
 *  Input: AliCaloCellsQA analysis output.
 *
 *  TODO: cells time spectra: currently it is not put in use. Seems that time shape fitting is
 *        not trivial due to the presence of parasite peaks (one needs to remove them first?),
 *        and this is a separate issue to think about...
 *
 *  TODO: some PHOS-specific parts
 *
 *  Author: Olga Driga (SUBATECH)
 */


// ========================== Parameters section ==============================

// input
// char *infile = "CellQA_LHC11d_pass1_PHI7.root";
char *infile = "CellQA_LHC11d_pass1_AnyInt.root";

// supermodule colors
Color_t colorsSM[] = {0,2,3,4};

// cells to exclude from averages calculations and their number: bad cells can
// mess up the results, it is suggested to list (known and newly found with
// this code) problematic cell candidates explicitly

Int_t excells[] = {1,2};
Int_t nexc = 2;

// ====================== End of parameters section ===========================

static TFile *gInputFile;

void getCellsRunsQAPHOS()
{
  // Entry point for the analysis.

  gRandom->SetSeed(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);
  gStyle->SetPalette(1);
  gStyle->SetFillColor(10);

  Printf("Input: %s", infile);
  gInputFile = TFile::Open(infile);

  // You may wish to extract and draw a particular histogram from infile.
  // Here are the examples:
  //
  //   TH1* hNEventsProcessedPerRun = (TH1*) gInputFile->Get("hNEventsProcessedPerRun");
  //   hNEventsProcessedPerRun->Draw();
  //   return;
  //
  //   TH1* hNTimesInClusterElow = (TH1*) gInputFile->Get("run146686_hCellLocMaxNTimesInClusterElow");
  //   hNTimesInClusterElow->Draw();
  //   return;
  //
  //   TH1* hETotalClusterEhigh = (TH1*) gInputFile->Get("run146686_hCellLocMaxETotalClusterEhigh");
  //   hETotalClusterEhigh->Draw();
  //   return;

  // Draw a random cell spectrum;
  // 0.2-1 GeV -- fit region;
  // 3 GeV -- the spectrum is drawn up to this value, -1 = no limit;
  // last argument -- histogram name to process, possible values are:
  //    hCellAmplitude, hCellAmplitudeEHigh, hCellAmplitudeNonLocMax or fhCellAmplitudeEhighNonLocMax.
  DrawCell(1+gRandom->Integer(10752), 0.25, 1., 4., "hCellAmplitude");  // common cell region for EMCAL2010 and for PHOS

  // Draw a random cell time spectrum
  // DrawCellTime(1+gRandom->Integer(10752));


  /* RUN NUMBER SELECTION SECTION
   *
   * NOTE: at any point below runs are sorted in chronological order.
   */

  // array with run numbers and their number
  Int_t runNumbers[10000];
  Int_t nruns = 0;

  // First, fill run numbers ...
  GetRunNumbers(nruns, runNumbers);
  Printf("Total number of runs: %i", nruns);

  // ... draw events distribution ...
  // (the last argument is number of bins in this distribution)
  DrawRunsDistribution(nruns, runNumbers, 100);

  // ... and exclude runs with number of events < 1k.
  ExcludeSmallRuns(nruns, runNumbers, 100);

  // You may wish to exclude particular runs:
  //   Int_t runs2Exclude[] = {111222,333444,555666};
  //   ExcludeRunNumbers(nruns, runNumbers, 3, runs2Exclude);

  Printf("Number of runs to be analysed: %i", nruns);

  // Finally, print a nice table with run index / run number / number of events.
  PrintRunNumbers(nruns, runNumbers);


  /* PER RUN BAD CELLS SEARCHING CRITERIA
   *
   * Four primary criteria on a per run basis:
   *   1 and 2: number of times cell was a local maximum in a cluster at low/high energies;
   *   3 and 4: total cluster energy for a local maximum cell at low/high energies.
   */

  // Extract the histograms with dead/noisy cell candidates per run.
  // For each of the four criteria, for each run and for each cell:
  //   1) calculate cell factor = [cell value]/[average over cells];
  //   2) mark cell as dead (candidate) if factor <= factorDead (3rd argument below);
  //   3) mark cell as noisy (candidate) if factor >= factorNoisy (4th argument below).
  //
  // Factor thresholds are quite wide by default:
  // low energy criteria are not very sensitive to them,
  // while high energy criteria are very sensitive due to limited statistics.
  //
  // The function below also draws histograms with factor distributions in all the runs
  // for all the cells. It may help to take the decision about dead/noisy factor thresholds.
  // The last two arguments -- number of bins and maximum X axis value for such histograms.
  //
  // Output: hBadCellMap[4] -- bad cell candidates per run for each of the four criteria;
  // X axis -- cell absId, Y axis -- run index, content: 0=not bad, 1=noisy, -1=dead.
  //

  TH2** hBadCellMap = FindDeadNoisyCellsPerRun(nruns, runNumbers, 0.05, 2.5, 200, 10.);
  
  // Look for noisy channels only:
  // TH2** hBadCellMap = FindDeadNoisyCellsPerRun(nruns, runNumbers, -1.0, 3.0, 200, 20.);

  // Look for dead channels only:
  // TH2** hBadCellMap = FindDeadNoisyCellsPerRun(nruns, runNumbers, 0.05, 1.e9, 200, 10.);

  // Print cell numbers suggested for exclusion from averages calculation;
  // see excells[] array in the parameters section in the beginning of the macro.
  // It is highly suggested to run this macro several times and add the output
  // of this function to the parameters section.
  PrintCellsToExcludeFromAverages(hBadCellMap);

  // The results can be visualized (second argument -- canvas name):
  // either per run ...
  DrawDeadNoisyCellMap(hBadCellMap[0], "hMapRuns_NTimesInClusterElow");
  DrawDeadNoisyCellMap(hBadCellMap[1], "hMapRuns_NTimesInClusterEhigh");
//   DrawDeadNoisyCellMap(hBadCellMap[2], "hMapRuns_ETotalClusterElow");
//   DrawDeadNoisyCellMap(hBadCellMap[3], "hMapRuns_ETotalClusterEhigh");

  // ... or per number of events ...
//   DrawDeadNoisyCellMap((TH2*)ConvertMapRuns2Events(nruns,runNumbers,hBadCellMap[0]), "hMapEvents_NTimesInClusterElow");
//   DrawDeadNoisyCellMap((TH2*)ConvertMapRuns2Events(nruns,runNumbers,hBadCellMap[1]), "hMapEvents_NTimesInClusterEhigh");
//   DrawDeadNoisyCellMap((TH2*)ConvertMapRuns2Events(nruns,runNumbers,hBadCellMap[2]), "hMapEvents_ETotalClusterElow");
//   DrawDeadNoisyCellMap((TH2*)ConvertMapRuns2Events(nruns,runNumbers,hBadCellMap[3]), "hMapEvents_ETotalClusterEhigh");

  // ... or we can also draw the percent for each cell being dead/noisy
//   DrawDeadNoisyCellPercent(nruns, runNumbers, hBadCellMap[0], "hPercent_NTimesInClusterElow");
//   DrawDeadNoisyCellPercent(nruns, runNumbers, hBadCellMap[1], "hPercent_NTimesInClusterEhigh");
//   DrawDeadNoisyCellPercent(nruns, runNumbers, hBadCellMap[2], "hPercent_ETotalClusterElow");
//   DrawDeadNoisyCellPercent(nruns, runNumbers, hBadCellMap[3], "hPercent_ETotalClusterEhigh");

  // Our main criteria for analytical finding of bad cells is based on the following.
  // Note that factor distributions for NTimesInCluster and ETotalCluster are very
  // similar, both at low and at high energy. Note also that we can say the same
  // for dead/noisy cell visualizations above: they are very similar. This suggests
  // that NTimesInCluster and ETotalCluster give the same results. The criteria
  // NTimesInCluster, however, is more calibration-independent (though energy thresholds
  // are still calibration-dependent) and thus is more reliable and clear. Thus, we
  // limit ourselves to NTimesInClusterElow/NTimesInClusterEhigh criteria.
  // Now, you could probably noted that NTimesInClusterEhigh give more dead
  // cells than that at low energy. This is an expected statistical effect: we have
  // much smaller number of clusters at high energy. Consequently, we will not use dead
  // cell candidates at high energy.
  //
  // Finally, we combine candidates from low/high energies and produce one TH2
  // histogram which is the primary source of our analytical results.
  //
  TH2* hBadCellMapPrimary = CombineDeadNoisyElowEhigh(hBadCellMap[0], hBadCellMap[1]);

  // Note for PHOS: if you are not happy with NTimesInClusterEhigh results (due to a lack of statistics)
  // uncomment this line:
//   hBadCellMapPrimary = hBadCellMap[0];

  // Draw everything for combined
  DrawDeadNoisyCellMap(hBadCellMapPrimary, "Primary_hMapRuns");
  DrawDeadNoisyCellMap((TH2*)ConvertMapRuns2Events(nruns,runNumbers,hBadCellMapPrimary), "Primary_hMapEvents");
  DrawDeadNoisyCellPercent(nruns, runNumbers, hBadCellMapPrimary, "Primary_hPercent");

  // Print full information on cells which are dead/noisy;
  // kTRUE -- print also the percentage % dead/% noisy
  PrintDeadNoisyCells(hBadCellMapPrimary, 0.9, 1.);          // in 90-100% of runs
  // PrintDeadNoisyCells(hBadCellMapPrimary, 0.5, 0.9, kTRUE);  // in 50-90% of runs
  // PrintDeadNoisyCells(hBadCellMapPrimary, 0.3, 0.5);      // in 30-50% of runs
  PrintDeadNoisyCells(hBadCellMapPrimary, 0.1, 0.5);      // in 10-50% of runs

  // visualize dead/noisy cell map for EMCAL/PHOS; requires aliroot
  DrawOccupancy(nruns, runNumbers, hBadCellMapPrimary, "hDeadNoisyCellsOccupancy");

  // EMCAL: print full information on missing/noisy parts (e.g. RCUs); requires aliroot
//   PrintEMCALProblematicBlocks(nruns, runNumbers, hBadCellMapPrimary);


  /* RUNS QUALITY SECTION: CLUSTER AVERAGES + PI0 AVERAGES
   *
   */

  // First, draw cluster averages per run;
  //   1 = minimum number of cells in a cluster;
  //   0.3GeV = minimum cluster energy;
  //   -1     = maximum cluster energy = infinity (in fact, 20GeV in the default configuration).
  DrawClusterAveragesPerRun(nruns, runNumbers, 1, 0.3, -1);

  // Second, draw the same cluster averages per run, but corrected for detector acceptance
  DrawClusterAveragesPerRun(nruns, runNumbers, 1, 0.3, -1, hBadCellMapPrimary);

  // Draw random slices of the pi0 peak in one supermodule and in whole detector
  DrawPi0Slice(runNumbers[6],  1);
  DrawPi0Slice(runNumbers[6],  2);
  DrawPi0Slice(runNumbers[6],  3);
  DrawPi0Slice(runNumbers[gRandom->Integer(nruns)], -1);

  // Draw number of pi0s per event, pi0 mass and width
  DrawPi0Averages(nruns, runNumbers);

  // Draw pi0 values per event with pi0 photons both in the same supermodule
//   DrawPi0Averages(nruns, runNumbers, kTRUE);

  // Draw pi0 values per event with pi0 photons both in the same supermodule
  // + correct for supermodule acceptance (should not be treated to be 100% reliable!)
//   DrawPi0Averages(nruns, runNumbers, kTRUE, hBadCellMap[0]);

  // Draw pi0 values per event with pi0 
  // + correct for supermodule acceptance (should not be treated to be 100% reliable!)
  // DrawPi0Averages(nruns, runNumbers, kFALSE, hBadCellMap[0]);

  // Draw pi0 distributions (helps to take decision on bad runs);
  // first argument -- suffix from canvas title;
  // second argument -- number of bins in distributions
//   DrawPi0Distributions("", 100);
//   DrawPi0Distributions("SM1", 100)
//   DrawPi0Distributions("_sameSM", 100);
//   DrawPi0Distributions("_sameSM_corr4accep", 100);

  // !!!
  return;

  /* SHAPE ANALYSIS
   *
   * Lazy/curious boundary is here! -------------------------
   * Do not pass it if you do not fill curious enough! ;)
   *
   * The shape analysis below belongs to the main bad cell searching criteria.
   * However, it may require some extra work in order to make it usable.
   *
   * The idea behind is simple: fit the energy spectrum of a cell with
   * the function A*exp(-B*x)/x^2 (which proved to be a very good fit),
   * draw parameters A, B and chi2/ndf distributions and notice
   * cells which deviate too much from the averages.
   *
   * Huge statistics (>20M events) is necessary for this to work reliably.
   *
   * TODO: test on PHOS
   */

  // The analysis below defines a cell as being bad simply if it is outside
  // of some good region for any of the parameters A, B, chi2/ndf.
  // The regions are depicted by vertical orange lines. The problem is that these
  // regions are usually not automatically determined correctly.
  //
  // The function syntax is the following:
  //   0.1-1.0 GeV -- fitting region;
  //   hCellAmplitude -- the histogram which to take for processing,
  //                     the other possible choice is hCellAmplitudeNonLocMax;
  //   1000 -- number of bins is distributions.
  //
  // The next three groups of parameters are:
  //   <text label> <maximum value on distribution> <left edge of the good region> <right edge of the good region>
  //
  // It is left/right edges which usually require manual tuning.
  // -1 means to initialize a parameter automatically.
  //
  TestCellShapes(0.25, 1., "hCellAmplitude", 1000,
               // maxval / left margin / right margin
                 -1, -1,-1,     // fit A
                 -1, -1,-1,     // fit B
                 -1, -1,-1);    // fit chi2/ndf


  /* DISTRIBUTION ANALYSIS
   *
   * Test for bad cells by plotting a distribution among cells.
   *
   * It is especially useful in searching for miscalibrated cells.
   * The function parameters are similar to parameters from shape analysis section.
   */

  TestCellsTH1(nruns, runNumbers, "hCellLocMaxNTimesInClusterElow",
               "Number of times cell was local maximum, low energy, per event", "Entries",
               400, -1,-1,-1);   // nbins / maxval / left margin / right margin

  TestCellsTH1(nruns, runNumbers, "hCellLocMaxNTimesInClusterEhigh",
               "Number of times cell was local maximum, high energy, per event", "Entries",
               400, -1,-1,-1);   // nbins / maxval / left margin / right margin

//   TestCellsTH1(nruns, runNumbers, "hCellLocMaxETotalClusterElow",
//                "Total cluster energy for local maximum cell, low energy, per event", "Entries",
//                400, -1,-1,-1);   // nbins / maxval / left margin / right margin
//
//   TestCellsTH1(nruns, runNumbers, "hCellLocMaxETotalClusterEhigh",
//                "Total cluster energy for local maximum cell, high energy, per event", "Entries",
//                400, -1,-1,-1);   // nbins / maxval / left margin / right margin

  // Three more tests for bad cells:
  //  1) total deposited energy;
  //  2) total number of entries;
  //  3) average energy = [total deposited energy]/[total number of entries].
  //
  TestCellEandN(0.1, "hCellAmplitude", 1000,
             // maxval / left margin / right margin
                -1,-1,-1,    // cell E
                -1,-1,-1,    // cell N
                -1,-1,-1);   // cell E/N

  // the same at high energies
  TestCellEandN(0.1, "hCellAmplitudeEhigh", 1000,
             // maxval / left margin / right margin
                -1,-1,-1,    // cell E
                -1,-1,-1,    // cell N
                -1,-1,-1);   // cell E/N


  // The same as above, but for not a local maximum cells; require more statistics
//   TestCellsTH1(nruns, runNumbers, "hCellNonLocMaxNTimesInClusterElow",
//                "Number of times cell wasn't local maximum, low energy, per event", "Entries",
//                400, -1,-1,-1);   // nbins / maxval / left margin / right margin
//
//   TestCellsTH1(nruns, runNumbers, "hCellNonLocMaxNTimesInClusterEhigh",
//                "Number of times cell wasn't local maximum, high energy, per event", "Entries",
//                400, -1,-1,-1);   // nbins / maxval / left margin / right margin
//
//   TestCellsTH1(nruns, runNumbers, "hCellNonLocMaxETotalClusterElow",
//                "Total cluster energy for not local maximum cell, low energy, per event", "Entries",
//                400, -1,-1,-1);   // nbins / maxval / left margin / right margin
//
//   TestCellsTH1(nruns, runNumbers, "hCellNonLocMaxETotalClusterEhigh",
//                "Total cluster energy for not local maximum cell, high energy, per event", "Entries",
//                400, -1,-1,-1);   // nbins / maxval / left margin / right margin
//
//   TestCellEandN(0.1, "hCellAmplitudeNonLocMax", 1000,
//              // maxval / left margin / right margin
//                 -1,-1,-1,    // cell E
//                 -1,-1,-1,    // cell N
//                 -1,-1,-1);   // cell E/N
//
//   TestCellEandN(0.1, "hCellAmplitudeEhighNonLocMax", 1000,
//              // maxval / left margin / right margin
//                 -1,-1,-1,    // cell E
//                 -1,-1,-1,    // cell N
//                 -1,-1,-1);   // cell E/N


  // TODO: cells time

  // The end ;)
}



/* BAD CELLS SEARCHING FUNCTIONS
 *
 */

//_________________________________________________________________________
TH2** FindDeadNoisyCellsPerRun(const Int_t nruns, Int_t runNumbers[],
                               Double_t factorDead = 0.01, Double_t factorNoisy = 4.,
                               Int_t fnbins = 200, Double_t fxmax = 10.)
{
  // Return bad cell candidate maps for the four criteria;
  // X axis -- cell absId, Y axis -- run index, content: 0=not bad, 1=noisy, -1=dead.
  //
  // For each run and each cell calculate factor = [cell value]/[averate over cells],
  // mark cell as noisy if factor >= factorNoisy, mark cell as dead if factor <= factorDead.
  //
  // fnbins, fxmax -- number of bins and X axis maximum value for factor distributions.
  //
  // Four criteria in order:
  char *hname[4] = {"hCellLocMaxNTimesInClusterElow", "hCellLocMaxNTimesInClusterEhigh",
                     "hCellLocMaxETotalClusterElow", "hCellLocMaxETotalClusterEhigh"};

  TH1*  hFactorDistr[4];
  TH2** hBadCellMap = new TH2*[4];

  // take one histogram to get binning parameters
  TH1* one = (TH1*) gInputFile->Get(Form("run%i_%s",runNumbers[0],hname[0]));
  Int_t  ncells = one->GetNbinsX();
  Double_t amin = one->GetXaxis()->GetXmin();
  Double_t amax = one->GetXaxis()->GetXmax();

  // find dead/noisy cell candidates
  for (Int_t k = 0; k < 4; k++) {
    hBadCellMap[k] = new TH2C(Form("hBadCellMap_%s_fdead=%.3f_fnoisy=%.1f",hname[k],factorDead,factorNoisy),
                              Form("Dead/noisy cell map (%s, fdead=%.3f, fnoisy=%.1f)",hname[k],factorDead,factorNoisy),
                                    ncells,amin,amax, nruns,0,nruns);
    hBadCellMap[k]->SetXTitle("AbsId");
    hBadCellMap[k]->SetYTitle("Run index");
    hBadCellMap[k]->SetTitleOffset(1.7,"Y");

    hFactorDistr[k] = new TH1F(Form("hFactorDistr_%s_fdead=%.3f_fnoisy=%.1f",
                                     hname[k],factorDead,factorNoisy), "", fnbins,0,fxmax);
    hFactorDistr[k]->SetTitle("Factor distributions in all runs");
    hFactorDistr[k]->SetXTitle("Factor");
    hFactorDistr[k]->SetYTitle("Entries");

    // run index loop
    for (Int_t ri = 0; ri < nruns; ri++) {
      TH1* one = (TH1*) gInputFile->Get(Form("run%i_%s", runNumbers[ri], hname[k]));
      if (one == 0) {
	Error("FindDeadNoisyCellsPerRun","Run %d does contain histogram %s",runNumbers[ri],hname[k]);
	continue;
      }
      // calculate average
      Double_t av = 0;  // average
      Int_t count = 0;  // counted cells
      for (Int_t c = 1; c <= ncells; c++) {
        // do not count cells with zero content
        if (one->GetBinContent(c) == 0) continue;
        // do not count cells listed in the parameters section in the beginning of the macro
        if (IsCellMarked4Exclude(one->GetBinLowEdge(c))) continue;
        count++;
        av += one->GetBinContent(c);
      }

      // division by zero checks
      if (count == 0) {
        Warning("FindDeadNoisyCellsPerRun", Form("No cells counted at ri=%i",ri));
        continue;
      }
      av /= count;

      if (av == 0) {
        Warning("FindDeadNoisyCellsPerRun", Form("Average is zero at ri=%i",ri));
        continue;
      }

      // find dead/noisy candidates
      for (Int_t c = 1; c <= ncells; c++) {
        Double_t fac = one->GetBinContent(c)/av;
        hFactorDistr[k]->Fill(fac);

        if (fac <= factorDead)
          hBadCellMap[k]->SetBinContent(c, ri+1, -1);
        else if (fac >= factorNoisy)
          hBadCellMap[k]->SetBinContent(c, ri+1, 1);
      }

      delete one;
    } // run index loop
  } // criteria loop


  // draw factor distributions ...
  TCanvas *c1 = new TCanvas("hFactorDistr", "hFactorDistr", 400,400);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.08);
  gPad->SetLogy();
  hFactorDistr[0]->SetLineColor(kBlue);
  hFactorDistr[1]->SetLineColor(kRed);
  hFactorDistr[2]->SetLineColor(kGreen);
  hFactorDistr[3]->SetLineColor(kOrange);
  hFactorDistr[0]->Draw();
  hFactorDistr[1]->Draw("same");
  hFactorDistr[2]->Draw("same");
  hFactorDistr[3]->Draw("same");

  // ... with legend
  TLegend *leg = new TLegend(0.45,0.65,0.90,0.85);
  leg->AddEntry(hFactorDistr[0], "NTimesInCluster, low energy","l");
  leg->AddEntry(hFactorDistr[2], "ETotalCluster, low energy","l");
  leg->AddEntry(hFactorDistr[1], "NTimesInCluster, high energy","l");
  leg->AddEntry(hFactorDistr[3], "EtotalCluster, high energy","l");
  leg->Draw("same");

  c1->Update();

  return hBadCellMap;
}

//_________________________________________________________________________
void PrintCellsToExcludeFromAverages(TH2** hBadCellMap)
{
  // Print cells suggested for exclusion from averages calculation

  Int_t ncells = hBadCellMap[0]->GetNbinsX();
  Int_t nruns  = hBadCellMap[0]->GetNbinsY();

  Int_t *suggested = new Int_t[ncells];
  memset(suggested, 0, ncells*sizeof(Int_t));

  for (Int_t c = 1; c <= ncells; c++)
    for (Int_t ri = 1; ri <= nruns; ri++)
      if (hBadCellMap[0]->GetBinContent(c, ri) != 0 || hBadCellMap[2]->GetBinContent(c, ri) != 0 ||
          hBadCellMap[1]->GetBinContent(c, ri)  > 0 || hBadCellMap[3]->GetBinContent(c, ri)  > 0) // NOTE: dead not counted
        suggested[c-1]++;

  printf("Suggested cells to switch off in averages calculations (copai 'n paste!):\n");
  printf("Int_t excells[] = {");
  Int_t n = 0;
  for (Int_t c = 1; c <= ncells; c++)
    if (suggested[c-1] >= 0.4*nruns) {// 40% of runs threshold
      printf("%s%i", n == 0 ? "" : ",", c);
      n++;
    }
  printf("};\nInt_t nexc = %i;\n\n",n);
}

//_________________________________________________________________________
Bool_t IsCellMarked4Exclude(Int_t absId)
{
  // Return true if a cell is in excells[] array

  static TH1* hExclCells = NULL;

  // one time initialization
  if (!hExclCells) {
    hExclCells = new TH1F("hExclCells", "", 20000,0,20000);
    for (Int_t c = 0; c < nexc; c++)
      hExclCells->SetBinContent(hExclCells->FindBin(excells[c]), 1);
  }

  return (hExclCells->GetBinContent(hExclCells->FindBin(absId)) > 0 ? kTRUE : kFALSE);
}

//_________________________________________________________________________
void DrawDeadNoisyCellMap(TH2* hmap, char* cname)
{
  // Visualize dead/noisy cell map;
  // cname -- canvas name.

  // Define only 3 color: dead=blue, noisy=red, normal=white
  Int_t color_index[]={kBlue,kWhite,kRed};
  gStyle->SetPalette(3,color_index); 

  TCanvas *c1 = new TCanvas(cname, cname, 0,0,1000,600);
  c1->Divide(1,3);
  TH2* hDummy[3];

  // Use hardware numeration of modules (2,3,4)
  // instead of offline numeration (3,2,1)
  for (Int_t iM=1; iM<=3; iM++) {
    c1->cd(4-iM);
    gPad->SetLeftMargin(0.04);
    gPad->SetRightMargin(0.02);
    gPad->SetTopMargin(0.10);
    gPad->SetBottomMargin(0.13);
    hDummy[iM-1] = (TH2*) hmap->Clone(Form("hDummy_%s_mod%d",cname,iM));
    hDummy[iM-1]->SetTitle(Form("Module %d",5-iM));
    hDummy[iM-1]->SetLabelSize(0.08,"XY");
    hDummy[iM-1]->SetTitleSize(0.08,"XY");
    hDummy[iM-1]->SetTitleOffset(0.8,"X");
    hDummy[iM-1]->SetTitleOffset(0.25,"Y");
    hDummy[iM-1]->SetTickLength(0.01,"Y");
    hDummy[iM-1]->SetNdivisions(520,"X");
    hDummy[iM-1]->SetAxisRange(3584*(iM-1)+1, 3584*iM, "X");
    hDummy[iM-1]->Draw("col");
  }

  c1->Update();
}

// //_________________________________________________________________________
// void DrawDeadNoisyCellMap(TH2* hmap, char* cname)
// {
//   // Visualize dead/noisy cell map;
//   // cname -- canvas name.

//   TCanvas *c1 = new TCanvas(cname, cname, 0,0,600,600);
//   gPad->SetLeftMargin(0.14);
//   gPad->SetRightMargin(0.06);

//   // draw dummy to initialize axis ranges
//   TH2* hDummy = (TH2*) hmap->Clone(Form("hDummy_%s",cname));
//   hDummy->Reset();
//   hDummy->Draw();

//   for (Int_t c = 1; c <= hmap->GetNbinsX(); c++) //cell number
//     for (Int_t ri = 1; ri <= hmap->GetNbinsY(); ri++) { //run index
//       Double_t stat = hmap->GetBinContent(c, ri); // cell status
//       if (stat == 0) continue;

//       Double_t x  = hmap->GetXaxis()->GetBinCenter(c);
//       Double_t y1 = hmap->GetYaxis()->GetBinLowEdge(ri);
//       Double_t y2 = hmap->GetYaxis()->GetBinUpEdge(ri);

//       // draw a line; FIXME: what is a better choice?
//       TLine* line = new TLine(x,y1,x,y2);
//       line->SetLineWidth(1);
//       if (stat > 0) line->SetLineColor(kRed); // noisy cell
//       else line->SetLineColor(kBlue); // dead cell
//       line->Draw();
//     }

//   c1->Update();
// }

//_________________________________________________________________________
void DrawDeadNoisyCellPercent(Int_t nruns, Int_t runNumbers[], TH2* hmap, char* cname)
{
  // Show percent of runs/events for each cell being problematic;
  // cname -- canvas name.

  // binning parameters
  Int_t  ncells = hmap->GetNbinsX();
  Double_t amin = hmap->GetXaxis()->GetXmin();
  Double_t amax = hmap->GetXaxis()->GetXmax();

  // number of times cell was dead/noisy;
  Int_t *ndeadRuns    = new Int_t[ncells];
  Int_t *nnoisyRuns   = new Int_t[ncells];
  Double_t *ndeadEvents  = new Double_t[ncells];
  Double_t *nnoisyEvents = new Double_t[ncells];

  // fill arrays above
  for (Int_t c = 0; c < ncells; c++) {
    ndeadRuns[c] = 0;
    nnoisyRuns[c] = 0;
    ndeadEvents[c] = 0;
    nnoisyEvents[c] = 0;

    for (Int_t ri = 0; ri < nruns; ri++) {
      Double_t stat = hmap->GetBinContent(c+1, ri+1); // cell status
      Int_t nevents = GetNumberOfEvents(runNumbers[ri]);

      if (stat > 0) {
        nnoisyRuns[c]++;
        nnoisyEvents[c] += nevents;
      }
      else if (stat < 0) {
        ndeadRuns[c]++;
        ndeadEvents[c] += nevents;
      }
    } // run index loop
  } // cell loop

  // total number of events
  Double_t ntotal = GetTotalNumberOfEvents(nruns, runNumbers);

  TCanvas *c1 = new TCanvas(cname, cname, 0,0,600,600);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.06);

  // draw dummy histogram to initialize canvas
  TH1* hDummy = new TH1F(Form("hDummy_%s",cname), hmap->GetTitle(), ncells,amin,amax);
  hDummy->SetAxisRange(0,100, "Y");
  hDummy->SetXTitle("AbsId");
  hDummy->SetYTitle("Percent");
  hDummy->Draw();

  // to fill legend
  TLine *line1 = NULL;
  TLine *line2 = NULL;
  TLine *line3 = NULL;
  TLine *line4 = NULL;

  // draw results, FIXME: is where a better way?
  for (Int_t c = 0; c < ncells; c++) {
    Double_t x = hmap->GetXaxis()->GetBinCenter(c+1);
    Double_t y1 = 100.*ndeadEvents[c]/ntotal;
    Double_t y2 = 100.*ndeadRuns[c]/nruns;

    // events, dead bar
    if (ndeadEvents[c] > 0) {
      line1 = new TLine(x, 0, x, y1);
      line1->SetLineWidth(1);
      line1->SetLineColor(kBlue);
      line1->Draw();
    }

    // events, noisy bar
    if (nnoisyEvents[c] > 0) {
      line2 = new TLine(x, y1, x, y1 + 100.*nnoisyEvents[c]/ntotal);
      line2->SetLineWidth(1);
      line2->SetLineColor(kRed);
      line2->Draw();
    }

    // runs, dead bar
    if (ndeadRuns[c] > 0) {
      line3 = new TLine(x, 0, x, y2);
      line3->SetLineWidth(1);
      line3->SetLineStyle(7);
      line3->SetLineColor(kBlack);
      line3->Draw();
    }

    // runs, noisy bar
    if (nnoisyRuns[c] > 0) {
      line4 = new TLine(x, y2, x, y2 + 100.*nnoisyRuns[c]/nruns);
      line4->SetLineWidth(1);
      line4->SetLineStyle(7);
      line4->SetLineColor(kOrange);
      line4->Draw();
    }
  }

  // legend
  TLegend *leg = new TLegend(0.65,0.65,0.92,0.75);
  if (line1) leg->AddEntry(line1, "% of events, dead","l");
  if (line2) leg->AddEntry(line2, "% of events, noisy","l");
  if (line3) leg->AddEntry(line3, "% of runs, dead","l");
  if (line4) leg->AddEntry(line4, "% of runs, noisy","l");
  leg->Draw("same");

  c1->Update();
}

//_________________________________________________________________________
TH1* ConvertMapRuns2Events(const Int_t nruns, Int_t runNumbers[], TH1* inhisto)
{
  // Returns a histogram in which run index axis is converted to number of events
  // by making a variable axis bin width.
  // If inhisto inherits from TH2, convert Y axis; convert X axis otherwise.

  // bin widths
  Double_t *nevents = new Double_t[nruns+1];
  nevents[0] = 0;
  for (Int_t ri = 0; ri < nruns; ri++)
    nevents[1+ri] = nevents[ri] + GetNumberOfEvents(runNumbers[ri]);

  TH1* outh = (TH1*) inhisto->Clone(Form("%s_Events",inhisto->GetName()));

  if (inhisto->InheritsFrom("TH2")) {
    outh->GetYaxis()->Set(nruns, nevents);
    outh->SetYTitle("Number of events");
    outh->SetTitleOffset(1.7,"Y");
  }
  else {// TH1 case
    outh->GetXaxis()->Set(nruns, nevents);
    outh->SetXTitle("Number of events");
  }

  return outh;
}

//_________________________________________________________________________
TH2* CombineDeadNoisyElowEhigh(TH2* hmapElow, TH2* hmapEhigh)
{
  // Combine two maps at low/high energy into one,
  // do not count dead map at high energy.
  // NOTE: if cell is dead at Elow and noisy and Ehigh, set content = 2.

  TH2* hmap_combined = (TH2*) hmapElow->Clone(Form("%s+%s",hmapElow->GetName(),hmapEhigh->GetName()));

  for (Int_t c = 1; c <= hmapElow->GetNbinsX(); c++)
    for (Int_t ri = 1; ri <= hmapElow->GetNbinsY(); ri++)
      if (hmapEhigh->GetBinContent(c, ri) > 0) {
        hmap_combined->SetBinContent(c, ri, 1);

        if (hmapElow->GetBinContent(c, ri) < 0)
          hmap_combined->SetBinContent(c, ri, 2);
      }

  hmap_combined->SetTitle("Dead/noisy cell map, combined");

  return hmap_combined;
}

//_________________________________________________________________________
void PrintDeadNoisyCells(TH2* hmap, Double_t percent1 = 0.95, Double_t percent2 = 1., Bool_t kPrintPercentage = kFALSE)
{
  // Print full information on dead/noisy cells which are present in
  // percent1-percent2 portion of runs (percent1 excluded, percent2 included).

  Int_t ncells = hmap->GetNbinsX();
  Int_t nruns  = hmap->GetNbinsY();

  Int_t nbad = 0;

  printf("Dead/noisy cells in >%.1f%% and <=%.1f%% of runs:", 100*percent1, 100*percent2);

  for (Int_t c = 1; c <= ncells; c++) {
    Int_t nrdead = 0; // count number of runs for the current cell
    Int_t nrnoisy = 0;

    for (Int_t ri = 1; ri <= nruns; ri++)
      if      (hmap->GetBinContent(c,ri) < 0) nrdead++;
      else if (hmap->GetBinContent(c,ri) > 0) nrnoisy++;

    if (nrdead+nrnoisy > percent1*nruns && nrdead+nrnoisy <= percent2*nruns) {
      printf(" %.0f", hmap->GetBinLowEdge(c));
      if (kPrintPercentage)
        printf("(%-2.0f/%-2.0f)", 100.*nrdead/nruns, 100.*nrnoisy/nruns);
      nbad++;
    }
  }

  printf(" (%i total)\n\n", nbad);
}

//_________________________________________________________________________
void DrawOccupancy(Int_t nruns, Int_t runNumbers[], TH2* hmap, char* cname)
{
  // Draw bad cell map for EMCAL or PHOS;
  // cname -- canvas name.

  gStyle->SetPalette(1);

  // guess detector
  if (hmap->GetNbinsX() % 1152 == 0)
    DrawEMCALOccupancy(nruns, runNumbers, hmap, cname);
  else
    DrawPHOSOccupancy(nruns, runNumbers, hmap, cname);
}

//_________________________________________________________________________
void DrawEMCALOccupancy(Int_t nruns, Int_t runNumbers[], TH2* hmap, char* cname)
{
  // Draw bad cell map for EMCAL;
  // cname -- canvas name.

  Int_t nmods = hmap->GetNbinsX()/1152;  // number of supermodules
  Int_t vsize = ceil(nmods/2.);          // vertical size in canvas
  Int_t lastSector = (nmods-1)/2;        // for pad number calculation

  TCanvas *c1 = new TCanvas(cname, cname, 800,200*vsize);
  c1->Divide(2, vsize);

  for (Int_t sm = 0; sm < nmods; sm++)
  {
    // top left is SM0, bottom right is SM9
    Int_t side = sm%2;
    Int_t sector = sm/2;
    c1->cd((lastSector - sector)*2 + side + 1);

    TH2* hSM = new TH2C(Form("hEMCALSM%i_%s",sm,cname), Form("SM%i",sm), 48,0,48, 24,0,24);
    hSM->SetXTitle("iEta");
    hSM->SetYTitle("iPhi");

    // loop over supermodule cells
    for (Int_t c = 0; c < 1152; c++) {
      Int_t absId = 1152 * sm + c;

      for (Int_t ri = 0; ri < nruns; ri++) {
        if (hmap->GetBinContent(hmap->FindBin(absId,ri)) == 0) continue;

        Int_t nModule, eta, phi;
        AbsId2EtaPhi(absId, nModule, eta, phi);
        hSM->Fill(eta,phi);
      }
    }

    hSM->Draw("colz");
  } // supermodule loop

  c1->Update();
}

//_________________________________________________________________________
void DrawPHOSOccupancy(Int_t nruns, Int_t runNumbers[], TH2* hmap, char* cname)
{
  // Draw bad cell map for PHOS;
  // cname -- canvas name.

  Int_t nmods = hmap->GetNbinsX()/3584;  // number of supermodules
  Int_t vsize = nmods;                   // vertical size in canvas

  TCanvas *c1 = new TCanvas(cname, cname, 64*10*vsize,56*10);
  c1->Divide(vsize,1);

  TFile *fBadMap = TFile::Open("PHOS_BadMap.root","recreate");
  for (Int_t sm = 1; sm <= nmods; sm++)
  {
    c1->cd(sm);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.15);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.10);

    TH2* hSM = new TH2C(Form("hPHOSSM%i_%s",sm,cname), Form("Module %i",5-sm), 64,0.5,64.5, 56,0.5,56.5);
    hSM->SetXTitle("x_{cell}");
    hSM->SetYTitle("z_{cell}      ");
    hSM->SetZTitle("N_{runs}  ");

    // loop over supermodule cells
    for (Int_t c = 1; c <= 3584; c++) {
      Int_t absId = 3584*(sm-1) + c;

      for (Int_t ri = 0; ri < nruns; ri++) {
        if (hmap->GetBinContent(hmap->FindBin(absId,ri)) == 0) continue;

        Int_t nModule, xCell, zCell;
        AbsId2EtaPhi(absId, nModule, xCell, zCell, 1);
        hSM->Fill(xCell,zCell);
      }
    }
    hSM->Write();
    hSM->DrawCopy("colz");
  } // supermodule loop

  fBadMap->Close();
  c1->Update();
}

//_________________________________________________________________________
void PrintEMCALProblematicBlocks(Int_t nruns, Int_t runNumbers[], TH2* hmap)
{
  // Print, on a per run basis, complete information about EMCAL missing
  // (or dead or noisy) blocks. Missing/noisy EMCAL atomic part is a 2x2
  // block (288 blocks per supermodule).

  // number of supermodules
  Int_t nmods = hmap->GetNbinsX()/1152;

  Printf("Problematic (missing/dead/noisy) 2x2 block numbers in EMCAL:");

  // run index loop
  for (Int_t ri = 0; ri < nruns; ri++) {
    Printf("Run %i:", runNumbers[ri]);

    // supermodule loop
    for (Int_t sm = 0; sm < nmods; sm++) {

      // will be filled with the number of missing cells (from 0 to 4)
      Int_t blk2x2[288];
      memset(blk2x2, 0, 288*sizeof(Int_t));

      for (Int_t c = 0; c < 1152; c++) {
        Int_t absId = 1152 * sm + c;

        // select problematic cells only
        if (hmap->GetBinContent(hmap->FindBin(absId,ri)) == 0) continue;

        Int_t nModule, eta, phi;
        AbsId2EtaPhi(c, nModule, eta, phi);
        blk2x2[nModule]++;
      }

      // calculate number of bad blocks
      Int_t nbad2x2 = 0;
      for (Int_t b = 0; b < 288; b++)
        if (blk2x2[b] == 4) nbad2x2++;

      // no bad
      if (nbad2x2 == 0) continue;

      printf("  SM%i:", sm);

      // whole supermodule
      if (nbad2x2 == 288) {
        printf(" missing the whole supermodule!\n");
        continue;
      }

      // RCUs
      if (nbad2x2 >= 144) {
        Int_t nRCU[2];
        nRCU[0] = 0;
        nRCU[1] = 0;

        for (Int_t b = 0; b < 288; b++)
          if (blk2x2[b] == 4) nRCU[GetEMCALRCUNumber(b)]++;

        if (nRCU[0] == 144) {
          printf(" RCU0");
          for (Int_t b = 0; b < 288; b++)
            if (blk2x2[b] == 4 && GetEMCALRCUNumber(b) == 0) blk2x2[b] = 0;
        }

        if (nRCU[1] == 144) {
          printf(" RCU1");
          for (Int_t b = 0; b < 288; b++)
            if (blk2x2[b] == 4 && GetEMCALRCUNumber(b) == 1) blk2x2[b] = 0;
        }

        nbad2x2 -= 144;
      }

      // the rest
      if (nbad2x2 > 0) {
        for (Int_t b = 0; b < 288; b++)
          if (blk2x2[b] == 4) printf(" %i", b);
        printf(" (%i)", nbad2x2);
      }

      printf("\n");

    } // supermodule loop
  } // run index loop

  Printf("");
}

//_________________________________________________________________________
Int_t GetEMCALRCUNumber(Int_t nModule)
{
  // Returns RCU number for a 2x2 block in EMCAL;
  // nModule -- block number (0-287).

  static TH1* hRCUs = NULL;

  // one-time initialization
  if (!hRCUs) {
    hRCUs = new TH1F("hRCU1", "", 288,0,288);

    Int_t RCU1[144] = {8,9,10,11,20,21,22,23,32,33,34,35,44,45,46,47,56,57,58,59,68,69,70,71,80,81,82,83,
             92,93,94,95,104,105,106,107,116,117,118,119,128,129,130,131,140,141,142,143,148,149,150,151,
             152,153,154,155,160,161,162,163,164,165,166,167,172,173,174,175,176,177,178,179,184,185,186,
             187,188,189,190,191,196,197,198,199,200,201,202,203,208,209,210,211,212,213,214,215,220,221,
             222,223,224,225,226,227,232,233,234,235,236,237,238,239,244,245,246,247,248,249,250,251,256,
             257,258,259,260,261,262,263,268,269,270,271,272,273,274,275,280,281,282,283,284,285,286,287};

    for (Int_t i = 0; i < 144; i++) hRCUs->SetBinContent(RCU1[i]+1, 1);
  }

  return hRCUs->GetBinContent(nModule+1);
}

//_________________________________________________________________________
void AbsId2EtaPhi(Int_t absId, Int_t &nModule, Int_t &eta, Int_t &phi, Int_t det = 0)
{
  // Converts cell absId --> (sm,eta,phi);
  //
  // nModule -- 2x2 block number for EMCAL;
  // det -- detector: 0/EMCAL, 1/PHOS.

  // EMCAL
  if (det == 0) {
    AliEMCALGeometry *geomEMCAL = AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");

    Int_t nSupMod, nIphi, nIeta;
    geomEMCAL->GetCellIndex(absId, nSupMod, nModule, nIphi, nIeta);
    geomEMCAL->GetCellPhiEtaIndexInSModule(nSupMod, nModule, nIphi, nIeta, phi, eta);
    return;
  }

  // PHOS
  else if (det == 1) {
    AliPHOSGeometry *geomPHOS = AliPHOSGeometry::GetInstance("IHEP");

    Int_t relid[4];
    geomPHOS->AbsToRelNumbering(absId, relid);
    //sm = relid[0];
    eta = relid[2];
    phi = relid[3];
    return;
  }

  // DCAL
  // not implemented
  //

  Error("AbsId2EtaPhi", "Wrong detector");
  abort();
}

//_________________________________________________________________________
void TestCellsTH1(Int_t nruns, Int_t runNumbers[], char *hname,
                  char* htitle = "", char* hytitle = "",
                  Int_t dnbins = 200, Double_t dmaxval = -1, Double_t goodmin = -1, Double_t goodmax = -1)
{
  // Test for bad cells by plotting a distribution of a TH1 histogram.
  // The histogram is obtained as a sum over runs of TH1 per run.
  //
  // hname -- histogram name to process;
  // htitle, hytitle -- histogram title and Y axis title;
  // dnbins -- number of bins in distribution;
  // dmaxval -- X axis maximum in distribution.
  // goodmin,goodmax -- the region outside which a cell is considered as bad.
  //
  // -1 value for maxval/goodmin/goodmax -- process a variable automatically.

  // initialize histogram
  TH1* histo = (TH1*) gInputFile->Get(Form("run%i_%s", runNumbers[0], hname));
  histo->SetName(hname);
  histo->SetTitle(htitle);
  histo->SetXTitle("AbsId");
  histo->SetYTitle(hytitle);

  // sum over runs
  for (Int_t i = 1; i < nruns; i++) {
    TH1* h = (TH1*) gInputFile->Get(Form("run%i_%s", runNumbers[i], hname));
    histo->Add(h);
    delete h;
  }

  histo->Scale(1./(Double_t)GetTotalNumberOfEvents(nruns, runNumbers));
  Process1(histo, Form("TestCellsTH1_%s",hname), dnbins, dmaxval, goodmin, goodmax);
}

//_________________________________________________________________________
void TestCellEandN(Double_t Emin = 0.1, char* hname = "hCellAmplitude", Int_t dnbins = 200,
                   Double_t maxval1 = -1, Double_t goodmin1 = -1, Double_t goodmax1 = -1,
                   Double_t maxval2 = -1, Double_t goodmin2 = -1, Double_t goodmax2 = -1,
                   Double_t maxval3 = -1, Double_t goodmin3 = -1, Double_t goodmax3 = -1)
{
  // Three more tests for bad cells:
  //  1) total deposited energy;
  //  2) total number of entries;
  //  3) average energy = [total deposited energy]/[total number of entries].
  //
  // Based on summary histograms. Possible choises:
  //   hCellAmplitude, hCellAmplitudeEhigh, hCellAmplitudeNonLocMax, hCellAmplitudeEhighNonLocMax
  //
  // Emin -- minimum cell amplitude to count;
  // hname -- name (in file) of TH2 histogram with cell amplitudes;
  // dnbins -- number of bins in distributions;
  // maxval[123] -- maximum values on distributions for the criteria 1),2),3) respectively;
  // goodmin[123],goodmax[123] -- the regions on distributions outside those a cell is considered as bad.

  // input; X axis -- absId numbers
  TH2 *hCellAmplitude = (TH2*) gInputFile->Get(hname);

  // binning parameters
  Int_t ncells = hCellAmplitude->GetNbinsX();
  Double_t amin = hCellAmplitude->GetXaxis()->GetXmin();
  Double_t amax = hCellAmplitude->GetXaxis()->GetXmax();

  TH1* hCellEtotal = new TH1F(Form("%s_hCellEtotal_E%.2f",hname,Emin),
                              Form("Total deposited energy, E > %.2f GeV",Emin), ncells,amin,amax);
  hCellEtotal->SetXTitle("AbsId");
  hCellEtotal->SetYTitle("Energy, GeV");

  TH1F *hCellNtotal = new TH1F(Form("%s_hCellNtotal_E%.2f",hname,Emin),
                               Form("Total number of entries, E > %.2f GeV",Emin), ncells,amin,amax);
  hCellNtotal->SetXTitle("AbsId");
  hCellNtotal->SetYTitle("Entries");

  TH1F *hCellEtoNtotal = new TH1F(Form("%s_hCellEtoNtotal_E%.2f",hname,Emin),
                                  Form("Average energy per hit, E > %.2f GeV",Emin), ncells,amin,amax);
  hCellEtoNtotal->SetXTitle("AbsId");
  hCellEtoNtotal->SetYTitle("Energy, GeV");

  // fill cells
  for (Int_t c = 1; c <= ncells; c++) {
    Double_t Esum = 0;
    Double_t Nsum = 0;

    for (Int_t i = 1; i <= hCellAmplitude->GetNbinsY(); i++) {
      Double_t E = hCellAmplitude->GetYaxis()->GetBinCenter(i);
      Double_t N = hCellAmplitude->GetBinContent(c, i);
      if (E < Emin) continue;
      Esum += E*N;
      Nsum += N;
    }

    hCellEtotal->SetBinContent(c, Esum);
    hCellNtotal->SetBinContent(c, Nsum);

    if (Nsum > 0.5)  // number of entries >= 1
      hCellEtoNtotal->SetBinContent(c, Esum/Nsum);
  }

  delete hCellAmplitude;

  Process1(hCellEtotal,    Form("%s_CellE", hname),   dnbins, maxval1, goodmin1, goodmax1);
  Process1(hCellNtotal,    Form("%s_CellN", hname),   dnbins, maxval2, goodmin2, goodmax2);
  Process1(hCellEtoNtotal, Form("%s_CellE/N", hname), dnbins, maxval3, goodmin3, goodmax3);
}

//_________________________________________________________________________
void TestCellShapes(Double_t fitEmin, Double_t fitEmax, char* hname = "hCellAmplitude", Int_t dnbins = 1000,
                    Double_t maxval1 = -1, Double_t goodmin1 = -1, Double_t goodmax1 = -1,
                    Double_t maxval2 = -1, Double_t goodmin2 = -1, Double_t goodmax2 = -1,
                    Double_t maxval3 = -1, Double_t goodmin3 = -1, Double_t goodmax3 = -1)
{
  // Test cells shape using fit function f(x)=A*exp(-B*x)/x^2.
  // Produce values per cell + distributions for A, B and chi2/ndf parameters.
  //
  // fitEmin, fitEmax -- fit range;
  // hname -- name (in file) of TH2 histogram with cell amplitudes;
  // dnbins -- number of bins in distributions;
  // maxval[123] -- maximum values on distributions for the criteria A, B and chi2/ndf respectively;
  // goodmin[123],goodmax[123] -- the regions on distributions outside those a cell is considered as bad.
  //
  // -1 value for maxval/goodmin/goodmax -- process a variable automatically.
  //
  // Note: numbers are optimized for EMCAL.
  // TODO: check for PHOS

  // input; X axis -- absId numbers
  TH2 *hCellAmplitude = (TH2*) gInputFile->Get(hname);

  // binning parameters
  Int_t  ncells = hCellAmplitude->GetNbinsX();
  Double_t amin = hCellAmplitude->GetXaxis()->GetXmin();
  Double_t amax = hCellAmplitude->GetXaxis()->GetXmax();

  // initialize histograms
  TH1 *hFitA = new TH1F(Form("hFitA_%s",hname),"Fit A value", ncells,amin,amax);
  hFitA->SetXTitle("AbsId");
  hFitA->SetYTitle("A");

  TH1 *hFitB = new TH1F(Form("hFitB_%s",hname),"Fit B value", ncells,amin,amax);
  hFitB->SetXTitle("AbsId");
  hFitB->SetYTitle("B");

  TH1 *hFitChi2Ndf = new TH1F(Form("hFitChi2Ndf_%s",hname),"Fit #chi^{2}/ndf value", ncells,amin,amax);
  hFitChi2Ndf->SetXTitle("AbsId");
  hFitChi2Ndf->SetYTitle("#chi^{2}/ndf");

  // total number of events; to estimate A value
  TH1* hNEventsProcessedPerRun = (TH1*) gInputFile->Get("hNEventsProcessedPerRun");
  Double_t ntotal = hNEventsProcessedPerRun->Integral(1, hNEventsProcessedPerRun->GetNbinsX());

  // fitting function
  TF1 *fit = new TF1("fit", "[0]*exp(-[1]*x)/x^2");
  fit->SetParLimits(0, ntotal*1e-8,ntotal*1e-4);
  fit->SetParLimits(1, 0.,30.);
  fit->SetParameter(0, ntotal*1e-6);
  fit->SetParameter(1, 1.5);

  for (Int_t i = 1; i <= ncells; i++) {
    TH1 *hCell = hCellAmplitude->ProjectionY("",i,i);
    if (hCell->GetEntries() == 0) continue;

    hCell->Fit(fit, "0LQEM", "", fitEmin, fitEmax);
    delete hCell;

    hFitA->SetBinContent(i, fit->GetParameter(0));
    hFitB->SetBinContent(i, fit->GetParameter(1));
    if (fit->GetNDF() != 0)
      hFitChi2Ndf->SetBinContent(i, fit->GetChisquare()/fit->GetNDF());
  }

  delete hCellAmplitude;

  // automatic parameters, if requested
  if (maxval1 < 0) maxval1 = 4e-6 * ntotal;
  if (maxval2 < 0) maxval2 = 10.;
  if (maxval3 < 0) maxval3 = 15.;

  Process1(hFitA,       Form("%s_FitA", hname),       dnbins, maxval1, goodmin1, goodmax1);
  Process1(hFitB,       Form("%s_FitB", hname),       dnbins, maxval2, goodmin2, goodmax2);
  Process1(hFitChi2Ndf, Form("%s_FitChi2ndf", hname), dnbins, maxval3, goodmin3, goodmax3);
}

//_________________________________________________________________________
void Process1(TH1* inhisto, char* label = "", Int_t dnbins = 200,
              Double_t dmaxval = -1, Double_t goodmin = -1, Double_t goodmax = -1)
{
  // Actual distribution analysis for a TH1 histogram:
  //  1) create a distribution for the input histogram;
  //  2) draw nicely;
  //  3) take a decision about bad cells.
  //
  // inhisto -- input histogram;
  // label -- text label for stdout;
  // dnbins -- number of bins in distribution;
  // goodmin,goodmax -- cells outside this region are considered as bad;
  // dmaxval -- maximum value on distribution histogram.
  // The later is required in cases where a bad cell kills the whole distribution:
  // limiting distribution maximum value solves the problem.

  if (dmaxval < 0)
    dmaxval = inhisto->GetMaximum()*1.01;  // 1.01 - to see the last bin

  TH1 *distrib = new TH1F(Form("%sDistr",inhisto->GetName()), "", dnbins, inhisto->GetMinimum(), dmaxval);
  distrib->SetXTitle(inhisto->GetYaxis()->GetTitle());
  distrib->SetYTitle("Entries");

  // fill distribution
  for (Int_t c = 1; c <= inhisto->GetNbinsX(); c++)
    distrib->Fill(inhisto->GetBinContent(c));

  // draw histogram + distribution
  TCanvas *c1 = new TCanvas(inhisto->GetName(),inhisto->GetName(), 800,400);
  c1->Divide(2,1);

  c1->cd(1);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.06);
  gPad->SetLogy();
  inhisto->SetTitleOffset(1.7,"Y");
  inhisto->Draw();

  c1->cd(2);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.10);
  gPad->SetLogy();
  distrib->Draw();

  // simple way to estimate the left margin for good cells region:
  // go from left to right and find the first bin with content 2,
  // then go from this bin right to left while bin content is nonzero
  if (goodmin < 0) {
    goodmin = distrib->GetXaxis()->GetXmin();

    for (Int_t i = 1; i <= distrib->GetNbinsX(); i++)
      if (distrib->GetBinContent(i) == 2) {
        while (i > 1 && distrib->GetBinContent(i-1) > 0) i--;
        goodmin = distrib->GetBinLowEdge(i);
        break;
      }
  }

  // the same automatic algorithm as above, but reflected
  if (goodmax < 0) {
    goodmax = distrib->GetXaxis()->GetXmax();

    for (Int_t i = distrib->GetNbinsX(); i >= 1; i--)
      if (distrib->GetBinContent(i) == 2) {
        while (i < distrib->GetNbinsX() && distrib->GetBinContent(i+1) > 0) i++;
        goodmax = distrib->GetXaxis()->GetBinUpEdge(i);
        break;
      }
  }

  // lines
  TLine *lline = new TLine(goodmin, 0, goodmin, distrib->GetMaximum());
  lline->SetLineColor(kOrange);
  lline->SetLineStyle(7);
  lline->Draw();

  TLine *rline = new TLine(goodmax, 0, goodmax, distrib->GetMaximum());
  rline->SetLineColor(kOrange);
  rline->SetLineStyle(7);
  rline->Draw();

  // legend
  TLegend *leg = new TLegend(0.60,0.82,0.9,0.88);
  leg->AddEntry(lline, "Good region boundary","l");
  leg->Draw("same");

  c1->Update();


  printf("Bad cells by criterum \"%s\":", label);

  // print bad cell numbers (outside goodmin,goodmax region)
  Int_t ntot = 0;
  for (Int_t c = 1; c <= inhisto->GetNbinsX(); c++)
    if (inhisto->GetBinContent(c) < goodmin || inhisto->GetBinContent(c) > goodmax) {
      printf(" %.0f", inhisto->GetBinLowEdge(c));
      ntot++;
    }

  printf(" (%i total)\n\n", ntot);
}

//_________________________________________________________________________
void DrawCell(Int_t absId, Double_t fitEmin = 0.3, Double_t fitEmax = 1.,
              Double_t Emax = -1, char* hname = "hCellAmplitude")
{
  // Draw one cell spectrum with a fit.
  //
  // fitEmin, fitEmax -- fit range;
  // Emax -- maximum value on X axis to show (-1 = no limit);
  // hname -- TH2 histogram name to process, possible values:
  //    "hCellAmplitude", "hCellAmplitudeEHigh", "hCellAmplitudeNonLocMax", "hCellAmplitudeEhighNonLocMax".

  // input; X axis -- absId numbers
  TH2* hCellAmplitude = (TH2*) gInputFile->Get(hname);

  Int_t bin = hCellAmplitude->GetXaxis()->FindBin(absId);
  TH1* hCell = hCellAmplitude->ProjectionY(Form("hCell%i_%s",absId,hname),bin,bin);
  hCell->SetXTitle("Energy, GeV");
  hCell->SetYTitle("Entries");
  hCell->SetTitle(Form("Cell %i, %s", absId, hname));
  if (Emax > 0) hCell->SetAxisRange(0, Emax);

  // draw spectrum
  TCanvas *c1 = new TCanvas(Form("hCell%i_%s",absId,hname), Form("hCell%i_%s",absId,hname), 400,400);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.08);
  gPad->SetBottomMargin(0.12);
  gPad->SetLogy();
  hCell->Draw();

  // total number of events
  TH1* hNEventsProcessedPerRun = (TH1*) gInputFile->Get("hNEventsProcessedPerRun");
  Double_t ntotal = hNEventsProcessedPerRun->Integral(1, hNEventsProcessedPerRun->GetNbinsX());

  // fit
  TF1 *fit = new TF1("fit", "[0]*exp(-[1]*x)/x^2");
  fit->SetLineColor(kRed);
  fit->SetLineWidth(2);
  fit->SetParName(0, "A");
  fit->SetParName(1, "B");
  fit->SetParLimits(0, ntotal*1e-8,ntotal*1e-4);
  fit->SetParLimits(1, 0.,30.);
  fit->SetParameter(0, ntotal*1e-6);
  fit->SetParameter(1, 1.);
  hCell->Fit(fit,"LQEM", "", fitEmin, fitEmax);

  // legend
  TLegend *leg = new TLegend(0.5,0.75,0.9,0.8);
  leg->AddEntry(fit, "A*exp(-B*x)/x^{2}","l");
  leg->Draw("same");

  c1->Update();
}

//_________________________________________________________________________
void DrawCellTime(Int_t absId)
{
  // Draw one cell time spectrum

  // input; X axis -- absId numbers
  TH2* hCellTime = (TH2*) gInputFile->Get("hCellTime");

  Int_t bin = hCellTime->GetXaxis()->FindBin(absId);
  TH1* hCell = hCellTime->ProjectionY(Form("hCellTime%i",absId),bin,bin);
  hCell->SetXTitle("Time, s");
  hCell->SetYTitle("Entries");
  hCell->SetTitle(Form("Cell %i time", absId));

  // draw spectrum
  TCanvas *c1 = new TCanvas(Form("hCellTime%i",absId), Form("hCellTime%i",absId), 400,400);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.10);
  gPad->SetBottomMargin(0.12);
  gPad->SetLogy();
  hCell->Draw();

  c1->Update();
}


/* RUN AVERAGES AND RELATED FUNCTIONS
 *
 */

//_________________________________________________________________________
void DrawClusterAveragesPerRun(Int_t nruns, Int_t runNumbers[], Int_t ncellsMin = 1,
                               Double_t eclusMin = 0.3, Double_t eclusMax = -1,
                               TH2* hmap = NULL)
{
  // Draws cluster averages per run vs run index and number of events.
  // NOTE: the averages are "smoothed" a little due to a finite bin width.
  //
  // ncellsMin -- minimum number of cells in cluster (>= 1);
  // eclusMin,eclusMax -- minimum/maximum cluster energy cut
  //                      (eclusMax = -1 means infinity);
  // hmap -- dead/noisy cell map, which is used for acceptance calculation due
  //   to switched off detector parts. Acceptance is used to correct the average
  //   number of clusters per event.

  // names suffix
  TString s(Form("_NC%i_Emin=%.2fGeV",ncellsMin,eclusMin));
  if (eclusMax > 0) s += TString(Form("_Emax=%.2fGeV",eclusMax)).Data();
  if (hmap) s += "_corr4accept";
  char *suff = s.Data();

  if (eclusMax < 0) eclusMax = 1e+5;

  // supermodule region
  Int_t SM1 = 0;
  Int_t SM2 = 10;
  while (SM1 <= SM2 && !gInputFile->Get(Form("run%i_hNCellsInClusterSM%i",runNumbers[0],SM1)) ) SM1++;
  while (SM2 >= SM1 && !gInputFile->Get(Form("run%i_hNCellsInClusterSM%i",runNumbers[0],SM2)) ) SM2--;

  // initialize histograms
  hAvECluster = new TH1F(Form("hAvECluster%s",suff), "Average cluster energy", nruns,0,nruns);
  SetRunLabel(hAvECluster,nruns,runNumbers,1);
  // hAvECluster->SetXTitle("Run index");
  hAvECluster->SetTitleOffset(0.3,"Y");
  hAvECluster->SetYTitle("Energy, GeV");
  hAvECluster->SetTickLength(0.01,"Y");
  hAvECluster->SetLabelSize(0.06,"XY");
  hAvECluster->SetTitleSize(0.06,"XY");

  hAvNCluster = new TH1F(Form("hAvNCluster%s",suff), "Average number of clusters per event", nruns,0,nruns);
  SetRunLabel(hAvNCluster,nruns,runNumbers,1);
  // hAvNCluster->SetXTitle("Run index");
  hAvNCluster->SetTitleOffset(0.3,"Y");
  hAvNCluster->SetYTitle("Number of clusters");
  hAvNCluster->SetTickLength(0.01,"Y");
  hAvNCluster->SetLabelSize(0.06,"XY");
  hAvNCluster->SetTitleSize(0.06,"XY");

  hAvNCellsInCluster = new TH1F(Form("hAvNCellsInCluster%s",suff), "Average number of cells in cluster", nruns,0,nruns);
  SetRunLabel(hAvNCellsInCluster,nruns,runNumbers,1);
  // hAvNCellsInCluster->SetXTitle("Run index");
  hAvNCellsInCluster->SetTitleOffset(0.3,"Y");
  hAvNCellsInCluster->SetYTitle("Number of cells");
  hAvNCellsInCluster->SetTickLength(0.01,"Y");
  hAvNCellsInCluster->SetLabelSize(0.06,"XY");
  hAvNCellsInCluster->SetTitleSize(0.06,"XY");

  // initialize per SM histograms
  TH1* hAvEClusterSM[10];
  TH1* hAvNClusterSM[10];
  TH1* hAvNCellsInClusterSM[10];

  for (Int_t sm = SM1; sm <= SM2; sm++) {
    hAvEClusterSM[sm] = new TH1F(Form("hAvEClusterSM%i%s",sm,suff),"", nruns,0,nruns);
    hAvNClusterSM[sm] = new TH1F(Form("hAvNClusterSM%i%s",sm,suff),"", nruns,0,nruns);
    hAvNCellsInClusterSM[sm] = new TH1F(Form("hAvNCellsInClusterSM%i%s",sm,suff),"", nruns,0,nruns);
  }

  // fill all the histograms per run index
  for (Int_t ri = 0; ri < nruns; ri++)
  {
    Int_t nevents = GetNumberOfEvents(runNumbers[ri]);

    // number of switched off supermodules
    Int_t noSM = 0;

    Double_t Eclus_total = 0;  // total cluster energy
    Double_t Nclus_total = 0;  // total number of clusters
    Double_t Ncells_total = 0; // total number of cells

    // supermodule loop
    for (Int_t sm = SM1; sm <= SM2; sm++) {
      TH2* hNCellsInClusterSM = (TH2*) gInputFile->Get(Form("run%i_hNCellsInClusterSM%i",runNumbers[ri],sm));
      if (hNCellsInClusterSM == 0) {
	Error("DrawClusterAveragesPerRun","Run %d does contain histogram %s",runNumbers[ri],Form("run%i_hNCellsInClusterSM%i",runNumbers[ri],sm));
	continue;
      }

      // the same as above, but per supermodule
      Double_t Eclus_totalSM = 0;
      Double_t Nclus_totalSM = 0;
      Double_t Ncells_totalSM = 0;

      // X axis -- cluster energy, Y axis -- number of cells
      for (Int_t x = 1; x <= hNCellsInClusterSM->GetNbinsX(); x++)
        for (Int_t y = 1+ncellsMin; y <= hNCellsInClusterSM->GetNbinsY(); y++) {//NOTE: bin 1 correspond to ncellsMin=0
          Double_t Eclus = hNCellsInClusterSM->GetXaxis()->GetBinCenter(x);
          Double_t Ncells = hNCellsInClusterSM->GetYaxis()->GetBinLowEdge(y);
          Double_t Nclus = hNCellsInClusterSM->GetBinContent(x,y);

          // cut on cluster energy
          if (Eclus < eclusMin || Eclus > eclusMax) continue;

          Eclus_totalSM += Eclus * Nclus;
          Nclus_totalSM += Nclus;
          Ncells_totalSM += Ncells * Nclus;
        }

        delete hNCellsInClusterSM;

      // correct for acceptance
      if (hmap) {
        Double_t accep = GetAcceptance(sm, hmap, ri);
        if (accep > 0) {
          Eclus_totalSM /= accep;
          Nclus_totalSM /= accep;
          Ncells_totalSM /= accep;
        }
        else noSM++;
      }

      Eclus_total += Eclus_totalSM;
      Nclus_total += Nclus_totalSM;
      Ncells_total += Ncells_totalSM;

      hAvNClusterSM[sm]->SetBinContent(ri+1, Nclus_totalSM/nevents);
      if (Nclus_totalSM > 0) hAvEClusterSM[sm]->SetBinContent(ri+1, Eclus_totalSM/Nclus_totalSM);
      if (Nclus_totalSM > 0) hAvNCellsInClusterSM[sm]->SetBinContent(ri+1, Ncells_totalSM/Nclus_totalSM);
    } // supermodule loop

    hAvNCluster->SetBinContent(ri+1, Nclus_total/nevents/(SM2-SM1+1-noSM));
    if (Nclus_total > 0) hAvECluster->SetBinContent(ri+1, Eclus_total/Nclus_total);
    if (Nclus_total > 0) hAvNCellsInCluster->SetBinContent(ri+1, Ncells_total/Nclus_total);
  } // run index loop


  /* Draw results vs run index
   */

  TCanvas *c1 = new TCanvas(Form("ClusterAveragesRuns%s",suff),
                            Form("ClusterAveragesRuns%s",suff), 1000,707);
  c1->Divide(1,3);

  // average cluster energy
  c1->cd(1);
  gPad->SetLeftMargin(0.04);
  gPad->SetRightMargin(0.02);
  gPad->SetTopMargin(0.10);
  gPad->SetBottomMargin(0.13);
  TLegend *leg = new TLegend(0.65,0.15,0.85,0.15+0.04*(SM2-SM1+1));

  hAvECluster->SetAxisRange(hAvECluster->GetMinimum()*0.5, hAvECluster->GetMaximum()*1.5,"Y");
  hAvECluster->SetLineWidth(2);
  hAvECluster->Draw();
  leg->AddEntry(hAvECluster, "All Modules","l");
  for (Int_t sm = SM1; sm <= SM2; sm++) {
    hAvEClusterSM[sm]->SetLineColor(colorsSM[sm]);
    hAvEClusterSM[sm]->SetLineWidth(1);
    hAvEClusterSM[sm]->Draw("same");
    leg->AddEntry(hAvEClusterSM[sm], Form("Module %i",5-sm),"l");
  }
  hAvECluster->Draw("same"); // to top
  leg->Draw("same");

  // average number of clusters
  c1->cd(2);
  gPad->SetLeftMargin(0.04);
  gPad->SetRightMargin(0.02);
  gPad->SetTopMargin(0.10);
  gPad->SetBottomMargin(0.13);
  leg = new TLegend(0.65,0.15,0.85,0.15+0.06*(SM2-SM1+1));

  hAvNCluster->SetAxisRange(0, hAvNCluster->GetMaximum()*1.5,"Y");
  hAvNCluster->SetLineWidth(2);
  hAvNCluster->Draw();
  leg->AddEntry(hAvNCluster, Form("(All Moduls)/%i",SM2-SM1+1),"l");
  for (Int_t sm = SM1; sm <= SM2; sm++) {
    hAvNClusterSM[sm]->SetLineColor(colorsSM[sm]);
    hAvNClusterSM[sm]->SetLineWidth(1);
    hAvNClusterSM[sm]->Draw("same");
    leg->AddEntry(hAvNClusterSM[sm], Form("Module %i",5-sm),"l");
  }
  hAvNCluster->Draw("same"); // to top
  leg->Draw("same");

  // average number of cells in cluster
  c1->cd(3);
  gPad->SetLeftMargin(0.04);
  gPad->SetRightMargin(0.02);
  gPad->SetTopMargin(0.10);
  gPad->SetBottomMargin(0.13);
  leg = new TLegend(0.65,0.15,0.85,0.15+0.04*(SM2-SM1+1));

  hAvNCellsInCluster->SetAxisRange(0, hAvNCellsInCluster->GetMaximum()*1.3,"Y");
  hAvNCellsInCluster->SetLineWidth(2);
  hAvNCellsInCluster->Draw();
  leg->AddEntry(hAvNCellsInCluster, "All Modules","l");
  for (Int_t sm = SM1; sm <= SM2; sm++) {
    hAvNCellsInClusterSM[sm]->SetLineColor(colorsSM[sm]);
    hAvNCellsInClusterSM[sm]->SetLineWidth(1);
    hAvNCellsInClusterSM[sm]->Draw("same");
    leg->AddEntry(hAvNCellsInClusterSM[sm], Form("Module %i",5-sm),"l");
  }
  hAvNCellsInCluster->Draw("same"); // to top
  leg->Draw("same");


  /* Draw the same vs number of events
   */

  TCanvas *c2 = new TCanvas(Form("ClusterAveragesEvents%s",suff),
                            Form("ClusterAveragesEvents%s",suff), 1000,707);
  c2->Divide(1,3);

  // average cluster energy
  c2->cd(1);
  gPad->SetLeftMargin(0.04);
  gPad->SetRightMargin(0.02);
  gPad->SetTopMargin(0.10);
  gPad->SetBottomMargin(0.13);
  leg = new TLegend(0.65,0.15,0.85,0.15+0.04*(SM2-SM1+1));

  TH1* hAvEClusterEvents = ConvertMapRuns2Events(nruns, runNumbers, hAvECluster);
  hAvEClusterEvents->Draw();
  leg->AddEntry(hAvEClusterEvents, "All Modules","l");
  for (Int_t sm = SM1; sm <= SM2; sm++) {
    TH1* hAvEClusterSMEvents = ConvertMapRuns2Events(nruns, runNumbers, hAvEClusterSM[sm]);
    hAvEClusterSMEvents->Draw("same");
    leg->AddEntry(hAvEClusterSMEvents, Form("Module %i",5-sm),"l");
  }
  hAvEClusterEvents->Draw("same"); // to top
  leg->Draw("same");

  // average number of clusters
  c2->cd(2);
  gPad->SetLeftMargin(0.04);
  gPad->SetRightMargin(0.02);
  gPad->SetTopMargin(0.10);
  gPad->SetBottomMargin(0.13);
  leg = new TLegend(0.65,0.15,0.85,0.15+0.04*(SM2-SM1+1));

  TH1* hAvNClusterEvents = ConvertMapRuns2Events(nruns, runNumbers, hAvNCluster);
  hAvNClusterEvents->Draw();
  leg->AddEntry(hAvNClusterEvents, Form("(All Modules)/%i",SM2-SM1+1),"l");
  for (Int_t sm = SM1; sm <= SM2; sm++) {
    TH1* hAvNClusterSMEvents = ConvertMapRuns2Events(nruns, runNumbers, hAvNClusterSM[sm]);
    hAvNClusterSMEvents->Draw("same");
    leg->AddEntry(hAvNClusterSMEvents, Form("Module %i",5-sm),"l");
  }
  hAvNClusterEvents->Draw("same"); // to top
  leg->Draw("same");

  // average number of cells in cluster
  c2->cd(3);
  gPad->SetLeftMargin(0.04);
  gPad->SetRightMargin(0.02);
  gPad->SetTopMargin(0.10);
  gPad->SetBottomMargin(0.13);
  leg = new TLegend(0.65,0.15,0.85,0.15+0.04*(SM2-SM1+1));

  TH1* hAvNCellsInClusterEvents = ConvertMapRuns2Events(nruns, runNumbers, hAvNCellsInCluster);
  hAvNCellsInClusterEvents->Draw();
  leg->AddEntry(hAvNCellsInClusterEvents, "All Modules","l");
  for (Int_t sm = SM1; sm <= SM2; sm++) {
    TH1* hAvNCellsInClusterSMEvents = ConvertMapRuns2Events(nruns, runNumbers, hAvNCellsInClusterSM[sm]);
    hAvNCellsInClusterSMEvents->Draw("same");
    leg->AddEntry(hAvNCellsInClusterSMEvents, Form("Module %i",5-sm),"l");
  }
  hAvNCellsInClusterEvents->Draw("same"); // to top
  leg->Draw("same");

  c1->Update();
  c2->Update();
}

//_________________________________________________________________________
Double_t GetAcceptance(Int_t sm, TH2* hmap, Int_t ri)
{
  // Returns [#cells - #dead]/#cells for a supermodule.
  // hmap -- dead/noisy cell map;
  // ri -- run index.

  // guess number of cells per supermodule
  Int_t nSM = 1152; // EMCAL
  if (hmap->GetXaxis()->GetXmin() == 1) {// PHOS
    nSM = 3584;
    sm--; // starts from 1, convenient from 0
  }

  // count dead cells
  Int_t ndead = 0;
  for (Int_t k = 1; k <= nSM; k++)
    if (hmap->GetBinContent(nSM*sm + k, ri+1) < 0)
      ndead++;

  return ((Double_t) (nSM - ndead))/nSM;
}

//_________________________________________________________________________
void DrawPi0Slice(Int_t run, Int_t sm = -1)
{
  // Draw the pi0 peak;
  // run,sm -- run number and supermodule to take;
  // sm < 0 -- draw for whole detector.

  TH1* h = NULL;
  if (sm >= 0) {//particular supermodule
    h = (TH1*) gInputFile->Get(Form("run%i_hPi0MassSM%iSM%i",run,sm,sm));
    if (h == 0) {
      Error("DrawPi0Slice","Run %d does contain histogram %s",run,Form("run%i_hPi0MassSM%iSM%i",run,sm,sm));
      continue;
    }
    h->SetName(Form("hPi0SliceSM%i_run%i",sm,run));
    h->SetTitle(Form("#pi^{0} in M%i, run %i, %.0fk events", sm, run, GetNumberOfEvents(run)/1e+3));
  }
  else {// whole detector
    for (Int_t sm1 = 0; sm1 < 10; sm1++)
      for (Int_t sm2 = sm1; sm2 < 10; sm2++) {
        TH1* one = (TH1*) gInputFile->Get(Form("run%i_hPi0MassSM%iSM%i",run,sm1,sm2));
        if (!one) continue;

        if (!h) h = one;
        else {
          h->Add(one);
          delete one;
        }
      }
    h->SetName(Form("hPi0Slice_run%i",run));
    h->SetTitle(Form("#pi^{0} in all modules, run %i, %.0fk events", run, GetNumberOfEvents(run)/1e+3));
  }

  h->SetXTitle("M_{#gamma#gamma}, GeV");
  h->SetYTitle("Entries");
  h->SetTitleOffset(1.7,"Y");

  TCanvas *c1;
  if (sm >= 0) c1 = new TCanvas(Form("hPi0SliceSM%i_run%i",sm,run),Form("hPi0SliceSM%i_run%i",sm,run), 400,400);
  else         c1 = new TCanvas(Form("hPi0Slice_run%i",run),Form("hPi0Slice_run%i",run), 400,400);

  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.06);
  h->Draw();

  Double_t nraw, enraw, mass, emass, sigma, esigma;
  FitPi0(h, nraw, enraw, mass, emass, sigma, esigma);

  // draw background
  TF1* fitfun = h->GetFunction("fitfun");
  if (fitfun) {
    Double_t emin, emax;
    fitfun->GetRange(emin, emax);

    backgr = new TF1("mypol2", "[0] + [1]*(x-0.135) + [2]*(x-0.135)^2", emin, emax);
    backgr->SetLineColor(kBlue);
    backgr->SetLineWidth(2);
    backgr->SetLineStyle(3);
    backgr->SetParameters(fitfun->GetParameter(3), fitfun->GetParameter(4), fitfun->GetParameter(5));
    backgr->Draw("same");
  }

  c1->Update();
}

//_________________________________________________________________________
void DrawPi0Averages(Int_t nruns, Int_t runNumbers[], Bool_t samesm = kFALSE, TH2* hmap = NULL)
{
  // Draw average number of pi0s per event, pi0 mass position and pi0 width per run index.
  // Errors are also drawn.
  //
  // samesm -- take only pi0s for which gammas were is same SM;
  // hmap -- if not NULL, do simple (area law based) acceptance correction.
  //
  // TODO: PHOS needs pi0 between SMs rather than in one SM

  // suffix to names
  char *suff = TString(Form("%s%s", samesm ? "_sameSM":"", hmap ? "_corr4accep":"")).Data();

  // supermodule region
  Int_t SM1 = 0;
  Int_t SM2 = 10;
  while (SM1 <= SM2 && !gInputFile->Get(Form("run%i_hPi0MassSM%iSM%i",runNumbers[0],SM1,SM1)) ) SM1++;
  while (SM2 >= SM1 && !gInputFile->Get(Form("run%i_hPi0MassSM%iSM%i",runNumbers[0],SM2,SM2)) ) SM2--;

  // initialize histograms for the entire detector
  TH1* hPi0Num = new TH1F(Form("hPi0Num%s",suff),"Average number of #pi^{0}s per event", nruns,0,nruns);
  SetRunLabel(hPi0Num,nruns,runNumbers,1);
  // hPi0Num->SetXTitle("Run index");
  hPi0Num->SetYTitle("Number of #pi^{0}s");

  TH1* hPi0Mass = new TH1F(Form("hPi0Mass%s",suff),"#pi^{0} mass position", nruns,0,nruns);
  SetRunLabel(hPi0Mass,nruns,runNumbers,1);
  // hPi0Mass->SetXTitle("Run index");
  hPi0Mass->SetYTitle("M_{#pi^{0}}, GeV");

  TH1* hPi0Sigma = new TH1F(Form("hPi0Sigma%s",suff),"#pi^{0} width", nruns,0,nruns);
  SetRunLabel(hPi0Sigma,nruns,runNumbers,1);
  // hPi0Sigma->SetXTitle("Run index");
  hPi0Sigma->SetYTitle("#sigma_{#pi^{0}}, GeV");

  // initialize histograms per SM
  TH1* hPi0NumSM[10];
  TH1* hPi0MassSM[10];
  TH1* hPi0SigmaSM[10];

  for (Int_t sm = SM1; sm <= SM2; sm++) {
    hPi0NumSM[sm] = new TH1F(Form("hPi0NumSM%i%s",sm,suff),"", nruns,0,nruns);
    hPi0MassSM[sm] = new TH1F(Form("hPi0MassSM%i%s",sm,suff),"", nruns,0,nruns);
    hPi0SigmaSM[sm] = new TH1F(Form("hPi0SigmaSM%i%s",sm,suff),"", nruns,0,nruns);
  }

  // run index loop
  for (Int_t ri = 0; ri < nruns; ri++)
  {
    fprintf(stderr, "\rDrawPi0Averages(): analysing run index %i/%i ...  ", ri, nruns-1);

    Int_t nevents = GetNumberOfEvents(runNumbers[ri]);

    // per SM histos
    for (Int_t sm = SM1; sm <= SM2; sm++) {
      TH1* h = (TH1*) gInputFile->Get(Form("run%i_hPi0MassSM%iSM%i",runNumbers[ri],sm,sm));
      if (h == 0) {
	Error("DrawPi0Averages","Run %d does contain histogram %s",runNumbers[ri],Form("run%i_hPi0MassSM%iSM%i",runNumbers[ri],sm,sm));
	continue;
      }

      // supermodule acceptance
      Double_t accep = 1.;
      if (hmap) accep = GetAcceptance(sm, hmap, ri);
      if (accep == 0) continue; // missing SM

      Double_t nraw, enraw, mass, emass, sigma, esigma;
      FitPi0(h, nraw, enraw, mass, emass, sigma, esigma);

      hPi0NumSM[sm]->SetBinContent(ri+1, nraw/nevents/accep);
      hPi0MassSM[sm]->SetBinContent(ri+1, mass);
      hPi0SigmaSM[sm]->SetBinContent(ri+1, sigma);

      hPi0NumSM[sm]->SetBinError(ri+1, enraw/nevents/accep);
      hPi0MassSM[sm]->SetBinError(ri+1, emass);
      hPi0SigmaSM[sm]->SetBinError(ri+1, esigma);

      delete h;
    } // supermodule loop


    /* fill for the entire detector
     */
    TH1* hsum = (TH1*) gInputFile->Get(Form("run%i_hPi0MassSM%iSM%i",runNumbers[ri],SM1,SM1));
    if (hsum == 0) {
      Error("DrawPi0Averages","Run %d does contain histogram %s",runNumbers[ri],Form("run%i_hPi0MassSM%iSM%i",runNumbers[ri],SM1,SM1));
      continue;
    }
    hsum->SetName("hSumTMP");

    // for acceptance correction
    Int_t noSM = 0;

    for (Int_t sm1 = SM1; sm1 <= SM2; sm1++)
      for (Int_t sm2 = sm1; sm2 <= SM2; sm2++) {
        if (sm1 == SM1 && sm2 == SM2) continue;
        if (samesm && sm1 != sm2) continue;

        TH1* h = (TH1*) gInputFile->Get(Form("run%i_hPi0MassSM%iSM%i",runNumbers[ri],sm1,sm2));
	if (h == 0) {
	  Error("DrawPi0Averages","Run %d does contain histogram %s",runNumbers[ri],Form("run%i_hPi0MassSM%iSM%i",runNumbers[ri],SM1,SM1));
	  continue;
	}

        // correct for SM acceptance
        if (hmap) {
          Double_t accep = GetAcceptance(sm1, hmap, ri);
          if (accep > 0) h->Scale(1./accep);
          else noSM++;
        }

        hsum->Add(h);
        delete h;
      }

    Double_t nraw, enraw, mass, emass, sigma, esigma;
    FitPi0(hsum, nraw, enraw, mass, emass, sigma, esigma);

    hPi0Num->SetBinContent(ri+1, nraw/nevents/(SM2-SM1+1-noSM));
    hPi0Mass->SetBinContent(ri+1, mass);
    hPi0Sigma->SetBinContent(ri+1, sigma);

    hPi0Num->SetBinError(ri+1, enraw/nevents/(SM2-SM1+1-noSM));
    hPi0Mass->SetBinError(ri+1, emass);
    hPi0Sigma->SetBinError(ri+1, esigma);

    delete hsum;
  } // run index loop

  fprintf(stderr, "\n");


  /* Draw results
   */

  // number of pi0s vs run index
  TCanvas *c1 = new TCanvas(Form("hPi0NumRuns%s",suff),Form("hPi0NumRuns%s",suff), 800,400);
  gPad->SetLeftMargin(0.08);
  gPad->SetRightMargin(0.03);
  gPad->SetTopMargin(0.12);
  TLegend *leg = new TLegend(0.75,0.15,0.89,0.15+0.06*(SM2-SM1+1));

  hPi0Num->SetAxisRange(0., hPi0Num->GetMaximum()*1.2, "Y");
  hPi0Num->SetTitleOffset(0.8, "Y");
  hPi0Num->SetLineWidth(2);
  hPi0Num->Draw();
  leg->AddEntry(hPi0Num, Form("(All Modules)/%i",SM2-SM1+1),"l");
  for (Int_t sm = SM1; sm <= SM2; sm++) {
    hPi0NumSM[sm]->SetLineColor(colorsSM[sm]);
    hPi0NumSM[sm]->SetLineWidth(2);
    hPi0NumSM[sm]->Draw("same");
    leg->AddEntry(hPi0NumSM[sm], Form("Module %i",5-sm),"l");
  }
  hPi0Num->Draw("same"); // to the top
  hPi0Num->Draw("hist,same");
  leg->Draw("same");


  // number of pi0s vs event count
  TCanvas *c2 = new TCanvas(Form("hPi0NumEvents%s",suff),Form("hPi0NumEvents%s",suff), 800,400);
  gPad->SetLeftMargin(0.08);
  gPad->SetRightMargin(0.06);
  gPad->SetTopMargin(0.12);
  leg = new TLegend(0.75,0.15,0.92,0.15+0.04*(SM2-SM1+1));

  TH1* hPi0NumEvents = ConvertMapRuns2Events(nruns, runNumbers, hPi0Num);
  hPi0NumEvents->Draw();
  leg->AddEntry(hPi0NumEvents, Form("(All Modules)/%i",SM2-SM1+1),"l");
  for (Int_t sm = SM1; sm <= SM2; sm++) {
    TH1* hPi0NumSMEvents = ConvertMapRuns2Events(nruns, runNumbers, hPi0NumSM[sm]);
    hPi0NumSMEvents->Draw("same");
    leg->AddEntry(hPi0NumSMEvents, Form("Module %i",5-sm),"l");
  }
  hPi0NumEvents->Draw("same"); // to the top
  hPi0NumEvents->Draw("hist,same");
  leg->Draw("same");


  // pi0 mass and width vs run index
  TCanvas *c3 = new TCanvas(Form("hPi0MassSigmaRuns%s",suff),Form("hPi0MassSigmaRuns%s",suff), 1000,707);
  c3->Divide(1,2);

  c3->cd(1);
  gPad->SetLeftMargin(0.06);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.10);
  gPad->SetBottomMargin(0.13);
  leg = new TLegend(0.75,0.15,0.95,0.15+0.06*(SM2-SM1+1));

  hPi0Mass->SetAxisRange(0.125, 0.145, "Y");
  hPi0Mass->SetTitleOffset(0.6, "Y");
  hPi0Mass->SetLineWidth(2);
  hPi0Mass->Draw();
  leg->AddEntry(hPi0Mass, "All Modules","l");
  for (Int_t sm = SM1; sm <= SM2; sm++) {
    hPi0MassSM[sm]->SetLineColor(colorsSM[sm]);
    hPi0MassSM[sm]->SetLineWidth(2);
    hPi0MassSM[sm]->Draw("same");
    leg->AddEntry(hPi0MassSM[sm], Form("Module %i",5-sm),"l");
  }
  hPi0Mass->Draw("same"); // to the top
  hPi0Mass->Draw("hist,same");
  leg->Draw("same");

  c3->cd(2);
  gPad->SetLeftMargin(0.06);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.10);
  gPad->SetBottomMargin(0.13);
  leg = new TLegend(0.75,0.15,0.95,0.15+0.06*(SM2-SM1+1));

  hPi0Sigma->SetAxisRange(0., hPi0Sigma->GetMaximum()*1.5, "Y");
  hPi0Sigma->SetTitleOffset(0.8, "Y");
  hPi0Sigma->SetLineWidth(2);
  hPi0Sigma->Draw();
  leg->AddEntry(hPi0Sigma, "All Modules","l");
  for (Int_t sm = SM1; sm <= SM2; sm++) {
    hPi0SigmaSM[sm]->SetLineColor(colorsSM[sm]);
    hPi0SigmaSM[sm]->SetLineWidth(2);
    hPi0SigmaSM[sm]->Draw("same");
    leg->AddEntry(hPi0SigmaSM[sm], Form("Module %i",5-sm),"l");
  }
  hPi0Sigma->Draw("same"); // to the top
  hPi0Sigma->Draw("hist,same");
  leg->Draw("same");


  // pi0 mass and width vs number of events
  TCanvas *c4 = new TCanvas(Form("hPi0MassSigmaEvents%s",suff),Form("hPi0MassSigmaEvents%s",suff), 800,400);
  c4->Divide(1,2);

  c4->cd(1);
  gPad->SetLeftMargin(0.16);
  gPad->SetRightMargin(0.08);
  leg = new TLegend(0.75,0.15,0.91,0.15+0.04*(SM2-SM1+1));

  TH1* hPi0MassEvents = ConvertMapRuns2Events(nruns, runNumbers, hPi0Mass);
  hPi0MassEvents->Draw();
  leg->AddEntry(hPi0MassEvents, "All Modules","l");
  for (Int_t sm = SM1; sm <= SM2; sm++) {
    TH1* hPi0MassSMEvents = ConvertMapRuns2Events(nruns, runNumbers, hPi0MassSM[sm]);
    hPi0MassSMEvents->Draw("same");
    leg->AddEntry(hPi0MassSMEvents, Form("Module %i",5-sm),"l");
  }
  hPi0MassEvents->Draw("same"); // to the top
  hPi0MassEvents->Draw("hist,same");
  leg->Draw("same");

  c4->cd(2);
  gPad->SetLeftMargin(0.16);
  gPad->SetRightMargin(0.08);
  leg = new TLegend(0.75,0.15,0.91,0.15+0.04*(SM2-SM1+1));

  TH1* hPi0SigmaEvents = ConvertMapRuns2Events(nruns, runNumbers, hPi0Sigma);
  hPi0SigmaEvents->Draw();
  leg->AddEntry(hPi0SigmaEvents, "All Modules","l");
  for (Int_t sm = SM1; sm <= SM2; sm++) {
    TH1* hPi0SigmaSMEvents = ConvertMapRuns2Events(nruns, runNumbers, hPi0SigmaSM[sm]);
    hPi0SigmaSMEvents->Draw("same");
    leg->AddEntry(hPi0SigmaSMEvents, Form("Module %i",5-sm),"l");
  }
  hPi0SigmaEvents->Draw("same"); // to the top
  hPi0SigmaEvents->Draw("hist,same");
  leg->Draw("same");


  c1->Update();
  c2->Update();
  c3->Update();
  c4->Update();
}

//_________________________________________________________________________
void DrawPi0Distributions(char *suff, Int_t nbins = 100)
{
  // Draw distributions for
  //    1) average number of pi0s per event;
  //    2) pi0 mass position;
  //    3) pi0 width.
  //
  // Must be called after DrawPi0Averages() because it
  // searches for the corresponding histograms by name.
  //
  // suff -- histograms suffix;
  // nbins -- number of bins in distributions.

  TCanvas *c1 = new TCanvas(Form("Pi0Distributions%s",suff),Form("Pi0Distributions%s",suff), 999,333);
  c1->Divide(3,1);

  // number of pi0s
  c1->cd(1)->SetLogy();
  gPad->SetLeftMargin(0.16);
  gPad->SetRightMargin(0.04);
  TH1* hPi0Num = (TH1*) gROOT->FindObject(Form("hPi0Num%s",suff));
  MakeDistribution(hPi0Num,nbins)->Draw();

  // pi0 mass
  c1->cd(2)->SetLogy();
  gPad->SetLeftMargin(0.16);
  gPad->SetRightMargin(0.04);
  TH1* hPi0Mass = (TH1*) gROOT->FindObject(Form("hPi0Mass%s",suff));
  MakeDistribution(hPi0Mass,nbins)->Draw();

  // pi0 width
  c1->cd(3)->SetLogy();
  gPad->SetLeftMargin(0.16);
  gPad->SetRightMargin(0.04);
  TH1* hPi0Sigma = (TH1*) gROOT->FindObject(Form("hPi0Sigma%s",suff));
  MakeDistribution(hPi0Sigma,nbins)->Draw();

  c1->Update();
}

//_________________________________________________________________________
TH1* MakeDistribution(TH1* inhisto, Int_t nbins)
{
  // Returns distribution for the input histogram.
  //
  // inhisto -- input histogram;
  // nbins -- number of bins in distribution.

  // initialize distribution; 1.01 -- to see the last bin
  TH1* distr = new TH1F(Form("%sDistr",inhisto->GetName()),"",
                        nbins, inhisto->GetMinimum(),inhisto->GetMaximum()*1.01);
  distr->SetXTitle(inhisto->GetYaxis()->GetTitle());
  distr->SetYTitle("Number of runs");

  // fill distribution
  for (Int_t i = 1; i <= inhisto->GetNbinsX(); i++)
    distr->Fill(inhisto->GetBinContent(i));

  return distr;
}

//_________________________________________________________________________
void FitPi0(TH1* h, Double_t &nraw, Double_t &enraw,
            Double_t &mass, Double_t &emass,
            Double_t &sigma, Double_t &esigma,
            Double_t emin = 0.05, Double_t emax = 0.3, Int_t rebin = 1)
{
  // Fits the pi0 peak with crystal ball + pol2,
  // fills number of pi0s, mass, width and their errors.

  nraw = enraw = 0;
  mass = emass = 0;
  sigma = esigma = 0;

  if (h->GetEntries() == 0) return NULL;

  if (rebin > 1) h->Rebin(rebin);

  // crystal ball parameters (fixed by hand for EMCAL); TODO: PHOS parameters?
  Double_t alpha = 1.1;  // alpha >= 0
  Double_t n = 2.;       // n > 1

  // CB tail parameters
  Double_t a = pow((n/alpha), n) * TMath::Exp(-alpha*alpha/2.);
  Double_t b = n/alpha - alpha;

  // integral value under crystal ball with amplitude = 1, sigma = 1
  // (will be sqrt(2pi) at alpha = infinity)
  Double_t nraw11 = a * pow(b+alpha, 1.-n)/(n-1.) + TMath::Sqrt(TMath::Pi()/2.) * TMath::Erfc(-alpha/TMath::Sqrt(2.));

  // signal (crystal ball);
  new TF1("cball", Form("(x-[1])/[2] > -%f ?                        \
                           [0]*exp(-(x-[1])*(x-[1])/(2*[2]*[2]))    \
                         : [0]*%f*(%f-(x-[1])/[2])^(-%f)", alpha, a, b, n));

  // background
  new TF1("mypol2", "[0] + [1]*(x-0.135) + [2]*(x-0.135)^2", emin, emax);

  // signal + background
  TF1 *fitfun = new TF1("fitfun", "cball + mypol2", emin, emax);
  fitfun->SetParNames("A", "M", "#sigma", "a_{0}", "a_{1}", "a_{2}");
  fitfun->SetLineColor(kRed);
  fitfun->SetLineWidth(2);

  // make a preliminary fit to estimate parameters
  TF1* ff = new TF1("fastfit", "gaus(0) + [3]");
  ff->SetParLimits(0, 5., h->GetMaximum()*1.5);
  ff->SetParLimits(1, 0.1, 0.2);
  ff->SetParLimits(2, 0.005,0.030);
  ff->SetParameters(h->GetMaximum()/3., 0.135, 0.010, 0.);
  h->Fit(ff, "0QEM", "", 0.105, 0.165);

  fitfun->SetParameters(ff->GetParameter(0), ff->GetParameter(1), ff->GetParameter(2), ff->GetParameter(3));
  h->Fit(fitfun,"QLEMR", "");

  // calculate pi0 values
  mass = fitfun->GetParameter(1);
  emass = fitfun->GetParError(1);

  sigma = fitfun->GetParameter(2);
  esigma = fitfun->GetParError(2);

  Double_t A = fitfun->GetParameter(0);
  Double_t eA = fitfun->GetParError(0);

  nraw = nraw11 * A * sigma / h->GetBinWidth(1);
  enraw = nraw * (eA/A + esigma/sigma);
}


/* COMMON FUNCTIONS FOR RUN NUMBERS EXTRACTION, VISUALIZATION AND FILTERING.
 *
 * Input: histogram hNEventsProcessedPerRun which is taken from gInputFile.
 */

//_________________________________________________________________________
void GetRunNumbers(Int_t &nruns, Int_t runNumbers[])
{
  // Fill runNumbers array with run numbers according to hNEventsProcessedPerRun
  // histogram: actually analysed runs have positive bin content.

  TH1* hNEventsProcessedPerRun = (TH1*) gInputFile->Get("hNEventsProcessedPerRun");

  // check underflow/overflow
  if (hNEventsProcessedPerRun->GetBinContent(0) > 0.5)
    Warning("GetRunNumbers", "some run numbers are in underflow; they will be absent in the list of runs");
  if (hNEventsProcessedPerRun->GetBinContent(hNEventsProcessedPerRun->GetNbinsX()+1) > 0.5)
    Warning("GetRunNumbers", "some run numbers are in overflow; they will be absent in the list of runs");

  nruns = 0;

  for (Int_t i = 1; i <= hNEventsProcessedPerRun->GetNbinsX(); i++)
    if (hNEventsProcessedPerRun->GetBinContent(i) > 0.5) {
      runNumbers[nruns] = hNEventsProcessedPerRun->GetBinLowEdge(i);
      nruns++;
    }
}

//_________________________________________________________________________
Int_t GetNumberOfEvents(Int_t run)
{
  // Return number of events in run;
  // run -- run number.

  TH1* hNEventsProcessedPerRun = (TH1*) gInputFile->Get("hNEventsProcessedPerRun");

  // round the number of events to avoid precision surprizes
  return TMath::Nint( hNEventsProcessedPerRun->GetBinContent(hNEventsProcessedPerRun->FindBin(run)) );
}

//_________________________________________________________________________
Long64_t GetTotalNumberOfEvents(Int_t nruns, Int_t runNumbers[])
{
  // Return total number of events in all runs

  Long64_t ntotal = 0;

  for (Int_t i = 0; i < nruns; i++)
    ntotal += GetNumberOfEvents(runNumbers[i]);

  return ntotal;
}

//_________________________________________________________________________
void DrawRunsDistribution(Int_t nruns, Int_t runNumbers[], Int_t dnbins = 100)
{
  // Draw events and events distribution;
  // dnbins -- number of bins in distribution.

  // initialize run index histogram ...
  TH1* hNEventsPerRunIndex = new TH1F("hNEventsPerRunIndex", "Number of processed events per run index", nruns,0,nruns);
  SetRunLabel(hNEventsPerRunIndex,nruns,runNumbers,1);
  // hNEventsPerRunIndex->SetXTitle("Run index");
  hNEventsPerRunIndex->SetYTitle("Number of events");

  // ... and fill it
  for (Int_t i = 0; i < nruns; i++)
    hNEventsPerRunIndex->SetBinContent(i+1, GetNumberOfEvents(runNumbers[i]));

  // initialize distribution; 1.01 - to see the last bin
  TH1* distrib = new TH1F("hNEventsPerRunIndexDistr", "", dnbins, 0, hNEventsPerRunIndex->GetMaximum()*1.01);
  distrib->SetXTitle("Number of events");
  distrib->SetYTitle("Number of runs");

  // fill distribution
  for (Int_t i = 1; i <= nruns; i++)
    distrib->Fill(hNEventsPerRunIndex->GetBinContent(i));

  // draw histogram + distribution
  TCanvas *c1 = new TCanvas("hNEventsPerRunIndex","hNEventsPerRunIndex", 800,600);
  c1->Divide(1,2);

  c1->cd(1);
  gPad->SetLeftMargin(0.06);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.10);
  gPad->SetBottomMargin(0.13);
  hNEventsPerRunIndex->SetTitleOffset(0.6,"Y");
  hNEventsPerRunIndex->SetTickLength(0.01,"Y");
  hNEventsPerRunIndex->Draw();

  c1->cd(2);
  gPad->SetLeftMargin(0.06);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.10);
  gPad->SetBottomMargin(0.13);
  distrib->SetTitleOffset(0.6,"Y");
  distrib->SetTickLength(0.01,"Y");
  distrib->Draw();

  c1->Update();
}

//_________________________________________________________________________
void ExcludeSmallRuns(Int_t &nruns, Int_t runNumbers[], Int_t nEventsMin = 100e+3)
{
  // Exclude runs with number of events < nEventsMin

  printf("Excluding runs (< %.0fk events):", nEventsMin/1000.);

  for (Int_t i = 0; i < nruns; i++)
    if (GetNumberOfEvents(runNumbers[i]) < nEventsMin) {
      printf(" %i", runNumbers[i]);
      nruns--;
      runNumbers[i] = runNumbers[nruns];
      i--;
    }

  printf("\n");

  SortArray(nruns, runNumbers);
}

//_________________________________________________________________________
void ExcludeRunNumbers(Int_t &nruns, Int_t runNumbers[], Int_t nexcl, Int_t runs2Exclude[])
{
  // Exclude particular runs.
  //
  // nexcl -- number of runs in runs2Exclude[];
  // runs2Exclude -- array with run numbers to exclude.

  for (Int_t i = 0; i < nruns; i++)
    for (Int_t e = 0; e < nexcl; e++)
      if (runNumbers[i] == runs2Exclude[e]) {
        Printf("Excluding run: %i", runs2Exclude[e]);
        nruns--;
        runNumbers[i] = runNumbers[nruns];
        i--;
        break;
      }

  SortArray(nruns, runNumbers);
}

//_________________________________________________________________________
void SortArray(const Int_t nruns, Int_t runNumbers[])
{
  // Sort an array;
  // used for runNumbers array sorting.

  Int_t indices[nruns];
  Int_t runNumbers_unsort[nruns];

  for (Int_t i = 0; i < nruns; i++)
    runNumbers_unsort[i] = runNumbers[i];

  TMath::Sort(nruns, runNumbers_unsort, indices, kFALSE);

  for (Int_t i = 0; i < nruns; i++)
    runNumbers[i] = runNumbers_unsort[indices[i]];
}

//_________________________________________________________________________
void PrintRunNumbers(Int_t nruns, Int_t runNumbers[])
{
  // Print a table run index / run number / number of events.

  Printf("");
  Printf("| index | run number |  nevents  |");
  Printf("|-------|------------|-----------|");

  for (Int_t i = 0; i < nruns; i++)
    Printf("|  %-4i |  %-9i |  %-8i |", i, runNumbers[i], GetNumberOfEvents(runNumbers[i]));

  Printf("| Events in total:    %-10lli |\n", GetTotalNumberOfEvents(nruns, runNumbers));
}

//------------------------------------------------------------------------
void SetRunLabel(TH1 *histo, Int_t nruns, Int_t *runNumbers, Int_t axis)
{
  TString runText;
  for (Int_t i=0; i<nruns; i++) {
    runText = Form("%d",runNumbers[i]);
    if      (axis == 1)
      histo->GetXaxis()->SetBinLabel(i+1,runText.Data());
    else if (axis == 2)
      histo->GetYaxis()->SetBinLabel(i+1,runText.Data());
  }
  histo->LabelsOption("v");
}
