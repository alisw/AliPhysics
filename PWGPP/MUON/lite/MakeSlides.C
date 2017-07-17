#if !defined(__CINT__) || defined(__MAKECINT__)

#include <map>
#include <vector>

#include <Riostream.h>

// ROOT includes
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TEnv.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TArrayI.h"
#include "TPRegexp.h"
#include "TRegexp.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH2.h"
#include "TLegend.h"
#include "TTree.h"
#endif

Int_t fNpages = 0;
const char* outPdf = "trending.pdf";

//_________________________________
TString GetTriggerShort ( TString trigName )
{
  TObjArray* arr = trigName.Tokenize("-");
  TString shortName = arr->At(0)->GetName();
  if ( arr->GetEntries() > 2 ) {
    shortName.Append(Form("-%s",arr->At(1)->GetName()));
  }
  delete arr;
  return shortName;
}

//_________________________________
TObjArray* GetRunList ( TString filename )
{
  TFile* file = TFile::Open(filename.Data());
  if ( ! file ) return NULL;

  TCanvas* can = static_cast<TCanvas*>(file->Get("Global/AllTriggers"));
  if ( ! can ) return NULL;

  TH2* histo = NULL;
  TIter next(can->GetListOfPrimitives());
  TObject* obj = NULL;

  while ( (obj = next()) ) {
    if ( obj->InheritsFrom(TH2::Class()) ) {
      histo = static_cast<TH2*>(obj);
      break;
    }
  }

  if ( ! histo ) return NULL;

  TAxis* axis = histo->GetXaxis();
  if ( axis->GetNbins() == 0 ) return NULL;

  TObjArray* runList = new TObjArray();
  runList->SetOwner();

  for ( Int_t ibin=1; ibin<=axis->GetNbins(); ibin++ ) {
    runList->Add(new TObjString(axis->GetBinLabel(ibin)));
  }

  delete file;

  return runList;
}

//_________________________________
Bool_t PrintToPdf ( TString objName, TString filename, Bool_t closeFile, Bool_t warnIfMissing  )
{
  Bool_t isOk = kFALSE;
  TFile* file = TFile::Open(filename.Data());
  if ( ! file ) return isOk;

  TObject* obj = file->FindObjectAny(objName.Data());
  if ( obj ) {
    TCanvas* can = 0x0;
    if ( obj->IsA() == TCanvas::Class() ) {
      can = static_cast<TCanvas*>(obj);
      can->Draw(); // This is needed or the print will not work
    }
    else {
      can = new TCanvas(Form("%s_can",obj->GetName()),obj->GetTitle(),600,600);
      TString drawOpt = "e";
      if ( obj->InheritsFrom(TH2::Class()) ) {
        drawOpt = "COLZ";
        TH2* histo = static_cast<TH2*>(obj);
        if ( histo->GetMinimum() == 0. && histo->GetMaximum() == 0. ) histo->SetMaximum(0.1);
      }
      obj->Draw(drawOpt.Data());
    }
    TString outFile(outPdf);
    if ( fNpages == 0 ) outFile.Append("(");
    else if ( closeFile ) outFile.Append(")");
    fNpages++;
    can->Print(outFile.Data());
    if ( obj != can ) delete obj;
    delete can; // We do not want the canvas to be drawn
    isOk = kTRUE;
  }
  else if ( warnIfMissing ) {
    printf("Warning: cannot find %s\n",objName.Data());
  }

  delete file;

  return isOk;
}

//_________________________________
void EscapeSpecialChars ( TString& str )
{
  str.ReplaceAll("\\","");
  TString specials = "_";
  TObjArray* specialList = specials.Tokenize(" ");
  for ( Int_t ichar=0; ichar<specialList->GetEntries(); ichar++ ) {
    TString currChar = static_cast<TObjString*>(specialList->At(ichar))->GetString();
    if ( str.Contains(currChar.Data()) ) str.ReplaceAll(currChar.Data(),Form("\\\%s",currChar.Data()));
  }
  delete specialList;
}

//_________________________________
void BeginFrame ( TString title, ofstream& outFile, TString label = "" )
{
  if ( ! outFile.is_open() ) return;
  outFile << endl;
  outFile << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  outFile << "\\begin{frame}";
  if ( ! label.IsNull() ) outFile << "\\label{" << label.Data() << "}";
  outFile << endl;
  outFile << " \\frametitle{" << title.Data() << "}" << endl;
}

//_________________________________
void EndFrame ( ofstream& outFile )
{
  if ( ! outFile.is_open() ) return;
  outFile << "\\end{frame}" << endl;
}

//_________________________________
void MakeDefaultItem ( ofstream& outFile, TString defaultItem = "" )
{
  if ( ! outFile.is_open() ) return;
  outFile << "\\begin{itemize}" << endl;
  outFile << " \\item " << defaultItem.Data() << endl;
  outFile << "\\end{itemize}" << endl;
}

//_________________________________
Bool_t MakeSingleFigureSlide ( TString objName, TString filename, TString title, ofstream& outFile, TString label = "", Bool_t warnIfMissing = kTRUE, Bool_t closePdf = kFALSE )
{
  Bool_t isOk = PrintToPdf(objName,filename,closePdf,warnIfMissing);
  if ( ! isOk ) {
    return kFALSE;
  }

  if ( ! outFile.is_open() ) return kTRUE;

  BeginFrame(title,outFile,label);
  outFile << " \\begin{columns}[onlytextwidth]" << endl;
  outFile << "  \\column{\\textwidth}" << endl;
  outFile << "  \\centering" << endl;
  outFile << "  \\includegraphics[width=0.98\\textwidth,height=0.92\\textheight,page=" << fNpages << "]{" << outPdf << "}" <<endl;
  outFile << " \\end{columns}" << endl;
  EndFrame(outFile);
  return kTRUE;
}

//_________________________________
Bool_t MakeTriggerSlide ( TString filename, ofstream &outFile )
{
  TString baseName[3] = {"bendPlane","nonBendPlane","bothPlanes"};
  Int_t colors[3] = {kBlack,kRed,kBlue};
  TH2* histo[3] = {0x0,0x0,0x0};
  TH1* histoTotal = 0x0;
  TFile* file = TFile::Open(filename.Data());
  for ( Int_t ieff=0; ieff<3; ieff++ ) {
    histo[ieff] = static_cast<TH2*>(file->FindObjectAny(Form("effEvolution%sChamber",baseName[ieff].Data())));
    histo[ieff]->SetDirectory(0);
  }
  histoTotal = static_cast<TH1*>(file->FindObjectAny("totalEffEvolution"));
  histoTotal->SetDirectory(0);
  delete file;

  TCanvas* can = new TCanvas("trigChEff","trigChEff",600,600);
  can->Divide(2,2,0,0);
  TLegend* leg = 0x0;
  TString drawOpt = "";
  for ( Int_t ich=1; ich<=4; ich++ ) {
    can->cd(ich);
    if ( ich == 4 ) leg = new TLegend(0.3,0.17,0.7,0.37);
    gPad->SetTicks(1,1);
    if ( ich > 2 ) gPad->SetBottomMargin(0.15);
    for ( Int_t ieff=0; ieff<3; ieff++ ) {
      TH1* proj = histo[ieff]->ProjectionX(Form("%s%i",histo[ieff]->GetName(),10+ich),ich,ich);
      proj->SetLineColor(colors[ieff]);
      proj->SetMarkerColor(colors[ieff]);
      proj->SetMarkerStyle(20+ieff);
      proj->SetMarkerSize(0.7);
      proj->SetStats(kFALSE);
      proj->SetTitle(Form("Chamber %i",10+ich));
      proj->GetYaxis()->SetRangeUser(0.9,1.01);
      proj->GetXaxis()->SetLabelSize(0.06);
      proj->LabelsOption("v");
      drawOpt = "e";
      if ( gPad->GetListOfPrimitives() != 0 ) drawOpt.Append("same");
      proj->Draw(drawOpt.Data());
      if ( leg ) leg->AddEntry(proj,baseName[ieff].Data(),"lp");
    }
    if ( leg ) leg->Draw();
  }
  for ( Int_t ieff=0; ieff<3; ieff++ ) delete histo[ieff];
  can->Print(outPdf);
  fNpages++;
  delete can;
  can = new TCanvas("totalTrigEff","totalTrifEff",600,600);
  can->SetLeftMargin(0.15);
  can->SetBottomMargin(0.15);
  histoTotal->GetXaxis()->SetLabelSize(0.06);
  histoTotal->Draw("e");
  can->Print(outPdf);
  fNpages++;
  delete can;

  if ( ! outFile.is_open() ) return kTRUE;

  BeginFrame("Trigger chamber efficiencies", outFile);
  outFile << " \\begin{columns}[onlytextwidth]" << endl;
  outFile << "  \\column{0.66\\textwidth}" << endl;
  outFile << "  \\centering" << endl;
  outFile << "  \\includegraphics[width=0.98\\textwidth,page=" << fNpages-1 << "]{" << outPdf << "}" << endl;
  outFile << "  \\column{0.34\\textwidth}" << endl;
  outFile << "  \\centering" << endl;
  outFile << "  \\includegraphics[width=0.98\\textwidth,page=" << fNpages << "]{" << outPdf << "}" << endl;
  outFile << " \\end{columns}" << endl;
  EndFrame(outFile);
  return kTRUE;
}

//_________________________________
Bool_t MakeTriggerRPCslide ( TString filename, ofstream &outFile, Bool_t outliers = kFALSE )
{
  for ( Int_t ich=0; ich<4; ich++ ) {
    TString currName = Form("effEvolutionbothPlanesRPCCh%i",11+ich);
    if ( outliers ) currName += "_outlier";
    PrintToPdf(currName,filename,kFALSE,kFALSE);
  }

  if ( ! outFile.is_open() ) return kTRUE;
  TString baseName = outliers ? "eff.-<eff.> for outliers" : "efficiency";
  BeginFrame(Form("Trigger chamber %s per RPC",baseName.Data()),outFile,outliers?"":"rpcEff");
  outFile << " \\begin{columns}[onlytextwidth]" << endl;
  outFile << "  \\column{\\textwidth}" << endl;
  outFile << "  \\centering" << endl;
  for ( Int_t ich=0; ich<4; ich++ ) {
    if ( ich%2 == 0 ) outFile << endl;
    outFile << "  \\includegraphics[width=0.37\\textwidth,page=" << fNpages-3+ich << "]{" << outPdf << "}" << endl;
  }
  outFile << " \\end{columns}" << endl;
  EndFrame(outFile);
  return kTRUE;
}

//_________________________________
Bool_t MakeTriggerRPCslides ( TString filename, ofstream &outFile, Bool_t fromTrackerTrack = kFALSE )
{
  TString trType = ( fromTrackerTrack ) ? "Tracker" : "Trigger";
  for ( Int_t ich=0; ich<4; ich++ ) {
    Int_t chamber = 11+ich;
    TString trigEffName = Form("%s_fromRPCEffCh%i",trType.Data(),chamber);
    Bool_t isOk = MakeSingleFigureSlide(trigEffName.Data(), filename.Data(), Form("MTR RPC efficiency for chamber %i (from %s track)",chamber, trType.Data()), outFile);
    if ( ! isOk ) return kFALSE;
  }
  return kTRUE;
}

//_________________________________
void MakeSummary ( TString period, ofstream &outFile )
{
  if ( ! outFile.is_open() ) return;
  BeginFrame("Summary I",outFile);
  outFile << "General informations" << endl;
  outFile << "\\begin{itemize}" << endl;
  outFile << " \\item Runs selected for MUON on ALICE logbook:" << endl;
  outFile << " \\begin{itemize}" << endl;
  outFile << "  \\item Run Type: PHYSICS" << endl;
  outFile << "  \\item Duration: at least 10 min" << endl;
  outFile << "  \\item GDC mStream Recording: Yes" << endl;
  outFile << "  \\item Period: " << period.Data() << endl;
  outFile << "  \\item Detectors: At least [ MUON\\_TRG \\& MUON\\_TRK ] as Readout" << endl;
  outFile << "  \\item Quality: globally GOOD and NOT BAD for readout Detectors" << endl;
  outFile << "  \\item Beam Mode: STABLE" << endl;
  outFile << " \\end{itemize}" << endl;
  outFile << "\\end{itemize}" << endl;
  outFile << endl;
  outFile << " \\vspace{5mm}" << endl;
  outFile << endl;
  outFile << "\\begin{columns}[onlytextwidth]" << endl;
  outFile << " \\column{\\textwidth}" << endl;
  outFile << " \\centering" << endl;
  outFile << " \\begin{tabular}{|l|lll|}" << endl;
  outFile << "  \\hline" << endl;
  outFile << "  & Total runs & CMUL & CMSL \\\\" << endl;
  outFile << "  \\hline" << endl;
  outFile << "  ALICE logbook & xx & xx & xx\\\\" << endl;
  outFile << "  Good from QA & xx & xx & xx\\\\" << endl;
  outFile << "  \\hline" << endl;
  outFile << " \\end{tabular}" << endl;
  outFile << "\\end{columns}" << endl;
  EndFrame(outFile);

  BeginFrame("Summary II",outFile);
  outFile << endl;
  outFile << "General:" << endl;
  MakeDefaultItem(outFile);
  outFile << endl;
  outFile << "MTR efficiency:" << endl;
  MakeDefaultItem(outFile);
  outFile << endl;
  outFile << "MCH and MUON data quality:" << endl;
  MakeDefaultItem(outFile);
  EndFrame(outFile);
}

//_________________________________
std::map<Int_t,std::vector<Double_t>> GetRunInfo ( TString evsQA )
{
  std::map<Int_t,std::vector<Double_t>> map;
  if ( gSystem->AccessPathName(evsQA.Data()) == 0 ) {
    TFile* file = TFile::Open(evsQA);
    TTree* tree = static_cast<TTree*>(file->Get("trending"));
    if ( tree ) {
      Int_t run, fill, bcs;
      Double_t mu;
      tree->SetBranchAddress("run",&run);
      tree->SetBranchAddress("fill",&fill);
      tree->SetBranchAddress("bcs",&bcs);
      tree->SetBranchAddress("mu",&mu);
      for ( Long64_t ientry=0; ientry<tree->GetEntries(); ientry++ ) {
        tree->GetEntry(ientry);
        auto search = map.find(run);
        if ( search != map.end() ) continue;
        auto vec = &(map[run]);
        vec->push_back((Double_t)fill);
        vec->push_back((Double_t)bcs);
        vec->push_back(mu);
      }
    }
    delete file;
  }
  return map;
}

//_________________________________
void MakeRunSummary ( ofstream &outFile, TString trackerQA, TString evsQA = "", ifstream* inFile = 0x0 )
{

  std::map<Int_t,std::vector<Double_t>> map = GetRunInfo(evsQA);

  TObjArray* runListArr = GetRunList(trackerQA);
  runListArr->Sort();

  Bool_t readSummary = ( inFile ) ? kTRUE : kFALSE;

  TString romanNum[10] = {"I","II","III","IV","V","VI","VII","VIII","IX","X"};

  Int_t nRuns = runListArr->GetEntries();
  Int_t nRunsPerPage = 40;
  Int_t nRunsPerColumn = nRunsPerPage/2;

  Int_t nPages = nRuns/nRunsPerPage;
  if ( nRuns%nRunsPerPage > 0 ) nPages++;

  Int_t irun = 0;
  Int_t readRun = -2, currRun = -1;

  Int_t previousFill = -1;
  Double_t previousMu = 0.;

  for ( Int_t ipage=0; ipage<nPages; ipage++ ) {
    TString title = "Run summary";
    if ( nPages > 1 ) title += Form(" (%s)",romanNum[ipage].Data());
    BeginFrame(title,outFile);
    outFile << " \\begin{columns}[onlytextwidth,T]" << endl;
    outFile << "  \\footnotesize" << endl;
    for ( Int_t icol=0; icol<2; icol++ ) {
      Bool_t needsHline = ( icol == 0 || irun < nRuns );
      outFile << "  \\column{0.5\\textwidth}" << endl;
      outFile << "  \\centering" << endl;
      outFile << "  \\begin{tabular}{|cp{0.63\\textwidth}|}" << endl;
      if ( needsHline ) outFile << "   \\hline" << endl;
      if ( nRuns == 0 ) {
        outFile << "   \\runTab[\\errorColor]{xxx}{xxx}" << endl;
      }
      else {
        while ( irun<nRuns ) {
          currRun = static_cast<TObjString*>(runListArr->UncheckedAt(irun++))->GetString().Atoi();
          Bool_t isNew = kTRUE;
          TString readLines = "", currLine = "";
          while ( readSummary ) {
            currLine.ReadLine(*inFile,kFALSE);
            if ( currLine.Contains("runTab") ) {
              TString sRun = currLine(TRegexp("{[0-9][0-9][0-9][0-9][0-9][0-9]}"));
              sRun.Remove(TString::kLeading,'{');
              sRun.Remove(TString::kTrailing,'}');
              readRun = sRun.Atoi();
              if ( readRun <= currRun ) readLines += currLine + "\n";
              if ( readRun == currRun ) {
                isNew = kFALSE;
                break;
              }
            }
            else if ( currLine.Contains("colorLegend") ) readSummary = kFALSE;
          }
          if ( isNew ) {
            auto search = map.find(currRun);
            TString info = "";
            if ( search != map.end() ) {
              auto vec = search->second;
              Int_t fill = (Int_t)vec[0];
              Double_t mu = vec[2];
              Double_t ratio = previousMu > 0. ? mu/previousMu : 10.;
              if ( fill != previousFill || TMath::Abs(1.-ratio) > 0.5 ) {
                previousFill = fill;
                previousMu = mu;
                info = Form("Fill %i, IB %i, mu %.3f",fill,(Int_t)vec[1],mu);
              }
            }
            outFile << "   \\runTab{" << currRun << "}{" << info.Data() << "}" << endl;
          }
          else outFile << readLines.Data();
          if ( irun%nRunsPerColumn == 0 ) break;
        }
      }

      if ( needsHline ) outFile << "   \\hline" << endl;
      if ( icol == 1 && ipage == nPages -1 ) {
        outFile << "   \\hline" << endl;
        outFile << "   \\colorLegend" << endl;
        outFile << "   \\hline" << endl;
      }
      outFile << "  \\end{tabular}" << endl;
      if ( icol == 0 ) outFile << endl;
      else {
        outFile << " \\end{columns}" << endl;
        EndFrame(outFile);
      }
    } // loop on columns
  } // loop on pages

  delete runListArr;
}

//_________________________________
void MakePreamble ( ofstream &outFile )
{
  if ( ! outFile.is_open() ) return;
  outFile << "\\documentclass[9pt,table]{beamer}" << endl;
  outFile << "\\mode<presentation>" << endl;
  outFile << "\\usepackage[T1]{fontenc}" << endl;
  outFile << "\\usepackage{lmodern}" << endl;
  outFile << "\\usepackage{textcomp}" << endl;
  outFile << "\\usepackage{amsmath}" << endl;
  outFile << "\\usepackage{color,graphicx}" << endl;
  outFile << "\\usepackage{colortbl}" << endl;
  outFile << "\\usepackage{multirow}" << endl;
  outFile << "\\usepackage{pifont}" << endl;
  outFile << "\\usepackage{wasysym}" << endl;
  outFile << "\\usepackage{appendixnumberbeamer}" << endl;
  outFile << "\\usepackage[absolute,overlay]{textpos}" << endl;
  outFile << "\\usetheme{Madrid}" << endl;
  outFile << "\\useoutertheme{shadow}" << endl;
  outFile << endl;
  outFile << "\\setbeamersize{text margin left=0.5cm, text margin right=0.5cm}" << endl;
  outFile << endl;
  outFile << "\\hypersetup{colorlinks,linkcolor=red,urlcolor=blue}" << endl;
  outFile << endl;
  outFile << "% Slightly change the template" << endl;
  outFile << "\\setbeamertemplate{navigation symbols}{} %suppress navigation symbols (bottom-right of frame)" << endl;

  outFile << "\\setbeamercolor*{author in head/foot}{parent=palette tertiary}" << endl;
  outFile << "\\setbeamercolor*{title in head/foot}{parent=palette secondary}" << endl;
  outFile << "\\setbeamercolor*{date in head/foot}{parent=palette primary}" << endl;
  outFile << "\\setbeamercolor*{section in head/foot}{parent=palette tertiary}" << endl;
  outFile << "\\setbeamercolor*{subsection in head/foot}{parent=palette primary}" << endl;
  outFile << "\\newcommand{\\changeFootline}[1]{" << endl;
  outFile << " \\setbeamertemplate{footline}{" << endl;
  outFile << "  \\hbox{%" << endl;
  outFile << "   \\begin{beamercolorbox}[wd=.2\\paperwidth,ht=2.25ex,dp=1ex,center]{author in head/foot}%" << endl;
  outFile << "    \\insertshortauthor%~~(\\insertshortinstitute)" << endl;
  outFile << "   \\end{beamercolorbox}%" << endl;
  outFile << "   \\begin{beamercolorbox}[wd=.6\\paperwidth,ht=2.25ex,dp=1ex,center]{title in head/foot}%" << endl;
  outFile << "    \\hypersetup{hidelinks}%" << endl;
  outFile << "    \\insertshorttitle" << endl;
  outFile << "    \\hspace*{2em}\\insertshortdate{}" << endl;
  outFile << "   \\end{beamercolorbox}%" << endl;
  outFile << "   \\begin{beamercolorbox}[wd=.2\\paperwidth,ht=2.25ex,dp=1ex,right]{date in head/foot}%" << endl;
  outFile << "    #1\\hspace*{2ex}" << endl;
  outFile << "  \\end{beamercolorbox}}%" << endl;
  outFile << " }}" << endl;
  outFile << "\\changeFootline{\\insertframenumber / \\inserttotalframenumber}" << endl;
  outFile << "\\setbeamertemplate{headline}{}" << endl;
  outFile << endl;
  outFile << endl;
  outFile << "\\newcommand{\\badForPassColor}{magenta!50!white}" << endl;
  outFile << "\\newcommand{\\errorColor}{red!50!white}" << endl;
//  outFile << "\\newcommand{\\newColor}{blue!20!white}" << endl;
  outFile << "\\newcommand{\\notInLogColor}{black!20!white}" << endl;
  outFile << "\\newcommand{\\pendingColor}{yellow!50!white}" << endl;
  outFile << "\\newcommand{\\warningColor}{orange!50!white}" << endl;
  outFile << "\\newcommand{\\colorLegend}{" << endl;
//  outFile << "  \\multicolumn{2}{|l|}{\\colorbox{\\newColor}{~~} = newly analyzed}\\\\" << endl;
  outFile << "  \\multicolumn{2}{|l|}{\\colorbox{\\notInLogColor}{~~} = non-selected from e-logbook}\\\\" << endl;
  outFile << "  \\multicolumn{2}{|l|}{\\colorbox{\\pendingColor}{~~} = pending statement}\\\\" << endl;
  outFile << "  \\multicolumn{2}{|l|}{\\colorbox{\\warningColor}{~~} = possible problem}\\\\" << endl;
  outFile << "  \\multicolumn{2}{|l|}{\\colorbox{\\errorColor}{~~} = problem spotted}\\\\" << endl;
  outFile << "  \\multicolumn{2}{|l|}{\\colorbox{\\badForPassColor}{~~} = bad for this pass}\\\\}" << endl;
  outFile << "\\newcommand{\\runTab}[3][white]{\\cellcolor{#1} #2 & \\cellcolor{#1} #3\\\\}" << endl;
  outFile << endl;
  outFile << "\\newcommand{\\pik}{\\ensuremath{\\pi\\mathrm{/K}}}" << endl;
  outFile << "\\newcommand{\\mum}{\\mbox{$\\mu {\\rm m}$}}" << endl;
  outFile << "\\newcommand{\\mom}{\\mbox{GeV$\\kern-0.15em /\\kern-0.12em c$}}" << endl;
  outFile << "\\newcommand{\\pt}{\\ensuremath{p_{\\mathrm{T}}}}" << endl;
  outFile << "\\newcommand{\\dd}{\\text{d}}" << endl;
  outFile << "\\newcommand{\\raa}{\\ensuremath{R_{AA}}}" << endl;
  outFile << "\\newcommand{\\un}[1]{\\protect\\detokenize{#1}}" << endl;
  outFile << endl;
}

//_________________________________
void BeginSlides ( TString period, TString pass, TString authors, ofstream &outFile )
{
  if ( ! outFile.is_open() ) return;
  TString authorsShort = "";
  Bool_t previousIsLetter = kFALSE;
  for ( Int_t ichar=0; ichar<authors.Length(); ichar++ ) {
    TString currChar = authors[ichar];
    currChar.ToUpper();
    Int_t currentIsLetter = currChar.IsAlpha();
    if ( currentIsLetter && ! previousIsLetter ) authorsShort += currChar + ".";
    if ( currChar == "," ) authorsShort += ",";
    previousIsLetter = currentIsLetter;
  }

  outFile << "\\title{Muon QA: " << period.Data() << " " << pass.Data() << "}" << endl;
  outFile << "\\author[" << authorsShort.Data() << "]{" << authors.Data() << "}" << endl;
  outFile << "\\date{\\today}" << endl;

  outFile << "\\begin{document}" << endl;
  outFile << "\\setlength{\\TPHorizModule}{1bp}" << endl;
  outFile << "\\setlength{\\TPVertModule}{1bp}" << endl;
  outFile << "\\textblockorigin{0bp}{0bp}" << endl;
  outFile << endl;
  outFile << "\\graphicspath{{images/}}" << endl;

  outFile << endl;
  outFile << "\\begin{frame}" << endl;
  outFile << " \\titlepage" << endl;
  outFile << "\\end{frame}" << endl;
}

//_________________________________
void EndSlides ( ofstream &outFile )
{
  if ( ! outFile.is_open() ) return;
  outFile << "\\end{document}" << endl;
  outFile.close();
}

//_________________________________
void StartAppendix ( ofstream &outFile )
{
  if ( ! outFile.is_open() ) return;
  outFile << endl;
  outFile << endl;
  outFile << "\\AtBeginSection[] % Do nothing for \\section*" << endl;
  outFile << "{" << endl;
  outFile << "  \\begin{frame}<beamer>" << endl;
  outFile << "   \\begin{beamercolorbox}[sep=8pt,center,shadow=true,rounded=true]{title}" << endl;
  outFile << "   \\usebeamerfont{title}\\insertsectionhead" << endl;
  outFile << "   \\end{beamercolorbox}" << endl;
  outFile << "  \\end{frame}" << endl;
  outFile << "}" << endl;
  outFile << endl;
  outFile << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  outFile << "\\appendix" << endl;
  outFile << "\\renewcommand{\\theframenumber}{A.\\arabic{framenumber}}" << endl;
  outFile << "\\changeFootline{\\theframenumber}" << endl;
  outFile << "\\hypersetup{hidelinks}" << endl;
  outFile << "\\section{\\huge Backup slides}" << endl;
}

//_________________________________
void WriteRunList ( TString trackerQA, TString outFilename = "runListQA.txt" )
{
  TObjArray* runListArr = GetRunList(trackerQA);

  ofstream outFile(outFilename);
  for ( Int_t irun=0; irun<runListArr->GetEntries(); irun++ ) {
    outFile << static_cast<TObjString*>(runListArr->At(irun))->GetString().Atoi() << endl;
  }
  outFile.close();
  delete runListArr;
}


//_________________________________
TObjArray* UpdateExisting ( TString texFilename, TString trackerQA, TString evsQA )
{
  TObjArray* trigList = NULL;
  TString backupFile = texFilename;
  backupFile.Append(".backup");
  printf("Copying existing file into %s\n",backupFile.Data());
  TFile::Cp(texFilename.Data(),backupFile);
  ofstream outFile(texFilename.Data(),std::ofstream::out | std::ofstream::trunc);
  ifstream inFile(backupFile.Data());
  TString currLine = "", frameTitle = "";
  Int_t nBlanks = 0;
  while ( ! inFile.eof() ) {
    currLine.ReadLine(inFile,kFALSE);
    // Avoid too many blank lines
    if ( currLine.IsNull() ) {
      nBlanks++;
      if ( nBlanks > 2 ) continue;
    }
    else nBlanks = 0;
    if ( currLine == "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" ) {
      frameTitle.ReadLine(inFile,kFALSE);
      currLine += Form("\n%s",frameTitle.Data());
      frameTitle.ReadLine(inFile,kFALSE);
      currLine += Form("\n%s",frameTitle.Data());
      if ( frameTitle.Contains("frametitle") && frameTitle.Contains("Run summary") ) {
        MakeRunSummary(outFile, trackerQA, evsQA, &inFile);
        while ( currLine != "\\end{frame}" ) {
          currLine.ReadLine(inFile);
        }
        continue;
      }
    }
    else if ( currLine.Contains("%TriggerList=") ) {
      TString triggers = currLine;
      triggers.Remove(0,triggers.Index("=")+1);
      trigList = triggers.Tokenize(",");
    }
    outFile << currLine.Data() << endl;
  }
  inFile.close();

  return trigList;
}

//_________________________________
void MakeSlides ( TString period, TString pass, TString triggerList, TString authors, TString trackerQA = "QA_muon_tracker.root", TString triggerQA = "QA_muon_trigger.root", TString evsQA = "trending_evs.root", TString texFilename = "muonQA.tex", TString outRunList = "" )
{

  TString fileNames[2] = {trackerQA,triggerQA};
  for ( Int_t ifile=0; ifile<2; ifile++ ) {
    if ( gSystem->AccessPathName(fileNames[ifile].Data()) ) {
      printf("Fatal: cannot open %s\n",fileNames[ifile].Data());
      return;
    }
  }

  TObjArray* trigList = 0x0;
  ofstream outFile;
  if ( gSystem->AccessPathName(texFilename.Data()) == 0 ) {
    printf("Output file %s already exists: updating it\n", texFilename.Data());
    trigList = UpdateExisting(texFilename, trackerQA, evsQA);
  }
  else {
    EscapeSpecialChars(period);
    EscapeSpecialChars(pass);

    trigList = triggerList.Tokenize(",");

    outFile.open(texFilename);
    outFile << "%TriggerList=" << triggerList.Data() << endl;
    MakePreamble(outFile);
    BeginSlides(period,pass,authors,outFile);

    MakeSummary(period,outFile);
    MakeRunSummary(outFile, trackerQA, evsQA);
  }

  MakeSingleFigureSlide("AllTriggers",trackerQA,"Number of events per trigger",outFile);
  MakeSingleFigureSlide("L2AQAoverSCALERS",triggerQA,"Reconstruction: reconstructed triggers in QA wrt L2A from OCDB scalers",outFile);

    // Bool_t isNew = MakeSingleFigureSlide("Trigger_fromChamberEff", triggerQA, "MTR chamber efficiency (from trigger track)", outFile);
  Bool_t isNew = MakeTriggerRPCslides(triggerQA, outFile);
  if ( ! isNew ) MakeTriggerSlide(triggerQA,outFile);

  for ( Int_t itrig=0; itrig<trigList->GetEntries(); itrig++ ) {
    TString currTrig = trigList->At(itrig)->GetName();
    TString shortTrig = GetTriggerShort(currTrig);
    MakeSingleFigureSlide(Form("RatioTrackTypes_cent0trigger%s",currTrig.Data()),trackerQA,Form("Muon tracks / event in %s events",shortTrig.Data()),outFile);
    MakeSingleFigureSlide(Form("RatioTrackTypes_cent3trigger%s",currTrig.Data()),trackerQA,Form("Muon tracks / event in %s events (central collisions)",shortTrig.Data()),outFile,"",kFALSE);
    MakeSingleFigureSlide(Form("TrackMult_cent0trigger%s",currTrig.Data()),trackerQA,Form("Muon tracker-trigger tracks / event in %s events",shortTrig.Data()),outFile);
    MakeSingleFigureSlide(Form("AsymMatchedtrigger%s",currTrig.Data()),trackerQA,Form("Charge asymmetry in %s events",shortTrig.Data()),outFile);
    MakeSingleFigureSlide(Form("BeamGasMatchedtrigger%s",currTrig.Data()),trackerQA,Form("Rel. num. of beam-gas tracks (id. by p$\\times$DCA cuts) in %s events",shortTrig.Data()),outFile,currTrig);
  }
  MakeSingleFigureSlide("cNClusters",trackerQA,"Average number of clusters per track and dispersion",outFile);    MakeSingleFigureSlide("cNClustersPerCh",trackerQA,"Average number of clusters per chamber",outFile,"","clustersPerChamber");

  if ( outFile.is_open() ) StartAppendix(outFile);
  MakeSingleFigureSlide("PhysSelCutOnCollTrigger",trackerQA,"Physics selection effects",outFile);
  MakeSingleFigureSlide("cClusterHitMapPerCh",trackerQA,"Average cluster position per chamber",outFile,"","clustersPosition");
  MakeSingleFigureSlide("cChi2",trackerQA,"Tracking quality",outFile);

  if ( isNew ) MakeTriggerRPCslides(triggerQA, outFile, kTRUE);
  else {
    MakeTriggerRPCslide(triggerQA,outFile);
    MakeTriggerRPCslide(triggerQA,outFile,kTRUE);
  }
  MakeSingleFigureSlide("cLptHpt",trackerQA,"Trigger \\pt\\ cut",outFile,"",kTRUE,kTRUE);

  if ( outFile.is_open() ) {
    BeginFrame("Hardware issues",outFile);
    outFile << "MUON Trigger" << endl;
    MakeDefaultItem(outFile);
    outFile << endl;
    outFile << "MUON tracker" << endl;
    MakeDefaultItem(outFile);
    EndFrame(outFile);

    EndSlides(outFile);
  }

  delete trigList;

  if ( ! outRunList.IsNull() ) WriteRunList(trackerQA, outRunList);
}
