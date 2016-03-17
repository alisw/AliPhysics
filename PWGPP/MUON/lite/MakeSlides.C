#if !defined(__CINT__) || defined(__MAKECINT__)

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
#endif

//_________________________________
TString PdfToTxt ( TString filename, Bool_t clear = kFALSE )
{
  TString convertedFilename = filename;
  convertedFilename.ReplaceAll(".pdf",".txt");

  if ( clear ) {
    if ( gSystem->AccessPathName(convertedFilename.Data()) == 0 ) {
      gSystem->Exec(Form("rm %s",convertedFilename.Data()));
    }
    return "";
  }

  if ( gSystem->AccessPathName(convertedFilename.Data()) != 0 ) {
    gSystem->Exec(Form("gs -dBATCH -dNOPAUSE -sDEVICE=txtwrite -sOutputFile=- %s | xargs -I %% > %s",filename.Data(),convertedFilename.Data()));
  }

  return convertedFilename;
}

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
void EscapeSpecialCharsForRegex ( TString& str )
{
  TString specials = "+ ( )";
  TObjArray* specialList = specials.Tokenize(" ");
  for ( Int_t ichar=0; ichar<specialList->GetEntries(); ichar++ ) {
    TString currChar = static_cast<TObjString*>(specialList->At(ichar))->GetString();
    if ( str.Contains(currChar.Data()) ) str.ReplaceAll(currChar.Data(),Form("\\%s",currChar.Data()));
  }
  delete specialList;
}


//_________________________________
Int_t GetPage ( TString pattern, TString filename, TString trigger = "", Bool_t warnIfMissing = kTRUE )
{
  TString convertedFilename = PdfToTxt(filename);

  ifstream inFile(convertedFilename.Data());
  if ( ! inFile.is_open() ) return -1;

  EscapeSpecialCharsForRegex(pattern);

  TObjArray* patternArr = pattern.Tokenize("&");
  if ( ! trigger.IsNull() ) {
    trigger.Prepend("(^|[ ]|/)");
    trigger.Append("([ ]|$)");
    patternArr->Add(new TObjString(trigger));
  }

  TString currLine = "", currToken = "";
  Int_t currPage = -1, foundPage = -1;
  while ( ! inFile.eof() ) {
    currToken.ReadToken(inFile);
    if ( currToken == "Page" || inFile.eof() ) {
      Bool_t isOk = kTRUE;
      for ( Int_t ipat=0; ipat<patternArr->GetEntries(); ipat++ ) {
        TString currPattern(static_cast<TObjString*>(patternArr->At(ipat))->GetString());
        TPRegexp re(currPattern.Data());
        if ( ! currLine.Contains(re) ) {
          isOk = kFALSE;
          break;
        }
      }
      if ( isOk ) {
        foundPage = currPage;
        break;
      }
      if ( ! inFile.eof() ) {
        currToken.ReadToken(inFile);
        currPage = currToken.Atoi();
        currLine = "";
      }
    }
    else currLine += currToken + " ";
  }
  inFile.close();
  delete patternArr;
  if ( foundPage < 0 && warnIfMissing ) printf("Warning: cannot find %s\n",pattern.Data());

  return foundPage;
}

//_________________________________
TString GetRunList ( TString filename )
{
  TString convertedFilename = PdfToTxt(filename);

  ifstream inFile(convertedFilename.Data());
  if ( ! inFile.is_open() ) return -1;

  TString runList = "", currToken = "";
  TString keyword = "RUN:";
  while ( ! inFile.eof() ) {
    currToken.ReadToken(inFile);
    if ( currToken.Contains(keyword.Data()) ) {
      currToken.ReplaceAll(keyword.Data(),"");
      if ( currToken.Contains(",") || currToken.IsDigit() ) {
        runList = currToken;
        break;
      }
    }
  }
  inFile.close();
  return runList;
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
  outFile << "\\end{frame}" << endl;
}

//_________________________________
void MakeDefaultItem ( ofstream& outFile, TString defaultItem = "" )
{
  outFile << "\\begin{itemize}" << endl;
  outFile << " \\item " << defaultItem.Data() << endl;
  outFile << "\\end{itemize}" << endl;
}

//_________________________________
Bool_t MakeSingleFigureSlide ( TString pattern, TString filename, TString title, ofstream &outFile, TString trigger = "", TString label = "", Bool_t warnIfMissing = kTRUE )
{
  Int_t pageNum = GetPage(pattern,filename,trigger,warnIfMissing);
  if ( pageNum<0 ) {
    return kFALSE;
  }

  BeginFrame(title,outFile,label);
  outFile << " \\begin{columns}[onlytextwidth]" << endl;
  outFile << "  \\column{\\textwidth}" << endl;
  outFile << "  \\centering" << endl;
  outFile << "  \\includegraphics[width=0.98\\textwidth,height=0.92\\textheight,page=" << pageNum << "]{" << gSystem->BaseName(filename.Data()) << "}" <<endl;
  outFile << " \\end{columns}" << endl;
  EndFrame(outFile);
  return kTRUE;
}

//_________________________________
Bool_t MakeTriggerSlide ( TString filename, ofstream &outFile )
{
  BeginFrame("Trigger chamber efficiencies", outFile);
  outFile << " \\begin{columns}[onlytextwidth]" << endl;
  outFile << "  \\column{0.66\\textwidth}" << endl;
  outFile << "  \\centering" << endl;
  for ( Int_t ich=0; ich<4; ich++ ) {
    if ( ich%2 == 0 ) outFile << endl;
    Int_t ipage = GetPage(Form("Trigger chamber efficiency vs run for chamber %i",11+ich),filename);
    outFile << "  \\includegraphics[width=0.48\\textwidth,page=" << ipage << "]{" << filename.Data() << "}" << endl;
  }
  outFile << "  \\column{0.34\\textwidth}" << endl;
  outFile << "  \\centering" << endl;
  Int_t ipage = GetPage("Multinomial probability",filename);
  outFile << "  \\includegraphics[width=0.98\\textwidth,page=" << ipage << "]{" << filename.Data() << "}" << endl;
  outFile << " \\end{columns}" << endl;
  EndFrame(outFile);
  return kTRUE;
}

//_________________________________
Bool_t MakeTriggerRPCslide ( TString filename, ofstream &outFile, Bool_t outliers = kFALSE )
{
  TString baseName = outliers ? "eff.-<eff.> for outliers" : "efficiency";
  BeginFrame(Form("Trigger chamber %s per RPC",baseName.Data()),outFile,outliers?"":"rpcEff");
  outFile << " \\begin{columns}[onlytextwidth]" << endl;
  outFile << "  \\column{\\textwidth}" << endl;
  outFile << "  \\centering" << endl;
  for ( Int_t ich=0; ich<4; ich++ ) {
    if ( ich%2 == 0 ) outFile << endl;
    Int_t ipage = GetPage(Form("Trigger chamber %s vs run for chamber %i&RPC",baseName.Data(),11+ich),filename);
    outFile << "  \\includegraphics[width=0.37\\textwidth,page=" << ipage << "]{" << filename.Data() << "}" << endl;
  }
  outFile << " \\end{columns}" << endl;
  EndFrame(outFile);
  return kTRUE;
}

//_________________________________
void MakeSummary ( TString period, ofstream &outFile, TString trackerQA )
{
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
  MakeDefaultItem(outFile,"More than xx\\% efficiency in trigger chambers, stable");
  outFile << endl;
  outFile << "MCH and MUON data quality:" << endl;
  MakeDefaultItem(outFile);
  EndFrame(outFile);

  TString runList = GetRunList(trackerQA);
  TObjArray* runListArr = runList.Tokenize(",");
  runListArr->Sort();

  TString romanNum[10] = {"I","II","III","IV","V","VI","VII","VIII","IX","X"};

  Int_t nRuns = runListArr->GetEntries();
  Int_t nRunsPerPage = 40;
  Int_t nRunsPerColumn = nRunsPerPage/2;

  Int_t nPages = nRuns/nRunsPerPage;
  if ( nRuns%nRunsPerPage > 0 ) nPages++;

  Int_t irun = 0;

  for ( Int_t ipage=0; ipage<nPages; ipage++ ) {
    TString title = "Run summary";
    if ( nPages > 1 ) title += Form(" (%s)",romanNum[ipage].Data());
    BeginFrame(title,outFile);
    outFile << " \\begin{columns}[onlytextwidth,T]" << endl;
    outFile << "  \\footnotesize" << endl;
    outFile << "  \\column{0.5\\textwidth}" << endl;
    outFile << "  \\centering" << endl;
    outFile << "  \\begin{tabular}{|cp{0.63\\textwidth}|}" << endl;
    outFile << "   \\hline" << endl;
    if ( nRuns == 0 ) {
      outFile << "   \\runTab[\\errorColor]{xxx}{xxx}" << endl;
    }
    else {
      while ( irun<nRuns ) {
        outFile << "   \\runTab{" << static_cast<TObjString*>(runListArr->At(irun++))->GetString().Atoi() << "}{}" << endl;
        if ( irun%nRunsPerColumn == 0 ) break;
      }
    }

    outFile << "   \\hline" << endl;
    outFile << "  \\end{tabular}" << endl;
    outFile << endl;
    outFile << "  \\column{0.5\\textwidth}" << endl;
    outFile << "  \\begin{tabular}{|cp{0.63\\textwidth}|}" << endl;
    Bool_t hasRuns = ( irun < nRuns );
    if ( hasRuns ) outFile << "   \\hline" << endl;
    while ( irun<nRuns ) {
      outFile << "   \\runTab{" << static_cast<TObjString*>(runListArr->At(irun++))->GetString().Atoi() << "}{}" << endl;
      if ( irun%nRunsPerColumn == 0 ) break;
    }
    if ( hasRuns ) outFile << "   \\hline" << endl;
    if ( ipage == nPages -1 ) {
      outFile << "   \\hline" << endl;
      outFile << "   \\colorLegend" << endl;
      outFile << "   \\hline" << endl;
    }
    outFile << "  \\end{tabular}" << endl;
    outFile << " \\end{columns}" << endl;
    EndFrame(outFile);
  }

  delete runListArr;
}

//_________________________________
void MakePreamble ( ofstream &outFile )
{
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
  outFile << "\\changeFootline{\\insertframenumber{} / \\inserttotalframenumber}" << endl;
  outFile << "\\setbeamertemplate{headline}{}" << endl;
  outFile << endl;
  outFile << endl;
  outFile << "\\newcommand{\\badForPassColor}{magenta!50!white}" << endl;
  outFile << "\\newcommand{\\errorColor}{red!50!white}" << endl;
  outFile << "\\newcommand{\\newColor}{blue!20!white}" << endl;
  outFile << "\\newcommand{\\notInLogColor}{black!20!white}" << endl;
  outFile << "\\newcommand{\\pendingColor}{yellow!50!white}" << endl;
  outFile << "\\newcommand{\\warningColor}{orange!50!white}" << endl;
  outFile << "\\newcommand{\\colorLegend}{" << endl;
  outFile << "  \\multicolumn{2}{|l|}{\\colorbox{\\newColor}{~~} = newly analyzed}\\\\" << endl;
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
  outFile << "\\end{document}" << endl;
  outFile.close();
}

//_________________________________
void StartAppendix ( ofstream &outFile )
{
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
  outFile << "\\changeFootline{A.\\insertframenumber{}}" << endl;
  outFile << "\\hypersetup{hidelinks}" << endl;
  outFile << "\\section{\\huge Backup slides}" << endl;
}

//_________________________________
void WriteRunList ( TString trackerQA, TString outFilename = "runListQA.txt" )
{
  TString runList = GetRunList(trackerQA);
  TObjArray* runListArr = runList.Tokenize(",");

  ofstream outFile(outFilename);
  for ( Int_t irun=0; irun<runListArr->GetEntries(); irun++ ) {
    outFile << static_cast<TObjString*>(runListArr->At(irun))->GetString().Atoi() << endl;
  }
  outFile.close();
  delete runListArr;

  PdfToTxt(trackerQA,kTRUE);
}

//_________________________________
void MakeSlides ( TString period, TString pass, TString triggerList, TString authors, TString trackerQA = "QA_muon_tracker.pdf", TString triggerQA = "QA_muon_trigger.pdf", TString  texFilename = "muonQA.tex", TString outRunList = "" )
{
  if ( gSystem->AccessPathName(texFilename.Data()) == 0 ) {
    printf("Output file %s already exists\nPlease remove it!\n", texFilename.Data());
    return;
  }

  TString hasGs = gSystem->GetFromPipe("which gs");
  if ( hasGs.IsNull() ) {
    printf("The macro selects the pdf page with gs, but the program was not found on this machine. Sorry, but slides cannot be automatically generated on this machine.\n");
    return;
  }

  EscapeSpecialChars(period);
  EscapeSpecialChars(pass);

  TObjArray* trigList = triggerList.Tokenize(",");

  ofstream outFile(texFilename);
  outFile << "%TriggerList=" << triggerList.Data() << endl;
  MakePreamble(outFile);
  BeginSlides(period,pass,authors,outFile);

  MakeSummary(period,outFile,trackerQA);

  MakeSingleFigureSlide("Selections: RUN",trackerQA,"Number of events per trigger",outFile);
  MakeSingleFigureSlide("L2A from QA",triggerQA,"Reconstruction: reconstructed triggers in QA wrt L2A from OCDB scalers",outFile);
  MakeTriggerSlide(triggerQA,outFile);

  for ( Int_t itrig=0; itrig<trigList->GetEntries(); itrig++ ) {
    TString currTrig = trigList->At(itrig)->GetName();
    TString shortTrig = GetTriggerShort(currTrig);
    MakeSingleFigureSlide("Number of Tracks",trackerQA,Form("Muon tracks / event in %s events",shortTrig.Data()),outFile,currTrig);
    MakeSingleFigureSlide("Number of Tracks&for high mult.",trackerQA,Form("Muon tracks / event in %s events (central collisions)",shortTrig.Data()),outFile,currTrig,"",kFALSE);
    MakeSingleFigureSlide("Sum of trigger tracks (matched + trigger-only) / # events in",trackerQA,Form("Muon tracker-trigger tracks / event in %s events",shortTrig.Data()),outFile,currTrig);
    MakeSingleFigureSlide("Matched tracks charge asymmetry for&with acc. cuts",trackerQA,Form("Charge asymmetry in %s events",shortTrig.Data()),outFile,currTrig);
    MakeSingleFigureSlide("Identified beam-gas tracks (pxDCA cuts) in matched tracks for",trackerQA,Form("Rel. num. of beam-gas tracks (id. by p$\\times$DCA cuts) in %s events",shortTrig.Data()),outFile,currTrig);
  }
  MakeSingleFigureSlide("averaged number of associated clusters or of the number of chamber hit per track",trackerQA,"Average number of clusters per track and dispersion",outFile);
  MakeSingleFigureSlide("averaged number of clusters in chamber i per track",trackerQA,"Average number of clusters per chamber",outFile,"","clustersPerChamber");

  StartAppendix(outFile);
  MakeSingleFigureSlide("Physics Selection Cut on selected triggers:",trackerQA,"Physics selection effects",outFile);
  MakeSingleFigureSlide("<X> of clusters - associated to a track - in chamber i",trackerQA,"Average cluster position per chamber",outFile);
  MakeSingleFigureSlide("averaged normalized",trackerQA,"Tracking quality",outFile);

  MakeTriggerRPCslide(triggerQA,outFile);
  MakeTriggerRPCslide(triggerQA,outFile,kTRUE);
  MakeSingleFigureSlide("Trigger Lpt cut per run",trackerQA,"Trigger \\pt\\ cut",outFile);

  BeginFrame("Hardware issues",outFile);
  outFile << "MUON Trigger" << endl;
  MakeDefaultItem(outFile);
  outFile << endl;
  outFile << "MUON tracker" << endl;
  MakeDefaultItem(outFile);
  EndFrame(outFile);

  EndSlides(outFile);

  delete trigList;

  if ( ! outRunList.IsNull() ) WriteRunList(trackerQA, outRunList);

  // Clean converted txt files
  TString filenames[2] = {trackerQA, triggerQA};
  for ( Int_t ifile=0; ifile<2; ifile++ ) {
    PdfToTxt(filenames[ifile],kTRUE);
  }
}
