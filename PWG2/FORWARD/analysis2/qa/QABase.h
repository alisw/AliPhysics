/**
 * @file   QABase.h
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Nov 17 11:58:11 2011
 * 
 * @brief  Base class for QA active classes
 * 
 * @ingroup pwg2_forward_qa_scripts
 * 
 */
#ifndef QABASE_H
#define QABASE_H
#ifndef __CINT__
# include <TTree.h>
# include <TFile.h>
# include <TError.h>
# include <TCanvas.h>
# include <TSystem.h>
# include <fstream>
# include <TLatex.h>
# include <TStyle.h>
# include "QARing.h"
#else
class TTree;
class TFile;
class QARing;
class Global;
class TCanvas;
#endif

/**
 * Base class for active QA classes.  This manages the I/O files, like
 * the tree file, the LaTeX file, and the storage file 
 * 
 * @ingroup pwg2_forward_qa_scripts
 */
struct QABase 
{
  /** 
   * Constructor
   * 
   * @param single  If true, only process one file. 
   */
  QABase(Bool_t single=false) 
    : fFMD1i(0),
      fFMD2i(0),
      fFMD2o(0),
      fFMD3i(0),
      fFMD3o(0),
      fGlobal(0),
      fTree(0), 
      fOutput(0), 
      fStore(0),
      fTeX(0), 
      fTeXName(""),
      fToDelete(""),
      fCanvas(0), 
      fSingle(single)
  {}
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  QABase(const QABase& o)
    : fFMD1i(o.fFMD1i),
      fFMD2i(o.fFMD2i),
      fFMD2o(o.fFMD2o),
      fFMD3i(o.fFMD3i),
      fFMD3o(o.fFMD3o),
      fGlobal(o.fGlobal),
      fTree(o.fTree), 
      fOutput(o.fOutput),
      fStore(o.fStore),
      fTeX(o.fTeX), 
      fTeXName(o.fTeXName),
      fToDelete(o.fToDelete),
      fCanvas(o.fCanvas),
      fSingle(o.fSingle)
  {}
  /** 
   * Assignment operator
   * 
   * @return Reference to this object
   */
  QABase& operator=(const QABase&) { return *this; }
  /** 
   * Desctructor 
   */
  virtual ~QABase()
  {
    if (fFMD1i)  { delete fFMD1i; }
    if (fFMD2i)  { delete fFMD2i; }
    if (fFMD2o)  { delete fFMD2o; }
    if (fFMD3i)  { delete fFMD3i; }
    if (fFMD3o)  { delete fFMD3o; }
    if (fGlobal) { delete fGlobal; } 
    if (fTree)   { delete fTree; } 
    if (fOutput) { delete fOutput; }
    if (fStore)  { delete fStore; }
    if (fTeX)    { fTeX->close(); fTeX = 0; }
  }
  /** 
   * The name of the TTree output file
   * 
   * @return Output file name 
   */
  const char* OutputName() const { return "forward_trend.root"; }
  /** 
   * Initialize 
   * 
   * @param read If true initialise for reading tree file
   *
   * @return True on success
   */
  Bool_t Init(bool read=false)
  {
    fToDelete = "";
    if (!fSingle) { 
      fOutput = new TFile(OutputName(), (read ? "READ" : "RECREATE"));
      if (!fOutput) { 
	Error("Init", "Failed to open output file");
	return false;
      }
      if (read) fTree   = static_cast<TTree*>(fOutput->Get("T"));
      else      fTree   = new TTree("T", "T");
    }
    if (read) fGlobal = Global::SetBranch(fTree);
    else      fGlobal = Global::MakeBranch(fTree);
      
    fFMD1i->Init(fTree, read);
    fFMD2i->Init(fTree, read);
    fFMD2o->Init(fTree, read);
    fFMD3i->Init(fTree, read);
    fFMD3o->Init(fTree, read);

    return true;
  }
  /** 
   * Make a canvas, LaTeX file, and storage file 
   * 
   * @param title Title of canvas and files
   */
  void MakeCanvas(const char* title)
  {
    fTeX = new std::ofstream(Form("%s.tex", fTeXName.Data()));
    
    *fTeX << "\\documentclass[landscape,a4paper,12pt]{article}\n"
	  << "\\usepackage[margin=2cm,a4paper]{geometry}\n"
	  << "\\usepackage{graphicx}\n"
	  << "\\title{{\\Huge\\bf " << title << "}}\n"
	  << "\\author{{\\LARGE FMD Team}}\n"
	  << "\\date{{\\Large \\today}}\n"
	  << "\\begin{document}\n"
      	  << "\\thispagestyle{empty}\n"
	  << "\\maketitle" << std::endl;

    gStyle->SetPalette(1);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
    gStyle->SetTitleW(.4);
    gStyle->SetTitleH(.1);
    // gStyle->SetTitleColor(0);
    gStyle->SetTitleStyle(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleX(.6);

    fCanvas = new TCanvas("qa", title, 900, 700);
    fCanvas->SetFillColor(0);
    fCanvas->SetBorderSize(0);
    fCanvas->SetLeftMargin(0.15);
    fCanvas->SetRightMargin(0.02);
    fCanvas->SetTopMargin(0.10);
    fCanvas->SetBottomMargin(0.10);

    fStore = TFile::Open(Form("%s.root", fTeXName.Data()), "RECREATE");
    if (!fStore) 
      Warning("MakeCanvas", "Failed to make store %s.root", fTeXName.Data());
    
    fToDelete = "";
  }
  /** 
   * Put a title on the canvas.  Also clears page in LaTeX file
   * 
   * @param title Title
   */
  void CanvasTitle(const char* title)
  {
    if (!fCanvas) return;

    fCanvas->cd();
    fCanvas->Clear();
    fCanvas->SetBorderSize(0);
    fCanvas->SetLeftMargin(0.15);
    fCanvas->SetRightMargin(0.02);
    fCanvas->SetTopMargin(0.10);
    fCanvas->SetBottomMargin(0.10);

    *fTeX << "\\clearpage\n"
	  << "%% -----------------------------------------------------"
	  << std::endl;

    // Put title on top 
    TLatex* topText = new TLatex(.5, .99, title);
    topText->SetTextAlign(23);
    topText->SetTextSize(.038);
    topText->SetTextFont(42);
    topText->SetTextColor(kBlue+3);
    topText->SetNDC();
    topText->Draw();
  }
  /** 
   * Print the canvas to PNG and include it into LaTeX file
   * 
   * @param pngName Base name of PNG
   * @param runNo   Run number to append
   */
  void PrintCanvas(const char* pngName, UInt_t runNo)
  {
    TString s(Form("%s_%09d", pngName, runNo));
    PrintCanvas(s.Data());
  }
  /** 
   * Print the canvas to PNG and include it into LaTeX file.  Stores
   * canvas in storage file 
   * 
   * @param pngName Base name of PNG
   */
  void PrintCanvas(const char* pngName)
  {
    fCanvas->SaveAs(Form("%s.png", pngName));
    *fTeX << "\\begin{center}\n"
	  << "\\includegraphics[keepaspectratio,height=\\textheight]{"
	  << pngName << "}\n" 
	  << "\\end{center}" << std::endl;
    fToDelete.Append(Form(" %s.png", pngName));
    TDirectory* d = gDirectory;
    fStore->cd();
    fCanvas->Write();
    d->cd();
  }
  /** 
   * Close the LaTeX and storage files. Runs PDFLaTeX on LaTeX
   * file. Makes sure that the files are group writable.
   * 
   * @param deletePNGs  If true, delete intermident PNG files 
   */
  void Close(bool deletePNGs=true)
  {
    if (fTeX) {
      *fTeX << "\\end{document}\n"
	    << "%% EOF" << std::endl;
      fTeX->close();
      fTeX = 0;
      
      const char* base = fTeXName.Data();
      gSystem->Exec(Form("pdflatex %s.tex > /dev/null 2>&1", base));
      TString cmd(Form("rm -f %s.log %s.aux %s.tex %s", 
		       base, base, base, 
		       deletePNGs ? fToDelete.Data() : ""));
      gSystem->Exec(cmd.Data());
      gSystem->Exec(Form("chmod g+rw %s.pdf %s", base, 
			 deletePNGs ? "" : fToDelete.Data()));
      Info("Close", "PDF file %s.pdf generated", base);
    }
    if (fStore) {
      fStore->Write();
      fStore->Close();
      gSystem->Exec(Form("chmod g+rw %s.root", fTeXName.Data()));
    }
  }

  // --- Members -----------------------------------------------------
  QARing*        fFMD1i;	// Pointer to ring object
  QARing*        fFMD2i;	// Pointer to ring object
  QARing*        fFMD2o;	// Pointer to ring object
  QARing*        fFMD3i;	// Pointer to ring object
  QARing*        fFMD3o;	// Pointer to ring object
  Global*        fGlobal;	// Pointer to global run object
  TTree*         fTree;		// Pointer to tree object
  TFile*         fOutput;	// Pointer to tree file 
  TFile*         fStore;	// Pointer to storage file
  std::ofstream* fTeX;		// pointer to LaTeX stream
  TString        fTeXName;	// Base name of LaTeX file 
  TString        fToDelete;	// List of files to possibly delete
  TCanvas*       fCanvas;	// Pointer to canvas object
  Bool_t         fSingle;	// Whether we're processing one run only
};

#endif


// Local Variables:
//  mode: C++
// End:
