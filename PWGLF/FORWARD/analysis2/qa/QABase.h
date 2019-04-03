/**
 * @file   QABase.h
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Nov 17 11:58:11 2011
 * 
 * @brief  Base class for QA active classes
 * 
 * @ingroup pwglf_forward_qa_scripts
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
 * @ingroup pwglf_forward_qa_scripts
 */
struct QABase 
{
  QABase(const TString& dataType, 
	 Int_t          prodYear, 
	 const TString& period, 
	 const TString& pass) 
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
      fHtml(0),
      fTeXName("index"),
      fToDelete(""),
      fCanvas(0), 
      fOutputName("trending.root"),
      fDataType(dataType), 
      fYear(prodYear), 
      fPeriod(period), 
      fPass(pass)
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
      fHtml(o.fHtml),
      fTeXName(o.fTeXName),
      fToDelete(o.fToDelete),
      fCanvas(o.fCanvas),
      fOutputName(o.fOutputName),
      fDataType(o.fDataType), 
      fYear(o.fYear), 
      fPeriod(o.fPeriod), 
      fPass(o.fPass)
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
    if (fHtml)   { fHtml->close(); fHtml = 0; }
  }
  /**
   * Set the output file name 
   *
   * @param name Name of output (tree) file 
   */
  void SetOutputName(const char* name) { fOutputName = name; }
  /** 
   * The name of the TTree output file
   * 
   * @return Output file name 
   */
  const char* OutputName() const { return fOutputName.Data(); }
  /** 
   * Make a tree 
   * 
   * @param read If true, read from file 
   * 
   * @return True on success 
   */
  virtual Bool_t MakeTree(bool read)
  {
    fOutput = new TFile(OutputName(), (read ? "READ" : "RECREATE"));
    if (!fOutput) { 
      Error("MakeTree", "Failed to open output file");
      return false;
    }
    if (read) fTree   = static_cast<TTree*>(fOutput->Get("T"));
    else      fTree   = new TTree("T", "T");
    if (!fTree) { 
      Error("MakeTree", "No tree defined!");
      return false;
    }

    return true;
  }
  /** 
   * Initialize 
   * 
   * @param read If true initialise for reading tree file
   *
   * @return True on success
   */
  Bool_t Init(bool read=false)
  {
    if (!MakeTree(read)) return false;

    if (read) fGlobal = Global::SetBranch(fTree);
    else      fGlobal = Global::MakeBranch(fTree);

    fToDelete = "";
      
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
    if (fCanvas) { 
      delete fCanvas;
      fCanvas = 0;
    }
    const char* base = fTeXName.Data();
    gSystem->Exec(Form("rm -f %s.tex %s.html %s.root %s.pdf", 
		       base, base, base, base));

    fTeX = new std::ofstream(Form("%s.tex", base));
    fHtml = new std::ofstream(Form("%s.html", base));
    
    TString texTitle(title);
    texTitle.ReplaceAll("_", "-");

    const char* can = "http://aliqafmd.web.cern.ch/aliqafmd/";
    *fTeX << "\\documentclass[landscape,a4paper,12pt]{article}\n"
	  << "\\nonstopmode\n"
	  << "\\usepackage[margin=2cm,a4paper]{geometry}\n"
	  << "\\usepackage{graphicx}\n"
	  << "\\title{{\\Huge\\bf " << texTitle << "}}\n"
	  << "\\author{{\\LARGE FMD Team}}\n"
	  << "\\date{{\\Large \\today}}\n"
	  << "\\begin{document}\n"
      	  << "\\thispagestyle{empty}\n"
	  << "\\maketitle" << std::endl;

    *fHtml << "<!DOCTYPE html>\n"
	   << "<html>\n"
	   << " <head>\n"
	   << "  <title>QA information - "  << title << "</title>\n"
	   << "  <link rel='stylesheet' href='" << can << "style.css'>\n" 
	   << "  <link rel='stylesheet' href='style.css'>\n" 
	   << "  <link rel='shortcut icon' href='" << can
	   << "fmd_favicon.png' type='image/x-png'>\n"
	   << "  <link rel='shortcut icon' href='fmd_favicon.png' "
	   << "type='image/x-png'>\n"
	   << "  <script>\n"
	   << "  function browseRootFile(file) {\n"
           << "   var o = document.location.origin;\n"
	   << "   var b = o;\n"
           << "   if (b.search(/aliqafmd.web.cern.ch/)) b += '/aliqafmd/';\n"
           << "   b += 'jsRoot/';\n"
	   << "   var p = o + document.location.pathname;\n"
	   << "   p=p.replace(/" << base << ".html/,'');\n"
	   << "   if (p[p.lenght-1] != '/') p += '/';\n"
	   << "   var u=encodeURIComponent(p+file)+'&menu=no';\n"
	   << "   window.open(b + '?url=' + u, '_blank',\n"
	   << "               'location=no,menubar=no,status=no,titlebar=no');\n"
	   << "  }\n"
	   << "  </script>\n"
	   << " </head>\n"
	   << "<body>\n" 
	   << " <h1>" << title << "</h1>\n"
	   << " <table>" 
	   << std::endl;

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

    TString tit(title);
    tit.ReplaceAll("#LT", "&lt;");
    tit.ReplaceAll("#GT", "&gt;");
    tit.ReplaceAll("#Delta", "&Delta;");
    tit.ReplaceAll("#xi", "&xi;");
    tit.ReplaceAll("#sigma", "&sigma;");
    tit.ReplaceAll("#chi", "&chi;");
    tit.ReplaceAll("#nu", "&nu;");
    tit.ReplaceAll("_{p}", "<sub>p</sub>");
    tit.ReplaceAll("_{z}", "<sub>z</sub>");
    tit.ReplaceAll("^{2}", "<sup>2</sup>");
    *fHtml << "<tr><td>" << tit << "</td>" << std::flush;

    PutCanvasTitle(title);
  }
  /** 
   * Put a title on the canvas
   * 
   * @param title Title 
   */
  void PutCanvasTitle(const char* title)
  {
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
   */
  void PrintCanvas(const char* pngName, UInt_t /*runNo*/)
  {
    // TString s(Form("%s_%09d", pngName, runNo));
    TString s(pngName);
    PrintCanvas(s.Data());
  }
  /** 
   * Print the canvas to PNG and include it into LaTeX file.  Stores
   * canvas in storage file 
   * 
   * @param pngName Base name of PNG
   * @param areas   Areas to print 
   */
  void PrintCanvas(const char* pngName, TCollection* areas=0)
  {
    gSystem->Exec(Form("rm -f %s.png %s.html", pngName, pngName));
    fCanvas->SaveAs(Form("%s.png", pngName));
    *fTeX << "\\begin{center}\n"
	  << "\\includegraphics[keepaspectratio,height=\\textheight]{"
	  << pngName << "}\n" 
	  << "\\end{center}" << std::endl;
    *fHtml << "<td><a href='" << pngName << ".html'>Plot</a></td></tr>" 
	   << std::endl;

    std::ofstream img(Form("%s.html", pngName));
    img << "<!DOCTYPE html>\n"
	<< "<html>\n"
	<< " <head>\n"
	<< "  <title>" << pngName << "</title>\n"
	<< "  <link rel='stylesheet' href='style.css'></link>\n"
	<< "  <link rel='shortcut icon' href='fmd_favicon.png' "
	<< "type='image/x-png'>\n"
	<< " </head>\n"
	<< " <body>\n"
	<< "  <h1>" << pngName << "</h1>\n"
	<< "  <div id=\"imap\">\n"
	<< "    <img src=\"" << pngName << ".png\">\n";
    if (areas) { 
      TIter next(areas);
      TObject* o = 0;
      while ((o = next())) 
	img << "       " << o->GetName() << "\n";
    }
    img << "  </div>\n" << std::endl;
    WriteImageFooter(img, pngName);
    img << " </body>\n" 
	<< "</html>" << std::endl;
    img.close();
    gSystem->Exec(Form("chmod g+rw %s.html", pngName));

    fToDelete.Append(Form(" %s.png", pngName));
    TDirectory* d = gDirectory;
    fStore->cd();
    fCanvas->Write();
    d->cd();
  }
  /** 
   * Write out image footer 
   * 
   * @param o Output stream
   */
  virtual void WriteImageFooter(std::ostream& o, const char* /*pngName*/)
  {
    TDatime now;
    o << "<div class='back'>\n"
      << "<a href='" << fTeXName << ".html'>Back</a>\n"
      << "</div>\n"
      << "<div class='change'>\n"
      << "  Last update: " << now.AsString() << "\n"
      << "</div>" << std::endl;
  }
  /** 
   * Close the LaTeX and storage files. Runs PDFLaTeX on LaTeX
   * file. Makes sure that the files are group writable.
   * 
   * @param deletePNGs  If true, delete intermident PNG files 
   */
  void Close(bool deletePNGs=true)
  {
    const char* base = fTeXName.Data();
    if (fTeX) {
      *fTeX << "\\end{document}\n"
	    << "%% EOF" << std::endl;
      fTeX->close();
      fTeX = 0;
      
      gSystem->Exec(Form("pdflatex %s.tex > /dev/null 2>&1", base));
      Info("Close", "PDF file %s.pdf generated", base);
    }
    if (fHtml) {
      TDatime now;
      *fHtml << "</table>" << std::endl;
      WriteLinks();
      WriteFooter();
      *fHtml << "</body></html>" << std::endl;
      fHtml->close();
      fHtml = 0;
      gSystem->Exec(Form("chmod g+rw %s.html", fTeXName.Data()));
    }
    if (fStore) {
      fStore->Write();
      fStore->Close();
      gSystem->Exec(Form("chmod 664 %s.root", fTeXName.Data()));
    }
    TString cmd(Form("rm -f %s.log %s.aux %s.tex %s", 
		     base, base, base, 
		     deletePNGs ? fToDelete.Data() : ""));
    gSystem->Exec(cmd.Data());
    gSystem->Exec(Form("chmod g+rw %s.pdf %s", base, 
		       deletePNGs ? "" : fToDelete.Data()));
  }
  /** 
   * Output links 
   * 
   */
  virtual void WriteLinks() 
  {
    const char* jsRoot = "http://cern.ch/aliqafmd/jsRoot/";
    *fHtml << "<h3>Collection of plots</h3>\n" 
	   << "<ul>\n"
	   << "  <li><a href='" << fTeXName << ".pdf'>PDF</a></li>\n"
	   << "  <li><a href='" << fTeXName << ".root'>ROOT</a>\n"
	   << "    <button onclick='browseRootFile(\"" 
	   << fTeXName << ".root\")'>browse online</button></li>\n"
	   << "</ul>" << std::endl;
    if (fPeriod.IsNull()) return;
    Bool_t isMC = (fDataType.EqualTo("sim", TString::kIgnoreCase) || 
		   fPass.BeginsWith("passMC", TString::kIgnoreCase));
    *fHtml << "<ul>\n"
	   << " <li><a href='https://alimonitor.cern.ch/" 
	   <<  (isMC ? "job_details.jsp" : "production/raw.jsp") 
	   << "?jt_field1=" << fPeriod << "'>Production(s)</a></li>\n"
	   << "</ul>" << std::endl;
  }
  /** 
   * Write full job footer 
   * 
   */
  virtual void WriteFooter() 
  {
    TDatime now;
    *fHtml << "<div class='back'>\n"
	   << "<a href='index.html'>Back</a>\n"
	   << "</div>\n"
	   << "<div class='change'>\n"
	   << "  Last update: " << now.AsString() << "\n"
	   << "</div>" << std::endl;
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
  std::ofstream* fHtml;		// pointer to HTML stream
  TString        fTeXName;	// Base name of LaTeX file 
  TString        fToDelete;	// List of files to possibly delete
  TCanvas*       fCanvas;	// Pointer to canvas object
  TString        fOutputName;   // Output tree file name 
  TString        fDataType;	// Data type       
  Int_t          fYear;         // Production year
  TString        fPeriod;       // Period identifier 
  TString        fPass;         // Pass identifier 
};

#endif


// Local Variables:
//  mode: C++
// End:
