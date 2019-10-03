#ifndef __CINT__
# include <TFile.h>
# include <TString.h>
# include <TError.h>
# include <TSystem.h>
# include <fstream>
#else
class TFile;
class TString;
#endif


/** 
 * Extract information on a production. 
 * 
 */
struct ExtractProd
{
  TString fTmp;
  TString fPath;
  TString fPass;
  TString fRuns;
  Bool_t  fDebug;
  /** 
   * Constructor 
   */
  ExtractProd() 
    : fTmp(""), fPath(""), fPass(""), fRuns(""), fDebug(false)
  {
  }
  /** 
   * Run the code.  Query MonALisa for information.  
   * 
   * @param name      Production identifier 
   * @param mc        If true, assume simulations 
   * @param minSize   Least size of runs to use 
   * 
   * @return true on success. 
   */  
  Bool_t Run(const TString& name, Bool_t mc=false, ULong_t minSize=1000000)
  {
    fTmp = "prod";
    gSystem->TempFileName(fTmp);
    gSystem->Unlink(fTmp);
    gSystem->mkdir(fTmp);
    if (fDebug) Info("", "TMP=%s", fTmp.Data());

    fPath = "";
    fPass = "";
    fRuns = "";
    TString job;
    if (!GetJobUrl(name, mc, job)) return false;
    
    if (!GetRuns(job, mc, minSize)) return false;

    if (!fPass.IsNull()) fPass.Append("/");
    Printf("alien://%s?pattern=%s*/AliESDs.root&runs=%s", 
	   fPath.Data(), fPass.Data(), fRuns.Data());

    gSystem->Unlink(fTmp);
    return true;
  }
  /** 
   * Get the job url. 
   * 
   * @param name  Production name  
   * @param mc    Should be true for MC 
   * @param url   On return, the job url 
   * 
   * @return true on success 
   */
  Bool_t GetJobUrl(const TString& name, Bool_t mc, TString& url)
  {
    url = "";
    TString index("raw.jsp");
    if (!Download((mc ? "job_details.jsp" : "production/raw.jsp"), index))
	return false;
    
    std::ifstream in(index.Data());
    TString line;
    TString tgt(Form("<td class=\"table_row\">%s</td>", name.Data()));
    do { 
      line.ReadLine(in);
      if (!line.Contains(tgt)) continue;
      line.ReadLine(in);
      Int_t first = line.Index("href=\"");
      Int_t last  = line.Index("\"", first+7);
      url = line(first+6,last-first-6);
      break;
      
    } while (!in.eof());
    in.close();
    
    if (url.IsNull()) { 
      Error("GetJobUrl", "Production %s not found", name.Data());
      return false;
    }

    return true;
  }
  /** 
   * Get list of runs associated with production
   *  
   * @param url     The production URL 
   * @param mc      True of MC
   * @param minSize Least size of runs to use 
   * 
   * @return true on success
   */
  Bool_t GetRuns(const TString& url, Bool_t mc, ULong_t minSize)
  {
    TString index("job");
    if (!Download(Form("%s%s", (mc ? "" : "raw/"), url.Data()), index)) 
      return false;

    std::ifstream in(index.Data());
    TString tgt1(mc ? "window.open" : "runDetails");
    TString line  = "";
    do { 
      line.ReadLine(in);
      if (!line.Contains(tgt1)) continue;
      Int_t   first = -1;
      Int_t   last  = -1;
      if (!mc) { 
	first = line.Index(tgt1);
	last  = line.Index(")", first+tgt1.Length()+1);
      }
      else { 
	Int_t tmp = line.Index(">");
	first = line.Index(">", tmp+1);
	last  = line.Index("<", first);
      }
      if (first == kNPOS || last == kNPOS) { 
	Error("GetDir", "Failed to get directory from %s", line.Data());
	return false;
      }
      first += (mc ? 1 : tgt1.Length()+1);
      last  -= first;
      TString srun  = line(first, last);
      ULong_t runNo = srun.Atoll();
      if (fDebug) Info("", "Got run %lu (%s)", runNo, srun.Data());

      if (!GetSize(in, runNo, mc, minSize)) continue;
      if (!GetDir(in, runNo, mc))           continue;

      if (!fRuns.IsNull()) fRuns.Append(",");
      fRuns.Append(Form("%lu", runNo));
    } while (!in.eof());
    in.close();
    if (fRuns.IsNull()) return false;

    return true;
  }
  /** 
   * Get the size of a given run 
   * 
   * @param in       Input stream 
   * @param runNo    Run number to search for 
   * @param mc       True for simulations 
   * @param minSize  Least size 
   * 
   * @return true on success 
   */
  Bool_t GetSize(std::istream& in, ULong_t runNo, 
		 Bool_t mc, ULong_t minSize=100000)
  {
    TString line;
    TString tgt2(mc ? "table_row_right" : "ESDs size");
    Int_t   cnt = 0;
    do {
      line.ReadLine(in);
      if (!line.Contains(tgt2)) continue;
      cnt++;
      if (mc && cnt < 3) continue;
      if (!mc) line.ReadLine(in);
      if (fDebug) Info("", line);

      TString ssiz;
      if (mc) { 
	Int_t first       = line.Index(">");
	Int_t last        = line.Index("<",first+1);
	if (first == kNPOS || last == kNPOS) { 
	  Error("GetDir", "Failed to get directory from %s", line.Data());
	  return false;
	}
	ssiz = line(first+1, last-first-1);
      }
      else {
	for (Int_t i = 0; i < line.Length(); i++) { 
	  if (line[i] == '<') break;
	  if (line[i] == ' ' || line[i] == '\t' || line[i] == ',') continue;
	  ssiz.Append(line[i]);
	}
      }
      Long_t size = ssiz.Atoll();
      if (fDebug) Info("", "Got run %lu %lu" , runNo, size);
      if (size < 0) {
	Error("GetSize", "Failed to extract size for run %lu", runNo);
	return false;
      }
      if (ULong_t(size) < minSize) {
	Warning("GetSize","Run %lu does not have enough events %lu",runNo,size);
	return false;
      }
      break;
    } while (!in.eof());
    return true;
  }
  /** 
   * Get a directory 
   * 
   * @param in    Input stream
   * @param runNo The run number 
   * @param mc    True for MC 
   * 
   * @return true on success 
   */
  Bool_t GetDir(std::istream& in, ULong_t runNo, Bool_t mc)
  {
    TString line;
    TString tgt3("/catalogue/index.jsp");
    do { 
      line.ReadLine(in);
      // Info("", "line=%s", line.Data());
      if (!line.Contains(tgt3)) continue;
      if (fDebug) Info("", line);
      Int_t tmp         = mc ? line.Index(">")+1 : 0;
      Int_t first       = line.Index(">", tmp);
      Int_t last        = line.Index("<",first+1);
      if (first == kNPOS || last == kNPOS) { 
	Error("GetDir", "Failed to get directory from %s", line.Data());
	return false;
      }
      
      TString dir = line(first+1,last-first-1);
	
      if (fDebug) Info("", "Got run %lu %s", runNo, dir.Data());
      TString path, pass;
      if (!GetPathPass(dir, runNo, path, pass)) return false;
      
      if (fDebug) Info("", "Got run %lu %s %s", runNo,path.Data(),pass.Data());

      if      (fPath.IsNull()) fPath = path;
      else if (!fPath.EqualTo(path)) { 
	Warning("GetDir", "Run %lu location %s not %s", 
	      runNo, path.Data(), fPath.Data());
	return false;
      }

      if      (fPass.IsNull()) fPass = pass;
      else if (!fPass.EqualTo(pass)) { 
	Warning("GetDir", "Run %lu pass %s not %s", 
	      runNo, pass.Data(), fPass.Data());
	return false;
      }
      break;
    } while (!in.eof());
    return true;
  }
  /** 
   * Get the pass from the path 
   * 
   * @param dir   Directory
   * @param run   Run number 
   * @param path  On return, the path 
   * @param pass  On return, the pass 
   * 
   * @return true on success 
   */
  Bool_t GetPathPass(const TString& dir, ULong_t run, 
		     TString& path, TString& pass) 
  {
    Int_t first = dir.Index(Form("%lu", run));
    Int_t last  = dir.Index("/", first);
    if (last == kNPOS) last = dir.Length();
    if (first == kNPOS) { 
      Error("GetPathPass", "Run number %lu not in path %s", run, dir.Data());
      return false;
    }
    while (dir[first-1] == '0') first--;
    path = dir(0, first);
    pass = dir(last+1,dir.Length()-last-1);
    return true;
  }
  /** 
   * Download a file from monalisa 
   * 
   * @param url URL to download 
   * @param out On return, the content of the file 
   * 
   * @return true on success 
   */
  Bool_t Download(const TString& url, TString& out)
  {
    const TString base("http://alimonitor.cern.ch");
    out = Form("%s/%s", fTmp.Data(), out.Data());
    TString cmd(Form("wget -q \"%s/%s\" -O \"%s\"", 
		     base.Data(), url.Data(), out.Data()));
    if (fDebug) Info("Download", "%s", cmd.Data());
    Int_t ret = gSystem->Exec(cmd.Data());
    if (ret != 0)  {
      Error("Download", "Failed to get %s/%s -> %s", 
	    base.Data(), url.Data(), out.Data());
      return false;
    }
    return true;
  }
};

//
// EOF
// 
