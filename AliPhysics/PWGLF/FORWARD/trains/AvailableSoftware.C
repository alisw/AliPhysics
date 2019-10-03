/**
 * @file   AvailableSoftware.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Tue Oct 16 17:54:11 2012
 * 
 * @brief  Find available packages
 * 
 * @ingroup pwglf_forward_trains_util
 */

#ifndef AVAILABLESOFTWARE_C
#define AVAILABLESOFTWARE_C
#ifndef __CINT__
# include <TString.h>
# include <TSystem.h>
# include <TError.h>
# include <TObjArray.h>
# include <TList.h>
# include <TPRegexp.h>
# include <TObjString.h>
# include <fstream>
#else
class TString;
class TList;
#endif

/**
 * Railway code to find available packages on Grid
 * 
 *
 * @ingroup pwglf_forward_trains_util
 */
struct AvailableSoftware
{
  static const char* GetField(TObject* o)
  {
    static TString r;
    if (!o) return "";
    TString s(o->GetName());
    TString t = s.Strip(TString::kBoth,' ');
    r = t.Strip(TString::kBoth, '\t');
    return r.Data();
  }
  /** 
   * Get the full map of packages to dependencies
   * 
   * @return the list
   */
  static TList* GetMap()
  {
    static TList l;
    if (l.GetEntries() > 0) return &l;
    TString c("wget -q http://alimonitor.cern.ch/packages/ -O - | "
	      "sed -n -e '/<tr class=table_row/,/<\\/tr>/ p' | "
	      "sed -n -e '/<td/,/<\\/td>/ p' | "
	      "sed -e '/\\/*td.*>/d' | "
	      "sed -e 's/<a.*>\\(.*\\)<\\/a>/\\1/'");
    TString    values = gSystem->GetFromPipe(c);
    TObjArray* tokens = values.Tokenize("\n");
    Int_t      n      = tokens->GetEntries();
    std::ofstream out("values");
    out << values << std::endl;
    out.close();
    // Printf("%s", values.Data());
    // tokens->ls();
    for (Int_t i = 0; i < n; i += 1) { // 2-3 lines per package
      TObject* opack = tokens->At(i+0);
      TString pack = GetField(opack);
      if (pack.IsNull()) {
	i += 2;
	continue;
      }
      TObjString* odeps = static_cast<TObjString*>(tokens->At(i+1));
      TObjString* oblank= static_cast<TObjString*>(tokens->At(i+2));
      TObjString* oavail= static_cast<TObjString*>(tokens->At(i+3));
      TObjString* odate = static_cast<TObjString*>(tokens->At(i+4));
      i += 4;
      if (oblank) {
	TString blk = GetField(oblank);
	if (!blk.IsNull())
	  ::Warning("", "blanck is not blanck: \"%s\"", oblank->GetName());
      }
      if (oavail) {
	TString avail = GetField(oavail);
	if (!avail.EqualTo("Available")) continue;
      }
      
      // TString& pack = opack->String();
      pack.ReplaceAll("VO_ALICE@", "");
      
      // TString& deps = odeps->String();
      TString deps = GetField(odeps);
      deps.ReplaceAll("VO_ALICE@", "");
      if (!(pack.BeginsWith("AliPhysics") ||
	    pack.BeginsWith("AliRoot") ||
	    pack.BeginsWith("ROOT"))) continue;
      if (pack.Contains(".post_install")) continue;

      l.Add(new TNamed(pack, deps));
      // Printf("%-30s: %s", pack.Data(), deps.Data());
    }
    l.Sort();
    // l.ls();
    tokens->Delete();
    delete tokens;
    
    return &l;
  }
  /** 
   * Get a package 
   * 
   * @param name   Package Name  
   * @param query  Query.  Either a specific version or some combination of 
   *
   * - last: the newest
   * - regular: No special tags 
   * - release: Only release tags 
   * - analysis: Only analysis tags
   * 
   * @return Pointer to package or null
   */
  static TObject* GetPackage(const TString& name,
			    const TString& query)
  {
    TList* l = GetMap();
    Bool_t list = (query.Contains("list",       TString::kIgnoreCase) || 
		   query.Contains("help",       TString::kIgnoreCase));
    Bool_t last = (query.Contains("last",       TString::kIgnoreCase) || 
		   query.Contains("newest",     TString::kIgnoreCase));
    Bool_t nots = (query.Contains("nonspecial", TString::kIgnoreCase) ||
		   query.Contains("regular",    TString::kIgnoreCase) ||
		   query.Contains("normal",     TString::kIgnoreCase) ||
		   query.Contains("standard",   TString::kIgnoreCase));
    Bool_t rele = (query.Contains("release",    TString::kIgnoreCase));
    Bool_t anat = (query.Contains("analysis",   TString::kIgnoreCase));
    if (!query.IsNull() && (!list && !last && !nots && !rele && !anat)) {
      Info("GetPackage", "%s=%s already specified, leaving that",
	   name.Data(), query.Data());
      TObject* o = l->FindObject(Form("%s::%s",name.Data(),query.Data()));
      return o;
    }

    TString relPat;
    relPat.Form("%s::v[0-9]-[0-9]{2}+-(Rev-|)[0-9]{2}", name.Data());
    if (name.Contains("AliPhysics")) relPat.Append("-[0-9]{2}");
    if (name.Contains("ROOT"))       relPat.Append("(-alice[0-9]*|)");
    relPat.Append("(-[0-9]+|)$");
    TPRegexp pRele(relPat);
    TPRegexp pAnat(Form("%s::vAN-[0-9]{8}(-[0-9]+|)$", name.Data()));
    TString  vers(Form("%s::%s", name.Data(), query.Data()));

    if (list) {
      TString qual;
      if (nots) qual.Append("regular ");
      if (rele) qual.Append("release ");
      if (anat) qual.Append("analysis ");
      
      Printf("Available %sversion of %s", qual.Data(), name.Data());
    }
    TIter    prev(l, kIterBackward);
    TObject* o = 0;
    TObject* r = 0;
    Bool_t   m = false;
    while ((o = prev())) {
      TString n(o->GetName());
      if (!n.BeginsWith(name)) {
	if (m) break;
	continue;
      }

      // We found the package 
      m = true;

      if (last || list) {
	Bool_t isRele = pRele.MatchB(n);
	Bool_t isAnat = pAnat.MatchB(n);
	Bool_t isSpec = !(isRele || isAnat);
	if (nots && isSpec)  continue;
	if (anat && !isAnat) continue;
	if (rele && !isRele) continue;
	if (list) {
	  n.ReplaceAll(Form("%s::", name.Data()), "");
	  Printf("\t%s", n.Data());
	  continue;
	}
	r = o;
	break;
      }
      if (!vers.EqualTo(n)) continue;
      r = o;
      break;
    }
    if (!r && !list) 
      Warning("GetPackage", "No match found for %s", vers.Data());
    return r;
  }
  static Bool_t GetVer(TObject* pack, const TString& name, TString& ret)
  {
    if (!pack) {
      // if (!ret.IsNull())    ret = "";
      return false;
    }
    
    ret = pack->GetName();
    ret.ReplaceAll(Form("%s::", name.Data()), "");
    return true;
  }
    
  /** 
   * Get the dependencies 
   * 
   * @param pack   Package 
   * @param which  Which dependency
   * @param ret    Return version of dependency
   * 
   * @return true on success
   */
  static Bool_t GetDep(TObject* pack, const TString& which, TString& ret)
  {
    if (!pack) return false;
    TString deps(pack->GetTitle());
    TObjArray* tokens = deps.Tokenize(",");
    TIter next(tokens);
    TObjString* s = 0;
    while ((s = static_cast<TObjString*>(next()))) {
      TString t = s->String().Strip(TString::kBoth);
      if (t.BeginsWith(which)) {
	ret = t;
	ret.ReplaceAll(Form("%s::",which.Data()), "");
	break;
      }
    }
    tokens->Delete();
    delete tokens;
    return !(ret.IsNull());
  }
  static Bool_t Check(TString& aliphysics,
		      TString& aliroot,
		      TString& root)
  {
    // Figure out what to do.  
    Bool_t show = (aliphysics.Contains("list", TString::kIgnoreCase) ||
		   aliroot.Contains("list", TString::kIgnoreCase) ||
		   root.Contains(   "list", TString::kIgnoreCase));    

    TObject* foundPhysics = GetPackage("AliPhysics", aliphysics);
    GetVer(foundPhysics, "AliPhysics", aliphysics);
    GetDep(foundPhysics, "AliRoot", aliroot);

    TObject* foundAliRoot = GetPackage("AliRoot", aliroot);
    GetVer(foundAliRoot, "AliRoot", aliroot);
    GetDep(foundAliRoot, "ROOT", root);

    TObject* foundRoot = GetPackage("ROOT", root);
    GetVer(foundRoot, "ROOT", root);

    if (show) return false;

    if (aliphysics.IsNull() ||
	aliroot.IsNull()    ||
	root.IsNull()) {
      Warning("Check", "Missing packages (%s,%s,%s)",
	      aliphysics.Data(), aliroot.Data(), root.Data());
      return false;
    }
    return true;
  }
  static void Test(const TString& phy,
		   const TString& ali=TString(),
		   const TString& roo=TString())
  {
    TString aliphysics(Form("list,%s",  phy.Data()));
    TString aliroot   (Form("list,%s",  ali.Data()));
    TString root      (Form("list,%s",  roo.Data()));
    Printf("Checking with AliPhysics=%s AliROOT=%s ROOT=%s",
	   phy.Data(), ali.Data(), roo.Data());
    AvailableSoftware::Check(aliphysics,aliroot, root);

    aliphysics = Form("last,%s",phy.Data());
    aliroot    = Form("last,%s",ali.Data());
    root       = Form("last,%s",roo.Data());
    if (AvailableSoftware::Check(aliphysics,aliroot, root))
      Printf("Got AliPhysics=%s AliROOT=%s ROOT=%s",
	     aliphysics.Data(), aliroot.Data(), root.Data());
  }
  
  static void Test()
  {
    Printf("All available");
    AvailableSoftware::Test("");
    
    Printf("All regular");
    AvailableSoftware::Test("regular","regular","regular");

    Printf("All releases");
    AvailableSoftware::Test("release","release","release");

    Printf("All analysis tags");
    AvailableSoftware::Test("analysis","analysis","analysis");
  }
};
#endif
