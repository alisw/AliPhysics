/** 
 * @defgroup pwglf_forward_trains_util Utilities for Train setups
 *
 * @ingroup pwglf_forward_trains
 */
/**
 * @file   OutputUtilities.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Tue Oct 16 17:55:32 2012
 * 
 * @brief  Special output handling
 * 
 * @ingroup pwglf_forward_trains_util
 */
#ifndef TREEOUTPUTHELPER_C
#define TREEOUTPUTHELPER_C
#ifndef __CINT__
# include <TString.h>
# include <TError.h>
# include <AliAnalysisManager.h>
# include <AliAnalysisDataContainer.h>
# include <AliVEventHandler.h>
# include <TSystem.h>
// Below for finding free port number
# ifdef R__UNIX
#  include <sys/types.h>
#  include <sys/socket.h>
#  include <netinet/in.h>
#  include <netinet/ip.h>
#  include <unistd.h>
# endif
#else
class TString;
#endif

// ===================================================================
/**
 * Special output handling - data sets and remote storage
 *
 * @ingroup pwglf_forward_trains_util
 */
struct OutputUtilities
{
  /** 
   * Register output data set 
   * 
   * @param dsname Data set name 
   * 
   * @return true on success
   */
  static Bool_t RegisterDataset(const TString& dsname)
  {
    // Get the manager
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();

    // If we are asked to make a data-set, get the output handler and
    // common output container.
    AliVEventHandler*         handler = mgr->GetOutputEventHandler();
    if (!handler) return true;

    // Get the container 
    AliAnalysisDataContainer* cont    = mgr->GetCommonOutputContainer();
    if (!cont) { 
      Warning("OutputUtilities::RegisterDataset", 
	      "No common output container defined");
      return false;
    }

    // Make the name 
    TString nme(dsname);
    if (nme.IsNull()) nme = mgr->GetName();
    if (nme.IsNull()) { 
      Error("OutputUtilities::RegisterDataset", "No data set name specified");
      return false;
    }

    // Flag for data-set creation
    cont->SetRegisterDataset(true);

    handler->SetOutputFileName(nme);
    // cont->SetFileName(nme);

    TString base(handler->GetOutputFileName());
    base.ReplaceAll(".root","");
    Info("OutputUtilities::RegisterDataset", 
	 "Will register tree output AODs (%s%s) as dataset",
	 base.Data(), cont->GetTitle());

    return true;
  }
  /** 
   * Get the name of the registered data set
   * 
   * 
   * @return Name of the registered data set
   */
  static TString RegisteredDataset()
  {
    TString ret;

    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    AliVEventHandler*   oh  = mgr->GetOutputEventHandler();
    if (!oh) { 
      Warning("OutputUtilities::GetOutputDataSet", 
	      "No outout event handler defined");
      return ret;
    }
    AliAnalysisDataContainer* co  = mgr->GetCommonOutputContainer();
    if (!co) { 
      Warning("OutputUtilities::GetOutputDataSet", 
	      "No common output container defined");
      return ret;
    }
    if (!co->IsRegisterDataset()) { 
      Info("OutputUtilities::GetOutputDataSet", 
	   "Common output is not registered as dataset");
      return ret;
    }
    ret = oh->GetOutputFileName();
    // ret.ReplaceAll("TTree", "");
    ret.ReplaceAll(".root", "");
    // ret.Append(co->GetTitle());

    return ret;
  }
  static Int_t FindPort()
  {
#ifdef R_UNIX
    int sd = socket(AF_INET,SOCK_STREAM,0);
    if (sd < 0) {
      Warning("FindPort", "Failed to make socket");
      return -1;
    }

    // Binding to port 0 will give us back a free port.  The kernel
    // does not reuse port numbers immediately
    struct sockaddr_in addr;
    addr.sin_family        = AF_INET;
    addr.sin_addr.s_addr   = INADDR_ANY;
    addr.sin_port          = 0;
    int bd = bind(sd, (struct sockaddr*)&addr, sizeof(addr));
    if (bd != 0) {
      Warning("FindPort", "Failed to bind socket to port 0");
      close(sd);
      return -1;
    }

    // Get the address on the bound socket 
    struct sockaddr_in radd;
    socklen_t ladd = sizeof(radd);
    int nr = getsockname(sd, (struct sockaddr*)&radd, &ladd);
    if (nr != 0) {
      Warning("FindPort", "Failed get socket port");
      close(sd);
      return -1;
    }

    int port = radd.sin_port;
    close (sd);

    return port;
#else
    Warning("FindPort", "don't know how to do that on your system");
    return -1;
#endif
  }
  /** 
   * Start a unique XRootd server and return it's access URL 
   * 
   * @param url On return, the access url 
   * 
   * @return true if successful, false otherwise 
   */
  static Bool_t StartXrootd(TString& url)
  {
    url = "";
    Int_t port = FindPort();
    if (port < 0) return false;

    // Get host, current directory, and user name for unique name
    TString host(gSystem->HostName());
    TString dir(gSystem->WorkingDirectory());
    TString name(gSystem->GetUserInfo()->fUser.Data());

    // Form the command line.  Note, we put the PID file one level up,
    // so we know where to look for it. Otherwise it would be put in a
    // sub-directory based on the name of the server.  Since we later
    // on don't know the name of the server we wouldn't now where to
    // look for the PID file
    TString exec;
    exec.Form("xrootd -p %d -l xrd.log -s ../xrd.pid -b -n %s %s",
	      port, name.Data(), dir.Data());
    Info("StartXrootd", "Starting XRootD to serve %s on port %d",
	 dir.Data(), port);
    Info("StartXrootd", "%s", exec.Data());
    int ret = gSystem->Exec(exec);
    if (ret != 0) {
      Warning("StartXrootd", "Failed to start XRootd server");
      return false;
    }
    
    // Form the access URL
    url = Form("root://%s@%s:%d/%s",
	       name.Data(), host.Data(), port, dir.Data());
    Info("StartXrootd", "Access URL is \"%s\"", url.Data());

    return true;
  }
  /** 
   * Stop a previously started Xrootd server 
   * 
   * @return true if stopped, false otherwise 
   */
  static Bool_t StopXrootd()
  {
    std::ifstream pidFile("xrd.pid");
    if (!pidFile) return false;

    TString s; s.ReadFile(pidFile);
    pidFile.close();
    gSystem->Unlink("xrd.pid");
    
    if (s.IsNull()) return false;

    Info("StopXrootd", "Stopping XRootd server (pid: %s)", s.Data());
    return gSystem->Exec(Form("kill -9 %s", s.Data())) == 0;
  }
  /** 
   * Register special putput storage 
   * 
   * @param url Url (root://host/full_path)
   * 
   * @return true on success
   */
  static Bool_t RegisterStorage(const TString& url)
  {
    if (url.IsNull()) { 
      Error("OutputUtilities::RegisterStorage", "No storage URI specified");
      return false;
    }

    // Get the manager
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    
    // Get the container 
    AliAnalysisDataContainer* cont    = mgr->GetCommonOutputContainer();
    if (!cont) { 
      Warning("OutputUtilities::RegisterStorage", 
	      "No common output container defined");
      return false;
    }

    TString u(url);
    if (u.EqualTo("auto")) {
      if (!StartXrootd(u)) {
	Warning("OutputUtilities::RegisterStorage",
		"Couldn't start the XRootD server");
	return false;
      }
    }

    cont->SetSpecialOutput();
    mgr->SetSpecialOutputLocation(u);

    return true;
  }
};
#endif
//
// EOF
//
