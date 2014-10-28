/**
 * Configuration of AOD pass.  The base class VirtualAODCfg is
 * declared in AOD.C
 */
struct AODCfg : public VirtualAODCfg
{
  /** 
   * @{ 
   * @name Plug-in settings 
   * Settings that make sense when using the Alien plugin
   */
  /** @return Connect to CDB */
  virtual Bool_t UseCDBconnect() const { return true; }
  /** @return use physics selection */
  virtual Bool_t UsePhysicsSelection() const { return kTRUE; }
  /** @return use tender wagon */
  virtual Bool_t UseTender() const { return kFALSE; }
  /** @return centrality */
  virtual Bool_t UseCentrality() const { return kTRUE; }
  /** @return use V0 correction in tender */
  virtual Bool_t UseV0tender() const { return kFALSE; }
  /** @return activate debugging */
  virtual Bool_t UseDBG() const { return kTRUE; }
  /** @return use MC info */
  virtual Bool_t UseMC() const { return kTRUE; }
  /** @return use Kinematics filter */
  virtual Bool_t UseKFILTER() const { return kTRUE; }
  /** @return use track references */
  virtual Bool_t UseTR() const { return kTRUE; }
  /** @return do not change */
  virtual Bool_t UseCORRFW() const { return kFALSE; }
  /** @return use AOD tags */
  virtual Bool_t UseAODTAGS() const { return kFALSE; }
  /** @return use sys info */
  virtual Bool_t UseSysInfo() const { return kFALSE; }
  /* @} */
  
  /** 
   * @{ 
   * @name Modules 
   *  Analysis modules to be included. Some may not be yet fully implemented.
   */
  /** @return Analysis produces an AOD or dAOD's */
  virtual Bool_t UseAODhandler() const { return true; }
  /** @return ESD to AOD filter (barrel + muon tracks) */
  virtual Bool_t UseESDfilter() const { return true; }
  /** @return Use Muon train  */
  virtual Bool_t UsePWGMuonTrain() const { return false; }
  /** @return Task that copies only muon events */
  virtual Bool_t UseMUONcopyAOD() const { return false; }
  /** @return Jet analysis (PWG4) */
  virtual Bool_t UseJETAN() const { return false; }
  /** @return Jet delta AODs */
  virtual Bool_t UseJETANdelta() const { return false; }
  /** @return Vertexing HF task (PWG3) */
  virtual Bool_t UsePWGHFvertexing() const { return false; }
  /** @return JPSI filtering (PWG3) */
  virtual Bool_t UsePWGDQJPSIfilter() const { return false; }
  /** @return D0->2 hadrons (PWG3) */
  virtual Bool_t UsePWGHFd2h() const { return false; }
  /** @return PID response */
  virtual Bool_t UsePIDResponse() const { return false; }
  /** @return Forward mult task (PWGLF) */
  virtual Bool_t UsePWGLFForward() const { return true; }
  /* @} */
};

/** 
 * Creating our configuration 
 * 
 */
void AODConfig()
{
  Info("AODConfig", "Creating configuration object");
  // MUST create the global object "aodCfg" here!
  aodCfg = new AODCfg();
}
// 
// EOF
// 
