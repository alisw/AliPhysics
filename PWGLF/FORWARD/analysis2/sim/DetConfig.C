/**
 * @file   DetConfig.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Oct 15 13:19:34 2014
 * 
 * @brief  Particular setup of detectors to use. 
 * 
 * This is used by Simulate.C, Reconstruct.C, QA.C, and AOD.C
 */
/** 
 * Particular set-up of detectors to use.  By default everything but
 * ACORDE is on.  Here, we turn off a few other detectors.
 */
struct DetCfg : public VirtualDetCfg
{
  // virtual Bool_t UseEMCAL()  const { return false; }
  // virtual Bool_t UseMUON()   const { return false; }
  // virtual Bool_t UsePHOS()   const { return false; }
  // virtual Bool_t UseHMPID()  const { return false; }
  // virtual Bool_t UseTOF()    const { return false; }
  // virtual Bool_t UseTRD()    const { return false; }
  const char* GeometrySource() const { return "*OCDB*"; }
};

/** 
 * Create the detector configuration 
 * 
 * @return Pointer to global @c detCfg
 */
void DetConfig()
{
  Info("DetConfig", "Creating detector configuration");
  detCfg = new DetCfg;
  detCfg->Print();
}
// 
// EOF
//

