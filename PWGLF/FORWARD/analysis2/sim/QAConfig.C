/**
 * @file   QAConfig.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Oct 15 13:27:27 2014
 * 
 * @brief  Configuratin of QA pass 
 */
/** 
 * Configuration of which QA tasks to run. Base class is declared in QA.C 
 */
struct QACfg : public VirtualQACfg
{
  /** @return  */
  virtual Bool_t DoCDBconnect()  const { return true; }
  /** @return  */
  virtual Bool_t DoEventStat()   const { return true; }
  /** @return  */
  virtual Bool_t DoCentrality()  const { return true; }
  /** @return  */
  virtual Bool_t DoQAsym()       const { return false; }
  /** @return  there is a 2nd file */
  virtual Bool_t DoVZERO()       const { return true; }
  /** @return  */
  virtual Bool_t DoVZEROPbPb()   const { return false; }
  /** @return  */
  virtual Bool_t DoVertex()      const { return true; }
  /** @return  needs RP    */
  virtual Bool_t DoSPD()         const { return true; }
  /** @return  */
  virtual Bool_t DoTPC()         const { return true; }
  /** @return  */
  virtual Bool_t DoHLT()         const { return true; }
  /** @return  needs RP */
  virtual Bool_t DoSDD()         const { return true; }
  /** @return  */
  virtual Bool_t DoSSDdEdx()     const { return true; }
  /** @return  */
  virtual Bool_t DoTRD()         const { return true; }
  /** @return  */
  virtual Bool_t DoITS()         const { return true; }
  /** @return  */
  virtual Bool_t DoITSsaTracks() const { return true; }
  /** @return  */
  virtual Bool_t DoITSalign()    const { return true; }
  /** @return  */
  virtual Bool_t DoCALO()        const { return true; }
  /** @return  */
  virtual Bool_t DoMUONTrig()    const { return true; }
  /** @return  */
  virtual Bool_t DoImpParRes()   const { return true; }
  /** @return  */
  virtual Bool_t DoMUON()        const { return true; }
  /** @return  */
  virtual Bool_t DoTOF()         const { return true; }
  /** @return  */
  virtual Bool_t DoHMPID()       const { return true; }
  /** @return  */
  virtual Bool_t DoT0()          const { return true; }
  /** @return  */
  virtual Bool_t DoZDC()         const { return true; }
  /** @return  */
  virtual Bool_t DoPIDResponse() const { return true; }
  /** @return  */
  virtual Bool_t DoPIDqa()       const { return true; }
  /** @return  */
  virtual Bool_t DoFWD()         const { return true; }
  /** @return  */
  virtual Bool_t DoPHOS()        const { return true; }
  /** @return  */
  virtual Bool_t DoPHOSTrig()    const { return true; }
  /** @return  */
  virtual Bool_t DoEMCAL()       const { return false; }
  /** @return  */
  virtual Bool_t DoFBFqa()       const { return true; }
  /** @return  NEEDS geometry */
  virtual Bool_t DoMUONEff()     const { return false; }
  /** @return  NEEDS MCtruth  */
  virtual Bool_t DoV0()          const { return false; }
  /** @return Get Debug level */
  virtual Int_t DebugLevel() const { return 3; }
};

/** 
 * Create our QA configuration 
 * 
 */
void QAConfig()
{
  Info("QAConfig", "Creating configuration object");
  qaCfg = new QACfg();
}

// 
// EOF
// 
