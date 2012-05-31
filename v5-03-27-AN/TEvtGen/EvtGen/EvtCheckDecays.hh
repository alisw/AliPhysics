//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed 
//      for the BaBar collaboration.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Module: EvtCheckDecays
//
// Description:  Holds code to conduct various checks on the 
//      EvtDecayTable::decaytable()
//
// Modification history:
//      Abi Soffer             Nov 29, 2007, created
//
//------------------------------------------------------------------------

#ifndef EVTCHECKDECAYS
#define EVTCHECKDECAYS

class EvtId;

class EvtCheckDecays {
public:
// check CP conservation in the decay BRs, daughters, and models:
  static void checkConj(bool compareArguments = false);  

  // Used by checkConj() to identify self-conjugate particles:
  static bool selfConj(const EvtId & id);

private:
};

#endif
