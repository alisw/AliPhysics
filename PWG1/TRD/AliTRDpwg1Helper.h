#ifndef ALITRDPWG1HELPER_H
#define ALITRDPWG1HELPER_H

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Helper class for PWG1 TRD train                                       //
//                                                                        //
//  Authors:                                                              //
//    Markus Fasel <M.Fasel@gsi.de>                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

class AliTRDpwg1Helper{
public:
  enum ETRDinfoGenOutSlots {
    kEventInfo     = 1
    ,kTracksBarrel  = 2
    ,kTracksSA      = 3
    ,kTracksKink    = 4
    ,kV0List        = 5
    ,kMonitor       = 6
    ,kNOutSlots     = 7
  };

  enum ETRDrecoTasks{
    kCheckESD      = 0
    ,kInfoGen       = 1
    ,kCheckDET      = 2
    ,kEfficiency    = 3
    ,kResolution    = 4
    ,kCheckPID      = 5
    ,kV0Monitor     = 6
    ,kCalibration   = 7
    ,kEfficiencyMC  = 8
    ,kAlignment     = 9
    ,kPIDRefMaker   =10
    ,kClErrParam    =11
    ,kMultiplicity  =12
  };

  enum{
    kNTRDQATASKS = 7,
    kNTRDCALIBTASKS = 6,
    kNTRDTASKS = kNTRDQATASKS + kNTRDCALIBTASKS
  };

  AliTRDpwg1Helper();
  ~AliTRDpwg1Helper();

  static Int_t  GetTaskIndex(const Char_t *name);
  static Bool_t HasReadMCData(Char_t *opt);
  static Bool_t HasReadFriendData(Char_t *opt);
  static const Char_t * TaskOpt(Int_t itask) {return fgkTRDtaskOpt[itask];}
  static const Char_t * TaskClassName(Int_t itask) {return fgkTRDtaskClassName[itask];}

  static void   MergeProd(const Char_t *mark, const Char_t *files, Int_t nBatch);
  static Int_t  ParseOptions(Char_t *trd);

private:
  AliTRDpwg1Helper(const AliTRDpwg1Helper& ref);
  const AliTRDpwg1Helper& operator=(const AliTRDpwg1Helper& ref);

  static const Char_t * fgkTRDtaskOpt[kNTRDTASKS + 1];  //! task options
  static const Char_t * fgkTRDtaskClassName[kNTRDTASKS];//! task class name
};

#endif
