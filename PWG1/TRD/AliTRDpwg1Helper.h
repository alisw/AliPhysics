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
    ,kTracksBarrel
    ,kTracksSA
    ,kTracksKink
    ,kV0List 
    ,kMonitor
    ,kNOutSlots
  };

  enum ETRDrecoTasks{
    kCheckESD = 0
    ,kInfoGen
    ,kCheckDET
    ,kEfficiency
    ,kResolution
    ,kCheckPID
    ,kV0Monitor
    ,kCheckTRK
    ,kCalibration
    ,kEfficiencyMC
    ,kAlignment
    ,kPIDRefMaker
    ,kClErrParam
    ,kMultiplicity
  };

  enum{
    kNTRDQATASKS = 8,
    kNTRDCALIBTASKS = 6,
    kNTRDTASKS = kNTRDQATASKS + kNTRDCALIBTASKS
  };

  AliTRDpwg1Helper();
  ~AliTRDpwg1Helper();

  static Bool_t DoTask(Int_t idx, Int_t map);
  static Int_t  GetTaskIndex(const Char_t *name);
  static Bool_t HasReadMCData(const Char_t *opt);
  static Bool_t HasReadFriendData(const Char_t *opt);
  static const Char_t * TaskOpt(Int_t itask) {return fgkTRDtaskOpt[itask];}
  static const Char_t * TaskClassName(Int_t itask) {return fgkTRDtaskClassName[itask];}

  static const Char_t*  Basename(const char* filepath);
  static const Char_t*  Dirname(const char* filepath);
  static Int_t  MergeBatch(const Char_t *mark, const Char_t *files, const Int_t nfiles=20, const Int_t first=0, Bool_t kSVN=kTRUE, Bool_t kCLEAR=kFALSE);
  static void   MergeProd(const Char_t *mark, const Char_t *files, const Int_t nBatch=20, Int_t level=0);
  static Int_t  ParseOptions(const Char_t *trd);

private:
  AliTRDpwg1Helper(const AliTRDpwg1Helper& ref);
  const AliTRDpwg1Helper& operator=(const AliTRDpwg1Helper& ref);
  static const Char_t * fgkTRDtaskOpt[kNTRDTASKS + 1];  //! task options
  static const Char_t * fgkTRDtaskClassName[kNTRDTASKS];//! task class name
};

#endif
