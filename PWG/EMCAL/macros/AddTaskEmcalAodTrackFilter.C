    static AliEmcalAodTrackFilterTask* AddTaskEmcalAodTrackFilter(    
                                  const char *name         = "FilterTracks",
                                  const char *inname       = "tracks",
                                  const char *runperiod    = "", 
                                  const char *taskName     = "AliEmcalAodTrackFilterTask",
                                  const char *suffix       = "" ){
                return AliEmcalAodTrackFilterTask::AddTaskEmcalAodTrackFilter(name,inname,runperiod,taskName,suffix);
                                                                }
