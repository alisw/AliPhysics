AliEmcalRejectMCBackground* AddTaskRejectMCBackground(const char *outname = "mcparticlesbgrej", const Int_t signalRejection = 0, const Int_t debug = 0){
  return AliEmcalRejectMCBackground::AddTaskRejectMCBackground(outname, signalRejection, debug);
}
