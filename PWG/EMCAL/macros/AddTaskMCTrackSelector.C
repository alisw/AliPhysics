AliEmcalMCTrackSelector* AddTaskMCTrackSelector(
  const char *outname    = "mcparticles",
  const char *nTrackCont = "usedefault",
  Bool_t      nk         = kFALSE,
  Bool_t      ch         = kFALSE,
  Double_t    etamax     = 1,
  Bool_t      physPrim   = kTRUE
)
{
  return AliEmcalMCTrackSelector::AddTaskMCTrackSelector(outname, nTrackCont, nk, ch, etamax, physPrim);
}
