// $Header$
// Cause ROOT to properly crash and dump core on SigSEGV.

{
 gSystem->IgnoreSignal(kSigSegmentationViolation, true);
}
