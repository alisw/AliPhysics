/*! \page READMEunfolding Unfolding

__Page under construction__

Unfolding is generally performed via RooUnfold.

# FAQ

1. Do I need to normalize the response before passing to RooUnfold?

   No. RooUnfold will normalize the response when appropriate for the given method. For 1D Bayes, the response is normalized via the code
   block listed below (either in the `RooUnfoldBayes::train(...)` or `RooUnfoldBayes::unfold(...)`, depending on the particular version.

   ~~~{.cxx}
   TMatrixD PEjCi(_ne,_nc), PEjCiEff(_ne,_nc);
   for (Int_t i = 0 ; i < _nc ; i++) {
       if (_nCi[i] <= 0.0) { _efficiencyCi[i] = 0.0; continue; }
       Double_t eff = 0.0;
       for (Int_t j = 0 ; j < _ne ; j++) {
           Double_t response = _Nji(j,i) / _nCi[i];
           PEjCi(j,i) = PEjCiEff(j,i) = response;  // efficiency of detecting the cause Ci in Effect Ej
           eff += response;
       }
       _efficiencyCi[i] = eff;
       Double_t effinv = eff > 0.0 ? 1.0/eff : 0.0;   // reset PEjCiEff if eff=0
       for (Int_t j = 0 ; j < _ne ; j++) PEjCiEff(j,i) *= effinv;
   }
   ~~~

   For SVD, the response should not be normalized.

2. Do I need to explicitly correct for the kinematic efficiency?

   If you provide a fully efficient prior as a separate 1D histogram to the RooUnfoldResponse, then you don't need to correct for the kinematic
   efficiency. RooUnfold will take care of it. However, if you don't provide the fully efficient prior (or you are in the 2D case where you
   cannot provide it), then you do need to correct for the kinematic efficiency after unfolding.

*/
