# Usage of the macro

## Requirements
Before running the macro make sure to have:
- the output file obtained running the task on the data,
- the output file obtained running the task on the general purpose MC,
- the output file obtained running the task on the Nuclei injected MC,
- a working _ROOT6_ installation.

## Usage
By executing``RunPlotting.C`` the raw primary signal of (anti-)proton and (anti-)deuteron is computed by taking the methods from the following files:
1. ``Functions.C`` defines all needed fit functions.
2. ``RawParticleSpectra.C``  Extracts the signal from the data information.
3. ``Secondary.C``  Computes the correction due to secondary and material particles.
4. ``Summary.C``  Generates the raw primary spectra (ratios) and produces plots
