README PT resolution correction:
================================

1. Start from filtered high pt trees
    at GSI (example,PbPb:) /lustre/nyx/alice/DETdata/triggeredESD/alice/data/2015/LHC15o/000245145/pass1_pidfix/AOD/007/FilterEvents_Trees.root
    on grid (eample,XeXe:) alien:///alice/data/2017/LHC17n/000280235/pass1/MergedTrees/001/FilterEvents_Trees.root

    and reate a list of those files, see list_tree_pbpb.txt (gsi) or list_tree_xexe.txt (grid) for examples 

2. Extract the normalized pT resoultion with the macro getscaledcov.C
   params: listname = "list_tree.txt"   Filename of the input list
           outfile = "scaledcov.root"   Filename of the output file (existing file is overwritten)
           maxfiles = -1                Max Number of files from list to process (-1 means all)
           tag = ""                     An extra tag to identify this file later

3. Run the macro pts.C to generate the resolution correction factors
   compile this macro because otherwise its very slow!
   params: power = 0.                   exponent of the power law used as input pt spectrum, ignored if a "ptdist" function is provided
           infile = "scaledcov.root"    Filename of inputfile with the scaled pt resolutions and the scaling function
           outfile = "ptrescorr.root"   output file (existing file is overwritten)
           ntracks = 1e8                number of tracks to simulate
           resfactor = 1.               scale the resolution from cov. matrix by this additional factor
           minpt = 7.                   minimal pt to conisder
           extrasmear = 0.              extra smearing of the resolution added in quadrature
           meanshifta = 0.              MEAN q/pt shift TPC A-Side
           sigmashifta = 0.             VARIANCE q/pt shift TPC A-Side
           meanshiftc = 0.              MEAN q/pt shift TPC C-Side
           sigmashiftc = 0.             VARIANCE q/pt shift TPC C-Side
           TF1* ptdist = 0.             provide a TF1 function here that describes the pt spectrum, in this case "power" is ignored
           tag = ""                     An extra tag to identify this file later

TODO: 
- check difference mc resolution vs. cov. matrix
- for low pT part use esd with fine binning (if prev. is ok)