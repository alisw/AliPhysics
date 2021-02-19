# How to Use the Ntuplizer

## Run the Ntuplizer

On a system with the [ALICE software](https://alisw.github.io/alibuild/) installed, get the code by

```bash
git pull https://github.com/yslai/ntuple-gj
```

Set the `alien_CLOSE_SE` (e.g. to `ALICE::LBL::EOS` for users at LBL and UC Berkeley), which will be honored for storing the ntuple output.

Then for a configuration file, e.g. `lhc16k_dijet_skim_1run.yaml`, execute the Grid job as:

```bash
./macros/runNTGJ.C lhc16k_dijet_skim_1run.yaml full
```

to submit onto the ALICE grid. Remove the argument `full` (or replace it with `test`) causes a small subset of the job to be executed locally. The file `macros/runNTGJ.C` is both a ROOT CINT macro, and also a valid shell script. The command shown above will automatically run as `root -l -b -q ./macros/runNTGJ.C("lhc16k_dijet_skim_1run.yaml","full")`.

### Configuration File

The configuration file is in a YAML-like format, and is parsed using a parser that is implemented in a custom parser that is compatible with ROOT 5&rsquo; CINT. The key difference is multiple keys with the same name, e.g. repeated `package: ...` is permissible (and being exploited).

## Run the Ntuple Merger

Do not attempt to merge ntuples on the Grid. They will inevitably fail as ALICE always call `hadd` to merge, which is in-memory. Instead, another macro is provided for disk-resident merging. For example, to merge `rx/lhc16h2a_bis-246994-1/001/AnalysisResults.root`, `rx/lhc16h2a_bis-246994-1/002/AnalysisResults.root`, etc., execute:

```bash
./macros/MergeNtuple.C "rx/lhc16h2a_bis-246994-1/*/*.root" merged_output.root
```

Note to quote the wildcard, in order to prevent its expansion by the shell before calling the macro (otherwise the maximum shell command line length is quickly reached). The wildcard is then expanded (macro internally) via a separate shell, and is not limited to ROOT&rsquo;s `TChain` expansion (which does not allow merging across multiple directories).

`MergeNtuple.C` will prefer ROOT 6, if it is available in the local `aliBuild` environment. It is also possible to compile `MergeNtuple.C` via ROOT ACLiC.
