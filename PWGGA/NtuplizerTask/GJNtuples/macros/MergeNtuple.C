#if 0
set +o posix; function join_by { local d=$1; shift; echo -n "$1";
    shift; printf "%s" "${@/#/$d}"; }
root=root; root6="$ROOT_BASEDIR/root6-1/bin/root";
[[ -x "$root6" ]] && root="$root6"; exec $root -l -b -q \
    "$0($(join_by \",\" \"$*\" | /bin/sed \
s/\"\\\([0-9]\\+\\\)\"/\\1/g\;s/^\"\"\$//))"; exit 0
#endif

// Low memory complexity ntuple merger. Adapted from CMS Heavy Ion's
// ntuple merger (https://github.com/richard-cms/hiForestMerging),
// with performance optimization and compatibility with ROOT 5.x/CINT.

#include <TSystem.h>
#include <TChain.h>
#include <TFile.h>
#include <linux/limits.h>

void chain_add_glob(TChain &chain, const char *pattern,
                    int nentries = -1)
{
    const TString command = TString("ls -1 ") + pattern;
    FILE *pipe = gSystem->OpenPipe(command, "r");

    if (pipe == NULL) {
        return;
    }

    char line[PATH_MAX + 1];

    while (fgets(line, PATH_MAX + 1, pipe) != NULL) {
        line[PATH_MAX] = '\0';
        line[strcspn(line, "\r\n")] = '\0';
        // Passing Long64_t does not appear to be handled gracefully
        // by ROOT
        if (nentries == -1) {
            chain.Add(line);
        }
        else {
            chain.Add(line, nentries);
        }
    }
    gSystem->ClosePipe(pipe);
}

void search_ntuple(TObjArray &dir_name, TObjArray &dir_tree_name,
                   const char *glob)
{
    TChain test_chain("AliAnalysisTaskNTGJ/_tree_event");

    fprintf(stderr, "Search for ntuples... ");
    // Maximum one event as to avoid loading all files
    chain_add_glob(test_chain, glob, 1);
    if (!(test_chain.GetEntries() > 0)) {
        fprintf(stderr, "error: no entries found, abort.\n");
        gSystem->Exit(1);
    }

    TFile *test_file = test_chain.GetFile();
    const TList *key_list_1 = test_file->GetListOfKeys();

    for (Int_t i = 0; i < key_list_1->GetEntries(); i++) {
        const TDirectoryFile *dir_file = (TDirectoryFile *)
            test_file->Get(key_list_1->At(i)->GetName());

        if (strcmp(dir_file->ClassName(), "TDirectoryFile") != 0) {
            continue;
        }

        const TList *key_list_2 = dir_file->GetListOfKeys();

        for (Int_t j = 0; j < key_list_2->GetEntries(); j++) {
            const TString n = dir_file->GetName() + TString("/") +
                key_list_2->At(j)->GetName();
            const TTree *t = (TTree *)test_file->Get(n);

            if (t != NULL &&
                (strcmp(t->ClassName(), "TTree") == 0 ||
                 strcmp(t->ClassName(), "TNtuple") == 0) &&
                ((dir_tree_name.GetEntries() == 0) ||
                 // skip duplicate tree entries
                 n != *(TString *)dir_tree_name.Last())) {
                dir_name.AddLast((TObject *)
                                 (new TString(dir_file->GetName())));
                dir_tree_name.AddLast((TObject *)(new TString(n)));
            }
        }
    }

    fprintf(stderr, "done.\n");
}

void merge_ntuple(const char *output_filename, TObjArray &dir_name,
                  TObjArray &dir_tree_name, const char *glob)
{
    TObjArray chain;

    for (Int_t i = 0; i < dir_tree_name.GetEntries(); i++) {
        const TString *dtn = (TString *)dir_tree_name.At(i);

        chain.AddLast((TObject *)(new TChain(dtn->Data())));
        chain_add_glob(*((TChain *)chain.Last()), glob);
        fprintf(stderr, "Created merged ntuple %d: %s\n", i,
                dtn->Data());
    }

    TFile output_file(output_filename, "RECREATE");

    for (Int_t i = 0; i < dir_tree_name.GetEntries(); i++) {
        output_file.cd();

        const TString *dn = (TString *)dir_name.At(i);
        const TString *dtn = (TString *)dir_tree_name.At(i);

        fprintf(stderr, "Processing %s... ", dtn->Data());

        if (i == 0 || *dtn != *(TString *)dir_tree_name.At(i - 1)) {
            (output_file.mkdir(*dn))->cd();
        }
        else {
            output_file.cd(*dn);
        }
        ((TChain *)chain[i])->Merge(&output_file, 0, "keep");

        fprintf(stderr, "done.\n");

        delete (TChain *)chain[i];
        delete dn;
        delete dtn;
    }
    output_file.Close();

    fprintf(stderr, "Output %s produced.\n", output_filename);
}

void MergeNtuple(
    const char *glob = NULL, const char *output_filename = NULL,
    const bool extended_glob = true)
{
    if (glob == NULL || output_filename == NULL) {
        fprintf(stderr, "Usage: MergeNtuple <glob> "
                "<output_filename>\nNote to quote <glob>, e.g. "
                "\"lhc15o/*/AnalysisResults.root\"\n");
        gSystem->Exit(1);
    }

    fprintf(stderr, "Merging %s => %s\n", glob, output_filename);

    TObjArray dir_name;
    TObjArray dir_tree_name;

    search_ntuple(dir_name, dir_tree_name, glob);
    merge_ntuple(output_filename, dir_name, dir_tree_name, glob);
    gSystem->Exit(0);
}
