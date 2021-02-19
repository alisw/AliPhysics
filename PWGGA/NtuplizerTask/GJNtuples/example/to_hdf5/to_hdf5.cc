#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <H5Cpp.h>

// This is chosen to be the CPU L2 cache size, which should exceed 512
// kB for many years now
#ifndef HDF5_DEFAULT_CACHE
#define HDF5_DEFAULT_CACHE (512 * 1024)
#endif // HDF5_DEFAULT_CACHE

#ifndef HDF5_USE_DEFLATE
#define HDF5_USE_DEFLATE
#endif // HDF5_USE_DEFLATE

// Rank of the tensor written in HDF5, rank 3 being (index_event,
// index_track, index_properties)
#define RANK 3

UInt_t find_ntrack_max(char *argv_first[], char *argv_last[])
{
    UInt_t ntrack_max = 0;

    for (char **p = argv_first; p != argv_last; p++) {
        // Cautious opening of the TTree, capturing all modes of
        // failure, and keep the TDirectoryFile (to be deleted later)
        // to avoid memory leak
        TFile *file = TFile::Open(*p);

        if (file == NULL) {
            continue;
        }

        TDirectoryFile *df = dynamic_cast<TDirectoryFile *>
            (file->Get("AliAnalysisTaskNTGJ"));

        if (df == NULL) {
            continue;
        }

        TTree *_tree_event = dynamic_cast<TTree *>
            (df->Get("_tree_event"));

        if (_tree_event == NULL) {
            continue;
        }

        // Get the maximum of the "ntrack"

        UInt_t ntrack;

        _tree_event->SetBranchAddress("ntrack", &ntrack);

        for (Long64_t i = 0; i < _tree_event->GetEntries(); i++) {
            _tree_event->GetEntry(i);
            ntrack_max = std::max(ntrack_max, ntrack);
        }

        // Fully delete everything

        _tree_event->Delete();
        delete df;
        file->Close();
        delete file;
    }

    return ntrack_max;
}

void write_track(H5::DataSet &data_set, hsize_t *offset,
                 const hsize_t *dim_extend, const UInt_t ntrack_max,
                 char *argv_first[], char *argv_last[])
{
    for (char **p = argv_first; p != argv_last; p++) {
        TFile *file = TFile::Open(*p);

        if (file == NULL) {
            continue;
        }

        TDirectoryFile *df = dynamic_cast<TDirectoryFile *>
            (file->Get("AliAnalysisTaskNTGJ"));

        if (df == NULL) {
            continue;
        }

        TTree *_tree_event = dynamic_cast<TTree *>
            (df->Get("_tree_event"));

        if (_tree_event == NULL) {
            continue;
        }

        UInt_t ntrack;
        std::vector<Float_t> track_e(ntrack_max, NAN);
        std::vector<Float_t> track_pt(ntrack_max, NAN);
        std::vector<Float_t> track_eta(ntrack_max, NAN);
        std::vector<Float_t> track_phi(ntrack_max, NAN);
        std::vector<UChar_t> track_quality(ntrack_max, NAN);
        std::vector<Float_t> track_eta_emcal(ntrack_max, NAN);
        std::vector<Float_t> track_phi_emcal(ntrack_max, NAN);

        _tree_event->SetBranchAddress("ntrack", &ntrack);
        _tree_event->SetBranchAddress("track_e", &track_e[0]);
        _tree_event->SetBranchAddress("track_pt", &track_pt[0]);
        _tree_event->SetBranchAddress("track_eta", &track_eta[0]);
        _tree_event->SetBranchAddress("track_phi", &track_phi[0]);
        _tree_event->SetBranchAddress("track_quality",
                                      &track_quality[0]);
        _tree_event->SetBranchAddress("track_eta_emcal",
                                      &track_eta_emcal[0]);
        _tree_event->SetBranchAddress("track_phi_emcal",
                                      &track_phi_emcal[0]);

        for (Long64_t i = 0; i < _tree_event->GetEntries(); i++) {
            _tree_event->GetEntry(i);

            std::vector<float> data(ntrack_max * 7, NAN);

            for (Long64_t j = 0; j < ntrack; j++) {
                // Note HDF5 is always row-major (C-like)
                data[j * 7 + 0] = track_e[j];
                data[j * 7 + 1] = track_pt[j];
                data[j * 7 + 2] = track_eta[j];
                data[j * 7 + 3] = track_phi[j];
                data[j * 7 + 4] = track_quality[j];
                data[j * 7 + 5] = track_eta_emcal[j];
                data[j * 7 + 6] = track_phi_emcal[j];
            }

            if (offset == 0) {
                // Writing the first event. The data space is already
                // created with space for one event (see when
                // file.createDataSet() was called)
                data_set.write(&data[0], H5::PredType::NATIVE_FLOAT);
            }
            else {
                // The new, extended-by-1 dimension
                const hsize_t dim_extended[RANK] = {
                    offset[0] + dim_extend[0], dim_extend[1],
                    dim_extend[2]
                };

                // Extend to the new dimension
                data_set.extend(dim_extended);

                // Select the hyperslab that only encompass the
                // difference from extending the data space (i.e. the
                // new event, but offsetted at the existing event)
                H5::DataSpace file_space = data_set.getSpace();

                file_space.selectHyperslab(
                    H5S_SELECT_SET, dim_extend, offset);

                // The memory space is the difference only (i.e. also
                // the new event, but at offset 0)
                H5::DataSpace memory_space(RANK, dim_extend, NULL);

                // Write the data from memory space to file space
                data_set.write(&data[0], H5::PredType::NATIVE_FLOAT,
                               memory_space, file_space);

            }
            offset[0]++;

            if (offset[0] % 1000 == 0) {
                fprintf(stderr, "%s:%d: %llu / %lld\n", __FILE__,
                        __LINE__, offset[0],
                        _tree_event->GetEntries());
            }
        }

        _tree_event->Delete();
        delete df;
        file->Close();
        delete file;
    }
}

int main(int argc, char *argv[])
{
    if (argc < 2) {
        exit(EXIT_FAILURE);
    }

    UInt_t ntrack_max = 1927 ; // = find_ntrack_max(argv + 1, argv + argc - 1);

    fprintf(stderr, "%s:%d: ntrack_max = %u\n", __FILE__, __LINE__, ntrack_max);

    // Access mode H5F_ACC_TRUNC truncates any existing file, while
    // not throwing any exception (unlike H5F_ACC_RDWR)
    H5::H5File file(argv[argc - 1], H5F_ACC_TRUNC);
    // How many properties per track is written
    static const size_t row_size = 7;
    // The tensor dimension increment for each new event
    hsize_t dim_extend[RANK] = { 1, ntrack_max, row_size };
    // The maximum tensor dimension, for unlimited number of events
    hsize_t dim_max[RANK] = { H5S_UNLIMITED, ntrack_max, row_size };
    // The extensible HDF5 data space
	H5::DataSpace data_space(RANK, dim_extend, dim_max);

    // To enable zlib compression (there will be many NANs) and
    // efficient chunking (splitting of the tensor into contingous
    // hyperslabs), a HDF5 property list is needed
    H5::DSetCreatPropList property = H5::DSetCreatPropList();

#ifdef HDF5_USE_DEFLATE
    // Check for zlib (deflate) availability and enable only if
    // present
    if (!H5Zfilter_avail(H5Z_FILTER_DEFLATE)) {
        fprintf(stderr, "%s:%d: warning: deflate filter not "
                "available\n", __FILE__, __LINE__);
    }
    else {
        unsigned int filter_info;

        H5Zget_filter_info(H5Z_FILTER_DEFLATE, &filter_info);
        if (!(filter_info & H5Z_FILTER_CONFIG_ENCODE_ENABLED)) {
            fprintf(stderr, "%s:%d: warning: deflate filter not "
                    "available for encoding\n", __FILE__, __LINE__);
        }
        else {
            property.setDeflate(1);
        }
    }
#endif // HDF5_USE_DEFLATE

    // Activate chunking, while observing the HDF5_DEFAULT_CACHE being
    // the CPU L2 cache size
    hsize_t dim_chunk[RANK] = {
        std::max(static_cast<unsigned long long>(1),
                 HDF5_DEFAULT_CACHE /
                 std::max(static_cast<unsigned long long>(1),
                          dim_extend[1] * dim_extend[2] *
                          sizeof(float))),
        dim_extend[1],
        dim_extend[2]
    };

    property.setChunk(RANK, dim_chunk);

    // Create the data set, which will have space for the first event
    H5::DataSet data_set =
        file.createDataSet("track", H5::PredType::NATIVE_FLOAT,
                           data_space, property);
    hsize_t offset[RANK] = {0, 0, 0};

    write_track(data_set, offset, dim_extend, ntrack_max,
                argv + 1, argv + argc - 1);

    file.close();

    return EXIT_SUCCESS;
}
