#include <string>
#include <unistd.h>
#include <dlfcn.h>

/////////////////////////////////////////////////////////////////////

#ifdef HAVE_CUDA
#include <cuda_runtime_api.h>
#include <cublas_v2.h>

#define CASE_MSG(s)	{ case s: msg = #s; break; }

#define CHECK_CUDA(c) {                                           \
    cudaError_t status = (c);                                     \
    if (status != cudaSuccess) {                                  \
        const char *msg = NULL;                                   \
        switch (status) {                                         \
        CASE_MSG(cudaSuccess);                                    \
        CASE_MSG(cudaErrorMissingConfiguration);                  \
        CASE_MSG(cudaErrorMemoryAllocation);                      \
        CASE_MSG(cudaErrorInitializationError);                   \
        CASE_MSG(cudaErrorLaunchFailure);                         \
        CASE_MSG(cudaErrorPriorLaunchFailure);                    \
        CASE_MSG(cudaErrorLaunchTimeout);                         \
        CASE_MSG(cudaErrorLaunchOutOfResources);                  \
        CASE_MSG(cudaErrorInvalidDeviceFunction);                 \
        CASE_MSG(cudaErrorInvalidConfiguration);                  \
        CASE_MSG(cudaErrorInvalidDevice);                         \
        CASE_MSG(cudaErrorInvalidValue);                          \
        CASE_MSG(cudaErrorInvalidPitchValue);                     \
        CASE_MSG(cudaErrorInvalidSymbol);                         \
        CASE_MSG(cudaErrorMapBufferObjectFailed);                 \
        CASE_MSG(cudaErrorUnmapBufferObjectFailed);               \
        CASE_MSG(cudaErrorInvalidHostPointer);                    \
        CASE_MSG(cudaErrorInvalidDevicePointer);                  \
        CASE_MSG(cudaErrorInvalidTexture);                        \
        CASE_MSG(cudaErrorInvalidTextureBinding);                 \
        CASE_MSG(cudaErrorInvalidChannelDescriptor);              \
        CASE_MSG(cudaErrorInvalidMemcpyDirection);                \
        CASE_MSG(cudaErrorAddressOfConstant);                     \
        CASE_MSG(cudaErrorTextureFetchFailed);                    \
        CASE_MSG(cudaErrorTextureNotBound);                       \
        CASE_MSG(cudaErrorSynchronizationError);                  \
        CASE_MSG(cudaErrorInvalidFilterSetting);                  \
        CASE_MSG(cudaErrorInvalidNormSetting);                    \
        CASE_MSG(cudaErrorMixedDeviceExecution);                  \
        CASE_MSG(cudaErrorCudartUnloading);                       \
        CASE_MSG(cudaErrorUnknown);                               \
        CASE_MSG(cudaErrorNotYetImplemented);                     \
        CASE_MSG(cudaErrorMemoryValueTooLarge);                   \
        CASE_MSG(cudaErrorInvalidResourceHandle);                 \
        CASE_MSG(cudaErrorNotReady);                              \
        CASE_MSG(cudaErrorInsufficientDriver);                    \
        CASE_MSG(cudaErrorSetOnActiveProcess);                    \
        CASE_MSG(cudaErrorInvalidSurface);                        \
        CASE_MSG(cudaErrorNoDevice);                              \
        CASE_MSG(cudaErrorECCUncorrectable);                      \
        CASE_MSG(cudaErrorSharedObjectSymbolNotFound);            \
        CASE_MSG(cudaErrorSharedObjectInitFailed);                \
        CASE_MSG(cudaErrorUnsupportedLimit);                      \
        CASE_MSG(cudaErrorDuplicateVariableName);                 \
        CASE_MSG(cudaErrorDuplicateTextureName);                  \
        CASE_MSG(cudaErrorDuplicateSurfaceName);                  \
        CASE_MSG(cudaErrorDevicesUnavailable);                    \
        CASE_MSG(cudaErrorInvalidKernelImage);                    \
        CASE_MSG(cudaErrorNoKernelImageForDevice);                \
        CASE_MSG(cudaErrorIncompatibleDriverContext);             \
        CASE_MSG(cudaErrorPeerAccessAlreadyEnabled);              \
        CASE_MSG(cudaErrorPeerAccessNotEnabled);                  \
        CASE_MSG(cudaErrorDeviceAlreadyInUse);                    \
        CASE_MSG(cudaErrorProfilerDisabled);                      \
        CASE_MSG(cudaErrorProfilerNotInitialized);                \
        CASE_MSG(cudaErrorProfilerAlreadyStarted);                \
        CASE_MSG(cudaErrorProfilerAlreadyStopped);                \
        CASE_MSG(cudaErrorAssert);                                \
        CASE_MSG(cudaErrorTooManyPeers);                          \
        CASE_MSG(cudaErrorHostMemoryAlreadyRegistered);           \
        CASE_MSG(cudaErrorHostMemoryNotRegistered);               \
        CASE_MSG(cudaErrorOperatingSystem);                       \
        CASE_MSG(cudaErrorPeerAccessUnsupported);                 \
        CASE_MSG(cudaErrorLaunchMaxDepthExceeded);                \
        CASE_MSG(cudaErrorLaunchFileScopedTex);                   \
        CASE_MSG(cudaErrorLaunchFileScopedSurf);                  \
        CASE_MSG(cudaErrorSyncDepthExceeded);                     \
        CASE_MSG(cudaErrorLaunchPendingCountExceeded);            \
        CASE_MSG(cudaErrorNotPermitted);                          \
        CASE_MSG(cudaErrorNotSupported);                          \
        CASE_MSG(cudaErrorHardwareStackError);                    \
        CASE_MSG(cudaErrorIllegalInstruction);                    \
        CASE_MSG(cudaErrorMisalignedAddress);                     \
        CASE_MSG(cudaErrorInvalidAddressSpace);                   \
        CASE_MSG(cudaErrorInvalidPc);                             \
        CASE_MSG(cudaErrorIllegalAddress);                        \
        CASE_MSG(cudaErrorInvalidPtx);                            \
        CASE_MSG(cudaErrorInvalidGraphicsContext);                \
        CASE_MSG(cudaErrorStartupFailure);                        \
        CASE_MSG(cudaErrorApiFailureBase);                        \
        default: break;                                           \
        }                                                         \
        if (msg == NULL) {                                        \
            fprintf(stderr, "%s:%d: CUDA error (%s)\n",           \
                    __FILE__, __LINE__, msg);                     \
        }                                                         \
        else {                                                    \
            fprintf(stderr, "%s:%d: CUDA error (%d)\n",           \
                    __FILE__, __LINE__, status);                  \
        }                                                         \
    }                                                             \
}

#define CHECK_CUBLAS(c) {                                         \
    cublasStatus_t status = (c);                                  \
    if (status != CUBLAS_STATUS_SUCCESS) {                        \
        const char *msg = NULL;                                   \
        switch (status) {                                         \
        CASE_MSG(CUBLAS_STATUS_SUCCESS);                          \
        CASE_MSG(CUBLAS_STATUS_NOT_INITIALIZED);                  \
        CASE_MSG(CUBLAS_STATUS_ALLOC_FAILED);                     \
        CASE_MSG(CUBLAS_STATUS_INVALID_VALUE);                    \
        CASE_MSG(CUBLAS_STATUS_ARCH_MISMATCH);                    \
        CASE_MSG(CUBLAS_STATUS_MAPPING_ERROR);                    \
        CASE_MSG(CUBLAS_STATUS_EXECUTION_FAILED);                 \
        CASE_MSG(CUBLAS_STATUS_INTERNAL_ERROR);                   \
        CASE_MSG(CUBLAS_STATUS_NOT_SUPPORTED);                    \
        CASE_MSG(CUBLAS_STATUS_LICENSE_ERROR);                    \
        default: break;                                           \
        }                                                         \
        if (msg == NULL) {                                        \
            fprintf(stderr, "%s:%d: cuBLAS error (%s)\n",         \
                    __FILE__, __LINE__, msg);                     \
        }                                                         \
        else {                                                    \
            fprintf(stderr, "%s:%d: cuBLAS error (%d)\n",         \
                    __FILE__, __LINE__, status);                  \
        }                                                         \
    }                                                             \
}

extern "C" {
    void vsmul_cuda(int, const float *, const float *, float *);
    void inner3f_cuda(int, const float *, const float *,
                      const float *, float *);
}

#endif // HAVE_CUDA

/////////////////////////////////////////////////////////////////////

#define CHECK_DLOPEN2(filename, flag)                \
    if (dlopen(filename, flag) == NULL) {            \
        if (dlopen("./" filename, flag) == NULL) {   \
            fprintf(stderr, "%s:%d: %s\n", __FILE__, \
                    __LINE__, dlerror());            \
        }                                            \
    }

#define CHECK_HANDLE(h)                          \
    if ((h) == NULL) {                           \
        fprintf(stderr, "%s:%d: %s\n", __FILE__, \
                __LINE__, dlerror());            \
    }

#define CHECK_DLOPEN3(handle, filename, flag)               \
    if ((handle = dlopen(filename, flag)) == NULL) {        \
        CHECK_HANDLE(handle = dlopen("./" filename, flag)); \
    }

#define SYMLINK_SO(basename)                                 \
    if (access(basename "_so", F_OK) != -1 &&                \
        access(basename ".so", F_OK) == -1) {                \
        if (symlink(basename "_so", basename ".so") == -1) { \
            perror("symlink");                               \
        }                                                    \
    }

/////////////////////////////////////////////////////////////////////

namespace {

    class cblas_t {
    public:
#ifndef _MKL_TYPES_H_
        typedef int MKL_INT;
#endif // _MKL_TYPES_H_
        enum CBLAS_LAYOUT {
            CblasRowMajor  = 101,
            CblasColMajor  = 102
        };
        enum CBLAS_TRANSPOSE {
            CblasNoTrans   = 111,
            CblasTrans     = 112,
            CblasConjTrans = 113
        };
        enum CBLAS_UPLO {
            CblasUpper     = 121,
            CblasLower     = 122
        };
        enum CBLAS_SIDE {
            CblasLeft      = 141,
            CblasRight     = 142
        };
        void (*_get_version_string)(char *, int);
        void (*_set_num_threads_local)(int);
        void (*_set_num_threads)(int);
        int (*_get_max_threads)(void);
        void (*_set_num_stripes)(int);
        int (*_get_num_stripes)(void);
        void (*_domain_set_num_threads)(int);
        int (*_domain_get_max_threads)(void);
        void (*_set_dynamic)(int);
        int (*_get_dynamic)(void);
        float (*_sdot)(const MKL_INT, const float *, const MKL_INT,
                       const float *, const MKL_INT);
        void (*_scopy)(const MKL_INT, const float *, const MKL_INT,
                       float *, const MKL_INT);
        void (*_sscal)(const MKL_INT, const float, float *,
                       const MKL_INT);
        void (*_sgemv)(const CBLAS_LAYOUT, const CBLAS_TRANSPOSE,
                       const MKL_INT, const MKL_INT, const float,
                       const float *, const MKL_INT, const float *,
                       const MKL_INT, const float, float *,
                       const MKL_INT);
        void (*_ssymv)(const CBLAS_LAYOUT, const CBLAS_UPLO,
                       const MKL_INT, const float, const float *,
                       const MKL_INT, const float *, const MKL_INT,
                       const float, float *, const MKL_INT);
        void (*_sger)(const CBLAS_LAYOUT, const MKL_INT,
                      const MKL_INT, const float, const float *,
                      const MKL_INT, const float *, const MKL_INT,
                      float *, const MKL_INT);
        void (*_sgemm)(const CBLAS_LAYOUT, const CBLAS_TRANSPOSE,
                       const CBLAS_TRANSPOSE, const MKL_INT,
                       const MKL_INT, const MKL_INT, const float,
                       const float *, const MKL_INT, const float *,
                       const MKL_INT, const float, float *,
                       const MKL_INT);
        void (*_ssymm)(const CBLAS_LAYOUT, const CBLAS_SIDE,
                       const CBLAS_UPLO, const MKL_INT,
                       const MKL_INT, const float, const float *,
                       const MKL_INT, const float *, const MKL_INT,
                       const float, float *, const MKL_INT);
        void (*_vsmul)(const MKL_INT, const float [], const float [],
                       float []);
#ifdef HAVE_CUDA
        cublasHandle_t _cublas_handle;
        float *_device_const_0;
        float *_device_const_1;
#endif // HAVE_CUDA
        void initialize_mkl(void)
        {
            SYMLINK_SO("libmkl_avx");
            SYMLINK_SO("libmkl_avx2");
            SYMLINK_SO("libmkl_def");
            SYMLINK_SO("libmkl_vml_avx");
            SYMLINK_SO("libmkl_vml_avx2");
            SYMLINK_SO("libmkl_vml_def");

            // The (non-"single") dynamic loading method for Intel
            // Math Kernel Library
            SYMLINK_SO("libiomp5");
            CHECK_DLOPEN2("libiomp5.so", RTLD_NOW | RTLD_GLOBAL);
            // libmkl_core.so has circular dependencies and has to be
            // loaded lazily
            SYMLINK_SO("libmkl_core");
            CHECK_DLOPEN2("libmkl_core.so", RTLD_LAZY | RTLD_GLOBAL);
            SYMLINK_SO("libmkl_intel_thread");
            CHECK_DLOPEN2("libmkl_intel_thread.so",
                          RTLD_NOW | RTLD_GLOBAL);
            SYMLINK_SO("libmkl_core");
            CHECK_DLOPEN2("libmkl_core.so", RTLD_NOW | RTLD_GLOBAL);

            SYMLINK_SO("libmkl_intel_lp64");

            void *mkl_intel_lp64;

            CHECK_DLOPEN3(mkl_intel_lp64, "libmkl_intel_lp64.so",
                          RTLD_NOW | RTLD_GLOBAL);

            _get_version_string =
                reinterpret_cast<void (*)(char *, int)>
                (dlsym(mkl_intel_lp64, "MKL_Get_Version_String"));
            CHECK_HANDLE(_get_version_string);
            _set_num_threads_local =
                reinterpret_cast<void (*)(int)>
                (dlsym(mkl_intel_lp64, "MKL_Set_Num_Threads_Local"));
            CHECK_HANDLE(_set_num_threads_local);
            _set_num_threads =
                reinterpret_cast<void (*)(int)>
                (dlsym(mkl_intel_lp64, "MKL_Set_Num_Threads"));
            CHECK_HANDLE(_set_num_threads);
            _get_max_threads =
                reinterpret_cast<int (*)(void)>
                (dlsym(mkl_intel_lp64, "MKL_Get_Max_Threads"));
            CHECK_HANDLE(_get_max_threads);
            _set_num_stripes =
                reinterpret_cast<void (*)(int)>
                (dlsym(mkl_intel_lp64, "MKL_Set_Num_Stripes"));
            CHECK_HANDLE(_set_num_stripes);
            _get_num_stripes =
                reinterpret_cast<int (*)(void)>
                (dlsym(mkl_intel_lp64, "MKL_Get_Num_Stripes"));
            CHECK_HANDLE(_get_num_stripes);
            _domain_set_num_threads =
                reinterpret_cast<void (*)(int)>
                (dlsym(mkl_intel_lp64, "MKL_Domain_Set_Num_Threads"));
            CHECK_HANDLE(_set_num_threads);
            _domain_get_max_threads =
                reinterpret_cast<int (*)(void)>
                (dlsym(mkl_intel_lp64, "MKL_Domain_Get_Max_Threads"));
            CHECK_HANDLE(_domain_get_max_threads);
            _set_dynamic =
                reinterpret_cast<void (*)(int)>
                (dlsym(mkl_intel_lp64, "MKL_Set_Dynamic"));
            CHECK_HANDLE(_set_dynamic);
            _get_dynamic =
                reinterpret_cast<int (*)(void)>
                (dlsym(mkl_intel_lp64, "MKL_Get_Dynamic"));
            CHECK_HANDLE(_get_dynamic);

            _sdot =
                reinterpret_cast<
                float (*)(const MKL_INT, const float *,
                          const MKL_INT, const float *,
                          const MKL_INT)>
                (dlsym(mkl_intel_lp64, "cblas_sdot"));
            CHECK_HANDLE(_sdot);
            _scopy =
                reinterpret_cast<
                void (*)(const MKL_INT, const float *,
                         const MKL_INT, float *, const MKL_INT)>
                (dlsym(mkl_intel_lp64, "cblas_scopy"));
            CHECK_HANDLE(_scopy);
            _sscal =
                reinterpret_cast<
                void (*)(const MKL_INT, const float, float *,
                         const MKL_INT)>
                (dlsym(mkl_intel_lp64, "cblas_sscal"));
            CHECK_HANDLE(_sscal);
            _sgemv =
                reinterpret_cast<
                void (*)(const CBLAS_LAYOUT, const CBLAS_TRANSPOSE,
                         const MKL_INT, const MKL_INT, const float,
                         const float *, const MKL_INT, const float *,
                         const MKL_INT, const float, float *,
                         const MKL_INT)>
                (dlsym(mkl_intel_lp64, "cblas_sgemv"));
            CHECK_HANDLE(_sgemv);
            _ssymv =
                reinterpret_cast<
                void (*)(const CBLAS_LAYOUT, const CBLAS_UPLO,
                         const MKL_INT, const float, const float *,
                         const MKL_INT, const float *, const MKL_INT,
                         const float, float *, const MKL_INT)>
                (dlsym(mkl_intel_lp64, "cblas_ssymv"));
            CHECK_HANDLE(_ssymv);
            _sger =
                reinterpret_cast<
                void (*)(const CBLAS_LAYOUT, const MKL_INT,
                         const MKL_INT, const float, const float *,
                         const MKL_INT, const float *, const MKL_INT,
                         float *, const MKL_INT)>
                (dlsym(mkl_intel_lp64, "cblas_sger"));
            CHECK_HANDLE(_sger);
            _sgemm =
                reinterpret_cast<
                void (*)(const CBLAS_LAYOUT, const CBLAS_TRANSPOSE,
                         const CBLAS_TRANSPOSE, const MKL_INT,
                         const MKL_INT, const MKL_INT, const float,
                         const float *, const MKL_INT, const float *,
                         const MKL_INT, const float, float *,
                         const MKL_INT)>
                (dlsym(mkl_intel_lp64, "cblas_sgemm"));
            CHECK_HANDLE(_sgemm);
            _ssymm =
                reinterpret_cast<
                void (*)(const CBLAS_LAYOUT, const CBLAS_SIDE,
                         const CBLAS_UPLO, const MKL_INT,
                         const MKL_INT, const float, const float *,
                         const MKL_INT, const float *, const MKL_INT,
                         const float, float *, const MKL_INT)>
                (dlsym(mkl_intel_lp64, "cblas_ssymm"));
            CHECK_HANDLE(_ssymm);
            _vsmul =
                reinterpret_cast<
                void (*)(const MKL_INT, const float [],
                         const float [], float [])>
                (dlsym(mkl_intel_lp64, "vsMul"));
            CHECK_HANDLE(_vsmul);
        }
#ifdef HAVE_CUDA
        void initialize_cuda(void)
        {
            CHECK_CUBLAS(cublasCreate(&_cublas_handle));
            CHECK_CUBLAS(cublasSetPointerMode(
                _cublas_handle, CUBLAS_POINTER_MODE_DEVICE));
            CHECK_CUDA(cudaMalloc(reinterpret_cast<void **>
                                  (&_device_const_0),
                                  sizeof(float)));
            CHECK_CUDA(cudaMalloc(reinterpret_cast<void **>
                                  (&_device_const_1),
                                  sizeof(float)));

            const float host_const_0 = 0;
            const float host_const_1 = 1;

            CHECK_CUBLAS(cublasSetVector(1, sizeof(float),
                                         &host_const_0, 1,
                                         _device_const_0, 1));
            CHECK_CUBLAS(cublasSetVector(1, sizeof(float),
                                         &host_const_1, 1,
                                         _device_const_1, 1));
        }
#endif // HAVE_CUDA
        cblas_t(void)
            : _get_version_string(NULL),
              _set_num_threads_local(NULL), _set_num_threads(NULL),
              _get_max_threads(NULL), _set_num_stripes(NULL),
              _get_num_stripes(NULL), _domain_set_num_threads(NULL),
              _domain_get_max_threads(NULL),
              _set_dynamic(NULL), _get_dynamic(NULL),
              _sdot(NULL), _scopy(NULL), _sscal(NULL), _sgemv(NULL),
              _ssymv(NULL), _sger(NULL), _sgemm(NULL), _ssymm(NULL),
              _vsmul(NULL)
        {
#ifdef HAVE_CUDA
            initialize_cuda();
#else // HAVE_CUDA
            initialize_mkl();
#endif // HAVE_CUDA
        }
        bool is_cublas(void) const
        {
#ifdef HAVE_CUDA
            return _cublas_handle != NULL;
#else // HAVE_CUDA
            return false;
#endif // HAVE_CUDA
        }
        std::string version_str(void)
        {
            // The maximum buffer size MKL_Get_Version_String can
            // handle appears to be 4k, larger buffer will cause
            // MKL_Get_Version_String to silently return (and leaving
            // an empty string)
            char buf[4096] = { '\0' };

            if (_get_version_string != NULL) {
                (*_get_version_string)(buf, 4096);
            }
            return std::string(buf);
        }
        void sdot(size_t n, const float *x, long int incx,
                  const float *y, long int incy, float *result)
        {
#ifdef HAVE_CUDA
            CHECK_CUBLAS(cublasSdot(_cublas_handle,
                                    n, x, incx, y, incy, result));
#else // HAVE_CUDA
            *result = _sdot(n, x, incx, y, incy);
#endif // HAVE_CUDA
        }
        void scopy(size_t n, const float *x, float *y)
        {
#ifdef HAVE_CUDA
            CHECK_CUBLAS(cublasScopy(_cublas_handle, n, x, 1, y, 1));
#else // HAVE_CUDA
            _scopy(n, x, 1, y, 1);
#endif // HAVE_CUDA
        }
        void sscal(size_t n, const float *sa, float *sx, long int incx)
        {
#ifdef HAVE_CUDA
            // This code path should not happen, vs. cublasSdgmm
            CHECK_CUBLAS(cublasSscal(_cublas_handle,
                                     n, sa, sx, incx));
#else // HAVE_CUDA
            _sscal(n, *sa, sx, incx);
#endif // HAVE_CUDA
        }
        void sgemv(bool nota, size_t n, const float *a,
                   const float *x, float *y)
        {
#ifdef HAVE_CUDA
            CHECK_CUBLAS(cublasSgemv(_cublas_handle,
                                     nota ?
                                     CUBLAS_OP_N : CUBLAS_OP_T,
                                     n, n, _device_const_1, a, n,
                                     x, 1, _device_const_0, y, 1));
#else // HAVE_CUDA
            _sgemv(cblas_t::CblasColMajor,
                   nota ?
                   cblas_t::CblasNoTrans :
                   cblas_t::CblasTrans,
                   n, n, 1.0F, a, n, x, 1, 0.0F, y, 1);
#endif // HAVE_CUDA
        }
        void ssymv(size_t n, const float *a, const float *x,
                   float *y)
        {
#ifdef HAVE_CUDA
            CHECK_CUBLAS(cublasSsymv(_cublas_handle,
                                     CUBLAS_FILL_MODE_UPPER,
                                     n, _device_const_1, a, n,
                                     x, 1, _device_const_0, y, 1));
#else // HAVE_CUDA
            _ssymv(cblas_t::CblasColMajor,
                   cblas_t::CblasUpper,
                   n, 1.0F, a, n, x, 1, 0.0F, y, 1);
#endif // HAVE_CUDA
        }
        void sger(size_t n, const float *x, const float *y,
                  float *a)
        {
#ifdef HAVE_CUDA
            CHECK_CUBLAS(cublasSger(_cublas_handle,
                                    n, n, _device_const_1, x, 1,
                                    y, 1, a, n));
#else // HAVE_CUDA
            _sger(cblas_t::CblasColMajor,
                  n, n, 1.0F, x, 1, y, 1, a, n);
#endif // HAVE_CUDA
        }
        void sgemm(bool transa, bool transb, size_t n,
                   const float *a, const float *b, float *c)
        {
#ifdef HAVE_CUDA
            CHECK_CUBLAS(cublasSgemm(_cublas_handle,
                                     transa ?
                                     CUBLAS_OP_T : CUBLAS_OP_N,
                                     transb ?
                                     CUBLAS_OP_T : CUBLAS_OP_N,
                                     n, n, n, _device_const_1, a, n,
                                     b, n, _device_const_0, c, n));
#else // HAVE_CUDA
            _sgemm(cblas_t::CblasColMajor,
                   transa ?
                   cblas_t::CblasTrans : cblas_t::CblasNoTrans,
                   transb ?
                   cblas_t::CblasTrans : cblas_t::CblasNoTrans,
                   n, n, n, 1.0F, a, n, b, n, 0.0F, c, n);
#endif // HAVE_CUDA
        }
        void ssymm(bool sider, size_t n, const float *a,
                   const float *b, float *c)
        {
#ifdef HAVE_CUDA
            CHECK_CUBLAS(cublasSsymm(_cublas_handle,
                                     sider ?
                                     CUBLAS_SIDE_RIGHT :
                                     CUBLAS_SIDE_LEFT,
                                     CUBLAS_FILL_MODE_UPPER,
                                     n, n, _device_const_1, a, n,
                                     b, n, _device_const_0, c, n));
#else // HAVE_CUDA
            _ssymm(cblas_t::CblasColMajor,
                   sider ? cblas_t::CblasRight : cblas_t::CblasLeft,
                   cblas_t::CblasUpper,
                   n, n, 1.0F, a, n, b, n, 0.0F, c, n);
#endif // HAVE_CUDA
        }
        void sdgmm(bool side_right, size_t n,
                   const float *a, const float *x, float *c)
        {
#ifdef HAVE_CUDA
            CHECK_CUBLAS(cublasSdgmm(_cublas_handle,
                                     side_right ?
                                     CUBLAS_SIDE_RIGHT :
                                     CUBLAS_SIDE_LEFT,
                                     n, n, a, n, x, 1, c, n));
#else // HAVE_CUDA
            scopy(std::pow(n, 2), a, c);
            if (side_right) {
                // a,ab->ab
                for (size_t i = 0; i < n; i++) {
                    sscal(n, x + i, c + i * n, 1);
                }
            }
            else {
                // b,ab->ab
                for (size_t i = 0; i < n; i++) {
                    sscal(n, x + i, c + i, n);
                }
            }
#endif // HAVE_CUDA
        }
        void vsmul(size_t n, const float a[], const float b[],
                   float r[])
        {
#ifdef HAVE_CUDA
            vsmul_cuda(n, a, b, r);
#else // HAVE_CUDA
            _vsmul(n, a, b, r);
#endif // HAVE_CUDA
        }

        /////////////////////////////////////////////////////////////

        bool test_avx(bool *have_fma3 = NULL,
                      bool *have_avx512f = NULL)
        {
            unsigned int eax = 1;
            unsigned int cpuid[4];

            __asm__ __volatile__ (
                "xchg   %%ebx, %%edi\n\t"
                "cpuid\n\t"
                "xchg   %%ebx, %%edi\n"
                : "=a" (cpuid[0]), "=D" (cpuid[1]), "=c" (cpuid[2]),
                  "=d" (cpuid[3])
                : "0" (eax));

            unsigned int ecx = 0;
            unsigned int edx;

            __asm__ __volatile__ (
                "xgetbv\n"
                : "=a" (eax), "=d"(edx)
                : "c" (ecx));

            static const unsigned int mask_osxsave_avx =
                (1 << 27) | (1 << 28);
            static const unsigned int mask_ymm_state =
                (1 << 1) | (1 << 2);

            const bool have_avx =
                (cpuid[2] & mask_osxsave_avx) == mask_osxsave_avx &&
                (eax & mask_ymm_state) == mask_ymm_state;

            if (have_fma3 != NULL) {
                static const unsigned int mask_fma3 = 1 << 12;

                *have_fma3 = (cpuid[2] & mask_fma3) == mask_fma3;
            }

            if (have_avx512f != NULL) {
                eax = 7;
                unsigned int ecx = 0;

                __asm__ __volatile__ (
                    "xchg   %%ebx, %%edi\n\t"
                    "cpuid\n\t"
                    "xchg   %%ebx, %%edi\n"
                    : "=a" (cpuid[0]), "=D" (cpuid[1]),
                      "=c" (cpuid[2]), "=d" (cpuid[3])
                    : "0" (eax), "2" (ecx));

                static const unsigned int mask_avx512f = 1 << 16;

                *have_avx512f =
                    (cpuid[1] & mask_avx512f) == mask_avx512f;
            }

            return have_avx;
        }

#define INNER3F_INIT_SSE                        \
    "xorps          %%xmm4, %%xmm4\n\t"         \
    "xorps          %%xmm5, %%xmm5\n\t"         \
    "xorps          %%xmm6, %%xmm6\n\t"         \
    "xorps          %%xmm7, %%xmm7\n"

#define INNER3F_INIT_AVX                        \
    "vxorps         %%xmm4, %%xmm4, %%xmm4\n\t" \
    "vxorps         %%xmm5, %%xmm5, %%xmm5\n\t" \
    "vxorps         %%xmm6, %%xmm6, %%xmm6\n\t" \
    "vxorps         %%xmm7, %%xmm7, %%xmm7\n"

#define INNER3F_BEGIN_LOOP                      \
    "1:\n\t"

#define INNER3F_LOAD_MUL_SSE                    \
    "movups         (%4,%0,4), %%xmm8\n\t"      \
    "movups         (%5,%0,4), %%xmm12\n\t"     \
    "movups         16(%4,%0,4), %%xmm9\n\t"    \
    "movups         16(%5,%0,4), %%xmm13\n\t"   \
    "movups         32(%4,%0,4), %%xmm10\n\t"   \
    "movups         32(%5,%0,4), %%xmm14\n\t"   \
    "movups         48(%4,%0,4), %%xmm11\n\t"   \
    "movups         48(%5,%0,4), %%xmm15\n\t"   \
    "mulps          %%xmm12, %%xmm8\n\t"        \
    "mulps          %%xmm13, %%xmm9\n\t"        \
    "mulps          %%xmm14, %%xmm10\n\t"       \
    "mulps          %%xmm15, %%xmm11\n\t"

#define INNER3F_LOAD_MUL_AVX                       \
    "vmovups        (%4,%0,4), %%ymm8\n\t"         \
    "vmovups        (%5,%0,4), %%ymm12\n\t"        \
    "vmovups        32(%4,%0,4), %%ymm9\n\t"       \
    "vmovups        32(%5,%0,4), %%ymm13\n\t"      \
    "vmovups        64(%4,%0,4), %%ymm10\n\t"      \
    "vmovups        64(%5,%0,4), %%ymm14\n\t"      \
    "vmovups        96(%4,%0,4), %%ymm11\n\t"      \
    "vmovups        96(%5,%0,4), %%ymm15\n\t"      \
    "vmulps         %%ymm12, %%ymm8, %%ymm8\n\t"   \
    "vmulps         %%ymm13, %%ymm9, %%ymm9\n\t"   \
    "vmulps         %%ymm14, %%ymm10, %%ymm10\n\t" \
    "vmulps         %%ymm15, %%ymm11, %%ymm11\n\t"

#define INNER3F_LOAD_MUL_AVX512F                   \
    "vmovups        (%4,%0,4), %%zmm8\n\t"         \
    "vmovups        (%5,%0,4), %%zmm12\n\t"        \
    "vmovups        64(%4,%0,4), %%zmm9\n\t"       \
    "vmovups        64(%5,%0,4), %%zmm13\n\t"      \
    "vmovups        128(%4,%0,4), %%zmm10\n\t"     \
    "vmovups        128(%5,%0,4), %%zmm14\n\t"     \
    "vmovups        192(%4,%0,4), %%zmm11\n\t"     \
    "vmovups        192(%5,%0,4), %%zmm15\n\t"     \
    "vmulps         %%zmm12, %%zmm8, %%zmm8\n\t"   \
    "vmulps         %%zmm13, %%zmm9, %%zmm9\n\t"   \
    "vmulps         %%zmm14, %%zmm10, %%zmm10\n\t" \
    "vmulps         %%zmm15, %%zmm11, %%zmm11\n\t"

#define INNER3F_MADD_SSE                        \
    "mulps          (%6,%0,4), %%xmm8\n\t"      \
    "mulps          16(%6,%0,4), %%xmm9\n\t"    \
    "mulps          32(%6,%0,4), %%xmm10\n\t"   \
    "mulps          48(%6,%0,4), %%xmm11\n\t"   \
    "addps          %%xmm8, %%xmm4\n\t"         \
    "addps          %%xmm9, %%xmm5\n\t"         \
    "addps          %%xmm10, %%xmm6\n\t"        \
    "addps          %%xmm11, %%xmm7\n\t"

#define INNER3F_MADD_AVX                                \
    "vmulps         (%6,%0,4), %%ymm8, %%ymm8\n\t"      \
    "vmulps         32(%6,%0,4), %%ymm9, %%ymm9\n\t"    \
    "vmulps         64(%6,%0,4), %%ymm10, %%ymm10\n\t"  \
    "vmulps         96(%6,%0,4), %%ymm11, %%ymm11\n\t"  \
    "vaddps         %%ymm8, %%ymm4, %%ymm4\n\t"         \
    "vaddps         %%ymm9, %%ymm5, %%ymm5\n\t"         \
    "vaddps         %%ymm10, %%ymm6, %%ymm6\n\t"        \
    "vaddps         %%ymm11, %%ymm7, %%ymm7\n\t"

#define INNER3F_MADD_AVX_FMA3                           \
    "vfmadd231ps    (%6,%0,4), %%ymm8, %%ymm4\n\t"      \
    "vfmadd231ps    32(%6,%0,4), %%ymm9, %%ymm5\n\t"    \
    "vfmadd231ps    64(%6,%0,4), %%ymm10, %%ymm6\n\t"   \
    "vfmadd231ps    96(%6,%0,4), %%ymm11, %%ymm7\n\t"

#define INNER3F_MADD_AVX512F                            \
    "vfmadd231ps    (%6,%0,4), %%zmm8, %%zmm4\n\t"      \
    "vfmadd231ps    64(%6,%0,4), %%zmm9, %%zmm5\n\t"    \
    "vfmadd231ps    128(%6,%0,4), %%zmm10, %%zmm6\n\t"  \
    "vfmadd231ps    192(%6,%0,4), %%zmm11, %%zmm7\n\t"

#define INNER3F_END_LOOP_SSE                    \
    "addq           $16, %0\n\t"                \
    "subq           $16, %1\n\t"                \
    "jnz            1b\n\t"

#define INNER3F_END_LOOP_AVX                    \
    "addq           $32, %0\n\t"                \
    "subq           $32, %1\n\t"                \
    "jnz            1b\n\t"

#define INNER3F_END_LOOP_AVX512F                \
    "addq           $64, %0\n\t"                \
    "subq           $64, %1\n\t"                \
    "jnz            1b\n\t"

#define INNER3F_HADD_STORE_SSE                  \
    "addps          %%xmm5, %%xmm4\n\t"         \
    "addps          %%xmm7, %%xmm6\n\t"         \
    "addps          %%xmm6, %%xmm4\n\t"         \
    "haddps         %%xmm4, %%xmm4\n\t"         \
    "haddps         %%xmm4, %%xmm4\n\t"         \
    "movss          %%xmm4, (%7)"

#define INNER3F_HADD_STORE_AVX                      \
    "vaddps         %%ymm5, %%ymm4, %%ymm4\n\t"     \
    "vaddps         %%ymm7, %%ymm6, %%ymm6\n\t"     \
    "vaddps         %%ymm6, %%ymm4, %%ymm4\n\t"     \
    "vextractf128   $1, %%ymm4, %%xmm5\n\t"         \
    "vaddps         %%xmm5, %%xmm4, %%xmm4\n\t"     \
    "vhaddps        %%xmm4, %%xmm4, %%xmm4\n\t"     \
    "vhaddps        %%xmm4, %%xmm4, %%xmm4\n\t"     \
    "vmovss         %%xmm4, (%7)\n\t"               \
    "vzeroupper"

#define INNER3F_HADD_STORE_AVX512F                  \
    "vaddps         %%zmm5, %%zmm4, %%zmm4\n\t"     \
    "vaddps         %%zmm7, %%zmm6, %%zmm6\n\t"     \
    "vaddps         %%zmm6, %%zmm4, %%zmm4\n\t"     \
    "vextractf64x4  $1, %%zmm4, %%ymm5\n\t"         \
    "vaddps         %%ymm5, %%ymm4, %%ymm4\n\t"     \
    "vextractf128   $1, %%ymm4, %%xmm5\n\t"         \
    "vaddps         %%xmm5, %%xmm4, %%xmm4\n\t"     \
    "vhaddps        %%xmm4, %%xmm4, %%xmm4\n\t"     \
    "vhaddps        %%xmm4, %%xmm4, %%xmm4\n\t"     \
    "vmovss         %%xmm4, (%7)"

#define INNER3F_OPERAND_CLOBBER                                     \
    : "=r" (i), "=r" (n)                                            \
    : "0" (i), "1" (n), "r" (a), "r" (b), "r" (c), "r" (dot)        \
    : "cc", "%xmm4", "%xmm5", "%xmm6", "%xmm7", "%xmm8", "%xmm9",   \
      "%xmm10", "%xmm11", "%xmm12", "%xmm13", "%xmm14", "%xmm15",   \
      "memory"

        void inner3f_sse_16(size_t n, const float *a, const float *b,
                            const float *c, float *dot)
        {
            register size_t i = 0;

            __asm__ __volatile__ (
                INNER3F_INIT_SSE
                INNER3F_BEGIN_LOOP
                INNER3F_LOAD_MUL_SSE
                INNER3F_MADD_SSE
                INNER3F_END_LOOP_SSE
                INNER3F_HADD_STORE_SSE
                INNER3F_OPERAND_CLOBBER);
        }

        void inner3f_avx_32(size_t n, const float *a, const float *b,
                            const float *c, float *dot)
        {
            register size_t i = 0;

            __asm__ __volatile__ (
                INNER3F_INIT_AVX
                INNER3F_BEGIN_LOOP
                INNER3F_LOAD_MUL_AVX
                INNER3F_MADD_AVX
                INNER3F_END_LOOP_AVX
                INNER3F_HADD_STORE_AVX
                INNER3F_OPERAND_CLOBBER);
        }

        void inner3f_avx_fma3_32(size_t n, const float *a,
                                 const float *b, const float *c,
                                 float *dot)
        {
            register size_t i = 0;

            __asm__ __volatile__ (
                INNER3F_INIT_AVX
                INNER3F_BEGIN_LOOP
                INNER3F_LOAD_MUL_AVX
                INNER3F_MADD_AVX_FMA3
                INNER3F_END_LOOP_AVX
                INNER3F_HADD_STORE_AVX
                INNER3F_OPERAND_CLOBBER);
        }

        void inner3f_avx512f_64(size_t n, const float *a,
                                const float *b, const float *c,
                                float *dot)
        {
            register size_t i = 0;

            __asm__ __volatile__ (
                INNER3F_INIT_AVX
                INNER3F_BEGIN_LOOP
                INNER3F_LOAD_MUL_AVX512F
                INNER3F_MADD_AVX512F
                INNER3F_END_LOOP_AVX512F
                INNER3F_HADD_STORE_AVX512F
                INNER3F_OPERAND_CLOBBER);
        }

        void inner3f(size_t n, const float *a, const float *b,
                     const float *c, float *dot)
        {
#ifdef HAVE_CUDA
            inner3f_cuda(n, a, b, c, dot);
#else // HAVE_CUDA
            static bool have_fma3;
            static bool have_avx512f;
            static bool have_avx =
                test_avx(&have_fma3, &have_avx512f);
            size_t i;

            *dot = 0;
            if (have_avx512f) {
                size_t n1 = n & (size_t)(-64);

                if (n1 != 0) {
                    inner3f_avx512f_64(n1, a, b, c, dot);
                }
                i = n1;
                n1 = (n - n1) & (size_t)(-32);

                if (n1 != 0) {
                    float dot_32 = 0;

                    inner3f_avx_fma3_32(n1, a + i, b + i, c + i,
                                        &dot_32);
                    *dot += dot_32;
                }
                i += n1;
            }
            else if (have_avx) {
                const size_t n1 = n & (size_t)(-32);

                if (n1 != 0) {
                    if (have_fma3) {
                        inner3f_avx_fma3_32(n1, a, b, c, dot);
                    }
                    else {
                        inner3f_avx_32(n1, a, b, c, dot);
                    }
                }
                i = n1;
            }
            else {
                size_t nu = (4 - ((size_t)c >> 2)) & 0x3;

                for (i = 0; i < nu; i++) {
                    *dot += a[i] * b[i] * c[i];
                }

                const size_t n1 = (n - nu) & (size_t)(-16);
                float dot_a;

                if (n1 != 0) {
                    inner3f_sse_16(n1, a + nu, b + nu, c + nu,
                                   &dot_a);
                    *dot += dot_a;
                }
                i = n1;
            }

            for (; i < n; i++) {
                *dot += a[i] * b[i] * c[i];
            }
#endif // HAVE_CUDA
        }
    };
}
