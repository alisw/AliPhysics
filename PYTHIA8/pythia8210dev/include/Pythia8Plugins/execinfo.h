// execinfo.h
//
// Dummy routines that implement execinfo, see
// http://www.gnu.org/software/libc/manual/html_node/Backtraces.html
// This file has to be in the include path to compile Pythia under Cygwin,
// but otherwise has no function and can be removed.

int backtrace (void **buffer, int size) { return 0;}

char ** backtrace_symbols (void *const *buffer, int size) { return 0;}

void backtrace_symbols_fd (void *const *buffer, int size, int fd) {}
