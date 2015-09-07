
// blame: Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch
// some of it is more of less based on czmq.h

//  --------------------------------------------------------------------------
//  Attach a socket to zero or more endpoints. If endpoints is not null,
//  parses as list of ZeroMQ endpoints, separated by commas, and prefixed by
//  '@' (to bind the socket) or '>' (to attach the socket). Returns 0 if all
//  endpoints were valid, or -1 if there was a syntax error. If the endpoint
//  does not start with '@' or '>', the serverish argument defines whether
//  it is used to bind (serverish = true) or connect (serverish = false).

int alizmq_attach (void *self, const char *endpoints, bool serverish=false);

