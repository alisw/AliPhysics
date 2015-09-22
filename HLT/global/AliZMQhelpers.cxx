#ifndef __AliZMQhelpers__
#define __AliZMQhelpers__

// blame: Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch
//see header file for details
//

#include "AliZMQhelpers.h"

#include "zmq.h"
#include <cstring>
#include <cassert>

int alizmq_attach (void *self, const char *endpoints, bool serverish)
{
    assert (self);
    if (!endpoints)
        return 0;

    //  We hold each individual endpoint here
    char endpoint [256];
    while (*endpoints) {
        const char *delimiter = strchr (endpoints, ',');
        if (!delimiter)
            delimiter = endpoints + strlen (endpoints);
        if (delimiter - endpoints > 255)
            return -1;
        memcpy (endpoint, endpoints, delimiter - endpoints);
        endpoint [delimiter - endpoints] = 0;

        int rc;
        if (endpoint [0] == '@')
            rc = zmq_bind (self, endpoint + 1);
        else
        if (endpoint [0] == '>')
            rc = zmq_connect (self, endpoint + 1);
        else
        if (serverish)
            rc = zmq_bind (self, endpoint);
        else
            rc = zmq_connect (self, endpoint);

        if (rc == -1)
            return -1;          //  Bad endpoint syntax

        if (*delimiter == 0)
            break;
        endpoints = delimiter + 1;
    }
    return 0;
}

int alizmq_msg_recv(std::map<std::string,std::string>& message, void* socket, int flags)
{
  int rc = 0;
  while (true)
  {
    zmq_msg_t requestTopicMsg;
    rc = zmq_msg_init(&requestTopicMsg);
    rc = zmq_msg_recv(&requestTopicMsg, socket, 0);
    if (!zmq_msg_more(&requestTopicMsg)) break;

    zmq_msg_t requestMsg;
    rc = zmq_msg_init(&requestMsg);
    rc = zmq_msg_recv(&requestMsg, socket, 0);
    if (!zmq_msg_more(&requestMsg)) break;

    std::string requestBody;
    std::string requestTopic;
    requestTopic.assign(static_cast<char*>(zmq_msg_data(&requestTopicMsg)), zmq_msg_size(&requestTopicMsg));
    requestBody.assign(static_cast<char*>(zmq_msg_data(&requestMsg)), zmq_msg_size(&requestMsg));

    message[requestTopic] = requestBody;

    zmq_msg_close(&requestTopicMsg);
    zmq_msg_close(&requestMsg);
  }
  return 0;
}
#endif
