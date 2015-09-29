#!/usr/bin/python
import zmq,sys

mode="SUB"
endpoint=">tcp://localhost:60202"
requestHeader="INFO"
requestBody=""
if len(sys.argv)==1:
    print "Usage:"
    print sys.argv[0]+" PULL \"@tcp://localhost:60202\""
    print sys.argv[0]+" SUB \">tcp://localhost:60202\" \"optional subscription\""
    print sys.argv[0]+" REQ \">tcp://localhost:60202\" INFO \"-reset someoption=value\""
    print sys.argv[0]+" REP \">tcp://localhost:60202\" INFO \"-reset someoption=value\""
    quit()
if len(sys.argv)>1:
    mode=sys.argv[1]
if len(sys.argv)>2:
    endpoint=sys.argv[2]
if len(sys.argv)>3:
    requestHeader=sys.argv[3]
if len(sys.argv)>4:
    requestBody=sys.argv[4]

print "mode: "+mode+" endpoint: "+endpoint

#  Prepare our context and sockets
context = zmq.Context()

if mode=="PUSH":
    socket = context.socket(zmq.PUSH)
elif mode=="PULL":
    socket = context.socket(zmq.PULL)
elif mode=="SUB":
    socket = context.socket(zmq.SUB)
    socket.setsockopt(zmq.SUBSCRIBE, requestHeader)
elif mode=="PUB":
    socket = context.socket(zmq.PUB)
elif mode=="REQ":
    socket = context.socket(zmq.REQ)
elif mode=="REP":
    socket = context.socket(zmq.REP)
else:
    print "not a valid mode"
    quit()

if endpoint[0]=='>':
    print "connect to: "+endpoint[1:]
    socket.connect(str(endpoint[1:]))
elif endpoint[0]=='@':
    print "bind to: "+endpoint[1:]
    socket.bind(str(endpoint[1:]))

endless=True
if mode=="REQ" or mode=="PUSH" or mode=="PUB":
    endless=False
    socket.send(requestHeader,zmq.SNDMORE)
    socket.send(requestBody)
    print("sent: "+requestHeader+" "+requestBody)

while True:
    if mode=="SUB" or mode=="PULL" or mode=="REP" or mode=="REQ":
        msg = socket.recv_multipart();
        print "_________________________"
        i=0;
        for message in msg:
            if i==0:
                print "topic: "+str(message)
                i=1
            elif i==1:
                print "  messagesize: "+str(sys.getsizeof(message))
                i=0
    if mode=="REP":
        socket.send(requestHeader,zmq.SNDMORE)
        socket.send(requestBody)
        print "  sent back: "+requestHeader+" "+requestBody
    if not endless:
        break
