ZMQESDserver reads AliESDs.root files and sends the (selected) data in a columnar format via a ZeroMQ socket. The data consists of a multi-part O2 data block message, each payload contains an array filled with a single variable, i.e. they are arrays if either int or float.

Any O2 device can connect to the zeromq socket and access the data.

Can work in push mode, or as a more traditional server replying to requests.
examples:

push mode, pushing to a PULL socket on localhost:2323
```
ZMQESDserver in="$(ls */AliESDs.root)" out=PUSH+tcp://localhost:2323
```

one can check how the data looks like by starting a monitoring tool:
```
ZMQmonitor.py PULL@tcp://*:2323
```

server mode listening on tcp port 2324 for connections from anywhere:
```
ZMQESDserver in="$(ls */AliESDs.root)" out=REP@tcp://*:2324
```

... and the monitoring of the client side:
```
ZMQmonitor.py REQ+tcp://localhost:2324
```

run ```ZMQESDserver``` (without args) to see the list of options.

