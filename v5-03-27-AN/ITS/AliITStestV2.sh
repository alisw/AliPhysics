#!/bin/sh


# run TPC
aliroot -q -b "$ALICE_ROOT/TPC/AliTPCtest.C"

# run ITS
aliroot -q -b "$ALICE_ROOT/ITS/AliITStestV2.C"


