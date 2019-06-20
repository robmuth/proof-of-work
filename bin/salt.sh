#!/bin/bash

# Read 32 bytes (256 bits) from /dev/urandom as hex
# and capitalizing letters in hex
dd if=/dev/urandom bs=1 count=32 2>/dev/null | xxd -p -c 64 | tr a-z A-Z

