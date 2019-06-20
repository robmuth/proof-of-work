#!/usr/bin/env python3

import sys
import hashlib
import random
import time
from bitstring import BitArray

from zokrates.eddsa import PrivateKey, PublicKey, Point
from zokrates.field import FQ
from zokrates.utils import to_bytes

def terminalSign(employeeHash,ts=time.time()):
    """Add timestamp and sign employee message"""
    # Prepare mst to sign
    timestamp = int(ts)
    timestampBin = to_bytes(timestamp)
    terminalMsg = hashlib.sha256(employeeHash + timestampBin).digest()

    # Pad with 256bit since implemented eddsa expects a 512bit message
    terminalMsg = terminalMsg + to_bytes(int(0))

    # sk = PrivateKey.from_rand()
    # Seeded for debug purpose
    key = FQ(1997011358982923168928344992199991480689546837621580239342656433234255379025)
    sk = PrivateKey(key)
    sig = sk.sign(terminalMsg)

    pk = PublicKey.from_private(sk)

    return {'pk': pk, 'sig': sig, 'msg': terminalMsg, 'timestamp': timestamp}

def write_input_for_zokrates_cli(pk, sig, msg, timestamp):
    "Writes the input arguments for ZoKrates to file."
    sig_R, sig_S = sig
    args = [sig_R.x, sig_R.y, sig_S, pk.p.x.n, pk.p.y.n]
    args = " ".join(map(str, args))

    M0 = msg.hex()[:64]
    M1 = msg.hex()[64:]
    b0 = BitArray(to_bytes(int(M0, 16))).bin
    b1 = BitArray(to_bytes(int(M1, 16))).bin
    args = args + " " + " ".join(b0 + b1)

    args = args + " " + " ".join(map(str,[timestamp]))

    print(args)

if __name__ == "__main__":
    employeeHash = int('0x{}'.format(sys.argv[1]), 16)

    terminalSignature = terminalSign(to_bytes(employeeHash))

    write_input_for_zokrates_cli(terminalSignature["pk"], terminalSignature["sig"], terminalSignature["msg"], terminalSignature["timestamp"])
