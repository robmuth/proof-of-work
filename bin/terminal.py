import hashlib
import random
import time
from bitstring import BitArray

from zokrates.eddsa import PrivateKey, PublicKey, Point
from zokrates.field import FQ

def generateNonce(length=8):
    """Generate pseudorandom number."""
    return ''.join([str(random.randint(0, 9)) for i in range(length)])

def generateEmployeeHash(eid):
    """Generate employee hash as concat of eid and a nonce."""
    eidBin = int(eid).to_bytes(32, "big")
    nonce = generateNonce()
    nonceBin = int(nonce).to_bytes(32, "big")

    employeeHash = hashlib.sha256(eidBin + nonceBin).digest()
    return {'eid': eid, 'nonce': nonce ,'employeeHash': employeeHash}




def terminalSign(employeeHash,ts=time.time()):
    """Add timestamp and sign employee message"""
    # Prepare mst to sign
    timestamp = int(ts)
    timestampBin = timestamp.to_bytes(32,"big")
    terminalMsg = hashlib.sha256(employeeHash + timestampBin).digest()

    # Pad with 256bit since implemented eddsa expects a 512bit message
    terminalMsg = terminalMsg + int(0).to_bytes(32, "big")

    # sk = PrivateKey.from_rand()
    # Seeded for debug purpose
    key = FQ(1997011358982923168928344992199991480689546837621580239342656433234255379025)
    sk = PrivateKey(key)
    sig = sk.sign(terminalMsg)

    pk = PublicKey.from_private(sk)

    return {'pk': pk, 'sig': sig, 'msg': terminalMsg, 'timestamp': timestamp}

def to_bytes(*args):
    "Returns byte representation for objects used in this module."
    result = b""
    for M in args:
        if isinstance(M, Point):
            result += to_bytes(M.x)
            # result += to_bytes(M.y)
        elif isinstance(M, FQ):
            result += to_bytes(M.n)
        elif isinstance(M, int):
            result += M.to_bytes(32, "big")
        elif isinstance(M, BitArray):
            result += M.tobytes()
        elif isinstance(M, bytes):
            result += M
        elif isinstance(M, (list, tuple)):
            result += b"".join(to_bytes(_) for _ in M)
        else:
            raise TypeError("Bad type for M: " + str(type(M)))
    return result


def write_input_for_zokrates_cli(pk, sig, msg, path, timestamp, eid, nonce):
    "Writes the input arguments for ZoKrates to file."
    sig_R, sig_S = sig
    args = [sig_R.x, sig_R.y, sig_S, pk.p.x.n, pk.p.y.n]
    args = " ".join(map(str, args))

    M0 = msg.hex()[:64]
    M1 = msg.hex()[64:]
    b0 = BitArray(int(M0, 16).to_bytes(32, "big")).bin
    b1 = BitArray(int(M1, 16).to_bytes(32, "big")).bin
    args = args + " " + " ".join(b0 + b1)

    args = args + " " + " ".join(map(str,[timestamp, eid, nonce]))



    with open(path, "w+") as file:
        for l in args:
            file.write(l)



if __name__ == "__main__":



    eid = 1

    employeeHash = generateEmployeeHash(eid)

    terminalSignature = terminalSign(employeeHash['employeeHash'])


    path = './zokrates_args'
    write_input_for_zokrates_cli(terminalSignature["pk"], terminalSignature["sig"], terminalSignature["msg"], path, terminalSignature["timestamp"], employeeHash["eid"], employeeHash["nonce"])
