// This file is LGPL3 Licensed

/**
 * @title Elliptic curve operations on twist points for alt_bn128
 * @author Mustafa Al-Bassam (mus@musalbas.com)
 */
library BN256G2 {
    uint256 internal constant FIELD_MODULUS = 0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47;
    uint256 internal constant TWISTBX = 0x2b149d40ceb8aaae81be18991be06ac3b5b4c5e559dbefa33267e6dc24a138e5;
    uint256 internal constant TWISTBY = 0x9713b03af0fed4cd2cafadeed8fdf4a74fa084e52d1852e4a2bd0685c315d2;
    uint internal constant PTXX = 0;
    uint internal constant PTXY = 1;
    uint internal constant PTYX = 2;
    uint internal constant PTYY = 3;
    uint internal constant PTZX = 4;
    uint internal constant PTZY = 5;

    /**
     * @notice Add two twist points
     * @param pt1xx Coefficient 1 of x on point 1
     * @param pt1xy Coefficient 2 of x on point 1
     * @param pt1yx Coefficient 1 of y on point 1
     * @param pt1yy Coefficient 2 of y on point 1
     * @param pt2xx Coefficient 1 of x on point 2
     * @param pt2xy Coefficient 2 of x on point 2
     * @param pt2yx Coefficient 1 of y on point 2
     * @param pt2yy Coefficient 2 of y on point 2
     * @return (pt3xx, pt3xy, pt3yx, pt3yy)
     */
    function ECTwistAdd(
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy,
        uint256 pt2xx, uint256 pt2xy,
        uint256 pt2yx, uint256 pt2yy
    ) public pure returns (
        uint256, uint256,
        uint256, uint256
    ) {
        if (
            pt1xx == 0 && pt1xy == 0 &&
            pt1yx == 0 && pt1yy == 0
        ) {
            if (!(
                pt2xx == 0 && pt2xy == 0 &&
                pt2yx == 0 && pt2yy == 0
            )) {
                assert(_isOnCurve(
                    pt2xx, pt2xy,
                    pt2yx, pt2yy
                ));
            }
            return (
                pt2xx, pt2xy,
                pt2yx, pt2yy
            );
        } else if (
            pt2xx == 0 && pt2xy == 0 &&
            pt2yx == 0 && pt2yy == 0
        ) {
            assert(_isOnCurve(
                pt1xx, pt1xy,
                pt1yx, pt1yy
            ));
            return (
                pt1xx, pt1xy,
                pt1yx, pt1yy
            );
        }

        assert(_isOnCurve(
            pt1xx, pt1xy,
            pt1yx, pt1yy
        ));
        assert(_isOnCurve(
            pt2xx, pt2xy,
            pt2yx, pt2yy
        ));

        uint256[6] memory pt3 = _ECTwistAddJacobian(
            pt1xx, pt1xy,
            pt1yx, pt1yy,
            1,     0,
            pt2xx, pt2xy,
            pt2yx, pt2yy,
            1,     0
        );

        return _fromJacobian(
            pt3[PTXX], pt3[PTXY],
            pt3[PTYX], pt3[PTYY],
            pt3[PTZX], pt3[PTZY]
        );
    }

    /**
     * @notice Multiply a twist point by a scalar
     * @param s     Scalar to multiply by
     * @param pt1xx Coefficient 1 of x
     * @param pt1xy Coefficient 2 of x
     * @param pt1yx Coefficient 1 of y
     * @param pt1yy Coefficient 2 of y
     * @return (pt2xx, pt2xy, pt2yx, pt2yy)
     */
    function ECTwistMul(
        uint256 s,
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy
    ) public pure returns (
        uint256, uint256,
        uint256, uint256
    ) {
        uint256 pt1zx = 1;
        if (
            pt1xx == 0 && pt1xy == 0 &&
            pt1yx == 0 && pt1yy == 0
        ) {
            pt1xx = 1;
            pt1yx = 1;
            pt1zx = 0;
        } else {
            assert(_isOnCurve(
                pt1xx, pt1xy,
                pt1yx, pt1yy
            ));
        }

        uint256[6] memory pt2 = _ECTwistMulJacobian(
            s,
            pt1xx, pt1xy,
            pt1yx, pt1yy,
            pt1zx, 0
        );

        return _fromJacobian(
            pt2[PTXX], pt2[PTXY],
            pt2[PTYX], pt2[PTYY],
            pt2[PTZX], pt2[PTZY]
        );
    }

    /**
     * @notice Get the field modulus
     * @return The field modulus
     */
    function GetFieldModulus() public pure returns (uint256) {
        return FIELD_MODULUS;
    }

    function submod(uint256 a, uint256 b, uint256 n) internal pure returns (uint256) {
        return addmod(a, n - b, n);
    }

    function _FQ2Mul(
        uint256 xx, uint256 xy,
        uint256 yx, uint256 yy
    ) internal pure returns(uint256, uint256) {
        return (
            submod(mulmod(xx, yx, FIELD_MODULUS), mulmod(xy, yy, FIELD_MODULUS), FIELD_MODULUS),
            addmod(mulmod(xx, yy, FIELD_MODULUS), mulmod(xy, yx, FIELD_MODULUS), FIELD_MODULUS)
        );
    }

    function _FQ2Muc(
        uint256 xx, uint256 xy,
        uint256 c
    ) internal pure returns(uint256, uint256) {
        return (
            mulmod(xx, c, FIELD_MODULUS),
            mulmod(xy, c, FIELD_MODULUS)
        );
    }

    function _FQ2Add(
        uint256 xx, uint256 xy,
        uint256 yx, uint256 yy
    ) internal pure returns(uint256, uint256) {
        return (
            addmod(xx, yx, FIELD_MODULUS),
            addmod(xy, yy, FIELD_MODULUS)
        );
    }

    function _FQ2Sub(
        uint256 xx, uint256 xy,
        uint256 yx, uint256 yy
    ) internal pure returns(uint256 rx, uint256 ry) {
        return (
            submod(xx, yx, FIELD_MODULUS),
            submod(xy, yy, FIELD_MODULUS)
        );
    }

    function _FQ2Div(
        uint256 xx, uint256 xy,
        uint256 yx, uint256 yy
    ) internal pure returns(uint256, uint256) {
        (yx, yy) = _FQ2Inv(yx, yy);
        return _FQ2Mul(xx, xy, yx, yy);
    }

    function _FQ2Inv(uint256 x, uint256 y) internal pure returns(uint256, uint256) {
        uint256 inv = _modInv(addmod(mulmod(y, y, FIELD_MODULUS), mulmod(x, x, FIELD_MODULUS), FIELD_MODULUS), FIELD_MODULUS);
        return (
            mulmod(x, inv, FIELD_MODULUS),
            FIELD_MODULUS - mulmod(y, inv, FIELD_MODULUS)
        );
    }

    function _isOnCurve(
        uint256 xx, uint256 xy,
        uint256 yx, uint256 yy
    ) internal pure returns (bool) {
        uint256 yyx;
        uint256 yyy;
        uint256 xxxx;
        uint256 xxxy;
        (yyx, yyy) = _FQ2Mul(yx, yy, yx, yy);
        (xxxx, xxxy) = _FQ2Mul(xx, xy, xx, xy);
        (xxxx, xxxy) = _FQ2Mul(xxxx, xxxy, xx, xy);
        (yyx, yyy) = _FQ2Sub(yyx, yyy, xxxx, xxxy);
        (yyx, yyy) = _FQ2Sub(yyx, yyy, TWISTBX, TWISTBY);
        return yyx == 0 && yyy == 0;
    }

    function _modInv(uint256 a, uint256 n) internal pure returns(uint256 t) {
        t = 0;
        uint256 newT = 1;
        uint256 r = n;
        uint256 newR = a;
        uint256 q;
        while (newR != 0) {
            q = r / newR;
            (t, newT) = (newT, submod(t, mulmod(q, newT, n), n));
            (r, newR) = (newR, r - q * newR);
        }
    }

    function _fromJacobian(
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy,
        uint256 pt1zx, uint256 pt1zy
    ) internal pure returns (
        uint256 pt2xx, uint256 pt2xy,
        uint256 pt2yx, uint256 pt2yy
    ) {
        uint256 invzx;
        uint256 invzy;
        (invzx, invzy) = _FQ2Inv(pt1zx, pt1zy);
        (pt2xx, pt2xy) = _FQ2Mul(pt1xx, pt1xy, invzx, invzy);
        (pt2yx, pt2yy) = _FQ2Mul(pt1yx, pt1yy, invzx, invzy);
    }

    function _ECTwistAddJacobian(
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy,
        uint256 pt1zx, uint256 pt1zy,
        uint256 pt2xx, uint256 pt2xy,
        uint256 pt2yx, uint256 pt2yy,
        uint256 pt2zx, uint256 pt2zy) internal pure returns (uint256[6] memory pt3) {
            if (pt1zx == 0 && pt1zy == 0) {
                (
                    pt3[PTXX], pt3[PTXY],
                    pt3[PTYX], pt3[PTYY],
                    pt3[PTZX], pt3[PTZY]
                ) = (
                    pt2xx, pt2xy,
                    pt2yx, pt2yy,
                    pt2zx, pt2zy
                );
                return pt3;
            } else if (pt2zx == 0 && pt2zy == 0) {
                (
                    pt3[PTXX], pt3[PTXY],
                    pt3[PTYX], pt3[PTYY],
                    pt3[PTZX], pt3[PTZY]
                ) = (
                    pt1xx, pt1xy,
                    pt1yx, pt1yy,
                    pt1zx, pt1zy
                );
                return pt3;
            }

            (pt2yx,     pt2yy)     = _FQ2Mul(pt2yx, pt2yy, pt1zx, pt1zy); // U1 = y2 * z1
            (pt3[PTYX], pt3[PTYY]) = _FQ2Mul(pt1yx, pt1yy, pt2zx, pt2zy); // U2 = y1 * z2
            (pt2xx,     pt2xy)     = _FQ2Mul(pt2xx, pt2xy, pt1zx, pt1zy); // V1 = x2 * z1
            (pt3[PTZX], pt3[PTZY]) = _FQ2Mul(pt1xx, pt1xy, pt2zx, pt2zy); // V2 = x1 * z2

            if (pt2xx == pt3[PTZX] && pt2xy == pt3[PTZY]) {
                if (pt2yx == pt3[PTYX] && pt2yy == pt3[PTYY]) {
                    (
                        pt3[PTXX], pt3[PTXY],
                        pt3[PTYX], pt3[PTYY],
                        pt3[PTZX], pt3[PTZY]
                    ) = _ECTwistDoubleJacobian(pt1xx, pt1xy, pt1yx, pt1yy, pt1zx, pt1zy);
                    return pt3;
                }
                (
                    pt3[PTXX], pt3[PTXY],
                    pt3[PTYX], pt3[PTYY],
                    pt3[PTZX], pt3[PTZY]
                ) = (
                    1, 0,
                    1, 0,
                    0, 0
                );
                return pt3;
            }

            (pt2zx,     pt2zy)     = _FQ2Mul(pt1zx, pt1zy, pt2zx,     pt2zy);     // W = z1 * z2
            (pt1xx,     pt1xy)     = _FQ2Sub(pt2yx, pt2yy, pt3[PTYX], pt3[PTYY]); // U = U1 - U2
            (pt1yx,     pt1yy)     = _FQ2Sub(pt2xx, pt2xy, pt3[PTZX], pt3[PTZY]); // V = V1 - V2
            (pt1zx,     pt1zy)     = _FQ2Mul(pt1yx, pt1yy, pt1yx,     pt1yy);     // V_squared = V * V
            (pt2yx,     pt2yy)     = _FQ2Mul(pt1zx, pt1zy, pt3[PTZX], pt3[PTZY]); // V_squared_times_V2 = V_squared * V2
            (pt1zx,     pt1zy)     = _FQ2Mul(pt1zx, pt1zy, pt1yx,     pt1yy);     // V_cubed = V * V_squared
            (pt3[PTZX], pt3[PTZY]) = _FQ2Mul(pt1zx, pt1zy, pt2zx,     pt2zy);     // newz = V_cubed * W
            (pt2xx,     pt2xy)     = _FQ2Mul(pt1xx, pt1xy, pt1xx,     pt1xy);     // U * U
            (pt2xx,     pt2xy)     = _FQ2Mul(pt2xx, pt2xy, pt2zx,     pt2zy);     // U * U * W
            (pt2xx,     pt2xy)     = _FQ2Sub(pt2xx, pt2xy, pt1zx,     pt1zy);     // U * U * W - V_cubed
            (pt2zx,     pt2zy)     = _FQ2Muc(pt2yx, pt2yy, 2);                    // 2 * V_squared_times_V2
            (pt2xx,     pt2xy)     = _FQ2Sub(pt2xx, pt2xy, pt2zx,     pt2zy);     // A = U * U * W - V_cubed - 2 * V_squared_times_V2
            (pt3[PTXX], pt3[PTXY]) = _FQ2Mul(pt1yx, pt1yy, pt2xx,     pt2xy);     // newx = V * A
            (pt1yx,     pt1yy)     = _FQ2Sub(pt2yx, pt2yy, pt2xx,     pt2xy);     // V_squared_times_V2 - A
            (pt1yx,     pt1yy)     = _FQ2Mul(pt1xx, pt1xy, pt1yx,     pt1yy);     // U * (V_squared_times_V2 - A)
            (pt1xx,     pt1xy)     = _FQ2Mul(pt1zx, pt1zy, pt3[PTYX], pt3[PTYY]); // V_cubed * U2
            (pt3[PTYX], pt3[PTYY]) = _FQ2Sub(pt1yx, pt1yy, pt1xx,     pt1xy);     // newy = U * (V_squared_times_V2 - A) - V_cubed * U2
    }

    function _ECTwistDoubleJacobian(
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy,
        uint256 pt1zx, uint256 pt1zy
    ) internal pure returns(
        uint256 pt2xx, uint256 pt2xy,
        uint256 pt2yx, uint256 pt2yy,
        uint256 pt2zx, uint256 pt2zy
    ) {
        (pt2xx, pt2xy) = _FQ2Muc(pt1xx, pt1xy, 3);            // 3 * x
        (pt2xx, pt2xy) = _FQ2Mul(pt2xx, pt2xy, pt1xx, pt1xy); // W = 3 * x * x
        (pt1zx, pt1zy) = _FQ2Mul(pt1yx, pt1yy, pt1zx, pt1zy); // S = y * z
        (pt2yx, pt2yy) = _FQ2Mul(pt1xx, pt1xy, pt1yx, pt1yy); // x * y
        (pt2yx, pt2yy) = _FQ2Mul(pt2yx, pt2yy, pt1zx, pt1zy); // B = x * y * S
        (pt1xx, pt1xy) = _FQ2Mul(pt2xx, pt2xy, pt2xx, pt2xy); // W * W
        (pt2zx, pt2zy) = _FQ2Muc(pt2yx, pt2yy, 8);            // 8 * B
        (pt1xx, pt1xy) = _FQ2Sub(pt1xx, pt1xy, pt2zx, pt2zy); // H = W * W - 8 * B
        (pt2zx, pt2zy) = _FQ2Mul(pt1zx, pt1zy, pt1zx, pt1zy); // S_squared = S * S
        (pt2yx, pt2yy) = _FQ2Muc(pt2yx, pt2yy, 4);            // 4 * B
        (pt2yx, pt2yy) = _FQ2Sub(pt2yx, pt2yy, pt1xx, pt1xy); // 4 * B - H
        (pt2yx, pt2yy) = _FQ2Mul(pt2yx, pt2yy, pt2xx, pt2xy); // W * (4 * B - H)
        (pt2xx, pt2xy) = _FQ2Muc(pt1yx, pt1yy, 8);            // 8 * y
        (pt2xx, pt2xy) = _FQ2Mul(pt2xx, pt2xy, pt1yx, pt1yy); // 8 * y * y
        (pt2xx, pt2xy) = _FQ2Mul(pt2xx, pt2xy, pt2zx, pt2zy); // 8 * y * y * S_squared
        (pt2yx, pt2yy) = _FQ2Sub(pt2yx, pt2yy, pt2xx, pt2xy); // newy = W * (4 * B - H) - 8 * y * y * S_squared
        (pt2xx, pt2xy) = _FQ2Muc(pt1xx, pt1xy, 2);            // 2 * H
        (pt2xx, pt2xy) = _FQ2Mul(pt2xx, pt2xy, pt1zx, pt1zy); // newx = 2 * H * S
        (pt2zx, pt2zy) = _FQ2Mul(pt1zx, pt1zy, pt2zx, pt2zy); // S * S_squared
        (pt2zx, pt2zy) = _FQ2Muc(pt2zx, pt2zy, 8);            // newz = 8 * S * S_squared
    }

    function _ECTwistMulJacobian(
        uint256 d,
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy,
        uint256 pt1zx, uint256 pt1zy
    ) internal pure returns(uint256[6] memory pt2) {
        while (d != 0) {
            if ((d & 1) != 0) {
                pt2 = _ECTwistAddJacobian(
                    pt2[PTXX], pt2[PTXY],
                    pt2[PTYX], pt2[PTYY],
                    pt2[PTZX], pt2[PTZY],
                    pt1xx, pt1xy,
                    pt1yx, pt1yy,
                    pt1zx, pt1zy);
            }
            (
                pt1xx, pt1xy,
                pt1yx, pt1yy,
                pt1zx, pt1zy
            ) = _ECTwistDoubleJacobian(
                pt1xx, pt1xy,
                pt1yx, pt1yy,
                pt1zx, pt1zy
            );

            d = d / 2;
        }
    }
}


// This file is MIT Licensed.
//
// Copyright 2017 Christian Reitwiessner
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

pragma solidity ^0.5.0;
library Pairing {
    struct G1Point {
        uint X;
        uint Y;
    }
    // Encoding of field elements is: X[0] * z + X[1]
    struct G2Point {
        uint[2] X;
        uint[2] Y;
    }
    /// @return the generator of G1
    function P1() pure internal returns (G1Point memory) {
        return G1Point(1, 2);
    }
    /// @return the generator of G2
    function P2() pure internal returns (G2Point memory) {
        return G2Point(
            [11559732032986387107991004021392285783925812861821192530917403151452391805634,
             10857046999023057135944570762232829481370756359578518086990519993285655852781],
            [4082367875863433681332203403145435568316851327593401208105741076214120093531,
             8495653923123431417604973247489272438418190587263600148770280649306958101930]
        );
    }
    /// @return the negation of p, i.e. p.addition(p.negate()) should be zero.
    function negate(G1Point memory p) pure internal returns (G1Point memory) {
        // The prime q in the base field F_q for G1
        uint q = 21888242871839275222246405745257275088696311157297823662689037894645226208583;
        if (p.X == 0 && p.Y == 0)
            return G1Point(0, 0);
        return G1Point(p.X, q - (p.Y % q));
    }
    /// @return the sum of two points of G1
    function addition(G1Point memory p1, G1Point memory p2) internal returns (G1Point memory r) {
        uint[4] memory input;
        input[0] = p1.X;
        input[1] = p1.Y;
        input[2] = p2.X;
        input[3] = p2.Y;
        bool success;
        assembly {
            success := call(sub(gas, 2000), 6, 0, input, 0xc0, r, 0x60)
            // Use "invalid" to make gas estimation work
            switch success case 0 { invalid() }
        }
        require(success);
    }
    /// @return the sum of two points of G2
    function addition(G2Point memory p1, G2Point memory p2) internal pure returns (G2Point memory r) {
        (r.X[1], r.X[0], r.Y[1], r.Y[0]) = BN256G2.ECTwistAdd(p1.X[1],p1.X[0],p1.Y[1],p1.Y[0],p2.X[1],p2.X[0],p2.Y[1],p2.Y[0]);
    }
    /// @return the product of a point on G1 and a scalar, i.e.
    /// p == p.scalar_mul(1) and p.addition(p) == p.scalar_mul(2) for all points p.
    function scalar_mul(G1Point memory p, uint s) internal returns (G1Point memory r) {
        uint[3] memory input;
        input[0] = p.X;
        input[1] = p.Y;
        input[2] = s;
        bool success;
        assembly {
            success := call(sub(gas, 2000), 7, 0, input, 0x80, r, 0x60)
            // Use "invalid" to make gas estimation work
            switch success case 0 { invalid() }
        }
        require (success);
    }
    /// @return the result of computing the pairing check
    /// e(p1[0], p2[0]) *  .... * e(p1[n], p2[n]) == 1
    /// For example pairing([P1(), P1().negate()], [P2(), P2()]) should
    /// return true.
    function pairing(G1Point[] memory p1, G2Point[] memory p2) internal returns (bool) {
        require(p1.length == p2.length);
        uint elements = p1.length;
        uint inputSize = elements * 6;
        uint[] memory input = new uint[](inputSize);
        for (uint i = 0; i < elements; i++)
        {
            input[i * 6 + 0] = p1[i].X;
            input[i * 6 + 1] = p1[i].Y;
            input[i * 6 + 2] = p2[i].X[0];
            input[i * 6 + 3] = p2[i].X[1];
            input[i * 6 + 4] = p2[i].Y[0];
            input[i * 6 + 5] = p2[i].Y[1];
        }
        uint[1] memory out;
        bool success;
        assembly {
            success := call(sub(gas, 2000), 8, 0, add(input, 0x20), mul(inputSize, 0x20), out, 0x20)
            // Use "invalid" to make gas estimation work
            switch success case 0 { invalid() }
        }
        require(success);
        return out[0] != 0;
    }
    /// Convenience method for a pairing check for two pairs.
    function pairingProd2(G1Point memory a1, G2Point memory a2, G1Point memory b1, G2Point memory b2) internal returns (bool) {
        G1Point[] memory p1 = new G1Point[](2);
        G2Point[] memory p2 = new G2Point[](2);
        p1[0] = a1;
        p1[1] = b1;
        p2[0] = a2;
        p2[1] = b2;
        return pairing(p1, p2);
    }
    /// Convenience method for a pairing check for three pairs.
    function pairingProd3(
            G1Point memory a1, G2Point memory a2,
            G1Point memory b1, G2Point memory b2,
            G1Point memory c1, G2Point memory c2
    ) internal returns (bool) {
        G1Point[] memory p1 = new G1Point[](3);
        G2Point[] memory p2 = new G2Point[](3);
        p1[0] = a1;
        p1[1] = b1;
        p1[2] = c1;
        p2[0] = a2;
        p2[1] = b2;
        p2[2] = c2;
        return pairing(p1, p2);
    }
    /// Convenience method for a pairing check for four pairs.
    function pairingProd4(
            G1Point memory a1, G2Point memory a2,
            G1Point memory b1, G2Point memory b2,
            G1Point memory c1, G2Point memory c2,
            G1Point memory d1, G2Point memory d2
    ) internal returns (bool) {
        G1Point[] memory p1 = new G1Point[](4);
        G2Point[] memory p2 = new G2Point[](4);
        p1[0] = a1;
        p1[1] = b1;
        p1[2] = c1;
        p1[3] = d1;
        p2[0] = a2;
        p2[1] = b2;
        p2[2] = c2;
        p2[3] = d2;
        return pairing(p1, p2);
    }
}

contract Verifier {
    using Pairing for *;
    struct VerifyingKey {
        Pairing.G1Point a;
        Pairing.G2Point b;
        Pairing.G2Point gamma;
        Pairing.G2Point delta;
        Pairing.G1Point[] gammaABC;
    }
    struct Proof {
        Pairing.G1Point A;
        Pairing.G2Point B;
        Pairing.G1Point C;
    }
    function verifyingKey() pure internal returns (VerifyingKey memory vk) {
        vk.a = Pairing.G1Point(uint256(0x034f45c3c618ecafc9dc1b76a2cb4789fb2d4e0ef4767a259d55ac6274027312), uint256(0x0d3f9fda463c29f28e281b0ecb216e3f3c747371a5b251e00665c280c3fad7b9));
        vk.b = Pairing.G2Point([uint256(0x26a88e13037d63b0d9978b23d8b32df832d82be6e74d0dea37b0d441b6ddba09), uint256(0x3016649ee46e0846b9bfadda88f9ad8331c1e68485c96f7f6d8c39dab0d57bc3)], [uint256(0x3040e29ab99254b2b622e35823484f9457a9db2206e7cfab341d51c07db504d2), uint256(0x1b798e8c54e94019efb165be5130e11df57cc696272fb59c522a25185a569f71)]);
        vk.gamma = Pairing.G2Point([uint256(0x046fb43949445e183eff0ca536888370c9f42a7e6ead52d6cb7f06603dd9eac1), uint256(0x25bd7310c8ee1f63c3e1778eb61a05973bb2c16cf2a9246b2973258e5383d8f6)], [uint256(0x06a3c1dfa04081ec6f789bdc8143da46d183c04438c3942db575aadc3b5967f2), uint256(0x04211e6c3f20a197c33778aec474e3704aa45916c7ae0de93d40960a5bd2425c)]);
        vk.delta = Pairing.G2Point([uint256(0x1b23380248699a26ff700e61179d69576947c487141dd46785620e4622dea551), uint256(0x03b806b0ac61960faf6ee00f7561c8c9e3ee18745d5873c151c523761540951b)], [uint256(0x0d56577e69493b17ccded0e63084cb5c363ebf94ad614ea2b8a43d6b59ba7fab), uint256(0x2066b35ab71573fab93b5d7b1ec6e1a85acaaa9898078af9c518d152301d51b9)]);
        vk.gammaABC = new Pairing.G1Point[](264);
        vk.gammaABC[0] = Pairing.G1Point(uint256(0x0128b4975b18307230a8402917568e94aa03b6348c4ee54079e7afef04724922), uint256(0x2a427b2be7b236d969873f13fffb0f5dfc7529566a8b9311cd4521c0db1d9448));
        vk.gammaABC[1] = Pairing.G1Point(uint256(0x2304db2fbe983cf3a527f28a4a6513e1dc03ff5a2e5d3703a7de392c97b6a23c), uint256(0x2c7f6d5690fbadf346a11332a9af8917698758333963f218e151623280df98b0));
        vk.gammaABC[2] = Pairing.G1Point(uint256(0x23eb1bd6afa6637a9b821ad9c9bd42b4e434405f7b68acf7c6ba23b3fde6f8f4), uint256(0x0645cfe9edd1eba40cbb7508e5e89c67dbe9007f2bfaa626851fac593a7b389b));
        vk.gammaABC[3] = Pairing.G1Point(uint256(0x27a551a09d27aeba5d702bd62ab1f13a0b10ea389688760c88f0d0ecd648ec91), uint256(0x1abd038f5c3240eefcd80ea100f8f2cecd055f91bce2d1a08511cfd0b8af585f));
        vk.gammaABC[4] = Pairing.G1Point(uint256(0x078b9f50ecad3b8e6eb42a03d8013d977b709fd41af12277f6969f6793836fff), uint256(0x0ab425bf3b56e79ae270ec773a53d0f973eb15f4400a64b3536910118d6741d8));
        vk.gammaABC[5] = Pairing.G1Point(uint256(0x1a14ebaee840157aeeaa9091923ab3da483737673de6b3c3e93c5ebd287a4a95), uint256(0x0d0e682a5da92a688d89cf3d3fb04f0b060cc56d00a5c8cdeb1261159af487ad));
        vk.gammaABC[6] = Pairing.G1Point(uint256(0x145cac89335ed1cfcb36c6408aa900317344ca8046f5e5688aa42e7905e452e8), uint256(0x09abbb17614af95a7ec3373be71af0ca608088167a58201fbd5e72bdd7d2ede2));
        vk.gammaABC[7] = Pairing.G1Point(uint256(0x0313dea1393f121df3a4d239fd35cebb8f345a687e94f361f4c48d8900d49f7d), uint256(0x11869bb636f3ef8ab3ea044849622c8341600ffdc045c9f29420b3488895c0d5));
        vk.gammaABC[8] = Pairing.G1Point(uint256(0x0e2964ad910d9ef37461413743bad131594e961f81a6ea0e8693125eb5d6f18e), uint256(0x2a538dd13a0d0dbbf4a244e71c96e655764eb2818cc5124bd4e3a5a44588c574));
        vk.gammaABC[9] = Pairing.G1Point(uint256(0x2bd5ade608fe7566b1047d4b8503481a589da8bb771fe81dea172cbfe25917c9), uint256(0x2e8191b5bf5b8256339f3c8d9f570d2b56c13fd8b5e8affee962d399ca40d565));
        vk.gammaABC[10] = Pairing.G1Point(uint256(0x270d533419b334f7c9f99d9374759faf60878215c4c43587edf51fedcf158947), uint256(0x0faf17fb75f68a80c31dec78fb2625526b755ab24eed69e680f819a5ab852ad1));
        vk.gammaABC[11] = Pairing.G1Point(uint256(0x059fec0d6d3c95bb445aadb24a88fdeb74f8fdee3b4a0c5a002ec3dd9947fa55), uint256(0x243f08236e7f72fdf90ba1adf874315286fd7d8ff3676210f17bda944583035f));
        vk.gammaABC[12] = Pairing.G1Point(uint256(0x0daaa956aa75d69a33f9df05ff4772acf6d17c4bb594fd669fd4a702bc1d016f), uint256(0x159270f4c5b567a57cd96701716e97f5193c0ff514c3ca800120f034125c98c0));
        vk.gammaABC[13] = Pairing.G1Point(uint256(0x057d6cf9e37a0457ebdebc1c515b3b458229e25b90fe3b72b9f59f7b34b06ff1), uint256(0x10a1059d8effef00b6547f7b4943ec382b721d712d5efc7d44702e100ffca7ac));
        vk.gammaABC[14] = Pairing.G1Point(uint256(0x1a902b7aaf1b68f4b6666d5f4a0ff165315a1304411f4914d78ee172a0e8c699), uint256(0x0a1512903ac63558e566c69126077529248656c03eb060bfceca07e7f67f3f23));
        vk.gammaABC[15] = Pairing.G1Point(uint256(0x1c620508e755fee5f78dfa4a7c28627ee9884cdec712c93bf2be4aba486b908e), uint256(0x2ca3b8e48a6dfb5cc11d90995ecb4114c67f62fd42774d67001deab3f1222f27));
        vk.gammaABC[16] = Pairing.G1Point(uint256(0x03250ec009e0e6907371e4997ea8b83a3ebc68a4be3ea4a62920386d71bc265c), uint256(0x155999d452760be2d0230fc309dde51d98abac8a0f4b86178fe39401a4d6f651));
        vk.gammaABC[17] = Pairing.G1Point(uint256(0x11225711f781da4093dcbd463cdb9bf471a2cba5994eb8edab772e881ea8826b), uint256(0x2fd2562fed4bfb1136b72a8842e49e976ef94c35a7a3bf09f66e6339a7e2321d));
        vk.gammaABC[18] = Pairing.G1Point(uint256(0x17cd123b278b679ac768f314de3ca4007d0d142a0f6317b773f3a952ab77be07), uint256(0x0a5b3bf4f261192bed4917c4a1a125a448d8565706717a69511cb50e89af19b3));
        vk.gammaABC[19] = Pairing.G1Point(uint256(0x128cc8d655ef39bcb24f0e4e4cdeb38538c06c273c806d2f4f81cc3a2ca10129), uint256(0x25bb072d7c281391a210c7f3c65244d9ff5c748e7f54b807c4cd20fb579be8bf));
        vk.gammaABC[20] = Pairing.G1Point(uint256(0x217d49a389b0e44a223237ee363cb79d8fffedf533fd0dc0ae5116a566c4d2e2), uint256(0x1e2c71ad3735d388e26a2e91bef08c52e3a42fe36a6ccf00dcd97e9607909f64));
        vk.gammaABC[21] = Pairing.G1Point(uint256(0x2c63922ea932f3ee9b69f0479609cf76e41b157f966d34ece6fab414ffc1d832), uint256(0x16a88a3f7690f853c406e89ba8b42cc1ab3cc9ed36a1135503cd915bd7494346));
        vk.gammaABC[22] = Pairing.G1Point(uint256(0x1b749547aca6f047b555e2f06a7e0cd21adc37b46b4fe863e30bba4c8a631ea1), uint256(0x2eab3f533739fd711bc9348a89e0b52a27f2fca6efa8434ef381d18a17470ad8));
        vk.gammaABC[23] = Pairing.G1Point(uint256(0x08107c0b37a3e2f28dfb5c2db1063a06fc4ddeb0e95a287d37ba83468821307d), uint256(0x26da7289c611531b40a7b278abf286d66088df3e3f7b4b4980736c44a42eb818));
        vk.gammaABC[24] = Pairing.G1Point(uint256(0x00827ba7be7679db69ff1bd364d1c59801f871fdcf88549150b537a2e9e0375c), uint256(0x1973468f1233c605c80b00358168aac8950ec571b99966bf161485267f567cd8));
        vk.gammaABC[25] = Pairing.G1Point(uint256(0x0560805bd258a23c4e1fc79ec3faf7a25d43b8d71bbd928c017f69f6cc751a6d), uint256(0x2bce742eb8229d7416c40a825ff0f32a55f2658c66ea4ea5a1c667c2dfca02f8));
        vk.gammaABC[26] = Pairing.G1Point(uint256(0x2d7e4f5d8586778beddc1aa6a77ae05c4c7b2b6b29a52350e925232fcd423e80), uint256(0x0503dc7bae518be8668a81084c317df96ddbc9b026a34ad1e6af0de76bf8ce17));
        vk.gammaABC[27] = Pairing.G1Point(uint256(0x145dac9849ecb59ec0b0201d0ab12cdecee70e311233e13b130603b0f0368b02), uint256(0x0cf63d824998f6c4faf16f7ffdf079f2c4f16ac01c769f2ced9340ca36c150d5));
        vk.gammaABC[28] = Pairing.G1Point(uint256(0x1142c40609e230d1457933628300e79f34ff845bb360ec791c99d3c6b8a94ad7), uint256(0x1bb36dc2c6620bcae7797355788571fbcfae3889738599f34ca114fd60f0987c));
        vk.gammaABC[29] = Pairing.G1Point(uint256(0x1e5141d47ebdc14a37ac316910fbeb4e02583f7270c49328f651f051b3f226c3), uint256(0x1e0949a5de21bdc7a2f8ac7de5f286f6f4961978e1937f5afcc1c1182d1d2d0a));
        vk.gammaABC[30] = Pairing.G1Point(uint256(0x28aeb2c21a6ff9cb2b5761641751c7b2c9f9b5018f2c518dcdf46e403081b4e4), uint256(0x0d4f37430de7e447c5a030b49b53615a5881d0c89d17f2c6476b0d8846afd255));
        vk.gammaABC[31] = Pairing.G1Point(uint256(0x0f9863a6a3be705d38dfd4a68bb50a94b9b218943f7616a9b513a543a204a4f8), uint256(0x0da2c1a9c55089a10f5491c2b2fc176d4b8d27e3cd3a774f42a34a32cbda0c83));
        vk.gammaABC[32] = Pairing.G1Point(uint256(0x188e7b815105f6d3d3a04e76b6ccedeef89f013fcf4a67abdbcbf71d7f9d528f), uint256(0x0ef16a3379323a4ddcc78ba39af7e1ac8eccefda56386db4b2bb056a9f2663e7));
        vk.gammaABC[33] = Pairing.G1Point(uint256(0x1c88ba243fcb9600afcbd12534afa26e35a60ffa58d4cafc2b9494efd9a65a5b), uint256(0x2bed0d3797b89d1778c2e72b48f26c04465bc1bcaf7d6668206b1f149991844d));
        vk.gammaABC[34] = Pairing.G1Point(uint256(0x219aa9b3fd342349923478bf3bdccc4ec2b44115f64b948b3c4424abc63592ac), uint256(0x0c675d4fc73f987ebb815701be55eb47097b0f0896745a0380ebc4304f870176));
        vk.gammaABC[35] = Pairing.G1Point(uint256(0x228db9b31713e5959d6c8e8efefe044afa35f7bcb5267eefcf780bdb73e4d5f4), uint256(0x235007afd468024627b386d82883b12f2293770ffe2d2cea9bfb5bdd90242c05));
        vk.gammaABC[36] = Pairing.G1Point(uint256(0x034afcc5ea8b8278d7c4edc768138c197df971eb526328dba493d679b4d7a6b8), uint256(0x16b7bbc0d324649911f49edd519d918502c6c40bd97a84e1f022fff866efc699));
        vk.gammaABC[37] = Pairing.G1Point(uint256(0x1ea7e7da0c804cb415d0cbd86963568569f98f5a5201f0994e4baf0b7cd544af), uint256(0x16d1dec005322ccf269d2198840b56a8d503d415fc4eacd6ba9141a62415b4fc));
        vk.gammaABC[38] = Pairing.G1Point(uint256(0x0ca84e004b73a94849ea2c65a9f86e60032a600b7ecf31893ee0ff5a30b834a2), uint256(0x3035943be6d371aedc62054ae33b5d1de24e6dea5b7f1651dd20c70017dd8c23));
        vk.gammaABC[39] = Pairing.G1Point(uint256(0x01546fdd5ce45424841909278fc5a71942dd3bfb46ae44c1c718eb15ebc413f5), uint256(0x239115d132cd26407d58d169e6f9ff2d3b7d04623a59a31fc81f6df91e16e27f));
        vk.gammaABC[40] = Pairing.G1Point(uint256(0x275af2050160f09d7684a218a7884f0bb00ef21fa16be8063e23144fd73efcde), uint256(0x0d1fc67fefdcfb6fb0b5aaf521281043b9f89d88c0de3463ee1f68e6717f85d5));
        vk.gammaABC[41] = Pairing.G1Point(uint256(0x29c9be08368b31e0131f65dce3855951f0b65541d8d7dd62aeafd99de9efbcfb), uint256(0x2788f25d74971d6f4b689cba3aa43c4f79e5da2871ff23f4ccd371595630df53));
        vk.gammaABC[42] = Pairing.G1Point(uint256(0x2c8dd35600d33c65f1f9eb9ad93ab5d7c23f8941638db3fd83e00ae9581e9fa1), uint256(0x0df7768e0be07e832487f5f8ed037e1296cf56cfcb1ec0395c247a89d5a52cea));
        vk.gammaABC[43] = Pairing.G1Point(uint256(0x15b141898ddd90cba223c2b78ed8c7f8405bd519f1f632faf5a0eab39de95728), uint256(0x081d84289da30376a817ad1379a1ece566753cbf9fc3f7ffbadfd0e4e8dcd663));
        vk.gammaABC[44] = Pairing.G1Point(uint256(0x17b36ad06b82ddadfd40d7865d7445ca1d3ad947023f135be315c050e2367c2d), uint256(0x2cb3709850a173ee28c4ae74da042f73f281d4c5387b7baa05f062c3c299b67e));
        vk.gammaABC[45] = Pairing.G1Point(uint256(0x09120337d878951e538a05a70222b602cf8b4a484fa2776d3ff2ff7fdbbde64f), uint256(0x0d7c31894966eb21528f4cc76774fc3557896deb0e4352cd20bab3dc69f3d5a3));
        vk.gammaABC[46] = Pairing.G1Point(uint256(0x24f701c8dae4579abffad681c1ff7cf35ea44c6d0d58c540e9b01806ac9c4881), uint256(0x135233abb306f407b9c0e0597d354b1821f8facfdf8e33aff16d7056c9bb82e9));
        vk.gammaABC[47] = Pairing.G1Point(uint256(0x09051f5f01b0620d1e0b50dea68b10b7f74a88ec116bd667052506b4c7425a2f), uint256(0x15aa8e80dc18c72915ce1c88029aeeb91d61a87e87c09870f882fa19930e5b90));
        vk.gammaABC[48] = Pairing.G1Point(uint256(0x1082d068ce0e28f4aba53851fc2a9c9175b0bcbc83375f49c964d37d5eb33f14), uint256(0x22b721fff710363b89c687d4656a42d1d18a906fce85d7d85a81a2969c1ba075));
        vk.gammaABC[49] = Pairing.G1Point(uint256(0x0c77c4320720bb5c3235e52cd1885dc3baae7fceea5d44afe43c4528e0596b55), uint256(0x280fe5ce21a9d4372c4d8a8f8398cd2a2183174f6540b6f1a73ab20a39b0ac99));
        vk.gammaABC[50] = Pairing.G1Point(uint256(0x19e18aa8278d7729c205bec74163400f8dedd25fba10c599ca9ca2a8b0f12794), uint256(0x2fc01de5d829ac8cd3b8c360dcac83b0b6890c86297fc09ae9aff617eb42c605));
        vk.gammaABC[51] = Pairing.G1Point(uint256(0x14a431f20c816bd21c2f6c860c9be9d180375a6bc3d82213514d7a7313d3deb9), uint256(0x2f981e1d2a6aff6f39120de88052b507d5bbf67a9ef0f1bee669162c527ab27e));
        vk.gammaABC[52] = Pairing.G1Point(uint256(0x2600de6cc0c577ba0e0ea8ad0d01e386e0ca5214035bb67cb171fc6f533cf280), uint256(0x220b50d509f3cc27d18f751504b136c49a081b898c99f5b913a07da741b20e78));
        vk.gammaABC[53] = Pairing.G1Point(uint256(0x283e55a58833c924a79df0598d035cda73b2993c9f32ae8a58965562f64cdde3), uint256(0x1d4893a34b8930f5584d4a22b93a765994f7d940f9fd0b9e140bc96fe0857531));
        vk.gammaABC[54] = Pairing.G1Point(uint256(0x24d9c9e99e52aeb297afe8c72ef931356b047c867b89fa55a66b91f1cd69e79f), uint256(0x04a1936c5002716acff6f99be44f2df0514c5bbff669157dd40ac736b9254907));
        vk.gammaABC[55] = Pairing.G1Point(uint256(0x10ddd740d41363a0437de35b9c6ccba72cd4216f3181f5f525257b27adc3c9d4), uint256(0x1448b5e00475d23d5a75f3d07caf101e3d1519fc7943f969d06b2428932bdb63));
        vk.gammaABC[56] = Pairing.G1Point(uint256(0x0c9244a9551831ea60c77b5718c2ff8ce2d6ce71432c0db23692e9236ed813cf), uint256(0x1f45dacc6de0664d9bf774582cee46be12b4bd64933f1bf4cae30ae811edeabe));
        vk.gammaABC[57] = Pairing.G1Point(uint256(0x0abf453835d25c76e06788f8d418eb37081c6cfe00b3305de86a4b6ee20f6e87), uint256(0x23e78bdf5f07a1f677f17185e271039dd753a865c5b7b4bb68e0ae30682cce87));
        vk.gammaABC[58] = Pairing.G1Point(uint256(0x0f83261890a0795fe3057e6b221b5f486b0f1d949e7e544fb79c47a2aa89b671), uint256(0x1bc67ac111790609715d98b6fce59db0ec9c8a71c09aee84f31572102a48d218));
        vk.gammaABC[59] = Pairing.G1Point(uint256(0x24181b89292bd25c411ac7431da03f55e0c15a140d60665ad1655bfd06989da2), uint256(0x0077e94531b1b4c1f3d07af798d000ff65ba85db0aa248e05b2792cd1cde5c68));
        vk.gammaABC[60] = Pairing.G1Point(uint256(0x1ae6d940e4ff817037b28c6fa44c89fb6157979a5e3826260f8f3f3c8c5ded65), uint256(0x0ad031cc4cb1535014490a4a64c93dd23e52b9bd79cbf12480480636d0d7c62a));
        vk.gammaABC[61] = Pairing.G1Point(uint256(0x122565440060c3829fcb2a71dfb5bf8e16eaa0bf676a371a88e7817bb2577a33), uint256(0x1a9c1919ca12b257a43a5c59e2b3d962cd30f9595753036205ec18369a499c37));
        vk.gammaABC[62] = Pairing.G1Point(uint256(0x12dc8c9f2c77e7865f27f4d9b2660bc14467a848d0b5c1374235f0f1629aa48f), uint256(0x06c84e337c115ab84680c571e0fce716c1ccb24190cc73b1f79f14d171e86df3));
        vk.gammaABC[63] = Pairing.G1Point(uint256(0x197714b1567104abbea2bc21cb38d1d6f5016e972589e51558798895a9d3087c), uint256(0x08fb3d36665b3ddae2126ea03d428feaa60f6602b104d1bbefc9d4e7668bac1e));
        vk.gammaABC[64] = Pairing.G1Point(uint256(0x27691d435e400923593ccaf35405c66a7d5fd19120ab85f704eb7ba347392e52), uint256(0x29a9083bd755356e4d1cbd1851817f70fd450e6e72ec7290efe02206cccd32fb));
        vk.gammaABC[65] = Pairing.G1Point(uint256(0x2a99bf4b9cf80f95e078d8ea35dcbe181a339967b5665747faa0eb6f73e7e490), uint256(0x1b46516124f00a9add84de60a0e9bab5f76954ae2559be4b2efb9f3d0de9c2d1));
        vk.gammaABC[66] = Pairing.G1Point(uint256(0x1e97bac51f62445447b0d31aa6e541c3db3e717a13f575f664be051f5247d83c), uint256(0x2192375fdf6b4b671f44e86f9b0f557836bcaafbe3e3c601b60fe8a824671d03));
        vk.gammaABC[67] = Pairing.G1Point(uint256(0x19102cd7dcb7f3f0df6d501735fe33d23666b85b40bd08fd4f0b6c0fc6a03bb4), uint256(0x1e5939cc8885cdb0878a110f1d80e2fe9d84038a2a20400e08817fdf60abca81));
        vk.gammaABC[68] = Pairing.G1Point(uint256(0x1b1dc2c8ad76eb98555aaba1a6662e80c03175e9a933680ac932f42d38c09c58), uint256(0x1a2d039cdbd3bc56cd54e21777be6bed5233e2a36143f7ddf911c93cefc66d43));
        vk.gammaABC[69] = Pairing.G1Point(uint256(0x2de07d8f0948949c9fc860428a34b1bf3af36b127fa73dd288db7ca27dee6f4c), uint256(0x1d30d7268c6db5439370f0d41bc1767fed2e9ebd83408e4998e1fe2a358be0b7));
        vk.gammaABC[70] = Pairing.G1Point(uint256(0x01b4c77535bff6c6bed33aedb3cc5e0b9df22eee7a9d8fd81160be20e084d4a9), uint256(0x23584cde1830e0e1684e242e45971cd4a21e42046914ca9085f0bd2f98ced759));
        vk.gammaABC[71] = Pairing.G1Point(uint256(0x0788c77a6f76b786719fca86f2f562f56fc1ae3d8b3a508ff14a4ea651e8b8f0), uint256(0x03e1936487b8991cbecbcf072bbc3eaed07ccd18abef0611fbd8027d2a7c7d03));
        vk.gammaABC[72] = Pairing.G1Point(uint256(0x131c410f337b002774d85afa568d02cc1b40178375323cb8449cf8fbe21311d3), uint256(0x141d6bf07edf5b0a535e8bb2ed7f083ea573330c3110e05b86765ec38cc6cebb));
        vk.gammaABC[73] = Pairing.G1Point(uint256(0x0e637276ecf8b4a47899bd1ee5f9fb9da4f4957c08267693bc8305d8849ca871), uint256(0x12f8be13a7e60183780e171f2c9d2e108b43145f4d7d831b5ba2456d5c23f07c));
        vk.gammaABC[74] = Pairing.G1Point(uint256(0x1a79450c9c2620e72171052afe3f82ac569095a1b50b6f6c0ac6fa07388494d1), uint256(0x252a67f549ae15dc1e6ada12cc6cf3c068174bd7fecaa899cd5b781d1019f088));
        vk.gammaABC[75] = Pairing.G1Point(uint256(0x032b655e4c18f12a1e20060bb9c5c56bde8afe444ab424f93688841fd1d0d9df), uint256(0x1fe0f5d9d845ba75a8b85f350e7c97aa8604c431d783cff4a7c34d59c4d3f0ea));
        vk.gammaABC[76] = Pairing.G1Point(uint256(0x156a213efa501d3d7dc32a5701d508c16adc5a508d9db11d5a32dffd5480e9de), uint256(0x29e68e7f44bc00ea4d2b8dd1d6aae21a23f09f7a0a42426bd2772effc9306978));
        vk.gammaABC[77] = Pairing.G1Point(uint256(0x12a8aeae0aa44a578bd125cd34eae91a97f5bd338b6bc479f8d304b1c309edb8), uint256(0x1804738ea81d7e03c8e2d43ff9ad7c915dd1031abdbd55f52d2277a5321d3cc4));
        vk.gammaABC[78] = Pairing.G1Point(uint256(0x0cf892fb3ca27933d0859a45986c25bf74e01048f399fa48d4681165a6209d30), uint256(0x19cbed43c2cb9605b88bbb6d6033c5d60dfc81de448ad3e68abffab667f16ed8));
        vk.gammaABC[79] = Pairing.G1Point(uint256(0x0cce39c2cd747c3b6ed424118fa9b9ee208c2fadba56159e1d737aabb4fefa74), uint256(0x0cd76018537259d53c61f307fdd957dbbf11b1ec7757ba86a8e1529e20c18d9b));
        vk.gammaABC[80] = Pairing.G1Point(uint256(0x2645131460054b5bc33cca4cafc4fe8a7481170a08fd2285b94e3bb57db21e02), uint256(0x27f13d33caa4accd5a7961e04079b8442a0557266b6fba268174b3592685c886));
        vk.gammaABC[81] = Pairing.G1Point(uint256(0x29e8247000b4be0c20966bb25c6c4ed9648a3213bd9b622506a6fdd4fd4bdb88), uint256(0x023331e01effe85d13f14c8205178d067985bacb7ee1be5c7e1e559955877828));
        vk.gammaABC[82] = Pairing.G1Point(uint256(0x1472a605d39d52b0b21f77e6d54ccd0265bc85717e9aea494faef893c04f0ef8), uint256(0x284586cc38680a53474194d43e8df95a81aef2559550cd45803c3e9bc8ed1fda));
        vk.gammaABC[83] = Pairing.G1Point(uint256(0x04a75aac6883057adf34ca3c8c11997a9b1376692e8a2904fe185d8d974d501b), uint256(0x27068061152c8a537a527ca4aea183efd660430710534ae087abd67fc60f5079));
        vk.gammaABC[84] = Pairing.G1Point(uint256(0x12994513b69cc7dff8e453e41126aea12131dd2335b1c5f1d65a9cf63b83970b), uint256(0x302275b3366e8bf4d13e9ef609f6d3dd57eedac6f2b355e02ce5b112d2e037eb));
        vk.gammaABC[85] = Pairing.G1Point(uint256(0x0939e4222dcc95e4d031dcb956d618116b8eb220b91cc060a15501a3ab25e916), uint256(0x01a2d08faf02b9d19f5cb2d39579e41947b93d3c69c902412139d901885f1f93));
        vk.gammaABC[86] = Pairing.G1Point(uint256(0x15ee7901185c5322063f6ea23cd5ccfe3e329f12bd9403185ae1e70799875268), uint256(0x1592d9303e7ef96f14c42917cd5555d5cfff1a01e2b83f1bc3df119f5d5d7ce5));
        vk.gammaABC[87] = Pairing.G1Point(uint256(0x15c94dcef9585cbdddfbe578c49895c5efb20fbf6d05293968357151c360b32b), uint256(0x11fa99cd039029c0d99a44c62403c8446a4c5a8ccb5280e11031087137b9cfdb));
        vk.gammaABC[88] = Pairing.G1Point(uint256(0x23bb8597894f2b4624cc7e69656ad8e5b65d63c68d2a7466c007fba5724238e4), uint256(0x123571ba523992de9468de6b20598fd3d1a4ebdd52057a7eeb9ba71eba055032));
        vk.gammaABC[89] = Pairing.G1Point(uint256(0x247d6c7a048b6800c925ae751c88bef277b2184c341298ed38b1ed74be69586f), uint256(0x171ad5f7dbec6894656368d3f274a1eb4c060f2ecd7e45b5a536b87434bb63c7));
        vk.gammaABC[90] = Pairing.G1Point(uint256(0x0e7fd48d778c2e18f20e1dff32d08cc9fe76f50d661f03c18d618d83c086d736), uint256(0x1025916d80d81a119b40162c61562e9c28e0b5c3b5f7237450b0a319d2e991cc));
        vk.gammaABC[91] = Pairing.G1Point(uint256(0x116ce3f62982fb291f3784a213dd9e8e8fc25500d6acb274a775bd62dd3eb6b4), uint256(0x272d3a00d64e8e0469c5bef0e553f11d43d8edcd7d63b12de59969b4eaf584d9));
        vk.gammaABC[92] = Pairing.G1Point(uint256(0x23d4c3ab11407adf5251fe8a72ad4ce2077c0440db9ca0d3af91618156be83d8), uint256(0x012f64c10a833f363e1d6dd82593b957f0e1c9abd6aa67265527ab69c9545f9e));
        vk.gammaABC[93] = Pairing.G1Point(uint256(0x0928ea50d2292d887ab43335ae1e5c6820b0661e8c8ac45ee8d52026a8d7bcb4), uint256(0x1897e96ac0bd813391eef5eaea5bc40be886a12038154af597b67b4686c199c3));
        vk.gammaABC[94] = Pairing.G1Point(uint256(0x20daf1dac0a62cb062a5c01dd88e5abc66d8efcd28a7dc3dc2805e9ccf4ff2ca), uint256(0x090a69baf63efb54bc9ee78577f72d220d4907b324acebbaacecaf971ea8c9a1));
        vk.gammaABC[95] = Pairing.G1Point(uint256(0x29cbc26bb03dd6e520eeff89b7e9e8fd73dfdbd6bf59995eca2fad35f12f244d), uint256(0x15a451a91b90a8fc8137933262a72edec57721e04390bb57583d0fdd84050b46));
        vk.gammaABC[96] = Pairing.G1Point(uint256(0x14dd85a19487d07e204e8b50f46d9cb507cd8cd9acbe07a84ee2af8400e4a39a), uint256(0x1df5cfcd6f88fcce46c3800db9407515fae702dc4635930dcd6ee4a44b9b8ea5));
        vk.gammaABC[97] = Pairing.G1Point(uint256(0x0a8b665fe9aaa01179397b60c4d4ec96d7832be5214290f2e6ca6347871cac60), uint256(0x2c6cfaa7b809fbe210f7a38e5a026b4289827afa36002adbfe733464e7de2fdb));
        vk.gammaABC[98] = Pairing.G1Point(uint256(0x097d8864ef1c4ece175e0c518dc64d6a2c1c848b9c28d6ce95023d71d98c6b20), uint256(0x0c33c6803b002cf35afc511753f59e2a5e31f89c971301c368ce2a18c3567b48));
        vk.gammaABC[99] = Pairing.G1Point(uint256(0x2ebf670dd278bc778da507137d6452d6853a66c607238d3295cb8369b5a33390), uint256(0x0b37dbddeb5470627ee50ce3239cf1a3d3a91932ef7981ecd8202daa3db03df9));
        vk.gammaABC[100] = Pairing.G1Point(uint256(0x240ce8753a4330ba410e985d9d793a0c49cb14e8e34b519e81a9a26196fc22b3), uint256(0x1d52e773e3c3dc2977324c9ec4d9b8956b7d6a6f030698aa6e894a97acb403dc));
        vk.gammaABC[101] = Pairing.G1Point(uint256(0x303222cdb5eda9d1526a35f8c5d108f271729573084bbee241f90a5283696f03), uint256(0x06c1b11f248abb09a93673f5347ea27f10cd006c7457cd74e7072290fa6de890));
        vk.gammaABC[102] = Pairing.G1Point(uint256(0x1eb9c37c8e268ecebcca3e3119820ec566180f708a03e6dd856449e274cfeff4), uint256(0x282989ffac5baa61fbe51095c77d97c7affc1bb2b3461308bcd8f47ffd1757a2));
        vk.gammaABC[103] = Pairing.G1Point(uint256(0x305c0ea4a139b666a07dfd5bbf7c68981f399f62e720d1663af0b5fac604edbf), uint256(0x27403b38e5ec90d285a3c2764088e44d37625d28b9d2f7fd4ecd8c58ea7d68c8));
        vk.gammaABC[104] = Pairing.G1Point(uint256(0x116eb1ba95752791c2952ea4369929e737501b15cd352361918c1efb550c4f38), uint256(0x086aeb18c093218426baf600b5b0b13398a808fc1c8806dda814e4106a26f9eb));
        vk.gammaABC[105] = Pairing.G1Point(uint256(0x1c449c2010f6026bf1c83a4ffcaf2a856a62fb15baf80cd3881f0b9f7b9f6968), uint256(0x0bf21431b133536279f48c5e5306c67267d148aa98bb8f30d352c20d83d8e368));
        vk.gammaABC[106] = Pairing.G1Point(uint256(0x14e85e07f849f1c2d4cfcd40abad6553f81569eb58b9eb1e4ef308797c800005), uint256(0x0ee80b88d46d5069960815e3b2df0d7e7fbd894350e463e70aaa4752f91c2a8f));
        vk.gammaABC[107] = Pairing.G1Point(uint256(0x11bc3a25715fe67f34e42e2188f9fae285e1b51f5b53e56bbc28739acac6903a), uint256(0x1bcb1add41ac995151f460c853e14ef190242f1b033d650c7118c9bb2f271190));
        vk.gammaABC[108] = Pairing.G1Point(uint256(0x00485aea5c6daadba063db577dde20b0dbfc030f7af33953ea91e608c3615729), uint256(0x0df5431478a91f777c5d21d12136e7d25350485df8abc96ffdfa1f1b5cabf6f7));
        vk.gammaABC[109] = Pairing.G1Point(uint256(0x2056407611f1fb5ed635e43ff50e5f047fbc8c3f2780307eb7650b8ea1a558e5), uint256(0x26babe7526661e679abc4d3f0181f4cf8a894557af0e85efad860fdbc2582d95));
        vk.gammaABC[110] = Pairing.G1Point(uint256(0x282de112653e91079f0470742fcba4137042ac8ec83f27e854e399439fee5ee1), uint256(0x1a6f364dddf6eee0142da334c2e515cd5f3915c05566a565750c60d72b656137));
        vk.gammaABC[111] = Pairing.G1Point(uint256(0x2dca340c4253c6e17e506d8ddbb6b0a6d093dcf8be1d7d6653fcc5f63c711222), uint256(0x136aaa5d59523ca3f8d35dbdec36356b90636af72e516ff26db9e0363764e4e6));
        vk.gammaABC[112] = Pairing.G1Point(uint256(0x0605dfcd071534bed94a5e4a01476e384b04bcb2cd1a62b6708f1fc04b8481b0), uint256(0x1a50c09c8e1459316bb49f3e8c488a875b3bd408c20379c30b82d6d922658ee1));
        vk.gammaABC[113] = Pairing.G1Point(uint256(0x2a9b1f7ff8eec6627c6e64265bd78078bdf57a6063292e03fc46c3c6c03d2be6), uint256(0x229cd4c1bf30d54ac7dcccb03336bd4475e5869f12e412b415175b22b1549ffe));
        vk.gammaABC[114] = Pairing.G1Point(uint256(0x00e73ae43e0910d34937e06b1b730f75cd588fbab5fb281f46b22b3b66c35afa), uint256(0x25c7af796934463b95622dd5d8d8aa6469b07dde4a450f54fda53fa6d0ecc69c));
        vk.gammaABC[115] = Pairing.G1Point(uint256(0x23f03f02b12ff1b9f0d7584a2cd95a7a8e58b1f09053dc8c5ef1117cb381ac06), uint256(0x2b55e3b96fc04e3b014649d67b69cc1458523e69020b181dc4c7e383c618a7bd));
        vk.gammaABC[116] = Pairing.G1Point(uint256(0x22838560e648968fc2fd167f79c508957deecb356d5797666e403324d7ed2d2d), uint256(0x1e4e6956dae72dd1e02b89a7485f21e9778b291b5d31cc5be7d3d34abb1c3b6e));
        vk.gammaABC[117] = Pairing.G1Point(uint256(0x2a233703ff3b0b617743eb0d69af6cf47654aff6cef403d49f086fd1ef6621c8), uint256(0x04e62889b89de667854fae48803c2ac0b62f16e2bf5f2d56f528019f1a514d1b));
        vk.gammaABC[118] = Pairing.G1Point(uint256(0x18470d9a3ffe63379239ca486d528c070ee64eae5359ef63c961cc9ef60f9239), uint256(0x20b94cbbcc6d3eaddc1603dc4df6d60f2b57b434b20ace0e7e6182b9dd97859c));
        vk.gammaABC[119] = Pairing.G1Point(uint256(0x08c1bfd2fb670b491b4e7c2cb4b3c1e70176e097a0f03ea2fbf4c5ee1cddcf99), uint256(0x10a05c1ea5e1ddab2b9671a88974feebaa143f2891e382d8d4786edaa3c3492c));
        vk.gammaABC[120] = Pairing.G1Point(uint256(0x2f87c5d4ac869d3d3f8ab0ec2c2df5702da0f3c18225a5cc715d982d56adab6f), uint256(0x0b23624b18be5ca2f3b8f10d75ee367303107c58c39b3d9d89ebcd3d924c1668));
        vk.gammaABC[121] = Pairing.G1Point(uint256(0x0f3279c50d78d86fcf3cd0c205e99126e00a4416dac8f6c94623e61b5916365e), uint256(0x2ad34756fd823f96fa85ffb3fc45c42b2fd8414f9ba00ef495e668ee293e5b7a));
        vk.gammaABC[122] = Pairing.G1Point(uint256(0x1b31a78cea70ffc90f64d610de2e7e93fa247d8d59bd3eab2599a83aeb46cfe9), uint256(0x1ef05e8b7ce8e03166af90117843a3c831437e90fef8d070a97812cfab4f2569));
        vk.gammaABC[123] = Pairing.G1Point(uint256(0x1b1ca43297e59067dcf2c20bde851c5993755215d6b1cb192c92d0c59d105a3d), uint256(0x1563831c217aae45b3ad9675901dceeefed96bb1997c9be7b00fc94794060390));
        vk.gammaABC[124] = Pairing.G1Point(uint256(0x28333e1c0d0ea94848d3ff3d51d59e44cadd6bfb2aa5aea675e85ce3181fc311), uint256(0x1f12bbd8179d1a80b9c6564009284df99b93857bba063a77a73f16ba08433566));
        vk.gammaABC[125] = Pairing.G1Point(uint256(0x02b59f2128b1cad69b3da77ff30ed211de21963dd1136a404baf22ab0a1ffc24), uint256(0x2c531cad5b4cfc91425f7d3d68aa8a1667381aa3d3c440693de7e9d6716e8667));
        vk.gammaABC[126] = Pairing.G1Point(uint256(0x10fbe50525ce345fd4a77911485349a090cd81eaac86523833ecc5b498dcfab2), uint256(0x244dbceeae7907f90b960695d05841c8b5805d4d025db8e56e5608ae3886597b));
        vk.gammaABC[127] = Pairing.G1Point(uint256(0x02288e241ef93ecb0615e0192b79f558a4ae2350cbd7970d4d9494a2a2e69ee9), uint256(0x11036aaa5b3968a0a342bf7ec9ca508ae1a3c23d666b756d9c7ee7fc761b6127));
        vk.gammaABC[128] = Pairing.G1Point(uint256(0x189cfc3552c8e4a4b3c8733e2d0c6a702a5d80cad3bf23497a90f8e58e41251a), uint256(0x0f06e55d5ed2ff45aac8d2e4df6964199b8c7f2198d0f058c49cbe233bc20978));
        vk.gammaABC[129] = Pairing.G1Point(uint256(0x0271710b8c97dd4baab7a8108a048fca4533fa9f057eaec5c009019f264391c6), uint256(0x1d9db51b5bc5ccf0fe705d5e91745fcacf8e0daf8dafaf007dc2b61e6df8d6f6));
        vk.gammaABC[130] = Pairing.G1Point(uint256(0x1a63ddde0779a264b9ce9aed91a32d69110e01265fb080de21c7050d6889cc09), uint256(0x1df5f150c5cc7953b24b2e5b1f10d331fdef864e7d00f0a576a8237af2c63e9f));
        vk.gammaABC[131] = Pairing.G1Point(uint256(0x1616b5d8b8170684636422c280a3ebd105eb767c22437e6607a9ef5308b8b9bd), uint256(0x2db9c350c9c636342dcd640475d1424987856251b66e2f7151b71cc668aae0f1));
        vk.gammaABC[132] = Pairing.G1Point(uint256(0x28d5d984cd816712f639c80ab85667d9aae7b4dd391ff427cf11c5880b606df1), uint256(0x0fb900d6c60e75c5efc5a9921614fc379cb42117ae5cbd294ba13a77dd62ff06));
        vk.gammaABC[133] = Pairing.G1Point(uint256(0x1762b01ca6100018de42f3c8fb93af793641145ff42b0454ee0d65fd1890b1c2), uint256(0x160ade7d57b57d5dd42ccc3453d9468fb6519ec3246b9457ccb17d4981530ba3));
        vk.gammaABC[134] = Pairing.G1Point(uint256(0x0cbef8047e8f5dac94e325868bdc5fa780842170761a3bab42dd41d3bd122725), uint256(0x2745a9fb7254c0afbb1e5a3118e7a236d7f0a6e6035d9714a00658fee544ca65));
        vk.gammaABC[135] = Pairing.G1Point(uint256(0x09d99693d80a93af9cefe458fc6ea8ad0c540fb0a4759b984446816167d8e0dd), uint256(0x0fae0aff68d880c9f9c703cbdb163a7c0cda3073ba15e3a3f6ef62b91f301053));
        vk.gammaABC[136] = Pairing.G1Point(uint256(0x11c2bfb02e6f9c26542a371e7a93c5f17989304dadfe4cdce6deabbf6cf6c907), uint256(0x04baab01df5cee8044a169fbb4c28b9ffeea6740b0ded6655fd8cbb33ad2d78e));
        vk.gammaABC[137] = Pairing.G1Point(uint256(0x1da5c9090a7b89e7e008c1ca9c9ef86efb51ccd8231cff16fa9130c26efa1abc), uint256(0x198c565e68f698bdefb03a2e906c27888fd531ab978c4610908a4e908ec8c9d6));
        vk.gammaABC[138] = Pairing.G1Point(uint256(0x2b83e2cc2f7898110cd354c8b4e4b70be3b70e60cb9a444add8e80ac5092c0aa), uint256(0x29383618d990f42c22cbc2f52c2faddb1ec383249b184e35cf01c78054c2250b));
        vk.gammaABC[139] = Pairing.G1Point(uint256(0x2e04d32e896b01d9a0510ad8b960db8237bf1101e1d096dc2c05eb5a9bffb6a3), uint256(0x08820e8b552795b27151e0038c61d503f313172484e26fe2db65dd104c47c8bb));
        vk.gammaABC[140] = Pairing.G1Point(uint256(0x0f0b2c3e7d42e1cf35a2d490e638a092290c3c8d49aba43de005102e3b1277a0), uint256(0x2bc3728c13cbf06a1a593f4d74cf7265becec0870e8a6c4aa573e876becca894));
        vk.gammaABC[141] = Pairing.G1Point(uint256(0x0a4b8efc5c7dc649ef8814d5ed7ddfdeb2963285aca90318f988fa1079a068c6), uint256(0x2885e9e532d4fdb78d9b1fbbd4b3db8681ac4dae188df7ecf7f9e28b260c0fe2));
        vk.gammaABC[142] = Pairing.G1Point(uint256(0x142a11a13f0f1140cc3c8061ae1e440a7717c235dca1db1a5fcff1678e80340f), uint256(0x0028e489a4d71419afde4bd132e22f2f5d749b522d24ca47357ab56b5756e5ac));
        vk.gammaABC[143] = Pairing.G1Point(uint256(0x25ee7444f23e0e2b6b653947d33f803182794a6d1173c7d03e2b96eb919f546a), uint256(0x30207b40c90afc6d460a076827fe420a510ba67a3e22eb1475b5b606a2b3c4f9));
        vk.gammaABC[144] = Pairing.G1Point(uint256(0x2e4213458608a0eaad01f3ddb171bf9ff3440dd968c04c37fc673c3c2a55f48f), uint256(0x1ea33e8af81b8b6378eb1e8f9519fec0389e2b53a8727080fb469135bc449d22));
        vk.gammaABC[145] = Pairing.G1Point(uint256(0x288ea856af66d444f11eabfb0776cb26e58a0fb08692d231a791671a591433d6), uint256(0x00ad868ea4f0814a3123336b2627d7ff16bb83e25b6d28bbfc64f49b6f1c3e15));
        vk.gammaABC[146] = Pairing.G1Point(uint256(0x2488de1d38f78acc819bc20913e2c3cf73b6354d77050148d36628d55ca8793a), uint256(0x1ce3df533b9cebe769608aae74ecf02857d8ab9cf68af839b7b1a40a45f4fc9c));
        vk.gammaABC[147] = Pairing.G1Point(uint256(0x16c809a06b04c41e7b77689844345e9df8e76d5fddeefd5715425f74e86fc5bc), uint256(0x12a8c7586ba11b4b9ce738bcb6fbd1151f473701c1fe25e46a8ed1df8235f513));
        vk.gammaABC[148] = Pairing.G1Point(uint256(0x2dd9330df472d5f9d7773cb575cc67f31bcf407e849202c6ce0c6a6416b63139), uint256(0x2206f777375132e6645f91edeef788e699cbea52310444753508f648412e03cf));
        vk.gammaABC[149] = Pairing.G1Point(uint256(0x0c9aaa4d64d736f4a38052ea04300d5174381cf886390fc81e7c312cd3a6940d), uint256(0x15323da1289702f79d83ae0c79c5394457a8cc570ade77a62fda8bde50712a52));
        vk.gammaABC[150] = Pairing.G1Point(uint256(0x154f60a91dd88a029e2199694c12193dac3918f192c083a8f908b4438e1a43f0), uint256(0x131ea0ae8994f147c8578b1f615be19c35a1a26ec40dcc0db8fcbfc8920982f2));
        vk.gammaABC[151] = Pairing.G1Point(uint256(0x0da5fbe0f6136aa4212a3985dcc8c43c9db252f7a9672c31b4a54aa09d9b8f4b), uint256(0x2f2ce254688c4179b1435094bd4e6fefc1ded4e2d5c562e3dbac62336c848b83));
        vk.gammaABC[152] = Pairing.G1Point(uint256(0x046b155bd59160ba6c8435d26ada38de0fe4b60ee02c9205c04e5a5d469e0181), uint256(0x115ab0b5d85bf0cb3a2c8428b992053e577b10fdaf8318c5725aa67db90ed1b8));
        vk.gammaABC[153] = Pairing.G1Point(uint256(0x1ba938caa171e1ae22dd0d6efaf35215bb8d5a953004f2b94b353fd2f6332466), uint256(0x05720a13e082a1109a29ad00a02ebb17e8959348cc625242238df0b54da69b30));
        vk.gammaABC[154] = Pairing.G1Point(uint256(0x2ddcab6fa67d879d15470e633ab26ab12f1450d5c19c5650aad2adb5ee08ff38), uint256(0x157c03b341722b60a68b1d2fe3414dda21e4261af7301d7f9d108b21b6e7653d));
        vk.gammaABC[155] = Pairing.G1Point(uint256(0x236b071ede9c3d89cda4ac6a879bd5382256f44066221bc0c68d7cfb95efbd56), uint256(0x2c3b5d0551f8584840df8b77e38fd8dc692f86f1516f020d191b46066ba9f56f));
        vk.gammaABC[156] = Pairing.G1Point(uint256(0x1396d6136f316972fdf6759f6708b91afc3bcc068536219929918e6f67f7161c), uint256(0x259d4fdd028ae2bb01b90b52f73eed5b9291509f938795399e66ec0b377fa5e9));
        vk.gammaABC[157] = Pairing.G1Point(uint256(0x27fc54f9685b792dbabcd5adc1d2fddca7a7590470925a9b820c2257032c20fb), uint256(0x1bbcb258dca40479f4fe9c5913d622e3e723aad31ee18432c198b69d1a0a221b));
        vk.gammaABC[158] = Pairing.G1Point(uint256(0x28c640349b6f9480d7300d637f7d5d2284665a8f31e97c71cb0ca0fab212d942), uint256(0x2355f74a588c91b4e6e720d267c811e5d515ea72e35a7715653abbc15aa1a5be));
        vk.gammaABC[159] = Pairing.G1Point(uint256(0x05ad51f33e2a7928fec8d12ecbaaad4f82cea47474476539985eccee8e3d0a1e), uint256(0x0443efde63215693fed5459754ca6c95d2dde2ca34e00078b9dfe6d07c50ac57));
        vk.gammaABC[160] = Pairing.G1Point(uint256(0x17688f6a4724bbeeb1aec67a696c2a8d0aa58539ad08de578badcc526a0464e0), uint256(0x284d79359453c1e124f8b79dd98f926bf698d4c0141734e8621ebe7c722d5140));
        vk.gammaABC[161] = Pairing.G1Point(uint256(0x101ca79e69789ea5648863538205e0e4322df1b9f68e89eb721b3decb1f9c0a1), uint256(0x26dd41a2d7ee691358abeaf675e795fa446fef249eee47f4aaf2fb0fcec2581b));
        vk.gammaABC[162] = Pairing.G1Point(uint256(0x150e48922d879f72bb309c7e71cfa5bcc270100cf4dcd5f5b8eeb13b96cffce4), uint256(0x2ba8a1ff814b3c65c860e7b06a856aaedeffe683ca6532d8685eb075b20e2ba8));
        vk.gammaABC[163] = Pairing.G1Point(uint256(0x301ce732a8e2ae43187d68cf3681e7d9014ab22068bed2867543115a26c28a3d), uint256(0x0a49396c8cd5ff7162d2e39b2474e4c1f13d1eea0869fa3fb7af959cad87869e));
        vk.gammaABC[164] = Pairing.G1Point(uint256(0x2ee6629084b2fdee2fbfa5f93ae5f2f1438f3c961b1bbef5b8897de7fd802e2f), uint256(0x17280327fd23c9c92ed133bb2208652955a0ba04f478e398f8552ed7b8de224f));
        vk.gammaABC[165] = Pairing.G1Point(uint256(0x01f308f2a1e3ce2c5bc5f37ced56b9a53806d598c9aaa5aa07683f067a1c7abf), uint256(0x1d5f068a2d234214619ee3808aefc9b27a8b8944fcc6f68ab07a49f5a302a343));
        vk.gammaABC[166] = Pairing.G1Point(uint256(0x05e44e4b64d39c554e29413a5f7c7febea4554ad895cc647948b54f6552838f4), uint256(0x0b117de6d56f7b47cd948e18da118cd890da596b0fbf77c7385364c6c73a292e));
        vk.gammaABC[167] = Pairing.G1Point(uint256(0x2c23b19ec95be93b4af7f7bbb06b7dc1eb19a0025d11023ded8d5b33697fa7de), uint256(0x2de7f7e2e302101c75ca55386bd1a02ed3b342c599cf7cc5d8b26ce99cd02608));
        vk.gammaABC[168] = Pairing.G1Point(uint256(0x2f27a445fe941012098c19f467c6e290f2d6ed27642c07af4f1076c5d25be73e), uint256(0x2c725a3220cdd20ce690c5672a9f4bab9fdbddc9e484b32aecfe4724c81dce97));
        vk.gammaABC[169] = Pairing.G1Point(uint256(0x0071971b20d4467d6fb7da0023fe4434cb3d3ce7f49a2d6e81ecf31dbe139349), uint256(0x2d6da7274299fe3e061753a6af5a83a833aad8f1e0393d4a88c939c5e3dc3b7b));
        vk.gammaABC[170] = Pairing.G1Point(uint256(0x292f9c3f7e06ec97c4e3be338d2dc1bcb53aa520ae62559f808ca4d5f23a42be), uint256(0x2a9055cf4b05a0bc7eb0aadef722783d6a7cebe298178a06aafe80e1133997c1));
        vk.gammaABC[171] = Pairing.G1Point(uint256(0x03f28e753873b5b58b743368f03156817d8e87cf4fe5b0265e39f261df3e9dc3), uint256(0x08eee7d95fa81a3d34ccd3e1854a4c12a8d7f962ff92e0c81e532cd54d5c52dc));
        vk.gammaABC[172] = Pairing.G1Point(uint256(0x209b714e8a7c03bed09f884a8364b8740150e03e76f9c1cfa3bd4a2c2ae357b9), uint256(0x1972c6c1a5c331b9bc892f01d782872e4f0963657ce9534931da2182413d391b));
        vk.gammaABC[173] = Pairing.G1Point(uint256(0x05c7667417556a3c1ded9197a831f6174e863af8220269a88d09d797b529b665), uint256(0x058deb70d6d1c74a74566af065dcfac4bb2c4ca5bb7da56508324d9736ad6e5f));
        vk.gammaABC[174] = Pairing.G1Point(uint256(0x05ddd87d13abda7ad09eedd7631a344c89a5f0db44a2e4f8d7fa4bab6d29c788), uint256(0x143fa5860a189fa6ddfdb10484590c6a894b1d860ab7118c95d42cf904973421));
        vk.gammaABC[175] = Pairing.G1Point(uint256(0x221be42e3c75eb478ffea68a4f9271e702e3cd62b0643c44327539f61dfc964f), uint256(0x2590295b45ee75efe3d52bdd7c26d498cc212394611b5aa89d54ba3990f12c54));
        vk.gammaABC[176] = Pairing.G1Point(uint256(0x0356775e6ecee3f3ce11dc4adb9f4a75d62f01f5e1ddb0d93e19849ac365da88), uint256(0x18adb227e831e407ed0dc4178ba36136af44cd053b02143fa3c94a453d727aa5));
        vk.gammaABC[177] = Pairing.G1Point(uint256(0x07d16c66aa16ccccbf00a23abc8f5a8d86fa63db3b53434a11178abb456acff8), uint256(0x1f77236de4d9a252d2a4e44effc4c9c2bbfbb9c0724450b764034e309e0188d7));
        vk.gammaABC[178] = Pairing.G1Point(uint256(0x18ce10c4e44b5d7c10aac5ec1173c596961dcafb43038de84f504c9c478102b5), uint256(0x120dd5501ab6afdd397e82eec0f2822d03fd06c8c93d4ace7211a80dcb39ed1f));
        vk.gammaABC[179] = Pairing.G1Point(uint256(0x256b137451c5efe32429e5466138113e813120e14710c46f55476f55b75ee3db), uint256(0x1309a7a7708142157ac10ac800b9334224bd060d12a7880aa53a8342b06c30bf));
        vk.gammaABC[180] = Pairing.G1Point(uint256(0x08d3607e22faa6067e161271e91c479aa6866dc865dd0a44023d5c8679bbd542), uint256(0x0e1374911be700375122e1c8d0c155a7e21f7a0b66ca9c53f008826a0998925d));
        vk.gammaABC[181] = Pairing.G1Point(uint256(0x16ce9acc1887d3217e72110bc407ca3d13a5ce64e112b03e48f39e7693bb49ad), uint256(0x28a004804ff2ad54865890640e8a1594d28b7e1734a2f81517a7be94e51e8258));
        vk.gammaABC[182] = Pairing.G1Point(uint256(0x0a2592faac20d7f19faa65a9d2379b4803c828e6af96e19f99b66a32a8829a80), uint256(0x20afdf8aeaf1e72ebd998a4272934c632b2e3ad6d33b067098da90d8d919a976));
        vk.gammaABC[183] = Pairing.G1Point(uint256(0x16146998495ab99d9022797084d4de6a95991fa071f6b7324837201f58e3c1c7), uint256(0x00e40abfcc2f1a142ea551e87035ef6ebf735de32b6ad633c24561a485e5bb8f));
        vk.gammaABC[184] = Pairing.G1Point(uint256(0x1fc0a4f72ccfbf5dcbacda522f6389a40faa9fca13478d547771744303bad00c), uint256(0x09f68ef3fdf6a0cc7f02d35247391084734c056b45395d49e146d6cfd6c8c097));
        vk.gammaABC[185] = Pairing.G1Point(uint256(0x0c48c9d10575f44e6213ec842f36e65f2298ffbec9c004059a2d5daf80a3f2de), uint256(0x05b804eaf1788cb0c9437eae515545bc152abbcd500e24c9375c01c9d6662cd4));
        vk.gammaABC[186] = Pairing.G1Point(uint256(0x2d03e9c03391ee3e089207b9843f06d5de193e0621bc14861ff029fec8f78033), uint256(0x0a518cf9254b614f92d7651d3ef7d9f00193bb3b3070062b1732c903d491cb17));
        vk.gammaABC[187] = Pairing.G1Point(uint256(0x296b42364a5ff124315fd32223f76f39b5825729277664d65ff618a86af2bdbd), uint256(0x1b6daaca4f026d3d98190b8ea5c2069c48e27a4622099b7b16f0c53ae489dbb5));
        vk.gammaABC[188] = Pairing.G1Point(uint256(0x164591a01a4562e01a55b604366d26dd65f70142a071f14d4791b2de477adf86), uint256(0x1f51d58e46a041eb015999c2e60f13b8680649ba18159b720ff7fca92239bb53));
        vk.gammaABC[189] = Pairing.G1Point(uint256(0x03019bc49d30ea560312452a156d75627e6f2e667f7865c435222750b31938cf), uint256(0x0e96bb1cb0b05813c5b219f4f243504f1c19af4a481b4c7c3810811cedaf123a));
        vk.gammaABC[190] = Pairing.G1Point(uint256(0x1328b75e441de3fdd18fe6282e7e552733eb82ee115e7f08028d3dd2c13f767d), uint256(0x0faaf7dae56a890dea5a165eca11e72de11f84ab2fe4c2bc04d8fbab09a01aaa));
        vk.gammaABC[191] = Pairing.G1Point(uint256(0x084968e7bbd98c11b6a0553556703468101f045ed5fbe71597f9c8063306cd98), uint256(0x2fe93c4e068959393b80d53210cc559b90277d005c31bcbb1acf4d75dde8a2ea));
        vk.gammaABC[192] = Pairing.G1Point(uint256(0x286e2ad55c57d43914642570287b0665ad4d67665d57624f7898cb366285a3e3), uint256(0x083ad59c1f6c1b3ee48fc57e922a78a5c539189b4f69217e75b1857aec017d4a));
        vk.gammaABC[193] = Pairing.G1Point(uint256(0x20b2d9d392f77636756295ebf0e0908afc557476da093349adacedc351d0f459), uint256(0x2ab495afebcc94d81389f47bb97e053f45ecef4850af4df2414381de4ee5ba8d));
        vk.gammaABC[194] = Pairing.G1Point(uint256(0x2ca42272ec4cd0924a512e98329a92f5224f797776982336e0320df1ca9d53d6), uint256(0x2fe6fe99636b6604ebf3418f15fab22e7fcd41ae72d7de1678209bdf7f63ab6d));
        vk.gammaABC[195] = Pairing.G1Point(uint256(0x2600e2721a5aa084a78fdf49a267fe3b00d03a21c9d405dcadbf50900f712b18), uint256(0x16122bb64d5c7139fabdb1b47ea82897b0705da182c2c19aeba0baf06def621b));
        vk.gammaABC[196] = Pairing.G1Point(uint256(0x2832567682fc4a2a0db771901b3ae71d6d805b65c8d5f497e828afa45a97af16), uint256(0x1dd91a4c936e2776d5f3759e7d261d29be2cd6f7ee7e1b22dce654c74791cec5));
        vk.gammaABC[197] = Pairing.G1Point(uint256(0x293ebf96fdca7e1ce2526aff60a5dfb508432f106943a858540f9357d180f60c), uint256(0x19aa9c64b0a2db270fec6c46569ebed3daf85d6dd653728a3af9e1d1891f832f));
        vk.gammaABC[198] = Pairing.G1Point(uint256(0x259a968e17d29a7cb10901958f8d329cb5f5bc52dbedb220b3167b91695d8b16), uint256(0x044b15cdf2c3bff7d9591e87b8d29d445e2c6490545df1d8299cac6a81c261de));
        vk.gammaABC[199] = Pairing.G1Point(uint256(0x2738fdffb8113cdb9bb49405df2d9c5e4b3ec7910b84641633ebb959981648f5), uint256(0x0e15f3eeb0225ed32b5b5e075394569b7c74f79367f5db6645b2d846c1c809c0));
        vk.gammaABC[200] = Pairing.G1Point(uint256(0x214164ff76c0d2834d68817bc55ef896af1bc25a7dbacb1d7881d85ed394583b), uint256(0x0e39c335a2eea16c5b9d7efb2e0663fc378831c0926d857c2ccb499aa070a61a));
        vk.gammaABC[201] = Pairing.G1Point(uint256(0x2b11e9569160c6f3234a5c8691fceb37ddf867a96dd400241091c9181cf61a5f), uint256(0x231719615ea431cbe92ed7ee635f35fda4259d7ba920ec372f15f438a9d60df6));
        vk.gammaABC[202] = Pairing.G1Point(uint256(0x2f6e2bb3357bbd9a47e73c3b1a7574a9b858fbde2559516c9ffda2f49d8d72cd), uint256(0x032d4645678fce85b281b5e733eac92f75c63c991410f96063740902e6d12ed7));
        vk.gammaABC[203] = Pairing.G1Point(uint256(0x07007e0dec272b2f2d4b53509218ca311a2806e85ed7fd7fb5c5bd2d86f5f61a), uint256(0x2754bc12b9d2c2c19a9bfe0bf306904501acb816a7e31c62eee90e152fa65ab9));
        vk.gammaABC[204] = Pairing.G1Point(uint256(0x070852e9395bd9a53faf6c40e57829f2fc6ed4eaed03828b1e66252eab2819ba), uint256(0x119fb05a86e4ea523f6e40dac414426af98e714e75145d9ecd8211aa9df77769));
        vk.gammaABC[205] = Pairing.G1Point(uint256(0x2d4a276e3cd92c884ffef08bca1dbf2fd475a14292d8d14c730950a867d89d74), uint256(0x1d6d1b24d77b9d23316afc32e8acca9590d0e82b2114837ccc10ff8924560137));
        vk.gammaABC[206] = Pairing.G1Point(uint256(0x026d1cd35aed87d5d04bda1414dd1853ed2e9ea8ab6e0a70765bf57a390e5160), uint256(0x02352ae94802ffa0efde8ce168104a7be72efd6e8ed80b8a45db6d4c42df9885));
        vk.gammaABC[207] = Pairing.G1Point(uint256(0x187029ecd7e5b66bc5309702190dd53a96d3ca38ebe154d63d65d9eed365ea3e), uint256(0x24f3decc26edcecbe6932ca51414d9373d2f135a8dbaf8d93ea2792354fe8938));
        vk.gammaABC[208] = Pairing.G1Point(uint256(0x13e400c553176325f6851c10ed040c8cf192f840f7524e2fe4160fe102af026a), uint256(0x0e20e4a3757519ab428750ffbd44bedb9efd697e2abec7c48f814b37fff8fc8c));
        vk.gammaABC[209] = Pairing.G1Point(uint256(0x30387df677d3118b44dbe9c909bded28b9dabd0fd48114da7e55b8c8e352a2f6), uint256(0x1691532e11068a17b042089ccfc1721c56d6977fe38c327ef45cc7b21cf5a93c));
        vk.gammaABC[210] = Pairing.G1Point(uint256(0x2709b970358e3c44ea749a845528827dd7172c4d1d5e7ea3002b241c03362b47), uint256(0x112f40d5f5a1acb2def8983d5ac7230f7e58e22e4d3d201bfe9c908c587292b5));
        vk.gammaABC[211] = Pairing.G1Point(uint256(0x115842cc1e2044e9e98c76316d568e36f60c2f29391584ac244fc184242a383e), uint256(0x16705b2e8d6b92ef1b72067a70373fa51374b6344495f1ac69343009ea8c63ee));
        vk.gammaABC[212] = Pairing.G1Point(uint256(0x2bcaa8a3ad1bc4be3d127412765461adda89599c2e06ec845cd8b6d2efbd6cf4), uint256(0x17b9169da17d22bef76479b8298bf3ae5976bfc345f9e0d1e2ccc85e6c65679a));
        vk.gammaABC[213] = Pairing.G1Point(uint256(0x06a6ef79b322f8b9a459162edeb7dafdc700fa7eaf05cc16c998c9d67d65c284), uint256(0x138f83a063b38ce196082b309096bb2d2e855ef9530138cedf0937a62f8ff3a7));
        vk.gammaABC[214] = Pairing.G1Point(uint256(0x29326bb8ac84d161a3c70766492c8fbe159750df7415101ad1afe42869fe18fb), uint256(0x2739a7be351e143b930d64f4b57f95070ecd1d7c34fcfc18895ac02ccd757c9b));
        vk.gammaABC[215] = Pairing.G1Point(uint256(0x269b9c3a14d5f644110ceb95d895a8c1935816f636d82308230e5d3a036756ac), uint256(0x2a6df9b5fab2d88721dd24fb4822566eddd692156eeb6f543c71458bff0e62d1));
        vk.gammaABC[216] = Pairing.G1Point(uint256(0x02cf1c2234aca51cb6deb94e12d096c88daabc5269bb888c873f05f98e71a9d1), uint256(0x0daeb1ec50297d5354a2302e3ad5d27e58074baf608f6c3c8bfaa22a7e3ef566));
        vk.gammaABC[217] = Pairing.G1Point(uint256(0x100ed1beb38ffe2d0427f7ff6b06f230483d3e4cfc3d10f774981baef7ba520f), uint256(0x0736e9acc40371ddfe4036d0e42ce812936322799e435106a5fd6902a23981ea));
        vk.gammaABC[218] = Pairing.G1Point(uint256(0x2c40a1f4ca9280a1837eb2b647d43479141fbeac7dc19cc92f25d8003c0d321b), uint256(0x21079526bcbb4a37d0cfed9ca26d29be5644b75563d2a6805aeb4d6d95b66413));
        vk.gammaABC[219] = Pairing.G1Point(uint256(0x1a3bf321178221b4a7d656c028e99299bf8e59aec4687769697418bfd33e5e3f), uint256(0x1d9948042431b8decbe9042a493535432583846df408e54eb107ead1db5dbc43));
        vk.gammaABC[220] = Pairing.G1Point(uint256(0x1c9531dbdd844925ef1a21c5ceec53fc464154332bb5e9cab2c855084834d265), uint256(0x0b9db4d118be43cc1e512f10992ae728460530c721336396af1d5dae6844d3d9));
        vk.gammaABC[221] = Pairing.G1Point(uint256(0x2888ad8d5c34ebb0bbf0d96521e2603a892a708577b7b198f4467f09b6be636f), uint256(0x02f8e661aa21d45eee90a2a4658aead4ef90a7cacc60d5521464f7edca16529e));
        vk.gammaABC[222] = Pairing.G1Point(uint256(0x29b3f735ea27c5955baa17eac5f97489315766d4eec9fe5b86510a476b29c32b), uint256(0x181f8061068d2eea1d5e7035f98f38236dabe35cc751d2fbe1ac4456a4ba116c));
        vk.gammaABC[223] = Pairing.G1Point(uint256(0x12557cc7d093e6ec266362d87b39122df3dea18fa62b07385fd629d05636c567), uint256(0x05b3478da872f0254ae1b4954347475b3e5c704d19897cf7f1c3c08ad62da831));
        vk.gammaABC[224] = Pairing.G1Point(uint256(0x2d7d015869b8f24526a38c3197f579910a5e80c0dc29a3a080153155ea4e8fa9), uint256(0x1b2da33235b4ab84fbc3801f2a690e0aab7d52759eaeeae558c0e4f45f0d489c));
        vk.gammaABC[225] = Pairing.G1Point(uint256(0x2fa2f4b4d3fdd5d5df86e3c41ca51d03d29347faa2f386b3e6a48743acbb1107), uint256(0x26504c71118f9b327dd79d01e6f594de74dd898399ce45a5d7da25a1c34df759));
        vk.gammaABC[226] = Pairing.G1Point(uint256(0x2af6c2c9a2e4080be61a7bb003cbdf24897d056cfdd6a882f12fc18b719bd051), uint256(0x19593b8f500fa4abfe05840f34c919ec3c181370e0b206439834005be99c839e));
        vk.gammaABC[227] = Pairing.G1Point(uint256(0x2b5fd3ef85ccf3ba7b535c9a27ff9d3562409843fcbf9df602a8902039ce390c), uint256(0x1b92e83ddd11c54d7af33bd0fe40d806470ce7a3e00035de22d39fefaa798125));
        vk.gammaABC[228] = Pairing.G1Point(uint256(0x0427799aa41edfe1ede504f3258cfabf5c3b92d73287281f735285411b380652), uint256(0x2a33bc31371962d2a3924f1d839f922042be094af5eaff778c41e9cba3daa441));
        vk.gammaABC[229] = Pairing.G1Point(uint256(0x2f2e5982d3c890d7edf5de4d665b3ba7b1f82d3a24acb1cea90f2a596f84a112), uint256(0x165530ec089dd63326b53387680d75ede454841ec7f8151444821f1e206f481f));
        vk.gammaABC[230] = Pairing.G1Point(uint256(0x1468e1b20ba9835bb5bc96e47633d5a1ce697953597875a78056eedfe424d343), uint256(0x031fabfe8e4c2ada438b84c80a0cc93cf4a493a09c2c8c2766bf11dec9031296));
        vk.gammaABC[231] = Pairing.G1Point(uint256(0x287a91515d48f10c81abdfb3b3e2ec1a119a08798ee295674d4b5f53f5d093ab), uint256(0x11e7684f7b36d89b4f94a901d53e5645894513ce8c60d868fe2125a920e820bd));
        vk.gammaABC[232] = Pairing.G1Point(uint256(0x0db7e3dca2d509d5feb6a689f2da2ae9b9075611323da878dd4ab1d951b4eff7), uint256(0x033a7e523ab66ae2ac863edc46673fa980da87ac4bcb036c7b26011a58e2d80a));
        vk.gammaABC[233] = Pairing.G1Point(uint256(0x1a5ecf2ccd2bb0df56c88c58b5e5ff240840add4868b2fac559bdafe0365a8cc), uint256(0x22d85b1afe5dabd8e6f0de84da9196703fdf419bf9590bcd030de15607ba30eb));
        vk.gammaABC[234] = Pairing.G1Point(uint256(0x0be508b14876897489d497dffacb0aeb028869d58893ec116049da2d5756ee58), uint256(0x1d94723df4139d850069466bb8f9c3da6ed50424a50804cd37d47f7f78367e8d));
        vk.gammaABC[235] = Pairing.G1Point(uint256(0x1b76ea03c35b51a8a8d2428b6fdfcd0384b421cea9123d8bb2b1f466c4832333), uint256(0x050b032b1548c0115cfea8d58b95b9b3bb9bc7dec9259f516cbb4a5ba94211a6));
        vk.gammaABC[236] = Pairing.G1Point(uint256(0x1511cfa9ba1a5741b5abad7a3456234cd6ed04ce8e4a5d5f4eae41aec003227f), uint256(0x0ba1fcd65ebcbedfc017d93f431568b893a26b291afee7ae3cc3f43db88a226b));
        vk.gammaABC[237] = Pairing.G1Point(uint256(0x03c3ab1be3d2b508214456bdc0fc131d06fc5f41d9ff65616066467d4da3d1fc), uint256(0x157ad8e1c6adc60e6f1284e7310af71bade57934a9809f087432b4a71c1631f2));
        vk.gammaABC[238] = Pairing.G1Point(uint256(0x2ec84e2c355d608f4531f0927b4b02dfa8ceaf644e7208fdcd4369e67037a39f), uint256(0x28d0293e7bf2f7fcbd98abcc9ecc185eae5397ec0f50aa0da434aabe6fc65da3));
        vk.gammaABC[239] = Pairing.G1Point(uint256(0x0131271325cfb0815557d0d69e3a642b52a86e929b850e82c8905ad320d98a25), uint256(0x22dfbbe073dfa7373f241c80c824927d4c60909f04b798c5fc3e59199f14cf2b));
        vk.gammaABC[240] = Pairing.G1Point(uint256(0x0e7a604c6f50e19f541891ca77faa19ae257a0730f833f733beb826b6318a2aa), uint256(0x01bb780af7aa567d62945b9594afe3cd4344e2c24cb1c5b33c2d417db6b51030));
        vk.gammaABC[241] = Pairing.G1Point(uint256(0x2a98faff4522349a56fe1f71aa61965ac5f2d39c02baef4198d57e24a01f15d2), uint256(0x24572e00676f851b2c4204fc21173816437bfba2fbe3febac44325f777182abd));
        vk.gammaABC[242] = Pairing.G1Point(uint256(0x280708205db0ea37436c14a7c14f259192d8c7d3943231021f450b57322043fb), uint256(0x18e22190acf325c8889326bab4759e968d723a558c7a44a08a71d749760a86a3));
        vk.gammaABC[243] = Pairing.G1Point(uint256(0x03bd20cb0cfb997a777ebc086ba3be8476f2e3cc8f3eda48a43a410fbbb919fc), uint256(0x28eb00e88fced70c4993802c5da8e7208562dddc8b6713210609af5b164caeaa));
        vk.gammaABC[244] = Pairing.G1Point(uint256(0x0c7cf4e79b05cac76e3d66c72771bc3976431202895530a4190b66231ec45ae5), uint256(0x177c06313f24a5f6c13b3d6cc2dca27727c5035038b8effda9bffdf210c3a726));
        vk.gammaABC[245] = Pairing.G1Point(uint256(0x245c04101321bafbaeadb8048542c5c1cfcbb10192afb1f92a3bcac5cd4fb9e3), uint256(0x0f291a70f1fafd9d093739f317e80b98816bcf7d1e3c5bf3a18d25ea1b1d6d82));
        vk.gammaABC[246] = Pairing.G1Point(uint256(0x0c7330b5d2ef7da8f5071a81ceab777c35c4ece978e4f877c861cc46cb8b715f), uint256(0x1e9d3cb70c8e7e8127b10abf6432c4d89e806527967f5808b722849b81462a8f));
        vk.gammaABC[247] = Pairing.G1Point(uint256(0x116c9972a22fb8685a31d47e97cd62cf02f9d2553d5400e03f7aac1dd09fb349), uint256(0x15eeedc11cc038bed5d53297a83f6d7750e00ba5a0eb82f27105fe7191117591));
        vk.gammaABC[248] = Pairing.G1Point(uint256(0x10139936b87a25e884f48a21e6d52a675d5fc7f12b16374eb008d3229ac9d04d), uint256(0x14650551ec13184b3becd293d5df65c012e56284057e0d3338ddf645b8c87241));
        vk.gammaABC[249] = Pairing.G1Point(uint256(0x121974cdf016f7656345f14dc438b74814b2afba8d8d862de6c5cc982aa1390b), uint256(0x0df5cf2d1bf7c94634e74879941d2c8f564c20ae496dfe4f5b96fd7e5675967c));
        vk.gammaABC[250] = Pairing.G1Point(uint256(0x19ec1ca5d58e1847e40165216aec1cf1e5b91600b03732dd0287d2fc8890fe26), uint256(0x05f068944d1297f96c7834ee9b1d9f444db3951b00e003a8edc2fb3115c93105));
        vk.gammaABC[251] = Pairing.G1Point(uint256(0x12ee68a4842f5f35773fe355475fc5d9fe5450b43a7a6294b020257ad2ef41dd), uint256(0x2637d7b2b1998a47fc0d5b7a8e55bd3230116d112c2167acf31038df6bc3dfe4));
        vk.gammaABC[252] = Pairing.G1Point(uint256(0x1db85e0d35dd363217f6006fb89090bca77636e7b6e9e7a8247cc79472f42f26), uint256(0x296cc6c136280a963db858c2ffb8da38240949fa8b670d933db197ccdb7b2b25));
        vk.gammaABC[253] = Pairing.G1Point(uint256(0x0f8920f21f7ace64de2556133766a3a52f0bdae91493e82d67f73f44f5df5962), uint256(0x00cf84154c0b9c0a2efbadf2b1827943c9cfc164389fe71aeb3054a640d044f1));
        vk.gammaABC[254] = Pairing.G1Point(uint256(0x273f0361effd99f78ff0f752bb08c6f41361e261396a381a8ab9c90ec6d97667), uint256(0x195f73ccbdbb645ac4645bf0a9ecdb75bb3dc85115d7a757095b2e610d29118f));
        vk.gammaABC[255] = Pairing.G1Point(uint256(0x14028586f06438e40a52c1c4781fce8d0f67035d5bb59136894d505d5440bd31), uint256(0x23ce6845d6e76677f1e5ec352535f2a296dd27574fd125160b80bfcc9d6037fe));
        vk.gammaABC[256] = Pairing.G1Point(uint256(0x2bbb3e22deb3a82c2e33d4c83d9583cf5a5b392182e107e7f8a8e22810c78871), uint256(0x2ea59d9934a899cae7ddaeb83142fabc67765b4d4dec32096bf1bdc1f2d24157));
        vk.gammaABC[257] = Pairing.G1Point(uint256(0x0286b627155bb2d53dfec9589a11c77523bce07b476e9be522b720ecb97c9d34), uint256(0x2cb54c5beecdb112e3ea96555761ee04a026f67c89e2bc5018773ac7de7cf085));
        vk.gammaABC[258] = Pairing.G1Point(uint256(0x032ce795d3e8f7776982b7ef621d2842623df268af128b8d401566f11d81cc15), uint256(0x128a6df625e4a4c98898bb6017cd43f5871a6560bd41cc09b48e4ee2dcd8ed8c));
        vk.gammaABC[259] = Pairing.G1Point(uint256(0x2619fd6cd42b2b896120c88a66705c3253892e5b6b3296e9eedf9072eda1555d), uint256(0x0c56e530dd77c9143c2dcd898de0925016bc079aa37d9c9b032e69a3e528edf6));
        vk.gammaABC[260] = Pairing.G1Point(uint256(0x18ce97998b0ca0b845c71ab21a9c6062cb18231c22463e6e7e9b2d1a4c3e8b7a), uint256(0x222030a6a0af680e25b074ead5f6ad412a63b6d9f7b260cb65491766c933a9e0));
        vk.gammaABC[261] = Pairing.G1Point(uint256(0x1fa7596da237bc1c5ba61ea94f5c88441282bb7f7cf12d2ca38d6061fb4c53f8), uint256(0x13e0237a226592060f2e30acd57cdc895414908d25f8b10c080215992f2c3ae9));
        vk.gammaABC[262] = Pairing.G1Point(uint256(0x221a8dced5444512721dd98d5110d9f9d0bc16449ec71b75a92c812023c0600c), uint256(0x04b8b03f4ab2bb85242c2dad99694f75870fae9e6744094d87866c0896b1d13a));
        vk.gammaABC[263] = Pairing.G1Point(uint256(0x17d08795733f13a152c8078a45a21f6cc977234063b2f74002a9da256a36d63c), uint256(0x03856723d4387872af472bb8e616c07f58fbd91178ec1ff9f4299ff18039fe60));
    }
    function verify(uint[] memory input, Proof memory proof) internal returns (uint) {
        VerifyingKey memory vk = verifyingKey();
        require(input.length + 1 == vk.gammaABC.length);
        // Compute the linear combination vk_x
        Pairing.G1Point memory vk_x = Pairing.G1Point(0, 0);
        for (uint i = 0; i < input.length; i++)
            vk_x = Pairing.addition(vk_x, Pairing.scalar_mul(vk.gammaABC[i + 1], input[i]));
        vk_x = Pairing.addition(vk_x, vk.gammaABC[0]);
        if(!Pairing.pairingProd4(
             proof.A, proof.B,
             Pairing.negate(vk_x), vk.gamma,
             Pairing.negate(proof.C), vk.delta,
             Pairing.negate(vk.a), vk.b)) return 1;
        return 0;
    }
    event Verified(string s);
    function verifyTx(
            uint[2] memory a,
            uint[2][2] memory b,
            uint[2] memory c,
            uint[263] memory input
        ) public returns (bool r) {
        Proof memory proof;
        proof.A = Pairing.G1Point(a[0], a[1]);
        proof.B = Pairing.G2Point([b[0][0], b[0][1]], [b[1][0], b[1][1]]);
        proof.C = Pairing.G1Point(c[0], c[1]);
        uint[] memory inputValues = new uint[](input.length);
        for(uint i = 0; i < input.length; i++){
            inputValues[i] = input[i];
        }
        if (verify(inputValues, proof) == 0) {
            emit Verified("Transaction successfully verified.");
            return true;
        } else {
            return false;
        }
    }
}
