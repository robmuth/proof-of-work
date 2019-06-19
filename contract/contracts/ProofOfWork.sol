pragma solidity >=0.5.0 <0.6.0;

import "./verifier.sol";

contract ProofOfWork {
	event Worked(uint secs);

	Verifier verifier;

	constructor(Verifier _verifier) public {
		verifier = _verifier;
	}

	function addHours(uint[2] memory a,
            uint[2][2] memory b,
            uint[2] memory c,
            uint[5] memory input) public returns (bool success) {
		if(verifier.verifyTx(a, b, c, input)) {
			success = true;
			emit Worked(input[input.length - 1]);
		} else {
			success = false;
		}
	}
}