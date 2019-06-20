pragma solidity >=0.5.0 <0.6.0;

import "./verifier.sol";

contract ProofOfWork {
	event Worked(
		uint employee,
		uint secs
	);

	// Zokrates Verifies
	Verifier private verifier;

	// Employers address
	address payable private employer;

	mapping(uint => uint) private workedSecs;

	mapping(uint => mapping(uint => bool)) private usedNonces;

	mapping(uint => bool) private acceptecPublicKeys;

	constructor(Verifier _verifier) public {
		verifier = _verifier;
		employer = msg.sender;
	}

	function addHours(uint[2] memory a, uint[2][2] memory b, uint[2] memory c, uint[7] memory input) public {
		// Check ZKP
		bool zkpVerify = verifier.verifyTx(a, b, c, input);
		require(zkpVerify, "ZKP failed.");

		// Get ZKP return values
		uint employeeId = input[2];
		uint pubKeyHash = input[input.length - 3];
		uint nonceHash = input[input.length - 2];
		uint workedInSecs = input[input.length - 1];

		// Check whether nonce was used before
		require(usedNonces[employeeId][nonceHash] != true, "Nonce already used.");
		usedNonces[employeeId][nonceHash] = true;

		// Check public key
		require(acceptecPublicKeys[pubKeyHash] == true, "Public key not accepted.");

		// Add working seconds
		workedSecs[employeeId] += workedInSecs;

		emit Worked(employeeId, workedInSecs);
	}

	function addPublicKey(uint _publicKey) public {
		require(msg.sender == employer, "Only the employer can add keys.");
		acceptecPublicKeys[_publicKey] = true;
	}

	function () external {
		revert();
	}

	function shutdown() public {
		require(msg.sender == employer);
		selfdestruct(employer);
	}
}