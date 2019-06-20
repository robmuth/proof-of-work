var Verifier = artifacts.require("./verifier.sol");
var ProofOfWork = artifacts.require("./ProofOfWork.sol");

module.exports = async (deployer, network, accounts) => {
 let deplyVerifier = await Verifier.deployed();

 return deployer.deploy(ProofOfWork, deplyVerifier.address)
 	.then((instance) => {
 		instance.addPublicKey("0x0a5ec76e95922ec59c8fde1a370de81f0b41b35982af6aa77aecede99ae1ee49", {from: accounts[0]});
 	});
};
