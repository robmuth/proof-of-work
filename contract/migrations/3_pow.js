var Verifier = artifacts.require("./verifier.sol");
var ProofOfWork = artifacts.require("./ProofOfWork.sol");

module.exports = async (deployer) => {
 let deplyVerifier = await Verifier.deployed();

 return deployer.deploy(ProofOfWork, deplyVerifier.address);
};
