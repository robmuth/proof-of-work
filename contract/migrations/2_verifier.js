var Verifier = artifacts.require("./verifier.sol");
var ProofOfWork = artifacts.require("./ProofOfWork.sol");

module.exports = function(deployer) {
 return deployer.deploy(Verifier).then((instance) => {
  	return deployer.deploy(ProofOfWork, instance.address);
  }).then(async (pow) => {
  	let p = {
        "proof": {
            "a": ["0x1eeb6f9745e1b6918a2a16875f1abd09d9fa06f625dc69b6f804ffccaef1beec", "0x25ca32746bf354a2209cca02a5f6382c829e89682ead16759e32362ae207f3d4"],
            "b": [["0x2facd377c516e2d734bd88f3779000d19544ada48a64ec6858355730f7279e90", "0x09918d5507cfde080b969c5658977e9e7bedd6f1bddefddef5b8ce00488e5f04"], ["0x1fbba8c2ea2dd8de9335fb0fd27dec2027d92996ad7a426ecfad1c9f22f59cf4", "0x06bf073b7c91ede4f6e4c14b0105c1ae1036a02f3ad3639aa988ed9e16299473"]],
            "c": ["0x09039de974a78906075748a73115c19815ab359da6d7abc6c466de42b3f2df27", "0x00f0fde125b4d21070932c3d50f5c9a4813ab991f92c953e76361d51e364d14d"]
        },
        "inputs": ["0x0000000000000000000000000000000000000000000000000000000000000001", "0x0000000000000000000000000000000000000000000000000000000001254ce5", "0x0000000000000000000000000000000000000000000000000000000000000001", "0x0000000000000000000000000000000000000000000000000000000000e18a75", "0x0000000000000000000000000000000000000000000000000000000000000d68"]
    };

    await pow.addHours(p.proof.a, p.proof.b, p.proof.c, p.inputs);
    let added = await pow.addHours.call(p.proof.a, p.proof.b, p.proof.c, p.inputs);
    console.log("Added: " + JSON.stringify(added));
    return pow;
  });
};
