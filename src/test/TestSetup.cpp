#include "TestSetup.h"

void TestSetup::runSetup(Parameter param) {
    srand(time(NULL));
	SetNumThreads(16);
	PrintUtils::parameter(param, "Serialize SecretKey");

    Ring ring(param.logN, param.logQ);
    SecretKey secretKey(ring);

    string path = genSecretKeyPath(param);

    SerializationUtils::writeSecretKey(secretKey, path);
}

string TestSetup::genSecretKeyPath(Parameter param) {
    string path = "data/sk";
    path += "_" + to_string(param.logN);
    path += "_" + to_string(param.logQ);
    path += ".txt";
    return path;
}