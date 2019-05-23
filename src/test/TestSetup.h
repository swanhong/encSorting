#ifndef TESTSETUP_H_
#define TESTSETUP_H_

#include "../../HEAAN/src/HEAAN.h"
#include "../PrintUtils.h"

class TestSetup {
private:
public:
    static void runSetup(Parameter param);
    static string genSecretKeyPath(Parameter param);
};

#endif // !TESTSETUP_H_