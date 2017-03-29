#ifndef MEM_USAGE
#define MEM_USAGE

#include "sys/sysinfo.h"
#include "sys/types.h"
#include "sys/times.h"
#include "sys/vtimes.h"

// extern static unsigned long long lastTotalUser, lastTotalUserLow, lastTotalSys, lastTotalIdle;
// extern static clock_t lastCPU, lastSysCPU, lastUserCPU;
// extern static int numProcessors;

int parseLine(char*);
int getValue();

void initCPUCurrent();
double CPUCurrent();
void initCPUCurrentProcess();
double	CPUCurrentProcess();
void GetMemoryUsage();


#endif