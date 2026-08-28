#pragma once
#include "hic_slice.h"

int subsampleMain(int argc, char* argv[]);
void dumpHbs(const std::string& input, const std::string& output, int32_t resolution, ContactFilter filter);
