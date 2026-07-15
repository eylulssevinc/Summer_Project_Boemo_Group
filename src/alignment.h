//----------------------------------------------------------
// Copyright 2019-2020 University of Oxford
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#ifndef ALIGN_H
#define ALIGN_H

#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <limits>
#include "../fast5/include/fast5.hpp"
#include <memory>
#include <utility>
#include "reads.h"

//which recovered tail source(s), if any, get printed for a read - see the
//"tailSignal" writeup in alignment.cpp for what Source 1 and Source 2 mean.
//NONE is the default/safe choice everywhere except explicit user request via
//DNAscent align's --tail-source flag; detect/trainCNN always pass NONE since
//they have no use for this signal.
enum class TailCaptureMode {
	NONE,
	SOURCE1,
	SOURCE2,
	BOTH
};

int align_main( int argc, char** argv );
void eventalign( DNAscent::read &, unsigned int, TailCaptureMode);
bool hasSoftClipAtReferenceThreePrime( bam1_t * );

#endif
