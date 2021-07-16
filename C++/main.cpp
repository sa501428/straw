/*
  The MIT License (MIT)

  Copyright (c) 2011-2016 Broad Institute, Aiden Lab

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
*/
#include <iostream>
#include <string>
#include "straw.h"
using namespace std;

int main(int argc, char *argv[])
{
    if (argc != 3) {
        cerr << "Incorrect arguments" << endl;
        cerr << "Usage: straw <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG> <binsize>" << endl;
        cerr << "Usage: straw <oe> <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG> <binsize>" << endl;
        exit(1);
    }
    string matrixType = "observed";
    string norm = "NONE";
    string fname = argv[1];
    string unit = "BP";
    string size = argv[2];
    int32_t binsize = stoi(size);
    vector<chromosome> chroms = getChromosomes(fname);

    for(int ci = 0; ci < chroms.size(); ci++){
        printf("%s\n", chroms[ci].name.c_str());
    }
    int starting_index = 0;
    /*
    for(int ci = starting_index; ci < chroms.size(); ci++){
        for(int cj = ci; cj < chroms.size(); cj++){
            vector<contactRecord> records = straw(matrixType, norm, fname, chroms[ci].name, chroms[cj].name, unit, binsize);
            for (int z = 0; z < records.size(); z++) {
                printf("%s\t%d\t%s\t%d\t%.14g\n", chroms[ci].name.c_str(), records[z].binX,
                       chroms[cj].name.c_str(), records[z].binY, records[z].counts);
            }
        }
    }
     *
     */
}
