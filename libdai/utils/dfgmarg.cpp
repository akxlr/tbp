/*
 * A simple command line interface to run TBP on .dfg files. See https://github.com/akxlr/tbp.
 *
 * Author: Andrew Wrigley, National University of Singapore and Australian National University
 */



#include <iostream>
#include <fstream>
#include <cstdlib>
#include <dai/alldai.h>
#include <dai/tbp.h>


using namespace std;
using namespace dai;

void show_usage() {
    cout << "Get approximate marginals for a decomposed factor graph using TBP" << endl;
    cout << "(see https://github.com/akxlr/tbp)." << endl << endl;
    cout << "Usage:" << endl;
    cout << "./dfgmarg graph.dfg K" << endl;
    cout << "    Run TBP on graph.dfg with sample size K" << endl;
    cout << "cat file.dfg | dfgmarg K" << endl;
    cout << "    Read graph from stdin and run TBP with sample size K" << endl;
}

bool is_int(char *str) {
    while (*str)
        if (!isdigit(*str++))
            return false;
    return true;
}

int main( int argc, char *argv[] ) {

    istream *in;
    ifstream ifn;
    int K;

    if ( argc == 3 ) {
        // An input file has been passed in the command line, read from this
        ifn.open( argv[1] );
        if( ifn.is_open() ) {
            in = &ifn;
            K = fromString<size_t>( argv[2] );
        } else {
            DAI_THROWE(CANNOT_READ_FILE, "Cannot read from file " + std::string(argv[1]));
        }
    } else if ( argc == 2) {
        // Read from stdin
        in = &cin;
        if ( is_int( argv[1] )) {
            K = fromString<size_t>( argv[1] );
        } else {
            show_usage();
            cout << endl << "ERROR: Supplied value for K is not an int: " << std::string(argv[1]) << endl << endl;
            return 1;
        }
    } else {
        show_usage();
        return 1;
    }

    // Get marginals
    vector<DFactor> dfg = read_dfg_from_stream( *in );
    vector<Eigen::VectorXd> marg = tbp_marg(dfg, K);

    // Output in .MAR format
    cout << marg.size();
    for (int i=0; i<marg.size(); ++i) {
        cout << " " << marg[i].size();
        for (int j=0; j<marg[i].size(); ++j) {
            cout << " " << marg[i][j];
        }
    }
    cout << endl;

    return 0;
}

