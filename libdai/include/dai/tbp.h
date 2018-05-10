#include <iostream>
#include <cfloat>
#include <sstream>
#include <dai/alldai.h>
#include <Eigen/Dense>

#ifndef __defined_tbp_h
#define __defined_tbp_h

namespace dai {

    std::vector<Eigen::VectorXd> tbp_marg(std::vector<DFactor> factors, int k);
    std::vector<dai::DFactor> read_dfg_from_stream(std::istream& is);

}

#endif



