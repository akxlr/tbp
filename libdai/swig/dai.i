/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  Copyright (c) 2006-2011, The libDAI authors. All rights reserved.
 *
 *  Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
 */


%module dai

#include <fstream>

%include "std_string.i"
%include "std_vector.i"
%include "exception.i"
%include "std_pair.i"
%template(IntVector) std::vector<size_t>;
//%include "std_set.i"  /* for python */
%include "../include/dai/dai_config.h"

%{
#include "../include/dai/alldai.h"

using namespace dai;
%}

// ************************************************************************************************
%include "../include/dai/util.h"

// ************************************************************************************************
%ignore dai::Var::label() const;
%ignore dai::Var::states() const;
%ignore operator<<(std::ostream&, const Var&);
%include "../include/dai/var.h"
%extend dai::Var {
    inline const char* __str__() const { return (*self).toString().c_str(); }  /* for python */
    inline std::string __str() const { return (*self).toString(); }  /* for octave */
};

// ************************************************************************************************
%ignore operator<<(std::ostream&, const SmallSet&);
%rename(__eq__) operator==(const SmallSet&, const SmallSet&); /* for python */
%rename(__ne__) operator!=(const SmallSet&, const SmallSet&); /* for python */
%rename(__lt__) operator<(const SmallSet&, const SmallSet&);  /* for python */
%include "../include/dai/smallset.h"
%template(SmallSetVar) dai::SmallSet< dai::Var >;
%extend dai::SmallSet<dai::Var> {
    inline const char* __str__() const { return (*self).toString().c_str(); }  /* for python */
    inline std::string __str() const { return (*self).toString(); }  /* for octave */
}

// ************************************************************************************************
%ignore operator<<(std::ostream&, const VarSet&);
%include "../include/dai/varset.h"
%extend dai::VarSet {
    inline void append(const dai::Var &v) { (*self) |= v; }   /* for python, octave */
    inline const char* __str__() const { return (*self).toString().c_str(); }  /* for python */
    inline std::string __str() const { return (*self).toString(); }  /* for octave */
};

// ************************************************************************************************
%ignore dai::TProb::operator[];
%include "../include/dai/prob.h"
%template(Prob) dai::TProb<dai::Real>;
%extend dai::TProb<dai::Real> {
    inline dai::Real __getitem__(int i) const {return (*self).get(i);} /* for python */
    inline void __setitem__(int i,dai::Real d) {(*self).set(i,d);}   /* for python */
    inline dai::Real __paren__(int i) const {return (*self).get(i);}     /* for octave */
    inline void __paren_asgn__(int i,dai::Real d) {(*self).set(i,d);}  /* for octave */
    inline const char* __str__() const { return (*self).toString().c_str(); }  /* for python */
    inline std::string __str() const { return (*self).toString(); }  /* for octave */
};

// ************************************************************************************************
%ignore dai::TFactor::operator[];
%include "../include/dai/factor.h"
%extend dai::TFactor<dai::Real> {
    inline dai::Real __getitem__(int i) const {return (*self).get(i);} /* for python */
    inline void __setitem__(int i,dai::Real d) {(*self).set(i,d);}   /* for python */
    inline dai::Real __paren__(int i) const {return (*self).get(i);}     /* for octave */
    inline void __paren_asgn__(int i,dai::Real d) {(*self).set(i,d);}  /* for octave */
    inline const char* __str__() const { return (*self).toString().c_str(); }  /* for python */
    inline std::string __str() const { return (*self).toString(); }  /* for octave */
};
%template(Factor) dai::TFactor<dai::Real>;
%inline %{
typedef std::vector<dai::Factor> VecFactor;
typedef std::vector<VecFactor> VecVecFactor;
%}
%template(VecFactor) std::vector<dai::Factor>;
%template(VecVecFactor) std::vector<VecFactor>;

// ************************************************************************************************
%ignore operator<<(std::ostream&, const GraphAL&);
%rename(toInt) dai::Neighbor::operator size_t() const;
%include "../include/dai/graph.h"
%extend dai::GraphAL {
    inline const char* __str__() const { return (*self).toString().c_str(); }  /* for python */
    inline std::string __str() const { return (*self).toString(); }  /* for octave */
}

// ************************************************************************************************
%ignore operator<<(std::ostream&, const BipartiteGraph&);
%include "../include/dai/bipgraph.h"
%extend dai::BipartiteGraph {
    inline const char* __str__() const { return (*self).toString().c_str(); }  /* for python */
    inline std::string __str() const { return (*self).toString(); }  /* for octave */
}

// ************************************************************************************************
%ignore operator<<(std::ostream&, const FactorGraph&);
%ignore operator>>(std::istream&,FactorGraph&);
%include "../include/dai/factorgraph.h"
%extend dai::FactorGraph {
    inline const char* __str__() const { return (*self).toString().c_str(); }  /* for python */
    inline std::string __str() const { return (*self).toString(); }  /* for octave */
    inline void printDot(const std::string& fname) const {
      std::ofstream o(fname.c_str(), std::ofstream::out);
      self->printDot(o);
      o.close();
    }
}

// ************************************************************************************************
%ignore operator<<(std::ostream&, const RegionGraph&);
%include "../include/dai/regiongraph.h"
%extend dai::RegionGraph {
    inline const char* __str__() const { return (*self).toString().c_str(); }  /* for python */
    inline std::string __str() const { return (*self).toString(); }  /* for octave */
}

// ************************************************************************************************
%template(PairIntBigInt) std::pair<size_t,dai::BigInt>;
%constant dai::greedyVariableElimination::eliminationCostFunction eliminationCost_MinNeighbors;
%constant dai::greedyVariableElimination::eliminationCostFunction eliminationCost_MinWeight;
%constant dai::greedyVariableElimination::eliminationCostFunction eliminationCost_MinFill;
%constant dai::greedyVariableElimination::eliminationCostFunction eliminationCost_WeightedMinFill;
%ignore dai::eliminationCost_MinNeighbors;
%ignore dai::eliminationCost_MinWeight;
%ignore dai::eliminationCost_MinFill;
%ignore dai::eliminationCost_WeightedMinFill;
%ignore operator<<(std::ostream&, const ClusterGraph&);
%include "../include/dai/clustergraph.h"

// ************************************************************************************************
//%ignore operator<<(std::ostream&, const CobwebGraph&);
//%include "../include/dai/cobwebgraph.h"
//%extend dai::CobwebGraph {
//    inline const char* __str__() const { return (*self).toString().c_str(); }  /* for python */
//    inline std::string __str() const { return (*self).toString(); }  /* for octave */
//}
//TODO fix problems with CobwebGraph

// ************************************************************************************************
%ignore operator<<(std::ostream&, const PropertySet&);
%ignore operator>>(std::istream&,PropertySet&);
%include "../include/dai/properties.h"
%extend dai::PropertySet {
    inline void __setitem__(char *name, char *val) {
        self->set(std::string(name), std::string(val));
    }
    inline const char* __str__() const { return (*self).toString().c_str(); }  /* for python */
    inline std::string __str() const { return (*self).toString(); }  /* for octave */
}

// ************************************************************************************************
%ignore dai::IndexFor::operator++;
%rename(toInt) dai::IndexFor::operator size_t() const;
%ignore dai::Permute::operator[];
%ignore dai::multifor::operator++;
%ignore dai::multifor::operator[];
%rename(toInt) dai::multifor::operator size_t() const;
%ignore dai::State::operator++;
%rename(toInt) dai::State::operator size_t() const;
%ignore dai::State::operator const std::map<Var,size_t>&;
%include "../include/dai/index.h"
%extend dai::IndexFor {
    inline void next() { return (*self)++; }
};
%extend dai::Permute {
    inline size_t __getitem__(int i) const {return (*self)[i];} /* for python */
    inline size_t __paren__(int i) const {return (*self)[i];}   /* for octave */
};
%extend dai::multifor {
    inline void next() { return (*self)++; }
    inline size_t __getitem__(int i) const {return (*self)[i];} /* for python */
    inline size_t __paren__(int i) const {return (*self)[i];}   /* for octave */
};
%extend dai::State {
    inline void next() { return (*self)++; }
};

// ************************************************************************************************
%include "../include/dai/daialg.h"
//TODO: why do the following lines not work?
//%template(DAIAlgFG) dai::DAIAlg<dai::FactorGraph>;
//%template(DAIAlgRG) dai::DAIAlg<dai::RegionGraph>;
//%template(DAIAlgCG) dai::DAIAlg<dai::CobwebGraph>;

// ************************************************************************************************
%include "../include/dai/alldai.h"

// ************************************************************************************************
%ignore dai::BP::operator=;
%include "../include/dai/bp.h"

// ************************************************************************************************
%include "../include/dai/fbp.h"

// ************************************************************************************************
%include "../include/dai/trwbp.h"

// ************************************************************************************************
%include "../include/dai/mf.h"

// ************************************************************************************************
%include "../include/dai/hak.h"

// ************************************************************************************************
%include "../include/dai/lc.h"

// ************************************************************************************************
%include "../include/dai/jtree.h"

// ************************************************************************************************
%ignore dai::TreeEP::operator=;
%include "../include/dai/treeep.h"

// ************************************************************************************************
%include "../include/dai/mr.h"

// ************************************************************************************************
%include "../include/dai/gibbs.h"

// ************************************************************************************************
%include "../include/dai/cbp.h"

// ************************************************************************************************
%include "../include/dai/decmap.h"

// ************************************************************************************************
%include "../include/dai/glc.h"
