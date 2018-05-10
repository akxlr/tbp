/*
 * This file is essentially an implementation of upstream factor.h with decomposed factors and approximate
 * multiplication instead of tables and exact multiplication.
 *  - DFACTOR_NOT_IMPL means the method hasn't been implemented yet.
 *  - DFACTOR_CANNOT_IMPL means the method can't be implemented due to the decomposed representation.
 * USED_BY_JTREE shows which functions are invoked when running libDAI's junction tree algorithm.
 *
 * Author: Andrew Wrigley, National University of Singapore and Australian National University
 */

#include <dai/dfactor.h>
#include <cassert>

namespace dai {

using std::cout;
using std::endl;

// Constructors

DFactor::DFactor (double p ) : var_labels(), vs(), weights(), factors() {}

/// Constructs factor depending on the variable \a v with uniform distribution
DFactor::DFactor( const Var &v ) {
    DFACTOR_NOT_IMPL;
}

/// Constructs factor depending on variables in \a vars with uniform distribution
DFactor::DFactor( const VarSet& vars ) {
    DFACTOR_NOT_IMPL;
}

/// Constructs factor depending on variables in \a vars with all values set to \a p
/* USED_BY_JTREE */
DFactor::DFactor( const VarSet& vars, double p ) : vs(vars), var_labels(), weights(1), factors(vars.size()) {

    // jtree.cpp line 205 creates empty factors, not sure why
    if (vars.size() > 0) {

        // Iteration order of vars might be undefined, so we get var_labels list first and then use that
        std::vector<int> n_states;
        for (Var v : vars) {
            var_labels.push_back(v.label());
            n_states.push_back(v.states());
        }
        weights << 1.0;
        factors = std::vector<Eigen::MatrixXd>(vars.size());
        factors[0] = Eigen::MatrixXd::Constant(n_states[0], 1, p);
        for (size_t i=1; i<var_labels.size(); ++i) {
            factors[i] = Eigen::MatrixXd::Ones(n_states[i], 1);
        }
        assert(factors.size() == var_labels.size());
    }
}

/// Constructs factor depending on variables in \a vars, copying the values from an array
/** \param vars contains the variables that the new factor should depend on.
    *  \param p Points to array of values to be added.
    */
DFactor::DFactor( const VarSet& vars, const double* p ) {
    DFACTOR_NOT_IMPL;
}

/// Constructs factor depending on variables in \a vars, copying the values from \a p
DFactor::DFactor( const VarSet& vars, const TProb<double> &p ) {
    DFACTOR_NOT_IMPL;
}

/// Constructs factor depending on variables in \a vars, permuting the values given in \a p accordingly
DFactor::DFactor( const std::vector<Var> &vars, const std::vector<double> &p ) {
    DFACTOR_NOT_IMPL;
}

DFactor::DFactor(const std::vector<int> v, const Eigen::VectorXd w, const std::vector<Eigen::MatrixXd> f) :
    vs(), var_labels(v), weights(w), factors(f) {

    // Initialise VarSet vs based on var_labels
    assert(factors.size() == var_labels.size());
    SmallSet<Var> vars_list;
    for (size_t i=0; i<var_labels.size(); ++i) {
        size_t label = v[i];
        size_t n_states = f[i].rows();
        vars_list.insert(Var(label, n_states));
    }
    vs = VarSet(vars_list);
}

    /// \name Get/set individual entries
    //@{
        /// Sets \a i 'th entry to \a val
        /* USED_BY_JTREE */ void DFactor::set( size_t i, double val ) {
            DFACTOR_CANNOT_IMPL;
        }

        /// Gets \a i 'th entry
        double DFactor::get( size_t i ) const {
            DFACTOR_NOT_IMPL;
        }
    //@}

    /// \name Queries
    //@{
        /// Returns constant reference to value vector
        const TProb<double>& DFactor::p() const {
            DFACTOR_NOT_IMPL;
        }

        /// Returns reference to value vector
        /* USED_BY_JTREE */
        TProb<double>& DFactor::p() {
             DFACTOR_CANNOT_IMPL;
        }

        /// Returns a copy of the \a i 'th entry of the value vector
        /* USED_BY_JTREE */
        double DFactor::operator[] (size_t i) const {
            // We can return a particular value, but not sure what multi-dimensional point the "ith entry" refers to
            // Perhaps we need an as_factor method, can then just invoke as normal
            DFACTOR_NOT_IMPL;
        }

        /// Returns constant reference to variable set (i.e., the variables on which the factor depends)
        /* USED_BY_JTREE */
        const VarSet& DFactor::vars() const {
            return vs;
        }

        /// Returns reference to variable set (i.e., the variables on which the factor depends)
        /* USED_BY_JTREE */
        VarSet& DFactor::vars() {
            return vs;
        }

        /// Returns the number of possible joint states of the variables on which the factor depends, \f$\prod_{l\in L} S_l\f$
        /** \note This is equal to the length of the value vector.
         */
        size_t DFactor::nrStates() const {
            DFACTOR_NOT_IMPL;
        }

        /// Returns the Shannon entropy of \c *this, \f$-\sum_i p_i \log p_i\f$
        double DFactor::entropy() const {
            DFACTOR_NOT_IMPL;
        }

        /// Returns maximum of all values
        double DFactor::max() const {
            DFACTOR_NOT_IMPL;
        }

        /// Returns minimum of all values
        double DFactor::min() const {
            DFACTOR_NOT_IMPL;
        }

        /// Returns sum of all values
        double DFactor::sum() const {
            DFACTOR_NOT_IMPL;
        }
        
        /// Returns sum of absolute values
        double DFactor::sumAbs() const {
            DFACTOR_NOT_IMPL;
        }

        /// Returns maximum absolute value of all values
        double DFactor::maxAbs() const {
            DFACTOR_NOT_IMPL;
        }

        /// Returns \c true if one or more values are NaN
        bool DFactor::hasNaNs() const {
            DFACTOR_NOT_IMPL;
        }

        /// Returns \c true if one or more values are negative
        bool DFactor::hasNegatives() const {
            DFACTOR_NOT_IMPL;
        }

        /// Returns strength of this factor (between variables \a i and \a j), as defined in eq. (52) of [\ref MoK07b]
//         double strength( const Var &i, const Var &j ) const;

        /// Comparison
        bool DFactor::operator==( const DFactor& y ) const {
            DFACTOR_NOT_IMPL;
        }

        /// Formats a factor as a string
        std::string DFactor::toString() const {
            DFACTOR_NOT_IMPL;
        }
    //@}

    /// \name Unary transformations
    //@{
        /// Returns negative of \c *this
        DFactor DFactor::operator- () const { 
            DFACTOR_NOT_IMPL;
        }

        /// Returns pointwise absolute value
        DFactor DFactor::abs() const {
            DFACTOR_NOT_IMPL;
        }

        /// Returns pointwise exponent
        DFactor DFactor::exp() const {
            DFACTOR_NOT_IMPL;
        }

        /// Returns pointwise logarithm
        /** If \a zero == \c true, uses <tt>log(0)==0</tt>; otherwise, <tt>log(0)==-Inf</tt>.
         */
        DFactor DFactor::log(bool zero) const {
            DFACTOR_NOT_IMPL;
        }

        /// Returns pointwise inverse
        /** If \a zero == \c true, uses <tt>1/0==0</tt>; otherwise, <tt>1/0==Inf</tt>.
         */
        DFactor DFactor::inverse(bool zero) const {
            DFACTOR_NOT_IMPL;
        }

        /// Returns normalized copy of \c *this, using the specified norm
        /** \throw NOT_NORMALIZABLE if the norm is zero
         */
        DFactor DFactor::normalized( ProbNormType norm) const {
            DFactor res(*this);
            res.normalize(norm);
            return res;
        }
    //@}

    /// \name Unary operations
    //@{
        /// Draws all values i.i.d. from a uniform distribution on [0,1)
        DFactor& DFactor::randomize() {
            DFACTOR_NOT_IMPL;
        }

        /// Sets all values to \f$1/n\f$ where \a n is the number of states
        DFactor& DFactor::setUniform() {
            DFACTOR_NOT_IMPL;
        }

        /// Applies absolute value pointwise
        DFactor& DFactor::takeAbs() {
            DFACTOR_NOT_IMPL;
        }

        /// Applies exponent pointwise
        DFactor& DFactor::takeExp() {
            DFACTOR_NOT_IMPL;
        }

        /// Applies logarithm pointwise
        /** If \a zero == \c true, uses <tt>log(0)==0</tt>; otherwise, <tt>log(0)==-Inf</tt>.
         */
        DFactor& DFactor::takeLog( bool zero) {
            DFACTOR_NOT_IMPL;
        }

        /// Normalizes factor using the specified norm
        /** \throw NOT_NORMALIZABLE if the norm is zero
         */
        /* USED_BY_JTREE */
        double DFactor::normalize( ProbNormType norm) {
            Eigen::VectorXd factor_prod = Eigen::VectorXd::Ones(weights.size());
            if( norm == dai::NORMPROB ) {
                // Calculate Z, the sum of all variable settings; this can be done efficiently due to factorisation
                for (int i=0; i<factors.size(); ++i) {
                    factor_prod.array() *= factors[i].colwise().sum().array();
                }
                double Z = weights.dot(factor_prod);
                if (Z <= 0) {
                    std::cerr << "WARNING: Z<=0, this probably indicates numerical underflow." << std::endl;
                }
                assert (Z >= 0);
                if (Z > DBL_MIN) {
                    weights /= Z;
                }

                return Z;
            }
            else if( norm == dai::NORMLINF ) {
                DFACTOR_NOT_IMPL;
            }
        }
    //@}

    /// \name Operations with scalars
    //@{
        /// Sets all values to \a x
        /* USED_BY_JTREE */
        DFactor& DFactor::fill (double x) {
            weights << 1.0;
            factors[0] = Eigen::MatrixXd::Constant(factors[0].rows(), 1, x);
            for (size_t i=1; i<vs.size(); ++i) {
                factors[i] = Eigen::MatrixXd::Ones(factors[i].rows(), 1);
            }
            return *this;
        }

        /// Adds scalar \a x to each value
        DFactor& DFactor::operator+= (double x) {
            DFACTOR_NOT_IMPL;
        }

        /// Subtracts scalar \a x from each value
        DFactor& DFactor::operator-= (double x) {
            DFACTOR_NOT_IMPL;
        }

        /// Multiplies each value with scalar \a x
        DFactor& DFactor::operator*= (double x) {
            DFACTOR_NOT_IMPL;
        }

        /// Divides each entry by scalar \a x
        DFactor& DFactor::operator/= (double x) {
            DFACTOR_NOT_IMPL;
        }

        /// Raises values to the power \a x
        DFactor& DFactor::operator^= (double x) {
            DFACTOR_NOT_IMPL;
        }
    //@}

    /// \name Transformations with scalars
    //@{
        /// Returns sum of \c *this and scalar \a x
        DFactor DFactor::operator+ (double x) const {
            DFACTOR_NOT_IMPL;
        }

        /// Returns difference of \c *this and scalar \a x
        DFactor DFactor::operator- (double x) const {
            DFACTOR_NOT_IMPL;
        }

        /// Returns product of \c *this with scalar \a x
        DFactor DFactor::operator* (double x) const {
            DFACTOR_NOT_IMPL;
        }

        /// Returns quotient of \c *this with scalar \a x
        DFactor DFactor::operator/ (double x) const {
            DFACTOR_NOT_IMPL;
        }

        /// Returns \c *this raised to the power \a x
        DFactor DFactor::operator^ (double x) const {
            DFACTOR_NOT_IMPL;
        }
    //@}

    /// \name Operations with other factors
    //@{
        /// Applies binary operation \a op on two factors, \c *this and \a g
        /** \tparam binOp Type of function object that accepts two arguments of type \a T and outputs a type \a T
         *  \param g Right operand
         *  \param op Operation of type \a binOp
         */
        /* USED_BY_JTREE */
        template<typename binOp> DFactor& DFactor::binaryOp( const DFactor &g, binOp op ) {
            DFACTOR_NOT_IMPL;
        }

        /// Adds \a g to \c *this
        /** The sum of two factors is defined as follows: if
         *  \f$f : \prod_{l\in L} X_l \to [0,\infty)\f$ and \f$g : \prod_{m\in M} X_m \to [0,\infty)\f$, then
         *  \f[f+g : \prod_{l\in L\cup M} X_l \to [0,\infty) : x \mapsto f(x_L) + g(x_M).\f]
         */
        DFactor& DFactor::operator+= (const DFactor& g) {
            DFACTOR_NOT_IMPL;
        }

        /// Subtracts \a g from \c *this
        /** The difference of two factors is defined as follows: if
         *  \f$f : \prod_{l\in L} X_l \to [0,\infty)\f$ and \f$g : \prod_{m\in M} X_m \to [0,\infty)\f$, then
         *  \f[f-g : \prod_{l\in L\cup M} X_l \to [0,\infty) : x \mapsto f(x_L) - g(x_M).\f]
         */
        DFactor& DFactor::operator-= (const DFactor& g) {
            DFACTOR_NOT_IMPL;
        }

        /// Multiplies \c *this with \a g
        /** The product of two factors is defined as follows: if
         *  \f$f : \prod_{l\in L} X_l \to [0,\infty)\f$ and \f$g : \prod_{m\in M} X_m \to [0,\infty)\f$, then
         *  \f[fg : \prod_{l\in L\cup M} X_l \to [0,\infty) : x \mapsto f(x_L) g(x_M).\f]
         */
        /* USED_BY_JTREE */
        DFactor& DFactor::operator*= (const DFactor& g) {

            // MAIN TBP METHOD: Multiply two factors using sampling.

            const int N_ARGS = 2;
            double n_terms_exact = weights.size() * g.weights.size();

            // Generate sample in the form of (indices -> freq) pairs, where indices is a list containing a term id from
            // each input DFactor, and freq is the count of that term combination. In the case of exact multiplication
            // (tbp_sample_size == 0), freq is instead the product of each combination of weights

            std::map<std::vector<int>, double> sample;
            if (n_terms_exact > EXACT_THRESHOLD) {

                // Initialise RNG
                std::random_device rd;
                std::mt19937 rng(rd());
                std::vector<std::discrete_distribution<int> > dists {
                    std::discrete_distribution<int>(weights.data(), weights.data() + weights.size()),
                    std::discrete_distribution<int>(g.weights.data(), g.weights.data() + g.weights.size())
                };

                // Generate sample
                int debug_total = 0;
                for (size_t i=0; i<tbp_sample_size; ++i) {
                    std::vector<int> ks(N_ARGS);
                    for (size_t j=0; j<N_ARGS; ++j) {
                        ks[j] = dists[j](rng);
                    }
                    sample[ks] += 1;
                }

            } else {
                for (size_t i=0; i<weights.size(); ++i) {
                    for (size_t j=0; j<g.weights.size(); ++j) {
                        sample[std::vector<int>{(int)i,(int)j}] = weights[i] * g.weights[j];
                    }
                }
            }

            // Initialise result vars, weights and factors
            std::map<int, int> var_sizes;
            for (DFactor df : {*this, g}) {
                for (size_t i=0; i<df.var_labels.size(); ++i) {
                    var_sizes[df.var_labels[i]] = df.factors[i].rows();
                }
            }

            Eigen::VectorXd res_weights(sample.size());
            std::vector<Eigen::MatrixXd> res_factors(var_sizes.size());
            std::vector<int> res_var_labels(var_sizes.size());
            std::map<int, int> var_indices;
            int i = 0;
            for (auto const &kv : var_sizes) {
                res_var_labels[i] = kv.first;
                var_indices[kv.first] = i;
                res_factors[i] = Eigen::MatrixXd::Ones(kv.second, sample.size());
                i++;
            }

            // Iterate over (ks -> freq) pairs in sample; each iteration corresponds to a single term in the result
            i = 0;
            for (auto const &term : sample) {
                // Sampled terms from each input DFactor: [k_1, k_2, ..., k_{N_ARGS}]
                std::vector<int> ks = term.first;
                // Frequency that this combination of sampled terms appears in sample, not an int because if using
                // exact multiplication we get floating point weights
                double freq = term.second;
                // Construct resulting term by multiplying sampled columns onto corresponding result columns according
                // to var ids
                for (size_t j=0; j<N_ARGS; ++j) {
                    for (size_t l=0; l<(j == 0 ? *this : g).var_labels.size(); ++l) {
                        // TODO: Is this the most efficient way to do in-place elementwise multiplication?
                        // http://stackoverflow.com/questions/25328694/eigen3-coefficient-wise-multiplication-in-place
                        res_factors[var_indices[(j == 0 ? *this : g).var_labels[l]]].col(i).array() *=
                            (j == 0 ? *this : g).factors[l].col(ks[j]).array();
                    }
                }
                res_weights[i] = freq;
                i++;
            }

            factors = res_factors;
            weights = res_weights;
            var_labels = res_var_labels;
            vs |= g.vs;

            reweight();

            return *this;
        }

        /// Divides \c *this by \a g (where division by zero yields zero)
        /** The quotient of two factors is defined as follows: if
         *  \f$f : \prod_{l\in L} X_l \to [0,\infty)\f$ and \f$g : \prod_{m\in M} X_m \to [0,\infty)\f$, then
         *  \f[\frac{f}{g} : \prod_{l\in L\cup M} X_l \to [0,\infty) : x \mapsto \frac{f(x_L)}{g(x_M)}.\f]
         */
        DFactor& DFactor::operator/= (const DFactor& g) {
            DFACTOR_NOT_IMPL;
        }
    //@}

    /// \name Transformations with other factors
    //@{
        /// Returns result of applying binary operation \a op on two factors, \c *this and \a g
        /** \tparam binOp Type of function object that accepts two arguments of type \a T and outputs a type \a T
         *  \param g Right operand
         *  \param op Operation of type \a binOp
         */
        /* USED_BY_JTREE */
        template<typename binOp> DFactor DFactor::binaryTr( const DFactor &g, binOp op ) const {
            DFACTOR_NOT_IMPL;
        }

        /// Returns sum of \c *this and \a g
        /** The sum of two factors is defined as follows: if
         *  \f$f : \prod_{l\in L} X_l \to [0,\infty)\f$ and \f$g : \prod_{m\in M} X_m \to [0,\infty)\f$, then
         *  \f[f+g : \prod_{l\in L\cup M} X_l \to [0,\infty) : x \mapsto f(x_L) + g(x_M).\f]
         */
        DFactor DFactor::operator+ (const DFactor& g) const {
            DFACTOR_NOT_IMPL;
        }

        /// Returns \c *this minus \a g
        /** The difference of two factors is defined as follows: if
         *  \f$f : \prod_{l\in L} X_l \to [0,\infty)\f$ and \f$g : \prod_{m\in M} X_m \to [0,\infty)\f$, then
         *  \f[f-g : \prod_{l\in L\cup M} X_l \to [0,\infty) : x \mapsto f(x_L) - g(x_M).\f]
         */
        DFactor DFactor::operator- (const DFactor& g) const {
            DFACTOR_NOT_IMPL;
        }

        /// Returns product of \c *this with \a g
        /** The product of two factors is defined as follows: if
         *  \f$f : \prod_{l\in L} X_l \to [0,\infty)\f$ and \f$g : \prod_{m\in M} X_m \to [0,\infty)\f$, then
         *  \f[fg : \prod_{l\in L\cup M} X_l \to [0,\infty) : x \mapsto f(x_L) g(x_M).\f]
         */
        DFactor DFactor::operator* (const DFactor& g) const {
            DFACTOR_NOT_IMPL;
        }

        /// Returns quotient of \c *this by \a f (where division by zero yields zero)
        /** The quotient of two factors is defined as follows: if
         *  \f$f : \prod_{l\in L} X_l \to [0,\infty)\f$ and \f$g : \prod_{m\in M} X_m \to [0,\infty)\f$, then
         *  \f[\frac{f}{g} : \prod_{l\in L\cup M} X_l \to [0,\infty) : x \mapsto \frac{f(x_L)}{g(x_M)}.\f]
         */
        /* USED_BY_JTREE */
        DFactor DFactor::operator/ (const DFactor& g) const {
            DFACTOR_CANNOT_IMPL;
        }
    //@}

    /// \name Miscellaneous operations
    //@{
        /// Returns a slice of \c *this, where the subset \a vars is in state \a varsState
        /** \pre \a vars sould be a subset of vars()
         *  \pre \a varsState < vars.nrStates()
         *
         *  The result is a factor that depends on the variables of *this except those in \a vars,
         *  obtained by setting the variables in \a vars to the joint state specified by the linear index
         *  \a varsState. Formally, if \c *this corresponds with the factor \f$f : \prod_{l\in L} X_l \to [0,\infty)\f$,
         *  \f$M \subset L\f$ corresponds with \a vars and \a varsState corresponds with a mapping \f$s\f$ that
         *  maps a variable \f$x_m\f$ with \f$m\in M\f$ to its state \f$s(x_m) \in X_m\f$, then the slice
         *  returned corresponds with the factor \f$g : \prod_{l \in L \setminus M} X_l \to [0,\infty)\f$
         *  defined by \f$g(\{x_l\}_{l\in L \setminus M}) = f(\{x_l\}_{l\in L \setminus M}, \{s(x_m)\}_{m\in M})\f$.
         */
//         DFactor slice( const VarSet& vars, size_t varsState ) const; 

        /// Embeds this factor in a larger VarSet
        /** \pre vars() should be a subset of \a vars 
         *
         *  If *this corresponds with \f$f : \prod_{l\in L} X_l \to [0,\infty)\f$ and \f$L \subset M\f$, then
         *  the embedded factor corresponds with \f$g : \prod_{m\in M} X_m \to [0,\infty) : x \mapsto f(x_L)\f$.
         */
        DFactor DFactor::embed(const VarSet & vars) const {
            DFACTOR_NOT_IMPL;
        }

        /// Returns marginal on \a vars, obtained by summing out all variables except those in \a vars, and normalizing the result if \a normed == \c true
//         DFactor marginal(const VarSet &vars, bool normed=true) const;

        /// Returns max-marginal on \a vars, obtained by maximizing all variables except those in \a vars, and normalizing the result if \a normed == \c true
//         DFactor maxMarginal(const VarSet &vars, bool normed=true) const;

        TFactor<Real> DFactor::asTFactor() const {
            TFactor<Real> res(vs, 0.0);
            for (size_t k=0; k<weights.size(); ++k) {
                for (size_t i=0; i<vs.nrStates(); ++i) {
                    std::map<Var, size_t> vals = calcState(vs, i);
                    double joint_p = 1;
                    for (size_t j=0; j<var_labels.size(); ++j) {
                        Var v;
                        for (size_t l=0; l<vs.size(); ++l) {
                            if (vs.elements()[l].label() == var_labels[j]) {
                                v = vs.elements()[l];
                                break;
                            }
                        }
                        joint_p *= factors[j](vals[v], k);
                    }
                    res.p().set(i, res.p().get(i) + weights[k] * joint_p);
                }
            }
            return res;
        }

    //@}

void DFactor::reweight() {

    // Note: This probably doesn't work if factors are negative

    // Possible reweightings:
    // "normonly": Normalise only
    // "maxnorm": Every term (rank one tensor) is divided by its max and corresponding weight multiplied by this, so
    //            that each term has max-norm 1. Then normalise.
    // "minvar": Reweighting to minimise sum-of-component variance. Then normalise.
    // Normalise here means make the total probability mass 1 (not the weight sum). We do this for numerical reasons,
    // and because the original factor.h code did it so it might be relied upon by junction tree (though probably not).
    // We never normalise the weight sum to be 1, since sampler in multiply method doesn't require normalised weights.

    // First, fix any negative weights by moving negative into the first factor
    for (int i=0; i<weights.size(); ++i) {
        if (weights[i] < 0) {
            weights[i] *= -1;
            factors[0].col(i).array() *= -1;
        }
    }

    if (REWEIGHT_METHOD == REWEIGHT_MAXNORM) {

        for (int i=0; i<factors.size(); ++i) {
//             factors[i].rowwise() /= factors[i].colwise().maxCoeff();
            // Divide each vector (column) by its maximum - broadcasting seems to works for += but not /=, so we use a loop instead
            for (int j=0; j<factors[i].cols(); ++j) {
                double max_coeff = factors[i].col(j).maxCoeff();
                if (max_coeff > DBL_MIN) {
                    factors[i].col(j) /= max_coeff;
                    weights[j] *= max_coeff;
                } else {
                    // It's possible to have nonzero doubles less than DBL_MIN (unnormalised, see IEEE floating point
                    // spec). These cause problems later, so we set to zero.
                    factors[i].col(j).setZero();
                    weights[j] = 0;
                }
            }
        }


    } else if (REWEIGHT_METHOD == REWEIGHT_MINVAR) {
        // Scale terms and weights in-place so the function is the same but the sum of cell variances is minimised when
        // sampling. We do not change terms other than scaling them.
        Eigen::VectorXd w_new = Eigen::VectorXd::Ones(weights.size());
        for (int i=0; i<factors.size(); ++i) {
            w_new.array() *= factors[i].colwise().squaredNorm().array();
        }
        w_new = weights.array() * w_new.array().sqrt();
        w_new.array() /= w_new.sum();

        for (int k=0; k<weights.size(); ++k) {
            // TODO choice of which factor absorbs the scaling term (i.e. 0 here) is currently arbitrary - either prove
            // it doesn't matter or change data structure to include an explicit scaling term
            factors[0].col(k).array() *= weights[k] / w_new[k];
            weights[k] = w_new[k];
        }
    }

    normalize(dai::NORMPROB);

}

DFactor DFactor::slice( const VarSet& vars, size_t varsState ) const {
    DFACTOR_NOT_IMPL;
}


/* USED_BY_JTREE */
DFactor DFactor::marginal(const VarSet &vars, bool normed) const {

    DFactor res(*this);
    res.vs = vs & vars;
    // / is set minus, operator/ is defined in smallset.h
    VarSet to_remove = vs / vars;
    // Marginalise out and remove each variable
    for (Var v : to_remove) {
        int var_index = std::find(res.var_labels.begin(), res.var_labels.end(), v.label()) - res.var_labels.begin();
        Eigen::VectorXd colsums = res.factors[var_index].colwise().sum();
        res.weights.array() *= colsums.array();

        res.var_labels.erase(res.var_labels.begin() + var_index);
        res.factors.erase(res.factors.begin() + var_index);
    }

    // If there is only one variable in resulting factor, we can simplify further to only use one term
    if (res.factors.size() == 1) {
        res.factors[0] = res.factors[0] * res.weights;
        // A vector with one element equal to 1
        Eigen::VectorXd one(1,1);
        one << 1.0;
        res.weights = one;
    }

    res.reweight();

    return res;
}


DFactor DFactor::maxMarginal(const VarSet &vars, bool normed) const {
    DFACTOR_NOT_IMPL;
}


double DFactor::strength( const Var &i, const Var &j ) const {
    DFACTOR_NOT_IMPL;
}


/// Writes a factor to an output stream
/** \relates DFactor
 */
std::ostream& operator<< (std::ostream& os, const DFactor& f) {
    std::string var_str;
    for (int var : f.var_labels) {
        var_str += std::to_string(var) + " ";
    }
    std::string weights_str;
    for (int i=0; i<f.weights.size(); ++i) {
        weights_str += std::to_string(f.weights[i]) + " ";
    }
    os << "[DFactor]\nVars (" << f.var_labels.size() << "): " << var_str << "\nWeights (" << f.weights.size() << "): "
       << weights_str << "\nFactors (" << f.factors.size() << "):\n";
    for (Eigen::MatrixXd mat : f.factors) {
        os << mat << endl << endl;
    }
    return os;
}

/// Returns distance between two factors \a f and \a g, according to the distance measure \a dt
/** \relates DFactor
 *  \pre f.vars() == g.vars()
 */
double dist( const DFactor &f, const DFactor &g, ProbDistType dt ) {
    DFACTOR_NOT_IMPL;
}


/// Returns the pointwise maximum of two factors
/** \relates DFactor
 *  \pre f.vars() == g.vars()
 */
DFactor max( const DFactor &f, const DFactor &g ) {
    DFACTOR_NOT_IMPL;
}


/// Returns the pointwise minimum of two factors
/** \relates DFactor
 *  \pre f.vars() == g.vars()
 */
DFactor min( const DFactor &f, const DFactor &g ) {
    DFACTOR_NOT_IMPL;
}


/// Calculates the mutual information between the two variables that \a f depends on, under the distribution given by \a f
/** \relates DFactor
 *  \pre f.vars().size() == 2
 */
double MutualInfo(const DFactor &f) {
    DFACTOR_NOT_IMPL;
}

} // end of namespace dai


