
/// \file
/// \brief Defines DFactor class which represents decomposed factors in probability distributions.


#ifndef __defined_libdai_dfactor_h
#define __defined_libdai_dfactor_h

#include <iostream>
#include <functional>
#include <cmath>
#include <dai/prob.h>
#include <dai/varset.h>
#include <dai/index.h>
#include <dai/util.h>

#include <initializer_list>
#include <Eigen/Dense>
#include <random>
#include <dai/dfactor_macros.h>
// Also includes realfactor.h, because DFactor has asTFactor() function
#include <dai/realfactor.h>

#define EXACT_THRESHOLD 1

// Possible reweightings:
// "NORMONLY": Normalise only
// "MAXNORM": Every term (rank one tensor) is divided by its max and corresponding weight multiplied by this, so that each term has max-norm 1. Then normalise.
// "MINVAR": Reweighting to minimise sum-of-component variance. Then normalise.
// Normalise here means make the total probability mass 1 (not the weight sum). We do this for numerical reasons.
// We never normalise the weight sum to be 1, since sampler in multiply method doesn't require normalised weights.

#define REWEIGHT_NORMONLY 0
#define REWEIGHT_MAXNORM 1
#define REWEIGHT_MINVAR 2

#define REWEIGHT_METHOD REWEIGHT_MAXNORM

namespace dai {

extern int tbp_sample_size;

class DFactor {
    public:
        std::vector<int> var_labels;
        VarSet vs;
        Eigen::VectorXd weights;
        std::vector<Eigen::MatrixXd> factors;

    public:
    /// \name Constructors and destructors
    //@{
        /// Constructs factor depending on no variables with value \a p
        /* USED_BY_JTREE */
        // TODO No idea what this is for; probably need to add a 1 to weights or factors
        DFactor (double p = 1 );

        /// Constructs factor depending on the variable \a v with uniform distribution
        DFactor( const Var &v );

        /// Constructs factor depending on variables in \a vars with uniform distribution
        DFactor( const VarSet& vars );

        /// Constructs factor depending on variables in \a vars with all values set to \a p
        /* USED_BY_JTREE */
        DFactor( const VarSet& vars, double p );

        /// Constructs factor depending on variables in \a vars, copying the values from a std::vector<>
        /** \tparam S Type of values of \a x
         *  \param vars contains the variables that the new factor should depend on.
         *  \param x Vector with values to be copied.
         */
        template<typename S>
        DFactor( const VarSet& vars, const std::vector<S> &x ) {
            DFACTOR_NOT_IMPL;
        }

        /// Constructs factor depending on variables in \a vars, copying the values from an array
        /** \param vars contains the variables that the new factor should depend on.
         *  \param p Points to array of values to be added.
         */
        DFactor( const VarSet& vars, const double* p );

        /// Constructs factor depending on variables in \a vars, copying the values from \a p
        DFactor( const VarSet& vars, const TProb<double> &p );

        /// Constructs factor depending on variables in \a vars, permuting the values given in \a p accordingly
        DFactor( const std::vector<Var> &vars, const std::vector<double> &p );

        DFactor(const std::vector<int> v, const Eigen::VectorXd w, const std::vector<Eigen::MatrixXd> f);

    //@}

    /// \name Get/set individual entries
    //@{
        /// Sets \a i 'th entry to \a val
        /* USED_BY_JTREE */
        void set( size_t i, double val );

        /// Gets \a i 'th entry
        double get( size_t i ) const;
    //@}

    /// \name Queries
    //@{
        /// Returns constant reference to value vector
        const TProb<double>& p() const;

        /// Returns reference to value vector
        /* USED_BY_JTREE */
        TProb<double>& p();

        /// Returns a copy of the \a i 'th entry of the value vector
        /* USED_BY_JTREE */
        double operator[] (size_t i) const;

        /// Returns constant reference to variable set (i.e., the variables on which the factor depends)
        /* USED_BY_JTREE */
        const VarSet& vars() const;

        /// Returns reference to variable set (i.e., the variables on which the factor depends)
        /* USED_BY_JTREE */
        VarSet& vars();

        /// Returns the number of possible joint states of the variables on which the factor depends, \f$\prod_{l\in L} S_l\f$
        /** \note This is equal to the length of the value vector.
         */
        size_t nrStates() const;

        /// Returns the Shannon entropy of \c *this, \f$-\sum_i p_i \log p_i\f$
        double entropy() const;

        /// Returns maximum of all values
        double max() const;

        /// Returns minimum of all values
        double min() const;

        /// Returns sum of all values
        double sum() const;
        
        /// Returns sum of absolute values
        double sumAbs() const;

        /// Returns maximum absolute value of all values
        double maxAbs() const;

        /// Returns \c true if one or more values are NaN
        bool hasNaNs() const;

        /// Returns \c true if one or more values are negative
        bool hasNegatives() const;

        /// Returns strength of this factor (between variables \a i and \a j), as defined in eq. (52) of [\ref MoK07b]
        double strength( const Var &i, const Var &j ) const;

        /// Comparison
        bool operator==( const DFactor& y ) const;

        /// Formats a factor as a string
        std::string toString() const;
    //@}

    /// \name Unary transformations
    //@{
        /// Returns negative of \c *this
        DFactor operator- () const;

        /// Returns pointwise absolute value
        DFactor abs() const;

        /// Returns pointwise exponent
        DFactor exp() const;

        /// Returns pointwise logarithm
        /** If \a zero == \c true, uses <tt>log(0)==0</tt>; otherwise, <tt>log(0)==-Inf</tt>.
         */
        DFactor log(bool zero=false) const;

        /// Returns pointwise inverse
        /** If \a zero == \c true, uses <tt>1/0==0</tt>; otherwise, <tt>1/0==Inf</tt>.
         */
        DFactor inverse(bool zero=true) const;

        /// Returns normalized copy of \c *this, using the specified norm
        /** \throw NOT_NORMALIZABLE if the norm is zero
         */
        DFactor normalized( ProbNormType norm=NORMPROB ) const;
    //@}

    /// \name Unary operations
    //@{
        /// Draws all values i.i.d. from a uniform distribution on [0,1)
        DFactor& randomize();

        /// Sets all values to \f$1/n\f$ where \a n is the number of states
        DFactor& setUniform();

        /// Applies absolute value pointwise
        DFactor& takeAbs();

        /// Applies exponent pointwise
        DFactor& takeExp();

        /// Applies logarithm pointwise
        /** If \a zero == \c true, uses <tt>log(0)==0</tt>; otherwise, <tt>log(0)==-Inf</tt>.
         */
        DFactor& takeLog( bool zero = false );

        /// Normalizes factor using the specified norm
        /** \throw NOT_NORMALIZABLE if the norm is zero
         */
        /* USED_BY_JTREE */
        double normalize( ProbNormType norm=NORMPROB );
    //@}

    /// \name Operations with scalars
    //@{
        /// Sets all values to \a x
        /* USED_BY_JTREE */
        DFactor& fill (double x);

        /// Adds scalar \a x to each value
        DFactor& operator+= (double x);

        /// Subtracts scalar \a x from each value
        DFactor& operator-= (double x);

        /// Multiplies each value with scalar \a x
        DFactor& operator*= (double x);

        /// Divides each entry by scalar \a x
        DFactor& operator/= (double x);

        /// Raises values to the power \a x
        DFactor& operator^= (double x);
    //@}

    /// \name Transformations with scalars
    //@{
        /// Returns sum of \c *this and scalar \a x
        DFactor operator+ (double x) const;

        /// Returns difference of \c *this and scalar \a x
        DFactor operator- (double x) const;

        /// Returns product of \c *this with scalar \a x
        DFactor operator* (double x) const;

        /// Returns quotient of \c *this with scalar \a x
        DFactor operator/ (double x) const;

        /// Returns \c *this raised to the power \a x
        DFactor operator^ (double x) const;
    //@}

    /// \name Operations with other factors
    //@{
        /// Applies binary operation \a op on two factors, \c *this and \a g
        /** \tparam binOp Type of function object that accepts two arguments of type \a T and outputs a type \a T
         *  \param g Right operand
         *  \param op Operation of type \a binOp
         */
        /* USED_BY_JTREE */
        template<typename binOp> DFactor& binaryOp( const DFactor &g, binOp op );

        /// Adds \a g to \c *this
        /** The sum of two factors is defined as follows: if
         *  \f$f : \prod_{l\in L} X_l \to [0,\infty)\f$ and \f$g : \prod_{m\in M} X_m \to [0,\infty)\f$, then
         *  \f[f+g : \prod_{l\in L\cup M} X_l \to [0,\infty) : x \mapsto f(x_L) + g(x_M).\f]
         */
        DFactor& operator+= (const DFactor& g);

        /// Subtracts \a g from \c *this
        /** The difference of two factors is defined as follows: if
         *  \f$f : \prod_{l\in L} X_l \to [0,\infty)\f$ and \f$g : \prod_{m\in M} X_m \to [0,\infty)\f$, then
         *  \f[f-g : \prod_{l\in L\cup M} X_l \to [0,\infty) : x \mapsto f(x_L) - g(x_M).\f]
         */
        DFactor& operator-= (const DFactor& g);

        /// Multiplies \c *this with \a g
        /** The product of two factors is defined as follows: if
         *  \f$f : \prod_{l\in L} X_l \to [0,\infty)\f$ and \f$g : \prod_{m\in M} X_m \to [0,\infty)\f$, then
         *  \f[fg : \prod_{l\in L\cup M} X_l \to [0,\infty) : x \mapsto f(x_L) g(x_M).\f]
         */
        /* USED_BY_JTREE */
        DFactor& operator*= (const DFactor& g);

        /// Divides \c *this by \a g (where division by zero yields zero)
        /** The quotient of two factors is defined as follows: if
         *  \f$f : \prod_{l\in L} X_l \to [0,\infty)\f$ and \f$g : \prod_{m\in M} X_m \to [0,\infty)\f$, then
         *  \f[\frac{f}{g} : \prod_{l\in L\cup M} X_l \to [0,\infty) : x \mapsto \frac{f(x_L)}{g(x_M)}.\f]
         */
        DFactor& operator/= (const DFactor& g);
    //@}

    /// \name Transformations with other factors
    //@{
        /// Returns result of applying binary operation \a op on two factors, \c *this and \a g
        /** \tparam binOp Type of function object that accepts two arguments of type \a T and outputs a type \a T
         *  \param g Right operand
         *  \param op Operation of type \a binOp
         */
        /* USED_BY_JTREE */
        template<typename binOp> DFactor binaryTr( const DFactor &g, binOp op ) const;

        /// Returns sum of \c *this and \a g
        /** The sum of two factors is defined as follows: if
         *  \f$f : \prod_{l\in L} X_l \to [0,\infty)\f$ and \f$g : \prod_{m\in M} X_m \to [0,\infty)\f$, then
         *  \f[f+g : \prod_{l\in L\cup M} X_l \to [0,\infty) : x \mapsto f(x_L) + g(x_M).\f]
         */
        DFactor operator+ (const DFactor& g) const;

        /// Returns \c *this minus \a g
        /** The difference of two factors is defined as follows: if
         *  \f$f : \prod_{l\in L} X_l \to [0,\infty)\f$ and \f$g : \prod_{m\in M} X_m \to [0,\infty)\f$, then
         *  \f[f-g : \prod_{l\in L\cup M} X_l \to [0,\infty) : x \mapsto f(x_L) - g(x_M).\f]
         */
        DFactor operator- (const DFactor& g) const;

        /// Returns product of \c *this with \a g
        /** The product of two factors is defined as follows: if
         *  \f$f : \prod_{l\in L} X_l \to [0,\infty)\f$ and \f$g : \prod_{m\in M} X_m \to [0,\infty)\f$, then
         *  \f[fg : \prod_{l\in L\cup M} X_l \to [0,\infty) : x \mapsto f(x_L) g(x_M).\f]
         */
        DFactor operator* (const DFactor& g) const;

        /// Returns quotient of \c *this by \a f (where division by zero yields zero)
        /** The quotient of two factors is defined as follows: if
         *  \f$f : \prod_{l\in L} X_l \to [0,\infty)\f$ and \f$g : \prod_{m\in M} X_m \to [0,\infty)\f$, then
         *  \f[\frac{f}{g} : \prod_{l\in L\cup M} X_l \to [0,\infty) : x \mapsto \frac{f(x_L)}{g(x_M)}.\f]
         */
        /* USED_BY_JTREE */
        DFactor operator/ (const DFactor& g) const;
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
        DFactor slice( const VarSet& vars, size_t varsState ) const;

        /// Embeds this factor in a larger VarSet
        /** \pre vars() should be a subset of \a vars 
         *
         *  If *this corresponds with \f$f : \prod_{l\in L} X_l \to [0,\infty)\f$ and \f$L \subset M\f$, then
         *  the embedded factor corresponds with \f$g : \prod_{m\in M} X_m \to [0,\infty) : x \mapsto f(x_L)\f$.
         */
        DFactor embed(const VarSet & vars) const;

        /// Returns marginal on \a vars, obtained by summing out all variables except those in \a vars, and normalizing the result if \a normed == \c true
        DFactor marginal(const VarSet &vars, bool normed=true) const;

        /// Returns max-marginal on \a vars, obtained by maximizing all variables except those in \a vars, and normalizing the result if \a normed == \c true
        DFactor maxMarginal(const VarSet &vars, bool normed=true) const;

        void reweight();

        TFactor<Real> asTFactor() const;

        double weightSum() const {
            return weights.sum();
        };

    //@}
};


/// Writes a factor to an output stream
/** \relates DFactor
 */
std::ostream& operator<< (std::ostream& os, const DFactor& f);

/// Returns distance between two factors \a f and \a g, according to the distance measure \a dt
/** \relates DFactor
 *  \pre f.vars() == g.vars()
 */
double dist( const DFactor &f, const DFactor &g, ProbDistType dt );

/// Returns the pointwise maximum of two factors
/** \relates DFactor
 *  \pre f.vars() == g.vars()
 */
DFactor max( const DFactor &f, const DFactor &g );

/// Returns the pointwise minimum of two factors
/** \relates DFactor
 *  \pre f.vars() == g.vars()
 */
DFactor min( const DFactor &f, const DFactor &g );

/// Calculates the mutual information between the two variables that \a f depends on, under the distribution given by \a f
/** \relates DFactor
 *  \pre f.vars().size() == 2
 */
double MutualInfo(const DFactor &f);



} // end of namespace dai


#endif
