#include <dai/tbp.h>
#include <assert.h>

namespace dai {

    using namespace std;

    int tbp_sample_size;

    vector<Eigen::VectorXd> tbp_marg(vector<DFactor> factors, int k) {

        tbp_sample_size = k;

        for (size_t i=0; i<factors.size(); ++i) {
            factors[i].reweight();
        }

        FactorGraph fg(factors);

        PropertySet opts;
        opts.set("verbose",(size_t)0);       // Verbosity (amount of output generated)

        string heuristic = "MAXCOMMONVARS";
//        string heuristic = "TBPMINWEIGHT";
//        string heuristic = "MINFILL";
//        string heuristic = "RANDOM";
        cerr << "Using elimintation heuristic: " << heuristic << endl;
        JTree jt(fg, opts("updates",string("SHSH"))("heuristic", heuristic));

        jt.init();
        jt.run();

        vector<Eigen::VectorXd> res(fg.nrVars());
        for (size_t i=0; i<fg.nrVars(); i++) {
            // If this assertion fails, it probably means the input graph had non-consecutively numbered variables
            // and we need to return a map of (label -> marginals) rather than a simple vector
            assert (fg.var(i).label() == i);
            TProb<Real> marg = jt.belief(fg.var(i)).asTFactor().p();
            res[i] = Eigen::VectorXd(marg.size());
            for (size_t j=0; j<marg.size(); ++j) {
                res[i][j] = marg[j];
            }
        }

        return res;
    }

    /// Reads a decomposed graph from an input stream
    std::vector<dai::DFactor> read_dfg_from_stream(std::istream& is) {

        // File format:

        // n_factors
        //
        // for each factor:
        //     n_terms
        //     <weights>
        //     n_variables
        //     <variable indices>
        //     <variable cardinalities>
        //     for each variable:
        //         n_nonzero
        //         1 0.5
        //         3 0.1
        //         4 0.1
        //         ...
        //         n_nonzero
        //         1 0.5
        //         3 0.1
        //         4 0.3
        //         ...

        // Increase this to debug
        int verbose = 0;

        std::vector<dai::DFactor> dfg;

        string line;

        while( (is.peek()) == '#' )
            getline(is,line);

        size_t nr_factors;
        is >> nr_factors;
        if( is.fail() )
            DAI_THROWE(INVALID_FACTORGRAPH_FILE,"Cannot read number of factors");

        getline(is,line);
        if( is.fail() || line.size() > 0 )
            DAI_THROWE(INVALID_FACTORGRAPH_FILE,"Expecting empty line");
        if( verbose >= 1 )
            cerr << "Reading " << nr_factors << " factors..." << endl;

        map<int,size_t> vardims;
        for( size_t I = 0; I < nr_factors; I++ ) {
            if (verbose >= 2)
                cerr << "Reading factor " << I << "..." << endl;

            // nr_terms
            size_t nr_terms;
            while ((is.peek()) == '#')
                getline(is, line);
            is >> nr_terms;
            if (verbose >= 2)
                cerr << "  nr_terms: " << nr_terms << endl;

            // weights
            vector<double> w;
            for (size_t mi = 0; mi < nr_terms; mi++) {
                double mi_weight;
                while ((is.peek()) == '#')
                    getline(is, line);
                is >> mi_weight;
                w.push_back(mi_weight);
            }
            Eigen::VectorXd weights = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(w.data(), w.size());
            if (verbose >= 2)
                cerr << "  weights: " << weights << endl;

            // nr_vars
            size_t nr_vars;
            while ((is.peek()) == '#')
                getline(is, line);
            is >> nr_vars;
            if (verbose >= 2)
                cerr << "  nr_vars: " << nr_vars << endl;

            // var labels
            vector<int> labels;
            for (size_t mi = 0; mi < nr_vars; mi++) {
                int mi_label;
                while ((is.peek()) == '#')
                    getline(is, line);
                is >> mi_label;
                labels.push_back(mi_label);
            }
            if (verbose >= 2)
                cerr << "  labels: " << labels << endl;

            // var cardinalities
            vector <size_t> dims;
            for (size_t mi = 0; mi < nr_vars; mi++) {
                size_t mi_dim;
                while ((is.peek()) == '#')
                    getline(is, line);
                is >> mi_dim;
                dims.push_back(mi_dim);
            }
            if (verbose >= 2)
                cerr << "  dimensions: " << dims << endl;

            // Factor matrices (one matrix per variable)

            std::vector <Eigen::MatrixXd> matrices;

            for (size_t mi = 0; mi < nr_vars; mi++) {

                if (verbose >= 2)
                    cerr << "  Reading matrix " << mi << endl;

                // Check cardinality of this variable is consistent with past occurrences
                map<int, size_t>::iterator vdi = vardims.find(labels[mi]);
                if (vdi != vardims.end()) {
                    if (vdi->second != dims[mi])
                        DAI_THROWE(INVALID_FACTORGRAPH_FILE,
                                   "Variable with label " + boost::lexical_cast<string>(labels[mi]) +
                                   " has inconsistent dimensions.");
                } else
                    vardims[labels[mi]] = dims[mi];

                // Read values into (<variable cardinality> * <n_terms>) matrix

                Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(dims[mi], nr_terms);
                size_t nr_nonzeros;
                while ((is.peek()) == '#')
                    getline(is, line);
                is >> nr_nonzeros;
                if (verbose >= 2)
                    cerr << "     nonzeroes: " << nr_nonzeros << endl;

                for (size_t k = 0; k < nr_nonzeros; k++) {
                    size_t li;
                    double val;
                    while ((is.peek()) == '#')
                        getline(is, line);
                    if (!(is >> li)) {
                        DAI_THROWE(INVALID_FACTORGRAPH_FILE,"Matrix index could not be parsed as int");
                    }
                    while ((is.peek()) == '#')
                        getline(is, line);
                    if (!(is >> val)) {
                        std::stringstream s;
                        s << "Matrix value could not be cast as double - check that matrix values in .dfg file are within "
                        "C++ double limits [" << DBL_MIN << ", " << DBL_MAX << "]";
                        DAI_THROWE(INVALID_FACTORGRAPH_FILE, s.str());
                    }

                    // The matrix we are trying to recover has dimensions (<variable cardinality> * <n_terms>)
                    // indexed in column-major order - i.e. we first read all the values for the first term, then the
                    // second term, etc. dims[mi] is the cardinality.

                    size_t col = li / dims[mi];
                    size_t row = li % dims[mi];
                    if (verbose >= 2)
                        cerr << "     index: " << li << " (" << row << ", " << col << ") value: " << val << endl;
                    mat(row, col) = val;


                }

                matrices.push_back(mat);

            }

            dfg.push_back(DFactor(labels, weights, matrices));
        }

        return dfg;


    }


}






