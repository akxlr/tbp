% [logZ,q,qv,qf,qmap,margs] = dai_jtree(psi,varsets,opts)
%
% INPUT:  psi        = linear cell array containing the factors
%                      (psi{i} should be a structure with a Member field
%                      and a P field).
%         varsets    = linear cell array containing varsets for which marginals
%                      are requested.
%         opts       = string of options.
% 
% OUTPUT: logZ       = logarithm of the partition sum.
%         q          = linear cell array containing all calculated marginals.
%         qv         = linear cell array containing all variable marginals.
%         qf         = linear cell array containing all factor marginals.
%         qmap       = linear array containing the MAP state.
%         margs      = linear cell array containing all requested marginals.
