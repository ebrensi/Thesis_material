function wt = tf_terms_weight(tf_terms)
% wt = tf_terms_weight(tf_terms)
%
% given terms of some expansion of the transfer function, compute a weight
% for each term in the sum, representing the relevance of that term.

wt = sum(abs(tf_terms),2);
wt = wt/sum(wt);

end
