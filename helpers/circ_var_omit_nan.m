function x = circ_var_omit_nan(these_oris)

these_oris_noNaN    = these_oris(~isnan(these_oris));

if ~isempty(these_oris_noNaN)
    x  = circ_var(these_oris_noNaN);
else
    x = NaN;
end