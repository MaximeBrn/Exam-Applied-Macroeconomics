function ds = dynamic_set_auxiliary_series(ds, params)
%
% Status : Computes Auxiliary variables of the dynamic model and returns a dseries
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

ds.AUX_ENDO_LAG_95_1=ds.E_ZEPS_PPI(-1);
ds.AUX_ENDO_LAG_95_2=ds.AUX_ENDO_LAG_95_1(-1);
ds.AUX_ENDO_LAG_95_3=ds.AUX_ENDO_LAG_95_2(-1);
end
