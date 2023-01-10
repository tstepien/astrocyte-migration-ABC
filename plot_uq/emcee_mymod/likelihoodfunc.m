function min_quantity = likelihoodfunc(param,filestring,growthfunction,...
    time,number_nodes,exp_px,exp_scale)
% min_quantity = likelihoodfunc(param)
%
% likelihood - minimization function
%
% input:
%   param = parameters (vector with Fk,kb,alpha) - note that Fk = F/k and
%           kb = k/b due to parameter identifiability
%
% output:
%   min_quantity = the quantity that is to be minimized in the parameter
%                  estimation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% forward model
[density,curve,time_simulation_end_early,centerdensincrease,g] ...
    = onelayer_par(param,filestring,growthfunction,time,...
    number_nodes,exp_px,exp_scale);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% calculate error
if time_simulation_end_early==0 && centerdensincrease==0
    min_quantity = data_errorfunction(density,curve,exp_scale,filestring,g);
else
    min_quantity = NaN;
end