clear variables
clc

pn = [18;13;11;9];
ninetype = {'bio','uni'};
P = length(pn);

% number of data points: 7 APC radii, 7 IPA radii, 1 time point = 15
numdata = 15;

param_val = zeros(P,18);
min_error = zeros(P,1);
AIC_val = zeros(P,1);
AICc_val = zeros(P,1);

nineiter = 1;

for j=1:P
    if pn(j)==9
        load(strcat('parameter_analysis/bestfitparam',num2str(pn(j)),...
            ninetype{nineiter},'.mat'));
        nineiter = nineiter + 1;
    else
        load(strcat('parameter_analysis/bestfitparam',num2str(pn(j)),...
            '.mat'));
    end
    min_error(j) = err_output;

    AIC_val(j) = 2*pn(j) - 2*log(min_error(j));
    AICc_val(j) = AIC_val(j) + 2*pn(j)*(pn(j)+1)/(numdata-pn(j)-1);
    
    clear ind
end

table(pn,min_error,AIC_val,AICc_val)