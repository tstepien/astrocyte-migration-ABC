clear variables
clc

pn = [18;13;11;9;9];
ninetype = {'bio','uni'};
P = length(pn);
pnstring = {'18';'13';'11';'9bio';'9uni'};


nineiter = 1;

threshold = [10^4 10 5 4];
T = length(threshold);
num_threshold = zeros(P,T);
ind = (1:5e5)';

%% load parameter distribution information
for j=1:P
    if pn(j)==9
        load(strcat('parameter_analysis/abc',num2str(pn(j)),...
            ninetype{nineiter},'_5e5.mat'));
        nineiter = nineiter + 1;
    else
        load(strcat('parameter_analysis/abc',num2str(pn(j)),...
            '_5e5.mat'));
    end

    for i=1:T
        ind_threshold = ind(err_tot < threshold(i));
        num_threshold(j,i) = length(ind_threshold);
        clear ind_threshold
    end
end

%% Bayes Factor calculation
bayesfactor = zeros(P,P,T);
for k=1:T
    for i=1:P
        for j=1:P
            if i==j
                bayesfactor(i,j,k) = 1;
            else
                bayesfactor(i,j,k) = num_threshold(i,k)/num_threshold(j,k);
            end
        end
    end
end

disp(bayesfactor)