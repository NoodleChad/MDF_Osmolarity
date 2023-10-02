clc
clear
cd('/Users/alexa/OneDrive - University of Toronto/Documents/PhD_UofT/Paper/Osm/'); % Your working directory
addpath('/Program Files/gurobi1002/win64/matlab/');

%% Load the Excel file
    T = readtable('Toy_model.xlsx'); % Your excel file
    S.rxns = T.Properties.VariableNames(2:end);
    S.mets = table2cell(T(1:end,1));
    S.S = table2array(T(1:end,2:end));
    S.S(isnan(S.S))=0;

%% Set-up the Gurobi problem
    model.lb(1:size(S.rxns,2)) = -10;
    model.ub(1:size(S.rxns,2)) = 10;
    model.lb(contains(S.rxns,"C_")) = 0;
    model.ub(contains(S.rxns,"C_")) = 10;
    model.lb(contains(S.rxns,"ln_")) = log(0.1);
    model.ub(contains(S.rxns,"ln_")) = log(1);
    model.lb(contains(S.rxns,"B")) = 0;
    model.ub(contains(S.rxns,"B")) = 10;
    model.ub(contains(S.rxns,"Bio")) = 10;
    model.lb(contains(S.rxns,"G")) = 1;
    model.ub(contains(S.rxns,"G")) = 10;
    model.lb(contains(S.rxns,"b_")) = -10;
    model.ub(contains(S.rxns,"b_")) = -1;
    model.lb(ismember(S.rxns,"K")) = 0;
    model.ub(ismember(S.rxns,"K")) = 10;
    model.lb(ismember(S.rxns,"M")) = 0;
    model.ub(ismember(S.rxns,"M")) = 1;
    model.A = sparse(S.S);
    model.obj(1:size(S.rxns,2)) = 0;
    model.rhs(1:size(S.mets,1)) = 0;
    sense(1:size(S.mets,1)) = {'='};
    sense(contains(S.mets,"_G")) = {'>'};
    sense(contains(S.mets,"LA_")) = {'>'};
    sense(contains(S.mets,"Osm")) = {'<'};
    model.sense = char(join(sense,""));
    vtype(1:size(S.rxns,2)) = {'C'};
    vtype(contains(S.rxns,"k_")) = {'B'};
    model.vtype = char(join(vtype,""));
    model.modelsense = 'min';
    model.varnames = S.rxns;
    model.mets = S.mets;

%% Find minimal K for each solutions in Gurobi
bio = [1 2 3] ;
model.obj(find(ismember(S.rxns,'K'))) = 1; % Smaller K allowing mu
var_results(1:size(S.rxns,2),1)=S.rxns'; % for data visualisation
for i = 1 : size(bio,2)
    model.lb(find(ismember(S.rxns,'Bio'))) = bio(i);
    gurobi_write(model, 'mip1.lp');
    params.outputflag = 0;
    result = gurobi(model, params);
    for j = 1:size(S.rxns,2)
        x(j,1) = {result.x(j)};
    end
    var_results(1:size(S.rxns,2),i+1)=x;
    result.objval
end