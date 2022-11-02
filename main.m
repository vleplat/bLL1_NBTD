%      Nonnegative block-term decomposition with the beta-divergence:     %
%            joint data fusion and blind spectral unmixing                %
%-------------------------------------------------------------------------%

% Copyright (c) 2022 Clemence Prevost, Valentin Leplat
% https://github.com/cprevost4/bLL1_NBTD
% Contact: clemence.prevost@univ-lille.fr

% This software reproduces the results from the preprint called:
% "Nonnegative block-term decomposition with the beta-divergence:
% joint data fusion and blind spectral unmixing" - C.Prévost, V. Leplat

% In order to run the demos, you will need to add to your MATLAB path:
% - Tensorlab 3.0: https://www.tensorlab.net
% - The HSMS data fusion toolbox: https://openremotesensing.net/knowledgebase/hyperspectral-and-multispectral-data-fusion/
% - The tensor HSR toolbox: https://drive.google.com/file/d/1dwLwwNTseMGIGhxS432BC3PREkJSMlo7/view?usp=share_link

%-------------------------------------------------------------------------%
%                              CONTENT
% - /data : contains data for degradation matrices and random vectors
% - /demos : contains demo files that produce tables and figures
% - /figures : where the tables and figures are saved
% - /src : contains helpful files to run the demos
%-------------------------------------------------------------------------%
%                               MENU
% You can launch a specified demo by typing its number. The resulting tables
% and figures produced will be stored in the figures folder.
%
% 1:  A demo with minimal requirements
% 2:  produces Fig. 2-3 and Table 1
% 3:  produces Fig. 4-5 and Table 2
% 4:  produces Fig. 6 and Table 3
%-------------------------------------------------------------------------%

list_demos = ["demo_minimal_requirements","demo_synthetic","demo_semireal",...
    "demo_semireal_semiblind"];

prompt = "Which file do you want to run ?";
num = input(prompt);
eval(list_demos(num));


