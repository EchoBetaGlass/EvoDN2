function autorun
warning off
%myFun - Description
%
% Syntax autorun
%
% This will be replaced when code added to mastercode
	Problem_name = 'Aman_dataset';
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	EvDtrain.in_index = [1:20];
	EvDtrain.out_index = [21];
	EvDtrain.generations = 500; % 10 max generations for evolution
	EvDtrain.NNet_str = [20 20 20];           %maximum number of nodes
	EvDtrain.Prey_popsize = 800;      %500 Initial popsize
	EvDtrain.no_Prey_preferred = 300; %500 Desired popsize
	EvDtrain.no_new_Prey = 50;       %500 new prey introduced every KillInterval
	EvDtrain.Predator_popsize = 100;  %Number of Predators
	EvDtrain.no_x = 30;              %lattice size (no of rows)
	EvDtrain.no_y = 30;              %lattice size (no of cols)
	EvDtrain.ploton = 1;            %set 0 for no plots or 1 for plots at every generation
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	parameters.EvDtrain = EvDtrain;
    parameters.EvDtrain.name = Problem_name;
    
	Train(Problem_name,parameters);
end