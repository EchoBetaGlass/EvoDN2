function PPGA(parameters)

%%Descrition: Use PPGA to create DNN models

%NNGA definitions
input_nodes = parameters.input_nodes; 	%Cell defining the input varables of 
										%various subnets
NNet_str = parameters.NNet_str; % %maximum number of nodes
P_omit_max = 0.99; %After ranking max value for random generated P_omit
P_omit_min = 0.3; %-''- min value 
wlow = -5; %Lower bound for randomly generated weights
whigh = 5; %Upper bound for randomly generated weights

%PP definitions
Prey_popsize = parameters.Prey_popsize;                %Initial popsize
no_Prey_preferred = parameters.no_Prey_preferred;      %Desired popsize
Predator_popsize = parameters.Predator_popsize;
no_new_Prey = parameters.no_new_Prey;
no_generations = parameters.generations;
P_move_Prey = 0.3;

%2D-lattice
no_x = parameters.no_x;
no_y = parameters.no_y;

%Killing of bad Prey
KillInterval = parameters.KillInterval;
maxrank = parameters.maxrank;
P_omit_match = P_omit_min;
F_bad = 1e6;%fitness assigned to Preys performing badly

Prey = create_population(Prey_popsize, input_nodes, NNet_str);

%Placement of Prey
for i = 1:length(Prey)
    [emptyx,emptyy] = find(lattice(2:no_x+1,2:no_y+1) == 0);
    if ~isempty(emptyx > 0)
        j = ceil(rand*(length(emptyx)));
        lattice(emptyx(j)+1,emptyy(j)+1) = i;
    end
end

%initialization of Predator_pop
Predators = (0:1/(Predator_popsize-1):1)'; %value to weight F1
%location of predators
for i = 1:length(Predators(:,1))
    [emptyx,emptyy] = find(lattice(2:no_x+1,2:no_y+1) == 0);
    if ~isempty(emptyx > 0)
        j = ceil(rand*(length(emptyx)));
        lattice(emptyx(j)+1,emptyy(j)+1) = -i; %Predator indentified with neg values
    end
end

MovePrey(0, 0, 0);

