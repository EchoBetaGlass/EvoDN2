function PPGA(parameters)

%%Descrition: Use PPGA to create DNN models

%NNGA definitions
NNet_str = parameters.NNet_str; %Cell defining the input varables and 
								%structure of various subnets
P_omit_max = 0.99;  %After ranking max value for random generated P_omit
P_omit_min = 0.3;   %-''- min value 

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

Prey = create_Population(Prey_popsize, NNet_str);

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

for num_Gen = 1:no_generations
    Pop_size = length(Prey);

    % Calculate Standard deviation in population as mutval here%

    % Movement of all population members
    for Pop_index = 1 : Pop_size
        if rand < P_move_prey
            %Identify location
            [xpos,ypos] = find(lattice(2:no_x+1,2:no_y+1) == i);
            xpos = xpos+1; ypos = ypos+1;
            for j = 1:10 %10 trial for free spot
                dx = round(rand*(3-eps)-1.5); %-1 to the left, 0 no move, 1 to the right
                dy = round(rand*(3-eps)-1.5); %-1 down, 0 no move, 1 up
                if lattice(xpos+dx,ypos+dy) == 0
                    lattice(xpos,ypos) = 0; %removing prey i
                    MovePrey(xpos+dx, ypos+dy, i); 
                    break
                end
            end
        end
    end

    % Breeding of all population members
    for Pop_index = 1 : Pop_size
        [xpos,ypos] = find(lattice(2:no_x+1,2:no_y+1) == i);
        xpos = xpos+1; ypos = ypos+1;
        moore = lattice(xpos-1:xpos+1,ypos-1:ypos+1);
        %Remove the i-prey
        moore(2,2) = 0;
        [matex,matey] = find(moore >= 1 & moore <= len);
        if ~isempty(matex)
            parent2 = ceil(rand*length(matex));
            parent2 = lattice(xpos-2+matex(parent2),ypos-2+matey(parent2));
            % Write the following function
            Offsprng = create_Offsprng(?);
            %Random placement of offsprings /10 trials
            for l = 1:2
                for j=1:10 %10 trial for free spot
                    xpos = round(rand(1,1)*(no_x-1)+1)+1;
                    ypos = round(rand(1,1)*(no_y-1)+1)+1;
                    if lattice(xpos,ypos)==0
                        Prey = [Prey Offsprng(l)];
                        lattice(xpos,ypos) = length(Prey{1}(:,1,1));
                        break
                    end
                end
            end
        end
        MovePrey(0, 0, 0);
    end 