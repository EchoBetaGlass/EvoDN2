function PP_NNGA(train_set_no)

global lattice no_x no_y
global generation
global parameters setno F_bad
global prey

setno = train_set_no;

%NNGA definitions
NNet_str = parameters.NNet_str; % %maximum number of nodes
noinnodes = parameters.noinnodes;
nooutnodes = parameters.nooutnodes;
P_omit_max = 0.99; %After ranking max value for random generated P_omit
P_omit_min = 0.3; %-''- min value 
wlow = -5; %Lower bound for randomly generated weights w
whigh = 5; %Upper bound for randomly generated weights w

%PP definitions
Prey_popsize = parameters.Prey_popsize;                %Initial popsize
no_Prey_preferred = parameters.no_Prey_preferred;      %Desired popsize
Predator_popsize = parameters.Predator_popsize;
no_new_Prey = parameters.no_new_Prey;
no_generations = parameters.generations;
P_move_prey = 0.3;
%2D-lattice
no_x = parameters.no_x;
no_y = parameters.no_y;
KillInterval = parameters.KillInterval;
maxrank = parameters.maxrank;

P_omit_match = P_omit_min;
F_bad = 1e6;        %fitness assigned to preys performing badly

%**************************************************
%   INITIALIZATION
lattice = zeros(no_x+2,no_y+2);
%initialization of Prey_pop
for subnet = 1:num_subnets
    NNet_str{subnet} = {length(input_nodes{subnet}) NNet_str};
    for layer = 1:(length(NNet_str)-1)
        Prey{subnet}{layer} = rand(Prey_popsize,NNet_str(layer)+1,NNet_str(layer+1)) ...
                    *(whigh-wlow)+wlow;
        %Eliminating some matches
        for i = 1:Prey_popsize
            for j = 1:NNet_str(layer)+1
                for k = 1:NNet_str(layer+1)
                    if rand(1,1) < P_omit_match
                        Prey{subnet}{layer}(i,j,k) = 0;                
                    end
                end
            end
        end
    end
%location of prey
for i = 1:length(Prey{1})
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

%**************************************************
%   GENERATIONS
for generation=1:no_generations
    Prey_new = [];
    len = length(Prey{1}(:,1,1));
    for layer = 1:(length(NNet_str)-1)
        mutval{layer} = squeeze(std(Prey{layer}));
    end
    % Move of prey /10 trials /1step
    for i = 1:len
        if rand < P_move_prey
            %Identify location
            [xpos,ypos] = find(lattice(2:no_x+1,2:no_y+1) == i);
            xpos = xpos+1; ypos = ypos+1;
            for j = 1:10 %10 trial for free spot
                dx = round(rand*(3-eps)-1.5); %-1 to the left, 0 no move, 1 to the right
                dy = round(rand*(3-eps)-1.5); %-1 down, 0 no move, 1 up
                if lattice(xpos+dx,ypos+dy) == 0
                    lattice(xpos,ypos) = 0; %removing prey i
                    MovePrey(xpos+dx, ypos+dy, i); break
                end
            end
        end
    end
    
    % Breeding of prey
    for i = 1:len
        [xpos,ypos] = find(lattice(2:no_x+1,2:no_y+1) == i);
        xpos = xpos+1; ypos = ypos+1;
        moore = lattice(xpos-1:xpos+1,ypos-1:ypos+1);
        %Remove the i-prey
        moore(2,2) = 0;
        [matex,matey] = find(moore >= 1 & moore <= len);
        if ~isempty(matex)
            parent2 = ceil(rand*length(matex));
            parent2 = lattice(xpos-2+matex(parent2),ypos-2+matey(parent2));
            Offsprng = create_offspring(Prey, i,parent2,mutval,generation,no_generations);
            %Random location of offspring /10 trials
            %Placing offspring
            for l = 1:2
                for j=1:10 %10 trial for free spot
                    xpos = round(rand(1,1)*(no_x-1)+1)+1;
                    ypos = round(rand(1,1)*(no_y-1)+1)+1;
                    if lattice(xpos,ypos)==0
                        for layer = 1:(length(NNet_str)-1)
                            Prey{layer} = [Prey{layer}; Offsprng{layer}(l,:,:)];
                        end
                        lattice(xpos,ypos) = length(Prey{1}(:,1,1));
                        break
                    end
                end
            end
        end
        MovePrey(0, 0, 0);
    end    

    [F1 F2 UW IC] = EvalF1F2(Prey); 
    [fonrank front] = NONDOM_SORT([F1 F2]);
    
    % Removing all except rank n
    if generation/KillInterval == round(generation/KillInterval)
        if generation < no_generations
            for layer = 1:(length(NNet_str)-1)
                Prey_new{layer} = rand(no_new_Prey,NNet_str(layer)+1,NNet_str(layer+1))*(whigh-wlow)+wlow;
                %Prey = rand(Prey_popsize,nonodes,noinnodes+1)*(whigh-wlow)+wlow;
                %Eliminating some matches
                for i = 1:Prey_popsize
                    for j = 1:NNet_str(layer)+1
                        for k = 1:NNet_str(layer+1)
                            if rand(1,1) < P_omit_match
                                Prey_new{layer}(i,j,k) = 0;                
                            end
                        end
                    end
                end
            end
        end
        indfr = find(fonrank > maxrank);
        F1(indfr) = F_bad+eps; F2(indfr) = F_bad+eps;
        [Prey F1 F2] = KillBadPrey(Prey, F1, F2);
        [fonrank front] = NONDOM_SORT([F1 F2]);
    end
    
    % Move of predators /killing
    crodit = CROW_SORT([F1 F2], front);
    F1 = (F1 - min(F1))./(max(F1) - min(F1));
    F2 = (F2 - min(F2))./(max(F2) - min(F2));
    PredMoves = floor((length(Prey{1}) - no_Prey_preferred)/Predator_popsize);
    fprintf('\nGeneration %i: Predatorpop %i PredMoves %i\n',generation,length(Predators(:,1)),PredMoves);
    fprintf('Preypop before: %i; Preypop after: ', length(Prey{1}))
    for i = 1:length(Predators)
        for k = 1:PredMoves
            [xpos, ypos] = find(lattice(2:no_x+1,2:no_y+1) == -i);
            xpos = xpos+1; ypos = ypos+1;
            [matex,matey] = find(lattice(xpos-1:xpos+1,ypos-1:ypos+1) > 0);
            if length(matex) > 1 %prey available
                rows=[];
                for t = 1:length(matex)
                    rows = [rows; lattice(matex(t)+xpos-2,matey(t)+ypos-2)];
                end
                f = (Predators(i)*F1(rows) + (1-Predators(i))*F2(rows)).*(front(rows)-1); % elitism for front = 1
                if length(unique(front(rows))) == 1
                    f = (f+1)./(1+crodit(rows));
                end
                [~,pos] = max(f);
                j = rows(pos(1)); lattice(xpos,ypos) = 0; %removing predator i
                [xpos,ypos] = find(lattice(2:no_x+1,2:no_y+1) == j);
                xpos = xpos+1; ypos = ypos+1;
                Prey = MovePredator(Prey, xpos, ypos, i, j);
                F1(j) = []; F2(j) = []; front(j) = []; crodit(j) = []; fonrank(j) = [];  %CHECK!!!!! => checked, works => see MovePredator
            else %Only move
                for j=1:10 %10 trial for free spot
                    dx=round(rand*(3-eps)-1.5); %-1 to the left, 0 no move, 1 to the right
                    dy=round(rand*(3-eps)-1.5); %-1 down, 0 no move, 1 up
                    if lattice(xpos+dx,ypos+dy)==0
                        lattice(xpos,ypos)=0; %removing predator i
                        Prey = MovePredator(Prey, xpos+dx, ypos+dy, i, inf);
                        break
                    end
                end
            end
        end
    end
    fprintf('%i\n' , length(Prey{1}))
    
    [F1 F2 UW IC] = EvalF1F2(Prey);
    PlotLattice
    PlotPareto(F1, F2, fonrank)
    
    % Placing new Prey
    for i = 1:length(Prey_new)
        for k = 1:10 %10 trial for free spot
            [emptyx,emptyy] = find(lattice(2:no_x+1,2:no_y+1) == 0);
            if ~isempty(emptyx > 0)
                j = ceil(rand*length(emptyx));
                lattice(emptyx(j)+1,emptyy(j)+1) = length(Prey{1}(:,1,1)) + 1;
                for layeer = 1:length(Prey)
                Prey{layer} = [Prey{layer}; Prey_new{layer}(i,:,:)];
                end
                break
            end
        end
    end
end

[F1 F2 UW IC] = EvalF1F2(Prey);
[fonrank, front] = NONDOM_SORT([F1 F2]);
PlotPareto(F1, F2, fonrank); figure(4)
for layer = 1:length(Prey)
    Prey{layer} = Prey{layer}(front == 1,:,:);
end
F1 = F1(front == 1); F2 = F2(front == 1);
UW = UW(front == 1); IC = IC(front == 1);
f = [(1:length(F1))' F1 F2]; f = sortrows(f, 3); f = sortrows(f, 2); i = 2;
while i ~= length(f(:,1))
    if f(i,2) == f(i-1,2), f(i,:) = []; else i = i+1; end
end
F1 = F1(f(:,1)); F2 = F2(f(:,1)); UW = UW(f(:,1));
IC = IC(f(:,1)); [~,i] = min(IC); hold on
parameters.dataset(setno).pareto.select = i;
plot(F2(i), F1(i), 'dr', 'LineWidth', 2)

P = [];
for i = 1:length(f(:,1))
    for layer = 1:(length(NNet_str)-1)
        P{layer}.w(i,:,:) = Prey{layer}(f(i,1),:,:);
    end
    info(i).W = UW(i).W; info(i).IC = IC(i);
end

parameters.dataset(setno).pareto.P = P;
parameters.dataset(setno).pareto.info = info;
parameters.dataset(setno).pareto.F1 = F1;
parameters.dataset(setno).pareto.F2 = F2;

end
%-------------------------------------------------
% Evaluate F1 and F2
%-------------------------------------------------
% function[F1 F2 UW IC] = EvalF1F2(Prey)
% global nonodes noinnodes F_bad

% for i = 1:length(Prey(:,1,1))
%     for j = 1:nonodes
%         for k = 1:noinnodes+1
%             w(j,k) = Prey(i,j,k);
%         end
%     end
%     [fval,W,InfoC] = NNevalnet(w);
%     F1(i) = fval;
%     if isnan(F1(i)) || isinf(F1(i))
%         F1(i) = F_bad+eps; F2(i) = F_bad+eps;
%     end
%     F2(i) = length(find(w(:,2:noinnodes+1))); %Does not include bias parameters!
%     UW(i).W = W;
%     IC(i) = InfoC;
%     %Should somehow be able to retrieve the different networks...
% end
% F1 = F1'; F2 = F2'; IC = IC'; UW = UW';
% end
function[F1 F2 UW IC] = EvalF1F2(Prey)
global nonodes noinnodes F_bad
F1 = zeros(length(Prey{1}),1);
F2 = F1; IC = F1;
for i = 1:length(Prey{1})
    w = {};
    for layer = 1:length(Prey)
        w{layer} = squeeze(Prey{layer}(i,:,:));  
    end
    [fval,complexity,W,InfoC] = DNevalnet(w);
    F1(i) = fval;
    if isnan(F1(i)) || isinf(F1(i))
        F1(i) = F_bad+eps; F2(i) = F_bad+eps;
    end
    %F2(i) = length(find(w(:,2:noinnodes+1))); %Does not include bias parameters!
    F2(i) = complexity;
    UW(i).W = W;
    IC(i) = InfoC;
    %Should somehow be able to retrieve the different networks...
end
end