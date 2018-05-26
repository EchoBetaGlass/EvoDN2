function Population = create_Population(Pop_size, Pop_str)
%create_Poulation - Creates cell array of Population of subnet NNs
%
% Syntax: Population = create_Poulation(Pop_size, Pop_str)
%
% Pop_size descries the population size
% Pop_str is a cell array with n cells. Each cell describes a subnet.
% Pop_str{i}{1} is a vector of the column indices of the input
% Pop_str{i}{2} is a vector describing the number of nodes in each layer
    
    P_omit_match = 0.3; % Probability with which a connection is randomly killed.
    whigh = 5;  %Max value for a connection
    wlow = -5;  %min value for a connection

    Population{Pop_size} = {};
    Num_subnets = length(Pop_str);
    for Pop_index = 1:Pop_size
        Population(Pop_index).structure = Pop_str;
        for s_index = 1:Num_subnets
            for layer = 1:length(Pop_str{s_index}{2})-1
                in = Pop_str{s_index}{2}(layer);    % number of input nodes
                out = Pop_str{s_index}{2}(layer+1); % number of output nodes
                % random genration of layer
                net = rand(in+1,out)*(whigh-wlow)+wlow;   % +1 for bias node
                % Killing some connections
                for i = 1:in+1
                    for j = 1:out
                        if rand(1,1) < P_omit_match
                            net(i,j) = 0;
                        end
                    end
                end
                Population(Pop_index).subnet{s_index}{layer} = net;
            end
        end
    end
end