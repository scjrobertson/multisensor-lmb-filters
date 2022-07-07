function [V, W] = generateMultisensorAssociationEvent(L, V, W)
%% Sizes
d = size(L);
m = d(1:end-1) - 1;
n = d(end);
S = length(m);
%% For each sensor
for s = 1:S
    %% For each object
    for i = 1:n
        k = V(i, s) + (V(i, s) == 0);
        u = [(V(i, :) + 1),  i];
        %% For each of sensor s' measurments
        for j = k:m(s)
            %% Update association vectors
            if ((W(j, s) == 0) || (W(j, s) == i))
                %% Determine sample probability
                % Generated
                u(s) = j + 1;
                q =  determineLinearIndex(u, d);
                % Missed
                u(s) = 1;
                r =  determineLinearIndex(u, d);
                % Sample probability
                P = 1 / (exp(L(r) - L(q)) + 1 ); % L is the loglikelihood of the event
                %% Sample
                if (rand() < P)
                    %% Object i generated measurement z_j^s
                    V(i, s) = j;
                    W(j, s) = i;
                    break;
                else
                    %% Object i did generate measurement z_j^s
                    V(i, s) = 0;
                    W(j, s) = 0;
                end
            end
        end
    end
end
end
%% Determine a linear index using an association vector
function ell = determineLinearIndex(u, d)
ell = u(1);
Pi = 1;
for i = 2:length(u)
    Pi = Pi * d(i-1);
    ell  = ell  + Pi * (u(i) - 1);
end
end