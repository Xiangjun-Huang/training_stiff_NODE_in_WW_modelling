%% Continuity check
% To conform the conservation law, the absolute result of Stoichiometric Matrix 
% (40 x 24) multiplying Composition Matrix (24 x 5) must be zero or less than 
% 1E-15.

Con_Mat = Stoi_Mat * Comp_Mat;
if ~isempty(Con_Mat(Con_Mat>1E-15))
    error("Continuity Check failed! Please check and correct your parameters input!");
end