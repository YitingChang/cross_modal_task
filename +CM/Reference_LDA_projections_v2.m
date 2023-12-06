load fisheriris
Mdl = fitcdiscr(meas,species, 'DiscrimType','quadratic')

figure, % nexttile
gscatter(meas(:,1), meas(:,2), species);
hold on
scatter(Mdl.Mu(:, 1), Mdl.Mu(:, 2), 100, 'xK', 'LineWidth', 1.5)
hold off
% nexttile
% gscatter(meas(:,1), meas(:,3), species);
% hold on
% scatter(Mdl.Mu(:, 1), Mdl.Mu(:, 3), 100, 'xK', 'LineWidth', 1.5)
% hold off
% legend('off')
% nexttile
% gscatter(meas(:,1), meas(:,4), species);
% hold on
% scatter(Mdl.Mu(:, 1), Mdl.Mu(:, 4), 100, 'xK', 'LineWidth', 1.5)
% hold off
% legend('off')
% nexttile
% gscatter(meas(:,2), meas(:,3), species);
% hold on
% scatter(Mdl.Mu(:, 2), Mdl.Mu(:, 3), 100, 'xK', 'LineWidth', 1.5)
% hold off
% legend('off')

V = Mdl.Mu; % positions of the centers
Comb = [1, 2; 1, 3; 2, 3]; % all possible combinations to form the axis with 3 groups
for ii = 1:3 % for each axis
    V1 = V(Comb(ii,2), :);
    V2 = V(Comb(ii,1), :);
    v = (V2-V1)./norm(V2-V1); % normalised vector for projection 
    
    for jj =1:size(meas, 1) % observations
        Q(jj, :) = dot(meas(jj, :)-V1,v)*v+V1;
        LD(jj, ii) = (Q(jj, 1)-V1(1))/(V2(1)-V1(1));
    end
end

figure
hold on
for ii = 1:3
    scatter3(LD(strcmp(species, 'setosa'), 1), LD(strcmp(species, 'setosa'), 2), LD(strcmp(species, 'setosa'), 3), 'r')
    scatter3(LD(strcmp(species, 'versicolor'), 1), LD(strcmp(species, 'versicolor'), 2), LD(strcmp(species, 'versicolor'), 3), 'g') 
    scatter3(LD(strcmp(species, 'virginica'), 1), LD(strcmp(species, 'virginica'), 2), LD(strcmp(species, 'virginica'), 3), 'b') 
end
hold off
legend({'setosa', 'versicolor', 'virginica' })
xlabel('LD1')
ylabel('LD2')
zlabel('LD3')
view(40,35)

% two predictors and two groups
Mdl = fitcdiscr(meas(1:100,1:2),species(1:100), 'DiscrimType','quadratic')
figure
gscatter(meas(1:100,1), meas(1:100,2), species(1:100));
hold on
scatter(Mdl.Mu(:, 1), Mdl.Mu(:, 2), 100, 'xK', 'LineWidth', 1.5)
hold off

V = Mdl.Mu; % positions of the centers
V1 = V(2,:);
V2 = V(1,:);
v = (V2-V1)./norm(V2-V1); % normalised vector for projection

for jj =1:size(meas(1:100), 1) % observations
    Q(jj, :) = dot(meas(jj, 1:2)-V1,v)*v+V1; % the projected position on the LD axis in the original coordinates 
    LD(jj, ii) = (Q(jj, 1)-V1(1))/(V2(1)-V1(1)); % the projected position in the LD coordinates 
end