% reference 
% https://stats.stackexchange.com/questions/111421/reproduce-linear-discriminant-analysis-projection-plot


% # data taken from your example
X = [-0.9437    -0.0433; -2.4165    -0.5211; -2.0249    -1.0120; ...
    -3.7482 0.2826; -3.3314 0.1653; -3.1927 0.0043; -2.2233 -0.8607; ...
    -3.1965 0.7736; -2.5039 0.2960; -4.4509 -0.3555];
G = [1 1 1 1 1 2 2 2 2 2];

% # overall mean
mu = mean(X);

% # loop over groups
for g=1:max(G)
    mus(g,:) = mean(X(G==g,:)); % # class means
    Ng(g) = length(find(G==g)); % # number of points per group
end

Sw = zeros(size(X,2)); % # within-class scatter matrix
Sb = zeros(size(X,2)); % # between-class scatter matrix
for g=1:max(G)
    Xg = bsxfun(@minus, X(G==g,:), mus(g,:)); % # centred group data
    Sw = Sw + transpose(Xg)*Xg;
    Sb = Sb + Ng(g)*(transpose(mus(g,:) - mu)*(mus(g,:) - mu));
end

St = transpose(bsxfun(@minus,X,mu)) * bsxfun(@minus,X,mu); % # total scatter matrix
assert(sum(sum((St-Sw-Sb).^2)) < 1e-10, 'Error: Sw + Sb ~= St')

% # LDA
[V,D] = eig(Sw\Sb);
[d, ind] = sort(diag(D));
Ds = D(ind,ind);
Vs = V(:,ind);

dataDim = size(X);
Z =  X*Vs;

figure;
colors = {'r' 'b'};
for g=1:max(G)
    scatter3(Z(G==g,1),Z(G==g,2),Z(G==g,3), colors{g}); hold on;
end
xlabel('LD1','FontSize',10);
ylabel('LD2','FontSize',10);
zlabel('LD3','FontSize',10);

            
            
% # projecting data points onto the first discriminant axis
Xcentred = bsxfun(@minus, X, mu);
Xprojected = Xcentred * U(:,1)*transpose(U(:,1));
Xprojected = bsxfun(@plus, Xprojected, mu);
        
        


% # preparing the figure
colors = [1 0 0; 0 0 1];
figure
hold on
axis([-5 0 -2.5 2.5])
axis square

% # plot discriminant axis
plot(mu(1) + U(1,1)*[-2 2], mu(2) + U(2,1)*[-2 2], 'k')
% # plot projection lines for each data point
for i=1:size(X,1)
    plot([X(i,1) Xprojected(i,1)], [X(i,2) Xprojected(i,2)], 'k--')
end
% # plot projected points
scatter(Xprojected(:,1), Xprojected(:,2), [], colors(G, :))
% # plot original points
scatter(X(:,1), X(:,2), [], colors(G, :), 'filled')






