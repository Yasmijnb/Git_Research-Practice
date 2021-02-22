clear all
close all

%% Retrieve the adjacency matrices
Ad1 = readmatrix('Adjacency_matrix_men.xlsx');
Ad2 = readmatrix('Adjacency_matrix_women.xlsx');
Ad3 = readmatrix('Adjacency_matrix_young.xlsx');
Ad4 = readmatrix('Adjacency_matrix_old.xlsx');
Ad5 = readmatrix('Adjacency_matrix_young_men.xlsx');
Ad6 = readmatrix('Adjacency_matrix_old_men.xlsx');
Ad7 = readmatrix('Adjacency_matrix_young_women.xlsx');
Ad8 = readmatrix('Adjacency_matrix_old_women.xlsx');

%% Make input for COVSCA
COVSCAinput = [Ad1 Ad2 Ad3 Ad4 Ad5 Ad6 Ad7 Ad8];

%% Input parameters
% Number of analyses
nanal = 1000; % Larger number will give more accurate resutls

% The number of loadings for each component
Q = [3 3]';
% The number of dimensions
L = length(Q);

%% Run COVSCA
[loadings, scores, fp, dys, func] = covsca(COVSCAinput, L, Q, 1, 1, nanal);

% Fit percentages
disp(fp) % Value between 0 and 100

%% Plot scores and Loadings
figure(1)
labels = {'all','all','young','old','young','old','young','old'};
plot(scores(1,1),scores(1,2),'b.', ...
    scores(2,1),scores(2,2),'r.', ...
    scores(3,1),scores(3,2),'k.', ...
    scores(4,1),scores(4,2),'k.', ...
    scores(5,1),scores(5,2),'b.', ...
    scores(6,1),scores(6,2),'b.', ...
    scores(7,1),scores(7,2),'r.', ...
    scores(8,1),scores(8,2),'r.', ...
    'MarkerSize',22);
text(scores(:,1),scores(:,2),labels,'VerticalAlignment','bottom','HorizontalAlignment','left')
leg = legend('men', 'women', 'all');
legtitle = get(leg, 'Title');
set(legtitle,'String','Gender')
xlabel('1st COVSCA component','FontSize',13);
ylabel('2nd COVSCA component','FontSize',13);
title('COVSCA scores (weights)','FontSize',16);

%% Plot scores and Loadings
figure(2)
subplot(4,1,1)
bar(loadings(:,1)');
ylabel('Loadings 1st comp','FontSize',13);
title('COVSCA loadings','FontSize',16);
subplot(4,1,2)
bar(loadings(:,2)');
ylabel('Loadings 2nd comp','FontSize',13);
subplot(4,1,3)
bar(loadings(:,3)');
ylabel('Loadings 3rd comp','FontSize',13);
subplot(4,1,4)
bar(loadings(:,4)');
xlabel('Variables','FontSize',13);
ylabel('Loadings 4th comp','FontSize',13);
