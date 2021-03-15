clear all
close all

%% Retrieve the adjacency matrices
Ad1 = readmatrix('\Adjacency_matrix_men.xlsx');
Ad2 = readmatrix('\Adjacency_matrix_women.xlsx');
% Ad3 = readmatrix('\Adjacency_matrix_young.xlsx');
% Ad4 = readmatrix('\Adjacency_matrix_old.xlsx');
Ad5 = readmatrix('\Adjacency_matrix_youngmen.xlsx');
Ad6 = readmatrix('\Adjacency_matrix_oldmen.xlsx');
Ad7 = readmatrix('\Adjacency_matrix_youngwomen.xlsx');
Ad8 = readmatrix('\Adjacency_matrix_oldwomen.xlsx');

%% Make input for COVSCA
COVSCAinput = [Ad1 Ad2 Ad5 Ad6 Ad7 Ad8];

%% Input parameters
% Number of analyses
nanal = 1000; % Larger number will give more accurate resutls

% The number of loadings for each component
Q = [2 2 2]';
% The number of dimensions
L = length(Q);

%% Run COVSCA
[loadings, scores, fp, dys, func] = covsca(COVSCAinput, L, Q, 1, 1, nanal);

% Fit percentages
disp(fp) % Value between 0 and 100

%% Plot scores and Loadings
figure(1)
labels = {'all','all','young','old','young','old'};
plot(scores(1,1),scores(1,2),'b.', ...
    scores(2,1),scores(2,2),'r.', ...
    scores(3,1),scores(3,2),'b.', ...
    scores(4,1),scores(4,2),'b.', ...
    scores(5,1),scores(5,2),'r.', ...
    scores(6,1),scores(6,2),'r.', ...
    'MarkerSize',22);
text(scores(:,1),scores(:,2),labels,'VerticalAlignment','bottom',...
    'HorizontalAlignment','left')
leg = legend('men', 'women');
legtitle = get(leg, 'Title');
set(legtitle,'String','Gender')
xlabel('1st COVSCA component','FontSize',13);
ylabel('2nd COVSCA component','FontSize',13);
title('COVSCA scores (weights)','FontSize',16);

%% Plot scores and Loadings

% Assign the variable names
variable_labels = {'Triglycerides, VLDL', 'Triglycerides, IDL', ...
    'Triglycerides, LDL', 'Triglycerides, HDL', 'Cholesterol, VLDL', ...
    'Cholesterol, IDL', ' Cholesterol, LDL', 'Cholesterol, HDL', ...
    'Free Cholesterol, VLDL', 'Free Cholesterol, IDL', ...
    'Free Cholesterol, LDL', 'Free Cholesterol, HDL', ...
    'Phospholipids, VLDL', 'Phospholipids, IDL', 'Phospholipids, LDL', ...
    'Phospholipids, HDL', 'Apo-A1, HDL', 'Apo-A2, HDL', 'Apo-B, VLDL', ...
    'Apo-B, IDL', ' Apo-B, LDL'};

% Show two loadings at a time
figure(2)
set(gcf, 'color', 'w');
title('COVSCA loadings','FontSize',16);
xlabel('Variables','FontSize',13);
% First
subplot(2,1,1)
bar(loadings(:,1)');
ylabel('Loadings 1st comp','FontSize',13);
set(gca, 'xtick', [1:21], 'xticklabel', variable_labels);
xtickangle(45);
% Second
subplot(2,1,2)
bar(loadings(:,2)');
ylabel('Loadings 2nd comp','FontSize',13);
set(gca, 'xtick', [1:21], 'xticklabel', variable_labels);
xtickangle(45);

% Show two loadings at a time
figure(3)
set(gcf, 'color', 'w');
title('COVSCA loadings','FontSize',16);
xlabel('Variables','FontSize',13);
% First
subplot(2,1,1)
bar(loadings(:,3)');
ylabel('Loadings 1st comp','FontSize',13);
set(gca, 'xtick', [1:21], 'xticklabel', variable_labels);
xtickangle(45);
% Second
subplot(2,1,2)
bar(loadings(:,4)');
ylabel('Loadings 2nd comp','FontSize',13);
set(gca, 'xtick', [1:21], 'xticklabel', variable_labels);
xtickangle(45);
