%% Confusion Pie Chart for Discord Usability Test
% Edit the values in 'counts' to match your data

clc; clear; close all;

% === Data: number of participants who reported confusion in each category ===
% For example, with 2 testers:
%   - 2 confused by pop-ups
%   - 2 confused by icons/layout
%   - 1 confused by terminology
%   - 1 confused by friend-adding format
counts = [2 2 1 1];

% Labels for each confusion category
labels = { ...
    'Pop-ups & ads', ...
    'Server icons & layout', ...
    'Terminology (threads, tags, etc.)', ...
    'Friend-adding format (username + tag)'};

% === Create the pie chart ===
figure;
p = pie(counts);

% Make labels show percentages + category names
% (This automatically uses the 'labels' text for each slice)
colormap(parula(numel(counts))); % you can change colormap if you want
title('Sources of User Confusion in Discord Usability Test');

% Add legend with category labels
legend(labels, 'Location', 'eastoutside');

% Improve readability
set(gca, 'FontSize', 12);
