%% 
% Cleaning from previous life.
clc
clear

% Set folders with source code.
addpath('./')
addpath('./gui')
addpath('./class')

%% Test classes
% IniDB class
db = db_mysql('ini_sources.mat');
conn = openConnection(db);
schemas = getSchemas(db);

dialog = choose_list_item(schemas,'ionosphere','Choose a schema from popup list:');
choice = dialog.getChoice()

%% Run GUI
GUI = main;
