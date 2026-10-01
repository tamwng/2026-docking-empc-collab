% setup.m — run once per VS Code session to load project paths
proj_root = fileparts(mfilename('fullpath'));
proj = openProject(fullfile(proj_root, '2026-docking-empc-collab.prj'));
cd(proj_root);   % openProject may change cwd; restore to project root

addpath(genpath(fullfile(proj_root, 'simulations')));
addpath(genpath(fullfile(proj_root, 'tests')));

disp('Project loaded.');