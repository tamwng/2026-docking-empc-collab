% setup.m — run once per VS Code session to load project paths
proj = openProject(fullfile(fileparts(mfilename('fullpath')), '2026-docking-empc-collab.prj'));
disp('Project loaded.');