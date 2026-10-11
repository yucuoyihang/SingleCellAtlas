% ENIGMA GitHub repository: https://github.com/FabricMC/Enigma
addpath(genpath(fullfile('ENIGMA', 'matlab')))

ct = readtable('Cortical.csv');
ct_d_fsa5 = parcel_to_surface(ct.t, 'aparc_fsa5');
f1 = figure,
    plot_cortical(ct_d_fsa5, ...
        'surface_name', 'fsa5', ...
        'label_text', 't-value', ...
        'color_range', [-2 2], ...
        'cmap', 'RdBu_r');
print(f1, 'Cortical', '-dpdf', '-r300')

sv = readtable('Subcortical.csv');

f2 = figure,
    plot_subcortical(sv.t, ...
        'color_range', [-2 2], ...
        'cmap', 'RdBu_r');
print(f2, 'Subcortical', '-dpdf', '-r300')