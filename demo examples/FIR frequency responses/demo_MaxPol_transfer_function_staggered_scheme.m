clear all
close all
clc

%%
current_directory = pwd;
cd ..
cd ..
addpath([cd, filesep, 'utilities'])
cd(current_directory)

%%  parameter initialization
disp('(1) lowpass')
disp('(2) fullband non-zero side-shift')
disp('(3) fullband zero side-shift (s=0)')
experiment_number = input('Please input the experiment number for demosntration: ');

switch experiment_number
    case 1
        sym_flag = true;
        tap_length_polynmial = 15;
        tap_length = tap_length_polynmial*ones(1,tap_length_polynmial);
        possible_cutoffs = [1: 2: 2*tap_length_polynmial-1];
        side_shifts = [0];
        n_ords = [0, 1, 2, 3, 4];
        col_bumbers = tap_length_polynmial;
    case 2
        sym_flag = false;
        tap_length = [5];
        possible_cutoffs = [2*tap_length-1];
        side_shifts = [0: tap_length];
        n_ords = [0, 1, 2, 3, 4];
        col_bumbers = tap_length+1;
    case 3
        sym_flag = false;
        tap_length = [1:15];
        possible_cutoffs = [2*tap_length-1];
        side_shifts = [0];
        n_ords = [0, 1, 2, 3, 4];
        col_bumbers = numel(tap_length);
end

%%  parameter for frequency response (fft) calculation
n_fft = 513;
stp = 2*pi/(n_fft);
omega = [-pi: stp: pi-stp];

%%  Color specification
col = colormap(hot);
col = col(1:end-20, :);
stp = length(col)/(col_bumbers+1);
col = col(round(1: stp: length(col)), :);
close all

iteration_n = 0;
for n_ord = n_ords
    iteration_n = iteration_n + 1;
    h_figure{iteration_n} = figure('rend','painters','pos',[n_ord*350, 50, 300, 300]);
    iteration = 0;
    iteration_s = 0;
    for s = side_shifts
        iteration_s = iteration_s + 1;
        iteration_P = 0;
        for p = possible_cutoffs(possible_cutoffs>=n_ord)
            iteration_P = iteration_P + 1;
            iteration = iteration + 1;
            l = tap_length(find(possible_cutoffs == p));
            %%  filter response
            [c{iteration_n, iteration}] = derivstag(l, p, s, n_ord, sym_flag);
            h_response{iteration_n, iteration} = fftshift(abs(fft((c{iteration_n, iteration}), n_fft)));
            legend_string{iteration_n, iteration} = ['l=', num2str(l), ', P=', num2str(p), ', s=', num2str(s)];
            figure(h_figure{iteration_n})
            plot(omega, h_response{iteration_n, iteration}, '-', 'LineWidth', 1, 'Color', col(iteration,:))
            hold on
        end
    end
    n_observation(iteration_n) = iteration;
end

%%
for iteration_n = 1: numel(n_ords)
    figure(h_figure{n_ords(iteration_n)+1})
    plot(omega, abs(omega).^n_ords(iteration_n), '--r')
    legend_string{iteration_n, n_observation(iteration_n)+1} = ['|(i\omega)^', num2str(n_ords(iteration_n)), '|'];
    range_plot = pi;
    axis([-range_plot, range_plot, 0, 1.025*(range_plot)^n_ords(iteration_n)])
    axis square
    ylabel(['|H^' num2str(n_ords(iteration_n)), '(\omega)|'])
    xlabel('\omega')
    title(['Derivative order n=', num2str(n_ords(iteration_n)), ', centralized'])
    set(gca, 'YScale', 'linear', 'FontSize', 8)
    if n_ords(iteration_n) == 0
        legend_location = 'South';
    else
        legend_location = 'north';
    end
    if iteration_n == numel(n_ords)
        legend(legend_string(iteration_n, 1: (n_observation(iteration_n)+1)), 'Location', legend_location)
    end
end
