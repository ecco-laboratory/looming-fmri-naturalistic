function activations = load_encoding_activations_allsubs(paths_activations, tr_duration)
    addpath('/home/data/eccolab/Code/GitHub/spm12'); % for spm_hrf convolution
    
    activations = [];
    disp('Loading encoding model activation timecourses')
    n_subjs = length(paths_activations);

    fprintf('Current subject:           ')
    for j=1:n_subjs
        fprintf('\b\b\b\b\b\b\b\b\b\b%03d of %03d', j, n_subjs)
        % then run
        for k=1:length(paths_activations{j})
            activations_this_run = readmatrix(paths_activations{j}{k});
            % Convolve features with HRF here
            % holy shit... matlab anonymous functions
            conv_activations = arrayfun(@(i) conv(double(activations_this_run(:, i)), spm_hrf(tr_duration)), 1:size(activations_this_run, 2), 'UniformOutput', false);
            conv_activations = cell2mat(conv_activations);
            % trim off the tail introduced by the convolution
            conv_activations = conv_activations(1:height(activations_this_run), :);
            % now after the run activations have been convolved we can concatenate onto the main
            % should end up with a 2D array where time x subject is on the rows and unit is on the cols
            activations = [activations; conv_activations];
        end
    end
    fprintf('\n')
end
