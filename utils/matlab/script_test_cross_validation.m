clear all;
clear java;

%% load data

% load flujet
% D = X(1:2:end, 1:2:end);
% D = D + (rand(size(D)) - 0.5);

% D = normrnd(32, 16, [sample_size, 1]);
% D = D(0.0 < D & D < 64);

% D = rand(sample_size, 1) * 64;

% D = randi(64, sample_size, 1);
% D = D + rand(sample_size, 1) - 0.5;
% D = D(0.0 < D & D < 64);

% D = poissrnd(32, sample_size, 1);
% D = D + rand(sample_size, 1) - 0.5;
% D = D(0.0 < D & D < 64);

% D = peaks;
% D = 64 * (D - min(D(:))) / (max(D(:)) - min(D(:)));

% Flash
fid = fopen('~/projects/flash-example/run/vn_128_2/dump/flash/feature_vector.rank_0.blk_0.t_1.rv_0.f.b3d', 'rb');
assert(0 < fid);
dim = fread(fid, 3, 'int32');
D = fread(fid, prod(dim), 'single');
fclose(fid);

%% record the final number of samples
sample_size = length(D(:));

%% plot the data
figure;
imagesc(D)
colormap(jet)
colorbar;
axis equal;
axis tight;
axis xy;
title('Input');

%% D_range of the data
D_range = [min(D(:)), max(D(:))];
% D_range(1) = 0; % floor(D_range(1));
% D_range(2) = 64; % ceil(D_range(2));

testing_bins = 2 .^ (3:16);
% testing_bins = 2 .^ (10:16);

%% H_range: range of the histogram
H_max_width = (D_range(2) - D_range(1))/testing_bins(1);
H_range(1) = D_range(1) - H_max_width/2;
H_range(2) = D_range(2) + H_max_width/2;
% H_range = D_range;

x_scale = 'log';
n_testing_bins = length(testing_bins);

%% apply cross-validation
figure;
n_cols = ceil(sqrt(n_testing_bins));
n_rows = ceil(n_testing_bins/n_cols);
% hold on;
cvs = [];
cv_noms = [];
cv_denoms = [];
for bi = 1:n_testing_bins
    n_bins = testing_bins(bi); 
	bin_width = (D_range(2) - D_range(1))/n_bins;
    
% 	bin_edges = H_range(1):bin_width:H_range(2);
	bin_edges = D_range(1)-bin_width/2:bin_width:D_range(2)+bin_width/2;
    bin_centers = (bin_edges(2:end) + bin_edges(1:end-1))/2;
    
    %% 
    D_hist = histc(D(:), bin_edges);
    D_hist = [D_hist(1:end-2); D_hist(end-1) + D_hist(end)];
    assert(sample_size == sum(D_hist));
    
    D_pmf = D_hist / sample_size;
    if( 1 == bi )
        y_max = max(D_pmf);
    end
 
    cv_nom = 2 - sum(D_pmf.^2) * (sample_size+1);
    cv_denom = (sample_size-1) * bin_width;
    cvs(end+1) = cv_nom / cv_denom;
    cv_noms(end+1) = cv_nom;
    cv_denoms(end+1) = cv_denom;
    
    %% plot the average histogram
    subplot(n_rows, n_cols, bi);
    hold on;
	bar(bin_centers, D_hist, 1.0, 'EdgeColor', 'none');
%     area(bin_edges, [D_pmf(:); D_pmf(end)], 'EdgeColor', 'none');
    title(sprintf('#Bins = %d, Cross-Validation = %.2e', n_bins, cvs(bi)));
    xlim(H_range);
    ylim([0, max(D_hist)]);
%     set(gca, 'YScale', 'log');
end

figure;
for ci = 1:3
    switch(ci)
        case 1
            subplot(3, 2, [1 2 3 4]);
            plot(testing_bins, cvs, '-o');
            title('Cross-Validation');
        case 2
            subplot(3, 2, 5);
            plot(testing_bins, cv_noms, '-o');
            title('Numerator');
        case 3
            subplot(3, 2, 6);
            plot(testing_bins, cv_denoms, '-o');
            title('Denominators');
    end
    set(gca, 'XScale', x_scale);
    set(gca, 'XTick', testing_bins);
    xlabel('# bins');
%     xlim([testing_bins(1), testing_bins(end)]);
    axis tight;
end

