% this script plot the vector field dumpped from the ITL wrapper

clear all;
clear java;

% dump_path = '~/projects/nek5000-example/eddy_uv/dump/eddy_uv';
% n_ranks = 8;
% n_block_per_rank  = 32;
% tstamps = 1:200:1001;
% rv = 0; % index of the random variable
% vec_rv = 0;

dump_path = '~/projects/nek5000-example/benard/dump/ray_dd';
n_ranks = 3;
n_block_per_rank  = 1;
tstamps = 1:10:201;
rv = 1; % index of the random variable
vec_rv = 1; % index of the rv for the vector field

color_map = colormap;
n_bins = 360;
max_entropy = 4;
% max_entropy = log2(n_bins);
n_tstamps = length(tstamps);

figure;
for ti = 1:n_tstamps
    ts = tstamps(ti);
    if( ti > 1 )
    	clf;
    end

    
    %% plot the vector as glyph and the coordinates as points
    X_all = []; % x coord.
    Y_all = []; % y coord.
    U_all = []; % u-component of the vector
    V_all = []; % v-component of the vector
    S_all = []; % random samples
    H_blocks = []; % entropy
    for ri = 1:n_ranks
        fid = fopen(sprintf('%s/ge.rank_%d.log', dump_path, ri - 1), 'rt'); % ri is convered from 1-based to 0-based index
        assert(fid > 0);

        ge = textscan(fid, 'Block %d RV %d Entropy %f');
        ge_rv = find(ge{2}==rv);
        block_entropies = ge{3};
        block_entropies = block_entropies(ge_rv);
%         block_entropies = block_entropies(end-n_block_per_rank+1:end);
        block_entropies = block_entropies((ti-1)*n_block_per_rank+(1:n_block_per_rank));
        n_blocks_per_rank = length(block_entropies(:));
    %     block_entropies = fread(fid, n_blocks_per_rank, 'single');
        fclose(fid);
%         H_blocks = [H_blocks; block_entropies(:)];

        for bi = 1:n_blocks_per_rank
            block_entropy = block_entropies(bi);
            color_entry = min(max(block_entropy / max_entropy * size(color_map, 1), 1), size(color_map, 1));
            block_color = interp1(1:size(color_map, 1), color_map, color_entry);

            %% read the dumpped geometry
            fid = fopen(sprintf('%s/geom.rank_%d.blk_%d.t_0.txt', dump_path, ri - 1, bi - 1), 'rt'); % ri is convered from 1-based to 0-based index
            assert(fid > 0);
            n_xs = fscanf(fid, '%d', 1);
            xs = textscan(fid, '%f %f', n_xs);
            xs = xs{1};
            x_range = fgets(fid);
            x_range = fgets(fid);
            n_ys = fscanf(fid, '%d', 1);
            ys = textscan(fid, '%f %f', n_ys);
            ys = ys{1};
            y_range = fgets(fid);
            y_range = fgets(fid);
            n_zs = fscanf(fid, '%d', 1);
            zs = textscan(fid, '%f %f', n_zs);
            zs = zs{1};
            z_range = fgets(fid);
            z_range = fgets(fid);
            fclose(fid);

            [X, Y, Z] = meshgrid(xs, ys, zs);

            if( vec_rv == rv )
                %% read the dumpped vector field
                fid = fopen(sprintf('%s/feature_vector.rank_%d.blk_%d.t_%d.rv_%d.7VECTOR3.b3d', dump_path, ri - 1, bi - 1, ts, vec_rv), 'rb'); % ri is convered from 1-based to 0-based index
                assert(fid > 0);
                dim = fread(fid, 3, 'int32');
                vec = fread(fid, 3 * prod(dim), 'single'); % 3 for VECTOR3
                fclose(fid);
                vec = reshape(vec, [3 dim(:)']); % the 1st '3' is because the element is VECTOR3
                U =     squeeze(vec(1, :, :))';
                V = squeeze(vec(2, :, :))';

                %% downsample the data to speed up the rendering
                down_step = 1;
                down_start = round(down_step/2);
                X_down = X(down_start:down_step:end, down_start:down_step:end);
                Y_down = Y(down_start:down_step:end, down_start:down_step:end);
                U_down = U(down_start:down_step:end, down_start:down_step:end);
                V_down = V(down_start:down_step:end, down_start:down_step:end);
                H_down = block_entropy * ones(size(V));

                %% plot the vector field
                hold on;
                quiver(X_down(:), Y_down(:), U_down(:), V_down(:), 'Color', 'black', 'LineWidth', 4.0);
                quiver(X_down(:), Y_down(:), U_down(:), V_down(:), 'Color', block_color, 'LineWidth', 2.0);

                X_all = [X_all; X_down];
                Y_all = [Y_all; Y_down];
                U_all = [U_all; U_down];
                V_all = [V_all; V_down];
            else
                %% read the dumpped random samples
                fid = fopen(sprintf('%s/feature_vector.rank_%d.block_%d.rv_%d.f.b3d', dump_path, ri - 1, bi - 1, rv), 'rb'); % ri is convered from 1-based to 0-based index
                assert(fid > 0);
                dim = fread(fid, 3, 'int32');
                vec = fread(fid, prod(dim), 'single'); % 3 for VECTOR3
                fclose(fid);
                S = reshape(vec, dim(:)'); % the 1st '3' is because the element is VECTOR3

                S_all = [S_all; S(:)]; % random samples
            end
        end
    end

    axis equal;
    axis off;
    title(sprintf('Entropy-colored Vectors (%d)', ts));
    
    pause(1);
end

% figure;
% hist(S_all(~isnan(S_all) & ~isinf(S_all)), 32);
% 
    % subplot(1, 2, 1);
    % title('Points');
    % 
    % subplot(1, 2, 2);
    % title('Vectors');

% %% plot the histogram
% figure;
% hist(H_blocks, 16);
% title('Entropy Histogram');
% xlabel('Entropy');
% ylabel('Frequencies');
% xlim([0, max_entropy]);
%     