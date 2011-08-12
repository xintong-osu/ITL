clear all;
clear java;

% nc_path = '/homes/leeten/projects/nek5000-example/benard/dump/ray_9';
% nc_filename_format = 'ray_9.rank_%d.nc';
% n_ranks = 3;

session_name = 'ray_dn';
dump_path = '~/projects/nek5000-example/benard/dump';  
nc_path = [dump_path '/' session_name];
nc_filename_format = [session_name '.rank_%d.nc'];
n_ranks = 3;

% nc_path = './eddy_uv';
% nc_filename_format = 'eddy_uv.rank_%d.nc';
% n_ranks = 8;

%% open the files
ncids = zeros(1, n_ranks);
for ri = 1:n_ranks
    rid = ri - 1; % from 0-based to 1-based
    ncids(ri) = netcdf.open(sprintf([nc_path, '/', nc_filename_format], rid), 'NC_NOWRITE');
end

%% parse the data from the 1st file
[n_dims, n_vars, n_gatts, tdimid] = netcdf.inq(ncids(1));

dim_keys = cell(1, n_dims);
dim_values = cell(1, n_dims);
for di = 1:n_dims
    dimid = di - 1;
    [dimname, dimlen] = netcdf.inqDim(ncids(1), dimid);
    dim_keys{di} = dimname;
    dim_values{di} = dimid;
end
dim_map = containers.Map(dim_keys, dim_values);

%% read the dimension and #blocks
[xdimid, xlen] = netcdf.inqDim(ncids(1), dim_map('x'));
[ydimid, ylen] = netcdf.inqDim(ncids(1), dim_map('y'));
[bdimid, n_blks] = netcdf.inqDim(ncids(1), dim_map('block'));


var_keys = cell(1, n_vars);
var_values = cell(1, n_vars);
for vi = 1:n_vars
    varid = vi - 1;
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncids(1),varid);
    var_keys{vi} = varname;
    var_values{vi} = varid;
end
var_map = containers.Map(var_keys, var_values);

%% get the time stamps
tstamps = netcdf.getVar(ncids(1), var_map('time'));
n_times = length(tstamps); % skip one time stamp since the 1st one represents the global information

c1d_x_init = cell(1, n_ranks);
c1d_y_init = cell(1, n_ranks);
for ri = 1:n_ranks
    c1d_x_init{ri} = netcdf.getVar(ncids(ri), var_map('x'), [0, 0, 0], [xlen, n_blks, 1]);
    c1d_y_init{ri} = netcdf.getVar(ncids(ri), var_map('y'), [0, 0, 0], [ylen, n_blks, 1]);
end

%% plot the data
figure;

%% prefetch the colormap
color_map = colormap;
n_bins = 16;
max_entropy = log2(n_bins);

% for each time step
for ti = 2:n_times  % skip the 1st time step since it is reserved for the default value
    tid = ti - 1;
    tstamp = tstamps(ti);
    
    %% clear the current frame
    if( tid > 1 )
    	clf;
    end
    
    xmin = +inf;
    ymin = +inf;
    xmax = -inf;
    ymax = -inf;
    % for each rank
    for ri = 1:n_ranks
        %% read the block-wise entropy
        %% NOTE: currently the log is still in the original format
        fid = fopen(sprintf('%s/%s/ge.rank_%d.log', dump_path, session_name, ri - 1), 'rt'); % ri is convered from 1-based to 0-based index
        assert(fid > 0);

        ge = textscan(fid, 'Block %d RV %d Entropy %f');
        ge_rv = find(ge{2}==1);
        blk_entropies = ge{3};
        blk_entropies = blk_entropies(ge_rv);
        blk_entropies = blk_entropies((tid-1)*n_blks + (1:n_blks));
        fclose(fid);
        
        %% read the XY coordinates
        x_time = netcdf.getVar(ncids(ri), var_map('x'), [0, 0, tid], [xlen, n_blks, 1]);
        if( x_time(1) > 1.0e36 )
            x_time = c1d_x_init{ri};
        end
        
        y_time = netcdf.getVar(ncids(ri), var_map('y'), [0, 0, tid], [ylen, n_blks, 1]);
        if( y_time(1) > 1.0e36 )
            y_time = c1d_y_init{ri};
        end
        
        xmin = min(xmin, min(x_time(:)));
        ymin = min(ymin, min(y_time(:)));
        xmax = max(xmax, max(x_time(:)));
        ymax = max(ymax, max(y_time(:)));
        % for each block
        for bi = 1:n_blks
            bid = bi - 1;
            % generate the coordinates in the given time step
            [X, Y] = meshgrid(x_time(:, bi), y_time(:, bi));

            %% read the U component at this time step
            U = netcdf.getVar(ncids(ri), var_map('data_u'), [0, 0, 0, 0, bid, tid], [xlen, ylen, 1, 1, 1, 1]);
            U = U'; % from row-major order to colume-major order
            
            %% read the V component at this time step
            V = netcdf.getVar(ncids(ri), var_map('data_v'), [0, 0, 0, 0, bid, tid], [xlen, ylen, 1, 1, 1, 1]);
            V = V'; % from row-major order to colume-major order

            %% read the temperature at this time step
            T = netcdf.getVar(ncids(ri), var_map('data_temperature'), [0, 0, 0, 0, bid, tid], [xlen, ylen, 1, 1, 1, 1]);
            T = T';
            
            %% map the block entropy to color
            blk_entropy = blk_entropies(bi);
            color_entry = min(max(blk_entropy / max_entropy * size(color_map, 1), 1), size(color_map, 1));
            blk_color = interp1(1:size(color_map, 1), color_map, color_entry);

            %% 
            subplot(2, 3, [1 2]);
            hold on;
            quiver(X, Y, U, V, 'Color', blk_color, 'LineWidth', 2.0);

            subplot(2, 3, [4 5]);
            hold on;
            imagesc(x_time(:, bi), y_time(:, bi), T, [0, 1]);
            
            %% read the random samples
            rs1 = netcdf.getVar(ncids(ri), var_map('rv_vec'), [0, 0, 0, 0, bid, tid], [xlen, ylen, 1, 1, 1, 1]);
            rs2 = netcdf.getVar(ncids(ri), var_map('rv_temperature'), [0, 0, 0, 0, bid, tid], [xlen, ylen, 1, 1, 1, 1]);
            
            subplot(2, 3, [3 6]);
            hold on;
            plot(rs1(:), rs2(:), '.');
        end
    end
    
    subplot(2, 3, [1 2]);
    axis equal;
    xlim([xmin, xmax]);
    ylim([ymin, ymax]);
    title(sprintf('Time %d', tstamp));
    colorbar;

    subplot(2, 3, [4 5]);
    axis equal;
    xlim([xmin, xmax]);
    ylim([ymin, ymax]);
    title('Temperature');
    colorbar;
    
    subplot(2, 3, [3 6]);
    axis equal;
    xlim([0, n_bins]);
    ylim([0, n_bins]);
    xlabel('RS (vec)');
    ylabel('RS (temp.)');
    title(regexprep(sprintf('Joint distribution (%s)', session_name), '_', '\\_'));

%     pause(1.0/15.0);
    saveas(gcf, sprintf([nc_path, '/', session_name, '.t_%d.png'], tstamp));
end

%% open the files
for ri = 1:n_ranks
    netcdf.close(ncids(ri));
end
