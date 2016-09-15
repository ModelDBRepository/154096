% Runs Crumiller code on a set of experimental results (i.e. different runs
% of an experiment, all held within subdirectories)

% This code descends into the subdirectories and runs Crumiller's calcinfo2
% on the default set of filenames (info-uniques.csv and info-repeats.csv).
% It also reads in the info-celldeaths.csv file to obtain times-of-death
% and to enable mapping between neurosim ('real') cell IDs and the cell IDs
% emerging from the Crumiller information code.

% Finally, it plots a scatter graph of information-contribution vs time of
% death for the given set of folders

% Usage: alzinfo(filepath[, makeclean])
% Where:   filepath is the path to a directory containing subdirs,
%            each with an info-uniques.csv and info-repeats.csv file
%          makeclean is an optional flag. Set true to force re-generation
%            of all info plots, or false to use existing info-processed.mat
%            files (if any). Defaults to false.

function alzinfo(filepath, varargin)
addpath('info-matlab') % Add info-matlab to path

% See if the 'make clean' flag was set
if length(varargin) > 0
    makeclean = varargin{1};
else
    makeclean = false;
end

% Use current directory if it contains info-uniques and info-repeats, otherwise
% search the subdirectories (this allows us to process either a single
% experiment or a whole directory of multiple experiments)
if exist(strcat(filepath, '/info-processed.mat')) | (exist(strcat(filepath, '/info-uniques.csv')) & exist(strcat(filepath, '/info-repeats.csv')))
    currentDir = filepath;
    subfolders={''};
else
    % Find all subdirectories under 'filepath' (each containing an
    % experimental run)
    d = dir(filepath);
    isub = [d(:).isdir]; %# returns logical vector
    subfolders = {d(isub).name}'
    subfolders(ismember(subfolders,{'.','..'})) = []; % Strip out '.' and '..'
end

allinfos = []; % Store all points for information contribution
alldeaths = []; % Store all points for time-of-death
allscales = []; % Store all points for scale factor at time of death
allpops = []; % Store all population labels
alldelproportions = {[], [], [], [], [], [], [], [], [], [], [], [], []}; % Cell struct containing 
% { { pop1del1 pop1del1 pop1del1; pop1del2 pop1del2 pop1del2 }; { pop2del1 pop2del1 pop2del1; pop2del1 pop2del1 pop2del1 } }
% Access using alldelproportions{pop}(run_no, :)

% Create sliding window for smoothed average plots
windowSize = 10;
h = ones(1,windowSize)/windowSize;

% Define colours for each population
colours = {'b', 'c', 'g', 'y', 'r', 'b', 'c', 'g', 'y', 'r', 'b', 'c', 'g'};

% Run info-matlab/calcinfo2 on the uniques/repeats files to generate the info-processed file
for experiment = 1:length(subfolders)
    currentDir = strcat(filepath, '/', subfolders{experiment});
    if exist(strcat(currentDir, '/info-celldeaths.csv'))

        cellidmap = []; % Make Crumiller info calc variables available
        info_cell = []; % outside the if-else-end block below
        
        if exist(strcat(currentDir, '/info-processed.mat')) & ~makeclean
            fprintf('Loading pre-processed info data from %s\n', currentDir)
            load(strcat(currentDir, '/info-processed.mat'))
        else
            fprintf('Calling Crumiller code for experiment %s\n', currentDir);
            % Call Crumiller info code
            % (This also saves processed data as 'info-processed.mat')
            [cellidmap,~,info_cell,~] = calcinfo2(currentDir);
        end
        
        % Load list of times-of-death
        deathtimescsv = csvread(strcat(currentDir, '/info-celldeaths.csv'));

        numcells = length(cellidmap);
        theseinfos = zeros(1,numcells); % Store points for this run's info contribution
        thesedeaths = zeros(1,numcells); % Store points for this run's cell death times
        thesescales = zeros(1,numcells); % Store points for this run's scale factors
        thesepops = zeros(1,numcells); % Store label of population for each point
        
        
        %% Process information-per-cell
        % Find cumulative sum of all values but bound final value to >= 0
        sum_info_cell = sum(info_cell'); % Sum info across the columns for each cell
        sum_info_cell(sum_info_cell < 0) = 0; % Bound all -ve values to zero
        
        % Find information per population
        cellscale = 1; % Unless we're running a bigger sim...
        poplabels = { 'E6' 'I6' 'I6L' 'E5B' 'E5R' 'I5' 'I5L' 'E4' 'I4' 'I4L' 'E2' 'I2' 'I2L'};
        popsizes = [59 25 13 17 65 25 13 30 20 14 150 25 13]*cellscale; % First element was 60, but cell 0 always seems to be missing
        %colours = [1 3 4 6 8 10 11 12 13 14 15 16 17]; % stretched rainbow
        %colours = [1 0 0 2 2.5 0 0 3 0 0 3.5 0 0]; % When only using E-populations
        
        % sum_info_cell now contains information-per-cell, but its indices
        % don't match up with actual cell IDs from the network simulation.
        % cellidmap provides us with the actual network cell ID for each
        % element in info_cell, so we can find that cell's time of death.
        for j = 1:length(sum_info_cell)
            theseinfos(j) = sum_info_cell(j);
            thesedeaths(j) = deathtimescsv((deathtimescsv(:,1) == cellidmap(j)), 2);
            thesescales(j) = deathtimescsv((deathtimescsv(:,1) == cellidmap(j)), 3);
            thesepops(j) = find(cellidmap(j) <= cumsum(popsizes), 1); % Find population to which this cell belongs (returns first index only)
        end
        thesedeaths = thesedeaths .* 1/1000/60/60/24; % Convert ms -> days
        
        
        %% Plot this run's info vs time scattergraph, coloured by population
        thisruninfofig = figure;
        xlabel('Time of death (days)');
        ylabel('Information contribution (bits/s)');
        hold on;
        
        thisrundelfig = figure;
        xlabel('Time of death (days)');
        ylabel('Proportion of cells deleted');
        hold on;
        
        plotnum = 0; % Keep track of current plot colour
        for pop = unique(thesepops)
            plotnum = plotnum + 1;
            % Extract only the points for this population
            thispopinfos = theseinfos(find(thesepops==pop));
            thispopdeaths = thesedeaths(find(thesepops==pop));
            
            % Sort data points into increasing-time order
            [sortedthispoptimes,i] = sort(thispopdeaths);
            sortedthispopinfos = thispopinfos(i);
           
            % For all times with > 1 corresponding death, replace each entry with avg
            [x,idx,bins] = unique(sortedthispoptimes); % x = death time, idx = unique indices locations, bin = death-time bin number for each cell
            averagedthispopinfos = zeros(length(idx),1);
            stdthispopinfos = zeros(length(idx),1);
            numdeleted = zeros(length(idx),1);
            
            % For each 'bin', get the mean, std, and total num deleted
            for bin = 1:length(idx)
                averagedthispopinfos(bin) = mean(sortedthispopinfos(find(bins==bin)));
                stdthispopinfos(bin) = std(sortedthispopinfos(find(bins==bin)));
                numdeleted(bin) = max(numdeleted) + length(sortedthispopinfos(find(bins==bin)));
            end
            numdeleted = numdeleted ./ popsizes(pop); % Normalise by population size
            % numdeleted gives us a y-axis for the line, but at arbitrary x points
            % So now we can interpolate across numdeleted to generate deletion-proportions at fixed x-intervals
            % which then allow us to find mean/std etc. for an average plot later.
            xi = (0:0.02:2); % Interpolate over ~100 points from 0->2 days
            intp = interp1q(x',numdeleted,xi');
            b = find(isnan(intp)); % Get indices of nans
            intp(b(b<length(intp)/2)) = 0; % Convert NaN->0 at start of line
            intp(b(b>length(intp)/2)) = 1; % Convert NaN->1 at end of line  intp(find(~isnan(intp), 1, 'last')); % Convert NaN->final value at end of line
            alldelproportions{plotnum} = vertcat(alldelproportions{plotnum}, num2cell(intp)');
            % { { pop1del1 pop1del1 pop1del1; pop1del2 pop1del2 pop1del2 }; { pop2del1 pop2del1 pop2del1; pop2del1 pop2del1 pop2del1 } }
            % Access using alldelproportions{pop}(run_no, :)
            
            % Draw the individual points
            figure(thisruninfofig);
            hold on;
            scgrph = scatter(thispopdeaths, thispopinfos, 60, colours{plotnum}, '.');
            set(get(get(scgrph,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off'); % Exclude points from legend
            
            % Draw the average line for this pop
            % Append mean of windowSize first data points to start of the
            % data, to eradicate the rise time.
            risetime = ones(1,windowSize)*mean(averagedthispopinfos(1:windowSize));
            filtereddata = filter(h, 1, [risetime averagedthispopinfos']);
            % Now plot only the part of the filtered data after the initial window size
            plot(x, filtereddata(windowSize+1:end), colours{plotnum});
            
            % Draw (on the other figure) the plot of deletion proportion
            figure(thisrundelfig);
            hold on;
            plot(x, numdeleted, colours{plotnum}, 'linewidth', 2);
        end
        
        % Plot smoothed moving-window average of info vs time
        % From http://stackoverflow.com/questions/1515977/how-to-smoothen-a-plot-in-matlab
        % Sort death times and info of dead cells into ascending order
        [sortedtimes,i] = sort(thesedeaths);
        sortedinfos = theseinfos(i);
        % For all times with > 1 corresponding death, replace each entry with avg
        [d,idx,bins] = unique(sortedtimes); % d = death time, idx = unique indices locations, bin = death-time bin number for each cell
        averagedinfos = zeros(length(idx),1);
        % For each 'bin', find the average info of cells
        for bin = 1:length(idx)
            averagedinfos(bin) = mean(sortedinfos(find(bins==bin)));
        end
        
        % Add average line to info plot.
        % Append mean of windowSize first data points to start of the
        % data, to eradicate the rise time.
        risetime = ones(1,windowSize)*mean(averagedinfos(1:windowSize));
        filtereddata = filter(h, 1, [risetime averagedinfos']);
        figure(thisruninfofig);
        hold on;
        plot(d, filtereddata(windowSize+1:end), 'k', 'linewidth', 2);
        
        % Add legend and save fig
        legend('E6', 'E5b', 'E5a', 'E4', 'E2/3', 'Average', 'Location', 'NorthWest');
        curylim = ylim;
        ylim([-2 curylim(2)]) % Set bottom of plot to -2, but set top to current value
        saveas(thisruninfofig, strcat(currentDir, '/infotime.fig'));
        print(thisruninfofig, '-painters', '-depsc', '-r900', strcat(currentDir, '/infotime'));
        
        % Add legend to the deletion plot too, and save
        figure(thisrundelfig);
        hold on;
        legend('E6', 'E5b', 'E5a', 'E4', 'E2/3', 'Location', 'NorthWest');
        saveas(thisrundelfig, strcat(currentDir, '/popdeltime.fig'));
        print(thisrundelfig, '-painters', '-depsc', '-r900', strcat(currentDir, '/popdeltime'));
        
        %% Plot info vs scalefactor scattergraph
        %figure;
        %plot(thesescales, theseinfos, 'x');
        %xlabel('Scale factor at time of death');
        %ylabel('Information contribution');
        % Save plot in individual run's directory
        %saveas(gcf, strcat(currentDir, '/infoscale.fig'));
        %print(gcf, '-painters', '-depsc', '-r300', strcat(currentDir, '/infoscale'));
        
        close all;
        
        %% Keep total record of info vs time-of-death points
        allinfos = [allinfos theseinfos];
        alldeaths = [alldeaths thesedeaths];
        allscales = [allscales thesescales];
        allpops = [allpops thesepops];
    else
        fprintf('No info-celldeaths.csv file present');
    end
end

%examine_data(filepath) % Run data analysis on the processed information measures


%% Plot average combination of all info vs time-of-death runs

%%% Process the data to find averages over all populations

% Sort data points into increasing-time order
[sortedalltimes,i] = sort(alldeaths);
sortedallinfos = allinfos(i);

% For all times with > 1 corresponding death, replace each entry with avg
[avgx,idx,bins] = unique(sortedalltimes); % x = death time, idx = unique indices locations, bin = death-time bin number for each cell
averagedallinfos = zeros(length(idx),1);
stdallinfos = zeros(length(idx),1);

% For each 'bin', get the mean and std
for bin = 1:length(idx)
    averagedallinfos(bin) = mean(sortedallinfos(find(bins==bin)));
    stdallinfos(bin) = std(sortedallinfos(find(bins==bin)));
end

% Plot line showing average over all populations
figure;
hold on;
risetime = ones(1,windowSize)*mean(averagedallinfos(1:windowSize));
filteredavgdata = filter(h, 1, [risetime averagedallinfos']);
boundedline(avgx, filteredavgdata(windowSize+1:end), stdallinfos, 'k');

% Set up rest of figure
xlabel('Time of death (days)');
ylabel('Information contribution (bits/s)');
curylim = ylim;
ylim([-0.2 curylim(2)]) % Set bottom of plot to 0, but set top to current value


%%% Add separate moving-window average for each population
plotnum = 0; % Keep track of current plot colour
for pop = unique(allpops)
    plotnum = plotnum + 1;
    % Extract only the points for this population
    thispopinfos = allinfos(find(allpops==pop));
    thispopdeaths = alldeaths(find(allpops==pop));
    
    % Sort data points into increasing-time order
    [sortedthispoptimes,i] = sort(thispopdeaths);
    sortedthispopinfos = thispopinfos(i);
    
    % For all times with > 1 corresponding death, replace each entry with avg
    [x,idx,bins] = unique(sortedthispoptimes); % x = death time, idx = unique indices locations, bin = death-time bin number for each cell
    averagedthispopinfos = zeros(length(idx),1);
    stdthispopinfos = zeros(length(idx),1);
    
    % For each 'bin', get the mean and std
    for bin = 1:length(idx)
        averagedthispopinfos(bin) = mean(sortedthispopinfos(find(bins==bin)));
        stdthispopinfos(bin) = std(sortedthispopinfos(find(bins==bin)));
    end
    
    % Draw the line
    %boundedline(x, filter(h, 1, averagedthispopinfos), stdthispopinfos, colourcelllist{pop}, 'alpha');
    risetime = ones(1,windowSize)*mean(averagedthispopinfos(1:windowSize));
    filtereddata = filter(h, 1, [risetime averagedthispopinfos']);
    toplot = filtereddata(windowSize+1:end);
    plot(x, toplot, colours{plotnum}, 'linewidth', 1);
end

% Re-draw average line (but not error patch) on top
toplot = filteredavgdata(windowSize+1:end);
plot(avgx, toplot, 'k', 'linewidth', 2);

% Save plot in parent directory
legend('Std.dev', 'Mean', 'E6', 'E5b', 'E5a', 'E4', 'E2/3', 'Location', 'NorthEast');
saveas(gcf, strcat(currentDir, '/../avginfotime.fig'));
%print(gcf, '-painters', '-depsc', '-r900', strcat(currentDir, '/../avginfotime'));
plot2svg(strcat(currentDir, '/../avginfotime.svg'), gcf, 'png'); % Save as SVG to get transparency
convcommand = sprintf('rsvg-convert -f pdf -w 500 -h 500 -o %s %s', strcat(currentDir, '/../avginfotime.pdf'), strcat(currentDir, '/../avginfotime.svg'));
system(convcommand);

%%% Add individual data points and save another figure (looks messy, not for printing!)
%scatter(alldeaths, allinfos, 60, allcolours, '.');
% Save plot in parent directory
%saveas(gcf, strcat(currentDir, '/../avginfotimepoints.fig'));
%print(gcf, '-painters', '-depsc', '-r900', strcat(currentDir, '/../avginfotimepoints'));


%% Plot average deletion proportion graphs over all the runs
% Store (time, proportion) tuples for each population in separate cell structs
%   (during processing in the loop above).
% For each experiment, append future tuples to the existing population cell struct
% Then sort the tuples on the time axis -- there should be plenty of overlap
% Finally, plot a boundedline for each population cell struct using mean,std

figure;
hold on;

for pop = 1:plotnum
    % Obtain set of deletion vs time data for this pop
    times = (0:0.02:2); % Fixed x-axis
    dels = cell2mat(alldelproportions{pop});
    zerocols = sum(dels)>0;
    dels = dels(:, zerocols); % Keep only non-zero columns
    times = times(:, zerocols);
    
    % Find the moving-window average
    %risetime = ones(1,windowSize)*mean(avgdels(1:windowSize));
    %filtereddata = filter(h, 1, [risetime avgdels']);
    %toplot = filtereddata(windowSize+1:end);
    [blplot, blpatch] = boundedline(times, mean(dels), std(dels), colours{pop}, 'alpha');
    set(get(get(blplot,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off'); % Exclude avg from legend
    set(get(get(blpatch,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off'); % Exclude std from legend
    plot(times, mean(dels), colours{pop}, 'linewidth', 2);
end
ylim([0 1]); % Bound axes to 0,1
xlabel('Time of death (days)');
ylabel('Proportion of cells deleted');
legend('E6', 'E5b', 'E5a', 'E4', 'E2/3', 'Location', 'NorthWest');
saveas(gcf, strcat(currentDir, '/../avgpopdeltime.fig'));
%print(gcf, '-painters', '-depsc', '-r900', strcat(currentDir, '/../avgpopdeltime'));
plot2svg(strcat(currentDir, '/../avgpopdeltime.svg'), gcf, 'png'); % Save as SVG to get transparency
convcommand = sprintf('rsvg-convert -f pdf -w 500 -h 500 -o %s %s', strcat(currentDir, '/../avgpopdeltime.pdf'), strcat(currentDir, '/../avgpopdeltime.svg'));
system(convcommand);

%% Plot combination of all info vs scalefactor runs
%clf;
%plot(allscales, allinfos, 'x');
%xlabel('Scale factor at time of death');
%ylabel('Information contribution');
% Save plot in parent directory
%saveas(gcf, strcat(currentDir, '/../avginfoscale.fig'));
%print(gcf, '-painters', '-dpdf', '-r300', strcat(currentDir, '/../avginfoscale'));

rmpath('info-matlab') % Housekeeping (remove info-matlab from path)
save(strcat(currentDir, '/../avginfo.mat')); % Save all data to file
end
